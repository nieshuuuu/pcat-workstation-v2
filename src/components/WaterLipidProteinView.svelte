<script lang="ts">
  /**
   * Whole-volume water/lipid/PROTEIN (3-material) decomposition viewer.
   *
   * Reproduces the delivered maps of wlp-decomposition/wlp_decomposition.jl:
   * every soft-tissue voxel of an axial slice decoded through the frozen poly2
   * calibration surface into (f_w, f_l, f_p), σ_f-weighted coupled Huber-TV,
   * jet 0..1 over a grayscale CT, σ_f in viridis, NaN (gated) transparent,
   * anterior-up.
   *
   * The surface is FROZEN (sim-calibrated at 70/150 keV), not self-calibrated —
   * see the caveat banner. Sibling of the 2-material Water/Lipid tab (its own
   * self-calibrated GLS view); the pericoronary ROI overlay also uses that GLS.
   *
   * Layout:
   *   ┌───────────────────────────────┬──────────────────────┐
   *   │                               │ Map: [f_w][f_l][f_p][σ]│
   *   │   axial slice canvas          │ ☐ CT only             │
   *   │   (CT + jet overlay)          │ colorbar              │
   *   │                               │ model readout + caveat │
   *   ├───────────────────────────────┴──────────────────────┤
   *   │  z = 179 / 339   [────●─────────]                     │
   *   └───────────────────────────────────────────────────────┘
   */
  import {
    runWaterLipidProtein,
    getWlpSlice,
    type WlpModel,
    type WlpSlice,
  } from '$lib/api';
  import { wlpStore } from '$lib/stores/wlpStore.svelte';
  import { jet, CT_LO, CT_HI } from '$lib/colormap';
  import { volumeStore } from '$lib/stores/volumeStore.svelte';
  import { saveSession } from '$lib/session';

  type WlMap = 'fw' | 'fl' | 'fp' | 'sf';

  /* ── State ──────────────────────────────────────────────── */

  // Single source of truth is wlpStore, so a model restored from a saved session
  // (or set by a Run) drives this view and is captured by session save.
  let model = $derived(wlpStore.model);
  let running = $state(false);
  let error = $state('');

  let slice = $state<WlpSlice | null>(null);
  let z = $state(0);
  let map = $state<WlMap>('fl');
  let ctOnly = $state(false);
  // TV smoothing strength (λ), on a median-normalized weight so it's noise-scale
  // invariant. 0 = raw per-voxel (salt-and-pepper), ~0.7 keeps rod/boundary shape
  // while cleaning flat tissue, ≳1.2 over-smooths (boundaries bleed, blotchy).
  let smoothing = $state(0.7);

  /** Latest-wins guard so slider scrubbing doesn't render a stale slice. */
  let fetchSeq = 0;
  /** Single-in-flight coalescing: at most one slice request is ever in flight.
   *  If z changes while one is running, we fetch exactly once more for the
   *  latest value when it returns. */
  let fetchInflight = false;
  let fetchPending = false;

  let canvasEl: HTMLCanvasElement | undefined = $state();
  // x,y = voxel coords; cx,cy = cursor position within the canvas pane (for the
  // floating tooltip). fw/fl/fp are the fractions, sf the std at (x,y).
  let hover = $state<{
    x: number; y: number; hu: number;
    fw: number; fl: number; fp: number; sf: number;
    cx: number; cy: number;
  } | null>(null);

  let nz = $derived(model?.dims[0] ?? 0);

  /* ── Colormaps ──────────────────────────────────────────── */

  // Viridis control points for the σ_f panel.
  const VIRIDIS: [number, number, number][] = [
    [68, 1, 84],
    [59, 82, 139],
    [33, 145, 140],
    [94, 201, 98],
    [253, 231, 37],
  ];
  function viridis(t: number): [number, number, number] {
    const u = Math.max(0, Math.min(1, t));
    const s = u * (VIRIDIS.length - 1);
    const i = Math.min(VIRIDIS.length - 2, Math.floor(s));
    const f = s - i;
    const a = VIRIDIS[i];
    const b = VIRIDIS[i + 1];
    return [
      Math.round(a[0] + f * (b[0] - a[0])),
      Math.round(a[1] + f * (b[1] - a[1])),
      Math.round(a[2] + f * (b[2] - a[2])),
    ];
  }

  /** Map a value (with the active map's colormap + range) to RGB. */
  function mapColor(v: number): [number, number, number] {
    if (map === 'sf') return viridis(v / sfCap);
    return jet(v); // fw / fl / fp are already in [0,1]
  }

  /* ── σ_f display cap (99th percentile of the current slice) ── */
  let sfCap = $state(0.5);
  function computeSfCap(s: WlpSlice): number {
    const vals: number[] = [];
    for (let i = 0; i < s.sf.length; i++) {
      const v = s.sf[i];
      if (Number.isFinite(v)) vals.push(v);
    }
    if (vals.length === 0) return 0.5;
    vals.sort((a, b) => a - b);
    const p99 = vals[Math.min(vals.length - 1, Math.floor(vals.length * 0.99))];
    return Math.max(0.05, p99);
  }

  /* ── Run / fetch ────────────────────────────────────────── */

  async function handleRun() {
    if (running) return;
    running = true;
    error = '';
    try {
      // model is derived from wlpStore; setting the store updates the view AND
      // makes the model part of the unified session save.
      wlpStore.set(await runWaterLipidProtein());
      if (volumeStore.dicomPath) {
        try {
          await saveSession(volumeStore.dicomPath);
        } catch (e) {
          console.error('auto-save after water/lipid/protein failed:', e);
        }
      }
    } catch (e) {
      error = e instanceof Error ? e.message : String(e);
      wlpStore.reset();
    } finally {
      running = false;
    }
  }

  async function fetchSlice() {
    if (!model) return;
    // Coalesce: never run two requests at once (see WL viewer note).
    if (fetchInflight) {
      fetchSeq++;
      fetchPending = true;
      return;
    }
    fetchInflight = true;
    const seq = ++fetchSeq;
    try {
      const s = await getWlpSlice(z, smoothing);
      if (seq === fetchSeq) {
        slice = s;
        sfCap = computeSfCap(s);
      }
    } catch (e) {
      if (seq === fetchSeq) error = e instanceof Error ? e.message : String(e);
    } finally {
      fetchInflight = false;
      if (fetchPending) {
        fetchPending = false;
        fetchSlice(); // chase the latest slider position
      }
    }
  }

  // When a new model appears (fresh Run, or one restored from a saved session),
  // jump to the mid-volume slice (heart). Track identity so unrelated reactive
  // ticks don't reset the user's slice position.
  let lastModelForZ: WlpModel | null = null;
  $effect(() => {
    if (model && model !== lastModelForZ) {
      z = Math.floor((model.dims[0] - 1) / 2);
    }
    lastModelForZ = model;
  });

  // Refetch when the model (run/re-run), slice index, or smoothing changes. The
  // selected map and CT-only toggle only re-render — no refetch.
  $effect(() => {
    void z;
    void smoothing;
    if (model) fetchSlice();
  });

  // Re-render whenever the decoded slice, the selected map, or the window changes.
  $effect(() => {
    void slice;
    void map;
    void ctOnly;
    void sfCap;
    render();
  });

  /* ── Rendering ──────────────────────────────────────────── */

  function valueAt(s: WlpSlice, i: number): number {
    if (map === 'fw') return s.fw[i];
    if (map === 'fl') return s.fl[i];
    if (map === 'fp') return s.fp[i];
    return s.sf[i];
  }

  function render() {
    const canvas = canvasEl;
    const s = slice;
    if (!canvas || !s) return;
    const { nx, ny } = s;

    // Offscreen image at native resolution, then scale to the display canvas.
    const off = document.createElement('canvas');
    off.width = nx;
    off.height = ny;
    const octx = off.getContext('2d');
    if (!octx) return;
    const img = octx.createImageData(nx, ny);

    const ctSpan = CT_HI - CT_LO;
    for (let i = 0; i < nx * ny; i++) {
      const hu = s.ct[i];
      const gray = Math.max(0, Math.min(255, Math.round(((hu - CT_LO) / ctSpan) * 255)));

      const v = ctOnly ? NaN : valueAt(s, i);
      let r: number, g: number, b: number;
      if (!ctOnly && Number.isFinite(v)) {
        [r, g, b] = mapColor(v);
      } else {
        r = g = b = gray;
      }
      img.data[i * 4] = r;
      img.data[i * 4 + 1] = g;
      img.data[i * 4 + 2] = b;
      img.data[i * 4 + 3] = 255;
    }
    octx.putImageData(img, 0, 0);

    // Fit the display canvas to its box, preserving aspect ratio. Row 0 is
    // anterior for standard head-first-supine axial data (anterior-up, no flip).
    const box = canvas.parentElement;
    const maxW = box ? box.clientWidth : nx;
    const maxH = box ? box.clientHeight : ny;
    const scale = Math.max(1, Math.min(maxW / nx, maxH / ny));
    canvas.width = Math.round(nx * scale);
    canvas.height = Math.round(ny * scale);

    const ctx = canvas.getContext('2d');
    if (!ctx) return;
    ctx.imageSmoothingEnabled = true;
    ctx.imageSmoothingQuality = 'high';
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    ctx.drawImage(off, 0, 0, nx, ny, 0, 0, canvas.width, canvas.height);
  }

  function handleWheel(e: WheelEvent) {
    if (!model) return;
    e.preventDefault();
    const next = z + (e.deltaY > 0 ? 1 : -1);
    z = Math.max(0, Math.min(nz - 1, next));
  }

  function handleMove(e: MouseEvent) {
    const canvas = canvasEl;
    const s = slice;
    if (!canvas || !s) {
      hover = null;
      return;
    }
    const rect = canvas.getBoundingClientRect();
    const px = Math.floor(((e.clientX - rect.left) / rect.width) * s.nx);
    const py = Math.floor(((e.clientY - rect.top) / rect.height) * s.ny);
    if (px < 0 || px >= s.nx || py < 0 || py >= s.ny) {
      hover = null;
      return;
    }
    const i = py * s.nx + px;
    const pane = canvas.parentElement;
    const prect = pane ? pane.getBoundingClientRect() : rect;
    hover = {
      x: px,
      y: py,
      hu: s.ct[i],
      fw: s.fw[i],
      fl: s.fl[i],
      fp: s.fp[i],
      sf: s.sf[i],
      cx: e.clientX - prect.left,
      cy: e.clientY - prect.top,
    };
  }

  /* ── Colorbar helpers ───────────────────────────────────── */

  let barRange = $derived<[number, number]>(map === 'sf' ? [0, sfCap] : [0, 1]);
  let barLabel = $derived(
    map === 'fw' ? 'f_w  (water fraction)'
      : map === 'fl' ? 'f_l  (lipid fraction)'
      : map === 'fp' ? 'f_p  (protein fraction)'
      : 'σ_f  (per-voxel std)',
  );
  // 64-stop gradient for the colorbar, using the active colormap.
  let barStops = $derived(
    Array.from({ length: 65 }, (_, k) => {
      const t = k / 64;
      const [r, g, b] = map === 'sf' ? viridis(t) : jet(t);
      return `rgb(${r},${g},${b}) ${(t * 100).toFixed(1)}%`;
    }).join(', '),
  );

  function fmt(n: number, d = 1): string {
    return Number.isFinite(n) ? n.toFixed(d) : '—';
  }
</script>

<div class="flex h-full w-full flex-col overflow-hidden bg-surface">
  <!-- Top bar -->
  <div class="flex shrink-0 items-center gap-3 border-b border-border bg-surface-secondary px-3 py-1.5">
    <span class="text-xs font-semibold text-text-primary">Water / Lipid / Protein (poly2 surface)</span>
    {#if model}
      <span class="text-[11px] text-text-secondary tabular-nums">
        {model.low_kev}+{model.high_kev} keV · {model.is_baked ? 'sim-calibrated (frozen)' : 'refit'}
      </span>
    {/if}
    <button
      class="ml-auto rounded bg-accent/15 px-3 py-1 text-xs font-medium text-accent hover:bg-accent/25 active:bg-accent/35 disabled:bg-surface-tertiary/40 disabled:text-text-secondary/70"
      onclick={handleRun}
      disabled={running}
      title="Register the frozen 70/150 keV poly2 surface and decompose the whole volume into water/lipid/protein"
    >
      {running ? 'Decomposing…' : model ? 'Re-run' : 'Run Water/Lipid/Protein Decomposition'}
    </button>
  </div>

  {#if error}
    <div class="shrink-0 border-b border-border bg-error/10 px-3 py-1 text-[11px] text-error">{error}</div>
  {/if}

  {#if !model}
    <div class="flex flex-1 items-center justify-center px-6 text-center">
      <div class="max-w-md text-xs leading-relaxed text-text-secondary">
        Load a dual-energy <span class="font-medium">70 keV + 150 keV</span> MonoPlus series, then
        <span class="font-medium text-accent">Run</span>. Every soft-tissue voxel is decoded through the
        frozen poly2 calibration surface into water / lipid / protein volume fractions (σ_f-weighted
        edge-preserving TV), rendered jet 0..1 over the CT.
        <br /><br />
        <span class="text-amber-500/90">The surface is calibrated on simulation (70/150 keV); on real
        scanner data the fractions are a first estimate — validate or refit per scanner.</span>
      </div>
    </div>
  {:else}
    <div class="flex min-h-0 flex-1 overflow-hidden">
      <!-- Canvas pane -->
      <!-- svelte-ignore a11y_no_static_element_interactions -->
      <div
        class="relative flex min-w-0 flex-1 items-center justify-center overflow-hidden bg-black"
        onwheel={handleWheel}
      >
        <canvas
          bind:this={canvasEl}
          class="max-h-full max-w-full"
          onmousemove={handleMove}
          onmouseleave={() => (hover = null)}
        ></canvas>

        <!-- Slice label -->
        <div class="pointer-events-none absolute left-2 top-2 rounded bg-black/50 px-2 py-0.5 text-[11px] text-yellow-300 tabular-nums">
          z = {z} / {nz - 1}
        </div>

        <!-- Hover readout: floating tooltip at the cursor with the estimated
             volume fractions. σ_f is a shared per-voxel uncertainty proxy. -->
        {#if hover}
          <div
            class="pointer-events-none absolute z-10 rounded bg-black/75 px-2 py-1 text-[10px] leading-snug text-white tabular-nums shadow-lg"
            style="left: {hover.cx + 14}px; top: {hover.cy + 14}px;"
          >
            {#if Number.isFinite(hover.fl)}
              <div>f_w&nbsp;&nbsp;{fmt(hover.fw, 2)}</div>
              <div>f_l&nbsp;&nbsp;{fmt(hover.fl, 2)}</div>
              <div>f_p&nbsp;&nbsp;{fmt(hover.fp, 2)}</div>
              <div class="text-white/55">σ_f {fmt(hover.sf, 2)} · HU {Math.round(hover.hu)} · ({hover.x}, {hover.y})</div>
            {:else}
              <div class="text-white/55">HU {Math.round(hover.hu)} · gated · ({hover.x}, {hover.y})</div>
            {/if}
          </div>
        {/if}
      </div>

      <!-- Controls + colorbar + model readout -->
      <div class="flex w-64 shrink-0 flex-col gap-3 overflow-y-auto border-l border-border bg-surface-secondary px-3 py-3">
        <!-- Map toggle -->
        <div>
          <div class="mb-1 text-[10px] font-semibold uppercase tracking-wider text-text-secondary/60">Map</div>
          <div class="flex gap-1">
            {#each [['fw', 'f_w'], ['fl', 'f_l'], ['fp', 'f_p'], ['sf', 'σ_f']] as opt}
              <button
                class="flex-1 rounded px-2 py-1 text-[11px] font-medium transition-colors {map === opt[0]
                  ? 'bg-accent text-white'
                  : 'bg-surface-tertiary text-text-secondary hover:text-text-primary'}"
                onclick={() => (map = opt[0] as WlMap)}
              >
                {opt[1]}
              </button>
            {/each}
          </div>
        </div>

        <!-- Smoothing (TV strength) -->
        <div>
          <div class="mb-1 flex items-center justify-between text-[10px] font-semibold uppercase tracking-wider text-text-secondary/60">
            <span>Smoothing (TV λ)</span>
            <span class="tabular-nums text-text-secondary">{smoothing === 0 ? 'raw' : smoothing.toFixed(1)}</span>
          </div>
          <input
            type="range"
            class="w-full accent-accent"
            min="0"
            max="3"
            step="0.1"
            bind:value={smoothing}
            title="0 = raw per-voxel; ~0.7 keeps boundaries while cleaning flat tissue; ≳1.2 over-smooths. Noise-scale invariant."
          />
        </div>

        <!-- CT only -->
        <label class="flex items-center gap-2 text-[11px] text-text-secondary">
          <input type="checkbox" bind:checked={ctOnly} class="accent-accent" />
          CT only (no overlay)
        </label>

        <!-- Colorbar -->
        {#if !ctOnly}
          <div>
            <div class="mb-1 text-[10px] text-text-secondary">{barLabel}</div>
            <div class="h-3 w-full rounded" style="background: linear-gradient(to right, {barStops});"></div>
            <div class="flex justify-between text-[10px] text-text-secondary tabular-nums">
              <span>{barRange[0].toFixed(map === 'sf' ? 2 : 1)}</span>
              <span>{barRange[1].toFixed(map === 'sf' ? 2 : 1)}</span>
            </div>
          </div>
        {/if}

        <!-- Sim-calibration caveat -->
        {#if model.is_baked}
          <div class="rounded border border-amber-500/30 bg-amber-500/5 px-2 py-1.5 text-[10px] leading-snug text-amber-500/90">
            <span class="font-semibold">Surface</span> frozen, sim-calibrated ({model.low_kev}/{model.high_kev} keV) —
            fractions are a first estimate, CCC is a sim number, refit per scanner.
            <span class="font-semibold">Noise</span> {model.noise_measured ? 'self-calibrated from this volume ✓' : 'sim (could not measure — too little soft tissue)'}.
          </div>
        {/if}

        <!-- Model readout (endpoints + noise) -->
        <div class="mt-1 border-t border-border pt-2 text-[10px] leading-relaxed text-text-secondary">
          <div class="mb-1 font-semibold uppercase tracking-wider text-text-secondary/60">poly2 surface model</div>
          <div class="grid grid-cols-[auto_1fr] gap-x-2 tabular-nums">
            <span>HU water</span><span class="text-text-primary">({fmt(model.hu_w[0])}, {fmt(model.hu_w[1])})</span>
            <span>HU lipid</span><span class="text-text-primary">({fmt(model.hu_l[0])}, {fmt(model.hu_l[1])})</span>
            <span>HU protein</span><span class="text-text-primary">({fmt(model.hu_p[0])}, {fmt(model.hu_p[1])})</span>
            <span>σ {model.low_kev}keV</span><span class="text-text-primary">{fmt(Math.sqrt(model.sigma_hu[0][0]))} {model.noise_measured ? '(meas)' : '(sim)'}</span>
            <span>σ {model.high_kev}keV</span><span class="text-text-primary">{fmt(Math.sqrt(model.sigma_hu[1][1]))} {model.noise_measured ? '(meas)' : '(sim)'}</span>
            <span>ρ</span><span class="text-text-primary">{model.rho.toFixed(3)}</span>
            <span>gate (150keV)</span><span class="text-text-primary">[{fmt(model.gate[0])}, {fmt(model.gate[1])}]</span>
          </div>
        </div>
      </div>
    </div>

    <!-- Slice slider -->
    <div class="flex shrink-0 items-center gap-3 border-t border-border bg-surface-secondary px-3 py-2">
      <span class="text-[11px] text-text-secondary tabular-nums">z {z}/{nz - 1}</span>
      <input
        type="range"
        class="flex-1 accent-accent"
        min="0"
        max={Math.max(0, nz - 1)}
        bind:value={z}
      />
    </div>
  {/if}
</div>

<script lang="ts">
  /**
   * Whole-volume noise-aware GLS water/lipid decomposition viewer.
   *
   * Reproduces fw_fl_maps_theolipid_57955439.png from wl-noise-aware-mmd:
   * self-calibrated from the patient's own fat/muscle, every soft-tissue voxel
   * of an axial slice projected onto the water→lipid line, jet f_w/f_l (0..1)
   * over a grayscale CT, σ_f in viridis, NaN (gated) transparent, anterior-up.
   *
   * Layout:
   *   ┌───────────────────────────────┬──────────────────────┐
   *   │                               │ Map: [f_w][f_l][σ_f]  │
   *   │   axial slice canvas          │ Anchor: [adi][theo]   │
   *   │   (CT + jet overlay)          │ ☐ CT only             │
   *   │                               │ colorbar              │
   *   │                               │ calibration readout   │
   *   ├───────────────────────────────┴──────────────────────┤
   *   │  z = 179 / 339   [────●─────────]                     │
   *   └───────────────────────────────────────────────────────┘
   */
  import {
    runWaterLipid,
    getWlSlice,
    type WlCalibration,
    type WlSlice,
    type WlAnchor,
  } from '$lib/api';

  type WlMap = 'fw' | 'fl' | 'sf';

  /* ── State ──────────────────────────────────────────────── */

  let calib = $state<WlCalibration | null>(null);
  let running = $state(false);
  let error = $state('');

  let slice = $state<WlSlice | null>(null);
  let z = $state(0);
  let map = $state<WlMap>('fw');
  // Default to the theoretical pure-lipid anchor so the maps match the
  // fw_fl_maps_theolipid reference; flip to Adipose for the fat-referenced scale.
  let anchor = $state<WlAnchor>('theoretical');
  let ctOnly = $state(false);

  /** Latest-wins guard so slider scrubbing doesn't render a stale slice. */
  let fetchSeq = 0;
  /** Single-in-flight coalescing: at most one slice request is ever in flight.
   *  If z/anchor change while one is running, we fetch exactly once more for
   *  the latest value when it returns — so a fast drag issues a handful of
   *  requests that converge on the final slice, not one IPC call per tick. */
  let fetchInflight = false;
  let fetchPending = false;

  let canvasEl: HTMLCanvasElement | undefined = $state();
  // x,y = voxel coords; cx,cy = cursor position within the canvas pane (for the
  // floating tooltip). fw/sf are the estimate + its standard deviation at (x,y).
  let hover = $state<{ x: number; y: number; hu: number; fw: number; sf: number; cx: number; cy: number } | null>(null);

  // CT window matching the reference figure's grayscale underlay.
  const CT_LO = -160;
  const CT_HI = 240;

  let nz = $derived(calib?.dims[0] ?? 0);

  /* ── Colormaps ──────────────────────────────────────────── */

  /** Classic jet (matplotlib): t∈[0,1] → blue→cyan→green→yellow→red. Matches
   *  the reference figure's f_w/f_l panels. */
  function jet(t: number): [number, number, number] {
    const u = Math.max(0, Math.min(1, t));
    const r = Math.max(0, Math.min(1, 1.5 - Math.abs(4 * u - 3)));
    const g = Math.max(0, Math.min(1, 1.5 - Math.abs(4 * u - 2)));
    const b = Math.max(0, Math.min(1, 1.5 - Math.abs(4 * u - 1)));
    return [Math.round(r * 255), Math.round(g * 255), Math.round(b * 255)];
  }

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
    return jet(v); // fw / fl are already in [0,1]
  }

  /* ── σ_f display cap (99th percentile of the current slice) ── */
  let sfCap = $state(0.5);
  function computeSfCap(s: WlSlice): number {
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
      calib = await runWaterLipid();
      z = Math.floor((calib.dims[0] - 1) / 2); // start at mid-volume (heart)
      // The fetch effect (tracks calib + z + anchor) loads the slice.
    } catch (e) {
      error = e instanceof Error ? e.message : String(e);
      calib = null;
    } finally {
      running = false;
    }
  }

  async function fetchSlice() {
    if (!calib) return;
    // Coalesce: never run two requests at once. Mark pending and return; the
    // in-flight call re-fires for the latest z/anchor when it settles.
    if (fetchInflight) {
      fetchPending = true;
      return;
    }
    fetchInflight = true;
    const seq = ++fetchSeq;
    try {
      const s = await getWlSlice(z, anchor);
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

  // Refetch when the calibration (run/re-run), slice index or anchor changes.
  // The selected map and CT-only toggle only re-render — no refetch.
  $effect(() => {
    void z;
    void anchor;
    if (calib) fetchSlice();
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

  function valueAt(s: WlSlice, i: number): number {
    if (map === 'fw') return s.fw[i];
    if (map === 'fl') return 1 - s.fw[i]; // f_l = 1 − f_w (NaN propagates)
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
    // anterior for standard head-first-supine axial data, so drawing row 0 at
    // the top keeps the spine down (radiological, anterior-up) — matching the
    // reference figure with no flip.
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
    if (!calib) return;
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
    // Position the tooltip relative to the canvas pane so it follows the cursor.
    const pane = canvas.parentElement;
    const prect = pane ? pane.getBoundingClientRect() : rect;
    hover = {
      x: px,
      y: py,
      hu: s.ct[i],
      fw: s.fw[i],
      sf: s.sf[i],
      cx: e.clientX - prect.left,
      cy: e.clientY - prect.top,
    };
  }

  /* ── Colorbar helpers ───────────────────────────────────── */

  let barRange = $derived<[number, number]>(map === 'sf' ? [0, sfCap] : [0, 1]);
  let barLabel = $derived(
    map === 'fw' ? 'f_w  (0 = lipid → 1 = water)'
      : map === 'fl' ? 'f_l  (0 = water → 1 = lipid)'
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
    <span class="text-xs font-semibold text-text-primary">Water / Lipid (noise-aware GLS)</span>
    {#if calib}
      <span class="text-[11px] text-text-secondary tabular-nums">
        {calib.low_kev}+{calib.high_kev} keV · self-calibrated · {(calib.fraction_kept * 100).toFixed(1)}% soft-tissue
      </span>
    {/if}
    <button
      class="ml-auto rounded bg-accent/15 px-3 py-1 text-xs font-medium text-accent hover:bg-accent/25 active:bg-accent/35 disabled:bg-surface-tertiary/40 disabled:text-text-secondary/70"
      onclick={handleRun}
      disabled={running}
      title="Self-calibrate from the patient's own fat/muscle and decompose the whole volume"
    >
      {running ? 'Decomposing…' : calib ? 'Re-run' : 'Run Water/Lipid Decomposition'}
    </button>
  </div>

  {#if error}
    <div class="shrink-0 border-b border-border bg-error/10 px-3 py-1 text-[11px] text-error">{error}</div>
  {/if}

  {#if !calib}
    <div class="flex flex-1 items-center justify-center px-6 text-center">
      <div class="max-w-md text-xs leading-relaxed text-text-secondary">
        Load a dual-energy (two-keV) series — e.g. MonoPlus 70 keV + 150 keV — then
        <span class="font-medium text-accent">Run</span>. The method self-calibrates from this
        patient's subcutaneous fat and muscle (no external phantom), then projects every
        soft-tissue voxel onto the water→lipid line. The whole-volume f_w / f_l maps render below,
        jet 0..1 over the CT.
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
             volume fractions ± standard deviation. -->
        {#if hover}
          <div
            class="pointer-events-none absolute z-10 rounded bg-black/75 px-2 py-1 text-[10px] leading-snug text-white tabular-nums shadow-lg"
            style="left: {hover.cx + 14}px; top: {hover.cy + 14}px;"
          >
            {#if Number.isFinite(hover.fw)}
              <div>f_w&nbsp;&nbsp;{fmt(hover.fw, 2)} ± {fmt(hover.sf, 2)}</div>
              <div>f_l&nbsp;&nbsp;{fmt(1 - hover.fw, 2)} ± {fmt(hover.sf, 2)}</div>
              <div class="text-white/55">HU {Math.round(hover.hu)} · ({hover.x}, {hover.y})</div>
            {:else}
              <div class="text-white/55">HU {Math.round(hover.hu)} · gated · ({hover.x}, {hover.y})</div>
            {/if}
          </div>
        {/if}
      </div>

      <!-- Controls + colorbar + calibration readout -->
      <div class="flex w-64 shrink-0 flex-col gap-3 overflow-y-auto border-l border-border bg-surface-secondary px-3 py-3">
        <!-- Map toggle -->
        <div>
          <div class="mb-1 text-[10px] font-semibold uppercase tracking-wider text-text-secondary/60">Map</div>
          <div class="flex gap-1">
            {#each [['fw', 'f_w'], ['fl', 'f_l'], ['sf', 'σ_f']] as opt}
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

        <!-- Anchor toggle -->
        <div>
          <div class="mb-1 text-[10px] font-semibold uppercase tracking-wider text-text-secondary/60">Lipid anchor</div>
          <div class="flex gap-1">
            {#each [['adipose', 'Adipose'], ['theoretical', 'Theoretical']] as opt}
              <button
                class="flex-1 rounded px-2 py-1 text-[11px] font-medium transition-colors {anchor === opt[0]
                  ? 'bg-accent text-white'
                  : 'bg-surface-tertiary text-text-secondary hover:text-text-primary'}"
                onclick={() => (anchor = opt[0] as WlAnchor)}
                title={opt[0] === 'adipose'
                  ? 'f_w=0 at this patient’s subcutaneous fat (fat-referenced)'
                  : 'f_l=1 ⇒ pure lipid (absolute scale, from the material LAC table)'}
              >
                {opt[1]}
              </button>
            {/each}
          </div>
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

        <!-- Calibration readout (mirrors calibration_*.toml) -->
        <div class="mt-1 border-t border-border pt-2 text-[10px] leading-relaxed text-text-secondary">
          <div class="mb-1 font-semibold uppercase tracking-wider text-text-secondary/60">Self-measured calibration</div>
          <div class="grid grid-cols-[auto_1fr] gap-x-2 tabular-nums">
            <span>HU_l fat</span><span class="text-text-primary">({fmt(calib.hu_l_adipose[0])}, {fmt(calib.hu_l_adipose[1])})</span>
            <span>HU_l theo</span><span class="text-text-primary">({fmt(calib.hu_l_theo[0])}, {fmt(calib.hu_l_theo[1])})</span>
            <span>HU muscle</span><span class="text-text-primary">({fmt(calib.hu_mus[0])}, {fmt(calib.hu_mus[1])})</span>
            <span>σ_fat</span><span class="text-text-primary">({fmt(calib.sigma_fat[0])}, {fmt(calib.sigma_fat[1])})</span>
            <span>σ {calib.low_kev}keV</span><span class="text-text-primary">{calib.ab_low[0].toFixed(4)}·HU+{fmt(calib.ab_low[1])}</span>
            <span>σ {calib.high_kev}keV</span><span class="text-text-primary">{calib.ab_high[0].toFixed(4)}·HU+{fmt(calib.ab_high[1])}</span>
            <span>ρ</span><span class="text-text-primary">{calib.rho.toFixed(3)}</span>
            <span>gate</span><span class="text-text-primary">[{fmt(calib.gate[0])}, {fmt(calib.gate[1])}]</span>
            <span>kept</span><span class="text-text-primary">{(calib.fraction_kept * 100).toFixed(1)}%</span>
            <span>fat / mus vox</span><span class="text-text-primary">{calib.n_fat.toLocaleString()} / {calib.n_mus.toLocaleString()}</span>
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

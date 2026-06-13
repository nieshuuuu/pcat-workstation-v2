<script lang="ts">
  /**
   * 3D surface plot (Plotly.js) for radial-angular material decomposition data.
   *
   * Displays the selected cross-section's material values on a (theta, r) grid.
   * The raw per-voxel GLS field is noisy, so the grid is denoised with a
   * NaN-aware Gaussian smoother (periodic in theta) and clamped to the physical
   * display range before plotting — turning the spiky lattice into a smooth
   * pericoronary "term structure" surface. The arc-length slider and the 2D
   * cross-section view share one index, so moving either moves both.
   */
  import { onMount, onDestroy } from 'svelte';
  import { MMD_MASS_MAX_MGML, type CrossSectionSurface } from '$lib/api';

  type Props = {
    surfaces: CrossSectionSurface[];
    selectedIndex: number;
    material: string;
    unit: string;
    onSliderChange: (index: number) => void;
    /** Absolute arc-length (mm) of the ostium along the centerline.
     *  Displayed arc = surface.arc_mm - arcOffsetMm. */
    arcOffsetMm?: number;
  };

  let {
    surfaces,
    selectedIndex,
    material,
    unit,
    onSliderChange,
    arcOffsetMm = 0,
  }: Props = $props();

  let plotDiv: HTMLDivElement | undefined = $state();
  let Plotly: typeof import('plotly.js-dist-min') | null = $state(null);
  let plotlyLoaded = $state(false);

  // Denoise the per-voxel GLS speckle in PHYSICAL units (derived per surface from
  // its angular/radial step), sized to PRESERVE the thin radial fat↔muscle
  // boundary that is the actual signal. Smoothing is mostly ANGULAR — averaging
  // around the ring at fixed radius cuts speckle without touching the radial
  // profile. The RADIAL kernel is deliberately tiny: enough to bridge a single
  // empty ring, NOT to smear the few-mm fat band (a wide radial blur is exactly
  // what erases the boundary we measure).
  const SMOOTH_THETA_FWHM_DEG = 35;
  const SMOOTH_R_FWHM_MM = 0.8;
  const FWHM_TO_SIGMA = 1 / 2.3548; // FWHM = 2·√(2 ln 2)·σ

  // Memoize the smoothed grid per (section, unit): revisiting a section — slider
  // scrub-back, camera rotate, resize — is then free; only a never-seen section
  // pays for the Gaussian. Reset when `surfaces` is replaced (a new MMD run).
  let smoothCache = new Map<string, (number | null)[][]>();
  let smoothCacheFor: CrossSectionSurface[] | null = null;

  // Coalesce rapid re-renders (slider drag) into one per animation frame.
  let rafId: number | null = null;

  onMount(async () => {
    const mod = await import('plotly.js-dist-min');
    Plotly = mod.default ?? mod;
    plotlyLoaded = true;
  });

  /** Material display labels. */
  function materialLabel(mat: string, u: string): string {
    if (mat === 'density') return 'Total Density (mg/mL)';
    const matName = mat.charAt(0).toUpperCase() + mat.slice(1);
    return u === 'fraction' ? `${matName} (vol %)` : `${matName} (mg/mL)`;
  }

  /** Physical display range. Fractions are clamped to [0, 100] vol% — the GLS
   *  estimator is intentionally unclamped (muscle/fibrous read f_w > 1, i.e.
   *  negative lipid, as a contamination signature), but those non-physical
   *  values are meaningless on a fat-fraction surface, so the viewer clamps
   *  them. Mass densities use a fixed [0, 1100] mg/mL scale. */
  function displayRange(u: string): [number, number] {
    return u === 'fraction' ? [0, 100] : [0, MMD_MASS_MAX_MGML];
  }

  /** 1D Gaussian kernel of the given sigma (radius = ceil(3σ)). */
  function gaussianKernel(sigma: number): number[] {
    const radius = Math.max(1, Math.ceil(3 * sigma));
    const k: number[] = [];
    for (let i = -radius; i <= radius; i++) {
      k.push(Math.exp(-(i * i) / (2 * sigma * sigma)));
    }
    return k;
  }

  /**
   * NaN-aware separable Gaussian smoothing of a [n_theta × n_radial] grid.
   * Theta (rows) wraps (periodic); r (cols) clamps at the edges. Gated cells
   * (NaN) contribute nothing and are filled from finite neighbours within
   * kernel reach (normalized convolution), so small holes close instead of
   * spiking. Returns finite values clamped to [lo, hi]; cells with no finite
   * neighbour stay NaN.
   */
  function smoothGrid(
    raw: number[],
    nTheta: number,
    nR: number,
    lo: number,
    hi: number,
    sigmaTheta: number,
    sigmaR: number,
  ): (number | null)[][] {
    // Pre-clamp finite values so a single −50 % outlier can't drag a whole
    // neighbourhood down before it is averaged.
    const v = new Float64Array(nTheta * nR);
    const ok = new Uint8Array(nTheta * nR);
    for (let i = 0; i < nTheta * nR; i++) {
      const x = raw[i];
      if (Number.isFinite(x)) {
        v[i] = Math.max(lo, Math.min(hi, x));
        ok[i] = 1;
      }
    }

    const kt = gaussianKernel(sigmaTheta);
    const kr = gaussianKernel(sigmaR);
    const rt = (kt.length - 1) / 2;
    const rr = (kr.length - 1) / 2;

    // Pass 1: smooth along theta (periodic wrap).
    const t1 = new Float64Array(nTheta * nR);
    const w1 = new Float64Array(nTheta * nR);
    for (let it = 0; it < nTheta; it++) {
      for (let ir = 0; ir < nR; ir++) {
        let acc = 0;
        let wsum = 0;
        for (let m = -rt; m <= rt; m++) {
          const jt = ((it + m) % nTheta + nTheta) % nTheta; // wrap
          const idx = jt * nR + ir;
          if (ok[idx]) {
            const w = kt[m + rt];
            acc += w * v[idx];
            wsum += w;
          }
        }
        const o = it * nR + ir;
        t1[o] = acc;
        w1[o] = wsum;
      }
    }

    // Pass 2: smooth the theta-smoothed result along r (clamped edges).
    const out: (number | null)[][] = [];
    for (let it = 0; it < nTheta; it++) {
      const row: (number | null)[] = [];
      for (let ir = 0; ir < nR; ir++) {
        let acc = 0;
        let wsum = 0;
        for (let m = -rr; m <= rr; m++) {
          let jr = ir + m;
          if (jr < 0) jr = 0;
          else if (jr >= nR) jr = nR - 1;
          const idx = it * nR + jr;
          if (w1[idx] > 0) {
            const w = kr[m + rr];
            // Re-weight by the theta-pass coverage so partially-gated columns
            // don't bias toward zero.
            acc += w * (t1[idx] / w1[idx]);
            wsum += w;
          }
        }
        row.push(wsum > 0 ? Math.max(lo, Math.min(hi, acc / wsum)) : null);
      }
      out.push(row);
    }
    return out;
  }

  /** Render / update the surface plot for the currently selected cross-section. */
  function renderPlot(P: typeof import('plotly.js-dist-min'), div: HTMLDivElement) {
    if (!surfaces || surfaces.length === 0 || selectedIndex < 0 || selectedIndex >= surfaces.length) {
      return;
    }

    const s = surfaces[selectedIndex];
    const [lo, hi] = displayRange(unit);

    // Smoothed grid, memoized per (section, unit). Cleared when surfaces change.
    if (smoothCacheFor !== surfaces) {
      smoothCache.clear();
      smoothCacheFor = surfaces;
    }
    const cacheKey = `${selectedIndex}|${unit}`;
    let z = smoothCache.get(cacheKey);
    if (!z) {
      // Convert fractions to vol% up front, then denoise + clamp on the smoother.
      const scaled = new Array<number>(s.surface.length);
      for (let i = 0; i < s.surface.length; i++) {
        const val = s.surface[i];
        scaled[i] = Number.isNaN(val) ? NaN : unit === 'fraction' ? val * 100 : val;
      }
      // Gaussian widths from THIS surface's own grid step → fixed physical FWHM,
      // independent of the backend sampling density.
      const dThetaDeg = s.n_theta > 0 ? 360 / s.n_theta : 1;
      const dRmm = s.n_radial > 1 ? Math.abs(s.r_mm[1] - s.r_mm[0]) : 1;
      const sigmaTheta = (SMOOTH_THETA_FWHM_DEG * FWHM_TO_SIGMA) / dThetaDeg;
      const sigmaR = (SMOOTH_R_FWHM_MM * FWHM_TO_SIGMA) / dRmm;
      z = smoothGrid(scaled, s.n_theta, s.n_radial, lo, hi, sigmaTheta, sigmaR);
      smoothCache.set(cacheKey, z);
    }

    const label = materialLabel(material, unit);
    const trace: Partial<Plotly.Data> = {
      type: 'surface' as const,
      x: s.r_mm,
      y: s.theta_deg,
      z,
      colorscale: 'Jet',
      cmin: lo,
      cmax: hi,
      showscale: true,
      colorbar: {
        title: { text: label, font: { size: 10, color: '#e5e5e7' } },
        tickfont: { size: 9, color: '#e5e5e7' },
        len: 0.6,
      },
      contours: {
        z: { show: true, usecolormap: true, width: 1, project: { z: false } },
      },
      hovertemplate: 'r=%{x:.1f} mm<br>theta=%{y:.0f} deg<br>value=%{z:.1f}<extra></extra>',
    };

    const layout: Partial<Plotly.Layout> = {
      paper_bgcolor: '#1c1c1e',
      plot_bgcolor: '#2c2c2e',
      font: { color: '#e5e5e7', size: 10 },
      margin: { l: 10, r: 10, t: 30, b: 10 },
      title: {
        text: `${label} — arc ${(s.arc_mm - arcOffsetMm).toFixed(1)} mm`,
        font: { size: 11, color: '#e5e5e7' },
      },
      scene: {
        xaxis: {
          title: { text: 'r (mm)', font: { size: 9 } },
          gridcolor: '#38383a',
          color: '#98989d',
        },
        yaxis: {
          title: { text: 'theta (deg)', font: { size: 9 } },
          gridcolor: '#38383a',
          color: '#98989d',
        },
        zaxis: {
          title: { text: label, font: { size: 9 } },
          gridcolor: '#38383a',
          color: '#98989d',
          range: [lo, hi],
        },
        bgcolor: '#2c2c2e',
      },
      autosize: true,
    };

    // `react` diffs against the existing plot and updates in place — far
    // cheaper than `newPlot`'s full teardown/rebuild, so dragging the arc
    // slider stays smooth.
    (P as any).react(div, [trace], layout, { responsive: true, displayModeBar: false });
  }

  // Plotly attaches a WebGL context to the plot div; if the div is removed
  // (surfaces emptied on a vessel switch / reload) without purging, the context
  // leaks and after enough cycles the webview refuses new ones and the surface
  // stops drawing. This action purges on the div's removal — covering both the
  // {#if} unmount and component teardown.
  function plotContainer(node: HTMLDivElement) {
    plotDiv = node;
    return {
      destroy() {
        if (Plotly) (Plotly as any).purge(node);
        if (plotDiv === node) plotDiv = undefined;
      },
    };
  }

  // Coalesce re-renders into one per frame (cancel any pending frame first), so
  // a fast slider drag smooths + redraws only the section the user lands on.
  function scheduleRender() {
    if (rafId != null) cancelAnimationFrame(rafId);
    rafId = requestAnimationFrame(() => {
      rafId = null;
      if (Plotly && plotDiv) renderPlot(Plotly, plotDiv);
    });
  }

  // Re-render when dependencies change.
  $effect(() => {
    if (!plotlyLoaded || !Plotly || !plotDiv) return;
    void surfaces;
    void selectedIndex;
    void material;
    void unit;
    scheduleRender();
  });

  onDestroy(() => {
    if (rafId != null) cancelAnimationFrame(rafId);
  });

  // True when the selected section's ring is entirely gated / outside the
  // volume, so the surface is all-NaN and Plotly would draw a blank — show a
  // note instead of a confusing empty plot.
  let currentEmpty = $derived(
    !!surfaces &&
      selectedIndex >= 0 &&
      selectedIndex < surfaces.length &&
      !surfaces[selectedIndex].surface.some((v) => Number.isFinite(v)),
  );
</script>

<div class="flex min-h-0 flex-1 flex-col gap-1.5 p-2">
  <!-- Plot area -->
  <div class="relative min-h-[16rem] w-full flex-1">
    {#if !plotlyLoaded}
      <div class="flex h-full items-center justify-center">
        <span class="text-xs text-text-secondary">Loading chart library...</span>
      </div>
    {:else if !surfaces || surfaces.length === 0}
      <div class="flex h-full items-center justify-center">
        <span class="text-xs text-text-secondary">Run MMD to generate surface data</span>
      </div>
    {:else}
      <div use:plotContainer class="h-full w-full"></div>
      {#if currentEmpty}
        <div class="pointer-events-none absolute inset-0 flex items-center justify-center">
          <span class="rounded bg-black/40 px-2 py-1 text-xs text-text-secondary">
            No decomposed tissue in this section's ring
          </span>
        </div>
      {/if}
    {/if}
  </div>

  <!-- Arc-length slider (shared with the 2D cross-section view) -->
  {#if surfaces && surfaces.length > 1}
    <div class="flex shrink-0 items-center gap-2 px-2">
      <span class="shrink-0 text-[10px] text-text-secondary">Arc</span>
      <input
        type="range"
        min="0"
        max={surfaces.length - 1}
        value={selectedIndex}
        oninput={(e) => onSliderChange(parseInt((e.target as HTMLInputElement).value))}
        class="h-1 flex-1 cursor-pointer appearance-none rounded-full bg-surface-tertiary accent-accent"
      />
      <span class="shrink-0 text-[10px] tabular-nums text-text-secondary">
        {surfaces[selectedIndex] !== undefined
          ? (surfaces[selectedIndex].arc_mm - arcOffsetMm).toFixed(1)
          : '—'} mm
      </span>
    </div>
  {/if}
</div>

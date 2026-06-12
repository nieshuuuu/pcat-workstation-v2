<script lang="ts">
  /**
   * 3D surface plot (Plotly.js) for radial-angular material decomposition data.
   *
   * Displays a surface plot of material values on a (theta, r) grid for the
   * selected cross-section. An arc-length slider selects which cross-section
   * to display.
   */
  import { onMount, tick } from 'svelte';
  import type { CrossSectionSurface } from '$lib/api';

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

  /** Render the surface plot for the currently selected cross-section. */
  function renderPlot(P: typeof import('plotly.js-dist-min'), div: HTMLDivElement) {
    if (!surfaces || surfaces.length === 0 || selectedIndex < 0 || selectedIndex >= surfaces.length) {
      return;
    }

    const s = surfaces[selectedIndex];

    // Display range per unit. Fractions are vol% in [0, 100]; the noise-aware
    // GLS estimate can read slightly outside (low-HU fat > 100 %, fibrous/wall
    // < 0 %), so clamp for a readable surface. Gated (non-soft-tissue) voxels
    // arrive as NaN and render as gaps. Raw values are still exported via CSV.
    const [zMin, zMax] = unit === 'fraction' ? [0, 100] : [0, 1000];

    // Build z matrix: [n_theta x n_radial], NaN -> null (gap), else clamped.
    const z: (number | null)[][] = [];
    for (let it = 0; it < s.n_theta; it++) {
      const row: (number | null)[] = [];
      for (let ir = 0; ir < s.n_radial; ir++) {
        const val = s.surface[it * s.n_radial + ir];
        if (isNaN(val)) {
          row.push(null);
          continue;
        }
        const scaled = unit === 'fraction' ? val * 100 : val;
        row.push(Math.max(zMin, Math.min(zMax, scaled)));
      }
      z.push(row);
    }

    const trace: Partial<Plotly.Data> = {
      type: 'surface' as const,
      x: s.r_mm,
      y: s.theta_deg,
      z: z,
      colorscale: 'Viridis',
      cmin: zMin,
      cmax: zMax,
      showscale: true,
      colorbar: {
        title: { text: materialLabel(material, unit), font: { size: 10, color: '#e5e5e7' } },
        tickfont: { size: 9, color: '#e5e5e7' },
        len: 0.6,
      },
      hovertemplate:
        'r=%{x:.1f} mm<br>theta=%{y:.0f} deg<br>value=%{z:.3f}<extra></extra>',
    };

    const layout: Partial<Plotly.Layout> = {
      paper_bgcolor: '#1c1c1e',
      plot_bgcolor: '#2c2c2e',
      font: { color: '#e5e5e7', size: 10 },
      margin: { l: 10, r: 10, t: 30, b: 10 },
      title: {
        text: `${materialLabel(material, unit)} — arc ${(s.arc_mm - arcOffsetMm).toFixed(1)} mm`,
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
          title: { text: materialLabel(material, unit), font: { size: 9 } },
          range: [zMin, zMax],
          gridcolor: '#38383a',
          color: '#98989d',
        },
        bgcolor: '#2c2c2e',
      },
      autosize: true,
    };

    (P as any).newPlot(div, [trace], layout, { responsive: true, displayModeBar: false });
  }

  // Re-render when dependencies change.
  $effect(() => {
    if (!plotlyLoaded || !Plotly || !plotDiv) return;
    void surfaces;
    void selectedIndex;
    void material;
    void unit;
    tick().then(() => {
      if (Plotly && plotDiv) renderPlot(Plotly, plotDiv);
    });
  });
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
      <div bind:this={plotDiv} class="h-full w-full"></div>
    {/if}
  </div>

  <!-- Arc-length slider -->
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

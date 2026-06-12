<script lang="ts">
  /**
   * FAI Analysis Dashboard — tabbed panel showing pipeline results.
   *
   * Tabs:
   *   1. Overview — per-vessel FAI summary cards
   *   2. Histograms — Plotly.js HU distribution charts
   *   3. Radial Profile — mean HU vs distance from vessel wall
   *   4. Angular — SVG ring cross-section with per-sector HU values
   */
  import { onMount, tick } from 'svelte';
  import { pipelineStore, type FaiStats } from '$lib/stores/pipelineStore.svelte';
  import { VESSEL_COLORS, type Vessel } from '$lib/stores/seedStore.svelte';

  const TABS = ['Overview', 'Histograms', 'Radial Profile', 'Angular'] as const;
  type Tab = (typeof TABS)[number];

  let activeTab = $state<Tab>('Overview');
  let chartDiv: HTMLDivElement | undefined = $state();
  let radialChartEl: HTMLDivElement | undefined = $state();
  let plotlyLoaded = $state(false);
  let Plotly: typeof import('plotly.js-dist-min') | null = $state(null);

  const vesselOrder: Vessel[] = ['RCA', 'LAD', 'LCx'];

  let vesselResults = $derived(
    pipelineStore.results
      ? vesselOrder
          .filter((v) => v in pipelineStore.results!)
          .map((v) => ({ vessel: v, stats: pipelineStore.results![v] }))
      : [],
  );

  onMount(async () => {
    const mod = await import('plotly.js-dist-min');
    Plotly = mod.default ?? mod;
    plotlyLoaded = true;
  });

  // Histogram effect
  $effect(() => {
    if (activeTab !== 'Histograms' || !plotlyLoaded || !Plotly) return;
    if (!pipelineStore.results || vesselResults.length === 0) return;
    tick().then(() => {
      if (!chartDiv || !Plotly) return;
      renderHistogram(Plotly, chartDiv, vesselResults);
    });
  });

  // Radial profile effect
  $effect(() => {
    if (activeTab !== 'Radial Profile' || !plotlyLoaded || !Plotly) return;
    if (!pipelineStore.results || vesselResults.length === 0) return;
    tick().then(() => {
      if (!radialChartEl || !Plotly) return;
      renderRadialProfile(Plotly, radialChartEl, vesselResults);
    });
  });

  function renderHistogram(
    P: typeof import('plotly.js-dist-min'),
    div: HTMLDivElement,
    data: { vessel: Vessel; stats: FaiStats }[],
  ) {
    const traces: Partial<Plotly.Data>[] = data.map(({ vessel, stats }) => ({
      x: stats.histogram_bins,
      y: stats.histogram_counts,
      type: 'bar' as const,
      name: vessel,
      marker: { color: VESSEL_COLORS[vessel], opacity: 0.7 },
    }));

    const shapes: Partial<Plotly.Shape>[] = [
      { type: 'line', x0: -190, x1: -190, y0: 0, y1: 1, yref: 'paper', line: { color: '#98989d', width: 1, dash: 'dash' } },
      { type: 'line', x0: -30, x1: -30, y0: 0, y1: 1, yref: 'paper', line: { color: '#98989d', width: 1, dash: 'dash' } },
      { type: 'line', x0: -70.1, x1: -70.1, y0: 0, y1: 1, yref: 'paper', line: { color: '#ff453a', width: 2, dash: 'dot' } },
    ];

    const layout: Partial<Plotly.Layout> = {
      paper_bgcolor: '#1c1c1e', plot_bgcolor: '#2c2c2e',
      font: { color: '#e5e5e7', size: 10 },
      margin: { l: 45, r: 15, t: 30, b: 40 },
      title: { text: 'HU Distribution (FAI Window)', font: { size: 12, color: '#e5e5e7' } },
      xaxis: { title: { text: 'HU', font: { size: 10 } }, gridcolor: '#38383a', zerolinecolor: '#38383a' },
      yaxis: { title: { text: 'Count', font: { size: 10 } }, gridcolor: '#38383a', zerolinecolor: '#38383a' },
      barmode: 'overlay',
      shapes: shapes as Plotly.Layout['shapes'],
      legend: { font: { color: '#e5e5e7', size: 10 }, bgcolor: 'rgba(0,0,0,0)' },
      autosize: true,
    };

    (P as any).newPlot(div, traces, layout, { responsive: true, displayModeBar: false });
  }

  function renderRadialProfile(
    P: typeof import('plotly.js-dist-min'),
    div: HTMLDivElement,
    data: { vessel: Vessel; stats: FaiStats }[],
  ) {
    const traces: Partial<Plotly.Data>[] = [];

    for (const { vessel, stats } of data) {
      const profile = stats.radial_profile;
      if (!profile) continue;

      const distances = profile.distances_mm ?? [];
      const meanHu = profile.mean_hu ?? [];
      const stdHu = profile.std_hu ?? [];
      const color = VESSEL_COLORS[vessel];

      // Filter out empty rings before plotting. Rust emits f64::NAN for an empty
      // ring, which serde serializes to JSON `null`; `isFinite(null)` is `true`
      // (null coerces to 0), so the `!= null` guard is required — without it a
      // null mean slips through and the band math `null + std` yields 0, spiking
      // the confidence band to the top of the chart. (Matches the null guards on
      // the sector table and the HU color scale below.)
      const valid = distances.map((d, i) => ({ d, m: meanHu[i], s: stdHu[i] }))
        .filter(v => v.m != null && isFinite(v.m));
      if (valid.length === 0) continue;

      const vd = valid.map(v => v.d);
      const vm = valid.map(v => v.m);
      const vs = valid.map(v => v.s);

      // Upper band
      traces.push({
        x: vd, y: vm.map((m, i) => m + (vs[i] ?? 0)),
        type: 'scatter' as const, mode: 'lines' as const,
        line: { width: 0 }, showlegend: false, hoverinfo: 'skip' as const,
      });
      // Lower band filled to upper
      traces.push({
        x: vd, y: vm.map((m, i) => m - (vs[i] ?? 0)),
        type: 'scatter' as const, mode: 'lines' as const,
        line: { width: 0 }, fill: 'tonexty' as const,
        fillcolor: color + '20', showlegend: false, hoverinfo: 'skip' as const,
      });
      // Mean line
      traces.push({
        x: vd, y: vm,
        type: 'scatter' as const, mode: 'lines+markers' as const,
        name: vessel, line: { color, width: 2 },
        marker: { size: 3, color },
      });
    }

    const shapes: Partial<Plotly.Shape>[] = [
      { type: 'rect', x0: 0, x1: 20, y0: -190, y1: -30, fillcolor: 'rgba(99,102,241,0.08)', line: { width: 0 } },
      { type: 'line', x0: 0, x1: 20, y0: -190, y1: -190, line: { color: '#98989d', width: 1, dash: 'dash' } },
      { type: 'line', x0: 0, x1: 20, y0: -30, y1: -30, line: { color: '#98989d', width: 1, dash: 'dash' } },
      { type: 'line', x0: 0, x1: 20, y0: -70.1, y1: -70.1, line: { color: '#ff453a', width: 2, dash: 'dot' } },
    ];

    const layout: Partial<Plotly.Layout> = {
      paper_bgcolor: '#1c1c1e', plot_bgcolor: '#2c2c2e',
      font: { color: '#e5e5e7', size: 10 },
      margin: { l: 50, r: 15, t: 30, b: 45 },
      title: { text: 'Radial Profile — Mean HU vs Distance from Wall', font: { size: 12, color: '#e5e5e7' } },
      xaxis: { title: { text: 'Distance from vessel wall (mm)', font: { size: 10 } }, range: [0, 20], gridcolor: '#38383a', zerolinecolor: '#38383a' },
      yaxis: { title: { text: 'Mean HU (fat range)', font: { size: 10 } }, range: [-100, -50], gridcolor: '#38383a', zerolinecolor: '#38383a' },
      shapes: shapes as Plotly.Layout['shapes'],
      legend: { font: { color: '#e5e5e7', size: 10 }, bgcolor: 'rgba(0,0,0,0)' },
      autosize: true,
    };

    (P as any).newPlot(div, traces, layout, { responsive: true, displayModeBar: false });
  }

  /** Map HU value to color: green (low/healthy) → yellow → red (high/inflamed). */
  function huToColor(hu: number | null): string {
    if (hu == null || !isFinite(hu)) return '#555';
    // Map -190..-30 to 0..1
    const t = Math.max(0, Math.min(1, (hu - (-190)) / ((-30) - (-190))));
    // Green → Yellow → Red
    const r = Math.round(t < 0.5 ? t * 2 * 255 : 255);
    const g = Math.round(t < 0.5 ? 255 : (1 - (t - 0.5) * 2) * 255);
    return `rgb(${r}, ${g}, 40)`;
  }
</script>

<div class="flex h-full w-full flex-col bg-surface-secondary">
  <!-- Tab bar -->
  <div class="flex shrink-0 border-b border-border">
    {#each TABS as tab}
      <button
        class="px-3 py-2 text-[11px] font-medium transition-colors"
        class:text-accent={activeTab === tab}
        class:border-b-2={activeTab === tab}
        class:border-accent={activeTab === tab}
        class:text-text-secondary={activeTab !== tab}
        class:hover:text-text-primary={activeTab !== tab}
        onclick={() => (activeTab = tab)}
      >
        {tab}
      </button>
    {/each}
  </div>

  <!-- Tab content -->
  <div class="min-h-0 flex-1 overflow-y-auto">
    {#if activeTab === 'Overview'}
      <div class="flex flex-col gap-2.5 p-3">
        {#if vesselResults.length === 0}
          <div class="flex flex-1 items-center justify-center py-8">
            <p class="text-xs text-text-secondary">No results available</p>
          </div>
        {:else}
          {#each vesselResults as { vessel, stats }}
            {@const isHigh = stats.fai_risk === 'HIGH'}
            <div class="rounded-lg border border-border bg-surface p-3">
              <div class="mb-2 flex items-center justify-between">
                <div class="flex items-center gap-2">
                  <span class="h-2.5 w-2.5 rounded-full" style="background-color: {VESSEL_COLORS[vessel]}"></span>
                  <span class="text-xs font-semibold text-text-primary">{vessel}</span>
                </div>
                <span
                  class="rounded px-1.5 py-0.5 text-[10px] font-bold"
                  style="background-color: {isHigh ? 'rgba(255,69,58,0.2)' : 'rgba(48,209,88,0.2)'}; color: {isHigh ? 'var(--color-error)' : 'var(--color-success)'}"
                >
                  {stats.fai_risk}
                </span>
              </div>
              <div class="mb-2">
                <span class="text-[10px] text-text-secondary">Mean HU</span>
                <p class="text-xl font-bold tabular-nums" style="color: {isHigh ? 'var(--color-error)' : 'var(--color-success)'}">
                  {stats.hu_mean.toFixed(1)}
                </p>
              </div>
              <div class="grid grid-cols-3 gap-2">
                <div><span class="text-[9px] text-text-secondary">Median HU</span><p class="text-[11px] tabular-nums text-text-primary">{stats.hu_median.toFixed(1)}</p></div>
                <div><span class="text-[9px] text-text-secondary">Std HU</span><p class="text-[11px] tabular-nums text-text-primary">{stats.hu_std.toFixed(1)}</p></div>
                <div><span class="text-[9px] text-text-secondary">Fat %</span><p class="text-[11px] tabular-nums text-text-primary">{(stats.fat_fraction * 100).toFixed(1)}%</p></div>
                <div><span class="text-[9px] text-text-secondary">VOI Voxels</span><p class="text-[11px] tabular-nums text-text-primary">{stats.n_voi_voxels.toLocaleString()}</p></div>
                <div><span class="text-[9px] text-text-secondary">Fat Voxels</span><p class="text-[11px] tabular-nums text-text-primary">{stats.n_fat_voxels.toLocaleString()}</p></div>
              </div>
              <div class="mt-2 border-t border-border/50 pt-1.5">
                <span class="text-[8px] text-text-secondary/50">CRISP-CT VOI: 1mm gap + 3mm ring | FAI: -190 to -30 HU | Risk: -70.1 HU | Segment: 0-40mm</span>
              </div>
            </div>
          {/each}
        {/if}
      </div>

    {:else if activeTab === 'Histograms'}
      <div class="flex h-full w-full items-center justify-center p-2">
        {#if !plotlyLoaded}
          <p class="text-xs text-text-secondary">Loading chart library...</p>
        {:else if vesselResults.length === 0}
          <p class="text-xs text-text-secondary">No histogram data available</p>
        {:else}
          <div bind:this={chartDiv} class="h-full w-full" style="min-height: 200px;"></div>
        {/if}
      </div>

    {:else if activeTab === 'Radial Profile'}
      <div class="h-full w-full p-2">
        {#if !plotlyLoaded}
          <div class="flex h-full items-center justify-center">
            <p class="text-xs text-text-secondary">Loading chart library...</p>
          </div>
        {:else if vesselResults.length === 0 || !vesselResults.some(({ stats }) => stats.radial_profile)}
          <div class="flex h-full items-center justify-center">
            <p class="text-xs text-text-secondary">No radial profile data available</p>
          </div>
        {:else}
          <div bind:this={radialChartEl} class="h-full w-full" style="min-height: 200px;"></div>
        {/if}
      </div>

    {:else if activeTab === 'Angular'}
      <!-- Angular: SVG ring cross-section with per-sector HU values -->
      <div class="flex flex-col items-center gap-6 p-4 overflow-y-auto">
        {#if vesselResults.length === 0 || !vesselResults.some(({ stats }) => stats.angular_asymmetry)}
          <div class="flex h-full w-full items-center justify-center">
            <p class="text-xs text-text-secondary">No angular data available</p>
          </div>
        {:else}
          {#each vesselResults as { vessel, stats }}
            {#if stats.angular_asymmetry}
              {@const sectors = stats.angular_asymmetry.sectors}
              {@const n = sectors.length}
              {@const cx = 200}
              {@const cy = 200}
              {@const rInner = 55}
              {@const rOuter = 140}
              {@const rLabel = 170}
              <div class="flex flex-col items-center gap-3">
                <div class="flex items-center gap-2">
                  <span class="h-3 w-3 rounded-full" style="background-color: {VESSEL_COLORS[vessel]}"></span>
                  <span class="text-sm font-semibold text-text-primary">{vessel}</span>
                </div>
                <svg viewBox="0 0 400 400" class="h-[340px] w-[340px]">
                  <!-- Lumen circle -->
                  <circle {cx} {cy} r={rInner} fill="#2c2c2e" stroke="#555" stroke-width="1.5" />
                  <text x={cx} y={cy} text-anchor="middle" dominant-baseline="central" fill="#98989d" font-size="13" font-weight="500">Lumen</text>

                  <!-- Sectors -->
                  {#each sectors as sector, i}
                    {@const a0 = (i / n) * 2 * Math.PI - Math.PI / 2}
                    {@const a1 = ((i + 1) / n) * 2 * Math.PI - Math.PI / 2}
                    {@const aMid = (a0 + a1) / 2}
                    {@const x1i = cx + rInner * Math.cos(a0)}
                    {@const y1i = cy + rInner * Math.sin(a0)}
                    {@const x1o = cx + rOuter * Math.cos(a0)}
                    {@const y1o = cy + rOuter * Math.sin(a0)}
                    {@const x2i = cx + rInner * Math.cos(a1)}
                    {@const y2i = cy + rInner * Math.sin(a1)}
                    {@const x2o = cx + rOuter * Math.cos(a1)}
                    {@const y2o = cy + rOuter * Math.sin(a1)}
                    {@const largeArc = (a1 - a0) > Math.PI ? 1 : 0}
                    {@const color = huToColor(sector.hu_mean)}
                    {@const txr = (rInner + rOuter) / 2}
                    {@const tx = cx + txr * Math.cos(aMid)}
                    {@const ty = cy + txr * Math.sin(aMid)}
                    {@const lx = cx + rLabel * Math.cos(aMid)}
                    {@const ly = cy + rLabel * Math.sin(aMid)}

                    <!-- Sector arc path -->
                    <path
                      d="M {x1i} {y1i} L {x1o} {y1o} A {rOuter} {rOuter} 0 {largeArc} 1 {x2o} {y2o} L {x2i} {y2i} A {rInner} {rInner} 0 {largeArc} 0 {x1i} {y1i}"
                      fill={color}
                      stroke="#1c1c1e"
                      stroke-width="1"
                      opacity="0.85"
                    />
                    <!-- HU value in sector -->
                    <text x={tx} y={ty} text-anchor="middle" dominant-baseline="central" fill="#000" font-size="13" font-weight="bold">
                      {sector.hu_mean != null && isFinite(sector.hu_mean) ? sector.hu_mean.toFixed(0) : '—'}
                    </text>
                    <!-- Sector label outside -->
                    <text x={lx} y={ly} text-anchor="middle" dominant-baseline="central" fill="#c0c0c4" font-size="11" font-weight="500">
                      {sector.label}
                    </text>
                  {/each}
                </svg>
                <!-- Color legend -->
                <div class="flex items-center gap-1 text-[9px] text-text-secondary">
                  <span class="inline-block h-2 w-6 rounded" style="background: linear-gradient(to right, rgb(0,255,40), rgb(255,255,40), rgb(255,0,40))"></span>
                  <span>-190</span>
                  <span class="mx-1">to</span>
                  <span>-30 HU</span>
                </div>
              </div>
            {/if}
          {/each}
        {/if}
      </div>

    {/if}
  </div>
</div>

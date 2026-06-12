<script lang="ts">
  /**
   * MMD Analysis View — full-screen analysis tab for multi-material
   * decomposition results.
   *
   * Layout:
   *   ┌────────────────────┬──────────────────────────┐
   *   │                    │  Surface Plot (Plotly 3D) │
   *   │  Cross-Section     │  [arc slider]             │
   *   │  (SnakeEditor)     ├──────────────────────────┤
   *   │                    │  (radial profile TBD)    │
   *   ├────────────────────┴──────────────────────────┤
   *   │ [Water|Lipid|Total ρ] [Vol%|Mass]             │
   *   ├───────────────────────────────────────────────┤
   *   │ CrossSectionStrip: [0mm] [2mm] [4mm*] ...     │
   *   ├───────────────────────────────────────────────┤
   *   │ Contour: [Init][Evolve][Reset][Accept] | MMD  │
   *   └───────────────────────────────────────────────┘
   */
  import {
    type AnnotationTarget,
    type CrossSectionSurface,
    type MmdSummary,
    generateAnnotationTargets,
    buildCprFrame,
    sampleSurfaces,
    runMmdOnRoi,
    saveAnnotations,
    exportMmdCsv,
    useVesselWallAsContour,
    getMmdOverlay,
  } from '$lib/api';
  import { volumeStore } from '$lib/stores/volumeStore.svelte';
  import { seedStore } from '$lib/stores/seedStore.svelte';

  import SnakeEditor from './SnakeEditor.svelte';
  import CrossSectionStrip from './CrossSectionStrip.svelte';
  import OverlaySelector from './OverlaySelector.svelte';
  import SurfacePlotPanel from './SurfacePlotPanel.svelte';

  type Props = {
    /** World-space centerline points for the vessel segment. */
    centerlineMm: [number, number, number][];
  };

  let { centerlineMm }: Props = $props();

  /* ── Derived from volume store ──────────────────────── */

  let dicomPath = $derived(volumeStore.current?.dicomPath ?? '');
  let patientName = $derived(volumeStore.current?.patientName ?? 'unknown');

  /* ── State ───────────────────────────────────────────── */

  let targets = $state<AnnotationTarget[]>([]);
  let selectedIndex = $state(0);
  let statusMap = $state<Record<number, 'pending' | 'in-progress' | 'done'>>({});
  let snakePoints = $state<Record<number, [number, number][]>>({});

  let material = $state('lipid');
  let unit = $state('fraction');

  let surfaces = $state<CrossSectionSurface[]>([]);

  let mmdSummary = $state<MmdSummary | null>(null);
  let mmdBusy = $state(false);
  let mmdError = $state('');

  /** Per-target flat material overlay (pixels×pixels) for the current
   *  material/unit. Lazily fetched after MMD runs when the user focuses a
   *  section, then cached until material/unit/mmdSummary change. */
  let overlayCache = $state<Record<number, number[]>>({});
  let currentOverlay = $derived<number[] | null>(
    material === 'ct' ? null : (overlayCache[selectedIndex] ?? null),
  );

  let loadingTargets = $state(false);
  let saveBusy = $state(false);
  let saveMsg = $state('');
  let exportBusy = $state(false);

  /* ── Derived ─────────────────────────────────────────── */

  let currentTarget = $derived(targets[selectedIndex] ?? null);
  let currentSnake = $derived(snakePoints[selectedIndex] ?? null);
  let contourCount = $derived(Object.keys(snakePoints).length);
  let totalCount = $derived(targets.length);

  // Arc-length offset: the first target sits at `start_arc_mm` (the ostium's
  // arc-length along the centerline). Subtracting this in the UI displays
  // distances as offsets from the ostium rather than from centerline point 0.
  let arcOffsetMm = $derived(targets[0]?.arc_mm ?? 0);

  /* ── Load targets on mount / centerline change ──────
   *
   * Fingerprint of the centerline we last loaded targets for. The parent
   * re-derives `centerlineMm` on every seedStore tick (returns a fresh
   * array even when the vessel path is identical), and Svelte tracks
   * dependencies by identity — without this guard every unrelated seed
   * tweak would fire `loadTargets()` and wipe an existing MMD result.
   * We skip when the centerline shape is unchanged: same point count +
   * same first and last waypoints. Only genuine edits trigger a reload. */
  let lastCenterlineKey = $state<string | null>(null);

  function centerlineKey(cl: [number, number, number][]): string {
    if (cl.length === 0) return '';
    const first = cl[0];
    const last = cl[cl.length - 1];
    return `${cl.length}|${first[0]},${first[1]},${first[2]}|${last[0]},${last[1]},${last[2]}`;
  }

  $effect(() => {
    if (centerlineMm.length < 2) return;
    // Include the ostium fraction in the key so SETTING/MOVING the ostium after
    // the centerline exists regenerates the cross-sections to start there
    // (otherwise the effect only fired on centerline geometry changes and the
    // sections kept starting at waypoint 0).
    const ostiumFrac = seedStore.getOstiumFraction(seedStore.activeVessel);
    const key = `${centerlineKey(centerlineMm)}|o:${ostiumFrac ?? 'none'}`;
    if (key === lastCenterlineKey) return;
    lastCenterlineKey = key;
    loadTargets();
  });

  /** Downsample a polyline to at most `maxPts` evenly-spaced points (matches
   *  CprView's frame build so both views produce an identical CPR frame). */
  function downsample(
    pts: [number, number, number][],
    maxPts: number,
  ): [number, number, number][] {
    if (pts.length <= maxPts) return pts;
    const step = (pts.length - 1) / (maxPts - 1);
    const result: [number, number, number][] = [];
    for (let i = 0; i < maxPts - 1; i++) result.push(pts[Math.round(i * step)]);
    result.push(pts[pts.length - 1]);
    return result;
  }

  async function loadTargets() {
    loadingTargets = true;
    try {
      // Rebuild the CPR frame from the CURRENT centerline first, so the
      // cross-sections AND the ostium offset are resolved against it rather than
      // a stale frame left over from a previous centerline (the cause of MMD
      // starting at point 0). Mirror CprView exactly: downsample to 100 points,
      // convert [x,y,z] (cornerstone) → [z,y,x] (patient), pixelsWide 768.
      const centerlineZyx = downsample(centerlineMm, 100).map(
        ([x, y, z]) => [z, y, x] as [number, number, number],
      );
      try {
        await buildCprFrame(centerlineZyx, 768);
      } catch (err) {
        console.warn('Failed to rebuild CPR frame for MMD:', err);
      }

      // Prefer the user-placed ostium marker over the first centerline waypoint
      // so MMD cross-sections start at the true coronary ostium (matches FAI).
      // seedStore returns [x, y, z] (cornerstone world); Rust pipeline expects
      // [z, y, x] — flip to match (same convention as pipelineStore.toZyx).
      const ostiumWorld = seedStore.getOstiumWorldPosForVessel(seedStore.activeVessel);
      const ostiumZyx: [number, number, number] | null = ostiumWorld
        ? [ostiumWorld[2], ostiumWorld[1], ostiumWorld[0]]
        : null;
      targets = await generateAnnotationTargets(centerlineMm, ostiumZyx);
      selectedIndex = 0;
      snakePoints = {};
      statusMap = {};
      surfaces = [];
      mmdSummary = null;
      mmdError = '';
      overlayCache = {};

      // Always auto-adopt the vessel wall — clean 16-point contours every
      // open. Prior saved annotations are reachable via the Save/Load flow
      // but never auto-applied; this avoids stale dense polygons from older
      // sessions leaking into a fresh open.
      try {
        const adopted = await useVesselWallAsContour({ all: true });
        const initSnake: Record<number, [number, number][]> = {};
        const initStatus: Record<number, 'pending' | 'in-progress' | 'done'> = {};
        for (const c of adopted) {
          initSnake[c.target_index] = c.points;
          initStatus[c.target_index] = 'done';
        }
        snakePoints = initSnake;
        statusMap = initStatus;

        // Auto-run so the cross-section heatmap + surface appear immediately,
        // reusing the whole-volume Water/Lipid calibration if it exists (no
        // manual click, no expensive re-calibration). The cross-section overlay
        // is computed from the same calibration, so the two views align.
        if (Object.keys(initSnake).length > 0) {
          handleRunMmd();
        }
      } catch (err) {
        console.warn('Auto-adopt vessel wall failed:', err);
      }
    } catch (err) {
      console.error('Failed to generate annotation targets:', err);
    } finally {
      loadingTargets = false;
    }
  }

  /* ── Handlers ─────────────────────────────────────────── */

  // Single shared cross-section index drives BOTH the 2D big window (SnakeEditor)
  // and the 3D surface plot, so moving either one moves the other. `targets`
  // and `surfaces` are 1:1 (every target auto-adopts a lumen contour → one
  // surface), so the same index addresses both.
  function handleSelect(index: number) {
    const max = Math.max(0, targets.length - 1);
    selectedIndex = Math.max(0, Math.min(max, index));
  }

  // `selectedIndex` is a TARGET index (drives the 2D editor + the overlay). The
  // 3D plot indexes its own `surfaces[]`, which the backend PACKS — a contourless
  // target yields no surface. Map between the two by arc-length so both views
  // always show the SAME physical section even when the counts differ; in the
  // common 1:1 case this is just the identity.
  function targetToSurfaceIndex(ti: number): number {
    if (surfaces.length === 0) return 0;
    const arc = targets[ti]?.arc_mm;
    if (arc == null) return Math.min(ti, surfaces.length - 1);
    let best = 0;
    let bestD = Infinity;
    for (let i = 0; i < surfaces.length; i++) {
      const d = Math.abs(surfaces[i].arc_mm - arc);
      if (d < bestD) { bestD = d; best = i; }
    }
    return best;
  }
  function surfaceToTargetIndex(si: number): number {
    const arc = surfaces[si]?.arc_mm;
    if (arc == null) return Math.min(si, Math.max(0, targets.length - 1));
    let best = 0;
    let bestD = Infinity;
    for (let i = 0; i < targets.length; i++) {
      const d = Math.abs(targets[i].arc_mm - arc);
      if (d < bestD) { bestD = d; best = i; }
    }
    return best;
  }
  let selectedSurfaceIndex = $derived(targetToSurfaceIndex(selectedIndex));

  async function handleSave() {
    if (saveBusy || !dicomPath) return;
    saveBusy = true;
    saveMsg = '';
    try {
      await saveAnnotations(dicomPath, centerlineMm);
      saveMsg = 'Saved';
      setTimeout(() => { saveMsg = ''; }, 2000);
    } catch (err) {
      saveMsg = 'Save failed';
      console.error('Save failed:', err);
    } finally {
      saveBusy = false;
    }
  }

  async function handleExportCsv() {
    if (exportBusy || !mmdSummary) return;
    exportBusy = true;
    try {
      const csv = await exportMmdCsv(patientName);
      const blob = new Blob([csv], { type: 'text/csv' });
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url;
      a.download = `${patientName}_mmd_surfaces.csv`;
      a.click();
      URL.revokeObjectURL(url);
    } catch (err) {
      console.error('CSV export failed:', err);
    } finally {
      exportBusy = false;
    }
  }

  function handleMaterialChange(m: string) {
    material = m;
    overlayCache = {}; // invalidate — overlay is material-specific
    if (mmdSummary) {
      refreshSurfaces();
    }
  }

  function handleUnitChange(u: string) {
    unit = u;
    overlayCache = {};
    if (mmdSummary) {
      refreshSurfaces();
    }
  }

  // Lazily fetch the material overlay for the currently-selected cross-section
  // whenever MMD has produced a result and the cache hasn't seen this target
  // under the active material/unit. Skip when the user has set material='ct'
  // (no overlay, plain HU grayscale).
  //
  // Single-in-flight coalescing: scrubbing the cross-section index fires this
  // effect on every step, but only one overlay request runs at a time. When it
  // settles we fetch the section the user actually landed on — so a fast scrub
  // issues a couple of requests, not one IPC call per intermediate section.
  let overlayInflight = false;
  function ensureOverlay() {
    if (!mmdSummary || material === 'ct') return;
    if (overlayInflight) return;
    const idx = selectedIndex;
    if (overlayCache[idx]) return;
    // Snapshot the material/unit this request is for. `overlayCache` is keyed by
    // section index ONLY and is cleared on every material/unit switch, so a
    // result that arrives after the user changed material must be discarded —
    // otherwise it would write water data into a lipid-keyed cache and the 2D
    // view would show the wrong decomposition.
    const reqMat = material;
    const reqUnit = unit;
    overlayInflight = true;
    getMmdOverlay(idx, reqMat, reqUnit)
      .then((data) => {
        if (reqMat === material && reqUnit === unit) {
          overlayCache = { ...overlayCache, [idx]: data };
        }
      })
      .catch((err) => {
        console.warn('Overlay fetch failed:', err);
      })
      .finally(() => {
        overlayInflight = false;
        // Re-fire only if the user moved (section/material/unit) during the
        // request, so we chase the latest selection without spinning on a
        // section whose fetch keeps failing.
        const moved = selectedIndex !== idx || material !== reqMat || unit !== reqUnit;
        if (moved && material !== 'ct' && !overlayCache[selectedIndex]) ensureOverlay();
      });
  }
  $effect(() => {
    void selectedIndex;
    void material;
    void unit;
    void mmdSummary;
    ensureOverlay();
  });

  async function handleRunMmd() {
    if (mmdBusy) return;
    mmdBusy = true;
    mmdError = '';
    overlayCache = {}; // stale now
    try {
      mmdSummary = await runMmdOnRoi('gls');
      await refreshSurfaces();
    } catch (err) {
      mmdError = err instanceof Error ? err.message : String(err);
      console.error('MMD failed:', err);
    } finally {
      mmdBusy = false;
    }
  }

  async function refreshSurfaces() {
    if (material === 'ct') return; // surfaces are only defined for decomposed materials
    try {
      surfaces = await sampleSurfaces(material, unit);
      // The 2D (targets[]) and 3D (surfaces[]) views are linked by arc-length
      // (see targetToSurfaceIndex), so a packed surfaces vector no longer desyncs
      // them — but a mismatch still means some sections lack a 3D surface, which
      // is worth surfacing.
      if (surfaces.length !== targets.length) {
        console.warn(
          `MMD surface/target count mismatch (${surfaces.length} vs ${targets.length}); ` +
            'views stay aligned by arc-length, but some sections have no 3D surface.',
        );
      }
      // selectedIndex is a TARGET index — clamp to the targets range (NOT
      // surfaces.length; the 3D plot is addressed via the arc-length mapping).
      selectedIndex = Math.min(selectedIndex, Math.max(0, targets.length - 1));
    } catch (err) {
      console.error('Failed to sample surfaces:', err);
    }
  }
</script>

<div class="flex h-full w-full flex-col overflow-hidden bg-surface">
  {#if loadingTargets}
    <div class="flex flex-1 items-center justify-center">
      <span class="text-xs text-text-secondary">Generating cross-sections...</span>
    </div>
  {:else if targets.length === 0}
    <div class="flex flex-1 items-center justify-center">
      <span class="text-xs text-text-secondary">No annotation targets. Ensure a centerline with 2+ points is selected.</span>
    </div>
  {:else}
    {#if !mmdSummary && !mmdBusy}
      <div class="shrink-0 border-b border-border bg-accent/5 px-3 py-1 text-[11px] text-text-secondary">
        Synced from the <span class="font-medium text-accent">Water/Lipid</span> decomposition — run that tab first if the maps are empty. The pericoronary lumen wall is detected automatically.
      </div>
    {/if}
    <!-- Main content: editor + surface plot side-by-side -->
    <div class="flex min-h-0 flex-1 overflow-hidden">
      <!-- Left: SnakeEditor for the selected cross-section -->
      <div class="flex w-1/2 shrink-0 flex-col overflow-hidden border-r border-border">
        {#if currentTarget}
          <SnakeEditor
            target={currentTarget}
            snakePoints={currentSnake}
            {arcOffsetMm}
            overlay={currentOverlay}
            {material}
            {unit}
            onStepTarget={(delta) => handleSelect(selectedIndex + delta)}
          />
        {/if}
      </div>

      <!-- Right: surface plot + radial profile placeholder -->
      <div class="flex min-h-0 min-w-0 flex-1 flex-col overflow-hidden">
        <SurfacePlotPanel
          {surfaces}
          selectedIndex={selectedSurfaceIndex}
          {material}
          {unit}
          onSliderChange={(si) => handleSelect(surfaceToTargetIndex(si))}
          {arcOffsetMm}
        />

        <!-- Radial profile placeholder -->
        <div class="shrink-0 border-t border-border px-3 py-2">
          <span class="text-[10px] text-text-secondary/60">Radial profile (coming soon)</span>
        </div>
      </div>
    </div>

    <!-- Cross-section strip (full width) -->
    <div class="shrink-0 border-t border-border bg-surface-secondary">
      <CrossSectionStrip
        {targets}
        {selectedIndex}
        {statusMap}
        onSelect={handleSelect}
        {arcOffsetMm}
      />
    </div>

    <!-- Bottom toolbar: overlay chips (left) + MMD actions (right) -->
    <div class="flex shrink-0 flex-wrap items-center gap-3 border-t border-border bg-surface-secondary px-3 py-1.5">
      <OverlaySelector
        {material}
        {unit}
        onMaterialChange={handleMaterialChange}
        onUnitChange={handleUnitChange}
      />

      <div class="ml-auto flex items-center gap-2">
        <button
          class="rounded bg-surface-tertiary px-3 py-1 text-xs font-medium text-text-primary hover:bg-surface-tertiary/80 active:bg-surface-tertiary/60 disabled:bg-surface-tertiary/40 disabled:text-text-secondary/70"
          onclick={handleSave}
          disabled={saveBusy || targets.length === 0}
          title="Save annotation state for this patient"
        >
          {saveBusy ? 'Saving...' : 'Save'}
        </button>

        <button
          class="rounded bg-surface-tertiary px-3 py-1 text-xs font-medium text-text-primary hover:bg-surface-tertiary/80 active:bg-surface-tertiary/60 disabled:bg-surface-tertiary/40 disabled:text-text-secondary/70"
          onclick={handleExportCsv}
          disabled={exportBusy || !mmdSummary}
          title="Export MMD surface data as CSV"
        >
          {exportBusy ? 'Exporting...' : 'Export CSV'}
        </button>

        <span class="text-[11px] tabular-nums text-text-secondary">
          {contourCount}/{totalCount}
        </span>

        {#if saveMsg}
          <span class="text-[10px] text-success">{saveMsg}</span>
        {/if}

        {#if mmdBusy}
          <span class="text-[11px] text-text-secondary">Syncing from Water/Lipid…</span>
        {:else if mmdSummary}
          <span class="text-[10px] text-success" title="Synced from the Water/Lipid calibration">
            Synced ({mmdSummary.n_voxels.toLocaleString()} voxels)
          </span>
        {/if}

        {#if mmdError}
          <span class="max-w-[28ch] truncate text-[10px] text-error" title={mmdError}>
            {mmdError.includes('calibration') ? 'Run the Water/Lipid tab first' : mmdError}
          </span>
        {/if}
      </div>
    </div>
  {/if}
</div>

<script lang="ts">
  /**
   * Root application shell.
   *
   * Layout:
   *   - Header toolbar (title + actions)
   *   - Tab bar (Editor | MMD Analysis)
   *   - Main area    (MprPanel or MmdAnalysisView fills remaining space)
   *   - Footer status bar (loading progress / ready state)
   */
  import MprPanel from './components/MprPanel.svelte';
  import MmdAnalysisView from './components/MmdAnalysisView.svelte';
  import WaterLipidView from './components/WaterLipidView.svelte';
  import PatientBrowser from './components/PatientBrowser.svelte';
  import SeedToolbar from './components/SeedToolbar.svelte';
  import HintLine from './components/HintLine.svelte';
  import ProgressOverlay from './components/ProgressOverlay.svelte';
  import {
    openDicomDialog,
    getRecentDicoms,
    loadSeeds,
    scanSeries,
    loadSeries,
    loadDualEnergy,
    loadPatientAll,
    setActiveVolume,
    setActiveVolumeMeta,
    onDicomLoadProgress,
    reuseLoadedVolume,
  } from '$lib/api';
  import type { LoadedSeriesDescriptor } from '$lib/api';
  import { loadSession, clearSession, saveSession } from '$lib/session';
  import { cache as cornerstoneCache } from '@cornerstonejs/core';
  import { buildVolume } from '$lib/cornerstone/volumeLoader';
  import { volumeStore } from '$lib/stores/volumeStore.svelte';
  import type { VolumeMetadata } from '$lib/stores/volumeStore.svelte';
  import { pipelineStore } from '$lib/stores/pipelineStore.svelte';
  import { seedStore, type Vessel } from '$lib/stores/seedStore.svelte';
  import { uiStore } from '$lib/stores/uiStore.svelte';
  import { navigateToWorldPos } from '$lib/navigation';
  import { mmdStore } from '$lib/stores/mmdStore.svelte';
  import { flagStore } from '$lib/stores/flagStore.svelte';
  import { derivePatientStatus, patientIdOf } from '$lib/patientStatus';

  /* ── Tab state ─────────────────────────────────────── */
  type AppTab = 'editor' | 'mmd' | 'wl';
  let activeTab = $state<AppTab>('editor');

  /** Centerline of the currently active vessel (for MmdAnalysisView). */
  let activeCenterlineMm = $derived(seedStore.activeVesselData.centerline ?? []);

  let errorMessage = $state('');
  let showFlagDialog = $state(false);
  let flagNoteDraft = $state('');
  let recentPaths = $state<string[]>([]);
  let showRecent = $state(false);
  let showPatientBrowser = $state(false);

  // Load recent paths on mount
  $effect(() => {
    getRecentDicoms().then((paths) => { recentPaths = paths; }).catch(() => {});
  });

  /** Auto-restore the full saved session (seeds + FAI + water/lipid) for a
   *  patient on open. Falls back to legacy seeds-only saves so patients saved
   *  before the unified session still restore their seeds. */
  async function autoRestoreSession(dicomPath: string) {
    // Clean slate per patient — never carry one patient's seeds / FAI / MMD /
    // water-lipid into the next. loadSession restores this patient's saved state.
    clearSession();
    try {
      const meta = await loadSession(dicomPath);
      if (!meta) {
        const seedsJson = await loadSeeds(dicomPath);
        if (seedsJson) seedStore.importJson(seedsJson);
      }
    } catch { /* no saved session/seeds for this patient */ }
    await flagStore.load(dicomPath);
  }

  /** Run FAI, then auto-persist the session so the patient list shows progress
   *  without a manual Save. Footer status updates live regardless. */
  async function runFaiAndSave() {
    await pipelineStore.run();
    if (pipelineStore.status === 'complete' && volumeStore.dicomPath) {
      try {
        await saveSession(volumeStore.dicomPath);
      } catch (e) {
        console.error('auto-save after FAI failed:', e);
      }
    }
  }

  // Clear stale FAI/pipeline results (and their overlay) when the centerline is
  // removed — deleting seeds, Escape, or switching patient. Without this the FAI
  // heatmap lingers on screen after the centerline it was computed from is gone.
  $effect(() => {
    const hasCenterline = Object.values(seedStore.vessels).some(
      (v) => v.centerline !== null && v.centerline.length >= 2,
    );
    if (!hasCenterline && pipelineStore.results !== null) {
      pipelineStore.reset();
    }
  });

  // ---- Keyboard shortcuts ----
  function handleKeydown(event: KeyboardEvent) {
    // Escape closes the fullscreen analysis overlay FIRST and stops there — it
    // must not also run the "clear vessel" shortcut below. Both fire on this one
    // window keydown, and the overlay's own stopPropagation can't prevent it, so
    // the guard lives here, at the single owner of the flag.
    if (event.key === 'Escape' && uiStore.analysisMaximized) {
      uiStore.analysisMaximized = false;
      return;
    }

    // Ignore if user is typing in an input/textarea
    const tag = (event.target as HTMLElement)?.tagName;
    if (tag === 'INPUT' || tag === 'TEXTAREA' || tag === 'SELECT') return;

    // Cmd+Shift+Z / Ctrl+Shift+Z: redo
    if ((event.ctrlKey || event.metaKey) && event.shiftKey && event.key === 'z') {
      event.preventDefault();
      seedStore.redo();
      return;
    }

    // Cmd+Z / Ctrl+Z: undo
    if ((event.ctrlKey || event.metaKey) && event.key === 'z') {
      event.preventDefault();
      seedStore.undo();
      return;
    }

    // Backspace / Delete: if seed selected, delete selected; else delete last
    if (event.key === 'Backspace' || event.key === 'Delete') {
      event.preventDefault();
      const selected = seedStore.selectedSeedIndex;
      if (selected !== null) {
        seedStore.removeSeed(selected);
      } else {
        const data = seedStore.activeVesselData;
        if (data.seeds.length > 0) {
          seedStore.removeSeed(data.seeds.length - 1);
        }
      }
      return;
    }

    // Escape: if seed selected, deselect; else clear vessel
    if (event.key === 'Escape') {
      if (seedStore.selectedSeedIndex !== null) {
        seedStore.deselectSeed();
      } else {
        seedStore.clearVessel(seedStore.activeVessel);
      }
      return;
    }

    // Arrow Left/Right: cycle through seeds
    if (event.key === 'ArrowLeft' || event.key === 'ArrowRight') {
      const data = seedStore.activeVesselData;
      if (data.seeds.length === 0) return;
      const current = seedStore.selectedSeedIndex;
      let next: number;
      if (current === null) {
        next = event.key === 'ArrowRight' ? 0 : data.seeds.length - 1;
      } else {
        next = event.key === 'ArrowRight'
          ? (current + 1) % data.seeds.length
          : (current - 1 + data.seeds.length) % data.seeds.length;
      }
      seedStore.selectSeed(next);
      navigateToWorldPos(data.seeds[next].position);
      return;
    }

    // 1/2/3: switch active vessel
    const vesselMap: Record<string, Vessel> = { '1': 'RCA', '2': 'LAD', '3': 'LCx' };
    if (vesselMap[event.key]) {
      seedStore.setActiveVessel(vesselMap[event.key]);
      return;
    }
  }

  /** Load DICOM from a specific folder path. */
  async function loadFromPath(path: string) {
    errorMessage = '';
    showRecent = false;

    // Idempotency guard: re-clicking the same patient (from recents, the
    // patient browser, or the same folder via the dialog) should be a no-op
    // rather than re-running ~5 s of SMB scan + decode. We also preserve
    // unsaved in-memory seed edits, which a full reload would wipe.
    //
    // Force-reload workflow: pick any other patient, then pick this one
    // again. The guard only short-circuits when a successful prior load is
    // still resident (volume + cornerstone id present) and no load is in
    // flight.
    if (
      volumeStore.current?.dicomPath === path
      && volumeStore.cornerstoneVolumeId !== null
      && !volumeStore.loading
    ) {
      showPatientBrowser = false;
      return;
    }

    let unlistenProgress: (() => void) | null = null;

    try {
      // Clear previous state
      seedStore.clearAll();
      volumeStore.clear();

      volumeStore.setLoading(true);
      volumeStore.setLoadProgress(0);

      // Subscribe to progress events. Decoding is the long leg; map its
      // `done / total` to 0-95% so the bar lands at 100% after buildVolume.
      unlistenProgress = await onDicomLoadProgress((p) => {
        if (p.phase === 'decoding' && p.total > 0) {
          volumeStore.setLoadProgress(Math.round((p.done / p.total) * 95));
        } else if (p.phase === 'done') {
          volumeStore.setLoadProgress(95);
        }
      });

      // 1. Header-only scan to discover series.
      const series = await scanSeries(path);
      if (series.length === 0) {
        throw new Error('No DICOM series found in folder.');
      }
      if (series.length > 1) {
        console.warn(
          `Folder has ${series.length} series; auto-selecting first:`,
          series.map((s) => `${s.description} (${s.num_slices} slices)`),
        );
      }
      const chosen = series[0];

      // Fast reload path: if both cornerstone (JS) and the Rust AppState already
      // hold this exact (path, uid), skip the decode+IPC entirely and just rewire
      // the store. Saves 30-70s on A→B→A workflows.
      const fastVolumeKey = `${path}::${chosen.uid}`;
      const fastCsId = `pcat:${fastVolumeKey}`;
      const cachedVolume = cornerstoneCache.getVolume(fastCsId);
      if (cachedVolume) {
        const cachedMeta = await reuseLoadedVolume(path, chosen.uid);
        if (cachedMeta) {
          // Both sides already have it. Rebuild frontend store state from metadata.
          const direction = computeDirectionMatrix(cachedMeta.orientation);
          const ipp = cachedMeta.image_position_patient;
          const storeMeta: VolumeMetadata = {
            volumeId: fastVolumeKey,
            shape: [cachedMeta.num_slices, cachedMeta.rows, cachedMeta.cols],
            spacing: [cachedMeta.slice_spacing, cachedMeta.pixel_spacing[0], cachedMeta.pixel_spacing[1]],
            origin: [cachedMeta.slice_positions_z[0] ?? ipp[2], ipp[1], ipp[0]],
            direction,
            windowCenter: cachedMeta.window_center,
            windowWidth: cachedMeta.window_width,
            patientName: cachedMeta.patient_name,
            studyDescription: cachedMeta.study_description,
            dicomPath: path,
          };
          volumeStore.set(storeMeta);
          volumeStore.setCornerstoneVolumeId(fastCsId);
          volumeStore.setLoadProgress(100);
          volumeStore.setLoading(false);
          await autoRestoreSession(path);
          getRecentDicoms().then((paths) => { recentPaths = paths; }).catch(() => {});
          return;
        }
      }

      // 2. Bulk load the chosen series (one binary IPC trip).
      const { metadata, voxels } = await loadSeries(path, chosen.uid);

      // 3. Build the cornerstone3D volume synchronously.
      // Stable key so cornerstone's cache.getVolume short-circuit fires on reload.
      const volumeKey = `${path}::${chosen.uid}`;
      const csId = buildVolume(volumeKey, metadata, voxels);

      // 4. Populate the legacy volumeStore shape.
      const direction = computeDirectionMatrix(metadata.orientation);
      const ipp = metadata.image_position_patient;
      const storeMeta: VolumeMetadata = {
        volumeId: volumeKey,
        shape: [metadata.num_slices, metadata.rows, metadata.cols],
        spacing: [metadata.slice_spacing, metadata.pixel_spacing[0], metadata.pixel_spacing[1]],
        // ZYX patient LPS mm; mirrors src-tauri bridge_into_state.
        origin: [metadata.slice_positions_z[0] ?? ipp[2], ipp[1], ipp[0]],
        direction,
        windowCenter: metadata.window_center,
        windowWidth: metadata.window_width,
        patientName: metadata.patient_name,
        studyDescription: metadata.study_description,
        dicomPath: path,
      };
      volumeStore.set(storeMeta);
      volumeStore.setCornerstoneVolumeId(csId);
      volumeStore.setLoadProgress(100);
      volumeStore.setLoading(false);

      // 5. Auto-load seeds for this patient.
      await autoRestoreSession(path);

      // 6. Refresh recent list.
      getRecentDicoms().then((paths) => { recentPaths = paths; }).catch(() => {});
    } catch (e) {
      volumeStore.setLoading(false);
      errorMessage = e instanceof Error ? e.message : String(e);
      console.error('Failed to load DICOM:', e);
    } finally {
      if (unlistenProgress) unlistenProgress();
    }
  }

  /**
   * Build a 3x3 direction matrix (row-major, 9 elements) from a 6-element
   * ImageOrientationPatient vector (row direction + column direction); the third
   * row is the cross product (slice normal).
   */
  function computeDirectionMatrix(orient: [number, number, number, number, number, number]): number[] {
    const row: [number, number, number] = [orient[0], orient[1], orient[2]];
    const col: [number, number, number] = [orient[3], orient[4], orient[5]];
    const normal: [number, number, number] = [
      row[1] * col[2] - row[2] * col[1],
      row[2] * col[0] - row[0] * col[2],
      row[0] * col[1] - row[1] * col[0],
    ];
    return [
      row[0], row[1], row[2],
      col[0], col[1], col[2],
      normal[0], normal[1], normal[2],
    ];
  }

  /** Open DICOM folder picker, then load. */
  async function handleOpenDicom() {
    const path = await openDicomDialog();
    if (!path) return;
    loadFromPath(path);
  }

  /** Load a low/high keV pair as a dual-energy volume. Triggered from the
   *  patient browser when the user picks a MonoPlus keV series that has a
   *  sibling at a different keV — the browser does the auto-pairing. */
  async function loadDualEnergyPair(lowDir: string, highDir: string) {
    errorMessage = '';
    let unlistenProgress: (() => void) | undefined;
    try {
      volumeStore.clear();
      volumeStore.setLoading(true);
      volumeStore.setLoadProgress(0);

      unlistenProgress = await onDicomLoadProgress((p) => {
        if (p.phase === 'decoding' && p.total > 0) {
          volumeStore.setLoadProgress(Math.round((p.done / p.total) * 95));
        } else if (p.phase === 'done') {
          volumeStore.setLoadProgress(95);
        }
      });

      // Backend parses keV + loads both series + populates state.dual_energy
      // + mirrors low into state.volume. Returns the framed low-energy bundle
      // so we can build the cornerstone volume without a second fetch.
      const { metadata, voxels } = await loadDualEnergy(lowDir, highDir);

      const volumeKey = `${lowDir}::dual_energy`;
      const csId = buildVolume(volumeKey, metadata, voxels);

      const direction = computeDirectionMatrix(metadata.orientation);
      const ipp = metadata.image_position_patient;
      const storeMeta: VolumeMetadata = {
        volumeId: volumeKey,
        shape: [metadata.num_slices, metadata.rows, metadata.cols],
        spacing: [metadata.slice_spacing, metadata.pixel_spacing[0], metadata.pixel_spacing[1]],
        origin: [metadata.slice_positions_z[0] ?? ipp[2], ipp[1], ipp[0]],
        direction,
        windowCenter: metadata.window_center,
        windowWidth: metadata.window_width,
        patientName: metadata.patient_name,
        studyDescription: metadata.study_description,
        dicomPath: lowDir,
      };
      volumeStore.set(storeMeta);
      volumeStore.setCornerstoneVolumeId(csId);
      volumeStore.setLoadProgress(100);
      volumeStore.setLoading(false);

      await autoRestoreSession(lowDir);

      getRecentDicoms().then((paths) => { recentPaths = paths; }).catch(() => {});
    } catch (e) {
      volumeStore.setLoading(false);
      errorMessage = e instanceof Error ? e.message : String(e);
      console.error('Failed to load dual-energy:', e);
    } finally {
      if (unlistenProgress) unlistenProgress();
    }
  }

  /** Load every DICOM series under a patient folder into the Rust cache.
   *  Hydrates the cornerstone volume for the chosen "active" series. After
   *  this, the volume switcher in the header swaps between siblings
   *  instantly (no re-decode). */
  async function loadPatientFolder(patientPath: string) {
    errorMessage = '';
    let unlistenProgress: (() => void) | undefined;
    try {
      volumeStore.clear();
      volumeStore.setLoading(true);
      volumeStore.setLoadProgress(0);
      volumeStore.setLoadMessage('');

      // Two-phase lazy load: header SCAN of all series (first 15%), then DECODE
      // of the active + dual-energy series (15-95%, the long leg). buildVolume
      // carries it to 100% afterward. Decode progress = coarse per-series index
      // ("patient_series") + fine per-slice ("decoding").
      const SCAN_PCT = 15;
      const DECODE_END = 95;
      let seriesIdx = 0;
      let seriesTotal = 0;
      const scaleToOuter = (withinSeries: number): number => {
        if (seriesTotal === 0) return SCAN_PCT;
        const span = (DECODE_END - SCAN_PCT) / seriesTotal;
        return Math.min(DECODE_END, Math.round(SCAN_PCT + seriesIdx * span + withinSeries * span));
      };

      unlistenProgress = await onDicomLoadProgress((p) => {
        if (p.phase === 'scanning' && p.total > 0) {
          volumeStore.setLoadProgress(Math.round((p.done / p.total) * SCAN_PCT));
          if (p.detail) {
            volumeStore.setLoadMessage(`Scanning ${p.done + 1}/${p.total}: ${p.detail}`);
          }
        } else if (p.phase === 'patient_series' && p.total > 0) {
          // A needed series is starting to decode.
          seriesIdx = p.done;
          seriesTotal = p.total;
          if (p.detail) {
            volumeStore.setLoadMessage(`Decoding ${p.done + 1}/${p.total}: ${p.detail}`);
          }
          volumeStore.setLoadProgress(scaleToOuter(0));
        } else if (p.phase === 'decoding' && p.total > 0 && seriesTotal > 0) {
          // Fine-grained slice decode within the current series.
          volumeStore.setLoadProgress(scaleToOuter(p.done / p.total));
        } else if (p.phase === 'done') {
          volumeStore.setLoadProgress(DECODE_END);
          volumeStore.setLoadMessage('Finalizing');
        }
      });

      const _t0 = performance.now();
      const result = await loadPatientAll(patientPath);
      const _t1 = performance.now();

      // Hydrate cornerstone for the active series via the cache-hit path.
      const active = result.series[result.active_index];
      const { metadata, voxels } = await setActiveVolume(active.path, active.uid);
      const _t2 = performance.now();

      const volumeKey = `${active.path}::${active.uid}`;
      const csId = buildVolume(volumeKey, metadata, voxels);
      const _t3 = performance.now();
      console.log(
        `[load-timing] loadPatientAll(scan+decode)=${(_t1 - _t0) | 0}ms · ` +
        `setActiveVolume(IPC ${(voxels.length * 2 / 1048576) | 0}MB)=${(_t2 - _t1) | 0}ms · ` +
        `buildVolume(cornerstone)=${(_t3 - _t2) | 0}ms`,
      );

      const direction = computeDirectionMatrix(metadata.orientation);
      const ipp = metadata.image_position_patient;
      const storeMeta: VolumeMetadata = {
        volumeId: volumeKey,
        shape: [metadata.num_slices, metadata.rows, metadata.cols],
        spacing: [metadata.slice_spacing, metadata.pixel_spacing[0], metadata.pixel_spacing[1]],
        origin: [metadata.slice_positions_z[0] ?? ipp[2], ipp[1], ipp[0]],
        direction,
        windowCenter: metadata.window_center,
        windowWidth: metadata.window_width,
        patientName: metadata.patient_name,
        studyDescription: metadata.study_description,
        dicomPath: active.path,
      };
      volumeStore.set(storeMeta);
      volumeStore.setCornerstoneVolumeId(csId);
      volumeStore.setLoaded(
        result.series.map((s) => ({
          name: s.name,
          path: s.path,
          uid: s.uid,
          seriesDescription: s.series_description,
          kev: s.kev,
          numSlices: s.num_slices,
          rows: s.rows,
          cols: s.cols,
        })),
      );
      volumeStore.setLoadProgress(100);
      volumeStore.setLoadMessage('');
      volumeStore.setLoading(false);

      await autoRestoreSession(active.path);

      if (result.failures.length > 0) {
        errorMessage = `Loaded ${result.series.length} series; skipped: ${result.failures.join('; ')}`;
      }

      getRecentDicoms().then((paths) => { recentPaths = paths; }).catch(() => {});
    } catch (e) {
      volumeStore.setLoading(false);
      volumeStore.setLoadMessage('');
      errorMessage = e instanceof Error ? e.message : String(e);
      console.error('Failed to load patient:', e);
    } finally {
      if (unlistenProgress) unlistenProgress();
    }
  }

  /** Switch the active volume to one already in the Rust cache. Re-hydrates
   *  the cornerstone3D volume (fast: no decode) and updates volumeStore. */
  async function switchToLoaded(entry: LoadedSeriesDescriptor | typeof volumeStore.loaded[number]) {
    if (volumeStore.loading) return;
    errorMessage = '';
    try {
      volumeStore.setLoading(true);
      volumeStore.setLoadProgress(0);

      const volumeKey = `${entry.path}::${entry.uid}`;
      const csVolumeId = `pcat:${volumeKey}`;

      // When cornerstone3D already holds this volume's pixels, skip the
      // 150–220 MB voxel IPC transfer entirely — fetch only the metadata and
      // re-bind the cached GPU volume. Switching between already-viewed series
      // (CCTA ⇄ CaScore ⇄ keV) is then a metadata round-trip, not a full
      // re-marshal of the volume. Only a genuine cornerstone miss pays for the
      // voxels.
      // `metadata` is the API series metadata (series_uid, num_slices, …),
      // distinct from the store's VolumeMetadata built below. Derive its type
      // from the API call to avoid the name collision with the store type.
      let metadata: Awaited<ReturnType<typeof setActiveVolumeMeta>>;
      let csId: string;
      if (cornerstoneCache.getVolume(csVolumeId)) {
        metadata = await setActiveVolumeMeta(entry.path, entry.uid);
        csId = csVolumeId;
      } else {
        const res = await setActiveVolume(entry.path, entry.uid);
        metadata = res.metadata;
        csId = buildVolume(volumeKey, metadata, res.voxels);
      }

      const direction = computeDirectionMatrix(metadata.orientation);
      const ipp = metadata.image_position_patient;
      const storeMeta: VolumeMetadata = {
        volumeId: volumeKey,
        shape: [metadata.num_slices, metadata.rows, metadata.cols],
        spacing: [metadata.slice_spacing, metadata.pixel_spacing[0], metadata.pixel_spacing[1]],
        origin: [metadata.slice_positions_z[0] ?? ipp[2], ipp[1], ipp[0]],
        direction,
        windowCenter: metadata.window_center,
        windowWidth: metadata.window_width,
        patientName: metadata.patient_name,
        studyDescription: metadata.study_description,
        dicomPath: entry.path,
      };
      volumeStore.set(storeMeta);
      volumeStore.setCornerstoneVolumeId(csId);
      volumeStore.setLoadProgress(100);
    } catch (e) {
      errorMessage = e instanceof Error ? e.message : String(e);
      console.error('Failed to switch volume:', e);
    } finally {
      volumeStore.setLoading(false);
    }
  }
</script>

<svelte:window onkeydown={handleKeydown} onclick={() => { showRecent = false; }} />

<div class="flex h-screen flex-col">
  <!-- ===== Header toolbar ===== -->
  <header
    class="flex h-11 shrink-0 items-center justify-between border-b border-border bg-surface-secondary px-4"
  >
    <div class="flex items-center gap-2">
      <h1 class="text-sm font-semibold tracking-wide text-text-primary">
        PCAT Workstation
      </h1>
      <span class="text-[11px] text-text-secondary">v2.0-dev</span>
    </div>

    <!-- Seed vessel selector (visible after volume load) -->
    {#if volumeStore.current}
      <div class="flex items-center gap-3">
        <SeedToolbar />
      </div>
    {/if}

    <!-- Volume switcher: shown only when a full patient has been loaded and
         multiple volumes are resident in the cache. Clicking a chip swaps
         the active volume (MPR/CPR/FAI re-bind) without another decode. -->
    {#if volumeStore.loaded.length > 1 && volumeStore.current}
      <div class="flex items-center gap-1 overflow-x-auto px-2" title="Active volume — click to switch">
        {#each volumeStore.loaded as entry (entry.path + entry.uid)}
          {@const isActive = volumeStore.current?.dicomPath === entry.path}
          <button
            class="shrink-0 rounded-full px-2.5 py-1 text-[10px] font-medium transition-colors
                   {isActive
                     ? 'bg-accent text-white'
                     : 'bg-surface-tertiary text-text-secondary hover:text-text-primary'}"
            onclick={() => !isActive && switchToLoaded(entry)}
            title={entry.seriesDescription || entry.name}
          >
            {entry.kev !== null ? `${entry.kev} keV` : entry.name}
          </button>
        {/each}
      </div>
    {/if}

    <div class="flex items-center gap-1.5">
      <!-- Pipeline action button -->
      {#if pipelineStore.status === 'complete'}
        <button
          class="rounded bg-accent/10 px-3 py-1 text-xs font-medium text-accent hover:bg-accent/20"
          onclick={() => { runFaiAndSave(); }}
          title="Re-run: centerline → contour extraction → CRISP-CT VOI (1mm gap + 3mm ring) → FAI stats"
        >
          Re-analyze
        </button>
      {:else if pipelineStore.canRun}
        <button
          class="rounded px-3 py-1 text-xs font-medium text-accent hover:bg-accent/10 active:bg-accent/20 disabled:opacity-40"
          onclick={() => runFaiAndSave()}
          disabled={pipelineStore.status === 'running'}
          title="Run FAI pipeline: centerline → contour extraction → CRISP-CT VOI (1mm gap + 3mm ring) → FAI stats"
        >
          {pipelineStore.status === 'running' ? 'Analyzing...' : 'Analyze'}
        </button>
      {/if}

      <button
        class="rounded px-3 py-1 text-xs font-medium text-accent hover:bg-accent/10 active:bg-accent/20 disabled:opacity-40"
        onclick={(e: MouseEvent) => { e.stopPropagation(); showPatientBrowser = true; }}
        disabled={volumeStore.loading}
        title="Browse patients in cohort directory"
      >
        Patients
      </button>

      <div class="relative flex items-center">
        <button
          class="rounded-l px-3 py-1 text-xs font-medium text-accent hover:bg-accent/10 active:bg-accent/20 disabled:opacity-40"
          onclick={handleOpenDicom}
          disabled={volumeStore.loading}
        >
          Open DICOM
        </button>
        {#if recentPaths.length > 0}
          <button
            class="rounded-r border-l border-border px-1.5 py-1 text-xs text-accent hover:bg-accent/10 disabled:opacity-40"
            onclick={(e: MouseEvent) => { e.stopPropagation(); showRecent = !showRecent; }}
            disabled={volumeStore.loading}
            title="Recent files"
          >
            &#9662;
          </button>
        {/if}

        <!-- Recent files dropdown -->
        {#if showRecent && recentPaths.length > 0}
          <!-- svelte-ignore a11y_no_static_element_interactions -->
          <div
            class="absolute right-0 top-full z-50 mt-1 max-h-64 w-80 overflow-y-auto rounded border border-border bg-surface-secondary shadow-lg"
          >
            <div class="px-3 py-1.5 text-[10px] font-semibold uppercase tracking-wider text-text-secondary/60">
              Recent
            </div>
            {#each recentPaths as rp}
              <button
                class="w-full px-3 py-1.5 text-left text-[11px] text-text-primary hover:bg-accent/10 truncate"
                onclick={() => loadFromPath(rp)}
                title={rp}
              >
                {rp.split('/').slice(-2).join('/')}
              </button>
            {/each}
          </div>
        {/if}
      </div>
    </div>
  </header>

  <!-- ===== Tab bar ===== -->
  {#if volumeStore.current}
    <nav class="flex shrink-0 items-center gap-1 border-b border-border bg-surface-secondary px-4">
      <button
        class="relative px-3 py-1.5 text-xs font-medium transition-colors {activeTab === 'editor'
          ? 'text-accent'
          : 'text-text-secondary hover:text-text-primary'}"
        onclick={() => { activeTab = 'editor'; }}
      >
        Editor
        {#if activeTab === 'editor'}
          <span class="absolute inset-x-0 bottom-0 h-[2px] bg-accent"></span>
        {/if}
      </button>
      <button
        class="relative px-3 py-1.5 text-xs font-medium transition-colors {activeTab === 'wl'
          ? 'text-accent'
          : 'text-text-secondary hover:text-text-primary'}"
        onclick={() => { activeTab = 'wl'; }}
        title="Whole-volume noise-aware GLS water/lipid decomposition — run this first"
      >
        Water/Lipid
        {#if activeTab === 'wl'}
          <span class="absolute inset-x-0 bottom-0 h-[2px] bg-accent"></span>
        {/if}
      </button>
      <button
        class="relative px-3 py-1.5 text-xs font-medium transition-colors {activeTab === 'mmd'
          ? 'text-accent'
          : 'text-text-secondary hover:text-text-primary'}"
        onclick={() => { activeTab = 'mmd'; }}
        title="Pericoronary cross-section + 3D surface, synced from the Water/Lipid decomposition"
      >
        MMD Analysis
        {#if activeTab === 'mmd'}
          <span class="absolute inset-x-0 bottom-0 h-[2px] bg-accent"></span>
        {/if}
      </button>
    </nav>
  {/if}

  <!-- ===== Main viewport area =====
       Both panes stay mounted; we toggle visibility instead of using {#if}
       so tab switches don't destroy component state (MMD results, cornerstone
       MIP viewport bindings, contour edits, etc.). Each pane is absolutely
       positioned inside the relative <main> so hiding one doesn't collapse
       the other's layout. `display: none` also stops offscreen work on
       hidden canvases. -->
  <main class="relative min-h-0 flex-1">
    {#if volumeStore.current}
      <div
        class="absolute inset-0 flex flex-col"
        class:hidden={activeTab !== 'editor'}
      >
        <MprPanel />
        <HintLine />
        {#if pipelineStore.status === 'running'}
          <ProgressOverlay />
        {/if}
      </div>
      <div
        class="absolute inset-0 flex flex-col"
        class:hidden={activeTab !== 'mmd'}
      >
        <MmdAnalysisView centerlineMm={activeCenterlineMm} />
      </div>
      <div
        class="absolute inset-0 flex flex-col"
        class:hidden={activeTab !== 'wl'}
      >
        <WaterLipidView />
      </div>
    {/if}
  </main>

  <!-- ===== Patient browser modal ===== -->
  {#if showPatientBrowser}
    <PatientBrowser
      currentPath={volumeStore.dicomPath}
      onSelect={(path) => { showPatientBrowser = false; loadFromPath(path); }}
      onSelectDualEnergy={(lowDir, highDir) => {
        showPatientBrowser = false;
        loadDualEnergyPair(lowDir, highDir);
      }}
      onSelectPatient={(path) => { showPatientBrowser = false; loadPatientFolder(path); }}
      onClose={() => { showPatientBrowser = false; }}
    />
  {/if}

  <!-- ===== Flag dialog ===== -->
  {#if showFlagDialog}
    <!-- svelte-ignore a11y_click_events_have_key_events -->
    <!-- svelte-ignore a11y_no_static_element_interactions -->
    <div
      class="absolute inset-0 z-50 flex items-center justify-center bg-black/40"
      onclick={(e) => { if (e.target === e.currentTarget) showFlagDialog = false; }}
    >
      <div class="w-80 rounded-lg border border-border bg-surface-secondary p-4 shadow-xl">
        <h3 class="mb-1 text-sm font-semibold text-text-primary">Flag this data</h3>
        <p class="mb-2 text-[11px] text-text-secondary">
          Mark this patient's data as problematic / unusable — e.g. vessel
          discontinuity, motion or recon artifact.
        </p>
        <textarea
          bind:value={flagNoteDraft}
          rows="3"
          placeholder="Optional note (what's wrong)…"
          class="w-full rounded border border-border bg-surface px-2 py-1 text-xs text-text-primary focus:border-accent focus:outline-none"
        ></textarea>
        <div class="mt-3 flex items-center justify-between gap-2">
          {#if flagStore.current?.flagged}
            <button
              class="rounded px-2 py-1 text-[11px] text-text-secondary hover:bg-surface-tertiary"
              onclick={async () => {
                if (volumeStore.dicomPath) await flagStore.set(volumeStore.dicomPath, false, '');
                showFlagDialog = false;
              }}
            >
              Remove flag
            </button>
          {:else}
            <span></span>
          {/if}
          <div class="flex gap-2">
            <button
              class="rounded px-2 py-1 text-[11px] text-text-secondary hover:bg-surface-tertiary"
              onclick={() => (showFlagDialog = false)}
            >
              Cancel
            </button>
            <button
              class="rounded bg-error/15 px-3 py-1 text-[11px] font-medium text-error hover:bg-error/25"
              onclick={async () => {
                if (volumeStore.dicomPath)
                  await flagStore.set(volumeStore.dicomPath, true, flagNoteDraft.trim());
                showFlagDialog = false;
              }}
            >
              Flag
            </button>
          </div>
        </div>
      </div>
    </div>
  {/if}

  <!-- ===== Footer status bar ===== -->
  <footer
    class="flex h-6 shrink-0 items-center justify-between border-t border-border bg-surface-secondary px-4"
  >
    <div class="flex items-center gap-2">
      {#if volumeStore.loading}
        <!-- Loading progress bar -->
        <div class="flex items-center gap-2">
          <div class="h-1.5 w-28 overflow-hidden rounded-full bg-surface-tertiary">
            <div
              class="h-full rounded-full bg-accent transition-all duration-150 ease-out"
              style="width: {volumeStore.loadProgress}%"
            ></div>
          </div>
          <span class="truncate text-[11px] text-text-secondary">
            {volumeStore.loadMessage || 'Loading volume...'} <span class="tabular-nums">{volumeStore.loadProgress}%</span>
          </span>
        </div>
      {:else if pipelineStore.status === 'error'}
        <span class="h-1.5 w-1.5 rounded-full bg-error"></span>
        <span class="truncate text-[11px] text-error">Analysis: {pipelineStore.error}</span>
      {:else if errorMessage}
        <span class="h-1.5 w-1.5 rounded-full bg-error"></span>
        <span class="truncate text-[11px] text-error">{errorMessage}</span>
      {:else if volumeStore.current}
        {@const st = derivePatientStatus({
          hasSeeds: (['LAD', 'LCx', 'RCA'] as const).some(
            (v) => seedStore.vessels[v].seeds.length > 0,
          ),
          hasFai: pipelineStore.results !== null,
          hasMmd: mmdStore.summary !== null,
        })}
        <span class="h-1.5 w-1.5 rounded-full {st === 'complete' ? 'bg-success' : st === 'in_progress' ? 'bg-warning' : 'bg-text-secondary/40'}"></span>
        <span class="text-[11px] font-medium text-text-primary">{patientIdOf(volumeStore.dicomPath)}</span>
        {#if volumeStore.current.studyDescription}
          <span class="truncate text-[11px] text-text-secondary">· {volumeStore.current.studyDescription}</span>
        {/if}
        <span
          class="rounded px-1.5 text-[10px] font-medium {st === 'complete'
            ? 'bg-success/15 text-success'
            : st === 'in_progress'
              ? 'bg-warning/15 text-warning'
              : 'bg-text-secondary/15 text-text-secondary'}"
        >
          {st === 'complete' ? 'complete' : st === 'in_progress' ? 'in progress' : 'not started'}
        </span>
        {#if flagStore.current?.flagged}
          <button
            class="rounded bg-error/15 px-1.5 text-[10px] font-medium text-error hover:bg-error/25"
            title={flagStore.current.note || 'Flagged: data problem / unusable'}
            onclick={() => { flagNoteDraft = flagStore.current?.note ?? ''; showFlagDialog = true; }}
          >
            ⚠ Flagged
          </button>
        {:else}
          <button
            class="rounded px-1.5 text-[10px] text-text-secondary hover:bg-surface-tertiary hover:text-text-primary"
            title="Flag this data as problematic / unusable"
            onclick={() => { flagNoteDraft = ''; showFlagDialog = true; }}
          >
            ⚑ Flag
          </button>
        {/if}
      {:else}
        <span class="h-1.5 w-1.5 rounded-full bg-text-secondary/40"></span>
        <span class="text-[11px] text-text-secondary">Ready</span>
      {/if}
    </div>
    <span class="text-[11px] text-text-secondary/60">Rust backend</span>
  </footer>
</div>

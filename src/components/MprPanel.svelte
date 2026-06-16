<script lang="ts">
  /**
   * 田字格 (2x2 grid) MPR layout.
   *
   * Orchestrates cornerstone3D: initialises the rendering engine, creates
   * three orthographic viewports (axial, coronal, sagittal), wires up tools,
   * and reacts to volume changes in the store.
   *
   * Grid positions:
   *   [0,0] Axial    [0,1] Coronal
   *   [1,0] Sagittal [1,1] ContextPanel
   */
  import { onMount } from 'svelte';
  import { tick } from 'svelte';
  import { Enums, setVolumesForViewports } from '@cornerstonejs/core';

  import { initCornerstone, getRenderingEngine } from '$lib/cornerstone/init';
  import { setupToolGroup } from '$lib/cornerstone/tools';
  import { registerNavigate } from '$lib/navigation';
  import { volumeStore } from '$lib/stores/volumeStore.svelte';

  import { seedStore } from '$lib/stores/seedStore.svelte';
  import { pipelineStore } from '$lib/stores/pipelineStore.svelte';

  import SliceViewport from './SliceViewport.svelte';
  import ContextPanel from './ContextPanel.svelte';

  // Viewport IDs
  const VP_AXIAL = 'vp-axial';
  const VP_CORONAL = 'vp-coronal';
  const VP_SAGITTAL = 'vp-sagittal';
  const VIEWPORT_IDS = [VP_AXIAL, VP_CORONAL, VP_SAGITTAL];

  // Bindable container elements from SliceViewport children
  let axialEl = $state<HTMLDivElement | null>(null);
  let coronalEl = $state<HTMLDivElement | null>(null);
  let sagittalEl = $state<HTMLDivElement | null>(null);

  let engineReady = $state(false);

  // Check if the active vessel has a computed centerline (requires 2+ seeds)
  let hasCenterline = $derived(
    seedStore.activeVesselData?.centerline !== null &&
    seedStore.activeVesselData?.centerline !== undefined &&
    (seedStore.activeVesselData?.centerline?.length ?? 0) >= 2,
  );

  // Derive the current workflow phase for the context panel
  // Priority: analysis > seeds (centerline) > dicom > empty
  let phase = $derived<'empty' | 'dicom' | 'seeds' | 'analysis'>(
    hasCenterline
      ? 'seeds'
      : volumeStore.current
        ? 'dicom'
        : 'empty',
  );

  // ---------- Initialise cornerstone3D on mount ----------
  onMount(async () => {
    await initCornerstone();

    // Ensure DOM elements are rendered and have layout dimensions
    await tick();

    const engine = getRenderingEngine();

    // Small delay to guarantee browser has painted and elements have size
    await new Promise((r) => setTimeout(r, 50));

    if (!axialEl || !coronalEl || !sagittalEl) {
      console.error('MprPanel: viewport container elements not available');
      return;
    }

    engine.setViewports([
      {
        viewportId: VP_AXIAL,
        type: Enums.ViewportType.ORTHOGRAPHIC,
        element: axialEl,
        defaultOptions: { orientation: Enums.OrientationAxis.AXIAL },
      },
      {
        viewportId: VP_CORONAL,
        type: Enums.ViewportType.ORTHOGRAPHIC,
        element: coronalEl,
        defaultOptions: { orientation: Enums.OrientationAxis.CORONAL },
      },
      {
        viewportId: VP_SAGITTAL,
        type: Enums.ViewportType.ORTHOGRAPHIC,
        element: sagittalEl,
        defaultOptions: { orientation: Enums.OrientationAxis.SAGITTAL },
      },
    ]);

    setupToolGroup(VIEWPORT_IDS);
    engineReady = true;

    // Register the imperative navigation function for the shared service
    registerNavigate(navigateToWorldPos);

    // ---------- Notify cornerstone of viewport element size changes ----------
    // Without this, opening dev tools / dragging the window edge causes the
    // viewport DOM to shrink while cornerstone keeps using the original canvas
    // size, so canvasToWorld returns wildly wrong world coords.
    const resizeObserver = new ResizeObserver(() => {
      try {
        const e = getRenderingEngine();
        e.resize(true, false); // immediate, don't reset camera
      } catch {
        // engine may not be ready yet — safe to ignore
      }
    });
    for (const el of [axialEl, coronalEl, sagittalEl]) {
      if (el) resizeObserver.observe(el);
    }
    // Cleanup on unmount
    (window as any).__pcatMprResizeObserver?.disconnect?.();
    (window as any).__pcatMprResizeObserver = resizeObserver;
  });

  // ---------- React to volume changes ----------
  $effect(() => {
    const csVolumeId = volumeStore.cornerstoneVolumeId;
    if (!csVolumeId || !engineReady) return;

    const engine = getRenderingEngine();
    const _tBind = performance.now();

    setVolumesForViewports(
      engine,
      [{ volumeId: csVolumeId }],
      VIEWPORT_IDS,
    ).then(() => {
      // Apply slab thickness to coronal/sagittal AFTER volume is bound.
      // Use the volume's native slice spacing (sz) as the slab — this is
      // the exact resolution of slice positions, so the slab covers
      // exactly one slice's worth of data. That's enough to mask the
      // discrete-voxel snap (~0.07mm subpixel offset) without averaging
      // unrelated anatomy. For NAEOTOM 0.35mm slices the slab is 0.35mm;
      // for thicker 0.8mm CCTA recons the slab is 0.8mm. Axial is left
      // thin since it's the native acquisition direction.
      try {
        // volumeStore.spacing is [sz, sy, sx] from the Rust loader.
        const sz = volumeStore.current?.spacing?.[0] ?? 1.0;
        const slabMm = Math.max(0.35, Math.min(sz, 1.0));
        for (const id of [VP_CORONAL, VP_SAGITTAL]) {
          const vp = engine.getViewport(id) as any;
          if (vp?.setSlabThickness) {
            vp.setSlabThickness(slabMm);
          }
          if (vp?.setBlendMode && (Enums as any).BlendModes?.AVERAGE_INTENSITY_BLEND !== undefined) {
            vp.setBlendMode((Enums as any).BlendModes.AVERAGE_INTENSITY_BLEND);
          }
        }
      } catch (e) {
        console.warn('MprPanel: failed to set slab', e);
      }
      engine.renderViewports(VIEWPORT_IDS);
      console.log(
        `[load-timing] MPR setVolumesForViewports + render (cornerstone GPU upload) = ${(performance.now() - _tBind) | 0}ms`,
      );
    });
  });

  // ---------- Navigate MPR views to selected seed position ----------
  /**
   * Safe navigation: only move the camera along the view's plane normal,
   * preserving orientation. Projects the target position onto the normal
   * and shifts focal + position by that distance.
   */
  function navigateToWorldPos(pos: [number, number, number]) {
    const engine = getRenderingEngine();
    for (const vpId of VIEWPORT_IDS) {
      const vp = engine.getViewport(vpId);
      if (!vp) continue;
      const camera = vp.getCamera();
      if (!camera.focalPoint || !camera.position || !camera.viewPlaneNormal) continue;
      const focal = [...camera.focalPoint] as [number, number, number];
      const position = [...camera.position] as [number, number, number];
      const normal = camera.viewPlaneNormal as [number, number, number];

      // Project the desired position onto the view plane normal
      // to find how far to move along the normal
      const dx = pos[0] - focal[0];
      const dy = pos[1] - focal[1];
      const dz = pos[2] - focal[2];
      const dist = dx * normal[0] + dy * normal[1] + dz * normal[2];

      // Move focal and position by same amount along normal
      focal[0] += dist * normal[0];
      focal[1] += dist * normal[1];
      focal[2] += dist * normal[2];
      position[0] += dist * normal[0];
      position[1] += dist * normal[1];
      position[2] += dist * normal[2];

      vp.setCamera({ ...camera, focalPoint: focal, position });
      vp.render();
    }
  }

  // Navigation is now handled imperatively via the shared navigation service.
  // Event handlers in SliceViewport and CprView call navigateToWorldPos directly.
</script>

<!--
  2x2 grid with 1px gap.
  The parent grid has bg-border so the gap colour is the border colour,
  while each cell has its own background — producing crisp 1px grid lines.
-->
<div class="grid h-full w-full grid-cols-2 grid-rows-2 gap-px bg-border">
  <!-- [0,0] Axial -->
  <SliceViewport
    orientation="axial"
    viewportId={VP_AXIAL}
    bind:containerEl={axialEl}
  />

  <!-- [0,1] Coronal -->
  <SliceViewport
    orientation="coronal"
    viewportId={VP_CORONAL}
    bind:containerEl={coronalEl}
  />

  <!-- [1,0] Sagittal -->
  <SliceViewport
    orientation="sagittal"
    viewportId={VP_SAGITTAL}
    bind:containerEl={sagittalEl}
  />

  <!-- [1,1] Context panel -->
  <ContextPanel {phase} />
</div>

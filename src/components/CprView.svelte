<script lang="ts">
  /**
   * Compound CPR view: straightened or stretched CPR image (left 70%) with
   * 3 cross-sections (right 30%) at needle positions A, B, C.
   *
   * Two-phase architecture:
   *   Phase 1: centerline changes -> build_cpr_frame (once, ~100ms)
   *   Phase 2: rotation/needle changes -> render_cpr_image + render_cross_sections (fast, ~10ms)
   *
   * All image IPC uses raw binary (tauri::ipc::Response) -- no base64.
   *
   * Layout:
   *   +-------------------------+----------+
   *   |                         |  A (xs)  |
   *   |   Straightened/Curved   |----------|
   *   |   3 needle lines A/B/C  |  B (xs)  |
   *   |                         |----------|
   *   |                         |  C (xs)  |
   *   +-------------------------+----------+
   *            70%                  30%
   */
  import { invoke } from '@tauri-apps/api/core';
  import { navigateToWorldPos } from '$lib/navigation';
  import { seedStore, VESSEL_COLORS } from '$lib/stores/seedStore.svelte';
  import {
    type CprProjectionInfo,
    worldToStraightenedCpr,
    worldToStretchedCpr,
    straightenedCprToWorld,
    stretchedCprToWorld,
  } from '$lib/cprProjection';
  import CrossSection from './CrossSection.svelte';
  import { volumeStore } from '$lib/stores/volumeStore.svelte';
  import { pipelineStore } from '$lib/stores/pipelineStore.svelte';

  // ---- Constants ----
  const CPR_WIDTH_MM = 25.0;

  // ---- Reactive state ----
  let cprCanvas: HTMLCanvasElement | undefined = $state();
  let rotationDeg = $state(0);
  // Initialize W/L from DICOM metadata (same as MPR views)
  let windowCenter = $state(volumeStore.current?.windowCenter ?? 40);
  let windowWidth = $state(volumeStore.current?.windowWidth ?? 400);

  // CPR mode: straightened (classic) vs stretched (Horos-style projected)
  let cprMode: 'straightened' | 'stretched' = $state('stretched');

  // FAI overlay toggle
  let showFaiOverlay = $state(false);

  // Track the FAI overlay to the pipeline status: auto-enable when analysis
  // completes, and auto-hide when it's reset/cleared (e.g. centerline deleted).
  // The effect only re-runs on a status change, so a manual toggle while
  // 'complete' is preserved.
  $effect(() => {
    showFaiOverlay = pipelineStore.status === 'complete';
  });

  // Needle B position as fraction (0..1); A and C are offset
  let needleBFraction = $state(0.5);
  let needleOffset = $state(0.05); // adjustable A-C spread (fraction of total arc)

  let needleAFraction = $derived(Math.max(0, needleBFraction - needleOffset));
  let needleCFraction = $derived(Math.min(1, needleBFraction + needleOffset));

  // Image dimensions from last CPR result
  let cprWidth = $state(768);
  let cprHeight = $state(384);
  let arclengths = $state<number[]>([]);
  let cprImageData = $state<Float32Array | null>(null);
  let loading = $state(false);

  // Phase 1: frame readiness
  let frameReady = $state(false);

  // Projection info for seed overlay (fetched after each render)
  let projectionInfo = $state<CprProjectionInfo | null>(null);

  // Seed dragging state on CPR canvas
  let draggingSeedIndex = $state<number | null>(null);
  let hoverSeedIndex = $state<number | null>(null);

  // Zoom state for CPR canvas
  let cprZoom = $state(1);
  let cprPanX = $state(0);
  let cprPanY = $state(0);

  // Current centerline from seed store
  let centerline = $derived(seedStore.activeVesselData?.centerline ?? null);

  // Ostium fraction for the active vessel (used in overlays + toolbar)
  let activeOstiumFrac = $derived(seedStore.getOstiumFraction(seedStore.activeVessel));

  // Absolute arc-length (mm) of the ostium along the centerline. Displayed
  // cross-section distances are computed relative to this so that lengths
  // match the MMD view (measured from the ostium, not centerline point 0).
  let ostiumArcMm = $derived.by<number | null>(() => {
    if (activeOstiumFrac === null || arclengths.length === 0) return null;
    const n = arclengths.length;
    const idx = Math.max(0, Math.min(n - 1, Math.round(activeOstiumFrac * (n - 1))));
    return arclengths[idx];
  });

  // Batch cross-section results (one per needle: A, B, C)
  type BatchCrossSectionItem = {
    imageData: Float32Array;
    pixels: number;
    arc_mm: number;
    vesselDiameterMm: number;
    /** Lumen boundary polygon in pixel coords, one [x, y] pair per vertex. */
    vesselWall: Float32Array;
  };
  let batchXsA = $state<BatchCrossSectionItem | null>(null);
  let batchXsB = $state<BatchCrossSectionItem | null>(null);
  let batchXsC = $state<BatchCrossSectionItem | null>(null);

  // ---- Helpers ----

  /** Downsample an array of points to at most maxPts, always keeping first and last. */
  function downsample(
    pts: [number, number, number][],
    maxPts: number,
  ): [number, number, number][] {
    if (pts.length <= maxPts) return pts;
    const step = (pts.length - 1) / (maxPts - 1);
    const result: [number, number, number][] = [];
    for (let i = 0; i < maxPts - 1; i++) {
      result.push(pts[Math.round(i * step)]);
    }
    result.push(pts[pts.length - 1]);
    return result;
  }

  /**
   * Decode the raw binary CPR response:
   *   [width: u32 LE][height: u32 LE][n_arclengths: u32 LE]
   *   [arclengths: n * f64 LE]
   *   [image: width*height * f32 LE]
   */
  function decodeCprBinary(buffer: ArrayBuffer): {
    width: number;
    height: number;
    arclengths: number[];
    image: Float32Array;
  } {
    const view = new DataView(buffer);
    const width = view.getUint32(0, true);
    const height = view.getUint32(4, true);
    const nArc = view.getUint32(8, true);

    const arcOffset = 12;
    const arclengths: number[] = [];
    for (let i = 0; i < nArc; i++) {
      arclengths.push(view.getFloat64(arcOffset + i * 8, true));
    }

    const imgOffset = arcOffset + nArc * 8;
    const image = new Float32Array(buffer, imgOffset, width * height);
    return { width, height, arclengths, image };
  }

  /**
   * Decode the raw binary cross-sections response:
   *   [n_sections: u32 LE]
   *   For each:
   *     [pixels: u32 LE][arc_mm: f64 LE][diameter_mm: f64 LE][n_wall: u32 LE]
   *     [image: pixels*pixels * f32 LE]
   *     [wall: n_wall * 2 * f32 LE]
   */
  function decodeCrossSectionsBinary(buffer: ArrayBuffer): BatchCrossSectionItem[] {
    const view = new DataView(buffer);
    const nSections = view.getUint32(0, true);
    const results: BatchCrossSectionItem[] = [];
    let offset = 4;

    for (let i = 0; i < nSections; i++) {
      const pixels = view.getUint32(offset, true);
      offset += 4;
      const arc_mm = view.getFloat64(offset, true);
      offset += 8;
      const vesselDiameterMm = view.getFloat64(offset, true);
      offset += 8;
      const nWall = view.getUint32(offset, true);
      offset += 4;

      const imgLen = pixels * pixels;
      const imageData = new Float32Array(buffer, offset, imgLen);
      offset += imgLen * 4;

      const vesselWall = new Float32Array(buffer, offset, nWall * 2);
      offset += nWall * 2 * 4;

      results.push({ imageData, pixels, arc_mm, vesselDiameterMm, vesselWall });
    }
    return results;
  }

  function renderCprToCanvas(
    cvs: HTMLCanvasElement,
    data: Float32Array,
    w: number,
    h: number,
    wc: number,
    ww: number,
  ) {
    cvs.width = w;
    cvs.height = h;
    const ctx = cvs.getContext('2d')!;
    const imgData = ctx.createImageData(w, h);
    const lo = wc - ww / 2;
    const range = ww;
    for (let i = 0; i < data.length; i++) {
      const raw = data[i];
      if (raw !== raw) {
        imgData.data[i * 4] = 0;
        imgData.data[i * 4 + 1] = 0;
        imgData.data[i * 4 + 2] = 0;
        imgData.data[i * 4 + 3] = 255;
        continue;
      }
      const gray = Math.max(0, Math.min(255, Math.round(((raw - lo) / range) * 255)));

      // FAI overlay: color fat-range pixels green→red (full color, no grayscale blend)
      if (showFaiOverlay && raw >= -190 && raw <= -30) {
        const t = (raw - (-190)) / ((-30) - (-190));
        const r = Math.round(t < 0.5 ? t * 2 * 255 : 255);
        const g = Math.round(t < 0.5 ? 255 : (1 - (t - 0.5) * 2) * 255);
        imgData.data[i * 4]     = r;
        imgData.data[i * 4 + 1] = g;
        imgData.data[i * 4 + 2] = 20;
      } else {
        imgData.data[i * 4]     = gray;
        imgData.data[i * 4 + 1] = gray;
        imgData.data[i * 4 + 2] = gray;
      }
      imgData.data[i * 4 + 3] = 255;
    }
    ctx.putImageData(imgData, 0, 0);
  }

  // ---- Centerline projection cache ----
  //
  // `worldToStretchedCpr` for every centerline point costs ~0.02 ms; with
  // 100 points that's 2 ms per `drawOverlays` call, and the polyline is
  // redrawn on every overlay effect (60 Hz during needle scroll). The
  // projection is a pure function of (projectionInfo, canvas w, canvas h)
  // — cache by identity so back-to-back repaints reuse the same array.

  let centerlineProjCache: Array<[number, number] | null> | null = null;
  let centerlineProjCacheKey: {
    info: CprProjectionInfo | null;
    w: number;
    h: number;
  } = { info: null, w: 0, h: 0 };

  function getCachedCenterlineProj(
    info: CprProjectionInfo,
    w: number,
    h: number,
  ): Array<[number, number] | null> {
    if (
      centerlineProjCache !== null
      && centerlineProjCacheKey.info === info
      && centerlineProjCacheKey.w === w
      && centerlineProjCacheKey.h === h
    ) {
      return centerlineProjCache;
    }
    const nPos = info.positions.length;
    const step = Math.max(1, Math.floor(nPos / 200));
    const out: Array<[number, number] | null> = [];
    for (let j = 0; j < nPos; j += step) {
      out.push(worldToStretchedCpr(info.positions[j], info, w, h));
    }
    centerlineProjCache = out;
    centerlineProjCacheKey = { info, w, h };
    return out;
  }

  /** Draw needle lines and arclength ticks as overlay. */
  function drawOverlays(cvs: HTMLCanvasElement) {
    const ctx = cvs.getContext('2d')!;
    const w = cvs.width;
    const h = cvs.height;

    // In stretched mode, needle lines are less meaningful (vessel is not straightened)
    // but we still draw them at the same fractional position for consistency.

    // Draw needle lines (vertical)
    const drawNeedle = (fraction: number, color: string, lineWidth: number) => {
      const x = Math.round(fraction * w);
      ctx.beginPath();
      ctx.strokeStyle = color;
      ctx.lineWidth = lineWidth;
      ctx.setLineDash([]);
      ctx.moveTo(x, 0);
      ctx.lineTo(x, h);
      ctx.stroke();
    };

    if (cprMode === 'straightened') {
      // Straightened mode: vertical needle lines
      // A (yellow, dashed)
      const xA = Math.round(needleAFraction * w);
      ctx.beginPath();
      ctx.strokeStyle = '#ffee00';
      ctx.lineWidth = 1;
      ctx.setLineDash([4, 3]);
      ctx.moveTo(xA, 0);
      ctx.lineTo(xA, h);
      ctx.stroke();
      ctx.setLineDash([]);

      // B (cyan, solid, thicker)
      drawNeedle(needleBFraction, '#00ffcc', 2);

      // C (yellow, dashed)
      const xC = Math.round(needleCFraction * w);
      ctx.beginPath();
      ctx.strokeStyle = '#ffee00';
      ctx.lineWidth = 1;
      ctx.setLineDash([4, 3]);
      ctx.moveTo(xC, 0);
      ctx.lineTo(xC, h);
      ctx.stroke();
      ctx.setLineDash([]);

      // Labels
      ctx.font = 'bold 11px -apple-system, sans-serif';
      ctx.textAlign = 'center';

      ctx.fillStyle = '#ffee00';
      ctx.fillText('A', xA, 14);

      ctx.fillStyle = '#00ffcc';
      ctx.fillText('B', Math.round(needleBFraction * w), 14);

      ctx.fillStyle = '#ffee00';
      ctx.fillText('C', xC, 14);
    } else if (projectionInfo) {
      // Stretched mode: draw needle lines perpendicular to vessel tangent
      const nPos = projectionInfo.positions.length;
      const drawCurvedNeedle = (frac: number, color: string, label: string, isDashed: boolean) => {
        const idx = Math.round(frac * (nPos - 1));
        const clampedIdx = Math.min(idx, nPos - 1);
        const pos = projectionInfo!.positions[clampedIdx];
        const projected = worldToStretchedCpr(pos, projectionInfo!, w, h);
        if (!projected) return;
        const [cx, cy] = projected;

        // Compute tangent direction from adjacent projected points
        const prevIdx = Math.max(0, clampedIdx - 1);
        const nextIdx = Math.min(nPos - 1, clampedIdx + 1);
        const prevPos = projectionInfo!.positions[prevIdx];
        const nextPos = projectionInfo!.positions[nextIdx];
        const prevProj = worldToStretchedCpr(prevPos, projectionInfo!, w, h);
        const nextProj = worldToStretchedCpr(nextPos, projectionInfo!, w, h);

        if (prevProj && nextProj) {
          const tx = nextProj[0] - prevProj[0];
          const ty = nextProj[1] - prevProj[1];
          const tlen = Math.sqrt(tx * tx + ty * ty);
          if (tlen > 0.1) {
            // Perpendicular to tangent (the needle line direction)
            const px = -ty / tlen;
            const py = tx / tlen;
            const lineLen = h * 0.4; // extend across a good portion of the view

            ctx.beginPath();
            ctx.strokeStyle = color;
            ctx.lineWidth = isDashed ? 1 : 1.5;
            ctx.setLineDash(isDashed ? [4, 3] : []);
            ctx.moveTo(cx - px * lineLen, cy - py * lineLen);
            ctx.lineTo(cx + px * lineLen, cy + py * lineLen);
            ctx.stroke();
            ctx.setLineDash([]);
          }
        }

        // Label
        ctx.font = 'bold 11px -apple-system, sans-serif';
        ctx.textAlign = 'center';
        ctx.fillStyle = color;
        ctx.fillText(label, cx, cy - 10);
      };

      drawCurvedNeedle(needleAFraction, '#ffee00', 'A', true);
      drawCurvedNeedle(needleBFraction, '#00ffcc', 'B', false);
      drawCurvedNeedle(needleCFraction, '#ffee00', 'C', true);
    }

    // --- Centerline: subtle horizontal dashed line at vertical midpoint ---
    const midY = Math.round(h / 2);
    ctx.beginPath();
    ctx.strokeStyle = 'rgba(255,255,255,0.15)';
    ctx.lineWidth = 0.5;
    ctx.setLineDash([6, 4]);
    ctx.moveTo(0, midY);
    ctx.lineTo(w, midY);
    ctx.stroke();
    ctx.setLineDash([]);

    // Arclength ticks every 10mm along the bottom (straightened mode only)
    if (cprMode === 'straightened' && arclengths.length > 0) {
      const totalArc = arclengths[arclengths.length - 1];
      ctx.font = '9px -apple-system, sans-serif';
      ctx.fillStyle = '#98989d';
      ctx.textAlign = 'center';
      ctx.strokeStyle = '#98989d';
      ctx.lineWidth = 0.5;

      for (let mm = 0; mm <= totalArc; mm += 10) {
        const frac = mm / totalArc;
        const x = Math.round(frac * w);
        // tick mark
        ctx.beginPath();
        ctx.moveTo(x, h - 10);
        ctx.lineTo(x, h - 2);
        ctx.stroke();
        // label
        ctx.fillText(`${mm}`, x, h - 12);
      }
    }

    // --- Ostium marker ---
    if (activeOstiumFrac !== null && cprMode === 'straightened') {
      const ox = Math.round(activeOstiumFrac * w);

      // Shaded region: proximal side (before ostium) is dimmed
      ctx.fillStyle = 'rgba(0,0,0,0.25)';
      ctx.fillRect(0, 0, ox, h);

      // Solid magenta vertical line
      ctx.beginPath();
      ctx.strokeStyle = '#ff00ff';
      ctx.lineWidth = 2;
      ctx.setLineDash([]);
      ctx.moveTo(ox, 0);
      ctx.lineTo(ox, h);
      ctx.stroke();

      // Label with background
      ctx.font = 'bold 11px -apple-system, sans-serif';
      ctx.textAlign = 'center';
      ctx.fillStyle = '#ff00ff';
      ctx.fillText('OSTIUM', ox, h - 6);
    } else if (activeOstiumFrac !== null && cprMode === 'stretched' && projectionInfo) {
      // In stretched mode, draw ostium marker at the projected position
      const nPos = projectionInfo.positions.length;
      const idx = Math.round(activeOstiumFrac * (nPos - 1));
      const clampedIdx = Math.min(idx, nPos - 1);
      const pos = projectionInfo.positions[clampedIdx];
      const projected = worldToStretchedCpr(pos, projectionInfo, w, h);
      if (projected) {
        const [cx, cy] = projected;

        // Diamond marker
        ctx.save();
        ctx.translate(cx, cy);
        ctx.rotate(Math.PI / 4);
        ctx.fillStyle = 'rgba(255,0,255,0.7)';
        ctx.fillRect(-5, -5, 10, 10);
        ctx.strokeStyle = 'white';
        ctx.lineWidth = 1.5;
        ctx.strokeRect(-5, -5, 10, 10);
        ctx.restore();

        // Label
        ctx.font = 'bold 11px -apple-system, sans-serif';
        ctx.textAlign = 'center';
        ctx.fillStyle = '#ff00ff';
        ctx.fillText('OSTIUM', cx, cy - 12);
      }
    }

    // --- Centerline polyline on CPR ---
    // The polyline is redrawn on every overlay repaint (which fires 60 Hz
    // during a needle scroll / drag). Recomputing 100+ `worldToStretchedCpr`
    // — each of which is an O(n_cols=128) nearest-segment search — burns
    // ~2 ms per paint. Cache the projected 2-D points and invalidate only
    // when the inputs that affect projection actually change.
    if (projectionInfo && cprMode === 'stretched') {
      const color = VESSEL_COLORS[seedStore.activeVessel];
      const projPts = getCachedCenterlineProj(projectionInfo, w, h);

      ctx.beginPath();
      ctx.strokeStyle = color;
      ctx.lineWidth = 1.5;
      ctx.globalAlpha = 0.5;
      let started = false;
      for (const p of projPts) {
        if (!p) { started = false; continue; }
        if (!started) { ctx.moveTo(p[0], p[1]); started = true; }
        else { ctx.lineTo(p[0], p[1]); }
      }
      ctx.stroke();
      ctx.globalAlpha = 1.0;
    }

    // --- Seed markers on CPR ---
    if (projectionInfo) {
      const vessel = seedStore.activeVessel;
      const data = seedStore.activeVesselData;
      const color = VESSEL_COLORS[vessel];
      const selectedIdx = seedStore.selectedSeedIndex;

      for (let i = 0; i < data.seeds.length; i++) {
        const seedPos = data.seeds[i].position;
        // Seeds are in [x,y,z] world coords; projection expects [z,y,x]
        const seedZyx: [number, number, number] = [seedPos[2], seedPos[1], seedPos[0]];

        const projected = cprMode === 'stretched'
          ? worldToStretchedCpr(seedZyx, projectionInfo, w, h)
          : worldToStraightenedCpr(seedZyx, projectionInfo, w, h);

        if (!projected) continue;
        const [cx, cy] = projected;

        const isSelected = i === selectedIdx;
        const isHovered = i === hoverSeedIndex;
        const radius = isSelected || isHovered ? 5 : 4;

        // Glow ring for selected/hovered
        if (isSelected) {
          ctx.beginPath();
          ctx.arc(cx, cy, radius + 3, 0, Math.PI * 2);
          ctx.strokeStyle = 'rgba(255,255,255,0.5)';
          ctx.lineWidth = 2;
          ctx.stroke();
        } else if (isHovered) {
          ctx.beginPath();
          ctx.arc(cx, cy, radius + 2, 0, Math.PI * 2);
          ctx.strokeStyle = 'rgba(255,255,255,0.3)';
          ctx.lineWidth = 1.5;
          ctx.stroke();
        }

        // Seed circle
        ctx.beginPath();
        ctx.arc(cx, cy, radius, 0, Math.PI * 2);
        ctx.fillStyle = color;
        ctx.fill();
        ctx.strokeStyle = 'white';
        ctx.lineWidth = 1.2;
        ctx.stroke();

        // Index label
        ctx.font = 'bold 11px -apple-system, sans-serif';
        ctx.textAlign = 'center';
        ctx.fillStyle = 'white';
        ctx.fillText(`${i}`, cx, cy - radius - 3);
      }
    }

    // Mode badge removed — toolbar already shows the mode.
  }

  // ---- Two-canvas compositor ----
  //
  // The old single-canvas repaint ran a 260k–520k pixel RGBA loop on every
  // mousemove (hover / seed drag / needle position change) because the
  // `$effect` that painted the image also depended on overlay state. Moving
  // the slider felt like ~30-60 ms per frame on the main thread.
  //
  // Now the HU→gray loop writes to a hidden buffer canvas once per image
  // change (W/L, FAI, new render). All overlay-only updates just
  // `drawImage(buffer)` + redraw overlays, which is ~1 ms regardless of
  // image size.

  /** Hidden buffer canvas that holds the rendered CPR image pixels. */
  let bufferCanvas: HTMLCanvasElement | null = null;
  /** True once `bufferCanvas` has image content matching current cprImageData. */
  let bufferReady = false;

  function ensureBuffer(w: number, h: number): HTMLCanvasElement {
    if (!bufferCanvas) bufferCanvas = document.createElement('canvas');
    if (bufferCanvas.width !== w) bufferCanvas.width = w;
    if (bufferCanvas.height !== h) bufferCanvas.height = h;
    return bufferCanvas;
  }

  /** Paint raw HU data through the window/level + FAI transform into the buffer. */
  function paintImageToBuffer() {
    if (!cprImageData) {
      bufferReady = false;
      return;
    }
    const buf = ensureBuffer(cprWidth, cprHeight);
    renderCprToCanvas(buf, cprImageData, cprWidth, cprHeight, windowCenter, windowWidth);
    bufferReady = true;
  }

  /** Fast path: blit the buffer to the visible canvas, then draw overlays. */
  function compositeToMain() {
    if (!cprCanvas) return;
    if (!bufferReady || !bufferCanvas) return;
    // Keep the visible canvas's pixel dims in sync with the buffer so the
    // blit is 1:1 (no sub-pixel blur, no implicit resampling).
    if (cprCanvas.width !== bufferCanvas.width) cprCanvas.width = bufferCanvas.width;
    if (cprCanvas.height !== bufferCanvas.height) cprCanvas.height = bufferCanvas.height;
    const ctx = cprCanvas.getContext('2d')!;
    ctx.drawImage(bufferCanvas, 0, 0);
    drawOverlays(cprCanvas);
  }

  // ---- Phase 1: Build frame when centerline changes ----

  let frameDebounce: ReturnType<typeof setTimeout> | undefined;
  let buildingFrame = false;

  $effect(() => {
    const cl = centerline;

    if (!cl || cl.length < 2) {
      cprImageData = null;
      frameReady = false;
      return;
    }

    // Centerline changed -- rebuild frame
    frameReady = false;
    loading = true;
    clearTimeout(frameDebounce);
    frameDebounce = setTimeout(async () => {
      if (buildingFrame) return;
      buildingFrame = true;
      try {
        // Downsample to max 100 points for smooth CPR while keeping IPC fast
        const sampled = downsample(cl, 100);
        // The spline centerline is in cornerstone3D world coords [x, y, z]
        // Rust expects [z, y, x]
        const centerlineZyx = sampled.map(
          ([x, y, z]) => [z, y, x] as [number, number, number],
        );

        await invoke('build_cpr_frame', {
          centerlineMm: centerlineZyx,
          pixelsWide: 768,
        });

        frameReady = true;
        // Trigger initial render
        renderCpr();
        renderCrossSections();
      } catch (e) {
        console.error('CprView: build_cpr_frame failed', e);
      } finally {
        buildingFrame = false;
        loading = false;
      }
    }, 200);

    return () => clearTimeout(frameDebounce);
  });

  // ---- Phase 2: Render when rotation or mode changes (uses cached frame) ----

  let renderDebounce: ReturnType<typeof setTimeout> | undefined;
  let renderingCpr = false;

  $effect(() => {
    // Track rotation and mode deps
    void rotationDeg;
    void cprMode;

    if (!frameReady) return;

    clearTimeout(renderDebounce);
    renderDebounce = setTimeout(() => {
      renderCpr();
      renderCrossSections();
    }, 50);

    return () => clearTimeout(renderDebounce);
  });

  async function renderCpr() {
    if (!frameReady || renderingCpr) return;
    renderingCpr = true;
    try {
      let buffer: ArrayBuffer;

      if (cprMode === 'stretched') {
        buffer = await invoke<ArrayBuffer>('render_stretched_cpr_image', {
          rotationDeg,
          widthMm: CPR_WIDTH_MM,
          pixelsWide: 512,
          pixelsHigh: 512,
          slabMm: 1.0,
        });
      } else {
        buffer = await invoke<ArrayBuffer>('render_cpr_image', {
          rotationDeg,
          widthMm: CPR_WIDTH_MM,
          pixelsHigh: 384,
          slabMm: 1.0,
        });
      }

      const decoded = decodeCprBinary(buffer);
      cprWidth = decoded.width;
      cprHeight = decoded.height;
      arclengths = decoded.arclengths;
      cprImageData = decoded.image;

      // Fetch projection info for seed overlay (lightweight JSON)
      fetchProjectionInfo();
    } catch (e) {
      console.error('CprView: render CPR failed', e);
    } finally {
      renderingCpr = false;
    }
  }

  /** Fetch projection info from Rust for seed overlay mapping. */
  async function fetchProjectionInfo() {
    if (!frameReady) {
      projectionInfo = null;
      return;
    }
    try {
      projectionInfo = await invoke<CprProjectionInfo>('get_cpr_projection_info', {
        rotationDeg,
        widthMm: CPR_WIDTH_MM,
        pixelsWide: 512,
        pixelsHigh: 512,
      });
    } catch (e) {
      console.error('CprView: get_cpr_projection_info failed', e);
      projectionInfo = null;
    }
  }

  // Image-level changes: raw HU data, W/L window, FAI colour overlay — these
  // all require re-running the per-pixel transform, so repaint the buffer
  // and composite.
  $effect(() => {
    void cprImageData;
    void windowCenter;
    void windowWidth;
    void showFaiOverlay;

    paintImageToBuffer();
    compositeToMain();
  });

  // Overlay-only changes: seed dots, needle lines, arclength ticks, hover
  // highlight, ostium marker, projection mode. None of these touch pixel
  // values, so we skip `paintImageToBuffer` entirely and just blit + draw
  // overlays (~1 ms vs the old ~30–60 ms).
  $effect(() => {
    void needleAFraction;
    void needleBFraction;
    void needleCFraction;
    void arclengths;
    void cprMode;
    void seedStore.activeVesselData;
    void seedStore.selectedSeedIndex;
    void activeOstiumFrac;
    void projectionInfo;
    void hoverSeedIndex;

    compositeToMain();
  });

  // ---- Cross-section computation (uses cached frame, raw binary) ----

  // Cross-section updates were on a 100 ms debounce — that means the three
  // panels sit frozen during a scroll and only update once the user stops.
  // Feels sluggish. Switch to leading-edge + in-flight-pipeline: fire
  // immediately, and if more needle/rotation changes arrive while the IPC
  // is still in flight, remember the latest one and fire once when the
  // previous completes. This gives real-time XS updates at whatever rate
  // the IPC round-trip allows (typically ~15–25 ms in release), without
  // queuing up stale requests.
  let xsPending = false;

  $effect(() => {
    // Track needle and rotation deps
    void needleAFraction;
    void needleBFraction;
    void needleCFraction;
    void rotationDeg;

    if (!frameReady) {
      batchXsA = null;
      batchXsB = null;
      batchXsC = null;
      return;
    }

    kickCrossSections();
  });

  function kickCrossSections() {
    if (!frameReady) return;
    if (computingXs) {
      // An IPC is already in flight; just note that inputs changed. The
      // in-flight call's `finally` will re-fire with the newest values.
      xsPending = true;
      return;
    }
    void renderCrossSections();
  }

  let computingXs = false;
  async function renderCrossSections() {
    if (!frameReady || computingXs) return;
    computingXs = true;
    xsPending = false;
    try {
      const buffer = await invoke<ArrayBuffer>('render_cross_sections', {
        positionFractions: [needleAFraction, needleBFraction, needleCFraction],
        rotationDeg,
        widthMm: 15.0,
        pixels: 128,
      });

      const results = decodeCrossSectionsBinary(buffer);
      if (results.length === 3) {
        batchXsA = results[0];
        batchXsB = results[1];
        batchXsC = results[2];
      }
    } catch (e) {
      console.error('CprView: render_cross_sections failed', e);
    } finally {
      computingXs = false;
      // If any needle/rotation change arrived while we were in flight,
      // re-fire now with the latest values. This is the "trailing" edge
      // of the leading-edge pipeline and keeps the XS panels fresh at
      // whatever rate the IPC round-trip allows.
      if (xsPending) {
        xsPending = false;
        void renderCrossSections();
      }
    }
  }

  // ---- Seed selection -> needle B navigation ----

  $effect(() => {
    const selectedIdx = seedStore.selectedSeedIndex;
    if (selectedIdx === null) return;
    const data = seedStore.activeVesselData;
    if (!data || selectedIdx >= data.seeds.length || !data.centerline) return;

    const seedPos = data.seeds[selectedIdx].position;
    const cl = data.centerline;
    if (cl.length < 2) return;

    // Find closest centerline point to the selected seed
    let minDist = Infinity;
    let closestIdx = 0;
    for (let i = 0; i < cl.length; i++) {
      const dx = cl[i][0] - seedPos[0];
      const dy = cl[i][1] - seedPos[1];
      const dz = cl[i][2] - seedPos[2];
      const d = dx * dx + dy * dy + dz * dz;
      if (d < minDist) {
        minDist = d;
        closestIdx = i;
      }
    }

    needleBFraction = closestIdx / (cl.length - 1);
  });

  // ---- Helpers: navigate MPR to needle B's current position ----

  /** Compute the world position at needle B and navigate MPR views there. */
  function navigateToNeedlePos() {
    const cl = seedStore.activeVesselData?.centerline;
    if (!cl || cl.length < 2) return;
    const idx = Math.round(needleBFraction * (cl.length - 1));
    const pos = cl[Math.min(idx, cl.length - 1)];
    navigateToWorldPos(pos);
  }

  // Mousemove / wheel fire *much* faster than cornerstone3D can re-render
  // the three MPR viewports. Each `navigateToWorldPos` call does
  // setCamera + render on all three, which on a CCTA volume is ~15–30 ms
  // of real WebGL work per call. Raw 60 Hz would ask for 1200 ms/sec of
  // MPR work, which is obviously impossible — the main thread stalls and
  // the UI drops frames.
  //
  // Throttle-with-leading-and-trailing: first scroll event fires
  // immediately so the MPR viewports start tracking; subsequent events
  // within the throttle interval are coalesced into a single trailing
  // call scheduled at the next interval boundary. This way the MPR
  // updates at ~7 fps during active scroll (continuous feedback, not
  // frozen), and always snaps to the final position on release.
  //
  // 140 ms ≈ 7 fps. Each navigate costs ~20 ms of WebGL, so this keeps
  // MPR overhead at ~15 % of the main thread while the user is actively
  // scrolling — enough headroom for the CPR overlay + cross-section
  // pipeline to stay smooth.
  const NAV_INTERVAL_MS = 140;
  let lastNavTime = 0;
  let navigateTimer: ReturnType<typeof setTimeout> | null = null;
  function scheduleNavigate() {
    const now = performance.now();
    const elapsed = now - lastNavTime;
    if (navigateTimer !== null) {
      clearTimeout(navigateTimer);
      navigateTimer = null;
    }
    if (elapsed >= NAV_INTERVAL_MS) {
      // Leading edge: fire immediately so MPR starts tracking as soon as
      // the user begins to scroll / drag.
      lastNavTime = now;
      navigateToNeedlePos();
    } else {
      // Trailing edge: schedule at the next interval boundary. This also
      // covers the "scroll ended; fire one final call with the true
      // final position" case.
      navigateTimer = setTimeout(() => {
        navigateTimer = null;
        lastNavTime = performance.now();
        navigateToNeedlePos();
      }, NAV_INTERVAL_MS - elapsed);
    }
  }

  // Precision touchpads fire wheel events at ~240 Hz. If we mutate
  // `needleBFraction` on every one, the overlay `$effect` fires 240×/sec
  // and `compositeToMain` (blit + draw overlays ≈ 5 ms) pegs the main
  // thread. Coalesce state updates to animation-frame rate: wheel/drag
  // handlers push into a non-reactive scratch variable, and we commit
  // once per frame. Downstream effects then fire at most 60 Hz.
  let pendingNeedleB: number | null = null;
  let needleBFrame: number | null = null;
  function commitNeedleB() {
    needleBFrame = null;
    if (pendingNeedleB === null) return;
    const frac = Math.max(0, Math.min(1, pendingNeedleB));
    pendingNeedleB = null;
    if (frac !== needleBFraction) {
      needleBFraction = frac;
      scheduleNavigate();
    }
  }
  function setNeedleBThrottled(fraction: number) {
    pendingNeedleB = fraction;
    if (needleBFrame !== null) return;
    needleBFrame = requestAnimationFrame(commitNeedleB);
  }

  // ---- Needle dragging ----

  let dragging = $state(false);

  /** Find the seed index nearest to the given canvas position, within hitRadius. */
  function findSeedAtCanvasPos(canvasX: number, canvasY: number, hitRadius: number): number | null {
    if (!projectionInfo || !cprCanvas) return null;
    const data = seedStore.activeVesselData;
    const w = cprCanvas.width;
    const h = cprCanvas.height;
    const rect = cprCanvas.getBoundingClientRect();
    const scaleX = rect.width / w;
    const scaleY = rect.height / h;

    let bestIdx: number | null = null;
    let bestDist = hitRadius * hitRadius;

    for (let i = 0; i < data.seeds.length; i++) {
      const seedPos = data.seeds[i].position;
      const seedZyx: [number, number, number] = [seedPos[2], seedPos[1], seedPos[0]];
      const projected = cprMode === 'stretched'
        ? worldToStretchedCpr(seedZyx, projectionInfo, w, h)
        : worldToStraightenedCpr(seedZyx, projectionInfo, w, h);
      if (!projected) continue;

      const dx = canvasX - projected[0] * scaleX;
      const dy = canvasY - projected[1] * scaleY;
      const d = dx * dx + dy * dy;
      if (d < bestDist) {
        bestDist = d;
        bestIdx = i;
      }
    }
    return bestIdx;
  }

  function onCanvasMouseDown(event: MouseEvent) {
    if (event.button !== 0 || !cprCanvas) return;

    // Shift+click sets ostium marker
    if (event.shiftKey) {
      const rect = cprCanvas.getBoundingClientRect();
      const x = event.clientX - rect.left;
      const fraction = Math.max(0, Math.min(1, x / rect.width));
      seedStore.setOstiumFraction(fraction);
      event.preventDefault();
      return;
    }

    const rect = cprCanvas.getBoundingClientRect();
    const x = event.clientX - rect.left;
    const y = event.clientY - rect.top;

    // Priority 2: Check if near a seed (8px hit radius)
    const seedIdx = findSeedAtCanvasPos(x, y, 8);
    if (seedIdx !== null) {
      seedStore.selectSeed(seedIdx);
      draggingSeedIndex = seedIdx;
      event.preventDefault();
      return;
    }

    if (cprMode === 'straightened') {
      const fraction = x / rect.width;
      // Check if click is near needle B (within ~8px)
      const bPixel = needleBFraction * rect.width;
      if (Math.abs(x - bPixel) < 8) {
        dragging = true;
        event.preventDefault();
      } else {
        needleBFraction = Math.max(0, Math.min(1, fraction));
        navigateToNeedlePos();
      }
    } else if (projectionInfo && cprCanvas) {
      // Stretched mode: find nearest centerline point to click position
      const canvasPixelX = (x / rect.width) * cprCanvas.width;
      const canvasPixelY = (y / rect.height) * cprCanvas.height;
      const nPos = projectionInfo.positions.length;
      const w = cprCanvas.width;
      const h = cprCanvas.height;

      let bestIdx = 0;
      let bestDist = Infinity;
      for (let j = 0; j < nPos; j++) {
        const projected = worldToStretchedCpr(projectionInfo.positions[j], projectionInfo, w, h);
        if (!projected) continue;
        const dx = canvasPixelX - projected[0];
        const dy = canvasPixelY - projected[1];
        const d = dx * dx + dy * dy;
        if (d < bestDist) {
          bestDist = d;
          bestIdx = j;
        }
      }
      needleBFraction = bestIdx / (nPos - 1);
      navigateToNeedlePos();
    }
  }

  function onCanvasMouseMove(event: MouseEvent) {
    if (!cprCanvas) return;
    const rect = cprCanvas.getBoundingClientRect();
    const x = event.clientX - rect.left;
    const y = event.clientY - rect.top;

    // Seed dragging takes priority
    if (draggingSeedIndex !== null && projectionInfo) {
      // Map canvas CSS position to canvas pixel position
      const canvasPixelX = (x / rect.width) * cprCanvas.width;
      const canvasPixelY = (y / rect.height) * cprCanvas.height;

      // Unproject to 3D
      const worldZyx = cprMode === 'stretched'
        ? stretchedCprToWorld(canvasPixelX, canvasPixelY, projectionInfo, cprCanvas.width, cprCanvas.height)
        : straightenedCprToWorld(canvasPixelX, canvasPixelY, projectionInfo, cprCanvas.width, cprCanvas.height);

      // Convert [z,y,x] back to [x,y,z] for seedStore
      const worldXyz: [number, number, number] = [worldZyx[2], worldZyx[1], worldZyx[0]];
      seedStore.moveSeed(draggingSeedIndex, worldXyz);
      return;
    }

    // Needle B dragging
    if (dragging) {
      if (cprMode === 'straightened') {
        setNeedleBThrottled(x / rect.width);
      }
      // In stretched mode, needle dragging isn't supported (use scroll instead)
      return;
    }

    // Hover detection for seed highlighting
    const seedIdx = findSeedAtCanvasPos(x, y, 8);
    if (seedIdx !== hoverSeedIndex) {
      hoverSeedIndex = seedIdx;
    }
  }

  function onCanvasMouseUp() {
    if (draggingSeedIndex !== null) {
      draggingSeedIndex = null;
      return;
    }
    if (dragging) {
      navigateToNeedlePos();
    }
    dragging = false;
  }

  /** Scroll/pinch on CPR canvas.
   *  Pinch (ctrlKey): zoom toward cursor.
   *  Scroll: move needle B position. */
  function onCanvasWheel(event: WheelEvent) {
    event.preventDefault();

    if (event.ctrlKey) {
      // Pinch-to-zoom toward cursor using transform-origin
      const zoomFactor = 1 - event.deltaY * 0.01;
      const newZoom = Math.max(1, Math.min(8, cprZoom * zoomFactor));
      if (cprCanvas) {
        const rect = cprCanvas.parentElement!.getBoundingClientRect();
        const cursorX = (event.clientX - rect.left) / rect.width;
        const cursorY = (event.clientY - rect.top) / rect.height;
        // Set transform-origin to cursor position (as offset from center)
        cprPanX = 0.5 - cursorX;
        cprPanY = 0.5 - cursorY;
      }
      cprZoom = newZoom;
    } else {
      // Scroll: move needle B. Accumulate into the pending value (not
      // `needleBFraction` directly) so the reactive graph only fires once
      // per animation frame, no matter how many wheel events arrive.
      const sensitivity = 0.0003;
      const delta = -event.deltaY * sensitivity;
      const base = pendingNeedleB ?? needleBFraction;
      setNeedleBThrottled(base + delta);
    }
  }

  // ---- W/L adjustment via right-drag on CPR canvas ----

  let wlDragging = $state(false);
  let wlStartX = $state(0);
  let wlStartY = $state(0);
  let wlStartCenter = $state(40);
  let wlStartWidth = $state(400);

  function onCanvasContextMenu(event: MouseEvent) {
    event.preventDefault();
  }

  function onCanvasRightDown(event: MouseEvent) {
    if (event.button !== 2) return;
    wlDragging = true;
    wlStartX = event.clientX;
    wlStartY = event.clientY;
    wlStartCenter = windowCenter;
    wlStartWidth = windowWidth;
    event.preventDefault();
  }

  function onCanvasRightMove(event: MouseEvent) {
    if (!wlDragging) return;
    const dx = event.clientX - wlStartX;
    const dy = event.clientY - wlStartY;
    windowWidth = Math.max(1, wlStartWidth + dx * 2);
    windowCenter = wlStartCenter - dy;
  }

  function onCanvasRightUp() {
    wlDragging = false;
  }

  // Combined mouse handlers
  function handleMouseDown(event: MouseEvent) {
    if (event.button === 0) onCanvasMouseDown(event);
    else if (event.button === 2) onCanvasRightDown(event);
  }

  function handleMouseMove(event: MouseEvent) {
    if (draggingSeedIndex !== null || dragging) onCanvasMouseMove(event);
    else if (wlDragging) onCanvasRightMove(event);
    else onCanvasMouseMove(event); // hover detection
  }

  function handleMouseUp(event: MouseEvent) {
    if (event.button === 0) onCanvasMouseUp();
    else if (event.button === 2) onCanvasRightUp();
  }
</script>

<!-- Attach global listeners for drag continuation outside canvas -->
<svelte:window
  onmousemove={handleMouseMove}
  onmouseup={handleMouseUp}
/>

<div class="flex h-full w-full flex-col bg-surface-secondary">
  <!-- Main layout: CPR + cross-sections -->
  <div class="flex min-h-0 flex-1">
    <!-- Left: CPR (70%) -->
    <div class="relative flex min-h-0 flex-[7] flex-col">
      <!-- CPR label -->
      <span
        class="pointer-events-none absolute left-2 top-1.5 z-10 text-[11px] font-semibold tracking-wider text-text-secondary/60"
      >
        CPR &mdash; {seedStore.activeVessel}
      </span>

      <!-- W/L display -->
      <span
        class="pointer-events-none absolute right-2 top-1.5 z-10 text-[10px] font-mono tabular-nums text-text-secondary/50"
      >
        C:{windowCenter} W:{windowWidth}
      </span>

      <!-- Loading indicator -->
      {#if loading}
        <div class="absolute inset-0 z-20 flex items-center justify-center bg-black/30">
          <span class="text-xs text-text-secondary">Computing CPR...</span>
        </div>
      {/if}

      {#if !centerline || centerline.length < 2}
        <div class="flex flex-1 items-center justify-center">
          <p class="text-xs text-text-secondary/60">
            Place seeds along the vessel (start in aorta, trace into coronary)
          </p>
        </div>
      {:else}
        <!-- svelte-ignore a11y_no_static_element_interactions -->
        <div
          class="min-h-0 flex-1 overflow-hidden relative"
          style:cursor={draggingSeedIndex !== null ? 'grabbing' : hoverSeedIndex !== null ? 'grab' : 'crosshair'}
          onmousedown={handleMouseDown}
          oncontextmenu={onCanvasContextMenu}
          onwheel={onCanvasWheel}
        >
          <canvas
            bind:this={cprCanvas}
            class="absolute inset-0"
            style="image-rendering: pixelated; width: 100%; height: 100%; transform: scale({cprZoom}); transform-origin: {(0.5 - cprPanX) * 100}% {(0.5 - cprPanY) * 100}%;"
          ></canvas>
        </div>
      {/if}
    </div>

    <!-- Right: 3 cross-sections (30%) -->
    <div
      class="flex min-h-0 flex-[3] flex-col gap-px border-l border-border bg-border"
    >
      {#if centerline && centerline.length >= 2}
        <div class="flex min-h-0 flex-1 flex-col bg-black">
          <CrossSection
            centerlineMm={centerline}
            positionFraction={needleAFraction}
            {rotationDeg}
            label="A"
            color="#ffee00"
            {windowCenter}
            {windowWidth}
            batchImageData={batchXsA?.imageData ?? null}
            batchPixels={batchXsA?.pixels ?? null}
            arcMmProp={batchXsA?.arc_mm ?? null}
            vesselDiameterMm={batchXsA?.vesselDiameterMm ?? null}
            vesselWall={batchXsA?.vesselWall ?? null}
            showFaiOverlay={showFaiOverlay}
            arcOffsetMm={ostiumArcMm}
          />
        </div>
        <div class="flex min-h-0 flex-1 flex-col bg-black">
          <CrossSection
            centerlineMm={centerline}
            positionFraction={needleBFraction}
            {rotationDeg}
            label="B"
            color="#00ffcc"
            {windowCenter}
            {windowWidth}
            batchImageData={batchXsB?.imageData ?? null}
            batchPixels={batchXsB?.pixels ?? null}
            arcMmProp={batchXsB?.arc_mm ?? null}
            vesselDiameterMm={batchXsB?.vesselDiameterMm ?? null}
            vesselWall={batchXsB?.vesselWall ?? null}
            showFaiOverlay={showFaiOverlay}
            arcOffsetMm={ostiumArcMm}
          />
        </div>
        <div class="flex min-h-0 flex-1 flex-col bg-black">
          <CrossSection
            centerlineMm={centerline}
            positionFraction={needleCFraction}
            {rotationDeg}
            label="C"
            color="#ffee00"
            {windowCenter}
            {windowWidth}
            batchImageData={batchXsC?.imageData ?? null}
            batchPixels={batchXsC?.pixels ?? null}
            arcMmProp={batchXsC?.arc_mm ?? null}
            vesselDiameterMm={batchXsC?.vesselDiameterMm ?? null}
            vesselWall={batchXsC?.vesselWall ?? null}
            showFaiOverlay={showFaiOverlay}
            arcOffsetMm={ostiumArcMm}
          />
        </div>
      {:else}
        <div class="flex flex-1 items-center justify-center bg-surface-secondary">
          <p class="text-[10px] text-text-secondary/40">No cross-sections</p>
        </div>
      {/if}
    </div>
  </div>

  <!-- Bottom toolbar: rotation slider + mode toggle -->
  <div
    class="flex shrink-0 flex-wrap items-center gap-x-3 gap-y-1 border-t border-border bg-surface-secondary px-3 py-1"
  >
    <label class="flex items-center gap-2 text-[10px] text-text-secondary">
      <span>Rot</span>
      <input
        type="range"
        min="0"
        max="360"
        step="1"
        bind:value={rotationDeg}
        class="h-1 w-24 cursor-pointer accent-accent"
        title="Rotate the cross-section viewing angle"
      />
      <span class="w-7 text-right tabular-nums">{rotationDeg}&deg;</span>
    </label>

    <span class="text-[10px] text-text-secondary/40">|</span>

    <!-- Stretched / Straightened toggle -->
    <button
      class="rounded px-1.5 py-0.5 text-[10px] font-medium transition-colors
        {cprMode === 'stretched'
          ? 'bg-accent/20 text-accent'
          : 'text-text-secondary/60 hover:text-text-secondary'}"
      onclick={() => { cprMode = 'stretched'; }}
    >
      Stretched
    </button>
    <button
      class="rounded px-1.5 py-0.5 text-[10px] font-medium transition-colors
        {cprMode === 'straightened'
          ? 'bg-accent/20 text-accent'
          : 'text-text-secondary/60 hover:text-text-secondary'}"
      onclick={() => { cprMode = 'straightened'; }}
    >
      Straightened
    </button>

    <button
      class="rounded px-1.5 py-0.5 text-[10px] font-medium transition-colors
        {showFaiOverlay
          ? 'bg-error/20 text-error'
          : 'text-text-secondary/60 hover:text-text-secondary'}"
      onclick={() => { showFaiOverlay = !showFaiOverlay; }}
      title="Fat Attenuation Index overlay: green = healthy fat, red = inflamed fat"
    >
      FAI
    </button>

    <span class="text-[10px] text-text-secondary/40">|</span>

    <span class="text-[10px] tabular-nums text-text-secondary/50">
      B: {(needleBFraction * 100).toFixed(0)}%
    </span>

    <label class="flex items-center gap-1 text-[10px] text-text-secondary">
      <span>Spread</span>
      <input
        type="range"
        min="0.01"
        max="0.20"
        step="0.01"
        bind:value={needleOffset}
        class="h-1 w-16 cursor-pointer accent-accent"
        title="Spacing between cross-section positions A and C"
      />
      <span class="w-6 text-right tabular-nums">{(needleOffset * 100).toFixed(0)}%</span>
    </label>

    <span class="text-[10px] text-text-secondary/40">|</span>

    <button
      class="rounded px-2 py-0.5 text-[10px] font-medium transition-colors"
      style={activeOstiumFrac !== null
        ? 'background-color: rgba(255,0,255,0.15); color: #ff00ff;'
        : 'color: #ff00ff;'}
      onclick={() => seedStore.setOstiumFraction(needleBFraction)}
      title="Mark where the coronary artery exits the aorta (proximal reference point for analysis)"
    >
      {activeOstiumFrac !== null ? 'Ostium Set' : 'Set Ostium'}
    </button>

    {#if activeOstiumFrac !== null}
      <span class="text-[10px] tabular-nums" style="color: #ff00ff;">
        {(activeOstiumFrac * 100).toFixed(0)}%
      </span>
      <button
        class="text-[10px] text-text-secondary/40 hover:text-error"
        onclick={() => seedStore.setOstiumFraction(null)}
        title="Clear ostium"
      >
        &times;
      </button>
    {/if}

    <span class="text-[10px] text-text-secondary/30">
      Scroll: needle &middot; Shift+click: ostium &middot; Drag: refine &middot; Right-drag: W/L
    </span>
  </div>
</div>

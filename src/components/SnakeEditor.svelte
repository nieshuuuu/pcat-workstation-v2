<script lang="ts">
  /**
   * Canvas overlay for editing active contour (snake) annotations on a
   * cross-section image.
   *
   * Renders the HU image as grayscale background, vessel wall as red dashed
   * polygon, init boundary as blue dashed polygon, and the snake contour as
   * a green polygon with draggable control points.
   */
  import type { AnnotationTarget } from '$lib/api';
  import {
    updateSnakePoints,
    addSnakePoint,
    useVesselWallAsContour,
  } from '$lib/api';

  type Props = {
    target: AnnotationTarget;
    targetIndex: number;
    /** Current snake contour points [x,y] in pixel coords */
    snakePoints: [number, number][] | null;
    onSnakeUpdate: (points: [number, number][]) => void;
    status: 'pending' | 'in-progress' | 'done';
    /** Absolute arc-length (mm) of the ostium along the centerline.
     *  Displayed arc = target.arc_mm - arcOffsetMm. */
    arcOffsetMm?: number;
    /** Optional material-decomposition overlay (flat pixels×pixels, values in
     *  volume fraction [0,1] or mass density mg/mL). Rendered as a rainbow
     *  colormap on top of the HU grayscale when non-null. NaN = skip voxel. */
    overlay?: number[] | null;
    /** Material label shown in the colorbar legend. */
    material?: string;
    /** Unit for the overlay: 'fraction' or 'mass'. Controls the colorbar
     *  range and label. */
    unit?: string;
    /** Scroll-wheel navigation: called with +1 to advance to the next
     *  cross-section, -1 to step back. Optional — omit to disable. */
    onStepTarget?: (delta: number) => void;
  };

  let {
    target,
    targetIndex,
    snakePoints,
    onSnakeUpdate,
    status,
    arcOffsetMm = 0,
    overlay = null,
    material = '',
    unit = 'fraction',
    onStepTarget,
  }: Props = $props();

  /* ── Canvas state ──────────────────────────────────────── */

  let canvasEl: HTMLCanvasElement | undefined = $state();
  let canvasSize = $state(512);
  let busy = $state(false);
  let addPointMode = $state(false);

  /* ── Drag state ────────────────────────────────────────── */

  let dragIndex = $state<number | null>(null);
  let hoverIndex = $state<number | null>(null);

  const WC = 40;
  const WW = 400;
  const CONTROL_POINT_RADIUS = 4; // visual radius in canvas px
  const HIT_RADIUS = 10; // mouse proximity threshold in canvas px

  /* ── Coordinate mapping ────────────────────────────────── */

  function pixelToCanvas(px: number): number {
    return (px / target.pixels) * canvasSize;
  }

  function canvasToPixel(cx: number): number {
    return (cx / canvasSize) * target.pixels;
  }

  /* ── Rendering ─────────────────────────────────────────── */

  /** Turbo-style rainbow colormap — good perceptual ordering, no ambiguous
   *  green band. Input t is clamped to [0, 1]. */
  function jetColor(t: number): [number, number, number] {
    const u = Math.max(0, Math.min(1, t));
    const r = Math.max(0, Math.min(1, 1.5 - Math.abs(4 * u - 3))) * 255;
    const g = Math.max(0, Math.min(1, 1.5 - Math.abs(4 * u - 2))) * 255;
    const b = Math.max(0, Math.min(1, 1.5 - Math.abs(4 * u - 1))) * 255;
    return [Math.round(r), Math.round(g), Math.round(b)];
  }

  /** Scale range for the overlay colormap: [0, 1] for volume fractions,
   *  [0, 1000] mg/mL (≈ water density) for mass densities. Fixed ranges so
   *  contiguous cross-sections share the same color scale. */
  function overlayRange(): [number, number] {
    if (unit === 'mass') return [0, 1000];
    return [0, 1];
  }

  function renderBackground(ctx: CanvasRenderingContext2D) {
    const srcSize = target.pixels;

    const srcCanvas = document.createElement('canvas');
    srcCanvas.width = srcSize;
    srcCanvas.height = srcSize;
    const srcCtx = srcCanvas.getContext('2d')!;
    const imgData = srcCtx.createImageData(srcSize, srcSize);

    const lo = WC - WW / 2;
    const range = WW;
    const [vMin, vMax] = overlayRange();
    const vSpan = vMax - vMin;

    for (let i = 0; i < target.image.length; i++) {
      const hu = target.image[i];
      const gray = Math.max(0, Math.min(255, Math.round(((hu - lo) / range) * 255)));

      const ov = overlay ? overlay[i] : NaN;
      if (overlay && !Number.isNaN(ov)) {
        // Material overlay inside the ROI: rainbow colormap.
        const t = vSpan > 0 ? (ov - vMin) / vSpan : 0;
        const [r, g, b] = jetColor(t);
        imgData.data[i * 4] = r;
        imgData.data[i * 4 + 1] = g;
        imgData.data[i * 4 + 2] = b;
      } else {
        imgData.data[i * 4] = gray;
        imgData.data[i * 4 + 1] = gray;
        imgData.data[i * 4 + 2] = gray;
      }
      imgData.data[i * 4 + 3] = 255;
    }

    srcCtx.putImageData(imgData, 0, 0);

    ctx.imageSmoothingEnabled = true;
    ctx.imageSmoothingQuality = 'high';
    ctx.drawImage(srcCanvas, 0, 0, srcSize, srcSize, 0, 0, canvasSize, canvasSize);
  }

  function drawClosedPolygon(
    ctx: CanvasRenderingContext2D,
    points: [number, number][],
    strokeColor: string,
    lineWidth: number,
    dashed: boolean,
  ) {
    if (points.length < 2) return;
    ctx.save();
    ctx.strokeStyle = strokeColor;
    ctx.lineWidth = lineWidth;
    if (dashed) ctx.setLineDash([6, 4]);
    else ctx.setLineDash([]);
    ctx.beginPath();
    const [x0, y0] = points[0];
    ctx.moveTo(pixelToCanvas(x0), pixelToCanvas(y0));
    for (let i = 1; i < points.length; i++) {
      ctx.lineTo(pixelToCanvas(points[i][0]), pixelToCanvas(points[i][1]));
    }
    ctx.closePath();
    ctx.stroke();
    ctx.restore();
  }

  function drawControlPoints(ctx: CanvasRenderingContext2D, points: [number, number][]) {
    for (let i = 0; i < points.length; i++) {
      const cx = pixelToCanvas(points[i][0]);
      const cy = pixelToCanvas(points[i][1]);
      const isHover = hoverIndex === i;
      const isDrag = dragIndex === i;

      ctx.beginPath();
      ctx.arc(cx, cy, CONTROL_POINT_RADIUS, 0, Math.PI * 2);
      ctx.fillStyle = isDrag ? '#ffffff' : isHover ? '#66ff66' : '#30d158';
      ctx.fill();
      ctx.strokeStyle = '#000000';
      ctx.lineWidth = 1;
      ctx.stroke();
    }
  }

  function render() {
    if (!canvasEl) return;
    const ctx = canvasEl.getContext('2d');
    if (!ctx) return;

    canvasEl.width = canvasSize;
    canvasEl.height = canvasSize;

    // 1. Background image
    renderBackground(ctx);

    // 2. Vessel wall (red dashed)
    if (target.vessel_wall.length > 0) {
      drawClosedPolygon(ctx, target.vessel_wall, '#ff453a', 1.5, true);
    }

    // 3. Snake or init boundary
    if (snakePoints && snakePoints.length > 0) {
      // Active snake contour (green solid)
      drawClosedPolygon(ctx, snakePoints, '#30d158', 2, false);
      drawControlPoints(ctx, snakePoints);
    } else if (target.init_boundary.length > 0) {
      // Init boundary (blue dashed)
      drawClosedPolygon(ctx, target.init_boundary, '#0a84ff', 1.5, true);
    }

    // 4. Add-point mode cursor indicator
    if (addPointMode) {
      ctx.save();
      ctx.fillStyle = 'rgba(255, 214, 10, 0.3)';
      ctx.fillRect(0, canvasSize - 24, canvasSize, 24);
      ctx.fillStyle = '#ffd60a';
      ctx.font = '11px -apple-system, sans-serif';
      ctx.textAlign = 'center';
      ctx.fillText('Click to add point (Esc to cancel)', canvasSize / 2, canvasSize - 8);
      ctx.restore();
    }
  }

  // Re-render when dependencies change.
  // Access reactive state inline to register as $effect dependencies.
  $effect(() => {
    void target.pixels;
    void target.image;
    void snakePoints;
    void hoverIndex;
    void dragIndex;
    void addPointMode;
    void canvasSize;
    void overlay;
    void unit;
    queueMicrotask(() => render());
  });

  /* ── Mouse interaction ─────────────────────────────────── */

  function getCanvasCoords(e: MouseEvent): [number, number] {
    if (!canvasEl) return [0, 0];
    const rect = canvasEl.getBoundingClientRect();
    const scaleX = canvasSize / rect.width;
    const scaleY = canvasSize / rect.height;
    return [
      (e.clientX - rect.left) * scaleX,
      (e.clientY - rect.top) * scaleY,
    ];
  }

  function findNearestPoint(canvasX: number, canvasY: number): number | null {
    if (!snakePoints) return null;
    let minDist = Infinity;
    let minIdx = -1;
    for (let i = 0; i < snakePoints.length; i++) {
      const cx = pixelToCanvas(snakePoints[i][0]);
      const cy = pixelToCanvas(snakePoints[i][1]);
      const dist = Math.hypot(canvasX - cx, canvasY - cy);
      if (dist < minDist) {
        minDist = dist;
        minIdx = i;
      }
    }
    return minDist < HIT_RADIUS ? minIdx : null;
  }

  function handleMouseDown(e: MouseEvent) {
    const [cx, cy] = getCanvasCoords(e);

    // Add-point mode: click adds a new point
    if (addPointMode) {
      const px = canvasToPixel(cx);
      const py = canvasToPixel(cy);
      handleAddPoint([px, py]);
      return;
    }

    // Check if clicking near a control point
    const idx = findNearestPoint(cx, cy);
    if (idx !== null) {
      dragIndex = idx;
      e.preventDefault();
    }
  }

  function handleMouseMove(e: MouseEvent) {
    const [cx, cy] = getCanvasCoords(e);

    if (dragIndex !== null && snakePoints) {
      // Dragging a control point
      const px = canvasToPixel(cx);
      const py = canvasToPixel(cy);
      const updated = snakePoints.map((p, i) =>
        i === dragIndex ? [px, py] as [number, number] : p,
      );
      onSnakeUpdate(updated);
    } else {
      // Hover detection
      hoverIndex = findNearestPoint(cx, cy);
    }
  }

  async function handleMouseUp() {
    if (dragIndex !== null && snakePoints) {
      // Sync dragged points to backend
      try {
        await updateSnakePoints(targetIndex, snakePoints);
      } catch (err) {
        console.error('Failed to sync snake points:', err);
      }
    }
    dragIndex = null;
  }

  function handleKeyDown(e: KeyboardEvent) {
    if (e.key === 'Escape' && addPointMode) {
      addPointMode = false;
      e.stopPropagation();
    }
  }

  /** Accumulated wheel deltaY — wheel events arrive as small fractional
   *  values on trackpads, so we stage them and fire a step when enough
   *  scroll distance has piled up. Threshold tuned so a normal trackpad
   *  swipe advances one cross-section. */
  let wheelAccum = 0;
  const WHEEL_STEP_THRESHOLD = 40;

  function handleWheel(e: WheelEvent) {
    if (!onStepTarget) return;
    // Don't hijack scroll while dragging a point.
    if (dragIndex !== null) return;
    e.preventDefault();
    wheelAccum += e.deltaY;
    while (wheelAccum >= WHEEL_STEP_THRESHOLD) {
      onStepTarget(1);
      wheelAccum -= WHEEL_STEP_THRESHOLD;
    }
    while (wheelAccum <= -WHEEL_STEP_THRESHOLD) {
      onStepTarget(-1);
      wheelAccum += WHEEL_STEP_THRESHOLD;
    }
  }

  /* ── Toolbar actions ───────────────────────────────────── */

  async function handleAddPoint(position: [number, number]) {
    if (busy) return;
    busy = true;
    addPointMode = false;
    try {
      // Backend inserts the point at the nearest edge and returns the new
      // polygon on the very next read; compute that locally to avoid an
      // extra fetch (addSnakePoint only returns the inserted index).
      await addSnakePoint(targetIndex, position);
      const inserted = insertPointLocal(snakePoints ?? [], position);
      onSnakeUpdate(inserted);
    } catch (err) {
      console.error('Add point failed:', err);
    } finally {
      busy = false;
    }
  }

  /** Mirror of `pcat_pipeline::active_contour::insert_control_point` —
   *  find the closest polygon edge and insert `position` after it. */
  function insertPointLocal(
    points: [number, number][],
    position: [number, number],
  ): [number, number][] {
    if (points.length < 2) return [...points, position];
    let best = 0;
    let bestDist = Infinity;
    for (let i = 0; i < points.length; i++) {
      const j = (i + 1) % points.length;
      const d = pointToSegmentDistance(position, points[i], points[j]);
      if (d < bestDist) {
        bestDist = d;
        best = i;
      }
    }
    const out = points.slice();
    out.splice(best + 1, 0, position);
    return out;
  }

  function pointToSegmentDistance(
    p: [number, number],
    a: [number, number],
    b: [number, number],
  ): number {
    const abx = b[0] - a[0];
    const aby = b[1] - a[1];
    const apx = p[0] - a[0];
    const apy = p[1] - a[1];
    const abSq = abx * abx + aby * aby;
    if (abSq < 1e-12) return Math.hypot(apx, apy);
    const t = Math.max(0, Math.min(1, (apx * abx + apy * aby) / abSq));
    const dx = p[0] - (a[0] + t * abx);
    const dy = p[1] - (a[1] + t * aby);
    return Math.hypot(dx, dy);
  }

  function handleAddPointMode() {
    addPointMode = !addPointMode;
  }

  /** Reset = re-adopt the auto-detected vessel wall (resampled). */
  async function handleReset() {
    if (busy) return;
    busy = true;
    try {
      const adopted = await useVesselWallAsContour({ targetIndex });
      if (adopted.length > 0) {
        onSnakeUpdate(adopted[0].points);
      }
    } catch (err) {
      console.error('Reset failed:', err);
    } finally {
      busy = false;
    }
  }
</script>

<svelte:window onkeydown={handleKeyDown} />

<div class="flex min-h-0 flex-1 flex-col overflow-hidden">
  <!-- Canvas area: square, fits whichever of width/height is smaller -->
  <div class="relative flex min-h-0 flex-1 items-center justify-center p-2">
    <div class="relative aspect-square h-full max-h-full w-auto max-w-full">
      <canvas
        bind:this={canvasEl}
        class="h-full w-full cursor-crosshair rounded"
        style="image-rendering: pixelated;"
        onmousedown={handleMouseDown}
        onmousemove={handleMouseMove}
        onmouseup={handleMouseUp}
        onmouseleave={() => { hoverIndex = null; if (dragIndex !== null) { dragIndex = null; } }}
        onwheel={handleWheel}
      ></canvas>

      <!-- Status badge overlay -->
      <div class="absolute right-2 top-2 flex items-center gap-1.5">
        {#if status === 'done'}
          <span class="rounded-full bg-success/20 px-2 py-0.5 text-[10px] font-medium text-success">
            Done
          </span>
        {:else if status === 'in-progress'}
          <span class="rounded-full bg-warning/20 px-2 py-0.5 text-[10px] font-medium text-warning">
            Editing
          </span>
        {:else}
          <span class="rounded-full bg-text-secondary/20 px-2 py-0.5 text-[10px] font-medium text-text-secondary">
            Pending
          </span>
        {/if}
      </div>

      <!-- Frame info overlay -->
      <div class="absolute left-2 top-2 rounded bg-black/40 px-1.5 py-0.5">
        <span class="text-[10px] tabular-nums text-text-primary">
          Frame {target.frame_index} | {(target.arc_mm - arcOffsetMm).toFixed(1)} mm
        </span>
      </div>

      <!-- MMD colorbar legend (only when overlay is showing) -->
      {#if overlay}
        <div class="pointer-events-none absolute bottom-2 right-2 flex items-end gap-1.5">
          <div class="flex flex-col items-end text-[9px] tabular-nums text-white drop-shadow">
            <span>{unit === 'mass' ? '1100' : '100%'}</span>
            <span class="flex-1"></span>
            <span>0</span>
          </div>
          <div
            class="h-20 w-2.5 rounded border border-white/40"
            style="background: linear-gradient(to top,
              rgb(128,  0,   0),
              rgb(255,  0,   0),
              rgb(255,128,   0),
              rgb(255,255,   0),
              rgb(128,255, 128),
              rgb(  0,255, 255),
              rgb(  0,128, 255),
              rgb(  0,  0, 255),
              rgb(  0,  0, 128));"
          ></div>
          <span class="text-[9px] font-medium text-white drop-shadow [writing-mode:vertical-rl] [transform:rotate(180deg)]">
            {material}{unit === 'mass' ? ' (mg/mL)' : ''}
          </span>
        </div>
      {/if}

      <!-- Loading overlay -->
      {#if busy}
        <div class="absolute inset-0 flex items-center justify-center rounded bg-black/30">
          <span class="text-xs text-text-primary">Processing...</span>
        </div>
      {/if}
    </div>
  </div>

  <!-- Toolbar -->
  <div class="flex shrink-0 items-center gap-1.5 border-t border-border bg-surface-secondary px-2 py-1.5">
    <button
      class="rounded px-2.5 py-1 text-xs font-medium transition-colors disabled:bg-surface-tertiary/40 disabled:text-text-secondary/70
             {addPointMode ? 'bg-warning/20 text-warning' : 'bg-accent/10 text-accent hover:bg-accent/20 active:bg-accent/30'}"
      onclick={handleAddPointMode}
      disabled={busy || !snakePoints}
      title="Click canvas to add a control point"
    >
      Add Point
    </button>
    <button
      class="rounded bg-surface-tertiary px-2.5 py-1 text-xs font-medium text-text-primary hover:bg-surface-tertiary/80 disabled:bg-surface-tertiary/40 disabled:text-text-secondary/70"
      onclick={handleReset}
      disabled={busy || target.vessel_wall.length === 0}
      title="Re-adopt the auto-detected vessel wall"
    >
      Reset
    </button>
  </div>
</div>

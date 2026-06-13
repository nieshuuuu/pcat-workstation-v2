<script lang="ts">
  /**
   * Cross-section material map for the MMD analysis view.
   *
   * Renders the HU image as the grayscale CT underlay with the material
   * decomposition (water/lipid fraction or mass density) as a jet overlay on
   * top — the same "jet over CT, gated voxels fall back to grayscale" standard
   * as the whole-volume Water/Lipid view, so the two read identically. The
   * auto-detected lumen wall is drawn as a non-interactive reference outline;
   * there is no manual contour editing (auto-detection is the single source).
   */
  import { MMD_MASS_MAX_MGML, type AnnotationTarget } from '$lib/api';
  import { jet, CT_LO, CT_HI } from '$lib/colormap';

  type Props = {
    target: AnnotationTarget;
    /** Auto-detected lumen contour [x,y] in pixel coords, drawn as a
     *  reference outline. Null hides it. */
    snakePoints: [number, number][] | null;
    /** Absolute arc-length (mm) of the ostium along the centerline.
     *  Displayed arc = target.arc_mm - arcOffsetMm. */
    arcOffsetMm?: number;
    /** Optional material-decomposition overlay (flat pixels×pixels, values in
     *  volume fraction [0,1] or mass density mg/mL). Rendered as a jet colormap
     *  over the CT when non-null. NaN = gated voxel → plain CT. */
    overlay?: number[] | null;
    /** Per-pixel 1σ uncertainty of `overlay`, in the same unit (from
     *  get_mmd_overlay), shown in the hover readout. Null when no overlay. */
    sigma?: number[] | null;
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
    snakePoints,
    arcOffsetMm = 0,
    overlay = null,
    sigma = null,
    material = '',
    unit = 'fraction',
    onStepTarget,
  }: Props = $props();

  /* ── Canvas state ──────────────────────────────────────── */

  let canvasEl: HTMLCanvasElement | undefined = $state();
  let canvasSize = $state(512);

  /* ── Hover readout ─────────────────────────────────────── */

  // Pixel under the cursor → material value ± σ + HU, mirroring the Water/Lipid
  // view. (col,row) is the cross-section pixel; cx,cy position the floating
  // tooltip within the canvas pane.
  let hover = $state<{
    col: number;
    row: number;
    hu: number;
    val: number;
    sig: number;
    cx: number;
    cy: number;
  } | null>(null);
  let matLabel = $derived(
    material ? material.charAt(0).toUpperCase() + material.slice(1) : 'Value',
  );
  // Density is a mass density regardless of the (disabled) unit toggle.
  let effUnit = $derived(material === 'density' ? 'mass' : unit);

  function handleMove(e: MouseEvent) {
    if (!canvasEl) return;
    const rect = canvasEl.getBoundingClientRect();
    const col = Math.floor(((e.clientX - rect.left) / rect.width) * target.pixels);
    const row = Math.floor(((e.clientY - rect.top) / rect.height) * target.pixels);
    if (row < 0 || row >= target.pixels || col < 0 || col >= target.pixels) {
      hover = null;
      return;
    }
    const i = row * target.pixels + col;
    hover = {
      col,
      row,
      hu: target.image[i],
      val: overlay ? overlay[i] : NaN,
      sig: sigma ? sigma[i] : NaN,
      cx: e.clientX - rect.left,
      cy: e.clientY - rect.top,
    };
  }

  /* ── Coordinate mapping ────────────────────────────────── */

  function pixelToCanvas(px: number): number {
    return (px / target.pixels) * canvasSize;
  }

  /** Scale range for the overlay colormap: [0, 1] for volume fractions,
   *  [0, 1000] mg/mL (≈ water density) for mass densities. Fixed ranges so
   *  contiguous cross-sections share one color scale (no per-section auto
   *  rescale), matching the Water/Lipid view's fixed [0,1]. */
  function overlayRange(): [number, number] {
    if (unit === 'mass') return [0, MMD_MASS_MAX_MGML];
    return [0, 1];
  }

  function renderBackground(ctx: CanvasRenderingContext2D) {
    const srcSize = target.pixels;

    const srcCanvas = document.createElement('canvas');
    srcCanvas.width = srcSize;
    srcCanvas.height = srcSize;
    const srcCtx = srcCanvas.getContext('2d')!;
    const imgData = srcCtx.createImageData(srcSize, srcSize);

    const ctSpan = CT_HI - CT_LO;
    const [vMin, vMax] = overlayRange();
    const vSpan = vMax - vMin;

    // Jet material color REPLACES the CT where the overlay is finite; gated
    // (NaN) voxels fall back to grayscale CT. Opaque, like the Water/Lipid view
    // — not an alpha blend — so the material map reads the same in both views.
    for (let i = 0; i < target.image.length; i++) {
      const hu = target.image[i];
      const gray = Math.max(0, Math.min(255, Math.round(((hu - CT_LO) / ctSpan) * 255)));

      const ov = overlay ? overlay[i] : NaN;
      let r: number, g: number, b: number;
      if (overlay && Number.isFinite(ov)) {
        const t = vSpan > 0 ? (ov - vMin) / vSpan : 0;
        [r, g, b] = jet(t);
      } else {
        r = g = b = gray;
      }
      imgData.data[i * 4] = r;
      imgData.data[i * 4 + 1] = g;
      imgData.data[i * 4 + 2] = b;
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

  function render() {
    if (!canvasEl) return;
    const ctx = canvasEl.getContext('2d');
    if (!ctx) return;

    canvasEl.width = canvasSize;
    canvasEl.height = canvasSize;

    // 1. CT + material overlay
    renderBackground(ctx);

    // 2. Auto-detected lumen wall (green outline, non-interactive). Prefer the
    //    adopted contour; fall back to the raw vessel wall, then the init
    //    boundary, purely as a visual reference for where the ring sits.
    if (snakePoints && snakePoints.length > 1) {
      drawClosedPolygon(ctx, snakePoints, '#30d158', 1.75, false);
    } else if (target.vessel_wall.length > 1) {
      drawClosedPolygon(ctx, target.vessel_wall, '#30d158', 1.5, true);
    } else if (target.init_boundary.length > 1) {
      drawClosedPolygon(ctx, target.init_boundary, '#0a84ff', 1.5, true);
    }
  }

  // Re-render when dependencies change.
  $effect(() => {
    void target.pixels;
    void target.image;
    void snakePoints;
    void canvasSize;
    void overlay;
    void unit;
    queueMicrotask(() => render());
  });

  /* ── Scroll-wheel cross-section navigation ─────────────── */

  /** Accumulated wheel deltaY — wheel events arrive as small fractional
   *  values on trackpads, so we stage them and fire a step when enough
   *  scroll distance has piled up. */
  let wheelAccum = 0;
  const WHEEL_STEP_THRESHOLD = 40;

  function handleWheel(e: WheelEvent) {
    if (!onStepTarget) return;
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

  /* ── Colorbar legend ───────────────────────────────────── */

  // 16-stop jet gradient (bottom = range min, top = range max), computed from
  // the same jet() the pixels use so legend and image never drift apart.
  let barStops = $derived(
    Array.from({ length: 17 }, (_, k) => {
      const t = k / 16;
      const [r, g, b] = jet(t);
      return `rgb(${r},${g},${b}) ${(t * 100).toFixed(1)}%`;
    }).join(', '),
  );
  let barMaxLabel = $derived(unit === 'mass' ? String(MMD_MASS_MAX_MGML) : '100%');
  let barUnitLabel = $derived(unit === 'mass' ? ` (mg/mL)` : ' (vol %)');
</script>

<div class="flex min-h-0 flex-1 flex-col overflow-hidden">
  <!-- Canvas area: square, fits whichever of width/height is smaller -->
  <div class="relative flex min-h-0 flex-1 items-center justify-center p-2">
    <div class="relative aspect-square h-full max-h-full w-auto max-w-full">
      <!-- svelte-ignore a11y_no_static_element_interactions -->
      <canvas
        bind:this={canvasEl}
        class="h-full w-full rounded {overlay ? 'cursor-crosshair' : ''}"
        style="image-rendering: pixelated;"
        onwheel={handleWheel}
        onmousemove={handleMove}
        onmouseleave={() => (hover = null)}
      ></canvas>

      <!-- Hover readout: value ± σ for the selected material/unit (mirrors the
           Water/Lipid view's tooltip). -->
      {#if hover}
        <div
          class="pointer-events-none absolute z-10 rounded bg-black/75 px-2 py-1 text-[10px] leading-snug text-white tabular-nums shadow-lg"
          style="left: {hover.cx + 14}px; top: {hover.cy + 14}px;"
        >
          {#if Number.isFinite(hover.val)}
            {#if effUnit === 'mass'}
              <div>{matLabel} {Math.round(hover.val)} ± {Math.round(hover.sig)} mg/mL</div>
            {:else}
              <div>{matLabel} {hover.val.toFixed(2)} ± {hover.sig.toFixed(2)}</div>
            {/if}
            <div class="text-white/55">HU {Math.round(hover.hu)} · ({hover.col}, {hover.row})</div>
          {:else}
            <div class="text-white/55">
              HU {Math.round(hover.hu)}{overlay ? ' · gated' : ''} · ({hover.col}, {hover.row})
            </div>
          {/if}
        </div>
      {/if}

      <!-- Frame info overlay -->
      <div class="absolute left-2 top-2 rounded bg-black/40 px-1.5 py-0.5">
        <span class="text-[10px] tabular-nums text-text-primary">
          Frame {target.frame_index} | {(target.arc_mm - arcOffsetMm).toFixed(1)} mm
        </span>
      </div>

      <!-- MMD colorbar legend (only when an overlay is showing) -->
      {#if overlay}
        <div class="pointer-events-none absolute bottom-2 right-2 flex items-end gap-1.5">
          <div class="flex flex-col items-end justify-between text-[9px] tabular-nums text-white drop-shadow">
            <span>{barMaxLabel}</span>
            <span>0</span>
          </div>
          <div
            class="h-20 w-2.5 rounded border border-white/40"
            style="background: linear-gradient(to top, {barStops});"
          ></div>
          <span class="text-[9px] font-medium text-white drop-shadow [writing-mode:vertical-rl] [transform:rotate(180deg)]">
            {material}{barUnitLabel}
          </span>
        </div>
      {/if}
    </div>
  </div>
</div>

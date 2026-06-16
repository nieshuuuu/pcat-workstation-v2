<script lang="ts">
  /**
   * Single cross-section canvas for a CPR needle position.
   *
   * The cross-section image, lumen polygon and diameter are all computed
   * in Rust (`pcat_pipeline::vessel_wall::compute_vessel_geometry`) and
   * delivered via `batchImageData` / `vesselWall` / `vesselDiameterMm`
   * from the parent. A legacy base64 fallback path is kept for robustness
   * when batch data is unavailable; in that fallback the image still
   * renders but no measurement overlay is shown (Rust is the sole source
   * of truth for the measurement).
   */
  import { invoke } from '@tauri-apps/api/core';

  type Props = {
    centerlineMm: [number, number, number][];
    positionFraction: number;
    rotationDeg: number;
    label: string;
    color: string;
    windowCenter?: number;
    windowWidth?: number;
    /** Pre-computed raw f32 image data from batch call. */
    batchImageData?: Float32Array | null;
    /** Pixel size of the batch image. */
    batchPixels?: number | null;
    /** Pre-computed arc-length in mm from batch call. */
    arcMmProp?: number | null;
    /** Rust-computed lumen diameter in mm (FWHM per-ray scan). */
    vesselDiameterMm?: number | null;
    /** Rust-computed lumen polygon: flat [x0,y0,x1,y1,...] in pixel coords. */
    vesselWall?: Float32Array | null;
    /** Whether to show FAI color overlay. */
    showFaiOverlay?: boolean;
    /** Absolute arc-length (mm) of the ostium along the centerline.
     *  Displayed arc = arc_mm - arcOffsetMm. */
    arcOffsetMm?: number | null;
  };

  let {
    centerlineMm,
    positionFraction,
    rotationDeg,
    label,
    color,
    windowCenter = 40,
    windowWidth = 400,
    batchImageData = null,
    batchPixels = null,
    arcMmProp = null,
    vesselDiameterMm = null,
    vesselWall = null,
    showFaiOverlay = false,
    arcOffsetMm = null,
  }: Props = $props();

  let canvas: HTMLCanvasElement | undefined = $state();
  let arcMm = $state<number | null>(null);
  let loading = $state(false);
  let pixels = $state(128);

  type CrossSectionResult = {
    image_base64: string;
    pixels: number;
    arc_mm: number;
  };

  function decodeBase64Float32(b64: string): Float32Array {
    const binary = atob(b64);
    const bytes = new Uint8Array(binary.length);
    for (let i = 0; i < binary.length; i++) bytes[i] = binary.charCodeAt(i);
    return new Float32Array(bytes.buffer);
  }

  function renderToCanvas(
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
    const hi = wc + ww / 2;
    for (let i = 0; i < data.length; i++) {
      const raw = data[i];
      const gray = Math.max(0, Math.min(255, Math.round(((raw - lo) / (hi - lo)) * 255)));
      if (showFaiOverlay && raw >= -190 && raw <= -30) {
        const t = (raw - (-190)) / ((-30) - (-190));
        const r = Math.round(t < 0.5 ? t * 2 * 255 : 255);
        const g = Math.round(t < 0.5 ? 255 : (1 - (t - 0.5) * 2) * 255);
        imgData.data[i * 4]     = r;
        imgData.data[i * 4 + 1] = g;
        imgData.data[i * 4 + 2] = 20;
      } else {
        imgData.data[i * 4] = gray;
        imgData.data[i * 4 + 1] = gray;
        imgData.data[i * 4 + 2] = gray;
      }
      imgData.data[i * 4 + 3] = 255;
    }
    ctx.putImageData(imgData, 0, 0);
  }

  /**
   * Horizontal / vertical caliper extents of the lumen polygon, in pixel
   * coords. Derived from the polygon's bounding box through its centroid —
   * a light visual echo of the diameter text.
   */
  type Caliper = {
    h: { left: number; right: number; y: number };
    v: { top: number; bottom: number; x: number };
  };

  function polygonCaliper(wall: Float32Array): Caliper | null {
    const n = wall.length / 2;
    if (n < 3) return null;
    let minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity;
    let sumX = 0, sumY = 0;
    for (let i = 0; i < n; i++) {
      const x = wall[2 * i];
      const y = wall[2 * i + 1];
      if (x < minX) minX = x;
      if (x > maxX) maxX = x;
      if (y < minY) minY = y;
      if (y > maxY) maxY = y;
      sumX += x;
      sumY += y;
    }
    return {
      h: { left: minX, right: maxX, y: sumY / n },
      v: { top: minY, bottom: maxY, x: sumX / n },
    };
  }

  let caliper = $derived.by<Caliper | null>(() =>
    vesselWall ? polygonCaliper(vesselWall) : null,
  );

  let wallPolygonPoints = $derived.by<string | null>(() => {
    if (!vesselWall || vesselWall.length < 6) return null;
    const n = vesselWall.length / 2;
    const parts: string[] = new Array(n);
    for (let i = 0; i < n; i++) {
      parts[i] = `${vesselWall[2 * i].toFixed(2)},${vesselWall[2 * i + 1].toFixed(2)}`;
    }
    return parts.join(' ');
  });

  // --- Render pre-computed batch data when provided ---
  $effect(() => {
    const imgData = batchImageData;
    const sz = batchPixels;
    const am = arcMmProp;
    const wc = windowCenter;
    const ww = windowWidth;

    if (imgData && sz && canvas) {
      pixels = sz;
      arcMm = am;
      renderToCanvas(canvas, imgData, sz, sz, wc, ww);
    }
  });

  // --- Fallback: invoke individually if no batch data (legacy base64 path) ---
  let debounceTimer: ReturnType<typeof setTimeout> | undefined;

  $effect(() => {
    // Only self-invoke if no batch data provided
    if (batchImageData) return;

    const cl = centerlineMm;
    const pf = positionFraction;
    const rd = rotationDeg;

    if (!cl || cl.length < 2) return;

    clearTimeout(debounceTimer);
    debounceTimer = setTimeout(async () => {
      loading = true;
      try {
        const centerlineZyx = cl.map(([x, y, z]) => [z, y, x] as [number, number, number]);

        const result = await invoke<CrossSectionResult>(
          'compute_cross_section_image',
          {
            centerlineMm: centerlineZyx,
            positionFraction: pf,
            rotationDeg: rd,
            widthMm: 15.0,
            pixels: 128,
          },
        );

        pixels = result.pixels;
        arcMm = result.arc_mm;

        if (canvas) {
          const data = decodeBase64Float32(result.image_base64);
          renderToCanvas(canvas, data, pixels, pixels, windowCenter, windowWidth);
        }
      } catch (e) {
        console.error(`CrossSection ${label}: computation failed`, e);
      } finally {
        loading = false;
      }
    }, 150);

    return () => clearTimeout(debounceTimer);
  });
</script>

<div class="relative flex flex-col items-center" style="min-height: 0; flex: 1;">
  <!-- Label badge -->
  <div
    class="absolute left-1.5 top-1.5 z-10 flex items-center gap-1"
  >
    <span
      class="inline-block h-3 w-3 rounded-sm text-center text-[10px] font-bold leading-3"
      style="background-color: {color}; color: #000;"
    >
      {label}
    </span>
    {#if arcMm !== null}
      <span class="text-[10px] tabular-nums text-text-secondary">
        {(arcMm - (arcOffsetMm ?? 0)).toFixed(1)} mm
      </span>
    {/if}
  </div>

  <!-- Loading indicator -->
  {#if loading}
    <div class="absolute inset-0 z-20 flex items-center justify-center bg-black/40">
      <span class="text-[10px] text-text-secondary">...</span>
    </div>
  {/if}

  <!-- Canvas -->
  <canvas
    bind:this={canvas}
    class="h-full w-full object-contain"
    width={pixels}
    height={pixels}
    style="image-rendering: pixelated;"
  ></canvas>

  <!-- Diameter measurement overlay -->
  {#if caliper && wallPolygonPoints && vesselDiameterMm !== null}
    <svg class="pointer-events-none absolute inset-0 h-full w-full" viewBox="0 0 {pixels} {pixels}" preserveAspectRatio="xMidYMid meet">
      <!-- Lumen boundary polygon (Rust FWHM per-ray) -->
      <polygon
        points={wallPolygonPoints}
        fill="#22d3ee"
        fill-opacity="0.1"
        stroke="#22d3ee"
        stroke-width="0.8"
        stroke-opacity="0.8"
      />

      <!-- Horizontal caliper -->
      <line x1={caliper.h.left} y1={caliper.h.y} x2={caliper.h.right} y2={caliper.h.y}
        stroke="#facc15" stroke-width="1" stroke-opacity="0.9" />
      <line x1={caliper.h.left} y1={caliper.h.y - 3} x2={caliper.h.left} y2={caliper.h.y + 3}
        stroke="#facc15" stroke-width="1" stroke-opacity="0.9" />
      <line x1={caliper.h.right} y1={caliper.h.y - 3} x2={caliper.h.right} y2={caliper.h.y + 3}
        stroke="#facc15" stroke-width="1" stroke-opacity="0.9" />

      <!-- Vertical caliper -->
      <line x1={caliper.v.x} y1={caliper.v.top} x2={caliper.v.x} y2={caliper.v.bottom}
        stroke="#facc15" stroke-width="1" stroke-opacity="0.9" />
      <line x1={caliper.v.x - 3} y1={caliper.v.top} x2={caliper.v.x + 3} y2={caliper.v.top}
        stroke="#facc15" stroke-width="1" stroke-opacity="0.9" />
      <line x1={caliper.v.x - 3} y1={caliper.v.bottom} x2={caliper.v.x + 3} y2={caliper.v.bottom}
        stroke="#facc15" stroke-width="1" stroke-opacity="0.9" />

      <!-- Diameter label -->
      <text x={caliper.h.right + 3} y={caliper.h.y - 2}
        fill="#facc15" font-size="9" font-family="-apple-system, sans-serif" font-weight="bold">
        {vesselDiameterMm.toFixed(1)} mm
      </text>
    </svg>
  {/if}
</div>

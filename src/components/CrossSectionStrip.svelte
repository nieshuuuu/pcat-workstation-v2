<script lang="ts">
  /**
   * Horizontal scrollable strip of cross-section thumbnails.
   *
   * Each thumbnail renders a 64x64 canvas from the annotation target's HU
   * image data. Clicking a thumbnail selects it for editing in SnakeEditor.
   * Status badges show annotation progress.
   */
  import type { AnnotationTarget } from '$lib/api';

  type Props = {
    targets: AnnotationTarget[];
    selectedIndex: number;
    /** Map of target index to annotation status */
    statusMap: Record<number, 'pending' | 'in-progress' | 'done'>;
    onSelect: (index: number) => void;
    /** Absolute arc-length (mm) of the ostium along the centerline.
     *  Displayed arc = target.arc_mm - arcOffsetMm. */
    arcOffsetMm?: number;
  };

  let { targets, selectedIndex, statusMap, onSelect, arcOffsetMm = 0 }: Props = $props();

  const THUMB_SIZE = 64;
  const WC = 40;
  const WW = 400;

  /** Map of target index -> bound canvas element. */
  let canvasRefs: Record<number, HTMLCanvasElement> = {};

  /** Svelte action to bind a canvas element by index. */
  function bindCanvas(node: HTMLCanvasElement, index: number) {
    canvasRefs[index] = node;
    renderThumbnail(node, targets[index]);
    return {
      update(newIndex: number) {
        delete canvasRefs[index];
        canvasRefs[newIndex] = node;
        if (targets[newIndex]) renderThumbnail(node, targets[newIndex]);
      },
      destroy() {
        delete canvasRefs[index];
      },
    };
  }

  /**
   * Render a cross-section HU image to a 64x64 canvas thumbnail.
   */
  function renderThumbnail(canvas: HTMLCanvasElement, tgt: AnnotationTarget) {
    if (!tgt) return;
    const ctx = canvas.getContext('2d');
    if (!ctx) return;

    const srcSize = tgt.pixels;
    canvas.width = THUMB_SIZE;
    canvas.height = THUMB_SIZE;

    // Create a full-resolution ImageData first, then draw scaled
    const srcCanvas = document.createElement('canvas');
    srcCanvas.width = srcSize;
    srcCanvas.height = srcSize;
    const srcCtx = srcCanvas.getContext('2d')!;
    const imgData = srcCtx.createImageData(srcSize, srcSize);

    const lo = WC - WW / 2;
    const hi = WC + WW / 2;

    for (let i = 0; i < tgt.image.length; i++) {
      const hu = tgt.image[i];
      const gray = Math.max(0, Math.min(255, Math.round(((hu - lo) / (hi - lo)) * 255)));
      imgData.data[i * 4] = gray;
      imgData.data[i * 4 + 1] = gray;
      imgData.data[i * 4 + 2] = gray;
      imgData.data[i * 4 + 3] = 255;
    }

    srcCtx.putImageData(imgData, 0, 0);

    // Draw scaled to thumbnail size
    ctx.imageSmoothingEnabled = true;
    ctx.imageSmoothingQuality = 'medium';
    ctx.drawImage(srcCanvas, 0, 0, srcSize, srcSize, 0, 0, THUMB_SIZE, THUMB_SIZE);
  }

  // Re-render all thumbnails when targets change
  $effect(() => {
    for (let i = 0; i < targets.length; i++) {
      const canvas = canvasRefs[i];
      if (canvas && targets[i]) {
        renderThumbnail(canvas, targets[i]);
      }
    }
  });
</script>

<div class="flex items-start gap-2 overflow-x-auto px-2 py-1.5" role="listbox" aria-label="Cross-section thumbnails">
  {#each targets as target, i}
    {@const status = statusMap[i] ?? 'pending'}
    {@const isSelected = i === selectedIndex}
    <button
      class="group flex shrink-0 flex-col items-center gap-1 rounded p-1 transition-colors
             {isSelected ? 'bg-accent/20' : 'hover:bg-surface-tertiary'}"
      role="option"
      aria-selected={isSelected}
      onclick={() => onSelect(i)}
      title="Frame {target.frame_index} — {(target.arc_mm - arcOffsetMm).toFixed(1)} mm"
    >
      <!-- Thumbnail canvas -->
      <div
        class="overflow-hidden rounded border-2 transition-colors
               {isSelected ? 'border-accent' : 'border-transparent group-hover:border-border'}"
      >
        <canvas
          width={THUMB_SIZE}
          height={THUMB_SIZE}
          style="width: {THUMB_SIZE}px; height: {THUMB_SIZE}px; image-rendering: pixelated;"
          use:bindCanvas={i}
        ></canvas>
      </div>

      <!-- Arc-length label (relative to ostium) -->
      <span class="text-[10px] tabular-nums text-text-secondary">
        {(target.arc_mm - arcOffsetMm).toFixed(0)}mm
      </span>

      <!-- Status badge -->
      <span class="flex items-center justify-center" title={status}>
        {#if status === 'done'}
          <svg class="h-3 w-3 text-success" viewBox="0 0 16 16" fill="currentColor">
            <path d="M8 0a8 8 0 1 1 0 16A8 8 0 0 1 8 0zm3.78 4.97a.75.75 0 0 0-1.06 0L7 8.69 5.28 6.97a.75.75 0 0 0-1.06 1.06l2.25 2.25a.75.75 0 0 0 1.06 0l4.25-4.25a.75.75 0 0 0 0-1.06z"/>
          </svg>
        {:else if status === 'in-progress'}
          <span class="inline-block h-2.5 w-2.5 rounded-full bg-warning"></span>
        {:else}
          <span class="inline-block h-2.5 w-2.5 rounded-full bg-text-secondary/40"></span>
        {/if}
      </span>
    </button>
  {/each}
</div>

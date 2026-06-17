<script lang="ts">
  /**
   * Faint angular (θ) reference overlaid on a cross-section image, matching the
   * radial-angular convention used by the MMD 3D surface plot: θ = 0° at the
   * RIGHT (east), increasing COUNTER-CLOCKWISE (90° top, 180° left, 270° bottom)
   * — i.e. image point = center + r·(cos θ, −sin θ). This lets you read which
   * physical location a θ on the surface plot corresponds to, and vice versa.
   *
   * Drawn as a non-interactive SVG with a square viewBox + `xMidYMid meet`, so it
   * lines up with the square, center-fit cross-section image regardless of the
   * surrounding box (same fit as an `object-contain` / `aspect-square` canvas).
   * Kept low-contrast and confined to the perimeter so it never obscures tissue.
   */
  let { color = '#7dd3fc' }: { color?: string } = $props();

  const C = 50; // center of the 100×100 viewBox
  const R_OUT = 49;
  const R_IN = 45.5;
  const R_LABEL = 40.5;

  /** Image-space point for an angle (deg) at radius r, in the θ=east / CCW frame. */
  function pt(deg: number, r: number): [number, number] {
    const t = (deg * Math.PI) / 180;
    return [C + r * Math.cos(t), C - r * Math.sin(t)];
  }

  const MINOR = [0, 45, 90, 135, 180, 225, 270, 315];
  const LABELED = new Set([0, 90, 180, 270]);
</script>

<svg
  class="pointer-events-none absolute inset-0 h-full w-full"
  viewBox="0 0 100 100"
  preserveAspectRatio="xMidYMid meet"
  aria-hidden="true"
>
  {#each MINOR as d}
    {@const a = pt(d, R_IN)}
    {@const b = pt(d, R_OUT)}
    <line
      x1={a[0]}
      y1={a[1]}
      x2={b[0]}
      y2={b[1]}
      stroke={color}
      stroke-width={LABELED.has(d) ? 0.6 : 0.4}
      opacity={LABELED.has(d) ? 0.5 : 0.3}
    />
  {/each}
  {#each [...LABELED] as d}
    {@const p = pt(d, R_LABEL)}
    <text
      x={p[0]}
      y={p[1]}
      fill={color}
      opacity="0.6"
      font-size="4.2"
      text-anchor="middle"
      dominant-baseline="central"
    >{d}°</text>
  {/each}
</svg>

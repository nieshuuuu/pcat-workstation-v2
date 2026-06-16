/**
 * CPR projection / unprojection utilities.
 *
 * Maps between 3-D world coordinates ([z, y, x] ordering, consistent with
 * the Rust backend) and 2-D canvas coordinates for both straightened and
 * stretched CPR views.
 */

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

export type CprProjectionInfo = {
  total_arc_mm: number;
  total_proj_arc_mm: number;
  half_width_mm: number;
  projection_normal: [number, number, number];
  mid_height_point: [number, number, number];
  dy_mm: number;
  pixels_wide: number;
  pixels_high: number;
  /**
   * Lookup table for `worldToStretchedCpr`, uniformly sampled in projected
   * arc-length. Its length is independent of the render resolution — the
   * backend caps it (~128) so per-seed projection stays cheap.
   */
  proj_col_pts: [number, number, number][];
  arclengths: number[];
  positions: [number, number, number][];
  normals: [number, number, number][];
};

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/** Dot product of two 3-component vectors. */
function dot3(a: [number, number, number], b: [number, number, number]): number {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

/**
 * Linear interpolation between two 3-component vectors.
 * Returns `a * (1 - t) + b * t`.
 */
function lerp3(
  a: [number, number, number],
  b: [number, number, number],
  t: number,
): [number, number, number] {
  return [
    a[0] + (b[0] - a[0]) * t,
    a[1] + (b[1] - a[1]) * t,
    a[2] + (b[2] - a[2]) * t,
  ];
}

/** Subtract two 3-component vectors: `a - b`. */
function sub3(
  a: [number, number, number],
  b: [number, number, number],
): [number, number, number] {
  return [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
}

/** Squared Euclidean distance between two 3-component vectors. */
function distSq3(a: [number, number, number], b: [number, number, number]): number {
  const dz = a[0] - b[0];
  const dy = a[1] - b[1];
  const dx = a[2] - b[2];
  return dz * dz + dy * dy + dx * dx;
}

// Small pixel margin for bounds checks.
const MARGIN_PX = 2;

// ---------------------------------------------------------------------------
// 1. World -> Straightened CPR
// ---------------------------------------------------------------------------

/**
 * Map a 3-D world position to straightened CPR canvas coordinates.
 *
 * Returns `[canvasX, canvasY]` or `null` when the projected point falls
 * outside the canvas (with a small margin).
 */
export function worldToStraightenedCpr(
  worldZyx: [number, number, number],
  info: CprProjectionInfo,
  canvasW: number,
  canvasH: number,
): [number, number] | null {
  const { positions, normals, half_width_mm } = info;
  const n = positions.length;
  if (n === 0) return null;

  // 1. Find closest centerline point.
  let bestIdx = 0;
  let bestDist = distSq3(worldZyx, positions[0]);
  for (let i = 1; i < n; i++) {
    const d = distSq3(worldZyx, positions[i]);
    if (d < bestDist) {
      bestDist = d;
      bestIdx = i;
    }
  }

  // 2. Column index -> canvasX.
  const canvasX = (bestIdx / (n - 1)) * canvasW;

  // 3. Lateral offset = projection of (world - centerline) onto the normal.
  const offset = sub3(worldZyx, positions[bestIdx]);
  const lateralOffset = dot3(offset, normals[bestIdx]);

  // 4. Map lateral offset -> canvasY.
  const canvasY = (0.5 - lateralOffset / (2 * half_width_mm)) * canvasH;

  // 5. Bounds check with margin.
  if (
    canvasX < -MARGIN_PX ||
    canvasX > canvasW + MARGIN_PX ||
    canvasY < -MARGIN_PX ||
    canvasY > canvasH + MARGIN_PX
  ) {
    return null;
  }

  return [canvasX, canvasY];
}

// ---------------------------------------------------------------------------
// 2. World -> Stretched CPR
// ---------------------------------------------------------------------------

/**
 * Map a 3-D world position to stretched CPR canvas coordinates.
 *
 * Returns `[canvasX, canvasY]` or `null` when the projected point falls
 * outside the canvas.
 */
export function worldToStretchedCpr(
  worldZyx: [number, number, number],
  info: CprProjectionInfo,
  canvasW: number,
  canvasH: number,
): [number, number] | null {
  const { projection_normal, mid_height_point, dy_mm, pixels_high, proj_col_pts } = info;
  const n_cols = proj_col_pts.length;
  if (n_cols < 2) return null;

  // 1. Depth along projection_normal (signed distance from mid-height plane).
  const offset = sub3(worldZyx, mid_height_point);
  const depth = dot3(offset, projection_normal);

  // 2. Project onto the mid-height plane.
  const worldProj: [number, number, number] = [
    worldZyx[0] - depth * projection_normal[0],
    worldZyx[1] - depth * projection_normal[1],
    worldZyx[2] - depth * projection_normal[2],
  ];

  // 3. Find the column by projecting worldProj onto each segment of proj_col_pts
  //    and choosing the segment with the smallest perpendicular distance.
  let bestFracIdx = 0;
  let bestPerpDistSq = Infinity;

  for (let i = 0; i < n_cols - 1; i++) {
    const a = proj_col_pts[i];
    const b = proj_col_pts[i + 1];
    const ab: [number, number, number] = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    const ap: [number, number, number] = [worldProj[0] - a[0], worldProj[1] - a[1], worldProj[2] - a[2]];
    const abLenSq = ab[0] * ab[0] + ab[1] * ab[1] + ab[2] * ab[2];
    if (abLenSq === 0) continue;
    const t = Math.max(0, Math.min(1, (ap[0] * ab[0] + ap[1] * ab[1] + ap[2] * ab[2]) / abLenSq));
    // Closest point on segment to worldProj
    const closest: [number, number, number] = [a[0] + t * ab[0], a[1] + t * ab[1], a[2] + t * ab[2]];
    const perpDistSq = distSq3(worldProj, closest);
    if (perpDistSq < bestPerpDistSq) {
      bestPerpDistSq = perpDistSq;
      bestFracIdx = i + t;
    }
  }

  const colFrac = bestFracIdx / (n_cols - 1);
  const canvasX = colFrac * canvasW;

  // 4. Vertical: invert the renderer's row formula.
  // renderer: y_offset_mm = (pixels_high / 2 - row_image) * dy_mm  =>  row_image = pixels_high / 2 - depth / dy_mm
  const rowImage = pixels_high / 2 - depth / dy_mm;
  const canvasY = (rowImage / pixels_high) * canvasH;

  // 5. Bounds check with margin.
  if (
    canvasX < -MARGIN_PX ||
    canvasX > canvasW + MARGIN_PX ||
    canvasY < -MARGIN_PX ||
    canvasY > canvasH + MARGIN_PX
  ) {
    return null;
  }

  return [canvasX, canvasY];
}

// ---------------------------------------------------------------------------
// 3. Straightened CPR -> World
// ---------------------------------------------------------------------------

/**
 * Map a straightened CPR canvas position back to 3-D world coordinates.
 *
 * Uses linear interpolation between centerline positions and normals for
 * fractional column indices.
 */
export function straightenedCprToWorld(
  canvasX: number,
  canvasY: number,
  info: CprProjectionInfo,
  canvasW: number,
  canvasH: number,
): [number, number, number] {
  const { positions, normals, half_width_mm } = info;
  const n = positions.length;

  // 1. Canvas X -> fractional index.
  const frac = canvasX / canvasW;
  const floatIdx = frac * (n - 1);
  const idx0 = Math.max(0, Math.min(n - 2, Math.floor(floatIdx)));
  const t = floatIdx - idx0;

  // 2. Interpolated centerline position.
  const pos = lerp3(positions[idx0], positions[idx0 + 1], t);

  // 3. Canvas Y -> lateral offset.
  const lateral = (0.5 - canvasY / canvasH) * 2 * half_width_mm;

  // 4. Interpolated normal.
  const normal = lerp3(normals[idx0], normals[idx0 + 1], t);

  // 5. Reconstruct world position = pos + lateral * normal.
  return [
    pos[0] + lateral * normal[0],
    pos[1] + lateral * normal[1],
    pos[2] + lateral * normal[2],
  ];
}

// ---------------------------------------------------------------------------
// 4. Stretched CPR -> World
// ---------------------------------------------------------------------------

/**
 * Map a stretched CPR canvas position back to 3-D world coordinates.
 */
export function stretchedCprToWorld(
  canvasX: number,
  canvasY: number,
  info: CprProjectionInfo,
  canvasW: number,
  canvasH: number,
): [number, number, number] {
  const { projection_normal, dy_mm, pixels_high, proj_col_pts } = info;
  const n_cols = proj_col_pts.length;

  // 1. Canvas X -> interpolated point on proj_col_pts.
  const colFrac = canvasX / canvasW;
  const floatIdx = colFrac * (n_cols - 1);
  const i0 = Math.max(0, Math.min(n_cols - 2, Math.floor(floatIdx)));
  const t = floatIdx - i0;
  const projPt = lerp3(proj_col_pts[i0], proj_col_pts[i0 + 1], t);

  // 2. Canvas Y -> depth offset along projection_normal.
  const rowImage = (canvasY / canvasH) * pixels_high;
  const y_offset_mm = (pixels_high / 2 - rowImage) * dy_mm;

  return [
    projPt[0] + y_offset_mm * projection_normal[0],
    projPt[1] + y_offset_mm * projection_normal[1],
    projPt[2] + y_offset_mm * projection_normal[2],
  ];
}

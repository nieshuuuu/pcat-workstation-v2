use ndarray::Array3;

use crate::cpr::CprFrame;

// ---------------------------------------------------------------------------
// Scanline polygon fill
// ---------------------------------------------------------------------------

/// Rasterize a filled polygon into a 2D grid using scanline fill.
///
/// `polygon` is a list of `[row, col]` points (floating point).
/// Returns a list of `(row, col)` integer coordinates inside the polygon.
fn scanline_fill(polygon: &[[f64; 2]], grid_rows: usize, grid_cols: usize) -> Vec<(usize, usize)> {
    if polygon.len() < 3 {
        return Vec::new();
    }

    // Find bounding box in row dimension
    let mut row_min = f64::INFINITY;
    let mut row_max = f64::NEG_INFINITY;
    for pt in polygon {
        row_min = row_min.min(pt[0]);
        row_max = row_max.max(pt[0]);
    }

    let row_start = (row_min.floor() as i64).max(0) as usize;
    let row_end = ((row_max.ceil() as i64).max(0) as usize).min(grid_rows.saturating_sub(1));

    let n = polygon.len();
    let mut result = Vec::new();

    for row in row_start..=row_end {
        let y = row as f64 + 0.5; // scanline at pixel center

        // Find all x-intersections of the scanline with polygon edges
        let mut intersections = Vec::new();
        for i in 0..n {
            let j = (i + 1) % n;
            let (y0, x0) = (polygon[i][0], polygon[i][1]);
            let (y1, x1) = (polygon[j][0], polygon[j][1]);

            // Skip horizontal edges
            if (y1 - y0).abs() < 1e-12 {
                continue;
            }

            // Check if scanline crosses this edge
            if (y < y0.min(y1)) || (y >= y0.max(y1)) {
                continue;
            }

            // Compute x-intersection
            let t = (y - y0) / (y1 - y0);
            let x = x0 + t * (x1 - x0);
            intersections.push(x);
        }

        intersections.sort_by(|a, b| a.partial_cmp(b).unwrap());

        // Fill between pairs of intersections
        for pair in intersections.chunks(2) {
            if pair.len() < 2 {
                break;
            }
            let col_start = (pair[0].ceil() as i64).max(0) as usize;
            let col_end = ((pair[1].floor() as i64).max(0) as usize).min(grid_cols.saturating_sub(1));

            for col in col_start..=col_end {
                result.push((row, col));
            }
        }
    }

    result
}

// ---------------------------------------------------------------------------
// Polygon radial dilation (used to build a ring ROI outside the lumen)
// ---------------------------------------------------------------------------

/// Dilate a closed polygon outward by `dilation_px` pixels along the radial
/// direction from its centroid. Each vertex moves to `centroid + (r + d)/r *
/// (vertex - centroid)`. Works well for roughly convex, centroid-enclosing
/// shapes like a coronary lumen; for non-convex shapes the dilation is
/// approximate but monotonic.
fn dilate_polygon_radial(polygon: &[[f64; 2]], dilation_px: f64) -> Vec<[f64; 2]> {
    if polygon.is_empty() || dilation_px <= 0.0 {
        return polygon.to_vec();
    }
    let n = polygon.len() as f64;
    let cx: f64 = polygon.iter().map(|p| p[0]).sum::<f64>() / n;
    let cy: f64 = polygon.iter().map(|p| p[1]).sum::<f64>() / n;

    polygon
        .iter()
        .map(|&[x, y]| {
            let dx = x - cx;
            let dy = y - cy;
            let r = (dx * dx + dy * dy).sqrt();
            if r < 1e-9 {
                [x, y]
            } else {
                let s = (r + dilation_px) / r;
                [cx + dx * s, cy + dy * s]
            }
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Contour interpolation
// ---------------------------------------------------------------------------

/// Linearly interpolate between two contours of equal length.
///
/// `t` in [0, 1]: t=0 returns `a`, t=1 returns `b`.
fn lerp_contour(a: &[[f64; 2]], b: &[[f64; 2]], t: f64) -> Vec<[f64; 2]> {
    debug_assert_eq!(a.len(), b.len(), "contours must have the same point count");
    a.iter()
        .zip(b.iter())
        .map(|(pa, pb)| {
            [
                pa[0] * (1.0 - t) + pb[0] * t,
                pa[1] * (1.0 - t) + pb[1] * t,
            ]
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Build a 3D boolean mask from snake contours interpolated along a Bishop frame.
///
/// Given a sparse set of annotated lumen contours (e.g. 20 cross-section
/// polygons), this function linearly interpolates between adjacent contours
/// at every frame column, maps each contour from cross-section pixel
/// coordinates to volume voxel coordinates via the Bishop frame, and
/// rasterizes a ring region into the output mask.
///
/// * `ring_outer_mm <= 0` — the mask is the interior of the contour (legacy
///   "lumen fill" mode).
/// * `ring_outer_mm > 0` — the mask is a ring *outside* the contour, from the
///   contour itself out to a radially-dilated copy `ring_outer_mm` farther
///   away. This is the pericoronary region used for PVAT / MMD analysis.
///
/// # Arguments
///
/// * `contours` — Annotated snake contours. Each entry is a `Vec<[f64; 2]>`
///   of `[x, y]` points in cross-section pixel coordinates.
/// * `frame` — Precomputed Bishop frame with positions/normals/binormals.
/// * `cross_section_indices` — Which frame column index each contour
///   corresponds to, length must equal `contours.len()`.
/// * `volume_dims` — `[depth, height, width]` of the output mask.
/// * `spacing` — `[sz, sy, sx]` voxel spacing in mm.
/// * `origin` — `[oz, oy, ox]` volume origin in mm.
/// * `cross_section_width_mm` — Physical width of cross-section images in mm.
/// * `cross_section_pixels` — Pixel dimension of cross-section images.
/// * `ring_outer_mm` — Outer radial extent of the ring in millimeters
///   (measured from each contour vertex). `<= 0` selects legacy lumen-fill.
///
/// # Panics
///
/// Panics if `contours` is empty, if `contours.len() != cross_section_indices.len()`,
/// if not all contours have the same number of points, or if `cross_section_indices`
/// is not strictly increasing.
pub fn build_3d_roi_mask(
    contours: &[Vec<[f64; 2]>],
    frame: &CprFrame,
    cross_section_indices: &[usize],
    volume_dims: [usize; 3],
    spacing: [f64; 3],
    origin: [f64; 3],
    direction: &[f64; 9],
    cross_section_width_mm: f64,
    cross_section_pixels: usize,
    ring_outer_mm: f64,
) -> Array3<bool> {
    assert!(!contours.is_empty(), "need at least one contour");
    assert_eq!(
        contours.len(),
        cross_section_indices.len(),
        "contours and cross_section_indices must have the same length"
    );

    let n_points = contours[0].len();
    for (i, c) in contours.iter().enumerate() {
        assert_eq!(
            c.len(),
            n_points,
            "contour {} has {} points but contour 0 has {}",
            i,
            c.len(),
            n_points,
        );
    }

    for k in 0..cross_section_indices.len().saturating_sub(1) {
        assert!(
            cross_section_indices[k] < cross_section_indices[k + 1],
            "cross_section_indices must be strictly increasing"
        );
    }

    let [depth, height, width] = volume_dims;
    let mut mask = Array3::<bool>::default((depth, height, width));

    let n_cols = frame.positions.len();
    if n_cols == 0 {
        return mask;
    }

    let center_px = cross_section_pixels as f64 / 2.0;
    // cross_section_width_mm is the half-width (center-to-edge); the full image
    // spans 2 * cross_section_width_mm across cross_section_pixels pixels.
    let scale = 2.0 * cross_section_width_mm / cross_section_pixels as f64; // mm per pixel
    let ring_px = if ring_outer_mm > 0.0 {
        ring_outer_mm / scale
    } else {
        0.0
    };

    let inv_spacing = [1.0 / spacing[0], 1.0 / spacing[1], 1.0 / spacing[2]];

    // For each frame column, determine which contour to use (interpolated or clamped).
    // Strategy:
    //   - Before first annotated index: use first contour
    //   - After last annotated index: use last contour
    //   - Between two annotated indices: linearly interpolate
    let n_contours = contours.len();

    for j in 0..n_cols {
        // Find which segment this column belongs to
        let contour_at_j: Vec<[f64; 2]> = if n_contours == 1 {
            // Single contour — use it everywhere
            contours[0].clone()
        } else if j <= cross_section_indices[0] {
            // Before or at first annotated index — clamp to first
            contours[0].clone()
        } else if j >= cross_section_indices[n_contours - 1] {
            // After or at last annotated index — clamp to last
            contours[n_contours - 1].clone()
        } else {
            // Find the bracketing pair
            let mut seg = 0;
            for k in 0..n_contours - 1 {
                if j >= cross_section_indices[k] && j < cross_section_indices[k + 1] {
                    seg = k;
                    break;
                }
            }
            let idx_a = cross_section_indices[seg];
            let idx_b = cross_section_indices[seg + 1];
            let t = (j - idx_a) as f64 / (idx_b - idx_a) as f64;
            lerp_contour(&contours[seg], &contours[seg + 1], t)
        };

        // Convert cross-section pixel coords → world mm → voxel indices
        // Cross-section coords: [x, y] where x is along binormal, y along normal
        // Center of cross-section image is at pixel (center_px, center_px)
        let pos = frame.positions[j]; // [z, y, x] in mm
        let normal = &frame.normals[j];
        let binormal = &frame.binormals[j];

        // Build polygon in voxel coordinates [vz, vy, vx] → we'll use [row=vz, col] but
        // actually the contour lives on a 2D cross-section plane embedded in 3D.
        // We need to project each contour point to 3D world mm, then to voxel coords,
        // then find which slice/row/col it lands on.
        //
        // Since the cross-section plane is generally oblique to the volume axes,
        // we rasterize by projecting the polygon into each affected volume slice.
        //
        // However, a simpler approach: convert each contour point to 3D voxel coords,
        // then for each "voxel slice" the polygon spans, intersect the polygon with
        // that slice plane. This is complex.
        //
        // Pragmatic approach: convert all contour points to voxel coordinates, then
        // for every point *inside* the polygon in the cross-section, mark the
        // corresponding voxel. We do this by rasterizing the polygon in
        // cross-section pixel space and mapping filled pixels to voxel coords.

        // Rasterize in cross-section pixel space. Convert [x, y] to [row, col]
        // = [y, x] for scanline fill.
        let inner_rc: Vec<[f64; 2]> = contour_at_j.iter().map(|pt| [pt[1], pt[0]]).collect();

        let filled: Vec<(usize, usize)> = if ring_px > 0.0 {
            // Ring mode: fill outer (dilated) polygon, subtract inner.
            let outer_xy = dilate_polygon_radial(&contour_at_j, ring_px);
            let outer_rc: Vec<[f64; 2]> =
                outer_xy.iter().map(|pt| [pt[1], pt[0]]).collect();
            let outer_pixels =
                scanline_fill(&outer_rc, cross_section_pixels, cross_section_pixels);
            let inner_pixels: std::collections::HashSet<(usize, usize)> =
                scanline_fill(&inner_rc, cross_section_pixels, cross_section_pixels)
                    .into_iter()
                    .collect();
            outer_pixels
                .into_iter()
                .filter(|p| !inner_pixels.contains(p))
                .collect()
        } else {
            // Lumen-fill mode (legacy): fill the contour interior.
            scanline_fill(&inner_rc, cross_section_pixels, cross_section_pixels)
        };

        for (row, col) in filled {
            // Cross-section pixel (col, row) = (x, y)
            let px_x = col as f64 + 0.5; // pixel center
            let px_y = row as f64 + 0.5;

            // Convert to mm offsets from centerline
            let offset_b = (px_x - center_px) * scale; // binormal offset
            let offset_n = (px_y - center_px) * scale; // normal offset

            // Map to world mm using Bishop frame
            // normal/binormal are nalgebra::Vector3<f64>
            let world_z = pos[0] + offset_n * normal[0] + offset_b * binormal[0];
            let world_y = pos[1] + offset_n * normal[1] + offset_b * binormal[1];
            let world_x = pos[2] + offset_n * normal[2] + offset_b * binormal[2];

            // Convert world mm to voxel indices (honors DICOM IOP).
            let [vzc, vyc, vxc] = crate::types::patient_to_voxel(
                [world_z, world_y, world_x],
                origin,
                inv_spacing,
                direction,
            );
            let vz = vzc.round() as i64;
            let vy = vyc.round() as i64;
            let vx = vxc.round() as i64;

            // Bounds check
            if vz >= 0
                && vy >= 0
                && vx >= 0
                && (vz as usize) < depth
                && (vy as usize) < height
                && (vx as usize) < width
            {
                mask[[vz as usize, vy as usize, vx as usize]] = true;
            }
        }
    }

    mask
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cpr::CprFrame;
    use nalgebra::Vector3;
    use std::f64::consts::PI;

    /// Helper: create a circular contour in cross-section pixel coordinates.
    fn make_circular_contour(
        cx: f64,
        cy: f64,
        radius_px: f64,
        n_points: usize,
    ) -> Vec<[f64; 2]> {
        (0..n_points)
            .map(|i| {
                let theta = 2.0 * PI * (i as f64) / (n_points as f64);
                [cx + radius_px * theta.cos(), cy + radius_px * theta.sin()]
            })
            .collect()
    }

    /// Helper: build a simple straight-line Bishop frame along the z-axis.
    ///
    /// Centerline from (z0, y0, x0) to (z0 + length, y0, x0) with `n_cols` samples.
    /// Normal = [0, 1, 0], Binormal = [0, 0, 1].
    fn make_z_axis_frame(z0: f64, y0: f64, x0: f64, length: f64, n_cols: usize) -> CprFrame {
        let mut positions = Vec::with_capacity(n_cols);
        let mut arclengths = Vec::with_capacity(n_cols);
        let tangents = vec![Vector3::new(1.0, 0.0, 0.0); n_cols];
        let normals = vec![Vector3::new(0.0, 1.0, 0.0); n_cols];
        let binormals = vec![Vector3::new(0.0, 0.0, 1.0); n_cols];

        for j in 0..n_cols {
            let s = length * (j as f64) / ((n_cols - 1) as f64);
            positions.push([z0 + s, y0, x0]);
            arclengths.push(s);
        }

        CprFrame {
            positions,
            tangents,
            normals,
            binormals,
            arclengths,
        }
    }

    // -- Test 1: Circular contours at 3 positions produce correct mask volume --

    #[test]
    fn test_circular_contours_mask_volume() {
        // Setup: volume 64x64x64 with 1mm spacing, origin at (0,0,0).
        // Centerline along z from z=10 to z=50, centered at y=32, x=32.
        let volume_dims = [64, 64, 64];
        let spacing = [1.0, 1.0, 1.0];
        let origin = [0.0, 0.0, 0.0];
        let n_cols = 41; // z = 10..50 inclusive
        let frame = make_z_axis_frame(10.0, 32.0, 32.0, 40.0, n_cols);

        // Cross-section: 64 pixels, half-width 32 mm (= full span 64 mm) → 1 mm/px.
        let cs_width_mm = 32.0;
        let cs_pixels = 64;
        let center = cs_pixels as f64 / 2.0; // 32.0

        // 3 circular contours of radius 5 pixels (= 5 mm), centered at (32, 32).
        let radius_px = 5.0;
        let n_points = 72;
        let contour = make_circular_contour(center, center, radius_px, n_points);

        let contours = vec![contour.clone(), contour.clone(), contour.clone()];
        let cross_section_indices = vec![0, 20, 40]; // beginning, middle, end

        let mask = build_3d_roi_mask(
            &contours,
            &frame,
            &cross_section_indices,
            volume_dims,
            spacing,
            origin,
            &crate::types::IDENTITY_DIRECTION,
            cs_width_mm,
            cs_pixels,
            0.0, // legacy lumen-fill mode
        );

        let total = mask.iter().filter(|&&v| v).count();

        // Expected: ~41 slices * pi * 5^2 ≈ 41 * 78.5 ≈ 3219 voxels
        // Allow generous tolerance for discrete rasterization
        let expected = (n_cols as f64) * PI * radius_px * radius_px;
        let ratio = total as f64 / expected;
        assert!(
            ratio > 0.7 && ratio < 1.4,
            "mask voxel count {} deviates too much from expected {:.0} (ratio {:.2})",
            total,
            expected,
            ratio,
        );

        // Center of vessel at (z=30, y=32, x=32) should be filled
        assert!(mask[[30, 32, 32]], "center of vessel should be filled");

        // Outside vessel at (z=30, y=32, x=40) should NOT be filled (8mm from center > 5mm radius)
        assert!(
            !mask[[30, 32, 40]],
            "point outside vessel should not be filled"
        );
    }

    // -- Test 2: Interpolation produces intermediate contours --

    #[test]
    fn test_interpolation_two_different_circles() {
        // Two contours: small circle (r=3) and large circle (r=9), at indices 0 and 20.
        // Midpoint at index 10 should produce ~r=6 circle.
        let volume_dims = [32, 64, 64];
        let spacing = [1.0, 1.0, 1.0];
        let origin = [0.0, 0.0, 0.0];
        let n_cols = 21;
        let frame = make_z_axis_frame(5.0, 32.0, 32.0, 20.0, n_cols);

        // Half-width 32 mm → full span 64 mm across 64 px → 1 mm/px.
        let cs_width_mm = 32.0;
        let cs_pixels = 64;
        let center = 32.0;

        let n_points = 72;
        let small = make_circular_contour(center, center, 3.0, n_points);
        let large = make_circular_contour(center, center, 9.0, n_points);

        let contours = vec![small, large];
        let cross_section_indices = vec![0, 20];

        let mask = build_3d_roi_mask(
            &contours,
            &frame,
            &cross_section_indices,
            volume_dims,
            spacing,
            origin,
            &crate::types::IDENTITY_DIRECTION,
            cs_width_mm,
            cs_pixels,
            0.0, // legacy lumen-fill mode
        );

        // Count voxels in the slice at z=5 (frame index 0, r=3) and z=15 (frame index 10, r≈6)
        let count_at_start = count_slice(&mask, 5);
        let count_at_mid = count_slice(&mask, 15);
        let count_at_end = count_slice(&mask, 25);

        // r=3 → area ≈ 28, r=6 → area ≈ 113, r=9 → area ≈ 254
        // Midpoint should be larger than start and smaller than end
        assert!(
            count_at_mid > count_at_start,
            "interpolated mid ({}) should be larger than start ({})",
            count_at_mid,
            count_at_start,
        );
        assert!(
            count_at_mid < count_at_end,
            "interpolated mid ({}) should be smaller than end ({})",
            count_at_mid,
            count_at_end,
        );

        // Midpoint area should be roughly pi * 6^2 ≈ 113
        let expected_mid = PI * 6.0 * 6.0;
        let ratio = count_at_mid as f64 / expected_mid;
        assert!(
            ratio > 0.6 && ratio < 1.5,
            "midpoint slice count {} deviates from expected {:.0} (ratio {:.2})",
            count_at_mid,
            expected_mid,
            ratio,
        );
    }

    // -- Test 3: Out-of-bounds voxels are handled gracefully --

    #[test]
    fn test_out_of_bounds_handling() {
        // Place contour near the edge so some points map outside volume bounds.
        let volume_dims = [32, 32, 32];
        let spacing = [1.0, 1.0, 1.0];
        let origin = [0.0, 0.0, 0.0];
        let n_cols = 11;
        // Centerline near corner: y=2, x=2 — the contour (r=5) will extend outside
        let frame = make_z_axis_frame(10.0, 2.0, 2.0, 10.0, n_cols);

        let cs_width_mm = 32.0;
        let cs_pixels = 32;
        let center = 16.0;

        let n_points = 72;
        let contour = make_circular_contour(center, center, 5.0 * 32.0 / 32.0, n_points);

        let contours = vec![contour.clone(), contour.clone()];
        let cross_section_indices = vec![0, 10];

        // Should not panic
        let mask = build_3d_roi_mask(
            &contours,
            &frame,
            &cross_section_indices,
            volume_dims,
            spacing,
            origin,
            &crate::types::IDENTITY_DIRECTION,
            cs_width_mm,
            cs_pixels,
            0.0, // legacy lumen-fill mode
        );

        // Some voxels should be filled (the portion inside the volume)
        let total = mask.iter().filter(|&&v| v).count();
        assert!(
            total > 0,
            "should have some voxels even when contour extends outside volume"
        );

        // But fewer than a full circle (pi*25 ≈ 78 per slice × 11 ≈ 863)
        // since part of the circle is clipped
        let full_circle_total = (n_cols as f64) * PI * 25.0;
        assert!(
            (total as f64) < full_circle_total,
            "clipped mask ({}) should have fewer voxels than full circle ({:.0})",
            total,
            full_circle_total,
        );
    }

    // -- Test 3b: ring mode excludes lumen and fills a band outside --

    #[test]
    fn test_ring_mode_excludes_lumen() {
        let volume_dims = [64, 64, 64];
        let spacing = [1.0, 1.0, 1.0];
        let origin = [0.0, 0.0, 0.0];
        let n_cols = 41;
        let frame = make_z_axis_frame(10.0, 32.0, 32.0, 40.0, n_cols);

        // Half-width 32 mm → full span 64 mm / 64 px → 1 mm/px.
        let cs_width_mm = 32.0;
        let cs_pixels = 64;
        let center = 32.0;

        // Lumen radius 5 px = 5 mm. Ring 3 mm outward → outer radius 8 mm.
        let radius_px = 5.0;
        let ring_outer_mm = 3.0;
        let n_points = 72;
        let contour = make_circular_contour(center, center, radius_px, n_points);
        let contours = vec![contour.clone(), contour.clone(), contour.clone()];
        let cross_section_indices = vec![0, 20, 40];

        let mask = build_3d_roi_mask(
            &contours,
            &frame,
            &cross_section_indices,
            volume_dims,
            spacing,
            origin,
            &crate::types::IDENTITY_DIRECTION,
            cs_width_mm,
            cs_pixels,
            ring_outer_mm,
        );

        // Lumen center (5 mm from centerline at y=32) must NOT be filled.
        assert!(!mask[[30, 32, 32]], "center of lumen should be empty in ring mode");
        // A point 7 mm from center along +y (within ring 5..8 mm) should be filled.
        assert!(mask[[30, 39, 32]], "point inside the ring should be filled");
        // A point 12 mm from center along +y (beyond outer ring) should NOT be filled.
        assert!(!mask[[30, 44, 32]], "point beyond the ring should be empty");

        // Total volume should be positive and consistent with ring annulus area.
        let total = mask.iter().filter(|&&v| v).count();
        let ring_area = std::f64::consts::PI * (8.0_f64.powi(2) - 5.0_f64.powi(2));
        let expected = (n_cols as f64) * ring_area;
        let ratio = total as f64 / expected;
        assert!(
            ratio > 0.6 && ratio < 1.5,
            "ring voxel count {} vs expected {:.0} (ratio {:.2})",
            total,
            expected,
            ratio,
        );
    }

    // -- Test 4: lerp_contour unit test --

    #[test]
    fn test_lerp_contour() {
        let a = vec![[0.0, 0.0], [10.0, 0.0]];
        let b = vec![[0.0, 10.0], [10.0, 10.0]];

        let mid = lerp_contour(&a, &b, 0.5);
        assert_eq!(mid.len(), 2);
        assert!((mid[0][0]).abs() < 1e-10);
        assert!((mid[0][1] - 5.0).abs() < 1e-10);
        assert!((mid[1][0] - 10.0).abs() < 1e-10);
        assert!((mid[1][1] - 5.0).abs() < 1e-10);

        // t=0 returns a
        let at_zero = lerp_contour(&a, &b, 0.0);
        assert!((at_zero[0][1]).abs() < 1e-10);

        // t=1 returns b
        let at_one = lerp_contour(&a, &b, 1.0);
        assert!((at_one[0][1] - 10.0).abs() < 1e-10);
    }

    // -- Test 5: scanline_fill unit test --

    #[test]
    fn test_scanline_fill_square() {
        // A square polygon from (2,2) to (6,6) in row,col coords
        let polygon = vec![
            [2.0, 2.0],
            [2.0, 6.0],
            [6.0, 6.0],
            [6.0, 2.0],
        ];

        let filled = scanline_fill(&polygon, 10, 10);
        // Should fill the interior region
        assert!(!filled.is_empty(), "square should have filled pixels");

        // All pixels should be within the polygon bounds
        for &(r, c) in &filled {
            assert!(r < 10 && c < 10, "pixel ({},{}) out of grid", r, c);
            assert!(
                r >= 2 && r <= 6 && c >= 2 && c <= 6,
                "pixel ({},{}) outside polygon extent [2..6]",
                r,
                c,
            );
        }

        // Pixel (3, 4) should definitely be filled (interior)
        assert!(
            filled.contains(&(3, 4)),
            "interior pixel (3,4) should be filled"
        );
    }

    /// Count the number of `true` voxels in a given z-slice.
    fn count_slice(mask: &Array3<bool>, z: usize) -> usize {
        let shape = mask.shape();
        if z >= shape[0] {
            return 0;
        }
        let mut count = 0;
        for y in 0..shape[1] {
            for x in 0..shape[2] {
                if mask[[z, y, x]] {
                    count += 1;
                }
            }
        }
        count
    }
}

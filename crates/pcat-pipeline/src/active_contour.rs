use ndarray::Array2;
use serde::{Deserialize, Serialize};

/// Parameters controlling snake (active contour) evolution.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SnakeParams {
    /// Elasticity weight (smoothness).
    pub alpha: f64,
    /// Bending weight (prevents sharp corners).
    pub beta: f64,
    /// Image force weight (attracts to edges).
    pub gamma: f64,
    /// Balloon pressure (expand/contract).
    pub balloon: f64,
    /// Time step for evolution.
    pub step_size: f64,
    /// Number of contour points.
    pub n_points: usize,
}

impl Default for SnakeParams {
    fn default() -> Self {
        Self {
            alpha: 0.3,
            beta: 0.1,
            gamma: 1.5,
            balloon: 0.3,
            step_size: 0.5,
            n_points: 72, // 5-degree spacing
        }
    }
}

/// Information about snake convergence after evolution.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConvergenceInfo {
    /// Number of iterations actually performed.
    pub iterations: usize,
    /// Largest point movement in the last iteration.
    pub max_displacement: f64,
    /// Whether the contour converged (max_displacement < 0.1 px).
    pub converged: bool,
}

/// Initialize a circular contour centered at `(cx, cy)` with the given `radius`.
///
/// Points are placed counter-clockwise at equal angular spacing.
pub fn init_circular_contour(cx: f64, cy: f64, radius: f64, n_points: usize) -> Vec<[f64; 2]> {
    let mut contour = Vec::with_capacity(n_points);
    for i in 0..n_points {
        let theta = 2.0 * std::f64::consts::PI * (i as f64) / (n_points as f64);
        contour.push([cx + radius * theta.cos(), cy + radius * theta.sin()]);
    }
    contour
}

/// Compute image gradient fields using a 3x3 Sobel operator.
///
/// Returns `(grad_x, grad_y)` — the gradient of the gradient magnitude.
/// This means the external force attracts the contour toward edges (high gradient
/// magnitude regions).
pub fn compute_gradient_field(image: &Array2<f32>) -> (Array2<f64>, Array2<f64>) {
    let (rows, cols) = image.dim();

    // Step 1: Sobel to get first-order image gradients gx, gy.
    let mut gx = Array2::<f64>::zeros((rows, cols));
    let mut gy = Array2::<f64>::zeros((rows, cols));

    // Sobel kernels:
    //  Kx = [[-1, 0, 1], [-2, 0, 2], [-1, 0, 1]]
    //  Ky = [[-1, -2, -1], [0, 0, 0], [1, 2, 1]]
    for r in 1..rows.saturating_sub(1) {
        for c in 1..cols.saturating_sub(1) {
            let tl = image[[r - 1, c - 1]] as f64;
            let tc = image[[r - 1, c]] as f64;
            let tr = image[[r - 1, c + 1]] as f64;
            let ml = image[[r, c - 1]] as f64;
            let mr = image[[r, c + 1]] as f64;
            let bl = image[[r + 1, c - 1]] as f64;
            let bc = image[[r + 1, c]] as f64;
            let br = image[[r + 1, c + 1]] as f64;

            gx[[r, c]] = -tl + tr - 2.0 * ml + 2.0 * mr - bl + br;
            gy[[r, c]] = -tl - 2.0 * tc - tr + bl + 2.0 * bc + br;
        }
    }

    // Step 2: Gradient magnitude.
    let mut gmag = Array2::<f64>::zeros((rows, cols));
    for r in 0..rows {
        for c in 0..cols {
            gmag[[r, c]] = (gx[[r, c]].powi(2) + gy[[r, c]].powi(2)).sqrt();
        }
    }

    // Step 3: Gradient OF gradient magnitude (what the external force follows).
    // Using the same Sobel operator on gmag.
    let mut grad_gm_x = Array2::<f64>::zeros((rows, cols));
    let mut grad_gm_y = Array2::<f64>::zeros((rows, cols));

    for r in 1..rows.saturating_sub(1) {
        for c in 1..cols.saturating_sub(1) {
            let tl = gmag[[r - 1, c - 1]];
            let tc = gmag[[r - 1, c]];
            let tr = gmag[[r - 1, c + 1]];
            let ml = gmag[[r, c - 1]];
            let mr = gmag[[r, c + 1]];
            let bl = gmag[[r + 1, c - 1]];
            let bc = gmag[[r + 1, c]];
            let br = gmag[[r + 1, c + 1]];

            grad_gm_x[[r, c]] = -tl + tr - 2.0 * ml + 2.0 * mr - bl + br;
            grad_gm_y[[r, c]] = -tl - 2.0 * tc - tr + bl + 2.0 * bc + br;
        }
    }

    // Step 4: Normalize the external-force field so its peak magnitude is 1.
    // Raw Sobel-of-gmag scales with HU range (~thousands for clinical data),
    // which multiplied by gamma * step_size explodes the snake. Normalizing
    // keeps forces in pixel-scale units independent of window/level.
    let mut max_mag = 0.0_f64;
    for r in 0..rows {
        for c in 0..cols {
            let m = (grad_gm_x[[r, c]].powi(2) + grad_gm_y[[r, c]].powi(2)).sqrt();
            if m > max_mag {
                max_mag = m;
            }
        }
    }
    if max_mag > 1e-12 {
        let s = 1.0 / max_mag;
        for r in 0..rows {
            for c in 0..cols {
                grad_gm_x[[r, c]] *= s;
                grad_gm_y[[r, c]] *= s;
            }
        }
    }

    (grad_gm_x, grad_gm_y)
}

/// Evolve a snake contour for `n_iterations`.
///
/// The contour evolves according to:
///   `v_new = v + dt * (F_internal + F_external + F_balloon)`
///
/// Internal forces enforce smoothness (elasticity) and prevents sharp corners (bending).
/// External forces attract the contour to image edges.
/// Balloon force modulated by tissue HU expands into fat and contracts elsewhere.
pub fn evolve_snake(
    contour: &mut Vec<[f64; 2]>,
    grad_x: &Array2<f64>,
    grad_y: &Array2<f64>,
    image: &Array2<f32>,
    params: &SnakeParams,
    n_iterations: usize,
) -> ConvergenceInfo {
    let (rows, cols) = image.dim();
    let max_x = (cols as f64) - 1.0;
    let max_y = (rows as f64) - 1.0;

    let n = contour.len();
    if n < 5 {
        return ConvergenceInfo {
            iterations: 0,
            max_displacement: 0.0,
            converged: true,
        };
    }

    let mut max_disp = 0.0;
    let mut iter_count = 0;

    for iter in 0..n_iterations {
        let normals = compute_normals(contour);
        let old_contour = contour.clone();
        max_disp = 0.0;

        for i in 0..n {
            // Modular neighbor indices for closed contour.
            let im2 = (i + n - 2) % n;
            let im1 = (i + n - 1) % n;
            let ip1 = (i + 1) % n;
            let ip2 = (i + 2) % n;

            // --- Internal forces ---

            // Elasticity: alpha * (v_{i-1} + v_{i+1} - 2*v_i)
            let f_elast_x =
                params.alpha * (old_contour[im1][0] + old_contour[ip1][0] - 2.0 * old_contour[i][0]);
            let f_elast_y =
                params.alpha * (old_contour[im1][1] + old_contour[ip1][1] - 2.0 * old_contour[i][1]);

            // Bending: -beta * (v_{i-2} - 4*v_{i-1} + 6*v_i - 4*v_{i+1} + v_{i+2})
            let f_bend_x = -params.beta
                * (old_contour[im2][0] - 4.0 * old_contour[im1][0] + 6.0 * old_contour[i][0]
                    - 4.0 * old_contour[ip1][0]
                    + old_contour[ip2][0]);
            let f_bend_y = -params.beta
                * (old_contour[im2][1] - 4.0 * old_contour[im1][1] + 6.0 * old_contour[i][1]
                    - 4.0 * old_contour[ip1][1]
                    + old_contour[ip2][1]);

            // --- External (image) forces ---
            // grad_x/grad_y store gradient of gradient magnitude.
            // Contour coordinates: [x, y] where x = column, y = row.
            let f_img_x = params.gamma * sample_field(grad_x, old_contour[i][0], old_contour[i][1]);
            let f_img_y = params.gamma * sample_field(grad_y, old_contour[i][0], old_contour[i][1]);

            // --- Balloon force ---
            // Sample original HU image at the contour point.
            let hu = sample_field_f32(image, old_contour[i][0], old_contour[i][1]);
            // Expand into fat tissue (-190 to -30 HU), contract elsewhere.
            let balloon_sign = if hu >= -190.0 && hu <= -30.0 {
                params.balloon
            } else {
                -params.balloon
            };
            let f_balloon_x = balloon_sign * normals[i][0];
            let f_balloon_y = balloon_sign * normals[i][1];

            // --- Total force and update ---
            let fx = f_elast_x + f_bend_x + f_img_x + f_balloon_x;
            let fy = f_elast_y + f_bend_y + f_img_y + f_balloon_y;

            contour[i][0] = (old_contour[i][0] + params.step_size * fx).clamp(0.0, max_x);
            contour[i][1] = (old_contour[i][1] + params.step_size * fy).clamp(0.0, max_y);

            let dx = contour[i][0] - old_contour[i][0];
            let dy = contour[i][1] - old_contour[i][1];
            let disp = (dx * dx + dy * dy).sqrt();
            if disp > max_disp {
                max_disp = disp;
            }
        }

        iter_count = iter + 1;

        if max_disp < 0.1 {
            return ConvergenceInfo {
                iterations: iter_count,
                max_displacement: max_disp,
                converged: true,
            };
        }
    }

    ConvergenceInfo {
        iterations: iter_count,
        max_displacement: max_disp,
        converged: max_disp < 0.1,
    }
}

/// Insert a new control point at the given position on the contour.
///
/// Finds the closest edge (segment between consecutive points) and inserts the point there.
/// Returns the index of the inserted point.
pub fn insert_control_point(contour: &mut Vec<[f64; 2]>, position: [f64; 2]) -> usize {
    let n = contour.len();
    if n < 2 {
        contour.push(position);
        return contour.len() - 1;
    }

    let mut best_seg = 0;
    let mut best_dist = f64::INFINITY;

    for i in 0..n {
        let j = (i + 1) % n;
        let dist = point_to_segment_distance(position, contour[i], contour[j]);
        if dist < best_dist {
            best_dist = dist;
            best_seg = i;
        }
    }

    // Insert after best_seg (i.e., between best_seg and best_seg+1).
    let insert_idx = best_seg + 1;
    contour.insert(insert_idx, position);
    insert_idx
}

/// Compute the outward normal at each contour point.
///
/// For a counter-clockwise contour, the outward normal at point `i` is the tangent
/// `v_{i+1} - v_{i-1}` rotated 90 degrees clockwise: `(ty, -tx)` normalized.
fn compute_normals(contour: &[[f64; 2]]) -> Vec<[f64; 2]> {
    let n = contour.len();
    let mut normals = Vec::with_capacity(n);

    for i in 0..n {
        let im1 = (i + n - 1) % n;
        let ip1 = (i + 1) % n;

        let tx = contour[ip1][0] - contour[im1][0];
        let ty = contour[ip1][1] - contour[im1][1];

        // Outward normal for CCW contour: rotate tangent 90 deg clockwise => (ty, -tx).
        let nx = ty;
        let ny = -tx;

        let mag = (nx * nx + ny * ny).sqrt();
        if mag > 1e-12 {
            normals.push([nx / mag, ny / mag]);
        } else {
            normals.push([0.0, 0.0]);
        }
    }

    normals
}

/// Bilinear interpolation to sample a 2D `f64` field at non-integer coordinates.
///
/// `x` is the column coordinate, `y` is the row coordinate.
fn sample_field(field: &Array2<f64>, x: f64, y: f64) -> f64 {
    let (rows, cols) = field.dim();
    if rows == 0 || cols == 0 {
        return 0.0;
    }

    let x = x.clamp(0.0, (cols as f64) - 1.0);
    let y = y.clamp(0.0, (rows as f64) - 1.0);

    let x0 = x.floor() as usize;
    let y0 = y.floor() as usize;
    let x1 = (x0 + 1).min(cols - 1);
    let y1 = (y0 + 1).min(rows - 1);

    let fx = x - x0 as f64;
    let fy = y - y0 as f64;

    let v00 = field[[y0, x0]];
    let v10 = field[[y0, x1]];
    let v01 = field[[y1, x0]];
    let v11 = field[[y1, x1]];

    v00 * (1.0 - fx) * (1.0 - fy) + v10 * fx * (1.0 - fy) + v01 * (1.0 - fx) * fy + v11 * fx * fy
}

/// Bilinear interpolation to sample a 2D `f32` field at non-integer coordinates.
fn sample_field_f32(field: &Array2<f32>, x: f64, y: f64) -> f64 {
    let (rows, cols) = field.dim();
    if rows == 0 || cols == 0 {
        return 0.0;
    }

    let x = x.clamp(0.0, (cols as f64) - 1.0);
    let y = y.clamp(0.0, (rows as f64) - 1.0);

    let x0 = x.floor() as usize;
    let y0 = y.floor() as usize;
    let x1 = (x0 + 1).min(cols - 1);
    let y1 = (y0 + 1).min(rows - 1);

    let fx = x - x0 as f64;
    let fy = y - y0 as f64;

    let v00 = field[[y0, x0]] as f64;
    let v10 = field[[y0, x1]] as f64;
    let v01 = field[[y1, x0]] as f64;
    let v11 = field[[y1, x1]] as f64;

    v00 * (1.0 - fx) * (1.0 - fy) + v10 * fx * (1.0 - fy) + v01 * (1.0 - fx) * fy + v11 * fx * fy
}

/// Compute the minimum distance from a point to a line segment.
fn point_to_segment_distance(p: [f64; 2], a: [f64; 2], b: [f64; 2]) -> f64 {
    let abx = b[0] - a[0];
    let aby = b[1] - a[1];
    let apx = p[0] - a[0];
    let apy = p[1] - a[1];

    let ab_sq = abx * abx + aby * aby;
    if ab_sq < 1e-12 {
        // a and b are effectively the same point.
        return (apx * apx + apy * apy).sqrt();
    }

    let t = ((apx * abx + apy * aby) / ab_sq).clamp(0.0, 1.0);
    let proj_x = a[0] + t * abx;
    let proj_y = a[1] + t * aby;

    let dx = p[0] - proj_x;
    let dy = p[1] - proj_y;
    (dx * dx + dy * dy).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_circular_contour_init() {
        let cx = 50.0;
        let cy = 50.0;
        let radius = 20.0;
        let n = 72;
        let contour = init_circular_contour(cx, cy, radius, n);

        assert_eq!(contour.len(), n);

        // Every point should be at the correct radius from center.
        for pt in &contour {
            let dx = pt[0] - cx;
            let dy = pt[1] - cy;
            let r = (dx * dx + dy * dy).sqrt();
            assert!(
                (r - radius).abs() < 1e-10,
                "Point {:?} has radius {r}, expected {radius}",
                pt
            );
        }

        // Centroid should be approximately (cx, cy).
        let mean_x: f64 = contour.iter().map(|p| p[0]).sum::<f64>() / n as f64;
        let mean_y: f64 = contour.iter().map(|p| p[1]).sum::<f64>() / n as f64;
        assert!((mean_x - cx).abs() < 1e-10);
        assert!((mean_y - cy).abs() < 1e-10);
    }

    #[test]
    fn test_gradient_field_known_edge() {
        // Create a 32x32 image with a vertical edge: left half = 0, right half = 100.
        let mut img = Array2::<f32>::zeros((32, 32));
        for r in 0..32 {
            for c in 16..32 {
                img[[r, c]] = 100.0;
            }
        }

        let (gx, gy) = compute_gradient_field(&img);

        // The gradient magnitude should peak at column ~16 (the edge).
        // The gradient-of-gradient-magnitude (returned fields) should point toward that edge.
        // At the edge itself, gx should be large in the interior rows.
        assert_eq!(gx.dim(), (32, 32));
        assert_eq!(gy.dim(), (32, 32));

        // Interior rows: the Sobel of the original image should give large gx near column 16.
        // The gradient of gradient magnitude should also peak near the edge.
        // Check that there is nonzero gradient activity near the edge.
        let edge_activity: f64 = (2..30)
            .map(|r| gx[[r, 15]].abs() + gx[[r, 16]].abs())
            .sum();
        assert!(
            edge_activity > 0.0,
            "Expected nonzero gradient-of-gmag near vertical edge"
        );

        // Far from the edge (col=5), gradient should be zero or negligible.
        let far_activity: f64 = (2..30).map(|r| gx[[r, 5]].abs()).sum();
        assert!(
            far_activity < edge_activity,
            "Far-from-edge gradient should be less than at-edge gradient"
        );
    }

    #[test]
    fn test_snake_convergence_on_circle() {
        // Create a 64x64 image with a bright circle (simulating an edge).
        let size = 64;
        let center = 32.0;
        let circle_radius = 15.0;
        let mut img = Array2::<f32>::zeros((size, size));

        for r in 0..size {
            for c in 0..size {
                let dx = c as f64 - center;
                let dy = r as f64 - center;
                let dist = (dx * dx + dy * dy).sqrt();
                if dist <= circle_radius {
                    // Inside the circle: set to fat-range HU so balloon expands.
                    img[[r, c]] = -100.0;
                } else {
                    img[[r, c]] = 0.0;
                }
            }
        }

        let (gx, gy) = compute_gradient_field(&img);

        // Init snake slightly smaller than the circle.
        let mut contour = init_circular_contour(center, center, circle_radius * 0.6, 72);

        let params = SnakeParams {
            alpha: 0.3,
            beta: 0.1,
            gamma: 1.5,
            balloon: 0.3,
            step_size: 0.5,
            n_points: 72,
        };

        let info = evolve_snake(&mut contour, &gx, &gy, &img, &params, 200);

        // The snake should have evolved outward toward the circle boundary.
        // Check that the average radius is closer to circle_radius than the initial.
        let avg_radius: f64 = contour
            .iter()
            .map(|p| {
                let dx = p[0] - center;
                let dy = p[1] - center;
                (dx * dx + dy * dy).sqrt()
            })
            .sum::<f64>()
            / contour.len() as f64;

        let initial_radius = circle_radius * 0.6;
        assert!(
            avg_radius > initial_radius,
            "Snake should have expanded: avg_radius={avg_radius}, initial={initial_radius}"
        );
        assert!(
            info.iterations > 0,
            "Should have done at least 1 iteration"
        );
    }

    #[test]
    fn test_insert_control_point() {
        // Create a square contour.
        let mut contour = vec![
            [0.0, 0.0],
            [10.0, 0.0],
            [10.0, 10.0],
            [0.0, 10.0],
        ];

        // Insert a point near the middle of the first edge (0,0)-(10,0).
        let idx = insert_control_point(&mut contour, [5.0, 0.5]);

        assert_eq!(contour.len(), 5);
        assert_eq!(idx, 1); // Should be inserted between index 0 and 1.
        assert_eq!(contour[idx], [5.0, 0.5]);
    }

    #[test]
    fn test_normals_circular_contour() {
        // For a CCW circular contour centered at origin, outward normals should point
        // radially outward (same direction as the point from center).
        let cx = 0.0;
        let cy = 0.0;
        let radius = 10.0;
        let n = 36;
        let contour = init_circular_contour(cx, cy, radius, n);
        let normals = compute_normals(&contour);

        assert_eq!(normals.len(), n);

        for i in 0..n {
            let px = contour[i][0] - cx;
            let py = contour[i][1] - cy;
            let pr = (px * px + py * py).sqrt();

            // Unit radial direction.
            let rx = px / pr;
            let ry = py / pr;

            // The normal should align with the radial direction (dot product ~ 1).
            let dot = normals[i][0] * rx + normals[i][1] * ry;
            assert!(
                dot > 0.9,
                "Normal at point {i} should point outward: dot={dot}, normal={:?}, radial=({rx},{ry})",
                normals[i]
            );

            // Normal should be unit length.
            let nmag = (normals[i][0].powi(2) + normals[i][1].powi(2)).sqrt();
            assert!(
                (nmag - 1.0).abs() < 1e-6,
                "Normal should be unit length, got {nmag}"
            );
        }
    }
}

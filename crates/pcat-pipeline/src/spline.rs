/// Natural cubic spline for 1D data parameterized by arc-length.
///
/// Coefficients: on interval [t_i, t_{i+1}], S(s) = a_i + b_i*(s-t_i) + c_i*(s-t_i)^2 + d_i*(s-t_i)^3
/// Uses natural boundary conditions (second derivative = 0 at endpoints).
pub struct CubicSpline1D {
    t: Vec<f64>,
    a: Vec<f64>,
    b: Vec<f64>,
    c: Vec<f64>,
    d: Vec<f64>,
}

impl CubicSpline1D {
    /// Fit a natural cubic spline through (t, y) data.
    /// `t` must be strictly increasing. Panics if lengths differ or < 2 points.
    pub fn new(t: &[f64], y: &[f64]) -> Self {
        let n = t.len();
        assert!(n >= 2, "need at least 2 data points");
        assert_eq!(t.len(), y.len(), "t and y must have same length");

        if n == 2 {
            // Linear segment — no curvature
            let h = t[1] - t[0];
            let slope = if h.abs() > 1e-15 {
                (y[1] - y[0]) / h
            } else {
                0.0
            };
            return Self {
                t: t.to_vec(),
                a: vec![y[0]],
                b: vec![slope],
                c: vec![0.0],
                d: vec![0.0],
            };
        }

        let m = n - 1; // number of intervals

        // Step 1: compute interval lengths h_i and slope differences
        let mut h = vec![0.0; m];
        for i in 0..m {
            h[i] = t[i + 1] - t[i];
        }

        // Step 2: Set up tridiagonal system for c coefficients
        // Natural spline: c[0] = 0, c[n-1] = 0
        // For i = 1..n-2:
        //   h[i-1]*c[i-1] + 2*(h[i-1]+h[i])*c[i] + h[i]*c[i+1]
        //     = 3*((y[i+1]-y[i])/h[i] - (y[i]-y[i-1])/h[i-1])
        let interior = n - 2; // number of interior c values to solve for

        let mut c_all = vec![0.0; n];

        if interior > 0 {
            // Build RHS
            let mut rhs = vec![0.0; interior];
            for i in 0..interior {
                let ii = i + 1; // index into original arrays
                rhs[i] = 3.0 * ((y[ii + 1] - y[ii]) / h[ii] - (y[ii] - y[ii - 1]) / h[ii - 1]);
            }

            // Tridiagonal solve (Thomas algorithm)
            // diagonal: 2*(h[i-1]+h[i]) for i=1..n-2
            // sub-diagonal: h[i] for i=1..n-3 (below main)
            // super-diagonal: h[i] for i=1..n-3 (above main)
            let mut diag = vec![0.0; interior];
            let mut upper = vec![0.0; interior]; // super-diagonal
            let mut lower = vec![0.0; interior]; // sub-diagonal

            for i in 0..interior {
                let ii = i + 1;
                diag[i] = 2.0 * (h[ii - 1] + h[ii]);
                if i + 1 < interior {
                    upper[i] = h[ii];
                    lower[i + 1] = h[ii];
                }
            }

            // Forward elimination
            let mut diag_m = diag.clone();
            let mut rhs_m = rhs.clone();
            for i in 1..interior {
                let factor = lower[i] / diag_m[i - 1];
                diag_m[i] -= factor * upper[i - 1];
                rhs_m[i] -= factor * rhs_m[i - 1];
            }

            // Back substitution
            let mut c_interior = vec![0.0; interior];
            c_interior[interior - 1] = rhs_m[interior - 1] / diag_m[interior - 1];
            for i in (0..interior - 1).rev() {
                c_interior[i] = (rhs_m[i] - upper[i] * c_interior[i + 1]) / diag_m[i];
            }

            for i in 0..interior {
                c_all[i + 1] = c_interior[i];
            }
        }

        // Step 3: compute a, b, d from c
        let mut a = vec![0.0; m];
        let mut b = vec![0.0; m];
        let mut c_coeff = vec![0.0; m];
        let mut d_coeff = vec![0.0; m];

        for i in 0..m {
            a[i] = y[i];
            c_coeff[i] = c_all[i];
            d_coeff[i] = (c_all[i + 1] - c_all[i]) / (3.0 * h[i]);
            b[i] = (y[i + 1] - y[i]) / h[i] - h[i] * (2.0 * c_all[i] + c_all[i + 1]) / 3.0;
        }

        Self {
            t: t.to_vec(),
            a,
            b,
            c: c_coeff,
            d: d_coeff,
        }
    }

    /// Find the interval index for parameter s (clamped to domain).
    #[inline]
    fn find_interval(&self, s: f64) -> (usize, f64) {
        let n_intervals = self.a.len();
        let s_clamped = s.clamp(self.t[0], self.t[self.t.len() - 1]);

        // Binary search for the right interval
        let mut lo = 0usize;
        let mut hi = n_intervals; // intervals indexed 0..n_intervals-1
        while lo < hi {
            let mid = (lo + hi) / 2;
            if self.t[mid + 1] < s_clamped {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        let i = lo.min(n_intervals - 1);
        let ds = s_clamped - self.t[i];
        (i, ds)
    }

    /// Evaluate S(s).
    pub fn eval(&self, s: f64) -> f64 {
        let (i, ds) = self.find_interval(s);
        self.a[i] + self.b[i] * ds + self.c[i] * ds * ds + self.d[i] * ds * ds * ds
    }

    /// Evaluate S'(s) — first derivative.
    pub fn deriv(&self, s: f64) -> f64 {
        let (i, ds) = self.find_interval(s);
        self.b[i] + 2.0 * self.c[i] * ds + 3.0 * self.d[i] * ds * ds
    }
}

/// 3D cubic spline through [z,y,x] points parameterized by cumulative arc-length.
pub struct CubicSpline3D {
    spline_z: CubicSpline1D,
    spline_y: CubicSpline1D,
    spline_x: CubicSpline1D,
    total_arc: f64,
}

impl CubicSpline3D {
    /// Fit a 3D cubic spline through a sequence of [z,y,x] points.
    /// Arc-length is computed as cumulative Euclidean distance between consecutive points.
    pub fn fit(points: &[[f64; 3]]) -> Self {
        assert!(points.len() >= 2, "need at least 2 points for spline");

        // Compute cumulative arc-length
        let mut arc = Vec::with_capacity(points.len());
        arc.push(0.0);
        for i in 1..points.len() {
            let dz = points[i][0] - points[i - 1][0];
            let dy = points[i][1] - points[i - 1][1];
            let dx = points[i][2] - points[i - 1][2];
            let d = (dz * dz + dy * dy + dx * dx).sqrt();
            arc.push(arc[i - 1] + d);
        }
        let total_arc = *arc.last().unwrap();

        // Extract per-dimension arrays
        let zs: Vec<f64> = points.iter().map(|p| p[0]).collect();
        let ys: Vec<f64> = points.iter().map(|p| p[1]).collect();
        let xs: Vec<f64> = points.iter().map(|p| p[2]).collect();

        Self {
            spline_z: CubicSpline1D::new(&arc, &zs),
            spline_y: CubicSpline1D::new(&arc, &ys),
            spline_x: CubicSpline1D::new(&arc, &xs),
            total_arc,
        }
    }

    /// Evaluate position [z,y,x] at arc-length s.
    pub fn eval(&self, s: f64) -> [f64; 3] {
        [
            self.spline_z.eval(s),
            self.spline_y.eval(s),
            self.spline_x.eval(s),
        ]
    }

    /// Analytic tangent vector at arc-length s, normalized.
    /// Returns [dz, dy, dx] / |[dz, dy, dx]|.
    pub fn tangent(&self, s: f64) -> [f64; 3] {
        let dz = self.spline_z.deriv(s);
        let dy = self.spline_y.deriv(s);
        let dx = self.spline_x.deriv(s);
        let norm = (dz * dz + dy * dy + dx * dx).sqrt();
        if norm < 1e-15 {
            // Degenerate — return arbitrary unit vector
            return [1.0, 0.0, 0.0];
        }
        [dz / norm, dy / norm, dx / norm]
    }

    /// Total arc-length of the fitted spline.
    pub fn total_arc(&self) -> f64 {
        self.total_arc
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cubic_spline_1d_linear() {
        // y = 2*t + 1, should be reproduced exactly
        let t = vec![0.0, 1.0, 2.0, 3.0];
        let y = vec![1.0, 3.0, 5.0, 7.0];
        let s = CubicSpline1D::new(&t, &y);

        assert!((s.eval(0.0) - 1.0).abs() < 1e-10);
        assert!((s.eval(1.5) - 4.0).abs() < 1e-10);
        assert!((s.eval(3.0) - 7.0).abs() < 1e-10);

        // Derivative should be 2.0 everywhere
        assert!((s.deriv(0.5) - 2.0).abs() < 1e-10);
        assert!((s.deriv(2.0) - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_cubic_spline_1d_two_points() {
        let t = vec![0.0, 5.0];
        let y = vec![10.0, 20.0];
        let s = CubicSpline1D::new(&t, &y);

        assert!((s.eval(0.0) - 10.0).abs() < 1e-10);
        assert!((s.eval(5.0) - 20.0).abs() < 1e-10);
        assert!((s.eval(2.5) - 15.0).abs() < 1e-10);
        assert!((s.deriv(2.5) - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_cubic_spline_1d_interpolates() {
        let t = vec![0.0, 1.0, 2.0, 3.0, 4.0];
        let y = vec![0.0, 1.0, 0.0, -1.0, 0.0];
        let s = CubicSpline1D::new(&t, &y);

        // Should pass through all data points
        for i in 0..t.len() {
            assert!(
                (s.eval(t[i]) - y[i]).abs() < 1e-10,
                "spline should interpolate at knot {}: got {} expected {}",
                i,
                s.eval(t[i]),
                y[i]
            );
        }
    }

    #[test]
    fn test_cubic_spline_3d_straight_line() {
        // Straight line along z-axis
        let pts: Vec<[f64; 3]> = (0..10).map(|i| [i as f64, 0.0, 0.0]).collect();
        let spline = CubicSpline3D::fit(&pts);

        assert!((spline.total_arc() - 9.0).abs() < 1e-10);

        // Midpoint
        let mid = spline.eval(4.5);
        assert!((mid[0] - 4.5).abs() < 1e-6);
        assert!(mid[1].abs() < 1e-6);
        assert!(mid[2].abs() < 1e-6);

        // Tangent should be [1, 0, 0]
        let t = spline.tangent(4.5);
        assert!((t[0] - 1.0).abs() < 1e-6);
        assert!(t[1].abs() < 1e-6);
        assert!(t[2].abs() < 1e-6);
    }

    #[test]
    fn test_cubic_spline_3d_endpoints() {
        let pts = vec![
            [0.0, 0.0, 0.0],
            [5.0, 3.0, 4.0],
            [10.0, 0.0, 0.0],
        ];
        let spline = CubicSpline3D::fit(&pts);

        let start = spline.eval(0.0);
        assert!((start[0] - 0.0).abs() < 1e-10);
        assert!((start[1] - 0.0).abs() < 1e-10);
        assert!((start[2] - 0.0).abs() < 1e-10);

        let end = spline.eval(spline.total_arc());
        assert!((end[0] - 10.0).abs() < 1e-10);
        assert!((end[1] - 0.0).abs() < 1e-10);
        assert!((end[2] - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_tangent_is_normalized() {
        let pts = vec![
            [0.0, 0.0, 0.0],
            [3.0, 4.0, 0.0],
            [6.0, 0.0, 5.0],
            [10.0, 2.0, 3.0],
        ];
        let spline = CubicSpline3D::fit(&pts);

        for i in 0..20 {
            let s = spline.total_arc() * (i as f64) / 19.0;
            let t = spline.tangent(s);
            let norm = (t[0] * t[0] + t[1] * t[1] + t[2] * t[2]).sqrt();
            assert!(
                (norm - 1.0).abs() < 1e-10,
                "tangent at s={} should be unit length, got norm={}",
                s,
                norm
            );
        }
    }
}

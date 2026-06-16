use ndarray::Array3;

/// Trilinear interpolation on a (Z, Y, X) f32 volume.
/// Coords are in voxel space (float). Returns NAN if out of bounds.
#[inline]
pub fn trilinear(vol: &Array3<f32>, z: f64, y: f64, x: f64) -> f32 {
    let (nz, ny, nx) = {
        let s = vol.shape();
        (s[0], s[1], s[2])
    };

    // Bounds check — need at least one voxel margin for interpolation
    if z < 0.0 || y < 0.0 || x < 0.0
        || z >= (nz - 1) as f64
        || y >= (ny - 1) as f64
        || x >= (nx - 1) as f64
    {
        return f32::NAN;
    }

    let z0 = z as usize;
    let y0 = y as usize;
    let x0 = x as usize;

    let z1 = z0 + 1;
    let y1 = y0 + 1;
    let x1 = x0 + 1;

    let zf = (z - z0 as f64) as f32;
    let yf = (y - y0 as f64) as f32;
    let xf = (x - x0 as f64) as f32;

    // Read 8 corner values
    let c000 = vol[[z0, y0, x0]];
    let c001 = vol[[z0, y0, x1]];
    let c010 = vol[[z0, y1, x0]];
    let c011 = vol[[z0, y1, x1]];
    let c100 = vol[[z1, y0, x0]];
    let c101 = vol[[z1, y0, x1]];
    let c110 = vol[[z1, y1, x0]];
    let c111 = vol[[z1, y1, x1]];

    // Bilinear on bottom face (z0)
    let c00 = c000 * (1.0 - xf) + c001 * xf;
    let c01 = c010 * (1.0 - xf) + c011 * xf;
    let c0 = c00 * (1.0 - yf) + c01 * yf;

    // Bilinear on top face (z1)
    let c10 = c100 * (1.0 - xf) + c101 * xf;
    let c11 = c110 * (1.0 - xf) + c111 * xf;
    let c1 = c10 * (1.0 - yf) + c11 * yf;

    // Linear between faces
    c0 * (1.0 - zf) + c1 * zf
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::Array3;

    #[test]
    fn test_trilinear_center() {
        // 2x2x2 volume with known values
        let mut vol = Array3::<f32>::zeros((2, 2, 2));
        vol[[0, 0, 0]] = 0.0;
        vol[[0, 0, 1]] = 1.0;
        vol[[0, 1, 0]] = 2.0;
        vol[[0, 1, 1]] = 3.0;
        vol[[1, 0, 0]] = 4.0;
        vol[[1, 0, 1]] = 5.0;
        vol[[1, 1, 0]] = 6.0;
        vol[[1, 1, 1]] = 7.0;

        // Center of the cube should be the average
        let val = trilinear(&vol, 0.5, 0.5, 0.5);
        assert!((val - 3.5).abs() < 1e-5);
    }

    #[test]
    fn test_trilinear_corners() {
        let mut vol = Array3::<f32>::zeros((3, 3, 3));
        vol[[0, 0, 0]] = 100.0;
        vol[[1, 1, 1]] = 200.0;

        assert!((trilinear(&vol, 0.0, 0.0, 0.0) - 100.0).abs() < 1e-5);
        assert!((trilinear(&vol, 1.0, 1.0, 1.0) - 200.0).abs() < 1e-5);
    }

    #[test]
    fn test_trilinear_out_of_bounds() {
        let vol = Array3::<f32>::zeros((3, 3, 3));
        assert!(trilinear(&vol, -0.1, 0.0, 0.0).is_nan());
        assert!(trilinear(&vol, 2.0, 0.0, 0.0).is_nan()); // nz-1 = 2, so >= 2.0 is OOB
        assert!(trilinear(&vol, 0.0, 0.0, 2.5).is_nan());
    }
}

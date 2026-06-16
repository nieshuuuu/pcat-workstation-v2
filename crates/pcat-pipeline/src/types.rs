use ndarray::Array3;
use std::sync::Arc;

/// CT volume loaded from DICOM, stored in Rust memory.
/// `data` is wrapped in `Arc` so commands can share it without
/// cloning ~300MB on every render call.
#[derive(Clone)]
pub struct LoadedVolume {
    pub data: Arc<Array3<f32>>,  // (Z, Y, X) HU values -- shared, not cloned
    pub spacing: [f64; 3],       // [sz, sy, sx] mm
    pub origin: [f64; 3],        // [oz, oy, ox] mm
    pub direction: [f64; 9],     // row-major 3x3
    pub window_center: f64,
    pub window_width: f64,
    pub patient_name: String,
    pub study_description: String,
}

/// Convert a patient-frame position (ZYX order, mm) into continuous voxel
/// indices (vz, vy, vx) using the volume's `origin`, `spacing`, and
/// `direction` (DICOM IOP).
///
/// `direction` is row-major 3×3 with rows `[iop_row, iop_col, slice_normal]` —
/// each row is a unit vector in patient XYZ. Because it is orthonormal, the
/// inverse is the transpose, which collapses to the three dot products below.
///
/// For axial CT with IOP `[1,0,0,0,1,0]` (the common case) this produces
/// bit-identical results to the previous component-wise `(s − o) / spacing`
/// because the matrix-vector product picks each component unchanged.
#[inline]
pub fn patient_to_voxel(
    sample_zyx: [f64; 3],
    origin: [f64; 3],
    inv_spacing: [f64; 3],
    direction: &[f64; 9],
) -> [f64; 3] {
    // Displacement in patient XYZ order (LPS), reversed from ZYX storage.
    let dx = sample_zyx[2] - origin[2];
    let dy = sample_zyx[1] - origin[1];
    let dz = sample_zyx[0] - origin[0];

    // Each row of `direction` is a unit IOP basis vector in patient XYZ.
    // Dot(row, [dx, dy, dz]) gives the displacement component along that
    // basis: row 0 → voxel-X axis, row 1 → voxel-Y axis, row 2 → voxel-Z axis.
    let vx_mm = direction[0] * dx + direction[1] * dy + direction[2] * dz;
    let vy_mm = direction[3] * dx + direction[4] * dy + direction[5] * dz;
    let vz_mm = direction[6] * dx + direction[7] * dy + direction[8] * dz;

    [vz_mm * inv_spacing[0], vy_mm * inv_spacing[1], vx_mm * inv_spacing[2]]
}

/// Identity direction matrix (axial CT default).  Useful for callers that do
/// not yet thread `direction` and want explicit-identity behavior.
pub const IDENTITY_DIRECTION: [f64; 9] = [
    1.0, 0.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 0.0, 1.0,
];

#[cfg(test)]
mod tests {
    use super::*;

    /// Identity direction + zero origin: voxel index = sample_zyx / spacing.
    #[test]
    fn axial_identity_origin_zero() {
        let inv_spacing = [1.0, 0.5, 0.25];
        let v = patient_to_voxel(
            [10.0, 20.0, 40.0],
            [0.0, 0.0, 0.0],
            inv_spacing,
            &IDENTITY_DIRECTION,
        );
        assert_eq!(v, [10.0, 10.0, 10.0]);
    }

    /// Identity direction with non-zero IPP: voxel index must subtract
    /// origin in *every* axis. The pre-fix code dropped X and Y, which would
    /// have returned [10, 100, 200] here instead of [0, 0, 0].
    #[test]
    fn axial_identity_full_ipp_subtracted() {
        let origin = [10.0, 50.0, 50.0]; // ZYX patient mm
        let inv_spacing = [1.0, 1.0, 1.0];
        let v = patient_to_voxel(
            [10.0, 50.0, 50.0],
            origin,
            inv_spacing,
            &IDENTITY_DIRECTION,
        );
        // sample == origin, so voxel index must be zero on every axis.
        assert_eq!(v, [0.0, 0.0, 0.0]);
    }

    /// Tilted direction: stepping +1mm along the IOP row direction must
    /// advance voxel-X by 1, even when that direction is not aligned with
    /// patient X. Pre-fix this silently mis-mapped.
    #[test]
    fn oblique_iop_indexes_along_basis() {
        // IOP row = +Y patient axis; IOP col = -X patient axis;
        // slice normal = +Z patient axis. (90° in-plane rotation.)
        let direction: [f64; 9] = [
            0.0, 1.0, 0.0,    // iop_row → patient +Y
            -1.0, 0.0, 0.0,   // iop_col → patient -X
            0.0, 0.0, 1.0,    // normal  → patient +Z
        ];
        let origin = [0.0, 0.0, 0.0];
        let inv_spacing = [1.0, 1.0, 1.0];
        // Step +1mm in patient Y. With IOP row = +Y, that should be voxel-X +1.
        let v = patient_to_voxel(
            [0.0, 1.0, 0.0],   // ZYX = (z=0, y=1, x=0) patient mm
            origin,
            inv_spacing,
            &direction,
        );
        assert!((v[2] - 1.0).abs() < 1e-12, "voxel x should be 1, got {v:?}");
        assert!(v[0].abs() < 1e-12, "voxel z should be 0, got {v:?}");
        assert!(v[1].abs() < 1e-12, "voxel y should be 0, got {v:?}");
    }

    /// Oblique IOP combined with non-zero origin: both bug fixes interact.
    #[test]
    fn oblique_iop_with_offset_origin() {
        let direction: [f64; 9] = [
            0.0, 1.0, 0.0,
            -1.0, 0.0, 0.0,
            0.0, 0.0, 1.0,
        ];
        let origin = [5.0, 7.0, 11.0]; // ZYX
        let inv_spacing = [1.0, 1.0, 1.0];
        // Walk +2mm along iop_row (patient +Y) from origin: sample_zyx = (5, 9, 11).
        let v = patient_to_voxel(
            [5.0, 9.0, 11.0],
            origin,
            inv_spacing,
            &direction,
        );
        assert!((v[2] - 2.0).abs() < 1e-12, "voxel x should be 2, got {v:?}");
        assert!(v[1].abs() < 1e-12 && v[0].abs() < 1e-12);
    }
}

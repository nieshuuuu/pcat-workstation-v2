//! Shared result type for material decomposition.
//!
//! Returned by the noise-aware water/lipid GLS solver
//! ([`super::water_lipid::decompose_volume_gls`]). Kept in its own module so it
//! is not coupled to any particular solver. Water/lipid only — the dropped
//! 3-material solvers' iodine/calcium fields are gone.

use ndarray::Array3;

/// Per-voxel water/lipid fractions and mass concentrations over a decomposed volume.
pub struct MmdResult {
    pub water_frac: Array3<f32>,
    pub lipid_frac: Array3<f32>,
    pub water_mass: Array3<f32>, // mg/mL
    pub lipid_mass: Array3<f32>, // mg/mL
    pub total_density: Array3<f32>,
    pub mask: Array3<bool>,
    pub iterations: usize,
    pub converged: bool,
}

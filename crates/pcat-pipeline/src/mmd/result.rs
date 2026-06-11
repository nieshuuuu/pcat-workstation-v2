//! Shared result type for material decomposition.
//!
//! Returned by the noise-aware water/lipid GLS solver
//! ([`super::water_lipid::decompose_volume_gls`]). Kept in its own module so it
//! is not coupled to any particular solver. The `iodine`/`calcium` fields are
//! retained in the shape for the overlay/surface plumbing but are zero under the
//! water/lipid model (no iodine or calcium basis).

use ndarray::Array3;

/// Per-voxel material fractions and mass concentrations over a decomposed volume.
pub struct MmdResult {
    pub water_frac: Array3<f32>,
    pub lipid_frac: Array3<f32>,
    /// Zero under the water/lipid model (no iodine basis).
    pub iodine_frac: Array3<f32>,
    /// Zero under the water/lipid model (no calcium basis).
    pub calcium_frac: Array3<f32>,
    pub water_mass: Array3<f32>,   // mg/mL
    pub lipid_mass: Array3<f32>,   // mg/mL
    pub iodine_mass: Array3<f32>,  // mg/mL (zero)
    pub calcium_mass: Array3<f32>, // mg/mL (zero)
    pub total_density: Array3<f32>,
    pub mask: Array3<bool>,
    pub iterations: usize,
    pub converged: bool,
}

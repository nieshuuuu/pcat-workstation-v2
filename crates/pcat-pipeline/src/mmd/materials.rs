//! Physical material constants for the water/lipid decomposition.
//!
//! Intrinsic mass densities (g/cm^3). Canonical home — `water_lipid.rs`
//! references `DENSITY_WATER` / `DENSITY_LIPID` to convert volume fractions to
//! mass concentrations (mg/mL = g/cm^3 × 1000). The former `MaterialLibrary` /
//! `Material` LAC machinery lived here for the 3-material direct/PWSQS solvers,
//! which were dropped (ill-conditioned at soft-tissue contrast).

pub(crate) const DENSITY_WATER: f64 = 1.000;
pub(crate) const DENSITY_LIPID: f64 = 0.950;

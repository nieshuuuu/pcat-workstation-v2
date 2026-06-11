use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Basis materials for multi-material decomposition.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Material {
    Water,
    Lipid,
    Iodine,
    Calcium,
}

impl Material {
    /// All four basis materials in canonical order.
    pub const ALL: [Material; 4] = [
        Material::Water,
        Material::Lipid,
        Material::Iodine,
        Material::Calcium,
    ];
}

/// Linear attenuation coefficients and physical properties for basis materials.
///
/// Stores LAC values (cm^-1) at a specific dual-energy keV pair for each
/// material, plus intrinsic mass densities. Provides helpers for building
/// the system matrix used in multi-material decomposition.
pub struct MaterialLibrary {
    /// Low energy in keV.
    pub low_kev: f64,
    /// High energy in keV.
    pub high_kev: f64,
    /// LAC in cm^-1 at [low_kev, high_kev] for each material.
    lac: HashMap<Material, [f64; 2]>,
    /// Intrinsic mass density in g/cm^3.
    density: HashMap<Material, f64>,
}

// ---------------------------------------------------------------------------
// Reference LAC data from NIST XCOM for interpolation.
// Energies in keV; LAC in cm^-1.
// Each row: (energy_keV, water, lipid, iodine, calcium_hydroxyapatite)
// ---------------------------------------------------------------------------

/// Reference energy-LAC table used for interpolation.
/// Columns: water, lipid (adipose), iodine (pure element basis), calcium (hydroxyapatite).
const LAC_TABLE: [(f64, [f64; 4]); 7] = [
    //  keV    Water    Lipid    Iodine   Calcium(HAp)
    (40.0, [0.2683, 0.2370, 7.580, 1.6100]),
    (60.0, [0.2059, 0.1830, 3.060, 0.7150]),
    (70.0, [0.1928, 0.1722, 1.943, 0.5730]),
    (80.0, [0.1837, 0.1647, 1.310, 0.4580]),
    (100.0, [0.1707, 0.1544, 0.6780, 0.3280]),
    (120.0, [0.1614, 0.1471, 0.5120, 0.2920]),
    (150.0, [0.1494, 0.1420, 0.5470, 0.3000]),
];

/// Intrinsic mass densities (g/cm^3). Canonical home — the water/lipid GLS
/// solver references `DENSITY_WATER`/`DENSITY_LIPID` rather than re-declaring them.
pub(crate) const DENSITY_WATER: f64 = 1.000;
pub(crate) const DENSITY_LIPID: f64 = 0.950;
const DENSITY_IODINE: f64 = 4.930;
const DENSITY_CALCIUM: f64 = 3.180;

/// Linearly interpolate a value from the LAC_TABLE for a given material index
/// at the requested energy. Clamps to the table bounds if energy is outside range.
fn interpolate_lac(energy_kev: f64, material_idx: usize) -> f64 {
    // Clamp to table range
    let first = LAC_TABLE[0];
    let last = LAC_TABLE[LAC_TABLE.len() - 1];

    if energy_kev <= first.0 {
        return first.1[material_idx];
    }
    if energy_kev >= last.0 {
        return last.1[material_idx];
    }

    // Find bracketing interval
    for i in 0..LAC_TABLE.len() - 1 {
        let (e0, vals0) = LAC_TABLE[i];
        let (e1, vals1) = LAC_TABLE[i + 1];
        if energy_kev >= e0 && energy_kev <= e1 {
            let t = (energy_kev - e0) / (e1 - e0);
            return vals0[material_idx] * (1.0 - t) + vals1[material_idx] * t;
        }
    }

    // Should not reach here, but fallback to nearest
    last.1[material_idx]
}

impl MaterialLibrary {
    /// Create a library for the NAEOTOM Alpha MonoPlus VMI+ 70/150 keV pair.
    ///
    /// LAC values are taken directly from the reference table (NIST XCOM
    /// approximations matched to VMI+ effective energies).
    pub fn naeotom_70_150() -> Self {
        let mut lac = HashMap::new();
        lac.insert(Material::Water, [0.1928, 0.1494]);
        lac.insert(Material::Lipid, [0.1722, 0.1420]);
        lac.insert(Material::Iodine, [1.943, 0.5470]);
        lac.insert(Material::Calcium, [0.5730, 0.3000]);

        let mut density = HashMap::new();
        density.insert(Material::Water, DENSITY_WATER);
        density.insert(Material::Lipid, DENSITY_LIPID);
        density.insert(Material::Iodine, DENSITY_IODINE);
        density.insert(Material::Calcium, DENSITY_CALCIUM);

        Self {
            low_kev: 70.0,
            high_kev: 150.0,
            lac,
            density,
        }
    }

    /// Create a library for an arbitrary keV pair.
    ///
    /// Uses linear interpolation between the reference energy points in
    /// [`LAC_TABLE`]. Energies outside the table range are clamped.
    pub fn new(low_kev: f64, high_kev: f64) -> Self {
        let mut lac = HashMap::new();
        for (idx, mat) in Material::ALL.iter().enumerate() {
            let low_lac = interpolate_lac(low_kev, idx);
            let high_lac = interpolate_lac(high_kev, idx);
            lac.insert(*mat, [low_lac, high_lac]);
        }

        let mut density = HashMap::new();
        density.insert(Material::Water, DENSITY_WATER);
        density.insert(Material::Lipid, DENSITY_LIPID);
        density.insert(Material::Iodine, DENSITY_IODINE);
        density.insert(Material::Calcium, DENSITY_CALCIUM);

        Self {
            low_kev,
            high_kev,
            lac,
            density,
        }
    }

    /// Get LAC for a material at the low energy (cm^-1).
    pub fn lac_low(&self, mat: Material) -> f64 {
        self.lac[&mat][0]
    }

    /// Get LAC for a material at the high energy (cm^-1).
    pub fn lac_high(&self, mat: Material) -> f64 {
        self.lac[&mat][1]
    }

    /// Get intrinsic mass density (g/cm^3).
    pub fn density(&self, mat: Material) -> f64 {
        self.density[&mat]
    }

    /// Get intrinsic mass density in mg/mL (= g/cm^3 x 1000).
    pub fn density_mg_ml(&self, mat: Material) -> f64 {
        self.density[&mat] * 1000.0
    }

    /// Enumerate all C(4,3) = 4 material triplets.
    ///
    /// Returns every combination of 3 materials chosen from the 4 basis
    /// materials, yielding exactly 4 triplets.
    pub fn triplets(&self) -> Vec<[Material; 3]> {
        let all = Material::ALL;
        let mut result = Vec::with_capacity(4);
        for i in 0..all.len() {
            for j in (i + 1)..all.len() {
                for k in (j + 1)..all.len() {
                    result.push([all[i], all[j], all[k]]);
                }
            }
        }
        result
    }

    /// Build the 2x3 system matrix A for a given material triplet.
    ///
    /// `A[i][j]` = LAC of material `j` at energy `i`, where:
    /// - row 0 = low energy
    /// - row 1 = high energy
    pub fn system_matrix(&self, triplet: &[Material; 3]) -> [[f64; 3]; 2] {
        [
            [
                self.lac_low(triplet[0]),
                self.lac_low(triplet[1]),
                self.lac_low(triplet[2]),
            ],
            [
                self.lac_high(triplet[0]),
                self.lac_high(triplet[1]),
                self.lac_high(triplet[2]),
            ],
        ]
    }

    /// Convert HU to LAC (cm^-1): mu = (HU / 1000 + 1) * mu_water.
    ///
    /// `energy_idx`: 0 for low energy, 1 for high energy.
    ///
    /// # Panics
    ///
    /// Panics if `energy_idx` is not 0 or 1.
    pub fn hu_to_lac(&self, hu: f64, energy_idx: usize) -> f64 {
        assert!(energy_idx <= 1, "energy_idx must be 0 (low) or 1 (high)");
        let mu_water = self.lac[&Material::Water][energy_idx];
        (hu / 1000.0 + 1.0) * mu_water
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-6;

    #[test]
    fn naeotom_70_150_lac_values() {
        let lib = MaterialLibrary::naeotom_70_150();

        assert_eq!(lib.low_kev, 70.0);
        assert_eq!(lib.high_kev, 150.0);

        // Water
        assert!((lib.lac_low(Material::Water) - 0.1928).abs() < EPS);
        assert!((lib.lac_high(Material::Water) - 0.1494).abs() < EPS);

        // Lipid
        assert!((lib.lac_low(Material::Lipid) - 0.1722).abs() < EPS);
        assert!((lib.lac_high(Material::Lipid) - 0.1420).abs() < EPS);

        // Iodine
        assert!((lib.lac_low(Material::Iodine) - 1.943).abs() < EPS);
        assert!((lib.lac_high(Material::Iodine) - 0.5470).abs() < EPS);

        // Calcium
        assert!((lib.lac_low(Material::Calcium) - 0.5730).abs() < EPS);
        assert!((lib.lac_high(Material::Calcium) - 0.3000).abs() < EPS);
    }

    #[test]
    fn triplets_returns_exactly_four() {
        let lib = MaterialLibrary::naeotom_70_150();
        let trips = lib.triplets();
        assert_eq!(trips.len(), 4);

        // Each triplet should have 3 distinct materials
        for t in &trips {
            assert_ne!(t[0], t[1]);
            assert_ne!(t[1], t[2]);
            assert_ne!(t[0], t[2]);
        }
    }

    #[test]
    fn system_matrix_correct() {
        let lib = MaterialLibrary::naeotom_70_150();
        let triplet = [Material::Water, Material::Iodine, Material::Calcium];
        let a = lib.system_matrix(&triplet);

        // Row 0 = low energy LACs
        assert!((a[0][0] - 0.1928).abs() < EPS); // water @ 70 keV
        assert!((a[0][1] - 1.943).abs() < EPS); // iodine @ 70 keV
        assert!((a[0][2] - 0.5730).abs() < EPS); // calcium @ 70 keV

        // Row 1 = high energy LACs
        assert!((a[1][0] - 0.1494).abs() < EPS); // water @ 150 keV
        assert!((a[1][1] - 0.5470).abs() < EPS); // iodine @ 150 keV
        assert!((a[1][2] - 0.3000).abs() < EPS); // calcium @ 150 keV
    }

    #[test]
    fn hu_to_lac_zero_gives_water() {
        let lib = MaterialLibrary::naeotom_70_150();

        // HU = 0 should give exactly mu_water
        let mu_low = lib.hu_to_lac(0.0, 0);
        assert!((mu_low - 0.1928).abs() < EPS);

        let mu_high = lib.hu_to_lac(0.0, 1);
        assert!((mu_high - 0.1494).abs() < EPS);
    }

    #[test]
    fn hu_to_lac_minus_1000_gives_zero() {
        let lib = MaterialLibrary::naeotom_70_150();

        // HU = -1000 (air) should give 0
        let mu = lib.hu_to_lac(-1000.0, 0);
        assert!(mu.abs() < EPS);
    }

    #[test]
    fn density_mg_ml_correct() {
        let lib = MaterialLibrary::naeotom_70_150();

        assert!((lib.density_mg_ml(Material::Water) - 1000.0).abs() < EPS);
        assert!((lib.density_mg_ml(Material::Lipid) - 950.0).abs() < EPS);
        assert!((lib.density_mg_ml(Material::Iodine) - 4930.0).abs() < EPS);
        assert!((lib.density_mg_ml(Material::Calcium) - 3180.0).abs() < EPS);
    }

    #[test]
    fn density_values() {
        let lib = MaterialLibrary::naeotom_70_150();

        assert!((lib.density(Material::Water) - 1.000).abs() < EPS);
        assert!((lib.density(Material::Lipid) - 0.950).abs() < EPS);
        assert!((lib.density(Material::Iodine) - 4.930).abs() < EPS);
        assert!((lib.density(Material::Calcium) - 3.180).abs() < EPS);
    }

    #[test]
    fn new_at_70_150_matches_naeotom() {
        // Using new() with 70/150 should match the table values exactly
        // since 70 and 150 are exact table entries.
        let lib = MaterialLibrary::new(70.0, 150.0);

        assert!((lib.lac_low(Material::Water) - 0.1928).abs() < EPS);
        assert!((lib.lac_high(Material::Water) - 0.1494).abs() < EPS);
        assert!((lib.lac_low(Material::Iodine) - 1.943).abs() < EPS);
        assert!((lib.lac_high(Material::Iodine) - 0.5470).abs() < EPS);
    }

    #[test]
    fn interpolation_midpoint() {
        // 65 keV is midpoint between 60 and 70 keV table entries.
        // Water: (0.2059 + 0.1928) / 2 = 0.19935
        let lib = MaterialLibrary::new(65.0, 150.0);
        let expected_water_65 = (0.2059 + 0.1928) / 2.0;
        assert!(
            (lib.lac_low(Material::Water) - expected_water_65).abs() < 1e-4,
            "got {}, expected {}",
            lib.lac_low(Material::Water),
            expected_water_65
        );
    }

    #[test]
    fn interpolation_clamped_below() {
        // Energy below table range should clamp to first entry
        let lib = MaterialLibrary::new(20.0, 150.0);
        assert!((lib.lac_low(Material::Water) - 0.2683).abs() < EPS);
    }

    #[test]
    fn interpolation_clamped_above() {
        // Energy above table range should clamp to last entry
        let lib = MaterialLibrary::new(70.0, 200.0);
        assert!((lib.lac_high(Material::Water) - 0.1494).abs() < EPS);
    }

    #[test]
    #[should_panic(expected = "energy_idx must be 0 (low) or 1 (high)")]
    fn hu_to_lac_invalid_index_panics() {
        let lib = MaterialLibrary::naeotom_70_150();
        lib.hu_to_lac(0.0, 2);
    }
}

mod materials;
pub mod result;
pub mod water_lipid;

pub use result::MmdResult;
pub use water_lipid::{
    decompose_slice, decompose_volume_gls, gls_fw, self_calibrate, sigma_cov, WlAnchor,
    WlCalibration,
};

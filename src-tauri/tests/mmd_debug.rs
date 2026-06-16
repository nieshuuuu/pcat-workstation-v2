use std::path::Path;
use pcat_pipeline::dicom_loader::load_dicom_directory;

// Manual diagnostic, not a CI test: it `.unwrap()`s hardcoded absolute paths to
// MonoPlus_{70,100,140,150}keV subdirs. Most patients (e.g. 57955439) only have
// 70/150, so the 100/140 unwraps panic. Run explicitly with `--ignored` when the
// full keV set is present locally.
#[test]
#[ignore = "manual diagnostic; needs all hardcoded MonoPlus keV subdirs present locally"]
fn debug_voxel_values() {
    let base = "/Users/shunie/Developer/PCAT/UCI NAEOTOM CCTA Data/57955439";
    let dirs = [
        format!("{}/MonoPlus_70keV", base),
        format!("{}/MonoPlus_100keV", base),
        format!("{}/MonoPlus_140keV", base),
        format!("{}/MonoPlus_150keV", base),
    ];

    let vols: Vec<_> = dirs.iter().map(|d| {
        load_dicom_directory(Path::new(d)).unwrap()
    }).collect();

    let shape = vols[0].data.raw_dim();
    let (nz, ny, nx) = (shape[0], shape[1], shape[2]);

    // NIST water LAC (cm⁻¹)
    let water_lac = [0.1937, 0.1707, 0.1538, 0.1505];

    // Center voxel
    let (z, y, x) = (nz/2, ny/2, nx/2);
    println!("Center voxel [{z},{y},{x}]:");
    for e in 0..4 {
        let hu = vols[e].data[[z, y, x]];
        let lac = (hu as f64 / 1000.0 + 1.0) * water_lac[e];
        println!("  E{e}: HU={hu:.1}, LAC={lac:.6} cm⁻¹");
    }

    // Find a fat voxel (HU ~ -80 at 70keV)
    'outer_fat: for zz in 0..nz {
        for yy in 0..ny {
            for xx in 0..nx {
                let hu70 = vols[0].data[[zz, yy, xx]];
                if hu70 > -90.0 && hu70 < -70.0 {
                    println!("\nFat voxel [{zz},{yy},{xx}]:");
                    for e in 0..4 {
                        let hu = vols[e].data[[zz, yy, xx]];
                        let lac = (hu as f64 / 1000.0 + 1.0) * water_lac[e];
                        println!("  E{e}: HU={hu:.1}, LAC={lac:.6} cm⁻¹");
                    }
                    break 'outer_fat;
                }
            }
        }
    }

    // Find an enhanced voxel (HU 300-400 at 70keV)
    'outer_enh: for zz in 0..nz {
        for yy in 0..ny {
            for xx in 0..nx {
                let hu70 = vols[0].data[[zz, yy, xx]];
                if hu70 > 300.0 && hu70 < 400.0 {
                    println!("\nEnhanced voxel [{zz},{yy},{xx}]:");
                    for e in 0..4 {
                        let hu = vols[e].data[[zz, yy, xx]];
                        let lac = (hu as f64 / 1000.0 + 1.0) * water_lac[e];
                        println!("  E{e}: HU={hu:.1}, LAC={lac:.6} cm⁻¹");
                    }
                    // Check iodine enhancement pattern (should drop with energy)
                    let hu70 = vols[0].data[[zz, yy, xx]] as f64;
                    let hu150 = vols[3].data[[zz, yy, xx]] as f64;
                    println!("  HU drop 70→150: {:.1} HU (iodine = large drop, bone = small drop)", hu70 - hu150);
                    break 'outer_enh;
                }
            }
        }
    }

    // Print basis material LACs for comparison
    let adipose_lac = [0.1785, 0.1604, 0.1454, 0.1425];
    let iodine_mu_rho = [5.0174, 1.9420, 0.8306, 0.6978];
    let conc = 0.05; // 50 mg/mL
    println!("\n=== Basis material LACs (cm⁻¹) ===");
    println!("  Energy   Water     Adipose   Iodine(50mg/mL)");
    for e in 0..4 {
        let i_lac = water_lac[e] + conc * iodine_mu_rho[e];
        println!("  {}keV:  {:.4}    {:.4}    {:.4}",
            [70, 100, 140, 150][e], water_lac[e], adipose_lac[e], i_lac);
    }
}

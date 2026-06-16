//! End-to-end sanity check for vessel_wall::compute_vessel_geometry.
//!
//! Loads a real DICOM series + the seeds JSON saved by the Tauri app,
//! rebuilds the CPR frame, and prints the measured vessel diameter at
//! uniform arc-length steps along the vessel. Intended to mimic what the
//! cross-section panels show in CprView, without needing the UI running.
//!
//! Run with:
//!
//!   cargo run -p pcat-pipeline --example measure_rca --release
//!
//! Edit the constants at the top of main() to point at a different patient
//! or vessel.

use std::path::Path;

use ndarray::Array3;
use pcat_pipeline::cpr::CprFrame;
use pcat_pipeline::dicom_load::load_series;
use pcat_pipeline::dicom_scan::scan_series;

#[tokio::main(flavor = "multi_thread")]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Override patient by setting `PATIENT=<id>` env var. Defaults to 512143294.
    let patient = std::env::var("PATIENT").unwrap_or_else(|_| "512143294".to_string());
    let dicom_dir_str = format!(
        "/Volumes/Molloilab/Shu Nie/UCI NAEOTOM CCTA Data/{}/MonoPlus_70keV",
        patient
    );
    let sanitize = |s: &str| -> String {
        s.chars()
            .map(|c| {
                if c.is_ascii_alphanumeric() || c == '-' || c == '_' || c == '.' {
                    c
                } else {
                    '_'
                }
            })
            .collect::<String>()
            .replace("..", "_")
    };
    let seed_filename = format!(
        "MonoPlus_70keV__{}.json",
        sanitize(&dicom_dir_str)
    );
    let seeds_path_owned = format!(
        "/Users/shunie/Library/Application Support/com.pcat.workstation/seeds/{}",
        seed_filename
    );
    let dicom_dir = Path::new(&dicom_dir_str);
    let seeds_path = Path::new(&seeds_path_owned);
    let vessel = "RCA";
    let width_mm = 15.0_f64;
    let pixels: usize = 128;
    let step_mm = 2.0_f64;

    // --- 1. Load seeds ---
    println!("seeds:  {}", seeds_path.display());
    let seeds_json = std::fs::read_to_string(seeds_path)?;
    let seeds_raw: serde_json::Value = serde_json::from_str(&seeds_json)?;
    let vessel_data = &seeds_raw["vessels"][vessel];
    let seeds_xyz: Vec<[f64; 3]> = vessel_data["seeds"]
        .as_array()
        .ok_or("no seeds array for this vessel")?
        .iter()
        .map(|p| {
            let a = p.as_array().expect("seed must be a 3-vector");
            [
                a[0].as_f64().expect("x"),
                a[1].as_f64().expect("y"),
                a[2].as_f64().expect("z"),
            ]
        })
        .collect();
    let ostium_fraction = vessel_data["ostiumFraction"].as_f64();
    println!(
        "{vessel}: {} seeds, ostium fraction = {:?}",
        seeds_xyz.len(),
        ostium_fraction
    );
    if seeds_xyz.len() < 2 {
        return Err("need at least 2 seeds to build a centerline".into());
    }

    // Seed positions in the JSON are patient-space [x, y, z] mm (the order the
    // frontend stores them). CprFrame expects [z, y, x] mm — matches what
    // CprView.svelte reorders before sending to the Rust backend.
    let seeds_zyx: Vec<[f64; 3]> = seeds_xyz.iter().map(|&[x, y, z]| [z, y, x]).collect();

    // --- 2. Load DICOM ---
    println!("dicom:  {}", dicom_dir.display());
    let descriptors = scan_series(dicom_dir).await?;
    let desc = descriptors
        .into_iter()
        .next()
        .ok_or("no DICOM series found in folder")?;
    println!(
        "series: {} ({} slices, {}x{})",
        desc.description, desc.num_slices, desc.rows, desc.cols
    );
    let uid = desc.uid.clone();
    let vol = load_series(dicom_dir, &uid, None).await?;

    let nz = vol.metadata.num_slices;
    let ny = vol.metadata.rows as usize;
    let nx = vol.metadata.cols as usize;
    let data_f32: Vec<f32> = vol.voxels_i16.iter().map(|&v| v as f32).collect();
    let volume_arr = Array3::from_shape_vec((nz, ny, nx), data_f32)?;

    let spacing = [
        vol.metadata.slice_spacing,
        vol.metadata.pixel_spacing[0],
        vol.metadata.pixel_spacing[1],
    ];
    let ipp = vol.metadata.image_position_patient;
    let origin = [
        vol.metadata
            .slice_positions_z
            .first()
            .copied()
            .unwrap_or(ipp[2]),
        ipp[1],
        ipp[0],
    ];

    let iop = vol.metadata.orientation;
    let iop_row = [iop[0], iop[1], iop[2]];
    let iop_col = [iop[3], iop[4], iop[5]];
    let normal = [
        iop_row[1] * iop_col[2] - iop_row[2] * iop_col[1],
        iop_row[2] * iop_col[0] - iop_row[0] * iop_col[2],
        iop_row[0] * iop_col[1] - iop_row[1] * iop_col[0],
    ];
    let direction = [
        iop_row[0], iop_row[1], iop_row[2], iop_col[0], iop_col[1], iop_col[2], normal[0],
        normal[1], normal[2],
    ];

    // --- 3. Build CPR frame ---
    let n_cols = 256;
    let frame = CprFrame::from_centerline(&seeds_zyx, n_cols);
    let total_arc = *frame.arclengths.last().expect("arclengths non-empty");
    println!("total arc-length: {:.1} mm", total_arc);
    if let Some(frac) = ostium_fraction {
        println!(
            "ostium at arc {:.1} mm (fraction {:.3})",
            frac * total_arc,
            frac
        );
    }

    // --- 4. Walk the centerline and measure diameter at each step ---
    let mm_per_px = 2.0 * width_mm / pixels as f64;
    let ostium_arc = ostium_fraction.map(|f| f * total_arc);
    let pixels_f = pixels as f64;

    println!();
    println!("arc_mm  arc_rel  D_mm   r_min  r_max  r_std  center_HU");
    println!("─────────────────────────────────────────────────────────");

    let mut arc = 0.0;
    while arc <= total_arc {
        let position_frac = (arc / total_arc).clamp(0.0, 1.0);
        let cs = frame.render_cross_section(
            &volume_arr,
            spacing,
            origin,
            &direction,
            position_frac,
            0.0,
            width_mm,
            pixels,
        );

        let centre_px = pixels / 2;
        let centre_hu = cs.image[centre_px * pixels + centre_px];

        let radii_mm: Vec<f64> = cs
            .vessel_wall
            .iter()
            .map(|&[x, y]| {
                let dx = x - pixels_f / 2.0;
                let dy = y - pixels_f / 2.0;
                (dx * dx + dy * dy).sqrt() * mm_per_px
            })
            .collect();
        let r_min = radii_mm.iter().cloned().fold(f64::INFINITY, f64::min);
        let r_max = radii_mm.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let r_mean = radii_mm.iter().sum::<f64>() / radii_mm.len() as f64;
        let r_var = radii_mm
            .iter()
            .map(|r| (r - r_mean).powi(2))
            .sum::<f64>()
            / radii_mm.len() as f64;
        let r_std = r_var.sqrt();

        let arc_rel = ostium_arc.map(|o| arc - o);
        let arc_rel_str = match arc_rel {
            Some(v) => format!("{:+6.1}", v),
            None => "   n/a".to_string(),
        };

        println!(
            "{:5.1}  {}  {:5.2}  {:5.2}  {:5.2}  {:5.2}  {:7.1}",
            arc, arc_rel_str, cs.vessel_diameter_mm, r_min, r_max, r_std, centre_hu,
        );
        arc += step_mm;
    }

    // --- 5. HU profile dump at a single arc position for eyeballing ---
    let probe_arc = 30.0_f64;
    let frac = (probe_arc / total_arc).clamp(0.0, 1.0);
    let cs = frame.render_cross_section(
        &volume_arr,
        spacing,
        origin,
        &direction,
        frac,
        0.0,
        width_mm,
        pixels,
    );

    println!();
    println!(
        "── HU profiles at arc {:.1} mm (D_algo = {:.2} mm) ──",
        probe_arc, cs.vessel_diameter_mm
    );
    println!("Columns: radius_mm, then HU along 0°/90°/180°/270° rays from image centre.");
    let centre = pixels as f64 / 2.0;
    let step_px = 0.5_f64;
    let n_radial = (centre / step_px) as usize;
    let dirs = [(1.0_f64, 0.0), (0.0, -1.0), (-1.0, 0.0), (0.0, 1.0)];

    // Collect per-ray profiles.
    let mut profiles: Vec<Vec<f32>> = dirs
        .iter()
        .map(|(cx, cy)| {
            (0..n_radial)
                .map(|ri| {
                    let r = ri as f64 * step_px;
                    let x = centre + r * cx;
                    let y = centre + r * cy;
                    let max_coord = (pixels as f64) - 1.0;
                    if x < 0.0 || y < 0.0 || x > max_coord || y > max_coord {
                        return f32::NAN;
                    }
                    let x0 = x.floor() as usize;
                    let y0 = y.floor() as usize;
                    let x1 = (x0 + 1).min(pixels - 1);
                    let y1 = (y0 + 1).min(pixels - 1);
                    let fx = (x - x0 as f64) as f32;
                    let fy = (y - y0 as f64) as f32;
                    let v00 = cs.image[y0 * pixels + x0];
                    let v10 = cs.image[y0 * pixels + x1];
                    let v01 = cs.image[y1 * pixels + x0];
                    let v11 = cs.image[y1 * pixels + x1];
                    v00 * (1.0 - fx) * (1.0 - fy)
                        + v10 * fx * (1.0 - fy)
                        + v01 * (1.0 - fx) * fy
                        + v11 * fx * fy
                })
                .collect()
        })
        .collect();

    // Print the first ~20 mm of each ray (covers any plausible coronary radius).
    let mm_per_px_probe = 2.0 * width_mm / pixels as f64;
    let max_r_mm = 8.0;
    let max_ri = (max_r_mm / (step_px * mm_per_px_probe)) as usize;
    println!("  r[mm]   east    north    west    south");
    for ri in 0..max_ri.min(n_radial) {
        let r_mm = ri as f64 * step_px * mm_per_px_probe;
        print!("  {:5.2}  ", r_mm);
        for p in &mut profiles {
            print!("{:7.1}  ", p[ri]);
        }
        println!();
    }

    Ok(())
}

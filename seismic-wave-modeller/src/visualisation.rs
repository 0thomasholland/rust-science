use ndarray::Array2;
use plotters::prelude::*;
use std::path::Path;
use colorcet::ColorMap;

pub struct WavefieldVisualiser {
    output_dir: String,
    width: u32,
    height: u32,
}

impl WavefieldVisualiser {
    pub fn new(output_dir: &str, width: u32, height: u32) -> Self {
        // Create output directory
        std::fs::create_dir_all(output_dir).unwrap();

        Self {
            output_dir: output_dir.to_string(),
            width,
            height,
        }
    }

    pub fn plot_field(
        &self,
        data: &Array2<f64>,
        timestep: usize,
        field_name: &str,
        time: f64,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let filename = format!("{}/{}_{:06}.png", self.output_dir, field_name, timestep);
        let root = BitMapBackend::new(&filename, (self.width, self.height)).into_drawing_area();

        root.fill(&WHITE)?;

        let (nx, nz) = data.dim();

        // Find min and max for color scale (symmetric around zero)
        let max_abs = data.iter().map(|&v| v.abs()).fold(0.0_f64, f64::max);

        let min_val = -max_abs;
        let max_val = max_abs;

        // Create title
        let title = format!("{} at t={:.4}s (step {})", field_name, time, timestep);

        let mut chart = ChartBuilder::on(&root)
            .caption(&title, ("sans-serif", 30))
            .margin(10)
            .x_label_area_size(40)
            .y_label_area_size(40)
            .build_cartesian_2d(0..nx, 0..nz)?;

        chart
            .configure_mesh()
            .x_desc("X (grid points)")
            .y_desc("Z (grid points)")
            .draw()?;

        // Build CET colormap once per frame.
        // Use "cet_l8" which corresponds to CET-L08.
        let cmap: ColorMap = "cet_l8".parse().expect("Failed to parse 'cet_l8' colormap. Check colorcet version/features.");
        let cmap_rgb: Vec<[u8; 3]> = cmap.get_rgb_int::<u8>();

        // Draw the field as colored rectangles
        for i in 0..nx {
            for k in 0..nz {
                let value = data[[i, k]];

                // Map value to color using the colorcet colormap
                let color = value_to_color_from_cmap(value, min_val, max_val, &cmap_rgb);

                chart.draw_series(std::iter::once(Rectangle::new(
                    [(i, k), (i + 1, k + 1)],
                    color.filled(),
                )))?;
            }
        }

        root.present()?;
        println!("Saved frame: {}", filename);

        Ok(())
    }
}

fn value_to_color_from_cmap(value: f64, min_val: f64, max_val: f64, cmap: &[[u8; 3]]) -> RGBColor {
    let normalized = if max_val > min_val {
        (value - min_val) / (max_val - min_val)
    } else {
        0.5
    }
    .clamp(0.0, 1.0);

    let idx = ((normalized * ((cmap.len() - 1) as f64)).round() as usize).min(cmap.len() - 1);
    let rgb = cmap[idx];
    RGBColor(rgb[0], rgb[1], rgb[2])
}
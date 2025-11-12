use ndarray::Array2;
use plotters::prelude::*;
use std::path::Path;

pub struct WavefieldVisualiser {
    output_dir: String,
    width: u32,
    height: u32,
    dx: f64,
    dz: f64,
    // Store as a boxed trait object
    gradient: Box<dyn colorgrad::Gradient>,
}

impl WavefieldVisualiser {
    pub fn new(output_dir: &str, width: u32, height: u32, dx: f64, dz: f64) -> Self {
        std::fs::create_dir_all(output_dir).unwrap();

        // Box the gradient when creating it
        let gradient = Box::new(colorgrad::preset::pu_or());

        Self {
            output_dir: output_dir.to_string(),
            width,
            height,
            dx,
            dz,
            gradient,
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
        let max_abs = data.iter().map(|&v| v.abs()).fold(0.0_f64, f64::max);

        // Clip to a fraction of max to make waves visible
        let clip_factor = 0.5; // Adjust this value between 0.01 and 0.5
        let min_val = -max_abs * clip_factor;
        let max_val = max_abs * clip_factor;

        // Calculate physical dimensions
        let x_max = (nx as f64) * self.dx;
        let z_max = (nz as f64) * self.dz;

        let title = format!("{} at t={:.4}s (step {})", field_name, time, timestep);
        let mut chart = ChartBuilder::on(&root)
            .caption(&title, ("sans-serif", 30))
            .margin(10)
            .x_label_area_size(40)
            .y_label_area_size(40)
            .build_cartesian_2d(0.0..x_max, 0.0..z_max)?;

        chart
            .configure_mesh()
            .x_desc("X (meters)")
            .y_desc("Z (meters)")
            .draw()?;

        for i in 0..nx {
            for k in 0..nz {
                let value = data[[i, k]];
                let color = self.value_to_color(value, min_val, max_val);

                // Convert grid indices to physical coordinates
                let x1 = (i as f64) * self.dx;
                let x2 = ((i + 1) as f64) * self.dx;
                let z1 = (k as f64) * self.dz;
                let z2 = ((k + 1) as f64) * self.dz;

                chart.draw_series(std::iter::once(Rectangle::new(
                    [(x1, z1), (x2, z2)],
                    color.filled(),
                )))?;
            }
        }

        root.present()?;
        println!("Saved frame: {}", filename);
        Ok(())
    }

    fn value_to_color(&self, value: f64, min_val: f64, max_val: f64) -> RGBColor {
        let normalized = if max_val > min_val {
            // Simple linear scaling with clamping
            (value - min_val) / (max_val - min_val)
        } else {
            0.5
        };
        let normalized = normalized.clamp(0.0, 1.0);
        let color_rgba = self.gradient.at(normalized as f32).to_rgba8();
        RGBColor(color_rgba[0], color_rgba[1], color_rgba[2])
    }
    pub fn plot_p_and_s_overlay(
        &self,
        p_wave: &Array2<f64>,
        s_wave: &Array2<f64>,
        timestep: usize,
        time: f64,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let filename = format!("{}/p_s_overlay_{:06}.png", self.output_dir, timestep);
        let root = BitMapBackend::new(&filename, (self.width, self.height)).into_drawing_area();
        root.fill(&WHITE)?;

        let (nx, nz) = p_wave.dim();

        // Calculate percentiles for each wave type independently
        let mut p_abs: Vec<f64> = p_wave.iter().map(|v| v.abs()).collect();
        p_abs.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let p_max = p_abs[((p_abs.len() as f64) * 0.98) as usize];

        let mut s_abs: Vec<f64> = s_wave.iter().map(|v| v.abs()).collect();
        s_abs.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let s_max = s_abs[((s_abs.len() as f64) * 0.98) as usize];

        let x_max = (nx as f64) * self.dx;
        let z_max = (nz as f64) * self.dz;

        let title = format!("P-waves (Red) & S-waves (Blue) at t={:.4}s", time);
        let mut chart = ChartBuilder::on(&root)
            .caption(&title, ("sans-serif", 30))
            .margin(10)
            .x_label_area_size(40)
            .y_label_area_size(40)
            .build_cartesian_2d(0.0..x_max, 0.0..z_max)?;

        chart
            .configure_mesh()
            .x_desc("X (meters)")
            .y_desc("Z (meters)")
            .draw()?;

        for i in 0..nx {
            for k in 0..nz {
                let p_val = p_wave[[i, k]];
                let s_val = s_wave[[i, k]];

                // Normalize each independently
                let p_norm = (p_val / p_max).clamp(-1.0, 1.0);
                let s_norm = (s_val / s_max).clamp(-1.0, 1.0);

                // Map to RGB: P-waves in red channel, S-waves in blue channel
                let r = ((p_norm.abs() * 255.0) as u8);
                let b = ((s_norm.abs() * 255.0) as u8);
                let g = 0; // Keep green at 0

                let color = RGBColor(r, g, b);

                let x1 = (i as f64) * self.dx;
                let x2 = ((i + 1) as f64) * self.dx;
                let z1 = (k as f64) * self.dz;
                let z2 = ((k + 1) as f64) * self.dz;

                chart.draw_series(std::iter::once(Rectangle::new(
                    [(x1, z1), (x2, z2)],
                    color.filled(),
                )))?;
            }
        }

        root.present()?;
        println!("Saved P-S overlay frame: {}", filename);
        Ok(())
    }
}

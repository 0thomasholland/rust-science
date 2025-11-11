use ndarray::Array2;
use plotters::prelude::*;
use std::path::Path;

pub struct WavefieldVisualiser {
    output_dir: String,
    width: u32,
    height: u32,
    // Store as a boxed trait object
    gradient: Box<dyn colorgrad::Gradient>,
}

impl WavefieldVisualiser {
    pub fn new(output_dir: &str, width: u32, height: u32) -> Self {
        std::fs::create_dir_all(output_dir).unwrap();

        // Box the gradient when creating it
        let gradient = Box::new(colorgrad::preset::rd_yl_bu()); 

        Self {
            output_dir: output_dir.to_string(),
            width,
            height,
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
        let min_val = -max_abs;
        let max_val = max_abs;

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

        for i in 0..nx {
            for k in 0..nz {
                let value = data[[i, k]];
                let color = self.value_to_color(value, min_val, max_val);
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

    fn value_to_color(&self, value: f64, min_val: f64, max_val: f64) -> RGBColor {
        let normalized = if max_val > min_val {
            (value - min_val) / (max_val - min_val)
        } else {
            0.5
        };
        let normalized = normalized.clamp(0.0, 1.0);
        let color_rgba = self.gradient.at(normalized as f32).to_rgba8();
        RGBColor(color_rgba[0], color_rgba[1], color_rgba[2])
    }
}
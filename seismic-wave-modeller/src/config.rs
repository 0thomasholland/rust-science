use anyhow::{anyhow, Result};
use serde::{Deserialize, Serialize};
use std::fs;

/// Grid configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GridConfig {
    pub nx: usize,
    pub nz: usize,
    pub dx: f64,
    pub dz: f64,
}

impl GridConfig {
    fn validate(&self) -> Result<()> {
        if self.nx == 0 || self.nz == 0 {
            return Err(anyhow!("Grid dimensions must be positive (nx={}, nz={})", self.nx, self.nz));
        }
        if self.dx <= 0.0 || self.dz <= 0.0 {
            return Err(anyhow!(
                "Grid spacing must be positive (dx={}, dz={})",
                self.dx,
                self.dz
            ));
        }
        Ok(())
    }
}

/// Homogeneous material properties
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MaterialConfig {
    pub vp: f64,   // P-wave velocity (m/s)
    pub vs: f64,   // S-wave velocity (m/s)
    pub rho: f64,  // Density (kg/m³)
}

impl MaterialConfig {
    fn validate(&self) -> Result<()> {
        if self.vp <= 0.0 || self.vs <= 0.0 || self.rho <= 0.0 {
            return Err(anyhow!(
                "Material properties must be positive (vp={}, vs={}, rho={})",
                self.vp,
                self.vs,
                self.rho
            ));
        }
        if self.vs > self.vp {
            return Err(anyhow!(
                "S-wave velocity must be less than P-wave velocity (vs={} > vp={})",
                self.vs,
                self.vp
            ));
        }
        Ok(())
    }
}

/// Simulation configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SimulationConfig {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub dt: Option<f64>,  // Optional: will be auto-computed from CFL if not provided
    pub total_time: f64,
    #[serde(default = "default_report_period")]
    pub report_period: usize,
    #[serde(default = "default_cfl_safety")]
    pub cfl_safety: f64,
}

fn default_report_period() -> usize {
    100
}

fn default_cfl_safety() -> f64 {
    0.5
}

impl SimulationConfig {
    fn validate(&self) -> Result<()> {
        if self.total_time <= 0.0 {
            return Err(anyhow!("total_time must be positive, got {}", self.total_time));
        }
        if self.cfl_safety <= 0.0 || self.cfl_safety > 1.0 {
            return Err(anyhow!(
                "cfl_safety must be in (0, 1], got {}",
                self.cfl_safety
            ));
        }
        Ok(())
    }

    /// Compute dt from CFL condition if not specified
    pub fn compute_dt_if_needed(&mut self, dx: f64, dz: f64, vp_max: f64) {
        if self.dt.is_none() {
            let min_spacing = dx.min(dz);
            self.dt = Some(self.cfl_safety * min_spacing / vp_max);
        }
    }

    /// Calculate number of timesteps given dt
    pub fn compute_nt(&self, dt: f64) -> usize {
        (self.total_time / dt).round() as usize
    }
}

/// Source configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SourceConfig {
    pub x: usize,
    pub z: usize,
    pub amplitude: f64,
    pub frequency: f64,
    #[serde(default)]
    pub trigger_time: f64,  // Default: 0.0
}

impl SourceConfig {
    fn validate(&self, nx: usize, nz: usize) -> Result<()> {
        if self.x >= nx || self.z >= nz {
            return Err(anyhow!(
                "Source position ({}, {}) is outside grid bounds ({}, {})",
                self.x,
                self.z,
                nx,
                nz
            ));
        }
        if self.amplitude <= 0.0 {
            return Err(anyhow!("Source amplitude must be positive, got {}", self.amplitude));
        }
        if self.frequency <= 0.0 {
            return Err(anyhow!("Source frequency must be positive, got {}", self.frequency));
        }
        if self.trigger_time < 0.0 {
            return Err(anyhow!(
                "Source trigger_time must be non-negative, got {}",
                self.trigger_time
            ));
        }
        Ok(())
    }
}

/// Visualization configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VisualizationConfig {
    #[serde(default = "default_video_length")]
    pub video_length: f64,
    #[serde(default = "default_fps")]
    pub fps: f64,
    #[serde(default = "default_field")]
    pub field: String,
    #[serde(default)]
    pub iteration: i64,
    #[serde(default = "default_image_width")]
    pub image_width: usize,
    #[serde(default = "default_image_height")]
    pub image_height: usize,
}

fn default_video_length() -> f64 {
    10.0
}

fn default_fps() -> f64 {
    30.0
}

fn default_field() -> String {
    "vmag".to_string()
}

fn default_image_width() -> usize {
    1200
}

fn default_image_height() -> usize {
    1000
}

impl VisualizationConfig {
    fn validate(&self) -> Result<()> {
        let valid_fields = [
            "vx", "vz", "vmag", "divergence", "curl", "sigma_xx", "sigma_zz", "sigma_xz",
        ];
        if !valid_fields.contains(&self.field.as_str()) {
            return Err(anyhow!(
                "Invalid field '{}'. Must be one of: {:?}",
                self.field,
                valid_fields
            ));
        }
        if self.video_length <= 0.0 {
            return Err(anyhow!("video_length must be positive, got {}", self.video_length));
        }
        if self.fps <= 0.0 {
            return Err(anyhow!("fps must be positive, got {}", self.fps));
        }
        if self.image_width == 0 || self.image_height == 0 {
            return Err(anyhow!(
                "Image dimensions must be positive (width={}, height={})",
                self.image_width,
                self.image_height
            ));
        }
        Ok(())
    }
}

/// Complete simulation configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Config {
    pub grid: GridConfig,
    pub materials: MaterialConfig,
    pub simulation: SimulationConfig,
    pub sources: Vec<SourceConfig>,
    pub visualization: VisualizationConfig,
}

impl Config {
    /// Load configuration from TOML file
    pub fn from_file(path: &str) -> Result<Self> {
        let content = fs::read_to_string(path)
            .map_err(|e| anyhow!("Failed to read config file '{}': {}", path, e))?;

        let mut config: Config = toml::from_str(&content)
            .map_err(|e| anyhow!("Failed to parse TOML config: {}", e))?;

        // Validate before returning
        config.validate()?;

        Ok(config)
    }

    /// Validate all configuration parameters
    pub fn validate(&mut self) -> Result<()> {
        self.grid.validate()?;
        self.materials.validate()?;
        self.simulation.validate()?;
        self.visualization.validate()?;

        // Compute dt from CFL if needed
        self.simulation
            .compute_dt_if_needed(self.grid.dx, self.grid.dz, self.materials.vp);

        // Validate sources against grid
        for source in &self.sources {
            source.validate(self.grid.nx, self.grid.nz)?;
        }

        if self.sources.is_empty() {
            return Err(anyhow!("At least one source must be defined"));
        }

        // Additional validation: check dt is not too small
        let dt = self.simulation.dt.unwrap();
        if dt < 1e-7 {
            eprintln!(
                "Warning: Computed dt is very small ({}), simulation may be slow",
                dt
            );
        }

        Ok(())
    }

    /// Print configuration summary
    pub fn print_summary(&self) {
        println!("=== Simulation Configuration ===");
        println!("Grid: {}x{} ({} x {} m)", self.grid.nx, self.grid.nz,
                 self.grid.nx as f64 * self.grid.dx,
                 self.grid.nz as f64 * self.grid.dz);
        println!(
            "Materials: Vp={} m/s, Vs={} m/s, ρ={} kg/m³",
            self.materials.vp, self.materials.vs, self.materials.rho
        );
        let dt = self.simulation.dt.unwrap();
        let nt = self.simulation.compute_nt(dt);
        println!(
            "Simulation: dt={} s, nt={}, total_time={} s",
            dt, nt, self.simulation.total_time
        );
        println!("Sources: {} source(s)", self.sources.len());
        for (i, src) in self.sources.iter().enumerate() {
            println!(
                "  Source {}: position ({}, {}), freq={} Hz, amp={}",
                i, src.x, src.z, src.frequency, src.amplitude
            );
            if src.trigger_time > 0.0 {
                println!("    Trigger delay: {} s", src.trigger_time);
            }
        }
        println!(
            "Visualization: {} at {} FPS for {} s, field={}",
            self.visualization.image_width, self.visualization.fps, self.visualization.video_length, self.visualization.field
        );
        println!("================================");
    }
}

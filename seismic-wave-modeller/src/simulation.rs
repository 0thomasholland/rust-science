use crate::grid::Grid;
use crate::materials::MaterialProperties;
use crate::visualisation::WavefieldVisualiser;
use crate::wavefield::Wavefield;
use ndarray::Zip;
use rayon::prelude::*;
use std::f64::consts::PI;

#[derive(Clone, Debug)]
pub struct Source {
    pub x: usize,
    pub z: usize,
    pub amplitude: f64,
    pub frequency: f64,
    triggered: bool,           // Add this
    trigger_time: Option<f64>, // Add this - when to trigger (in seconds)
}

impl Source {
    pub fn new(x: usize, z: usize, amplitude: f64, frequency: f64) -> Self {
        Self {
            x,
            z,
            amplitude,
            frequency,
            triggered: false,
            trigger_time: Some(0.0), // Trigger immediately by default
        }
    }

    pub fn with_delay(x: usize, z: usize, amplitude: f64, frequency: f64, delay: f64) -> Self {
        Self {
            x,
            z,
            amplitude,
            frequency,
            triggered: false,
            trigger_time: Some(delay),
        }
    }

    pub fn should_trigger(&self, current_time: f64) -> bool {
        !self.triggered && self.trigger_time.map_or(false, |t| current_time >= t)
    }

    pub fn mark_triggered(&mut self) {
        self.triggered = true;
    }
    pub fn ricker_wavelet(&self, t: f64) -> f64 {
        // Normalized Ricker wavelet for this source
        // Peak amplitude is frequency-independent when properly normalized
        let tau = t - self.trigger_time.unwrap_or(0.0);
        let arg = (PI * self.frequency * tau).powi(2);
        let wavelet = (1.0 - 2.0 * arg) * (-arg).exp();
        // Normalize by peak amplitude to maintain constant energy across frequencies
        // Peak of (1 - 2*arg) * exp(-arg) occurs at arg = 0.5, with value ~ 0.3677
        let normalization = 0.75 * PI * self.frequency * self.frequency;
        wavelet * normalization
    }
}

pub struct SimulationParams {
    pub dt: f64,              // Time step (seconds)
    pub nt: usize,            // Number of time steps
    pub report_period: usize, // How many times to report progress

    // Source parameters
    pub sources: Vec<Source>, // Time delay for source (seconds)
    pub cfl_safety: f64,      // CFL safety factor
}

impl SimulationParams {
    pub fn new(
        dt: f64,
        nt: usize,
        report_period: usize,
        source_i: usize,
        source_k: usize,
        source_f0: f64,
        amplitude: f64,
        cfl_safety: f64,
    ) -> Self {
        // Calculate appropriate t0 from f0
        // Typically t0 = 1.2 / f0 so the wavelet starts near zero
        let source_t0 = 1.2 / source_f0;

        Self {
            dt,
            nt,
            report_period,
            sources: vec![Source::with_delay(
                source_i, source_k, source_f0, amplitude, source_t0,
            )],
            cfl_safety,
        }
    }

    pub fn total_time(&self) -> f64 {
        // Total simulation time from number of time steps and time step size
        // nt * dt
        self.nt as f64 * self.dt
    }

    pub fn check_cfl(&self, dx: f64, dz: f64, vp_max: f64) -> bool {
        // CFL condition for 2D elastic waves with staggered grid:
        // dt <= CFL_safety × min(dx, dz) / vp_max
        // where CFL_safety is typically 0.5
        // returns true is safe, false if unstable

        let min_spacing = dx.min(dz);
        let dt_max = self.cfl_safety * min_spacing / vp_max;

        self.dt <= dt_max
    }

    pub fn compute_stable_dt(&self, dx: f64, dz: f64, vp_max: f64) -> f64 {
        // Helper function to compute a stable time step
        // Return CFL_safety × min(dx, dz) / vp_max

        self.cfl_safety * dx.min(dz) / vp_max
    }
}

pub struct Simulation {
    pub grid: Grid,
    pub materials: MaterialProperties,
    pub wavefield_current: Wavefield,
    pub wavefield_old: Wavefield,
    pub params: SimulationParams,
    current_timestep: usize,
}

impl Simulation {
    pub fn new(grid: Grid, materials: MaterialProperties, params: SimulationParams) -> Self {
        // Validate that the time step satisfies CFL

        let vp_max = Self::find_max_vp(&materials);
        if !params.check_cfl(grid.dx, grid.dz, vp_max) {
            panic!(
                "CFL condition violated! dt={}, max stable dt={}",
                params.dt,
                params.compute_stable_dt(grid.dx, grid.dz, vp_max)
            );
        }
        for (idx, source) in params.sources.iter().enumerate() {
            if source.x >= grid.nx || source.z >= grid.nz {
                // Changed from i,k to x,z
                panic!(
                    "Source {} at ({}, {}) is outside grid bounds ({}, {})",
                    idx,
                    source.x,
                    source.z,
                    grid.nx,
                    grid.nz // Changed from i,k to x,z
                );
            }
        }
        // Create wavefields
        let wavefield_current = Wavefield::new(grid.nx, grid.nz);
        let wavefield_old = Wavefield::new(grid.nx, grid.nz);

        Self {
            grid,
            materials,
            wavefield_current,
            wavefield_old,
            params,
            current_timestep: 0,
        }
    }

    pub fn current_time(&self) -> f64 {
        // Return current time based on timestep
        self.current_timestep as f64 * self.params.dt
    }

    pub fn is_finished(&self) -> bool {
        // Check if simulation has run all timesteps
        // current_timestep >= nt
        self.current_timestep >= self.params.nt
    }
    fn find_max_vp(materials: &MaterialProperties) -> f64 {
        // Find maximum P-wave velocity in the domain
        // Vp = sqrt((lambda + 2*mu) / rho)

        let (nx, nz) = materials.lambda.dim();
        let mut vp_max = 0.0;

        // TODO: Loop over all grid points and find max Vp
        // For each (i, k):
        //   vp = sqrt((lambda[[i,k]] + 2.0 * mu[[i,k]]) / rho[[i,k]])
        //   vp_max = vp_max.max(vp)
        for i in 0..nx {
            for k in 0..nz {
                let lambda = materials.lambda[[i, k]];
                let mu = materials.mu[[i, k]];
                let rho = materials.rho[[i, k]];
                let vp = ((lambda + 2.0 * mu) / rho).sqrt();
                if vp > vp_max {
                    vp_max = vp;
                }
            }
        }
        vp_max
    }

    fn update_stresses(&mut self) {
        let dx = self.grid.dx;
        let dz = self.grid.dz;
        let dt = self.params.dt;
        let nx = self.grid.nx;
        let nz = self.grid.nz;

        // Update sigma_xx and sigma_zz (normal stresses)
        // These live at (i, k) positions
        // Equations:
        // ∂σxx/∂t = (λ + 2μ) × ∂vx/∂x + λ × ∂vz/∂z
        // ∂σzz/∂t = λ × ∂vx/∂x + (λ + 2μ) × ∂vz/∂z

        for i in 1..nx - 1 {
            for k in 1..nz - 1 {
                // Get material properties at this point
                let lambda = self.materials.lambda[[i, k]];
                let mu = self.materials.mu[[i, k]];
                let lambda_plus_2mu = lambda + 2.0 * mu;

                // Compute velocity derivatives using staggered grid
                // ∂vx/∂x ≈ (vx[i,k] - vx[i-1,k]) / dx
                let dvx_dx = (self.wavefield_current.vx[[i, k]]
                    - self.wavefield_current.vx[[i - 1, k]])
                    / dx;

                // vz lives at (i, k+1/2), so for σxx at (i,k):
                // ∂vz/∂z ≈ (vz[i,k] - vz[i,k-1]) / dz
                let dvz_dz = (self.wavefield_current.vz[[i, k]]
                    - self.wavefield_current.vz[[i, k - 1]])
                    / dz;

                // Update stresses
                // sigma_xx_new = sigma_xx_old + dt * [(λ+2μ) * dvx_dx + λ * dvz_dz]
                // sigma_zz_new = sigma_zz_old + dt * [λ * dvx_dx + (λ+2μ) * dvz_dz]
                let sigma_xx_new = self.wavefield_current.sigma_xx[[i, k]]
                    + dt * (lambda_plus_2mu * dvx_dx + lambda * dvz_dz);
                let sigma_zz_new = self.wavefield_current.sigma_zz[[i, k]]
                    + dt * (lambda * dvx_dx + lambda_plus_2mu * dvz_dz);
                self.wavefield_current.sigma_xx[[i, k]] = sigma_xx_new;
                self.wavefield_current.sigma_zz[[i, k]] = sigma_zz_new
            }
        }

        // Update sigma_xz (shear stress)
        // This lives at (i+1/2, k+1/2) positions
        // Equation: ∂σxz/∂t = μ × (∂vx/∂z + ∂vz/∂x)

        for i in 0..nx - 1 {
            for k in 0..nz - 1 {
                // Get mu at shear stress position (already pre-computed)
                let mu = self.materials.mu_xz[[i, k]];

                // Compute velocity derivatives
                // vx lives at (i+1/2, k), we need ∂vx/∂z at (i+1/2, k+1/2)
                // ∂vx/∂z ≈ (vx[i,k+1] - vx[i,k]) / dz
                let dvx_dz = (self.wavefield_current.vx[[i, k + 1]]
                    - self.wavefield_current.vx[[i, k]])
                    / dz;

                // vz lives at (i, k+1/2), we need ∂vz/∂x at (i+1/2, k+1/2)
                // ∂vz/∂x ≈ (vz[i+1,k] - vz[i,k]) / dx
                let dvz_dx = (self.wavefield_current.vz[[i + 1, k]]
                    - self.wavefield_current.vz[[i, k]])
                    / dx;

                // Update shear stress
                // sigma_xz_new = sigma_xz_old + dt * μ * (dvx_dz + dvz_dx)
                let sigma_xz_new =
                    self.wavefield_current.sigma_xz[[i, k]] + dt * mu * (dvx_dz + dvz_dx);
                self.wavefield_current.sigma_xz[[i, k]] = sigma_xz_new;
            }
        }
    }

    fn update_velocities(&mut self) {
        let dx = self.grid.dx;
        let dz = self.grid.dz;
        let dt = self.params.dt;
        let nx = self.grid.nx;
        let nz = self.grid.nz;

        // Update vx (horizontal velocity)
        // Lives at (i+1/2, k) positions
        // Equation: ∂vx/∂t = (1/ρ) × (∂σxx/∂x + ∂σxz/∂z)

        for i in 1..nx - 1 {
            for k in 1..nz - 1 {
                // Get density at vx position (pre-computed)
                let rho = self.materials.rho_vx[[i, k]];

                // Compute stress derivatives
                // σxx lives at (i, k), we need ∂σxx/∂x at (i+1/2, k)
                // ∂σxx/∂x ≈ (σxx[i+1,k] - σxx[i,k]) / dx
                let dsigma_xx_dx = (self.wavefield_current.sigma_xx[[i + 1, k]]
                    - self.wavefield_current.sigma_xx[[i, k]])
                    / dx;

                // σxz lives at (i+1/2, k+1/2), we need ∂σxz/∂z at (i+1/2, k)
                // ∂σxz/∂z ≈ (σxz[i,k] - σxz[i,k-1]) / dz
                let dsigma_xz_dz = (self.wavefield_current.sigma_xz[[i, k]]
                    - self.wavefield_current.sigma_xz[[i, k - 1]])
                    / dz;

                // Update velocity
                // vx_new = vx_old + dt * (1/ρ) * (dsigma_xx_dx + dsigma_xz_dz)
                let vx_new = self.wavefield_current.vx[[i, k]]
                    + dt * (1.0 / rho) * (dsigma_xx_dx + dsigma_xz_dz);
                self.wavefield_current.vx[[i, k]] = vx_new;
            }
        }

        // Update vz (vertical velocity)
        // Lives at (i, k+1/2) positions
        // Equation: ∂vz/∂t = (1/ρ) × (∂σxz/∂x + ∂σzz/∂z)

        for i in 1..nx - 1 {
            for k in 1..nz - 1 {
                // Get density at vz position (pre-computed)
                let rho = self.materials.rho_vz[[i, k]];

                // Compute stress derivatives
                // σxz lives at (i+1/2, k+1/2), we need ∂σxz/∂x at (i, k+1/2)
                // ∂σxz/∂x ≈ (σxz[i,k] - σxz[i-1,k]) / dx
                let dsigma_xz_dx = (self.wavefield_current.sigma_xz[[i, k]]
                    - self.wavefield_current.sigma_xz[[i - 1, k]])
                    / dx;

                // σzz lives at (i, k), we need ∂σzz/∂z at (i, k+1/2)
                // ∂σzz/∂z ≈ (σzz[i,k+1] - σzz[i,k]) / dz
                let dsigma_zz_dz = (self.wavefield_current.sigma_zz[[i, k + 1]]
                    - self.wavefield_current.sigma_zz[[i, k]])
                    / dz;

                // Update velocity
                // vz_new = vz_old + dt * (1/ρ) * (dsigma_xz_dx + dsigma_zz_dz)
                let vz_new = self.wavefield_current.vz[[i, k]]
                    + dt * (1.0 / rho) * (dsigma_xz_dx + dsigma_zz_dz);
                self.wavefield_current.vz[[i, k]] = vz_new;
            }
        }
    }

    fn apply_sources(&mut self) {
        let t = self.current_time();

        for source in &mut self.params.sources {
            if source.should_trigger(t) {
                source.mark_triggered();
            }

            if !source.triggered {
                continue;
            }

            let tau = t - source.trigger_time.unwrap_or(0.0);

            // Only apply for ~3 periods of the dominant frequency
            let pulse_duration = 3.0 / source.frequency;
            if tau > pulse_duration {
                continue; // Stop applying after pulse is complete
            }

            let amplitude = source.ricker_wavelet(t);
            let scaled_amplitude = amplitude * source.amplitude;

            self.wavefield_current.sigma_xx[[source.x, source.z]] += scaled_amplitude;
            self.wavefield_current.sigma_zz[[source.x, source.z]] += scaled_amplitude;
        }
    }

    fn apply_boundary_conditions(&mut self) {
        let nx = self.grid.nx;
        let nz = self.grid.nz;

        // Rigid boundaries: set velocities to zero at edges
        // This is simple but causes reflections

        // Left and right boundaries
        for k in 0..nz {
            self.wavefield_current.vx[[0, k]] = 0.0;
            self.wavefield_current.vx[[nx - 1, k]] = 0.0;
            self.wavefield_current.vz[[0, k]] = 0.0;
            self.wavefield_current.vz[[nx - 1, k]] = 0.0;
        }

        // Top and bottom boundaries
        for i in 0..nx {
            self.wavefield_current.vx[[i, 0]] = 0.0;
            self.wavefield_current.vx[[i, nz - 1]] = 0.0;
            self.wavefield_current.vz[[i, 0]] = 0.0;
            self.wavefield_current.vz[[i, nz - 1]] = 0.0;
        }
    }

    pub fn run(&mut self) {
        // Main simulation loop
        println!("Starting simulation...");
        println!("Grid: {}x{}", self.grid.nx, self.grid.nz);
        println!("Time step: {:.6} s", self.params.dt);
        println!("Total time: {:.3} s", self.params.total_time());
        println!("Number of steps: {}", self.params.nt);

        while !self.is_finished() {
            self.step();

            // Print progress every 5%
            if self.current_timestep % (self.params.nt / self.params.report_period) == 0 {
                println!("Step {}/{}", self.current_timestep, self.params.nt);
            }
        }

        println!("Simulation complete!");
    }

    pub fn step(&mut self) {
        let current_time = self.current_time();

        // 1. Apply sources (inject energy) - only for sources that should trigger
        let should_apply = self.params.sources.iter_mut().any(|source| {
            if source.should_trigger(current_time) {
                source.mark_triggered();
                true
            } else {
                false
            }
        });

        if should_apply {
            self.apply_sources();
        }

        // 2. Update stresses from velocity gradients
        self.update_stresses_parallel();

        // 3. Update velocities from stress gradients
        self.update_velocities_parallel();

        // 4. Apply boundary conditions
        self.apply_boundary_conditions();

        // 5. Increment timestep
        self.current_timestep += 1;
    }

    pub fn step_serial(&mut self) {
        let current_time = self.current_time();

        // 1. Apply sources (inject energy) - only for sources that should trigger
        let should_apply = self.params.sources.iter_mut().any(|source| {
            if source.should_trigger(current_time) {
                source.mark_triggered();
                true
            } else {
                false
            }
        });

        if should_apply {
            self.apply_sources();
        }

        // 2. Update stresses from velocity gradients
        self.update_stresses();

        // 3. Update velocities from stress gradients
        self.update_velocities();

        // 4. Apply boundary conditions
        self.apply_boundary_conditions();

        // 5. Increment timestep
        self.current_timestep += 1;
    }

    /// Compute total kinetic energy in the domain
    /// KE = ∫ 0.5 * ρ * |v|² dV ≈ Σ 0.5 * ρ * (vx² + vz²) * dx * dz
    pub fn compute_kinetic_energy(&self) -> f64 {
        let dx = self.grid.dx;
        let dz = self.grid.dz;
        let cell_volume = dx * dz;
        let mut ke = 0.0;

        let (nx, nz) = self.wavefield_current.vx.dim();
        for i in 0..nx {
            for k in 0..nz {
                let vx = self.wavefield_current.vx[[i, k]];
                let vz = self.wavefield_current.vz[[i, k]];
                let v_squared = vx * vx + vz * vz;
                let rho = self.materials.rho[[i, k]];
                ke += 0.5 * rho * v_squared * cell_volume;
            }
        }
        ke
    }

    /// Compute total elastic potential energy in the domain
    /// PE = ∫ 0.5 * (σ:ε) dV where : is double dot product
    /// For 2D: PE ≈ Σ 0.5 * (σxx*exx + σzz*ezz + 2*σxz*exz) * dx * dz
    /// where exx = ∂vx/∂x, ezz = ∂vz/∂z, exz = 0.5*(∂vx/∂z + ∂vz/∂x)
    pub fn compute_potential_energy(&self) -> f64 {
        let dx = self.grid.dx;
        let dz = self.grid.dz;
        let cell_volume = dx * dz;
        let mut pe = 0.0;

        let (nx, nz) = self.wavefield_current.vx.dim();
        for i in 1..nx - 1 {
            for k in 1..nz - 1 {
                let sigma_xx = self.wavefield_current.sigma_xx[[i, k]];
                let sigma_zz = self.wavefield_current.sigma_zz[[i, k]];
                let sigma_xz = self.wavefield_current.sigma_xz[[i, k]];

                // Compute strain components using finite differences
                let exx = (self.wavefield_current.vx[[i, k]] - self.wavefield_current.vx[[i - 1, k]]) / dx;
                let ezz = (self.wavefield_current.vz[[i, k]] - self.wavefield_current.vz[[i, k - 1]]) / dz;
                let exz = 0.5 * (
                    (self.wavefield_current.vx[[i, k + 1]] - self.wavefield_current.vx[[i, k - 1]]) / (2.0 * dz)
                    + (self.wavefield_current.vz[[i + 1, k]] - self.wavefield_current.vz[[i - 1, k]]) / (2.0 * dx)
                );

                let strain_energy = sigma_xx * exx + sigma_zz * ezz + 2.0 * sigma_xz * exz;
                pe += 0.5 * strain_energy * cell_volume;
            }
        }
        pe
    }

    /// Compute and report total energy (kinetic + potential)
    pub fn report_energy(&self) {
        let ke = self.compute_kinetic_energy();
        let pe = self.compute_potential_energy();
        let total = ke + pe;
        println!(
            "t={:.4}s | KE={:.6e} | PE={:.6e} | Total={:.6e}",
            self.current_time(),
            ke,
            pe,
            total
        );
    }

    pub fn run_with_visualisation(
        &mut self,
        visualise_interval: usize,
        field: &str, // "vx", "vz", "vmag", "divergence", "curl"
        iteration: i64,
    ) {
        println!("Starting simulation with visualisation...");
        println!("Grid: {}x{}", self.grid.nx, self.grid.nz);
        println!("Time step: {:.6} s", self.params.dt);
        println!("Total time: {:.3} s", self.params.total_time());
        println!("Visualising '{}' every {} steps", field, visualise_interval);
        let output_directory = vec!["output".to_string(), iteration.to_string()].join("/");
        // Create visualiser
        let visualiser =
            WavefieldVisualiser::new(&output_directory, 1200, 1000, self.grid.dx, self.grid.dz);

        // Save initial state
        self.visualise_field(&visualiser, field);
        self.report_energy();

        while !self.is_finished() {
            self.step();

            // Visualise
            if self.current_timestep % visualise_interval == 0 {
                self.visualise_field(&visualiser, field);
            }

            // Print progress and energy every 100 steps
            if self.current_timestep % 100 == 0 {
                println!(
                    "Step {}/{} (t={:.4}s)",
                    self.current_timestep,
                    self.params.nt,
                    self.current_time()
                );
                self.report_energy();
            }
        }

        println!("Simulation complete!");
        println!("Frames saved to output/ directory");
    }

    fn visualise_field(&self, visualiser: &WavefieldVisualiser, field: &str) {
        let data = match field {
            "vx" => self.wavefield_current.vx.clone(),
            "vz" => self.wavefield_current.vz.clone(),
            "vmag" => self.wavefield_current.compute_velocity_magnitude(),
            "divergence" => self
                .wavefield_current
                .compute_divergence(self.grid.dx, self.grid.dz),
            "curl" => self
                .wavefield_current
                .compute_curl(self.grid.dx, self.grid.dz),
            "sigma_xx" => self.wavefield_current.sigma_xx.clone(),
            "sigma_zz" => self.wavefield_current.sigma_zz.clone(),
            "sigma_xz" => self.wavefield_current.sigma_xz.clone(),
            _ => panic!("Unknown field: {}", field),
        };

        if let Err(e) =
            visualiser.plot_field(&data, self.current_timestep, field, self.current_time())
        {
            eprintln!("Warning: Failed to visualise: {}", e);
        }
    }
    pub fn run_with_p_s_visualisation(&mut self, visualise_interval: usize, iteration: i64) {
        println!("Starting simulation with P-S wave visualisation...");
        println!("Grid: {}x{}", self.grid.nx, self.grid.nz);
        println!("Time step: {:.6} s", self.params.dt);
        println!("Total time: {:.3} s", self.params.total_time());
        println!(
            "Visualising P and S waves every {} steps",
            visualise_interval
        );

        let output_directory = format!("output/{}", iteration);

        // Create visualiser
        let visualiser =
            WavefieldVisualiser::new(&output_directory, 1200, 1000, self.grid.dx, self.grid.dz);

        // Save initial state
        let divergence = self
            .wavefield_current
            .compute_divergence(self.grid.dx, self.grid.dz);
        let curl = self
            .wavefield_current
            .compute_curl(self.grid.dx, self.grid.dz);
        visualiser
            .plot_p_and_s_overlay(&divergence, &curl, 0, 0.0)
            .unwrap();

        while !self.is_finished() {
            self.step();

            // Visualise
            if self.current_timestep % visualise_interval == 0 {
                let divergence = self
                    .wavefield_current
                    .compute_divergence(self.grid.dx, self.grid.dz);
                let curl = self
                    .wavefield_current
                    .compute_curl(self.grid.dx, self.grid.dz);

                if let Err(e) = visualiser.plot_p_and_s_overlay(
                    &divergence,
                    &curl,
                    self.current_timestep / visualise_interval,
                    self.current_time(),
                ) {
                    eprintln!("Warning: Failed to visualise: {}", e);
                }
            }

            // Print progress
            if self.current_timestep % 100 == 0 {
                println!(
                    "Step {}/{} (t={:.4}s)",
                    self.current_timestep,
                    self.params.nt,
                    self.current_time()
                );
            }
        }

        println!("Simulation complete!");
        println!("P-S overlay frames saved to {}/", output_directory);
    }

    fn update_stresses_parallel(&mut self) {
        let dx = self.grid.dx;
        let dz = self.grid.dz;
        let dt = self.params.dt;
        let nx = self.grid.nx;
        let nz = self.grid.nz;

        // For each row, process in parallel
        // This gives us good parallelization while maintaining safe access

        // Process normal stresses
        let indices: Vec<(usize, usize)> = (1..nx - 1)
            .flat_map(|i| (1..nz - 1).map(move |k| (i, k)))
            .collect();

        let updates: Vec<(usize, usize, f64, f64)> = indices
            .par_iter()
            .map(|&(i, k)| {
                let lambda = self.materials.lambda[[i, k]];
                let mu = self.materials.mu[[i, k]];
                let lambda_plus_2mu = lambda + 2.0 * mu;

                // Compute velocity derivatives
                let dvx_dx = (self.wavefield_current.vx[[i, k]]
                    - self.wavefield_current.vx[[i - 1, k]])
                    / dx;
                let dvz_dz = (self.wavefield_current.vz[[i, k]]
                    - self.wavefield_current.vz[[i, k - 1]])
                    / dz;

                let dsigma_xx = dt * (lambda_plus_2mu * dvx_dx + lambda * dvz_dz);
                let dsigma_zz = dt * (lambda * dvx_dx + lambda_plus_2mu * dvz_dz);

                (i, k, dsigma_xx, dsigma_zz)
            })
            .collect();

        // Apply updates
        for (i, k, dsigma_xx, dsigma_zz) in updates {
            self.wavefield_current.sigma_xx[[i, k]] += dsigma_xx;
            self.wavefield_current.sigma_zz[[i, k]] += dsigma_zz;
        }

        // Update shear stress (sigma_xz)
        let indices_xz: Vec<(usize, usize)> = (0..nx - 1)
            .flat_map(|i| (0..nz - 1).map(move |k| (i, k)))
            .collect();

        let updates_xz: Vec<(usize, usize, f64)> = indices_xz
            .par_iter()
            .map(|&(i, k)| {
                let mu = self.materials.mu_xz[[i, k]];

                let dvx_dz = (self.wavefield_current.vx[[i, k + 1]]
                    - self.wavefield_current.vx[[i, k]])
                    / dz;
                let dvz_dx = (self.wavefield_current.vz[[i + 1, k]]
                    - self.wavefield_current.vz[[i, k]])
                    / dx;

                let dsigma_xz = dt * mu * (dvx_dz + dvz_dx);

                (i, k, dsigma_xz)
            })
            .collect();

        for (i, k, dsigma_xz) in updates_xz {
            self.wavefield_current.sigma_xz[[i, k]] += dsigma_xz;
        }
    }

    fn update_velocities_parallel(&mut self) {
        let dx = self.grid.dx;
        let dz = self.grid.dz;
        let dt = self.params.dt;
        let nx = self.grid.nx;
        let nz = self.grid.nz;

        // Update vx
        let indices: Vec<(usize, usize)> = (1..nx - 1)
            .flat_map(|i| (1..nz - 1).map(move |k| (i, k)))
            .collect();

        let updates_vx: Vec<(usize, usize, f64)> = indices
            .par_iter()
            .map(|&(i, k)| {
                let rho = self.materials.rho_vx[[i, k]];

                let dsigma_xx_dx = (self.wavefield_current.sigma_xx[[i + 1, k]]
                    - self.wavefield_current.sigma_xx[[i, k]])
                    / dx;
                let dsigma_xz_dz = (self.wavefield_current.sigma_xz[[i, k]]
                    - self.wavefield_current.sigma_xz[[i, k - 1]])
                    / dz;

                let dvx = dt * (1.0 / rho) * (dsigma_xx_dx + dsigma_xz_dz);

                (i, k, dvx)
            })
            .collect();

        for (i, k, dvx) in updates_vx {
            self.wavefield_current.vx[[i, k]] += dvx;
        }

        // Update vz
        let updates_vz: Vec<(usize, usize, f64)> = indices
            .par_iter()
            .map(|&(i, k)| {
                let rho = self.materials.rho_vz[[i, k]];

                let dsigma_xz_dx = (self.wavefield_current.sigma_xz[[i, k]]
                    - self.wavefield_current.sigma_xz[[i - 1, k]])
                    / dx;
                let dsigma_zz_dz = (self.wavefield_current.sigma_zz[[i, k + 1]]
                    - self.wavefield_current.sigma_zz[[i, k]])
                    / dz;

                let dvz = dt * (1.0 / rho) * (dsigma_xz_dx + dsigma_zz_dz);

                (i, k, dvz)
            })
            .collect();

        for (i, k, dvz) in updates_vz {
            self.wavefield_current.vz[[i, k]] += dvz;
        }
    }
}

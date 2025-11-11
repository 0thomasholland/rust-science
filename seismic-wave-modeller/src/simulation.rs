use crate::grid::Grid;
use crate::materials::MaterialProperties;
use crate::visualisation::WavefieldVisualiser;
use crate::wavefield::Wavefield;
use ndarray::Zip;
use rayon::prelude::*;
use std::f64::consts::PI;

#[derive(Clone, Debug)]
pub struct Source {
    pub i: usize, // Grid position x
    pub k: usize, // Grid position z
    pub f0: f64,  // Peak frequency (Hz)
    pub t0: f64,  // Time delay (s)
}

impl Source {
    pub fn new(i: usize, k: usize, f0: f64) -> Self {
        // Auto-calculate t0 from f0
        let t0 = 1.2 / f0;
        Self { i, k, f0, t0 }
    }

    pub fn with_delay(i: usize, k: usize, f0: f64, t0: f64) -> Self {
        // If you want to manually specify t0
        Self { i, k, f0, t0 }
    }

    pub fn ricker_wavelet(&self, t: f64) -> f64 {
        // Ricker wavelet for this source
        let tau = t - self.t0;
        let arg = (PI * self.f0 * tau).powi(2);
        (1.0 - 2.0 * arg) * (-arg).exp()
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
        cfl_safety: f64,
    ) -> Self {
        // Calculate appropriate t0 from f0
        // Typically t0 = 1.2 / f0 so the wavelet starts near zero
        let source_t0 = 1.2 / source_f0;

        Self {
            dt,
            nt,
            report_period,
            sources: vec![Source::with_delay(source_i, source_k, source_f0, source_t0)],
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
            if source.i >= grid.nx || source.k >= grid.nz {
                panic!(
                    "Source {} at ({}, {}) is outside grid bounds ({}, {})",
                    idx, source.i, source.k, grid.nx, grid.nz
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

        for source in &self.params.sources {
            // Get wavelet amplitude at current time
            let amplitude = source.ricker_wavelet(t);

            // Apply as explosion source (isotropic expansion)
            // Add to normal stresses at source location
            // self.wavefield_current.sigma_xx[[source.i, source.k]] += amplitude;
            // self.wavefield_current.sigma_zz[[source.i, source.k]] += amplitude;
            // Note: We're adding directly to stress, not using +=
            // The amplitude needs to be scaled properly
            // A common scaling is: amplitude * (lambda + 2*mu)
            // This gives a physically meaningful force

            let lambda = self.materials.lambda[[source.i, source.k]];
            let mu = self.materials.mu[[source.i, source.k]];
            let scale = lambda + 2.0 * mu;
            let scaled_amplitude = amplitude * scale;
            self.wavefield_current.sigma_xx[[source.i, source.k]] += scaled_amplitude;
            self.wavefield_current.sigma_zz[[source.i, source.k]] += scaled_amplitude;
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
        // 1. Apply sources (inject energy)
        self.apply_sources();

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
        // 1. Apply sources (inject energy)
        self.apply_sources();

        // 2. Update stresses from velocity gradients
        self.update_stresses();

        // 3. Update velocities from stress gradients
        self.update_velocities();

        // 4. Apply boundary conditions
        self.apply_boundary_conditions();

        // 5. Increment timestep
        self.current_timestep += 1;
    }

    pub fn run_with_visualisation(
        &mut self,
        visualise_interval: usize,
        field: &str, // "vx", "vz", "vmag", "divergence", "curl"
    ) {
        println!("Starting simulation with visualisation...");
        println!("Grid: {}x{}", self.grid.nx, self.grid.nz);
        println!("Time step: {:.6} s", self.params.dt);
        println!("Total time: {:.3} s", self.params.total_time());
        println!("Visualising '{}' every {} steps", field, visualise_interval);

        // Create visualiser
        let visualiser = WavefieldVisualiser::new("output", 1200, 1000);

        // Save initial state
        self.visualise_field(&visualiser, field);

        while !self.is_finished() {
            self.step();

            // Visualise
            if self.current_timestep % visualise_interval == 0 {
                self.visualise_field(&visualiser, field);
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

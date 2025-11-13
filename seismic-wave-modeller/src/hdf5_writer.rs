use crate::simulation::Simulation;
use hdf5::{File, Result as H5Result};
use ndarray::Array2;

pub struct HDF5Writer {
    file: File,
    save_interval: usize,
    fields: Vec<String>,
}

impl HDF5Writer {
    pub fn new(filename: &str, save_interval: usize, fields: Vec<&str>) -> H5Result<Self> {
        // prepend output/ to filename
        let filename = format!("output/{}", filename);
        let file = File::create(filename)?;

        Ok(Self {
            file,
            save_interval,
            fields: fields.iter().map(|s| s.to_string()).collect(),
        })
    }

    pub fn write_metadata(&self, sim: &Simulation) -> H5Result<()> {
        // Create metadata group
        let meta = self.file.create_group("metadata")?;

        // Grid parameters
        meta.new_dataset::<usize>()
            .create("nx")?
            .write_scalar(&sim.grid.nx)?;
        meta.new_dataset::<usize>()
            .create("nz")?
            .write_scalar(&sim.grid.nz)?;
        meta.new_dataset::<f64>()
            .create("dx")?
            .write_scalar(&sim.grid.dx)?;
        meta.new_dataset::<f64>()
            .create("dz")?
            .write_scalar(&sim.grid.dz)?;

        // Simulation parameters
        meta.new_dataset::<f64>()
            .create("dt")?
            .write_scalar(&sim.params.dt)?;
        meta.new_dataset::<usize>()
            .create("nt")?
            .write_scalar(&sim.params.nt)?;
        meta.new_dataset::<usize>()
            .create("save_interval")?
            .write_scalar(&self.save_interval)?;
        meta.new_dataset::<f64>()
            .create("cfl_safety")?
            .write_scalar(&sim.params.cfl_safety)?;

        // Source information
        let n_sources = sim.params.sources.len();
        let mut source_x = vec![0usize; n_sources];
        let mut source_z = vec![0usize; n_sources];
        let mut source_freq = vec![0.0f64; n_sources];
        let mut source_amp = vec![0.0f64; n_sources];

        for (idx, source) in sim.params.sources.iter().enumerate() {
            source_x[idx] = source.x;
            source_z[idx] = source.z;
            source_freq[idx] = source.frequency;
            source_amp[idx] = source.amplitude;
        }

        meta.new_dataset::<usize>()
            .shape([n_sources])
            .create("source_x")?
            .write(&source_x)?;
        meta.new_dataset::<usize>()
            .shape([n_sources])
            .create("source_z")?
            .write(&source_z)?;
        meta.new_dataset::<f64>()
            .shape([n_sources])
            .create("source_frequency")?
            .write(&source_freq)?;
        meta.new_dataset::<f64>()
            .shape([n_sources])
            .create("source_amplitude")?
            .write(&source_amp)?;

        Ok(())
    }
    pub fn write_materials(&self, sim: &Simulation) -> H5Result<()> {
        let materials = self.file.create_group("materials")?;

        let (nx, nz) = sim.materials.lambda.dim();

        // Lambda - write as 1D
        let lambda_flat: Vec<f64> = sim.materials.lambda.iter().copied().collect();
        materials
            .new_dataset::<f64>()
            .shape([nx * nz]) // 1D shape
            .create("lambda")?
            .write(&lambda_flat)?;

        // Mu
        let mu_flat: Vec<f64> = sim.materials.mu.iter().copied().collect();
        materials
            .new_dataset::<f64>()
            .shape([nx * nz])
            .create("mu")?
            .write(&mu_flat)?;

        // Rho
        let rho_flat: Vec<f64> = sim.materials.rho.iter().copied().collect();
        materials
            .new_dataset::<f64>()
            .shape([nx * nz])
            .create("rho")?
            .write(&rho_flat)?;

        // Vp and Vs
        let mut vp = Array2::<f64>::zeros((nx, nz));
        let mut vs = Array2::<f64>::zeros((nx, nz));

        for i in 0..nx {
            for k in 0..nz {
                let lambda = sim.materials.lambda[[i, k]];
                let mu = sim.materials.mu[[i, k]];
                let rho = sim.materials.rho[[i, k]];
                vp[[i, k]] = ((lambda + 2.0 * mu) / rho).sqrt();
                vs[[i, k]] = (mu / rho).sqrt();
            }
        }

        let vp_flat: Vec<f64> = vp.iter().copied().collect();
        materials
            .new_dataset::<f64>()
            .shape([nx * nz])
            .create("vp")?
            .write(&vp_flat)?;

        let vs_flat: Vec<f64> = vs.iter().copied().collect();
        materials
            .new_dataset::<f64>()
            .shape([nx * nz])
            .create("vs")?
            .write(&vs_flat)?;

        Ok(())
    }

    pub fn write_snapshot(&self, sim: &Simulation, frame_number: usize) -> H5Result<()> {
        let (nx, nz) = (sim.grid.nx, sim.grid.nz);
        let timestep = sim.current_timestep;

        for field_name in &self.fields {
            let group_path = format!("fields/{}", field_name);
            let field_group = if self.file.group(&group_path).is_ok() {
                self.file.group(&group_path)?
            } else {
                self.file.create_group(&group_path)?
            };

            // Get the data
            let data = match field_name.as_str() {
                "vx" => sim.wavefield_current.vx.clone(),
                "vz" => sim.wavefield_current.vz.clone(),
                "sigma_xx" => sim.wavefield_current.sigma_xx.clone(),
                "sigma_zz" => sim.wavefield_current.sigma_zz.clone(),
                "sigma_xz" => sim.wavefield_current.sigma_xz.clone(),
                "vmag" => sim.wavefield_current.compute_velocity_magnitude(),
                "divergence" => sim
                    .wavefield_current
                    .compute_divergence(sim.grid.dx, sim.grid.dz),
                "curl" => sim.wavefield_current.compute_curl(sim.grid.dx, sim.grid.dz),
                _ => continue,
            };

            // Convert to flat vector (1D)
            let data_flat: Vec<f64> = data.iter().copied().collect();

            // Create dataset for this frame as 1D
            let dataset_name = format!("frame_{:06}", frame_number);
            field_group
                .new_dataset::<f64>()
                .shape([nx * nz]) // 1D shape
                .create(dataset_name.as_str())?
                .write(&data_flat)?;

            // Store metadata as attributes
            let dataset = field_group.dataset(dataset_name.as_str())?;
            dataset
                .new_attr::<usize>()
                .create("timestep")?
                .write_scalar(&timestep)?;
            dataset
                .new_attr::<f64>()
                .create("time")?
                .write_scalar(&sim.current_time())?;
        }

        Ok(())
    }
}

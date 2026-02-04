// TESTS //
mod config;
mod grid;
mod materials;
mod simulation;
mod visualisation;
mod wavefield;

use anyhow::Result;
use clap::Parser;
use config::Config;
use grid::Grid;
use materials::MaterialProperties;
use ndarray::Array2;
use simulation::{Simulation, SimulationParams, Source};

#[derive(Parser)]
#[command(name = "seismic-wave-modeller")]
#[command(about = "2D elastic wave propagation simulator", long_about = None)]
struct Cli {
    /// Path to TOML configuration file
    #[arg(short, long)]
    config: String,

    /// Override output directory iteration number
    #[arg(short, long)]
    iteration: Option<i64>,

    /// Verbose output
    #[arg(short, long)]
    verbose: bool,
}

fn run_simulation(config: Config, verbose: bool) -> Result<()> {
    // Print configuration summary if verbose
    if verbose {
        config.print_summary();
    }

    // Get simulation parameters
    let nx = config.grid.nx;
    let nz = config.grid.nz;
    let dx = config.grid.dx;
    let dz = config.grid.dz;
    let dt = config.simulation.dt.unwrap();
    let total_time = config.simulation.total_time;
    let nt = config.simulation.compute_nt(dt);

    // Build grid
    let grid = Grid::new(nx, nz, dx, dz);

    // Build homogeneous material properties
    let vp = Array2::from_elem((nx, nz), config.materials.vp);
    let vs = Array2::from_elem((nx, nz), config.materials.vs);
    let rho = Array2::from_elem((nx, nz), config.materials.rho);
    let materials = MaterialProperties::new(vp, vs, rho);

    // Build sources from configuration
    let sources: Vec<Source> = config
        .sources
        .into_iter()
        .map(|src| {
            if src.trigger_time > 0.0 {
                Source::with_delay(src.x, src.z, src.amplitude, src.frequency, src.trigger_time)
            } else {
                Source::new(src.x, src.z, src.amplitude, src.frequency)
            }
        })
        .collect();

    // Build simulation parameters
    let sim_params = SimulationParams {
        dt,
        nt,
        report_period: config.simulation.report_period,
        sources,
        cfl_safety: config.simulation.cfl_safety,
    };

    // Verify CFL condition
    sim_params.check_cfl(dx, dz, config.materials.vp);

    // Print summary
    println!("Running simulation with dt = {} seconds", dt);
    println!("Simulation length: {} seconds", total_time);
    println!("Total steps: {}", nt);
    println!("Visualization time: {} seconds", config.visualization.video_length);

    // Initialize simulation
    let mut sim = Simulation::new(grid, materials, sim_params);

    // Calculate visualization step rate
    let step_rate = (((nt as f64) / (config.visualization.video_length * config.visualization.fps)) as f64)
        .round() as usize;
    let step_rate = step_rate.max(1); // Ensure at least 1 step per frame

    if verbose {
        println!("Visualization step rate: {} steps per frame", step_rate);
    }

    // Run simulation with visualization
    sim.run_with_visualisation(step_rate, &config.visualization.field, config.visualization.iteration);

    Ok(())
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    // Load configuration from file
    let config = Config::from_file(&cli.config)?;

    // Override iteration if provided via CLI
    let mut config = config;
    if let Some(iteration) = cli.iteration {
        config.visualization.iteration = iteration;
    }

    // Run simulation
    run_simulation(config, cli.verbose)?;

    Ok(())
}


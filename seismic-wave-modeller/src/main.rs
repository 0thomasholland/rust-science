mod grid;
mod hdf5_writer;
mod materials;
mod simulation;
mod wavefield;

use grid::Grid;
use materials::MaterialProperties;
use ndarray::Array2;
use simulation::{Simulation, SimulationParams, Source};
use std::time::Instant;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let t0 = Instant::now();

    rayon::ThreadPoolBuilder::new()
        .num_threads(32) // All threads on 9950X
        .build_global()
        .unwrap();

    let nx = 800;
    let nz = 800;
    let dx = 5.0;
    let dz = 5.0;

    let grid = Grid::new(nx, nz, dx, dz);

    let vp = Array2::from_elem((nx, nz), 3000.0);
    let vs = Array2::from_elem((nx, nz), 1732.0);
    let rho = Array2::from_elem((nx, nz), 2500.0);
    let materials = MaterialProperties::new(vp, vs, rho);

    let vp_max = 3000.0;
    let cfl_safety = 0.4;

    // Calculate stable dt
    let dt_stable = SimulationParams::compute_stable_dt_static(dx, dz, vp_max, cfl_safety);

    println!("Stable dt: {:.6} s", dt_stable);

    // Create sources
    let sources = vec![
        Source::new(nx / 2, nz / 2, 25.0, 1.0), // Adjust based on your Source::new signature
    ];

    let simulated_time = 0.1; // seconds
    let number_steps = (simulated_time / dt_stable) as usize;
    let params = SimulationParams::new_multi_source(
        dt_stable,
        number_steps,
        100, // report_period
        sources,
        cfl_safety,
    );

    let mut sim = Simulation::new(grid, materials, params);

    let video_length = 5.0; // seconds
    let frame_rate = 30.0; // frames per second
    let total_frames = (video_length * frame_rate) as usize;
    let steps_per_frame = sim.params.nt / total_frames;
    // Run with HDF5 output
    println!("Total frames to save: {}", total_frames);
    println!(" (every {} steps)", steps_per_frame);
    println!();
    println!("Total simulation steps: {}", sim.params.nt);
    println!(
        "Expected simulated length: {:.2} seconds",
        sim.params.dt * (sim.params.nt as f64)
    );
    println!(
        "Expected video length: {:.2} seconds",
        (total_frames as f64) / frame_rate
    );
    let t1 = Instant::now();
    sim.run_with_hdf5(
        "11.h5",
        steps_per_frame, // Save every steps_per_frame steps
        vec!["vmag", "divergence", "curl"],
        10,
    )?;

    let t2 = Instant::now();

    println!("Initialization time: {:.2?}", t1.duration_since(t0));
    println!("Simulation time: {:.2?}", t2.duration_since(t1));
    println!("\nVisualize with: python3 visualize_hdf5.py");

    Ok(())
}

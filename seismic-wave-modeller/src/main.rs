// TESTS //
mod grid;
mod materials;
mod simulation;
mod visualisation;
mod wavefield;

use grid::Grid;
use materials::MaterialProperties;
use ndarray::Array2;
use simulation::{Simulation, SimulationParams, Source};
use wavefield::Wavefield;

fn run() {
    // Define grid
    let nx = 800;
    let nz = 800;
    let dx = 0.125; // meters
    let dz = 0.125;

    let grid = Grid::new(nx, nz, dx, dz);

    // Homogeneous material
    let vp = Array2::from_elem((nx, nz), 6000.0); // 6 km/s
    let vs = Array2::from_elem((nx, nz), 4000.0); // 4 km/s
    let rho = Array2::from_elem((nx, nz), 3000.0); // 3 g/cm³

    let materials = MaterialProperties::new(vp, vs, rho);
    let source0 = Source::new(nx / 2, nz / 2, 25.0, 10.0);
    let source1 = Source::new(nx / 4, nz / 2, 25.0, 10.0);
    let source2 = Source::new(3 * nx / 4, nz / 2, 25.0, 10.0);

    let dt = 0.00001; // 40 microseconds
    let length = 0.012; // simulated seconds
    let video_length = 10.0; // seconds
    let fps = 30.0; // frames per second
    let params = SimulationParams {
        dt: dt,
        nt: (length / dt).round() as usize,
        report_period: 100,
        sources: vec![source1, source2],
        cfl_safety: 0.5,
    };

    println!("Running simulation with dt = {} seconds", dt);
    println!("Simulation length: {} simulated seconds", length);
    println!("Total steps: {}", params.nt);
    println!("Visualisation time: {} seconds", video_length);
    // println!("Sampling rate: {} simulations steps per frame", (params.nt as f64 / video_length * fps) as f64);

    // Initialize simulation
    let mut sim = Simulation::new(grid, materials, params);
    let step_rate = (((length / dt) as f64 / (video_length * fps)) as f64).round() as usize; // 30 FPS type rounded then as usize
    sim.run_with_visualisation(step_rate, "vmag", 9);
}

fn safe_dt() {
    let dt = 0.00002;
    let mut nt = 2;
    let dx = 0.125;
    let dz = 0.125;
    let vp_max = 6000.0;
    let cfl_safety = 0.5;

    let params = SimulationParams {
        dt,
        nt,
        report_period: 1,
        sources: vec![],
        cfl_safety,
    };
    let mut safe_dt = params.compute_stable_dt(dx, dz, vp_max);
    println!("Computed stable dt: {}", safe_dt);

    // calculate sampling rate

    let video_length = 7.0;
    let fps = 30.0;
    let simulation_length = 0.012;
    safe_dt *= 1.0; // higher step resolution
    nt = (simulation_length / safe_dt) as usize;
    let step_rate = (((nt as f64) / (video_length * fps)) as f64).round() as usize;
    println!("Sampling rate (steps per frame): {}", step_rate);
    println!("Total steps: {}", nt);
    let actual_video_length = (nt as f64) / (step_rate as f64 * fps);
    println!(
        "Actual video length at {} FPS: {} seconds",
        fps, actual_video_length
    );
}

fn main() {
    run();
    // safe_dt();
}

fn safe_params() {
    //calcualtes safe simulation parameters
    let nx = 100;
    let nz = 100;
    let dx = 0.1; // meters
    let dz = 0.1;

    let grid = Grid::new(nx, nz, dx, dz);

    // Homogeneous material
    let vp = Array2::from_elem((nx, nz), 6000.0); // 6 km/s
    let vs = Array2::from_elem((nx, nz), 4000.0); // 4 km/s
    let rho = Array2::from_elem((nx, nz), 3000.0); // 3 g/cm³

    let materials = MaterialProperties::new(vp, vs, rho);

    let cfl_safety = 0.5;
    let parameters = SimulationParams {
        dt: 0.0, // Placeholder
        nt: 1000,
        report_period: 100,
        sources: vec![],
        cfl_safety,
    };
    let dt = parameters.compute_stable_dt(dx, dz, 6000.0); // Using max velocity
    println!("Calculated stable timestep: {} seconds", dt);
}

fn main() {
    run();
}

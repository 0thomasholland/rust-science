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
    let nx = 400;
    let nz = 400;
    let dx = 0.1; // meters
    let dz = 0.1;

    let grid = Grid::new(nx, nz, dx, dz);

    // Homogeneous material
    let vp = Array2::from_elem((nx, nz), 6000.0); // 6 km/s
    let vs = Array2::from_elem((nx, nz), 4000.0); // 4 km/s
    let rho = Array2::from_elem((nx, nz), 3000.0); // 3 g/cm³

    let materials = MaterialProperties::new(vp, vs, rho);
    let source = Source::new(nx / 2, nz / 2, 25.0);

    let dt = 0.000008;
    let length = 0.004;
    let params = SimulationParams {
        dt: dt, // We'll calculate proper value later
        nt: (length / dt) as usize,
        report_period: 100,
        sources: vec![source],
        cfl_safety: 0.5,
    };

    // Initialize simulation
    let mut sim = Simulation::new(grid, materials, params);

    // Run just a few steps
    sim.run_with_visualisation(10, "vmag");

    // "ffmpeg -framerate 30 -pattern_type glob -i 'output/2/vmag_*.png' -c:v libx264 -pix_fmt yuv420p output/2.mp4");
}

fn safe_params() { //calcualtes safe simulation parameters
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

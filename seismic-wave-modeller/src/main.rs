// TESTS //
mod grid;
mod materials;
mod simulation;
mod wavefield;

use grid::Grid;
use materials::MaterialProperties;
use ndarray::Array2;
use simulation::{Simulation, SimulationParams, Source};
use wavefield::Wavefield;

fn main() {
    // Define grid
    let nx = 200;
    let nz = 200;
    let dx = 1.0; // meters
    let dz = 1.0;

    let grid = Grid::new(nx, nz, dx, dz);

    // Homogeneous material
    let vp = Array2::from_elem((nx, nz), 6000.0); // 6 km/s
    let vs = Array2::from_elem((nx, nz), 4000.0); // 4 km/s
    let rho = Array2::from_elem((nx, nz), 3000.0); // 3 g/cm³

    let materials = MaterialProperties::new(vp, vs, rho);
    let source = Source::new(nx / 2, nz / 2, 25.0);

    let params = SimulationParams {
        dt: 0.00008, // We'll calculate proper value later
        nt: 1000,
        report_period: 100,
        sources: vec![source],
        cfl_safety: 0.5,
    };

    // Initialize simulation
    let mut sim = Simulation::new(grid, materials, params);

    // Run just a few steps
    for step in 0..100 {
        sim.step();
        let max_vx = sim
            .wavefield_current
            .vx
            .iter()
            .cloned()
            .fold(0.0_f64, f64::max);
        println!("Step {}: max vx = {:.6e}", step, max_vx);
    }
}

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

fn main() {
    // Define grid
    let nx = 400;
    let nz = 400;
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
        dt: 0.000008, // We'll calculate proper value later
        nt: 2000,
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

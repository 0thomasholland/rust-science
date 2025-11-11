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
    let dx = 1.0; // meters
    let dz = 1.0;

    let grid = Grid::new(nx, nz, dx, dz);

    // Homogeneous material
    let vp = Array2::from_elem((nx, nz), 6000.0); // 6 km/s
    let vs = Array2::from_elem((nx, nz), 4000.0); // 4 km/s
    let rho = Array2::from_elem((nx, nz), 3000.0); // 3 g/cm³

    let materials = MaterialProperties::new(vp, vs, rho);
    let source0 = Source::with_delay(nx / 2, nz / 2, 25.0,0.0);
    let source1 = Source::with_delay(nx / 4, nz / 2, 25.0, 0.002);
    let source2 = Source::with_delay(3* nx / 4, nz / 2, 15.0, 0.002);

    let dt = 0.00004; // 40 microseconds
    let length = 0.005; // simulated seconds
    let video_length = 5.0; // seconds
    let fps = 30.0; // frames per second
    let params = SimulationParams {
        dt: 0.000008, // We'll calculate proper value later
        nt: (length / dt) as usize,
        report_period: 100,
        sources: vec![source0, source1, source2],
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
    sim.run_with_visualisation(step_rate, "vmag");

    // "ffmpeg -framerate 30 -pattern_type glob -i 'output/2/vmag_*.png' -c:v libx264 -pix_fmt yuv420p output/2.mp4");

}


fn main() {
    run();
}

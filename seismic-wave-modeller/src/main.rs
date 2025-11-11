mod grid;
use grid::Grid;

// struct Simulation {
//     grid: Grid,
//     materials: MaterialProperties,
//     wavefield_current: Wavefield,
//     wavefield_old: Wavefield, // For time-stepping
//     params: SimulationParams,
//     current_time_step: usize,
// }

// impl Simulation {
//     // Constructor
//     // new(grid, materials, params) -> Self
//     // Initialize both wavefields to zero

//     // Method: Get current simulation time
//     // current_time(&self) -> f64

//     // Method: Swap wavefields (current becomes old, prepare for new current)
//     // swap_wavefields(&mut self)

//     // Placeholder for now - we'll implement these next:
//     // step(&mut self)  // Advance one time step
//     // run(&mut self)   // Main loop
//     // apply_source(&mut self)  // Add source term
// }

fn main() {
    let grid = Grid::new(100, 100, 10.0, 10.0);

    println!("Grid dimensions: {} x {}", grid.nx, grid.nz);
    println!("Grid spacing: {} m x {} m", grid.dx, grid.dz);
    println!("Domain size: {} m x {} m", grid.width(), grid.height());

    // Test coordinates
    println!(
        "Point (0, 0) is at ({}, {})",
        grid.x_coord(0),
        grid.z_coord(0)
    );
    println!(
        "Point (10, 10) is at ({}, {})",
        grid.x_coord(10),
        grid.z_coord(10)
    );

    // Test bounds checking
    println!("(50, 50) in bounds? {}", grid.in_bounds(50, 50));
    println!("(200, 50) in bounds? {}", grid.in_bounds(200, 50));
}

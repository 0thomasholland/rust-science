pub struct Grid {
    pub nx: usize, // Number of points in x direction
    pub nz: usize, // Number of points in z direction
    pub dx: f64,   // Grid spacing in x (meters)
    pub dz: f64,   // Grid spacing in z (meters)
}

impl Grid {
    pub fn new(nx: usize, nz: usize, dx: f64, dz: f64) -> Self {
        // Constructor
        Grid { nx, nz, dx, dz }
    }
    pub fn x_coord(&self, i: usize) -> f64 {
        // Convert grid index i to physical x coordinate
        self.dx * (i as f64)
    }
    pub fn z_coord(&self, k: usize) -> f64 {
        // Convert grid index k to physical z coordinate
        self.dz * (k as f64)
    }
    pub fn in_bounds(&self, i: usize, k: usize) -> bool {
        // Check if (i, k) is within grid bounds
        i < self.nx && k < self.nz
    }
    pub fn width(&self) -> f64 {
        // Total width of domain in x direction
        (self.nx - 1) as f64 * self.dx
    }

    pub fn height(&self) -> f64 {
        // Total height of domain in z direction
        (self.nz - 1) as f64 * self.dz
    }
}

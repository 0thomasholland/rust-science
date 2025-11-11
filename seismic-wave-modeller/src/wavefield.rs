struct Wavefield {
    // Five components
    vx: Array2<f64>,
    vz: Array2<f64>,
    sigma_xx: Array2<f64>,
    sigma_zz: Array2<f64>,
    sigma_xz: Array2<f64>,
}

impl Wavefield {
    // Constructor: Initialize all to zeros
    // new(nx: usize, nz: usize) -> Self
    
    // Method: Set all fields to zero (for resetting)
    // zero(&mut self)
    
    // Optional helper methods:
    // compute_divergence() -> Array2<f64>  // For visualization
    // compute_curl() -> Array2<f64>  // For visualization
}
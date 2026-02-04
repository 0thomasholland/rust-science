use ndarray::Array2;

pub struct Wavefield {
    // Five components
    pub vx: Array2<f64>,
    pub vz: Array2<f64>,
    pub sigma_xx: Array2<f64>,
    pub sigma_zz: Array2<f64>,
    pub sigma_xz: Array2<f64>,
}

impl Wavefield {
    pub fn new(nx: usize, nz: usize) -> Self {
        // Create all 5 arrays as zeros
        // Use Array2::zeros((nx, nz))
        Wavefield {
            vx: Array2::zeros((nx, nz)),
            vz: Array2::zeros((nx, nz)),
            sigma_xx: Array2::zeros((nx, nz)),
            sigma_zz: Array2::zeros((nx, nz)),
            sigma_xz: Array2::zeros((nx, nz)),
        }
    }

    pub fn zero(&mut self) {
        // Set all arrays back to zero
        // Use .fill(0.0) method on each array
        self.vx.fill(0.0);
        self.vz.fill(0.0);
        self.sigma_xx.fill(0.0);
        self.sigma_zz.fill(0.0);
        self.sigma_xz.fill(0.0);
    }

    pub fn compute_divergence(&self, dx: f64, dz: f64) -> Array2<f64> {
        let (nx, nz) = self.vx.dim();
        let mut div = Array2::<f64>::zeros((nx, nz));

        // Implement with centered differences
        for i in 1..nx - 1 {
            for k in 1..nz - 1 {
                let dvx_dx = (self.vx[[i + 1, k]] - self.vx[[i - 1, k]]) / (2.0 * dx);
                let dvz_dz = (self.vz[[i, k + 1]] - self.vz[[i, k - 1]]) / (2.0 * dz);
                div[[i, k]] = dvx_dx + dvz_dz;
            }
        }

        // Boundaries: use one-sided differences
        for k in 0..nz {
            // Left boundary (i=0)
            div[[0, k]] = (self.vx[[1, k]] - self.vx[[0, k]]) / dx
                + if k + 1 < nz {
                    (self.vz[[0, k + 1]] - self.vz[[0, k]]) / dz
                } else {
                    0.0
                };
            // Right boundary (i=nx-1)
            div[[nx - 1, k]] = (self.vx[[nx - 1, k]] - self.vx[[nx - 2, k]]) / dx
                + if k + 1 < nz {
                    (self.vz[[nx - 1, k + 1]] - self.vz[[nx - 1, k]]) / dz
                } else {
                    0.0
                };
        }
        for i in 0..nx {
            // Top boundary (k=0)
            div[[i, 0]] = if i + 1 < nx {
                (self.vx[[i + 1, 0]] - self.vx[[i, 0]]) / dx
            } else {
                0.0
            } + (self.vz[[i, 1]] - self.vz[[i, 0]]) / dz;
            // Bottom boundary (k=nz-1)
            div[[i, nz - 1]] = if i + 1 < nx {
                (self.vx[[i + 1, nz - 1]] - self.vx[[i, nz - 1]]) / dx
            } else {
                0.0
            } + (self.vz[[i, nz - 1]] - self.vz[[i, nz - 2]]) / dz;
        }
        div
    }

    pub fn compute_curl(&self, dx: f64, dz: f64) -> Array2<f64> {
        let (nx, nz) = self.vx.dim();
        let mut curl = Array2::<f64>::zeros((nx, nz));

        // For i in 1..nx-1 and k in 1..nz-1:
        //   dvx_dz = (vx[[i, k+1]] - vx[[i, k-1]]) / (2.0 * dz)
        //   dvz_dx = (vz[[i+1, k]] - vz[[i-1, k]]) / (2.0 * dx)
        //   curl[[i, k]] = dvx_dz - dvz_dx

        for i in 1..nx - 1 {
            for k in 1..nz - 1 {
                let dvx_dz = (self.vx[[i, k + 1]] - self.vx[[i, k - 1]]) / (2.0 * dz);
                let dvz_dx = (self.vz[[i + 1, k]] - self.vz[[i - 1, k]]) / (2.0 * dx);
                curl[[i, k]] = dvx_dz - dvz_dx;
            }
        }

        // Boundaries: use one-sided differences
        for k in 0..nz {
            // Left boundary (i=0)
            curl[[0, k]] = (self.vx[[0, k.min(nz - 1)]] - self.vx[[0, k.saturating_sub(1)]]) / dz
                - (self.vz[[1, k]] - self.vz[[0, k]]) / dx;
            // Right boundary (i=nx-1)
            curl[[nx - 1, k]] =
                (self.vx[[nx - 1, k.min(nz - 1)]] - self.vx[[nx - 1, k.saturating_sub(1)]]) / dz
                    - (self.vz[[nx - 1, k]] - self.vz[[nx - 2, k]]) / dx;
        }
        for i in 0..nx {
            // Top boundary (k=0)
            curl[[i, 0]] = (self.vx[[i, 1]] - self.vx[[i, 0]]) / dz
                - if i + 1 < nx {
                    (self.vz[[i + 1, 0]] - self.vz[[i, 0]]) / dx
                } else {
                    0.0
                };
            // Bottom boundary (k=nz-1)
            curl[[i, nz - 1]] = (self.vx[[i, nz - 1]] - self.vx[[i, nz - 2]]) / dz
                - if i + 1 < nx {
                    (self.vz[[i + 1, nz - 1]] - self.vz[[i, nz - 1]]) / dx
                } else {
                    0.0
                };
        }

        curl
    }

    pub fn compute_velocity_magnitude(&self) -> Array2<f64> {
        // |v| = sqrt(vx² + vz²)
        // Useful for visualization

        let (nx, nz) = self.vx.dim();
        let mut mag = Array2::<f64>::zeros((nx, nz));
        // For each point: mag[i][k] = sqrt(vx[i][k]² + vz[i][k]²)
        for i in 0..nx {
            for k in 0..nz {
                mag[[i, k]] = (self.vx[[i, k]].powi(2) + self.vz[[i, k]].powi(2)).sqrt();
            }
        }

        mag
    }
}

use ndarray::Array2;

pub struct MaterialProperties {
    // Primary storage at (i, k) positions - normal stress locations
    pub lambda: Array2<f64>,
    pub mu: Array2<f64>,
    pub rho: Array2<f64>,

    // Pre-computed averaged values at staggered positions
    // These are what you'll actually use in the time-stepping
    pub rho_vx: Array2<f64>, // Density at vx positions (i+1/2, k)
    pub rho_vz: Array2<f64>, // Density at vz positions (i, k+1/2)
    pub mu_xz: Array2<f64>,  // Mu at sigma_xz positions (i+1/2, k+1/2)
}

impl MaterialProperties {
    pub fn new(vp: Array2<f64>, vs: Array2<f64>, rho: Array2<f64>) -> Self {
        // Get dimensions
        let (nx, nz) = vp.dim();

        // Convert to lambda and mu

        // Create empty arrays for lambda and mu
        let mut lambda = Array2::<f64>::zeros((nx, nz));
        let mut mu = Array2::<f64>::zeros((nx, nz));

        // Fill lambda and mu arrays using the formulas
        for i in 0..nx {
            for k in 0..nz {
                let vp_ik = vp[[i, k]];
                let vs_ik = vs[[i, k]];
                let rho_ik = rho[[i, k]];

                mu[[i, k]] = rho_ik * vs_ik * vs_ik;
                lambda[[i, k]] = rho_ik * vp_ik * vp_ik - 2.0 * mu[[i, k]];
            }
        }

        // Pre-compute staggered values
        let rho_vx = Self::average_to_vx(&rho);
        let rho_vz = Self::average_to_vz(&rho);
        let mu_xz = Self::average_to_xz(&mu);

        Self {
            lambda,
            mu,
            rho,
            rho_vx,
            rho_vz,
            mu_xz,
        }
    }

    fn harmonic_mean(a: f64, b: f64) -> f64 {
        // Handle edge case: what if a or b are both zero?
        if a == 0.0 && b == 0.0 {
            0.0
        } else {
            2.0 * a * b / (a + b)
        }
    }

    fn average_to_vx(rho: &Array2<f64>) -> Array2<f64> {
        // Average to vx positions (i+1/2, k)
        let (nx, nz) = rho.dim();
        let mut result = Array2::<f64>::zeros((nx, nz));

        // For each (i, k):
        //   If i < nx-1: result[[i,k]] = harmonic_mean(rho[[i,k]], rho[[i+1,k]])
        //   Else: result[[i,k]] = rho[[i,k]] (boundary case)
        for i in 0..nx {
            for k in 0..nz {
                if i < nx - 1 {
                    result[[i, k]] = Self::harmonic_mean(rho[[i, k]], rho[[i + 1, k]]);
                } else {
                    result[[i, k]] = rho[[i, k]];
                }
            }
        }
        result
    }

    fn average_to_vz(rho: &Array2<f64>) -> Array2<f64> {
        // Average to vz positions (i, k+1/2)
        // TODO: Similar to average_to_vx but in z direction
    }

    fn average_to_xz(mu: &Array2<f64>) -> Array2<f64> {
        // Average to sigma_xz positions (i+1/2, k+1/2)
        // Need to average over 4 points
        let (nx, nz) = mu.dim();
        let mut result = Array2::<f64>::zeros((nx, nz));

        // TODO: Implement
        // For each (i, k):
        //   If i < nx-1 AND k < nz-1:
        //     Average mu from 4 corners: (i,k), (i+1,k), (i,k+1), (i+1,k+1)
        //     You can do harmonic mean of 4 values, or arithmetic mean is also okay
        //   Else: handle boundaries

        result
    }
}

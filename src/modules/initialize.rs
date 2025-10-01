use ndarray::prelude::{Array1, Array2};
use crate::modules::utils::spectral_field;
use crate::modules::consts::{PI,GRAV_ACCEL};


pub fn set_frequency(fmin: f64, fmax: f64, om: usize) -> (Array1<f64>, f64) {
    let dlnf: f64  = (fmax.ln() - fmin.ln())/ (om-1) as f64;

    let aux = Array1::from_shape_fn(om, |o| {
        (fmin.ln() + o as f64 *dlnf).exp()
    });
    (aux,dlnf)
}


pub fn wavenumber(istart: usize,
    iend  : usize,
    om    : usize,
    frequency: &Array2<f64>,
    water_depth: f64,
    surface_tension: f64,
    water_density: f64) -> Array2<f64> {
    
    let pi  = std::f64::consts::PI;
    
    
    let frequency_nondim = 2. * pi * frequency * (water_depth/GRAV_ACCEL).sqrt();
    let mut k: Array2<f64> = frequency_nondim.mapv(|x| x.powi(2));
    let surface_tension_nondim = surface_tension / (GRAV_ACCEL * water_density * water_depth.powi(2));

    for i in istart..iend{
        for o in 0..om{
        let mut count: usize = 0;
        let mut dk: f64 = 0.;
    
            // while dk.abs()>tol {
            loop {
                let t = (k[[i,o]]).tanh();
                
                dk = -(frequency_nondim[[i,o]].powi(2) -
                       k[[i,o]] * t * (1.+ surface_tension_nondim * (k[[i,o]]).powi(2)))
                       / ( 3. * surface_tension_nondim * (k[[i,o]]).powi(2) * t
                       + t
                       + k[[i,o]] * (1.+surface_tension_nondim*(k[[i,o]]).powi(2)) 
                       * (1.-t.powi(2)));
                k[[i,o]] -= dk;
                
                count += 1;    
                if count >= 1000 {
                    break;  // Escape if stuck
                }
                }
            k[[i,o]] = k[[i,o]]/water_depth;
            }
        }
    k
}


pub fn set_ustar_initial(wspd: & Array1::<f64>,
                        istart: usize,
                        iend:   usize) -> Array1::<f64>{

    let nx = iend - istart;
    let mut cd    = Array1::<f64>::ones(nx) * 1.2e-3;
    let mut ustar = Array1::<f64>::zeros(nx);

    for i in istart..iend{
        if wspd[i] > 11. {
            cd[i] = (0.49 + 0.065 * wspd[i]) * 1e-3;
        }
    }

    for i in istart..iend{
        ustar[i] = (cd[i]).sqrt() * wspd[i];
    }
    ustar
}


pub fn compute_velocities_and_adimensional_depth(
    om    : usize,
    istart: usize,
    iend  : usize,
    nghost: usize,
    sfct  : f64,
    rhow: f64,
    f: &Array1<f64>,
    k: &Array2<f64>,
    depth: &Array1<f64>,
    ) -> (Array2<f64>, Array2<f64>, Array2<f64>){
    
    let nx = iend - istart + 2 * nghost;
    let mut cp0 = spectral_field(nx,om);
    let mut cg0 = spectral_field(nx,om);
    let mut kd = spectral_field(nx,om);
    let g: f64 = 9.8;

    let twopi  = std::f64::consts::PI*2.0;

    for i in istart-1..iend+1 {
        for o in 0..om {
            kd[[i,o]]  = k[[i,o]] * depth[i];
            cp0[[i,o]] = twopi * f[o] / k[[i,o]];
            cg0[[i,o]] = cp0[[i,o]]*(0.5+k[[i,o]]*depth[i]/(2.*kd[[i,o]]).sinh()
                                    +sfct*k[[i,o]]*k[[i,o]]/(rhow*GRAV_ACCEL+sfct*k[[i,o]]*k[[i,o]]));
        }
    }
    return (kd, cp0, cg0)
}




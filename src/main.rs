mod modules;

use ndarray::prelude::{Array1, Array2, Axis,array};
use crate::modules::initialize as init;
use crate::modules::source_functions as sf;
use crate::modules::utils::ArrayExtend;
use crate::modules::utils as ut;
use crate::modules::physics as phy;
use crate::modules::diagnostics as diag;
use crate::modules::consts::{PI, GRAV_ACCEL};
use std::time::Instant;


fn main() {

    let fmax: f64 = 20.0;
    let fmin: f64 = 0.1;
    let om: usize = 50;
    let pi  = PI; //std::f64::consts::PI;
    
    let nx  : usize = 11;
    let istart : usize = 0;
    let iend   : usize = nx;

    let surface_tension: f64 = 0.074;
    let water_density: f64 = 1e3;
    let water_depth: f64 = 1e3;
    
    let x0 = ndarray::Array1::<f64>::linspace(0.,nx as f64 -1.,nx);    
    let (f0, dlnf) = init::set_frequency(fmin, fmax, om);
    let depth = Array2::<f64>::ones((nx,om)) * water_depth;

    let mut f = ndarray::Array2::<f64>::zeros((nx,om));
    let mut x = ndarray::Array2::<f64>::zeros((nx,om));

    pub const dt_panic: f32 = 1e-4;

    
    for i in 0..nx{
        for o in 0..om{
                f[[i,o]] = f0[o];
                x[[i,o]] = x0[i];
            }
    }
    
    let k = init::wavenumber(istart, iend, om, &f, 32., surface_tension, water_density);
    let cp = ut::phase_speed( &f, &k);
    let cg = ut::group_speed(surface_tension, &f, &k, &depth, water_density);
    let dk = 2.*pi* &f/&cg;
    
    
    
    let Fk_init = sf::wind_wave_balance(&sf::source_input(0.8, &f, &k, &cp), &f, &k);
    let wspd = Array1::<f64>::ones(nx) *15.;
    
    let start = Instant::now();
    let swh = phy::integrate(
        istart,
        iend,
        om,
        &x,
        &Fk_init,
        &f,
        &k,
        &cp,
        &cg,
        15.,
        60.,
        1.,
        1.,
        42.,
        0.1);
    let duration = start.elapsed();
    
    println!("Function took: {:?}", duration);
    println!("Function took: {} ms", duration.as_millis());
    println!("Function took: {} μs", duration.as_micros());

    println!("{:.2}", swh)

}
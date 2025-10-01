use ndarray::prelude::{Array2};
use crate::modules::utils::*;
use crate::modules::consts::{PI,GRAV_ACCEL};

pub fn source_input(
    wind_speed: f64,
    frequency: &Array2<f64>,
    wavenumber: &Array2<f64>,
    phase_speed:&Array2<f64>,
) -> Array2<f64> {

    let sheltering_coefficient: f64 = 0.11;
    let air_density: f64 = 1.2;
    let water_density: f64 = 1e3;
    let current: f64 = 0.;
    

    let wind_speed_relative = wind_speed - phase_speed - current;
    let mut s_in = sheltering_coefficient * &wind_speed_relative * &wind_speed_relative.abs_val()
                   * frequency * wavenumber ;
    s_in = s_in * air_density / water_density / GRAV_ACCEL * 2. * PI ;
    s_in
}

pub fn source_dissipation(spectrum_k: &Array2<f64>,
                      frequency: &Array2<f64>,
                      k: &Array2<f64>,
                      dk: &Array2<f64>) -> Array2<f64> {
    
    let dissipation_coefficient: f64 =42.;
    let dissipation_power: f64 =2.4;
    let mss_coefficient: f64 =120.;

    let omega = 2.*PI*frequency;
    
    let mut mss = Array2::<f64>::zeros((k.shape()[0], k.shape()[1]));
    if mss_coefficient >0. {
        let mss = mean_squared_slope_long(&spectrum_k, &k, &dk);
    }
    else {
        mss = Array2::<f64>::zeros((k.shape()[0], k.shape()[1]));
    }
    let mss_effect = (1. + mss_coefficient * &mss).mapv(|x| x.powi(2));

    let b_k = saturation_spectrum(&spectrum_k, &k);    
    dissipation_coefficient * omega * mss_effect * b_k.map(|x| x.powf(dissipation_power))
    
}


pub fn wind_wave_balance(source_input: &Array2<f64>,
                     frequency: &Array2<f64>,
                     wavenumber: &Array2<f64>) -> Array2<f64> {

let dissipation_coefficient: f64 = 42.;
let dissipation_power: f64 = 2.4;
let mss: f64 = 0.;
let mss_coefficient: f64 = 120.;

let omega = 2. * PI * frequency;
let mss_effect = (1. + mss_coefficient*mss).powi(2); 
let mut fk = (source_input 
          /(omega * dissipation_coefficient * mss_effect))
          .mapv(|x| x.powf(1. / dissipation_power))
          / wavenumber.mapv(|x| x.powi(4));

for val in fk.iter_mut() {
    if val.is_nan() {
        *val = 0.0;
    }
}

fk  
}


pub fn source_wave_interaction(istart: usize,
                           iend  : usize,
                           om    : usize,
                           spectrum_k: &Array2<f64>,
                           k: &Array2<f64>,
                           dk: &Array2<f64>) -> Array2<f64> {
    let snl_coefficient: f64=1.;
    let mut snl = Array2::<f64>::zeros((k.shape()[0], k.shape()[1]));

    for o in 0..om-1 {
        for i in istart..iend {
            snl[[i,o]] = snl_coefficient * (spectrum_k[[i,o+1]] - spectrum_k[[i,o]])/ dk[[i,o+1]];
        }
    }
    snl
}
    

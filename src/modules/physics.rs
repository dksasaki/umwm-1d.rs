use ndarray::prelude::{Array2};
use crate::modules::source_functions::*;
use crate::modules::utils::ArrayExtend;
use crate::modules::diagnostics as diag;
use crate::modules::consts::{PI};

pub fn advect(istart: usize,
          iend  : usize,
          om    : usize,
          q: &Array2<f64>,
          cp: &Array2<f64>,
          x: &Array2<f64>,
         ) -> Array2<f64> {
    let mut res = Array2::<f64>::zeros((q.shape()[0], q.shape()[1]));

    for i in istart..iend-1 {
        for o in 0..om {
            res[[i+1, o]] = (q[[i+1,o]] - q[[i,o]])/ (x[[i+1,o]] - x[[i,o]]);
            res[[0,o]] = q[[0,o]]/x[[1,o]] - x[[0,o]];
        }
    }
    &res * cp
}




pub fn integrate(
    istart: usize,
    iend:   usize,
    om:     usize,
    x: &Array2<f64>,
    fk_init: &Array2<f64>,
    f: &Array2<f64>,
    k: &Array2<f64>,
    cp: &Array2<f64>,
    cg: &Array2<f64>,
    wspd: f64,
    duration: f64,
    output_interval: f64,
    mss_coefficient: f64,
    snl_coefficient: f64,
    exp_growth_factor: f64) -> Array2<f64> {
    
    let num_time_steps:usize = (duration/output_interval) as usize;


    let dk = 2. * PI * f/cg;
    
    let mut fk = fk_init.clone();
    let mut swh = Array2::<f64>::zeros((num_time_steps, iend));
    
    // let mut count: usize = 0;
    for count in 0..num_time_steps {
        let mut elapsed: f64 = 0.;

        // while count < 4 {
        while elapsed < output_interval {
            let sin = source_input(wspd, &f, &k, &cp);
            let sds = source_dissipation(&fk, &f, &k, &dk);
            let snl = source_wave_interaction(istart, iend, om, &fk, &k, &dk);
            
            let aux = &sin -&sds;
            
            
            let dt = (exp_growth_factor / &aux.abs_val().max_val())
                .min(output_interval - elapsed);

            fk = &fk * (dt * &aux).mapv(|x| x.exp()) +  dt * (&snl - advect(istart, iend, om, &fk, &cg, &x));
            swh.row_mut(count).assign(&diag::significant_wave_height(&fk, &dk));
            
            elapsed += dt;
        }
        // count += 1;
    }
    swh
}
     
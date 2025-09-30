use ndarray::prelude::{Array2, Array1};


pub fn significant_wave_height(spectrum_k: &Array2<f64>, dk: &Array2<f64>) -> Array1<f64>{
    
    let istart: usize = 0;
    let iend  : usize = spectrum_k.shape()[0];
    let om    : usize = spectrum_k.shape()[1];
    
    let mut swh = Array1::<f64>::zeros(iend);

    for i in istart..iend {
        for o in 0..om{
            swh[[i]]+= spectrum_k[[i,o]] * dk[[i,o]];
        }
    }

    4. * swh.mapv(|x| x.sqrt())
        
}
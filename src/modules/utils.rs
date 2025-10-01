use ndarray::prelude::{Array1, Array2, Array3, Axis};
use ndarray::{ArrayBase, Data, Dimension, Array, array};
use crate::modules::consts::{PI,GRAV_ACCEL};
use num_traits::Float;

pub trait ArrayExtend<T,D> {
    fn min_val(&self) -> T;
    fn max_val(&self) -> T;
    fn abs_val(&self) -> Array<T,D>;
    fn exp_val(&self) -> Array<T,D>;
}

impl<S, D, T> ArrayExtend<T,D> for ArrayBase<S, D>
where
    S: Data<Elem = T>,
    D: Dimension,
    T: Float+ Copy //+ num_traits::Signed,
{
    fn min_val(&self) -> T {
        self.iter()
            .copied()
            .reduce(|a, b| if a < b { a } else { b })
            .expect("called min_val on an empty array!")
    }

    fn max_val(&self) -> T {
        self.iter()
            .copied()
            .reduce(|a, b| if a > b { a } else { b })
            .expect("called max_val on an empty array!")
    }

    fn abs_val(&self) -> Array<T,D>{
        self.mapv(|x| x.abs())
    }

    fn exp_val(&self) -> Array<T,D>{
        self.mapv(|x| x.exp())
    }
}




/// Computes sum(spectrum_k * k^2 * dk, axis=1) for 2D k and dk
pub fn weighted_sum_2d(spectrum_k: &Array2<f64>, k: &Array2<f64>, dk: &Array2<f64>) -> ndarray::Array1<f64> {
    let weighted = spectrum_k * &k.mapv(|x| x.powi(2)) * dk;
    weighted.sum_axis(Axis(1)) // sum along columns (axis=1)
}


pub fn cumsum_axis1(arr: &Array2<f64>) -> Array2<f64> {
    let mut out = arr.clone();
    for mut row in out.axis_iter_mut(Axis(0)) { // iterate over rows
        let mut sum = 0.0;
        for elem in row.iter_mut() {
            sum += *elem;
            *elem = sum;
        }
    } 
    out
}

pub fn weighted_cumsum(spectrum_k: &Array2<f64>, k: &Array2<f64>, dk: &Array2<f64>) -> Array2<f64> {
    // Elementwise multiply
    let weighted = spectrum_k * k.mapv(|x| x.powi(2)) * dk;

    // Cumulative sum along columns (axis=1)
    cumsum_axis1(&weighted)
}



pub fn phase_speed(frequency: &Array2<f64>,wavenumber: &Array2<f64>) -> Array2<f64> {
    2. *PI * frequency/wavenumber
}

pub fn group_speed(surface_tension: f64, frequency: &Array2<f64>, wavenumber: &Array2<f64>, water_depth: &Array2<f64>,
water_density: f64) -> Array2<f64> {
    let cp = phase_speed(&frequency, &wavenumber);
    let kd = wavenumber * water_depth;
    let sigma_k2 = surface_tension * wavenumber.mapv(|x| x.powi(2));
    cp * (0.5 + &kd / &kd.mapv(|x| x.sinh()) + &sigma_k2/ (water_density * GRAV_ACCEL * &sigma_k2))
    
}

                   
pub fn mean_squared_slope(spectrum_k: &Array2<f64>, k: &Array2<f64>, dk: &Array2<f64>) -> Array1<f64> {
    weighted_sum_2d(&spectrum_k, &k, &dk)
}

pub fn saturation_spectrum(spectrum_k: &Array2<f64>, k: &Array2<f64>) -> Array2<f64> {
    spectrum_k * k.mapv(|x| x.powi(4))
}

pub fn mean_squared_slope_long(spectrum_k: &Array2<f64>, k: &Array2<f64>, dk: &Array2<f64>) -> Array2<f64> {
    weighted_cumsum(&spectrum_k, &k, &dk)
}
                      
                   
pub fn spectral_field(nx: usize, om: usize) -> Array2<f64> {
    Array2::<f64>::zeros((nx, om))
}

pub fn spectral_field3d(nx: usize, pm: usize, om: usize) -> Array3<f64> {
    Array3::<f64>::zeros((nx, pm, om))
}
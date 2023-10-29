//// Implement derivatives for the Lane-Emden equation

pub fn get_yprime() -> &'static dyn Fn(f64,f64,f64)->f64 {
    fn yprime(_x:f64, _y:f64, z:f64) -> f64 {
        z
    }
    &yprime
}

// fn get_zprime(n:usize) -> &'static dyn Fn(f64,f64,f64)->f64 {
//     fn zprime(x:f64, y:f64, z:f64) -> f64 {
//         -1.0*y.powf(n as f64) - 2.0/x * z
//     }
//     &zprime
// }
pub fn get_zprime(n: f64) -> Box<dyn Fn(f64, f64, f64) -> f64> {
    Box::new(move |x, y, z| {
        -1.0 * y.powf(n) - 2.0 / x * z
    })
}
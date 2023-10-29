//// Rust module that does the solving
/// 
use crate::runge_kutta;
use crate::derivatives;


pub fn solve(
    x_init: f64,
    n: f64,
    h: f64,
    max_iter: u32
) -> (Vec<f64>,Vec<f64>) {
    let mut x_prev: f64 = x_init;
    let mut y_prev: f64 = 1.0;
    let mut z_prev: f64 = 0.0;
    let yprime = derivatives::get_yprime();
    let zprime = derivatives::get_zprime(n);
    let mut n_iter: u32 = 0;
    let mut xs: Vec<f64> = Vec::new();
    let mut ys: Vec<f64> = Vec::new();
    while (y_prev > 0.0) && (n_iter < max_iter) {
        n_iter += 1;
        xs.push(x_prev);
        ys.push(y_prev);
        let (x_next, y_next, z_next) = runge_kutta::get_next_xyz(
            yprime,
            &zprime,
            x_prev,
            y_prev,
            z_prev,
            h
        );
        x_prev = x_next;
        y_prev = y_next;
        z_prev = z_next;
    }
    (xs, ys)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_n0() {
        let (xs, ys) = solve(1e-3, 0.0, 0.01, 1000);
        assert!(xs.len() == ys.len());
    }
}
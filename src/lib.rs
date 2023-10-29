mod runge_kutta;
mod derivatives;
mod solve_poly;

use pyo3::prelude::*;
use pyo3::Python;
use pyo3::types::PyList as PyO3List;
// use std::marker::Tuple;

// /// Formats the sum of two numbers as string.
// #[pyfunction]
// fn sum_as_string(a: usize, b: usize) -> PyResult<String> {
//     Ok((a + b).to_string())
// }


#[pyfunction]
fn solve(
    py: Python,
    x_init:f64,
    n:f64,
    h:f64,
    max_iter:u32
) -> PyResult<(Py<PyO3List>,Py<PyO3List>)> {
    let (xs, ys): (Vec<f64>, Vec<f64>) = solve_poly::solve(x_init, n, h, max_iter);
    let x_list = PyO3List::new(py, xs);
    let y_list = PyO3List::new(py, ys);
    let result: PyResult<(Py<PyO3List>,Py<PyO3List>)> = 
        PyResult::Ok(
            (
                x_list.into_py(py),
                y_list.into_py(py)
            )
        );
    result
}



/// A Python module implemented in Rust.
#[pymodule]
fn polysolver_rust(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(solve, m)?)?;
    Ok(())
}

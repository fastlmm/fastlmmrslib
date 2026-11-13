use pyo3::prelude::*;

/// Return (upper_tail, Davies fault, trace). A failed calculation has
/// upper_tail=None, never a plausible-looking sentinel probability.
#[pyfunction]
#[pyo3(signature = (chi2val, coeffs, dof=None, noncentrality=None, sigma=0.0, lim=1_000_000, acc=1e-8))]
#[allow(clippy::too_many_arguments)]
fn qf(
    py: Python<'_>,
    chi2val: f64,
    coeffs: Vec<f64>,
    dof: Option<Vec<i32>>,
    noncentrality: Option<Vec<f64>>,
    sigma: f64,
    lim: usize,
    acc: f64,
) -> (Option<f64>, i32, Vec<f64>) {
    let dof = dof.unwrap_or_else(|| vec![1; coeffs.len()]);
    let noncentrality = noncentrality.unwrap_or_else(|| vec![0.0; coeffs.len()]);
    let out =
        py.detach(move || super::davies(chi2val, &coeffs, &dof, &noncentrality, sigma, acc, lim));
    (out.upper_tail, out.fault, out.trace.to_vec())
}

#[pymodule]
fn _core(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(qf, module)?)?;
    Ok(())
}

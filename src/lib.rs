//! Experimental translation of Robert Davies' QFC implementation, as maintained
//! by fastlmmclib. Retain NOTICE.Davies and the original third-party notice.
//! Numerical operation order and constants initially follow QFC.cpp. This is
//! not yet a production replacement or a claim of validated platform support.

// Match the legacy constant during differential qualification.
#[allow(clippy::approx_constant)]
const PI: f64 = 3.14159265358979;
const LOG28: f64 = 0.0866;

#[cfg(feature = "python")]
mod python;

#[derive(Clone, Debug, PartialEq)]
pub struct Outcome {
    /// None when Davies did not produce a probability; never fabricate a value.
    pub upper_tail: Option<f64>,
    /// Original Davies diagnostic: 0 success, 1 accuracy/work limit,
    /// 2 roundoff possibly significant, 3 invalid parameters, 4 search limit.
    pub fault: i32,
    pub trace: [f64; 7],
}

#[derive(Clone, Copy, Debug)]
enum SearchFailure {
    ParameterSearch,
    Accuracy,
}

struct Context<'a> {
    coefficients: &'a [f64],
    noncentrality: &'a [f64],
    dof: &'a [i32],
    order: Vec<usize>,
    sorted: bool,
    failed: bool,
    sigma_squared: f64,
    minimum: f64,
    maximum: f64,
    mean: f64,
    statistic: f64,
    integral: f64,
    error_sum: f64,
    count: usize,
    limit: usize,
    trace: [f64; 7],
}

fn square(x: f64) -> f64 {
    x * x
}

fn log1(x: f64, first: bool) -> f64 {
    // A NaN must not enter the convergence loop: NaN != NaN never terminates.
    // Callers turn a non-finite bound into an explicit search failure.
    if !x.is_finite() {
        return (1.0 + x).ln();
    }
    if x.abs() > 0.1 {
        return if first {
            (1.0 + x).ln()
        } else {
            (1.0 + x).ln() - x
        };
    }
    let mut y = x / (2.0 + x);
    let mut term = 2.0 * y * y * y;
    let mut k = 3.0;
    let mut sum = (if first { 2.0 } else { -x }) * y;
    y *= y;
    let mut next = sum + term / k;
    while next != sum {
        k += 2.0;
        term *= y;
        sum = next;
        next = sum + term / k;
    }
    sum
}

impl Context<'_> {
    fn counter(&mut self) -> Result<(), SearchFailure> {
        self.count += 1;
        if self.count > self.limit {
            Err(SearchFailure::ParameterSearch)
        } else {
            Ok(())
        }
    }

    fn error_bound(&mut self, mut u: f64) -> Result<(f64, f64), SearchFailure> {
        self.counter()?;
        let mut cutoff = u * self.sigma_squared;
        let mut sum = u * cutoff;
        u *= 2.0;
        for j in (0..self.coefficients.len()).rev() {
            let n = f64::from(self.dof[j]);
            let coefficient = self.coefficients[j];
            let nc = self.noncentrality[j];
            let x = u * coefficient;
            let y = 1.0 - x;
            cutoff += coefficient * (nc / y + n) / y;
            sum = sum + nc * square(x / y) + n * (square(x) / y + log1(-x, false));
        }
        let bound = (-0.5 * sum).exp();
        if !bound.is_finite() || !cutoff.is_finite() {
            return Err(SearchFailure::ParameterSearch);
        }
        Ok((bound, cutoff))
    }

    fn cutoff(&mut self, accuracy: f64, up: &mut f64) -> Result<f64, SearchFailure> {
        let mut u2 = *up;
        let mut u1 = 0.0;
        let mut c1 = self.mean;
        let rb = 2.0 * if u2 > 0.0 { self.maximum } else { self.minimum };
        let mut c2;
        loop {
            let (error, cutoff) = self.error_bound(u2 / (1.0 + u2 * rb))?;
            c2 = cutoff;
            if error <= accuracy {
                break;
            }
            u1 = u2;
            c1 = c2;
            u2 *= 2.0;
        }
        while (c1 - self.mean) / (c2 - self.mean) < 0.9 {
            let u = (u1 + u2) / 2.0;
            let (error, cutoff) = self.error_bound(u / (1.0 + u * rb))?;
            if error > accuracy {
                u1 = u;
                c1 = cutoff;
            } else {
                u2 = u;
                c2 = cutoff;
            }
        }
        *up = u2;
        Ok(c2)
    }

    fn truncation(&mut self, mut u: f64, tau_squared: f64) -> Result<f64, SearchFailure> {
        self.counter()?;
        let mut sum1 = 0.0;
        let mut prod2 = 0.0;
        let mut prod3 = 0.0;
        let mut s = 0.0;
        let sum2 = (self.sigma_squared + tau_squared) * square(u);
        let mut prod1 = 2.0 * sum2;
        u *= 2.0;
        for j in 0..self.coefficients.len() {
            let n = f64::from(self.dof[j]);
            let x = square(u * self.coefficients[j]);
            sum1 += self.noncentrality[j] * x / (1.0 + x);
            if x > 1.0 {
                prod2 += n * x.ln();
                prod3 += n * log1(x, true);
                s += n;
            } else {
                prod1 += n * log1(x, true);
            }
        }
        sum1 *= 0.5;
        prod2 += prod1;
        prod3 += prod1;
        let x = (-sum1 - 0.25 * prod2).exp() / PI;
        let y = (-sum1 - 0.25 * prod3).exp() / PI;
        let mut error1 = if s == 0.0 { 1.0 } else { x * 2.0 / s };
        let mut error2 = if prod3 > 1.0 { 2.5 * y } else { 1.0 };
        if error2 < error1 {
            error1 = error2;
        }
        let x = 0.5 * sum2;
        error2 = if x <= y { 1.0 } else { y / x };
        let bound = if error1 < error2 { error1 } else { error2 };
        if !bound.is_finite() {
            return Err(SearchFailure::Accuracy);
        }
        Ok(bound)
    }

    fn find_u(&mut self, ut: &mut f64, accuracy: f64) -> Result<(), SearchFailure> {
        let mut u = *ut / 4.0;
        if self.truncation(u, 0.0)? > accuracy {
            u = *ut;
            while self.truncation(u, 0.0)? > accuracy {
                *ut *= 4.0;
                u = *ut;
            }
        } else {
            *ut = u;
            u /= 4.0;
            while self.truncation(u, 0.0)? <= accuracy {
                *ut = u;
                u /= 4.0;
            }
        }
        for divisor in [2.0, 1.4, 1.2, 1.1] {
            u = *ut / divisor;
            if self.truncation(u, 0.0)? <= accuracy {
                *ut = u;
            }
        }
        Ok(())
    }

    fn integrate(&mut self, terms: usize, interval: f64, tau_squared: f64, main: bool) {
        let in_pi = interval / PI;
        for k in (0..=terms).rev() {
            let u = (k as f64 + 0.5) * interval;
            let mut sum1 = -2.0 * u * self.statistic;
            let mut sum2 = sum1.abs();
            let mut sum3 = -0.5 * self.sigma_squared * square(u);
            for j in (0..self.coefficients.len()).rev() {
                let n = f64::from(self.dof[j]);
                let x = 2.0 * self.coefficients[j] * u;
                let mut y = square(x);
                sum3 -= 0.25 * n * log1(y, true);
                y = self.noncentrality[j] * x / (1.0 + y);
                let z = n * x.atan() + y;
                sum1 += z;
                sum2 += z.abs();
                sum3 -= 0.5 * x * y;
            }
            let mut x = in_pi * sum3.exp() / u;
            if !main {
                x *= 1.0 - (-0.5 * tau_squared * square(u)).exp();
            }
            sum1 = (0.5 * sum1).sin() * x;
            sum2 *= 0.5 * x;
            self.integral += sum1;
            self.error_sum += sum2;
        }
    }

    fn convergence_factor(&mut self, x: f64) -> Result<f64, SearchFailure> {
        self.counter()?;
        if !self.sorted {
            let coefficients = self.coefficients;
            self.order
                .sort_by(|&a, &b| coefficients[b].abs().total_cmp(&coefficients[a].abs()));
            self.sorted = true;
        }
        let mut ax = x.abs();
        let sign = if x > 0.0 { 1.0 } else { -1.0 };
        let mut sum = 0.0;
        for j in (0..self.coefficients.len()).rev() {
            let t = self.order[j];
            if self.coefficients[t] * sign > 0.0 {
                let coefficient = self.coefficients[t].abs();
                let ax1 = ax - coefficient * (f64::from(self.dof[t]) + self.noncentrality[t]);
                let ax2 = coefficient / LOG28;
                if ax1 > ax2 {
                    ax = ax1;
                } else {
                    if ax > ax2 {
                        ax = ax2;
                    }
                    sum = (ax - ax1) / coefficient;
                    for k in (0..j).rev() {
                        let index = self.order[k];
                        sum += f64::from(self.dof[index]) + self.noncentrality[index];
                    }
                    break;
                }
            }
        }
        if sum > 100.0 {
            self.failed = true;
            Ok(1.0)
        } else {
            Ok(2.0_f64.powf(sum / 4.0) / (PI * square(ax)))
        }
    }

    fn run(&mut self, sigma: f64, accuracy: f64) -> Result<(Option<f64>, i32), SearchFailure> {
        let mut sd = self.sigma_squared;
        for j in 0..self.coefficients.len() {
            let coefficient = self.coefficients[j];
            let nc = self.noncentrality[j];
            let n = f64::from(self.dof[j]);
            sd += square(coefficient) * (2.0 * n + 4.0 * nc);
            self.mean += coefficient * (n + nc);
            if self.maximum < coefficient {
                self.maximum = coefficient;
            } else if self.minimum > coefficient {
                self.minimum = coefficient;
            }
        }
        if !sd.is_finite() || !self.mean.is_finite() {
            return Ok((None, 3));
        }
        if sd == 0.0 {
            // Retain QFC's degenerate-boundary convention during this initial
            // translation: at statistic == 0 it returns a p-value of one.
            return Ok((Some(if self.statistic > 0.0 { 0.0 } else { 1.0 }), 0));
        }
        if self.minimum == 0.0 && self.maximum == 0.0 && sigma == 0.0 {
            return Ok((None, 3));
        }
        // Exact support bounds also apply to zero-degree noncentral terms,
        // whose atom at zero can prevent the continuous-tail search converging.
        if sigma == 0.0 {
            if self.maximum <= 0.0 && self.statistic > 0.0 {
                return Ok((Some(0.0), 0));
            }
            if self.minimum >= 0.0 && self.statistic < 0.0 {
                return Ok((Some(1.0), 0));
            }
        }
        sd = sd.sqrt();
        let maximum_abs = if self.maximum < -self.minimum {
            -self.minimum
        } else {
            self.maximum
        };
        let mut ut = 16.0 / sd;
        let mut up = 4.5 / sd;
        let mut un = -up;
        let mut acc = accuracy;
        let mut remaining = self.limit as f64;
        self.find_u(&mut ut, 0.5 * acc)?;
        if self.statistic != 0.0 && maximum_abs > 0.07 * sd {
            let tau = 0.25 * acc / self.convergence_factor(self.statistic)?;
            if self.failed {
                self.failed = false;
            } else if self.truncation(ut, tau)? < 0.2 * acc {
                self.sigma_squared += tau;
                self.find_u(&mut ut, 0.25 * acc)?;
                self.trace[5] = tau.sqrt();
            }
        }
        self.trace[4] = ut;
        acc *= 0.5;
        loop {
            let d1 = self.cutoff(acc, &mut up)? - self.statistic;
            if d1 < 0.0 {
                return Ok((Some(0.0), 0));
            }
            let d2 = self.statistic - self.cutoff(acc, &mut un)?;
            if d2 < 0.0 {
                return Ok((Some(1.0), 0));
            }
            let interval = 2.0 * PI / if d1 > d2 { d1 } else { d2 };
            let required = ut / interval;
            let auxiliary = 3.0 / acc.sqrt();
            if required > auxiliary * 1.5 {
                if auxiliary > remaining {
                    return Ok((None, 1));
                }
                let terms = (auxiliary + 0.5).floor() as usize;
                let interval1 = ut / terms as f64;
                let x = 2.0 * PI / interval1;
                if x > self.statistic.abs() {
                    let tau = 0.33 * acc
                        / (1.1
                            * (self.convergence_factor(self.statistic - x)?
                                + self.convergence_factor(self.statistic + x)?));
                    if !self.failed {
                        acc *= 0.67;
                        self.integrate(terms, interval1, tau, false);
                        remaining -= auxiliary;
                        self.sigma_squared += tau;
                        self.trace[2] += 1.0;
                        self.trace[1] += terms as f64 + 1.0;
                        self.find_u(&mut ut, 0.25 * acc)?;
                        acc *= 0.75;
                        continue;
                    }
                }
            }
            self.trace[3] = interval;
            if required > remaining {
                return Ok((None, 1));
            }
            self.integrate((required + 0.5).floor() as usize, interval, 0.0, true);
            self.trace[2] += 1.0;
            self.trace[1] += (required + 0.5).floor() + 1.0;
            self.trace[0] = self.error_sum;
            let x = self.error_sum + accuracy / 10.0;
            let fault = if [1.0, 2.0, 4.0, 8.0]
                .iter()
                .any(|r| r * x == r * self.error_sum)
            {
                2
            } else {
                0
            };
            return Ok((Some(1.0 - (0.5 - self.integral)), fault));
        }
    }
}

/// Evaluate the upper tail of a quadratic form using the complete Davies
/// algorithm. Failure outcomes carry the original fault and diagnostic trace;
/// their probability is None, unlike the legacy wrapper's misleading 2.0.
///
/// This experimental API validates shape, finiteness, nonnegative degrees of
/// freedom/noncentrality, accuracy in (0, 1), and the positive work limit.
pub fn davies(
    statistic: f64,
    coefficients: &[f64],
    dof: &[i32],
    noncentrality: &[f64],
    sigma: f64,
    accuracy: f64,
    max_terms: usize,
) -> Outcome {
    let invalid = coefficients.len() != dof.len()
        || coefficients.len() != noncentrality.len()
        || !statistic.is_finite()
        || !sigma.is_finite()
        || !accuracy.is_finite()
        || accuracy <= 0.0
        || accuracy >= 1.0
        || max_terms == 0
        || max_terms > i32::MAX as usize
        || coefficients.iter().any(|x| !x.is_finite())
        || noncentrality.iter().any(|x| !x.is_finite() || *x < 0.0)
        || dof.iter().any(|x| *x < 0);
    if invalid {
        return Outcome {
            upper_tail: None,
            fault: 3,
            trace: [0.0; 7],
        };
    }
    let mut ctx = Context {
        coefficients,
        dof,
        noncentrality,
        order: (0..coefficients.len()).collect(),
        sorted: false,
        failed: false,
        sigma_squared: square(sigma),
        minimum: 0.0,
        maximum: 0.0,
        mean: 0.0,
        statistic,
        integral: 0.0,
        error_sum: 0.0,
        count: 0,
        limit: max_terms,
        trace: [0.0; 7],
    };
    let (probability, fault) = match ctx.run(sigma, accuracy) {
        Ok(result) => result,
        Err(SearchFailure::ParameterSearch) => (None, 4),
        Err(SearchFailure::Accuracy) => (None, 1),
    };
    ctx.trace[6] = ctx.count as f64;
    Outcome {
        upper_tail: probability,
        fault,
        trace: ctx.trace,
    }
}

#[cfg(test)]
mod tests {
    use super::davies;
    use rayon::prelude::*;

    #[test]
    fn nonfinite_search_bound_terminates() {
        // A zero-degree noncentral term can exhaust the truncation search.
        // NaN in the log-series convergence loop previously never terminated.
        let out = davies(
            11.498286033063403,
            &[1.4256991860196153],
            &[0],
            &[1.2791211646129712],
            0.0,
            1e-7,
            1_000_000,
        );
        assert_eq!(out.fault, 1);
        assert_eq!(out.upper_tail, None);
    }

    #[test]
    fn analytic_exponential_cases() {
        for scale in [0.5, 2.0] {
            for x in [0.0, 0.5, 3.0, 10.0] {
                let out = davies(x, &[scale], &[2], &[0.0], 0.0, 1e-8, 1_000_000);
                assert_eq!(out.fault, 0, "{out:?}");
                let expected = (-x / (2.0 * scale)).exp();
                assert!(
                    (out.upper_tail.unwrap() - expected).abs() <= 1e-8,
                    "{out:?}"
                );
            }
        }
    }

    #[test]
    fn analytic_unequal_unit_dof_cases() {
        for a in [1e-3, 1e-8] {
            let out = davies(
                20.0,
                &[a, a, 1.0, 1.0],
                &[1; 4],
                &[0.0; 4],
                0.0,
                1e-7,
                1_000_000,
            );
            assert_eq!(out.fault, 0);
            let expected = (-10.0_f64).exp() / (1.0 - a);
            assert!((out.upper_tail.unwrap() - expected).abs() < 1e-7, "{out:?}");
        }
    }

    #[test]
    fn invalid_parameters_have_no_probability() {
        for out in [
            davies(1.0, &[1.0], &[], &[0.0], 0.0, 1e-7, 100),
            davies(f64::NAN, &[1.0], &[1], &[0.0], 0.0, 1e-7, 100),
            davies(1.0, &[1.0], &[-1], &[0.0], 0.0, 1e-7, 100),
            davies(1.0, &[1.0], &[1], &[-1.0], 0.0, 1e-7, 100),
            davies(1.0, &[1.0], &[1], &[0.0], 0.0, 0.0, 100),
        ] {
            assert_eq!(out.fault, 3);
            assert_eq!(out.upper_tail, None);
        }
    }

    #[test]
    fn exhausted_search_has_no_probability() {
        let out = davies(1.0, &[1.0], &[1], &[0.0], 0.0, 1e-8, 1);
        assert_eq!(out.fault, 4);
        assert_eq!(out.upper_tail, None);
        assert!(out.trace[6] > 1.0);
    }

    #[test]
    fn calls_are_reentrant() {
        fn evaluate(index: usize) -> super::Outcome {
            let scale = 0.5 + (index % 7) as f64;
            let weights = [0.1 * scale, -0.3, 0.5, 0.8];
            let limit = if index % 6 == 5 { 1 } else { 1_000_000 };
            davies(
                2.0 + scale,
                &weights,
                &[1, 2, 3, 4],
                &[0.0, 0.2 * scale, 0.0, 1.0],
                0.5,
                1e-8,
                limit,
            )
        }
        let expected: Vec<_> = (0..48).map(evaluate).collect();
        assert!(expected.iter().any(|result| result.fault == 0));
        assert!(expected.iter().any(|result| result.fault == 4));
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(4)
            .build()
            .unwrap();
        let results: Vec<_> = pool.install(|| (0..48).into_par_iter().map(evaluate).collect());
        assert_eq!(results, expected);
    }
}

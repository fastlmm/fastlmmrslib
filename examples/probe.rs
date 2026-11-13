//! Line protocol for differential tests: n, statistic, sigma, limit, accuracy,
//! followed by n coefficients, n degrees of freedom, n noncentralities.
use std::error::Error;
use std::io::{self, BufRead};

fn next<T: std::str::FromStr>(
    words: &mut std::str::SplitWhitespace<'_>,
) -> Result<T, Box<dyn Error>>
where
    T::Err: Error + 'static,
{
    Ok(words.next().ok_or("missing field")?.parse()?)
}

fn main() -> Result<(), Box<dyn Error>> {
    for (index, line) in io::stdin().lock().lines().enumerate() {
        let line = line?;
        let mut words = line.split_whitespace();
        let n: usize = next(&mut words)?;
        let statistic = next(&mut words)?;
        let sigma = next(&mut words)?;
        let limit = next(&mut words)?;
        let accuracy = next(&mut words)?;
        let coefficients: Vec<f64> = (0..n).map(|_| next(&mut words)).collect::<Result<_, _>>()?;
        let dof: Vec<i32> = (0..n).map(|_| next(&mut words)).collect::<Result<_, _>>()?;
        let noncentrality: Vec<f64> = (0..n).map(|_| next(&mut words)).collect::<Result<_, _>>()?;
        if words.next().is_some() {
            return Err("unexpected extra field".into());
        }
        eprintln!("case {index}");
        let out = fastlmmrslib::davies(
            statistic,
            &coefficients,
            &dof,
            &noncentrality,
            sigma,
            accuracy,
            limit,
        );
        print!("{}", out.fault);
        if let Some(p) = out.upper_tail {
            print!(" {p:.17e}");
        } else {
            print!(" none");
        }
        for value in out.trace {
            print!(" {value:.17e}");
        }
        println!();
    }
    Ok(())
}

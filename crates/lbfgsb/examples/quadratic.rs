use ariadne_lbfgsb::{Bounds, Options, Solver};
use std::convert::Infallible;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut x = [3.0, -4.0];
    let lower = [-10.0; 2];
    let upper = [10.0; 2];
    let bounds = Bounds::new(&lower, &upper, x.len())?;
    let options = Options::new()
        .with_relative_function_tolerance(1e-12)?
        .with_projected_gradient_tolerance(1e-8)?;
    let mut solver = Solver::new(options);

    let report = solver.minimize(&mut x, bounds, |x, gradient| {
        gradient[0] = 2.0 * (x[0] - 1.0);
        gradient[1] = 2.0 * (x[1] + 2.0);
        Ok::<_, Infallible>((x[0] - 1.0).powi(2) + (x[1] + 2.0).powi(2))
    })?;

    println!(
        "x={x:?}, f={}, termination={:?}",
        report.value, report.termination
    );
    Ok(())
}

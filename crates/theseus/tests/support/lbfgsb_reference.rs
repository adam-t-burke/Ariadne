//! Used only in the baseline worktree by scripts/compare-optimizers.py.

use super::{Case, Outcome};
use std::convert::Infallible;

pub(super) fn solve(case: &Case) -> Outcome {
    use ariadne_lbfgsb::{Backend, Bounds, Control, Options, Solver};
    let options = Options::new()
        .with_backend(Backend::Faer)
        .unwrap()
        .with_history_size(case.history)
        .unwrap()
        .with_projected_gradient_tolerance(case.gradient_tolerance)
        .unwrap()
        .with_relative_function_tolerance(case.cost_tolerance)
        .unwrap();
    let mut x = case.start.clone();
    let report = Solver::new(options)
        .minimize_with_callback(
            &mut x,
            Bounds::new(&case.lower, &case.upper, case.start.len()).unwrap(),
            |x, g| Ok::<_, Infallible>(case.evaluate(x, g)),
            |iteration| {
                if case.gradient_tolerance == 0.0
                    && (iteration.stats.evaluations >= case.evaluation_limit
                        || iteration.projected_gradient_norm
                            <= 1e-10 * (1.0 + iteration.value.abs()))
                {
                    Control::Stop
                } else {
                    Control::Continue
                }
            },
        )
        .unwrap();
    Outcome {
        x,
        value: report.value,
        projected_gradient: report.projected_gradient_norm,
        iterations: report.stats.iterations,
        evaluations: report.stats.evaluations,
    }
}

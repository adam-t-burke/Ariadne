use super::{Case, Outcome};
use std::convert::Infallible;

impl basin::CostFunction for &Case {
    type Param = Vec<f64>;
    type Output = f64;
    type Error = Infallible;
    fn cost(&self, x: &Vec<f64>) -> Result<f64, Infallible> {
        Ok(self.evaluate(x, &mut vec![0.0; x.len()]))
    }
}
impl basin::Gradient for &Case {
    type Gradient = Vec<f64>;
    fn gradient(&self, x: &Vec<f64>) -> Result<Vec<f64>, Infallible> {
        Ok(self.cost_and_gradient(x)?.1)
    }
    fn cost_and_gradient(&self, x: &Vec<f64>) -> Result<(f64, Vec<f64>), Infallible> {
        let mut gradient = vec![0.0; x.len()];
        let value = self.evaluate(x, &mut gradient);
        Ok((value, gradient))
    }
}
impl basin::BoxConstraints for &Case {
    fn lower(&self) -> &Vec<f64> {
        &self.lower
    }
    fn upper(&self) -> &Vec<f64> {
        &self.upper
    }
}

fn projected_gradient(x: &[f64], g: &[f64], lower: &[f64], upper: &[f64]) -> f64 {
    x.iter()
        .zip(g)
        .enumerate()
        .map(|(i, (&x, &g))| {
            if g < 0.0 {
                g.max(x - upper[i]).abs()
            } else {
                g.min(x - lower[i]).abs()
            }
        })
        .fold(0.0, f64::max)
}

pub(super) fn solve(case: &Case) -> Outcome {
    use basin::{Executor, GradientState, LbfgsState, Lbfgsb, State, TerminationReason};
    let solver = Lbfgsb::new().with_absolute_projected_gradient_tolerance(
        (case.gradient_tolerance > 0.0).then_some(case.gradient_tolerance),
    );
    let mut previous: Option<f64> = None;
    let cost_tolerance = case.cost_tolerance;
    let custom_stop = case.gradient_tolerance == 0.0;
    let evaluation_limit = case.evaluation_limit;
    let lower = case.lower.clone();
    let upper = case.upper.clone();
    let result = Executor::new(
        case,
        solver,
        LbfgsState::new(case.start.clone(), case.history),
    )
    .max_iter(2000)
    .stop_when(move |state: &LbfgsState<Vec<f64>>| {
        let pg = projected_gradient(state.param(), state.gradient().unwrap(), &lower, &upper);
        let value = state.cost();
        let old = previous.replace(value);
        if state.iter() > 0
            && custom_stop
            && (state.cost_evals() as usize >= evaluation_limit
                || pg <= 1e-10 * (1.0 + value.abs()))
        {
            return Some(TerminationReason::UserRequested);
        }
        if cost_tolerance > 0.0
            && old.is_some_and(|old| {
                old - value <= cost_tolerance * old.abs().max(value.abs()).max(1.0)
            })
        {
            return Some(TerminationReason::RelativeCostTolerance);
        }
        None
    })
    .run()
    .unwrap();
    assert!(
        !result.reason.is_failure(),
        "{}: {:?}",
        case.name,
        result.reason
    );
    let pg = projected_gradient(
        result.param(),
        result.state.gradient().unwrap(),
        &case.lower,
        &case.upper,
    );
    Outcome {
        x: result.param().clone(),
        value: result.cost(),
        projected_gradient: pg,
        iterations: result.iter() as usize,
        evaluations: result.cost_evals() as usize,
    }
}

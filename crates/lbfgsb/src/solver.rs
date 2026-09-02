//! Callback-first solver API backed by the in-repo scalar engine.
#![allow(clippy::needless_range_loop)] // Preserve line-search arithmetic order.

use crate::core::{feasible_maximum_step, DcSearch, SearchConfig, SearchResult};
use crate::kernel::Kernel;
use crate::session::{DirectionError, Session};
#[cfg(feature = "benchmark-instrumentation")]
use crate::BenchmarkTimings;
use crate::{
    projected_gradient_norm, Backend, Bounds, Control, Convergence, EvaluationError, Failure,
    Iteration, Options, Report, SolveError, Stats, StopReason, Termination, Workspace,
};
#[cfg(feature = "benchmark-instrumentation")]
use std::time::Instant;

/// Stateful objective and accepted-iteration adapter.
///
/// Keeping both hooks on one object lets applications share mutable evaluation
/// state without reference counting or interior mutability.
pub trait SolveAdapter {
    /// Application error returned by either callback.
    type Error;

    /// Computes the objective value and overwrites every gradient entry.
    fn value_and_gradient(&mut self, x: &[f64], gradient: &mut [f64]) -> Result<f64, Self::Error>;

    /// Observes an accepted major iteration and optionally stops the solve.
    fn accepted_iteration(&mut self, _iteration: Iteration<'_>) -> Result<Control, Self::Error> {
        Ok(Control::Continue)
    }
}

struct ClosureAdapter<F, C> {
    value_and_gradient: F,
    callback: C,
}

impl<E, F, C> SolveAdapter for ClosureAdapter<F, C>
where
    F: FnMut(&[f64], &mut [f64]) -> Result<f64, E>,
    C: FnMut(Iteration<'_>) -> Result<Control, E>,
{
    type Error = E;

    fn value_and_gradient(&mut self, x: &[f64], gradient: &mut [f64]) -> Result<f64, Self::Error> {
        (self.value_and_gradient)(x, gradient)
    }

    fn accepted_iteration(&mut self, iteration: Iteration<'_>) -> Result<Control, Self::Error> {
        (self.callback)(iteration)
    }
}

/// An L-BFGS-B solver with reusable owned workspace.
#[derive(Debug)]
pub struct Solver {
    options: Options,
    workspace: Workspace,
    kernel: Kernel,
    #[cfg(feature = "benchmark-instrumentation")]
    instrumentation: BenchmarkTimings,
}

impl Default for Solver {
    fn default() -> Self {
        Self::new(Options::default())
    }
}

impl Solver {
    /// Creates a solver with empty reusable workspace.
    pub fn new(options: Options) -> Self {
        Self {
            options,
            workspace: Workspace::new(),
            kernel: Kernel::selected(options.backend()),
            #[cfg(feature = "benchmark-instrumentation")]
            instrumentation: BenchmarkTimings::default(),
        }
    }

    /// Creates a solver using an existing reusable workspace.
    pub fn with_workspace(options: Options, workspace: Workspace) -> Self {
        Self {
            options,
            workspace,
            kernel: Kernel::selected(options.backend()),
            #[cfg(feature = "benchmark-instrumentation")]
            instrumentation: BenchmarkTimings::default(),
        }
    }

    /// Returns the active solver configuration.
    pub fn options(&self) -> Options {
        self.options
    }

    /// Returns the solver's reusable workspace.
    pub fn workspace(&self) -> &Workspace {
        &self.workspace
    }

    /// Returns mutable access to the solver's reusable workspace.
    pub fn workspace_mut(&mut self) -> &mut Workspace {
        &mut self.workspace
    }

    /// Returns the backend resolved for this solver.
    pub fn resolved_backend(&self) -> Backend {
        self.kernel.public()
    }

    #[cfg(feature = "benchmark-instrumentation")]
    /// Returns timings from the most recent attempted solve.
    ///
    /// `total_nanoseconds` and `residual_nanoseconds` are finalized only when
    /// the solve returns `Ok`; phase counters may be partial after an error.
    pub fn benchmark_timings(&self) -> BenchmarkTimings {
        self.instrumentation
    }

    /// Minimizes an objective, updating `x` in place.
    ///
    /// The callback must overwrite every entry in `gradient` before returning.
    pub fn minimize<E, F>(
        &mut self,
        x: &mut [f64],
        bounds: Bounds<'_>,
        value_and_gradient: F,
    ) -> Result<Report, SolveError<E>>
    where
        F: FnMut(&[f64], &mut [f64]) -> Result<f64, E>,
    {
        self.minimize_with_callback(x, bounds, value_and_gradient, |_| Control::Continue)
    }

    /// Minimizes an objective and reports each accepted major iteration.
    pub fn minimize_with_callback<E, F, C>(
        &mut self,
        x: &mut [f64],
        bounds: Bounds<'_>,
        value_and_gradient: F,
        mut callback: C,
    ) -> Result<Report, SolveError<E>>
    where
        F: FnMut(&[f64], &mut [f64]) -> Result<f64, E>,
        C: FnMut(Iteration<'_>) -> Control,
    {
        self.minimize_with_fallible_callback(x, bounds, value_and_gradient, |iteration| {
            Ok(callback(iteration))
        })
    }

    /// Minimizes an objective with a fallible accepted-iteration callback.
    pub fn minimize_with_fallible_callback<E, F, C>(
        &mut self,
        x: &mut [f64],
        bounds: Bounds<'_>,
        value_and_gradient: F,
        callback: C,
    ) -> Result<Report, SolveError<E>>
    where
        F: FnMut(&[f64], &mut [f64]) -> Result<f64, E>,
        C: FnMut(Iteration<'_>) -> Result<Control, E>,
    {
        let mut adapter = ClosureAdapter {
            value_and_gradient,
            callback,
        };
        self.minimize_with_adapter(x, bounds, &mut adapter)
    }

    /// Minimizes an objective through one stateful application adapter.
    pub fn minimize_with_adapter<A>(
        &mut self,
        x: &mut [f64],
        bounds: Bounds<'_>,
        adapter: &mut A,
    ) -> Result<Report, SolveError<A::Error>>
    where
        A: SolveAdapter,
    {
        bounds.validate(x.len())?;
        if let Some(index) = x.iter().position(|value| !value.is_finite()) {
            return Err(SolveError::InvalidInitialValue { index });
        }
        for ((xi, &lower), &upper) in x.iter_mut().zip(bounds.lower()).zip(bounds.upper()) {
            *xi = xi.max(lower).min(upper);
        }

        let n = x.len();
        self.workspace
            .prepare(bounds.lower(), bounds.upper(), self.options.history_size())?;
        #[cfg(feature = "benchmark-instrumentation")]
        {
            self.instrumentation = BenchmarkTimings::default();
        }
        #[cfg(feature = "benchmark-instrumentation")]
        let solve_started = Instant::now();
        let mut session = Session::new(n, self.options.history_size(), self.kernel);
        let mut evaluations = 0usize;
        let mut iterations = 0usize;
        let mut line_search_probes = 0usize;
        #[cfg(feature = "benchmark-instrumentation")]
        let phase_start = Instant::now();
        let mut value = evaluate(
            x,
            &mut self.workspace.evaluation_gradient,
            adapter,
            &mut evaluations,
        )?;
        #[cfg(feature = "benchmark-instrumentation")]
        {
            self.instrumentation.evaluation_nanoseconds += elapsed_nanoseconds(phase_start);
        }
        let mut pg = projected_gradient_norm(x, &self.workspace.evaluation_gradient, bounds)
            .map_err(SolveError::Bounds)?;
        if self.options.projected_gradient_tolerance() > 0.0
            && pg <= self.options.projected_gradient_tolerance()
        {
            let report = report(
                value,
                pg,
                solve_stats(
                    &session,
                    x,
                    &self.workspace,
                    iterations,
                    evaluations,
                    line_search_probes,
                ),
                Termination::Converged(Convergence::ProjectedGradient),
            );
            self.workspace.commit_gradient()?;
            #[cfg(feature = "benchmark-instrumentation")]
            self.finish_instrumentation(solve_started);
            return Ok(report);
        }

        let termination = 'solve: loop {
            if iterations >= self.options.max_iterations() {
                break Termination::Stopped(StopReason::MaximumIterations);
            }
            if evaluations >= self.options.max_evaluations() {
                break Termination::Stopped(StopReason::MaximumEvaluations);
            }
            #[cfg(feature = "benchmark-instrumentation")]
            let phase_start = Instant::now();
            let direction = session.direction(x, &mut self.workspace);
            #[cfg(feature = "benchmark-instrumentation")]
            {
                self.instrumentation.direction_nanoseconds += elapsed_nanoseconds(phase_start);
            }
            if let Err(error) = direction {
                break Termination::Failed(match error {
                    DirectionError::NoDescent => Failure::NoDescentDirection,
                    DirectionError::Factorization => Failure::Factorization,
                    DirectionError::Numerical => Failure::Numerical,
                });
            }

            self.workspace.old_x.copy_from_slice(x);
            self.workspace
                .old_gradient
                .copy_from_slice(&self.workspace.evaluation_gradient);
            let old_value = value;
            let gdold = dot(&self.workspace.old_gradient, &self.workspace.direction);
            let maximum_step = feasible_maximum_step(
                &self.workspace.old_x,
                &self.workspace.direction,
                &self.workspace.lower,
                &self.workspace.upper,
                iterations,
            );
            let norm = dot(&self.workspace.direction, &self.workspace.direction).sqrt();
            let boxed = self
                .workspace
                .lower
                .iter()
                .zip(&self.workspace.upper)
                .all(|(&l, &u)| l.is_finite() && u.is_finite());
            let initial_step = if iterations == 0 && !boxed {
                (1.0 / norm).min(maximum_step)
            } else {
                1.0
            };
            let mut search = match DcSearch::start(
                value,
                gdold,
                initial_step,
                SearchConfig::lbfgsb(maximum_step),
            ) {
                Ok(search) => search,
                Err(_) if session.has_history() => {
                    session.reset_history();
                    continue;
                }
                Err(_) => break Termination::Failed(Failure::LineSearch),
            };
            let mut accepted = false;
            for _ in 0..self.options.max_line_search_evaluations() {
                if evaluations >= self.options.max_evaluations() {
                    x.copy_from_slice(&self.workspace.old_x);
                    value = old_value;
                    self.workspace
                        .evaluation_gradient
                        .copy_from_slice(&self.workspace.old_gradient);
                    break 'solve Termination::Stopped(StopReason::MaximumEvaluations);
                }
                let step = search.step();
                if step == 1.0 {
                    x.copy_from_slice(&self.workspace.cauchy);
                } else {
                    for i in 0..n {
                        x[i] = (self.workspace.old_x[i] + step * self.workspace.direction[i])
                            .max(self.workspace.lower[i])
                            .min(self.workspace.upper[i]);
                    }
                }
                #[cfg(feature = "benchmark-instrumentation")]
                let phase_start = Instant::now();
                let trial_value = match evaluate(
                    x,
                    &mut self.workspace.trial_gradient,
                    adapter,
                    &mut evaluations,
                ) {
                    Ok(value) => value,
                    Err(error) => {
                        x.copy_from_slice(&self.workspace.old_x);
                        self.workspace
                            .evaluation_gradient
                            .copy_from_slice(&self.workspace.old_gradient);
                        return Err(error);
                    }
                };
                #[cfg(feature = "benchmark-instrumentation")]
                {
                    self.instrumentation.evaluation_nanoseconds += elapsed_nanoseconds(phase_start);
                }
                line_search_probes += 1;
                let gd = dot(&self.workspace.trial_gradient, &self.workspace.direction);
                match search.advance(trial_value, gd) {
                    SearchResult::Evaluate(_) => {}
                    SearchResult::Converged(_) | SearchResult::Warning { .. } => {
                        value = trial_value;
                        self.workspace
                            .evaluation_gradient
                            .copy_from_slice(&self.workspace.trial_gradient);
                        accepted = true;
                        break;
                    }
                }
            }
            if !accepted {
                x.copy_from_slice(&self.workspace.old_x);
                value = old_value;
                self.workspace
                    .evaluation_gradient
                    .copy_from_slice(&self.workspace.old_gradient);
                if session.has_history() {
                    session.reset_history();
                    continue;
                }
                break Termination::Stopped(StopReason::MaximumLineSearchEvaluations);
            }

            iterations += 1;
            pg = projected_gradient_norm(x, &self.workspace.evaluation_gradient, bounds)
                .map_err(SolveError::Bounds)?;
            let stats = solve_stats(
                &session,
                x,
                &self.workspace,
                iterations,
                evaluations,
                line_search_probes,
            );
            #[cfg(feature = "benchmark-instrumentation")]
            let phase_start = Instant::now();
            let control = adapter
                .accepted_iteration(Iteration {
                    x,
                    value,
                    projected_gradient_norm: pg,
                    stats,
                })
                .map_err(SolveError::Callback)?;
            #[cfg(feature = "benchmark-instrumentation")]
            {
                self.instrumentation.observer_nanoseconds += elapsed_nanoseconds(phase_start);
            }
            if control == Control::Stop {
                break Termination::Stopped(StopReason::User);
            }
            if self.options.projected_gradient_tolerance() > 0.0
                && pg <= self.options.projected_gradient_tolerance()
            {
                break Termination::Converged(Convergence::ProjectedGradient);
            }
            let scale = old_value.abs().max(value.abs()).max(1.0);
            if self.options.relative_function_tolerance() > 0.0
                && old_value - value <= self.options.relative_function_tolerance() * scale
            {
                break Termination::Converged(Convergence::RelativeFunction);
            }
            if evaluations >= self.options.max_evaluations() {
                break Termination::Stopped(StopReason::MaximumEvaluations);
            }
            #[cfg(feature = "benchmark-instrumentation")]
            let phase_start = Instant::now();
            session.update(x, gdold, &mut self.workspace);
            #[cfg(feature = "benchmark-instrumentation")]
            {
                self.instrumentation.history_update_nanoseconds += elapsed_nanoseconds(phase_start);
            }
        };
        let report = report(
            value,
            pg,
            solve_stats(
                &session,
                x,
                &self.workspace,
                iterations,
                evaluations,
                line_search_probes,
            ),
            termination,
        );
        self.workspace.commit_gradient()?;
        #[cfg(feature = "benchmark-instrumentation")]
        self.finish_instrumentation(solve_started);
        Ok(report)
    }

    #[cfg(feature = "benchmark-instrumentation")]
    fn finish_instrumentation(&mut self, started: Instant) {
        self.instrumentation.total_nanoseconds = elapsed_nanoseconds(started);
        let attributed = self
            .instrumentation
            .evaluation_nanoseconds
            .saturating_add(self.instrumentation.direction_nanoseconds)
            .saturating_add(self.instrumentation.observer_nanoseconds)
            .saturating_add(self.instrumentation.history_update_nanoseconds);
        self.instrumentation.residual_nanoseconds = self
            .instrumentation
            .total_nanoseconds
            .saturating_sub(attributed);
    }
}

fn evaluate<A>(
    x: &[f64],
    gradient: &mut [f64],
    adapter: &mut A,
    evaluations: &mut usize,
) -> Result<f64, SolveError<A::Error>>
where
    A: SolveAdapter,
{
    *evaluations += 1;
    let value = adapter
        .value_and_gradient(x, gradient)
        .map_err(SolveError::Objective)?;
    if !value.is_finite() {
        return Err(SolveError::Evaluation(EvaluationError::NonFiniteValue));
    }
    if let Some(index) = gradient.iter().position(|value| !value.is_finite()) {
        return Err(SolveError::Evaluation(EvaluationError::NonFiniteGradient {
            index,
        }));
    }
    Ok(value)
}

fn report(
    value: f64,
    projected_gradient_norm: f64,
    stats: Stats,
    termination: Termination,
) -> Report {
    Report {
        value,
        projected_gradient_norm,
        termination,
        stats,
    }
}

fn solve_stats(
    session: &Session,
    x: &[f64],
    workspace: &Workspace,
    iterations: usize,
    evaluations: usize,
    line_search_probes: usize,
) -> Stats {
    let active_variables = x
        .iter()
        .zip(&workspace.evaluation_gradient)
        .zip(&workspace.lower)
        .zip(&workspace.upper)
        .filter(|(((x, gradient), lower), upper)| {
            **lower == **upper
                || (**x <= **lower && **gradient >= 0.0)
                || (**x >= **upper && **gradient <= 0.0)
        })
        .count();
    Stats {
        iterations,
        evaluations,
        line_search_probes,
        accepted_updates: session.accepted_updates(),
        skipped_updates: session.skipped_updates(),
        history_resets: session.history_resets(),
        active_variables,
        free_variables: x.len() - active_variables,
    }
}

fn dot(left: &[f64], right: &[f64]) -> f64 {
    left.iter().zip(right).map(|(&a, &b)| a * b).sum()
}

#[cfg(feature = "benchmark-instrumentation")]
fn elapsed_nanoseconds(start: Instant) -> u64 {
    start.elapsed().as_nanos().min(u64::MAX as u128) as u64
}

#[cfg(test)]
mod adapter_tests {
    use super::*;

    struct FailingObserver {
        evaluations: usize,
    }

    impl SolveAdapter for FailingObserver {
        type Error = &'static str;

        fn value_and_gradient(
            &mut self,
            x: &[f64],
            gradient: &mut [f64],
        ) -> Result<f64, Self::Error> {
            self.evaluations += 1;
            gradient[0] = 2.0 * (x[0] - 1.0);
            Ok((x[0] - 1.0).powi(2))
        }

        fn accepted_iteration(
            &mut self,
            _iteration: Iteration<'_>,
        ) -> Result<Control, Self::Error> {
            Err("observer failure")
        }
    }

    #[test]
    fn stateful_adapter_preserves_callback_error() {
        let mut x = [3.0];
        let bounds = Bounds::new(&[-10.0], &[10.0], 1).unwrap();
        let mut adapter = FailingObserver { evaluations: 0 };
        let mut solver = Solver::default();

        let error = solver
            .minimize_with_adapter(&mut x, bounds, &mut adapter)
            .unwrap_err();

        assert!(matches!(error, SolveError::Callback("observer failure")));
        assert!(adapter.evaluations >= 2);
    }
}

#[cfg(all(test, feature = "faer-backend"))]
mod backend_tests {
    use super::*;
    use std::convert::Infallible;

    fn trajectory(kernel: Kernel) -> (Vec<(f64, f64)>, Vec<f64>, Report) {
        const N: usize = 180;
        let mut x = vec![3.0; N];
        let mut lower = vec![-100.0; N];
        let upper = vec![100.0; N];
        for i in (0..N).step_by(2) {
            lower[i] = 1.0;
        }
        let options = Options::new()
            .with_history_size(10)
            .unwrap()
            .with_projected_gradient_tolerance(1.0e-5)
            .unwrap();
        let mut solver = Solver::new(options);
        solver.kernel = kernel;
        let mut values = Vec::new();
        let report = solver
            .minimize_with_callback(
                &mut x,
                Bounds::new(&lower, &upper, N).unwrap(),
                |x, gradient| {
                    let mut value = 0.25 * (x[0] - 1.0).powi(2);
                    for i in 1..N {
                        let residual = x[i] - x[i - 1].powi(2);
                        value += residual * residual;
                    }
                    value *= 4.0;
                    let mut residual = x[1] - x[0].powi(2);
                    gradient[0] = 2.0 * (x[0] - 1.0) - 16.0 * x[0] * residual;
                    for i in 1..N - 1 {
                        let previous = residual;
                        residual = x[i + 1] - x[i].powi(2);
                        gradient[i] = 8.0 * previous - 16.0 * x[i] * residual;
                    }
                    gradient[N - 1] = 8.0 * residual;
                    Ok::<_, Infallible>(value)
                },
                |iteration| {
                    values.push((iteration.value, iteration.projected_gradient_norm));
                    Control::Continue
                },
            )
            .unwrap();
        (values, x, report)
    }

    #[test]
    fn scalar_and_faer_match_full_driver_trajectory() {
        let (scalar, scalar_x, scalar_report) = trajectory(Kernel::Scalar);
        let (faer, faer_x, faer_report) = trajectory(Kernel::Faer);
        assert_eq!(faer.len(), scalar.len());
        assert_eq!(faer_report.stats, scalar_report.stats);
        assert_eq!(faer_report.termination, scalar_report.termination);
        for ((faer_value, faer_pg), (scalar_value, scalar_pg)) in faer.iter().zip(&scalar) {
            assert!((faer_value - scalar_value).abs() <= 1.0e-12 * scalar_value.abs().max(1.0));
            assert!((faer_pg - scalar_pg).abs() <= 1.0e-11 * scalar_pg.abs().max(1.0));
        }
        for (faer, scalar) in faer_x.iter().zip(&scalar_x) {
            assert!((faer - scalar).abs() <= 1.0e-12 * scalar.abs().max(1.0));
        }
    }
}

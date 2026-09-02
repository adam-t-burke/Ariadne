//! Typed progress and termination reports.

/// Why a solve terminated.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum Termination {
    /// A configured convergence test was satisfied.
    Converged(Convergence),
    /// A configured limit or user request stopped the solve.
    Stopped(StopReason),
    /// Numerical work failed before a convergence test was satisfied.
    Failed(Failure),
}

impl Termination {
    /// Returns `true` when a convergence test ended the solve.
    pub fn converged(self) -> bool {
        matches!(self, Self::Converged(_))
    }
}

/// A satisfied convergence test.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum Convergence {
    /// The projected-gradient infinity norm reached its tolerance.
    ProjectedGradient,
    /// Relative function reduction between accepted iterates reached its tolerance.
    RelativeFunction,
}

/// A non-failure stopping condition.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum StopReason {
    /// The accepted-major-iteration limit was reached.
    MaximumIterations,
    /// The objective/gradient evaluation limit was reached.
    MaximumEvaluations,
    /// The per-line-search evaluation limit was reached.
    MaximumLineSearchEvaluations,
    /// The accepted-iteration callback requested a stop.
    User,
}

/// A numerical failure reported by the scalar engine.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum Failure {
    /// The line search could not produce an acceptable step.
    LineSearch,
    /// The computed search direction was not a descent direction.
    NoDescentDirection,
    /// Compact-history factorization or solve failed after history recovery.
    Factorization,
    /// Invalid numerical state prevented further progress.
    Numerical,
}

/// Externally observable work performed by a solve.
///
/// `evaluations` includes the initial objective callback and all line-search
/// probes. `accepted_updates` and `skipped_updates` count curvature-history
/// update attempts; a failed compact factorization is skipped and also causes
/// a history reset. Active/free counts describe the final reported point.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
#[non_exhaustive]
pub struct Stats {
    /// Number of accepted major iterations.
    pub iterations: usize,
    /// Number of objective/gradient evaluations.
    pub evaluations: usize,
    /// Number of line-search trial evaluations.
    pub line_search_probes: usize,
    /// Number of accepted curvature-history updates.
    pub accepted_updates: usize,
    /// Number of skipped curvature-history updates.
    pub skipped_updates: usize,
    /// Number of times invalid compact history was discarded.
    pub history_resets: usize,
    /// Number of variables active at a bound at the final point.
    pub active_variables: usize,
    /// Number of free variables at the final point.
    pub free_variables: usize,
}

/// Coarse phase timings for development benchmarks.
///
/// Available only with `benchmark-instrumentation`; enabling it adds timer
/// reads around hot phases and is not intended for production measurements.
#[cfg(feature = "benchmark-instrumentation")]
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct BenchmarkTimings {
    /// Total wall-clock time inside the solve.
    pub total_nanoseconds: u64,
    /// Nanoseconds spent in objective/gradient callbacks.
    pub evaluation_nanoseconds: u64,
    /// Nanoseconds spent computing search directions.
    pub direction_nanoseconds: u64,
    /// Nanoseconds spent computing generalized Cauchy points.
    pub cauchy_nanoseconds: u64,
    /// Nanoseconds spent updating the free-variable partition.
    pub freev_nanoseconds: u64,
    /// Nanoseconds spent forming compact free-set matrices.
    pub formk_nanoseconds: u64,
    /// Nanoseconds spent forming the reduced gradient.
    pub cmprlb_nanoseconds: u64,
    /// Nanoseconds spent solving the subspace minimization problem.
    pub subsm_nanoseconds: u64,
    /// Nanoseconds spent in accepted-iteration callbacks.
    pub observer_nanoseconds: u64,
    /// Nanoseconds spent updating curvature history.
    pub history_update_nanoseconds: u64,
    /// Time not attributed to the instrumented callback and major kernels.
    ///
    /// This includes validation, projection, line-search bookkeeping, vector
    /// copies, report construction, and timer overhead.
    pub residual_nanoseconds: u64,
}

/// Summary of the final in-place solution.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub struct Report {
    /// Objective value at the final reported point.
    pub value: f64,
    /// Projected-gradient infinity norm at the final reported point.
    pub projected_gradient_norm: f64,
    /// Condition that ended the solve.
    pub termination: Termination,
    /// Work counters and final active-set sizes.
    pub stats: Stats,
}

/// Information supplied after an accepted major iteration.
#[derive(Debug, Clone, Copy)]
pub struct Iteration<'a> {
    /// Accepted variable values.
    pub x: &'a [f64],
    /// Objective value at the accepted point.
    pub value: f64,
    /// Projected-gradient infinity norm at the accepted point.
    pub projected_gradient_norm: f64,
    /// Cumulative work counters and current active-set sizes.
    pub stats: Stats,
}

/// Response from an accepted-iteration callback.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum Control {
    /// Continue optimization.
    Continue,
    /// Stop and return the current accepted point.
    Stop,
}

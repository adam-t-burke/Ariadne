//! L-BFGS and L-BFGS-B optimization through Basin.
//!
//! The evaluator shares the forward and adjoint solve across cost/gradient
//! requests. Accepted steps drive progress independently of line-search trials.

use crate::ffi::ProgressCallback;
use crate::gradients::value_and_gradient_cancellable;
use crate::linear_solver::ToleranceSchedule;
use crate::types::{
    FdmCache, OptimizationState, Problem, QParameterizationMode, SolverResult, TheseusError,
    VariableSupportKind,
};
use crate::variable_supports;
use basin::{
    BoxConstraints, CostFunction, Executor, Gradient, GradientState, LbfgsState, Lbfgsb,
    MoreThuente, Solver, State, StepOutcome, TerminationReason,
};
use std::cell::{Cell, RefCell};
use std::sync::atomic::{AtomicBool, Ordering};

const DIRECT_BOX_SCALE_EPS: f64 = 1e-12;
type OptimizerState = LbfgsState<Vec<f64>>;

struct Evaluation {
    cache: FdmCache,
    parameters: Vec<f64>,
    physical_gradient: Vec<f64>,
    gradient: Vec<f64>,
    x: Vec<f64>,
    value: f64,
    valid: bool,
    loss_trace: Vec<f64>,
    /// Linear-solver iterations spent by each evaluation (parallel to
    /// `loss_trace`; zeros on `Direct`).
    linear_solver_iterations: Vec<u32>,
}

struct FdmProblem<'a> {
    problem: &'a Problem,
    cancel_flag: &'a AtomicBool,
    anchor_scales: Vec<f64>,
    lower: Vec<f64>,
    upper: Vec<f64>,
    evaluation_lower: Vec<f64>,
    evaluation_upper: Vec<f64>,
    lower_indices: Vec<usize>,
    upper_indices: Vec<usize>,
    evaluation: RefCell<Evaluation>,
    // A callback may need an accepted point after a rejected trial evaluation.
    observer_cache: RefCell<Option<FdmCache>>,
    /// Relative-residual tolerance of the iterative linear solves, driven
    /// by the (projected) gradient norm of the last accepted iterate
    /// (`run_solver`) and read into `FdmCache::solve_tolerance` before each
    /// evaluation. Ignored by `Direct`.
    tolerance_schedule: Cell<ToleranceSchedule>,
}

impl<'a> FdmProblem<'a> {
    fn new(problem: &'a Problem, cancel_flag: &'a AtomicBool) -> Result<Self, TheseusError> {
        let n_lat = variable_supports::latent_dim(problem);
        let n = problem.topology.num_edges + n_lat;
        let bounded =
            problem.solver.q_parameterization_mode == QParameterizationMode::DirectBoxBounds;
        let anchor_scales = if bounded {
            anchor_optimizer_scales(problem)
        } else {
            vec![1.0; n_lat]
        };
        if anchor_scales.len() != n_lat {
            return Err(TheseusError::Shape(
                "anchor optimizer scale length mismatch".into(),
            ));
        }
        let (lower, upper) = if bounded {
            validate_direct_box_bounds(problem)?;
            direct_box_optimizer_bounds(problem, n_lat)
        } else {
            (vec![f64::NEG_INFINITY; n], vec![f64::INFINITY; n])
        };
        let (evaluation_lower, evaluation_upper) = parameter_bounds(problem);
        Ok(Self {
            problem,
            cancel_flag,
            anchor_scales,
            lower,
            upper,
            lower_indices: finite_indices(&evaluation_lower),
            upper_indices: finite_indices(&evaluation_upper),
            evaluation_lower,
            evaluation_upper,
            evaluation: RefCell::new(Evaluation {
                cache: FdmCache::new(problem)?,
                parameters: vec![0.0; n],
                physical_gradient: vec![0.0; n],
                gradient: vec![0.0; n],
                x: vec![0.0; n],
                value: 0.0,
                valid: false,
                loss_trace: Vec::new(),
                linear_solver_iterations: Vec::new(),
            }),
            observer_cache: RefCell::new(None),
            tolerance_schedule: Cell::new(ToleranceSchedule::new(
                problem.solver.iterative.tolerance,
            )),
        })
    }

    /// Feed the (projected) gradient norm of a newly accepted iterate to the
    /// tolerance schedule; the evaluations that follow use the resulting
    /// tolerance (`Direct` ignores it).
    fn observe_accepted_gradient(&self, gradient_norm: f64) {
        let mut schedule = self.tolerance_schedule.get();
        schedule.observe(gradient_norm);
        self.tolerance_schedule.set(schedule);
    }

    fn check_cancelled(&self) -> Result<(), TheseusError> {
        if self.cancel_flag.load(Ordering::Acquire) {
            Err(TheseusError::Cancelled)
        } else {
            Ok(())
        }
    }

    fn ensure_evaluated(&self, x: &[f64]) -> Result<(), TheseusError> {
        self.check_cancelled()?;
        let mut evaluation = self.evaluation.borrow_mut();
        if evaluation.valid && evaluation.x == x {
            return Ok(());
        }
        if x.len() != evaluation.x.len() || x.iter().any(|v| !v.is_finite()) {
            return Err(TheseusError::Solver(
                "invalid or non-finite optimizer parameters".into(),
            ));
        }
        // Invalidate before mutating the cache so a failed evaluation cannot be reused.
        evaluation.valid = false;
        let Evaluation {
            cache,
            parameters,
            physical_gradient,
            gradient,
            ..
        } = &mut *evaluation;
        fill_physical_parameters(self.problem, x, &self.anchor_scales, parameters);
        cache.solve_tolerance = self.tolerance_schedule.get().current();
        let iterations_before = cache.linear_solver_totals.iterations_total;
        let value = value_and_gradient_cancellable(
            cache,
            self.problem,
            parameters,
            physical_gradient,
            &self.evaluation_lower,
            &self.evaluation_upper,
            &self.lower_indices,
            &self.upper_indices,
            Some(self.cancel_flag),
        )?;
        let iterations = cache.linear_solver_totals.iterations_total - iterations_before;
        fill_scaled_gradient(
            self.problem,
            physical_gradient,
            &self.anchor_scales,
            gradient,
        );
        if !value.is_finite() || gradient.iter().any(|g| !g.is_finite()) {
            return Err(TheseusError::Solver(
                "value_and_gradient produced NaN or Inf".into(),
            ));
        }
        evaluation.x.copy_from_slice(x);
        evaluation.value = value;
        evaluation.valid = true;
        evaluation.loss_trace.push(value);
        evaluation
            .linear_solver_iterations
            .push(u32::try_from(iterations).unwrap_or(u32::MAX));
        self.check_cancelled()
    }

    fn xyz_for(&self, x: &[f64]) -> Result<Vec<f64>, TheseusError> {
        let evaluation = self.evaluation.borrow();
        if evaluation.valid && evaluation.x == x {
            return Ok(flatten_xyz(
                &evaluation.cache,
                self.problem.topology.num_nodes,
            ));
        }
        let mut theta = vec![0.0; x.len()];
        fill_physical_parameters(self.problem, x, &self.anchor_scales, &mut theta);
        let ne = self.problem.topology.num_edges;
        let anchors = variable_supports::map_latents_to_positions(self.problem, &theta[ne..])?;
        let mut cache = self.observer_cache.borrow_mut();
        if cache.is_none() {
            *cache = Some(FdmCache::new(self.problem)?);
        }
        let cache = cache.as_mut().unwrap();
        crate::fdm::solve_fdm_with_loads(cache, &theta[..ne], self.problem, &anchors, 1e-12)?;
        Ok(flatten_xyz(cache, self.problem.topology.num_nodes))
    }

    fn report(
        &self,
        iteration: usize,
        x: &[f64],
        value: f64,
        callback: Option<ProgressCallback>,
        frequency: usize,
    ) -> Result<(), TheseusError> {
        self.check_cancelled()?;
        if let Some(callback) =
            callback.filter(|_| should_report_major_iteration(iteration, frequency))
        {
            let xyz = self.xyz_for(x)?;
            let ne = self.problem.topology.num_edges;
            let keep_going = unsafe {
                callback(
                    iteration,
                    value,
                    xyz.as_ptr(),
                    self.problem.topology.num_nodes,
                    x.as_ptr(),
                    ne,
                )
            };
            if keep_going == 0 {
                return Err(TheseusError::Cancelled);
            }
        }
        self.check_cancelled()
    }

    fn finish(
        self,
        x: &[f64],
        iterations: usize,
        reason: TerminationReason,
        gradient_norm: f64,
        state: &mut OptimizationState,
    ) -> Result<SolverResult, TheseusError> {
        self.check_cancelled()?;
        let mut evaluation = self.evaluation.into_inner();
        fill_physical_parameters(
            self.problem,
            x,
            &self.anchor_scales,
            &mut evaluation.parameters,
        );
        let (q, latents) = unpack_parameters(self.problem, &evaluation.parameters);
        let anchors = variable_supports::map_latents_to_positions(self.problem, &latents)?;
        if !evaluation.valid || evaluation.x != x {
            crate::fdm::solve_fdm_with_loads(
                &mut evaluation.cache,
                &q,
                self.problem,
                &anchors,
                1e-12,
            )?;
        }
        if self.cancel_flag.load(Ordering::Acquire) {
            return Err(TheseusError::Cancelled);
        }
        let (converged, text) = termination_text(reason);
        let metric = if self.problem.solver.q_parameterization_mode
            == QParameterizationMode::DirectBoxBounds
        {
            "projected_gradient"
        } else {
            "gradient"
        };
        let mut termination_reason = format!(
            "{text}; iterations={iterations}; evaluations={}; {metric}={gradient_norm:.3e}",
            evaluation.loss_trace.len()
        );
        let linear_solver_totals = evaluation.cache.linear_solver_totals;
        // The direct path's termination string is part of its byte-identical
        // contract; only the iterative kinds append their solver summary.
        if evaluation.cache.linear_solver_kind.is_iterative() {
            termination_reason.push_str(&format!(
                "; linear solver: {}, {} solves, {} iterations",
                linear_solver_totals.backend,
                linear_solver_totals.solves,
                linear_solver_totals.iterations_total
            ));
        }
        state.force_densities = q.clone();
        state.variable_anchor_positions = anchors.clone();
        state.variable_anchor_latents = latents;
        state.iterations = iterations;
        state.loss_trace = evaluation.loss_trace.clone();
        Ok(SolverResult {
            q,
            anchor_positions: anchors,
            xyz: evaluation.cache.nf,
            member_lengths: evaluation.cache.member_lengths,
            member_forces: evaluation.cache.member_forces,
            reactions: evaluation.cache.reactions,
            cross_section_areas: evaluation.cache.cross_section_areas,
            loss_trace: evaluation.loss_trace,
            iterations,
            converged,
            termination_reason,
            linear_solver_iterations: evaluation.linear_solver_iterations,
            linear_solver_totals,
        })
    }
}

// Borrowing the evaluator lets callbacks inspect its cache without raw pointers.
impl CostFunction for &FdmProblem<'_> {
    type Param = Vec<f64>;
    type Output = f64;
    type Error = TheseusError;
    fn cost(&self, x: &Vec<f64>) -> Result<f64, TheseusError> {
        self.ensure_evaluated(x)?;
        Ok(self.evaluation.borrow().value)
    }
}
impl Gradient for &FdmProblem<'_> {
    type Gradient = Vec<f64>;
    fn gradient(&self, x: &Vec<f64>) -> Result<Vec<f64>, TheseusError> {
        self.ensure_evaluated(x)?;
        Ok(self.evaluation.borrow().gradient.clone())
    }
    fn cost_and_gradient(&self, x: &Vec<f64>) -> Result<(f64, Vec<f64>), TheseusError> {
        self.ensure_evaluated(x)?;
        let evaluation = self.evaluation.borrow();
        Ok((evaluation.value, evaluation.gradient.clone()))
    }
}
impl BoxConstraints for &FdmProblem<'_> {
    fn lower(&self) -> &Vec<f64> {
        &self.lower
    }
    fn upper(&self) -> &Vec<f64> {
        &self.upper
    }
}

#[derive(Clone, Copy)]
struct Tolerances {
    bounded: bool,
    gradient: f64,
    cost: f64,
}

impl Tolerances {
    fn new(problem: &Problem) -> Result<Self, TheseusError> {
        let bounded =
            problem.solver.q_parameterization_mode == QParameterizationMode::DirectBoxBounds;
        let gradient = problem.solver.absolute_tolerance;
        let cost = problem.solver.relative_tolerance;
        if !gradient.is_finite()
            || !cost.is_finite()
            || (!bounded && (gradient < 0.0 || cost < 0.0))
        {
            return Err(TheseusError::Solver(
                "optimizer tolerances must be finite and nonnegative".into(),
            ));
        }
        Ok(Self {
            bounded,
            gradient: gradient.max(0.0),
            cost: cost.max(0.0),
        })
    }

    fn check(&self, gradient: f64, value: f64, previous: Option<f64>) -> Option<TerminationReason> {
        if self.gradient > 0.0
            && if self.bounded {
                gradient <= self.gradient
            } else {
                gradient < self.gradient
            }
        {
            return Some(if self.bounded {
                TerminationReason::ProjectedGradientTolerance
            } else {
                TerminationReason::GradientTolerance
            });
        }
        let previous = previous?;
        // Theseus's box tolerance has a unit floor; Basin's relative check does not.
        if self.cost > 0.0
            && if self.bounded {
                previous - value <= self.cost * previous.abs().max(value.abs()).max(1.0)
            } else {
                (previous - value).abs() < self.cost
            }
        {
            Some(if self.bounded {
                TerminationReason::RelativeCostTolerance
            } else {
                TerminationReason::CostTolerance
            })
        } else {
            None
        }
    }
}

struct CompatibleSolver<S> {
    inner: S,
    tolerances: Tolerances,
    previous_cost: Option<f64>,
}

impl<P, S> Solver<P, OptimizerState> for CompatibleSolver<S>
where
    P: CostFunction<Param = Vec<f64>, Output = f64> + BoxConstraints,
    S: Solver<P, OptimizerState>,
{
    type Error = S::Error;
    fn init(
        &mut self,
        problem: &mut basin::Problem<P>,
        state: OptimizerState,
    ) -> Result<OptimizerState, Self::Error> {
        self.inner.init(problem, state)
    }
    fn next_iter(
        &mut self,
        problem: &mut basin::Problem<P>,
        state: OptimizerState,
    ) -> Result<(OptimizerState, Option<TerminationReason>), Self::Error> {
        self.previous_cost = Some(state.cost());
        self.inner.next_iter(problem, state)
    }
    fn reset_convergence(&mut self) {
        self.previous_cost = None;
        self.inner.reset_convergence();
    }
    fn check_convergence(
        &mut self,
        problem: &basin::Problem<P>,
        state: &OptimizerState,
    ) -> Option<TerminationReason> {
        let norm = gradient_norm(
            state,
            self.tolerances.bounded,
            problem.inner().lower(),
            problem.inner().upper(),
        );
        self.tolerances
            .check(norm, state.cost(), self.previous_cost)
    }
}

fn gradient_norm(state: &OptimizerState, bounded: bool, lower: &[f64], upper: &[f64]) -> f64 {
    try_gradient_norm(state, bounded, lower, upper)
        .expect("Basin initializes the gradient before publishing a state")
}

/// The convergence metric of `state`: the Euclidean gradient norm (soft
/// bounds) or the ∞-norm of the projected gradient (box bounds); `None`
/// while the state carries no gradient.
fn try_gradient_norm(
    state: &OptimizerState,
    bounded: bool,
    lower: &[f64],
    upper: &[f64],
) -> Option<f64> {
    let gradient = state.gradient()?;
    if !bounded {
        return Some(gradient.iter().map(|g| g * g).sum::<f64>().sqrt());
    }
    Some(
        state
            .param()
            .iter()
            .zip(gradient)
            .enumerate()
            .map(|(i, (&x, &g))| {
                if g < 0.0 {
                    g.max(x - upper[i]).abs()
                } else {
                    g.min(x - lower[i]).abs()
                }
            })
            .fold(0.0, f64::max),
    )
}

fn termination_text(reason: TerminationReason) -> (bool, &'static str) {
    match reason {
        TerminationReason::ProjectedGradientTolerance => {
            (true, "converged: projected gradient tolerance reached")
        }
        TerminationReason::GradientTolerance => (true, "converged: gradient tolerance reached"),
        TerminationReason::RelativeCostTolerance => (
            true,
            "converged: relative objective reduction tolerance reached",
        ),
        TerminationReason::CostTolerance => (
            true,
            "converged: absolute objective reduction tolerance reached",
        ),
        TerminationReason::MaxIter => (false, "stopped: maximum iterations reached"),
        TerminationReason::SolverFailed => (false, "failed: solver could not find a step"),
        _ => (false, "stopped: solver terminated without convergence"),
    }
}

fn should_report_major_iteration(iteration: usize, frequency: usize) -> bool {
    iteration == 1 || iteration.is_multiple_of(frequency.max(1))
}

fn flatten_xyz(cache: &FdmCache, nn: usize) -> Vec<f64> {
    (0..nn)
        .flat_map(|i| (0..3).map(move |d| cache.nf[[i, d]]))
        .collect()
}

fn run_solver<S>(
    fdm: &FdmProblem<'_>,
    solver: S,
    initial: Vec<f64>,
    tolerances: Tolerances,
    callback: Option<ProgressCallback>,
    frequency: usize,
) -> Result<(OptimizerState, TerminationReason), TheseusError>
where
    for<'a, 'b> S: Solver<&'a FdmProblem<'b>, OptimizerState, Error = TheseusError>,
{
    let max_iterations = if tolerances.bounded {
        fdm.problem.solver.max_iterations.max(1)
    } else {
        fdm.problem.solver.max_iterations
    };
    let solver = CompatibleSolver {
        inner: solver,
        tolerances,
        previous_cost: None,
    };
    let mut stepper = Executor::new(fdm, solver, LbfgsState::new(initial, 10))
        .max_iter(max_iterations as u64)
        .into_stepper()?;
    // The initial point is the first accepted iterate: its (projected)
    // gradient norm is the reference ‖g⁺_0‖ of the adaptive tolerance policy.
    let observe = |state: &OptimizerState| {
        if let Some(norm) = try_gradient_norm(state, tolerances.bounded, &fdm.lower, &fdm.upper) {
            fdm.observe_accepted_gradient(norm);
        }
    };
    observe(stepper.state());
    let mut previous_cost = None;
    let reason = loop {
        fdm.check_cancelled()?;
        let cost = stepper.state().cost();
        match stepper.step()? {
            StepOutcome::Continue => {
                previous_cost = Some(cost);
                let state = stepper.state();
                observe(state);
                fdm.report(
                    state.iter() as usize,
                    state.param(),
                    state.cost(),
                    callback,
                    frequency,
                )?;
            }
            StepOutcome::Stopped(mut reason) => {
                // The old bounded solver checks convergence before the iteration budget.
                if reason == TerminationReason::MaxIter && tolerances.bounded {
                    let state = stepper.state();
                    let norm = gradient_norm(state, true, &fdm.lower, &fdm.upper);
                    reason = tolerances
                        .check(norm, state.cost(), previous_cost)
                        .unwrap_or(reason);
                }
                break reason;
            }
        }
    };
    fdm.check_cancelled()?;
    Ok((stepper.into_state(), reason))
}

/// Optimize force densities and variable supports, reporting accepted iterations.
///
/// A zero callback return or a set cancellation flag aborts without publishing
/// a partial result to `state`.
pub fn optimize(
    problem: &Problem,
    state: &mut OptimizationState,
    progress_cb: Option<ProgressCallback>,
    report_freq: usize,
    cancel_flag: &AtomicBool,
) -> Result<SolverResult, TheseusError> {
    crate::objectives::validate_objectives(&problem.objectives)?;
    let tolerances = Tolerances::new(problem)?;
    let fdm = FdmProblem::new(problem, cancel_flag)?;
    fdm.check_cancelled()?;
    let (result, reason) = if tolerances.bounded {
        let initial = pack_direct_box_scaled(problem, state, &fdm.anchor_scales);
        let solver = Lbfgsb::new().with_absolute_projected_gradient_tolerance(None);
        run_solver(&fdm, solver, initial, tolerances, progress_cb, report_freq)?
    } else {
        let line_search = MoreThuente {
            ftol: 1e-4,
            gtol: 0.9,
            xtol: 1e-10,
            stpmin: f64::EPSILON.sqrt(),
            stpmax: f64::INFINITY,
            maxfev: u32::MAX,
            ..MoreThuente::new()
        };
        let solver = Lbfgsb::with_line_search(line_search).unbounded();
        run_solver(
            &fdm,
            solver,
            pack_parameters(problem, state),
            tolerances,
            progress_cb,
            report_freq,
        )?
    };
    let norm = gradient_norm(&result, tolerances.bounded, &fdm.lower, &fdm.upper);
    let x = if tolerances.bounded {
        result.param()
    } else {
        result.best_param()
    };
    fdm.finish(x, result.iter() as usize, reason, norm, state)
}

/// Pack q and latent support parameters into a single θ vector.
pub fn pack_parameters(problem: &Problem, state: &OptimizationState) -> Vec<f64> {
    let ne = problem.topology.num_edges;
    let n_lat = variable_supports::latent_dim(problem);
    let mut theta = Vec::with_capacity(ne + n_lat);
    for i in 0..ne {
        theta.push(state.force_densities[i]);
    }
    if n_lat > 0 {
        if state.variable_anchor_latents.len() == n_lat {
            theta.extend_from_slice(&state.variable_anchor_latents);
        } else {
            // Backward compatibility fallback for older state payloads.
            for i in 0..state.variable_anchor_positions.nrows() {
                theta.push(state.variable_anchor_positions[[i, 0]]);
                theta.push(state.variable_anchor_positions[[i, 1]]);
                theta.push(state.variable_anchor_positions[[i, 2]]);
            }
        }
    }
    theta
}

/// Unpack θ into q and latent support parameters.
pub fn unpack_parameters(problem: &Problem, theta: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let ne = problem.topology.num_edges;
    let q = theta[..ne].to_vec();
    let lat = theta[ne..].to_vec();
    (q, lat)
}

// ─────────────────────────────────────────────────────────────
//  Bound index precomputation
// ─────────────────────────────────────────────────────────────

fn parameter_bounds(problem: &Problem) -> (Vec<f64>, Vec<f64>) {
    let n_lat = variable_supports::latent_dim(problem);
    let ne = problem.topology.num_edges;
    let mut lb = vec![f64::NEG_INFINITY; ne];
    let mut ub = vec![f64::INFINITY; ne];
    if problem.solver.q_parameterization_mode == QParameterizationMode::DirectSoftBounds {
        lb.copy_from_slice(&problem.bounds.lower);
        ub.copy_from_slice(&problem.bounds.upper);
    }
    if n_lat > 0 {
        // Support feasibility is enforced via latent maps. Keep latent vars unbounded.
        lb.extend(vec![f64::NEG_INFINITY; n_lat]);
        ub.extend(vec![f64::INFINITY; n_lat]);
    }
    (lb, ub)
}

fn finite_indices(v: &[f64]) -> Vec<usize> {
    v.iter()
        .enumerate()
        .filter(|(_, &x)| x.is_finite())
        .map(|(i, _)| i)
        .collect()
}

fn anchor_lambda(problem: &Problem) -> f64 {
    problem
        .solver
        .anchor_saturation_lambda
        .abs()
        .max(DIRECT_BOX_SCALE_EPS)
}

fn anchor_optimizer_scales(problem: &Problem) -> Vec<f64> {
    let lambda = anchor_lambda(problem);
    if problem.anchors.variable_supports.is_empty() {
        return vec![lambda; variable_supports::latent_dim(problem)];
    }

    let mut scales = Vec::with_capacity(variable_supports::latent_dim(problem));
    for support in &problem.anchors.variable_supports {
        let support_lambda =
            if support.saturation_lambda.is_finite() && support.saturation_lambda > 0.0 {
                support.saturation_lambda
            } else {
                lambda
            };
        match &support.kind {
            VariableSupportKind::Sphere { radius } => {
                let scale = support_lambda * radius.abs().max(DIRECT_BOX_SCALE_EPS);
                scales.extend_from_slice(&[scale, scale, scale]);
            }
            VariableSupportKind::Roller { enabled, .. } => {
                for &is_enabled in enabled {
                    if is_enabled {
                        scales.push(1.0);
                    }
                }
            }
            VariableSupportKind::Rail { .. } | VariableSupportKind::NurbsCurve { .. } => {
                scales.push(1.0);
            }
            VariableSupportKind::NurbsSurface { .. } => {
                scales.extend_from_slice(&[1.0, 1.0]);
            }
        }
    }
    scales
}

fn validate_direct_box_bounds(problem: &Problem) -> Result<(), TheseusError> {
    let ne = problem.topology.num_edges;
    for i in 0..ne {
        let lb = problem.bounds.lower[i];
        let ub = problem.bounds.upper[i];
        if !lb.is_finite() || !ub.is_finite() {
            return Err(TheseusError::Shape(format!(
                "DirectBoxBounds model contract requires finite two-sided q bounds at edge {i}; \
                 use DirectSoftBounds for one-sided or unbounded q (got [{lb}, {ub}])"
            )));
        }
        let span = ub - lb;
        if !span.is_finite() || span <= DIRECT_BOX_SCALE_EPS {
            return Err(TheseusError::Shape(format!(
                "DirectBoxBounds model contract requires a non-fixed finite q interval wider than \
                 {DIRECT_BOX_SCALE_EPS:e} at edge {i}; use DirectSoftBounds for fixed q \
                 (got [{lb}, {ub}])"
            )));
        }
    }
    Ok(())
}

fn pack_direct_box_scaled(
    problem: &Problem,
    state: &OptimizationState,
    anchor_scales: &[f64],
) -> Vec<f64> {
    let ne = problem.topology.num_edges;
    let n_lat = variable_supports::latent_dim(problem);
    let mut x = Vec::with_capacity(ne + n_lat);
    for i in 0..ne {
        let lower = problem.bounds.lower[i];
        let upper = problem.bounds.upper[i];
        x.push(state.force_densities[i].clamp(lower, upper));
    }

    let latents: Vec<f64> = if state.variable_anchor_latents.len() == n_lat {
        state.variable_anchor_latents.clone()
    } else {
        state
            .variable_anchor_positions
            .iter()
            .copied()
            .take(n_lat)
            .collect()
    };
    for (i, &scale) in anchor_scales.iter().enumerate().take(n_lat) {
        let scale = scale.max(DIRECT_BOX_SCALE_EPS);
        x.push(latents.get(i).copied().unwrap_or(0.0) / scale);
    }
    x
}

#[cfg(test)]
fn direct_box_scaled_to_physical(problem: &Problem, x: &[f64], anchor_scales: &[f64]) -> Vec<f64> {
    let mut theta = vec![0.0; x.len()];
    fill_physical_parameters(problem, x, anchor_scales, &mut theta);
    theta
}

fn fill_physical_parameters(
    problem: &Problem,
    x: &[f64],
    anchor_scales: &[f64],
    theta: &mut [f64],
) {
    let ne = problem.topology.num_edges;
    theta[..ne].copy_from_slice(&x[..ne]);
    for i in 0..anchor_scales.len() {
        theta[ne + i] = x[ne + i] * anchor_scales[i];
    }
}

fn fill_scaled_gradient(
    problem: &Problem,
    grad_physical: &[f64],
    anchor_scales: &[f64],
    gradient: &mut [f64],
) {
    let ne = problem.topology.num_edges;
    gradient[..ne].copy_from_slice(&grad_physical[..ne]);
    for i in 0..anchor_scales.len() {
        gradient[ne + i] = grad_physical[ne + i] * anchor_scales[i];
    }
}

#[cfg(test)]
fn direct_box_scaled_gradient(
    problem: &Problem,
    grad_physical: &[f64],
    anchor_scales: &[f64],
) -> Vec<f64> {
    let mut gradient = vec![0.0; grad_physical.len()];
    fill_scaled_gradient(problem, grad_physical, anchor_scales, &mut gradient);
    gradient
}

fn direct_box_optimizer_bounds(problem: &Problem, n_lat: usize) -> (Vec<f64>, Vec<f64>) {
    let ne = problem.topology.num_edges;
    let mut lower = problem.bounds.lower[..ne].to_vec();
    let mut upper = problem.bounds.upper[..ne].to_vec();
    let (support_lower, support_upper) = variable_supports::direct_box_optimizer_bounds(problem);
    debug_assert_eq!(support_lower.len(), n_lat);
    debug_assert_eq!(support_upper.len(), n_lat);
    lower.extend(support_lower);
    upper.extend(support_upper);
    (lower, upper)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sparse::SparseColMatOwned;
    use crate::types::{AnchorInfo, Bounds, NetworkTopology, SolverOptions};
    use ndarray::Array2;
    use std::sync::Mutex;

    static REPORTED_ITERATIONS: Mutex<Vec<usize>> = Mutex::new(Vec::new());

    unsafe extern "C" fn record_progress(
        iteration: usize,
        _loss: f64,
        _xyz: *const f64,
        _num_nodes: usize,
        _q: *const f64,
        _num_edges: usize,
    ) -> u8 {
        REPORTED_ITERATIONS.lock().unwrap().push(iteration);
        1
    }

    unsafe extern "C" fn cancel_progress(
        _iteration: usize,
        _loss: f64,
        _xyz: *const f64,
        _num_nodes: usize,
        _q: *const f64,
        _num_edges: usize,
    ) -> u8 {
        0
    }

    fn tiny_problem(bounds: Bounds) -> Problem {
        let incidence =
            SparseColMatOwned::from_coo(1, 2, &[0usize, 0usize], &[0usize, 1usize], &[-1.0, 1.0])
                .unwrap();
        let free_node_indices = vec![0usize];
        let fixed_node_indices = vec![1usize];
        let free_incidence = incidence.extract_columns(&free_node_indices);
        let fixed_incidence = incidence.extract_columns(&fixed_node_indices);

        Problem {
            topology: NetworkTopology {
                incidence,
                free_incidence,
                fixed_incidence,
                num_edges: 1,
                num_nodes: 2,
                free_node_indices,
                fixed_node_indices,
            },
            free_node_loads: Array2::from_shape_vec((1, 3), vec![0.0, 0.0, -1.0]).unwrap(),
            fixed_node_positions: Array2::from_shape_vec((1, 3), vec![1.0, 0.0, 0.0]).unwrap(),
            anchors: AnchorInfo::all_fixed(
                Array2::from_shape_vec((1, 3), vec![1.0, 0.0, 0.0]).unwrap(),
            ),
            objectives: Vec::new(),
            bounds,
            solver: SolverOptions {
                q_parameterization_mode: QParameterizationMode::DirectBoxBounds,
                ..SolverOptions::default()
            },
            self_weight: None,
            pressure: None,
        }
    }

    #[test]
    fn reports_first_and_frequency_major_iterations() {
        let reported: Vec<usize> = (1..=10)
            .filter(|&iter| should_report_major_iteration(iter, 3))
            .collect();

        assert_eq!(reported, vec![1, 3, 6, 9]);
    }

    #[test]
    fn zero_frequency_reports_every_major_iteration() {
        let reported: Vec<usize> = (1..=4)
            .filter(|&iter| should_report_major_iteration(iter, 0))
            .collect();

        assert_eq!(reported, vec![1, 2, 3, 4]);
    }

    #[test]
    fn direct_box_q_is_physical() {
        let problem = tiny_problem(Bounds {
            lower: vec![2.0],
            upper: vec![10.0],
        });
        validate_direct_box_bounds(&problem).unwrap();
        let state = OptimizationState::new(vec![4.0], Array2::zeros((0, 3)));
        let packed = pack_direct_box_scaled(&problem, &state, &[]);
        let theta = direct_box_scaled_to_physical(&problem, &packed, &[]);

        assert_eq!(packed, vec![4.0]);
        assert_eq!(theta, vec![4.0]);
    }

    #[test]
    fn direct_box_pack_clamps_initial_q_to_bounds() {
        let problem = tiny_problem(Bounds {
            lower: vec![2.0],
            upper: vec![10.0],
        });
        let state = OptimizationState::new(vec![0.0], Array2::zeros((0, 3)));
        let packed = pack_direct_box_scaled(&problem, &state, &[]);

        assert_eq!(packed, vec![2.0]);
    }

    #[test]
    fn direct_box_optimizer_uses_physical_q_bounds() {
        let problem = tiny_problem(Bounds {
            lower: vec![2.0],
            upper: vec![10.0],
        });

        let (lower, upper) = direct_box_optimizer_bounds(&problem, 0);

        assert_eq!(lower, vec![2.0]);
        assert_eq!(upper, vec![10.0]);
    }

    #[test]
    fn direct_box_requires_finite_two_sided_q_bounds() {
        let problem = tiny_problem(Bounds {
            lower: vec![0.1],
            upper: vec![f64::INFINITY],
        });

        assert!(validate_direct_box_bounds(&problem).is_err());
    }

    #[test]
    fn direct_box_gradient_scales_q_and_anchor_blocks() {
        let problem = tiny_problem(Bounds {
            lower: vec![2.0],
            upper: vec![10.0],
        });
        let grad = direct_box_scaled_gradient(&problem, &[2.0, 3.0, -4.0], &[0.5, 2.0]);

        assert_eq!(grad, vec![2.0, 1.5, -8.0]);
    }

    #[test]
    fn progress_does_not_desynchronize_evaluation_cache() {
        let problem = tiny_problem(Bounds {
            lower: vec![0.1],
            upper: vec![10.0],
        });
        let cancel = AtomicBool::new(false);
        let fdm = FdmProblem::new(&problem, &cancel).unwrap();
        fdm.ensure_evaluated(&[2.0]).unwrap();
        let original = fdm.evaluation.borrow().cache.nf.clone();
        fdm.xyz_for(&[3.0]).unwrap();
        assert_eq!(fdm.evaluation.borrow().cache.nf, original);
        assert_eq!(fdm.evaluation.borrow().x, vec![2.0]);
        let _ = (&fdm).cost_and_gradient(&vec![2.0]).unwrap();
        assert_eq!(fdm.evaluation.borrow().loss_trace.len(), 1);
    }

    #[test]
    fn adapter_reports_only_configured_accepted_iterations_in_order() {
        REPORTED_ITERATIONS.lock().unwrap().clear();
        let problem = tiny_problem(Bounds {
            lower: vec![2.0],
            upper: vec![10.0],
        });
        let cancel = AtomicBool::new(false);
        let fdm = FdmProblem::new(&problem, &cancel).unwrap();
        let value = (&fdm).cost(&vec![2.0]).unwrap();
        for iteration in 1..=7 {
            fdm.report(iteration, &[2.0], value, Some(record_progress), 3)
                .unwrap();
        }
        assert_eq!(*REPORTED_ITERATIONS.lock().unwrap(), vec![1, 3, 6]);
    }

    #[test]
    fn adapter_preserves_progress_cancellation_error() {
        let problem = tiny_problem(Bounds {
            lower: vec![2.0],
            upper: vec![10.0],
        });
        let cancel = AtomicBool::new(false);
        let fdm = FdmProblem::new(&problem, &cancel).unwrap();
        let value = (&fdm).cost(&vec![2.0]).unwrap();
        assert!(matches!(
            fdm.report(1, &[2.0], value, Some(cancel_progress), 1),
            Err(TheseusError::Cancelled)
        ));
    }

    #[test]
    fn tolerance_checks_preserve_mode_specific_meaning_and_zero_behavior() {
        let boxed = Tolerances {
            bounded: true,
            gradient: 0.0,
            cost: 0.01,
        };
        assert_eq!(
            boxed.check(1.0, 0.001, Some(0.002)),
            Some(TerminationReason::RelativeCostTolerance)
        );
        assert_eq!(
            boxed.check(1.0, 100.0, Some(100.5)),
            Some(TerminationReason::RelativeCostTolerance)
        );
        let soft = Tolerances {
            bounded: false,
            ..boxed
        };
        assert_eq!(soft.check(1.0, 100.0, Some(100.5)), None);
        assert_eq!(
            soft.check(1.0, 100.0, Some(100.001)),
            Some(TerminationReason::CostTolerance)
        );
        for bounded in [true, false] {
            assert_eq!(
                Tolerances {
                    bounded,
                    gradient: 0.0,
                    cost: 0.0
                }
                .check(0.0, 0.0, Some(0.0)),
                None
            );
        }
    }

    #[test]
    fn solver_failure_is_not_convergence() {
        assert!(!termination_text(TerminationReason::SolverFailed).0);
        assert!(!termination_text(TerminationReason::MaxIter).0);
    }
}

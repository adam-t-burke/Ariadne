//! L-BFGS optimisation driver via the `argmin` crate.
//!
//! Wraps the hand-coded `value_and_gradient` into argmin's `CostFunction`
//! + `Gradient` traits, then runs L-BFGS with the user's solver options.
//!
//! Uses `Vec<f64>` as the argmin parameter type to avoid ndarray version
//! conflicts between our ndarray 0.16 and argmin-math's bundled ndarray.

use crate::ffi::ProgressCallback;
use crate::gradients::value_and_gradient;
use crate::types::{
    FdmCache, OptimizationState, Problem, QParameterizationMode, SolverResult, TheseusError,
    VariableSupportKind,
};
use crate::variable_supports;
use argmin::core::observers::{Observe, ObserverMode};
use argmin::core::{
    CostFunction, Error, Executor, Gradient, IterState, State, TerminationReason,
    TerminationStatus, KV,
};
use argmin::solver::linesearch::MoreThuenteLineSearch;
use argmin::solver::quasinewton::LBFGS;
use ariadne_lbfgsb::{
    Backend, Bounds as SolverBounds, Control, Convergence, Failure, Iteration, Options,
    SolveAdapter, SolveError, Solver, StopReason, Termination,
};
use std::cell::RefCell;
use std::rc::Rc;
use std::sync::atomic::{AtomicBool, Ordering};

const DIRECT_BOX_SCALE_EPS: f64 = 1e-12;
type EvaluationRecord = Rc<RefCell<Option<(Vec<f64>, f64, Vec<f64>)>>>;

// ─────────────────────────────────────────────────────────────
//  argmin problem wrapper
// ─────────────────────────────────────────────────────────────

/// Wraps the FDM problem + cache + barrier data so argmin can evaluate
/// cost and gradient.
///
/// `RefCell` is used for the cache because argmin's `CostFunction` /
/// `Gradient` traits take `&self`, but our forward solver mutates the cache.
/// The solver is single-threaded, so the borrow never actually conflicts,
/// but `RefCell` gives us debug-mode borrow checking for free.
///
/// **Evaluation cache**: argmin calls `cost(θ)` and `gradient(θ)` separately
/// at the same θ each iteration.  We cache the last `(θ, loss, grad)` so the
/// expensive forward + adjoint solve runs only once per unique θ.
struct FdmProblem<'a> {
    problem: &'a Problem,
    cache: Rc<RefCell<FdmCache>>,
    lb: Vec<f64>,
    ub: Vec<f64>,
    lb_idx: Vec<usize>,
    ub_idx: Vec<usize>,
    /// Cached (θ, loss, gradient) from the last evaluation.
    last_eval: EvaluationRecord,
    /// Loss value recorded at each unique evaluation.
    loss_trace: RefCell<Vec<f64>>,
    cancel_flag: &'a AtomicBool,
}

impl<'a> FdmProblem<'a> {
    fn check_cancelled(&self) -> Result<(), Error> {
        if self.cancel_flag.load(Ordering::Acquire) {
            return Err(Error::msg(TheseusError::Cancelled.to_string()));
        }
        Ok(())
    }

    /// Ensure the cache contains results for `theta`.
    /// If θ matches the cached value, this is a no-op.
    /// Otherwise, runs the full forward + adjoint solve.
    fn ensure_evaluated(&self, theta: &[f64]) -> Result<(), argmin::core::Error> {
        self.check_cancelled()?;
        {
            let cached = self.last_eval.borrow();
            if let Some((ref t, _, _)) = *cached {
                if t == theta {
                    return Ok(());
                }
            }
        }
        // Reject θ vectors containing NaN/Inf before attempting the solve
        if theta.iter().any(|v| !v.is_finite()) {
            return Err(argmin::core::Error::msg("theta contains NaN or Inf"));
        }

        // Cache miss — run the full solve
        let mut fdm_cache = self.cache.borrow_mut();
        let mut grad = vec![0.0; theta.len()];
        let val = value_and_gradient(
            &mut fdm_cache,
            self.problem,
            theta,
            &mut grad,
            &self.lb,
            &self.ub,
            &self.lb_idx,
            &self.ub_idx,
        )
        .map_err(|e| argmin::core::Error::msg(e.to_string()))?;

        // Guard against NaN/Inf in loss or gradient
        if !val.is_finite() || grad.iter().any(|g| !g.is_finite()) {
            return Err(argmin::core::Error::msg(
                "value_and_gradient produced NaN or Inf",
            ));
        }

        {
            let mut trace = self.loss_trace.borrow_mut();
            trace.push(val);
        }

        *self.last_eval.borrow_mut() = Some((theta.to_vec(), val, grad));
        self.check_cancelled()
    }
}

type LbfgsState = IterState<Vec<f64>, Vec<f64>, (), (), (), f64>;

/// Emits progress only after accepted L-BFGS iterations. Argmin calls
/// observers after `next_iter`, so this avoids More-Thuente trial states.
struct MajorIterationProgressObserver {
    problem: *const Problem,
    cancel_flag: *const AtomicBool,
    evaluation_cache: Rc<RefCell<FdmCache>>,
    observer_cache: RefCell<FdmCache>,
    last_eval: EvaluationRecord,
    callback: ProgressCallback,
    report_frequency: usize,
}

impl MajorIterationProgressObserver {
    fn problem(&self) -> Result<&Problem, Error> {
        unsafe {
            self.problem
                .as_ref()
                .ok_or_else(|| Error::msg("progress observer problem pointer is null"))
        }
    }

    fn check_cancelled(&self) -> Result<(), Error> {
        let cancelled = unsafe {
            self.cancel_flag
                .as_ref()
                .ok_or_else(|| Error::msg("progress observer cancellation pointer is null"))?
                .load(Ordering::Acquire)
        };
        if cancelled {
            Err(Error::msg(TheseusError::Cancelled.to_string()))
        } else {
            Ok(())
        }
    }

    fn xyz_for(&self, theta: &[f64]) -> Result<Vec<f64>, Error> {
        let problem = self.problem()?;
        let nn = problem.topology.num_nodes;
        let ne = problem.topology.num_edges;

        if let Some((cached_theta, _, _)) = self.last_eval.borrow().as_ref() {
            if cached_theta.as_slice() == theta {
                let cache = self.evaluation_cache.borrow();
                return Ok(flatten_xyz(&cache, nn));
            }
        }

        let q = theta[..ne].to_vec();
        let anchors = variable_supports::map_latents_to_positions(problem, &theta[ne..])
            .map_err(|e| Error::msg(e.to_string()))?;
        let mut cache = self.observer_cache.borrow_mut();
        crate::fdm::solve_fdm_with_loads(&mut cache, &q, problem, &anchors, 1e-12)
            .map_err(|e| Error::msg(e.to_string()))?;
        Ok(flatten_xyz(&cache, nn))
    }
}

impl Observe<LbfgsState> for MajorIterationProgressObserver {
    fn observe_iter(&mut self, state: &LbfgsState, _kv: &KV) -> Result<(), Error> {
        let major_iteration = state.get_iter() as usize + 1;
        if !should_report_major_iteration(major_iteration, self.report_frequency) {
            return Ok(());
        }

        let problem = self.problem()?;
        let ne = problem.topology.num_edges;
        let Some(theta) = state.get_param() else {
            // L-BFGS can return a terminated state without restoring `param`
            // when the inner line search exits early. That is not an accepted
            // major iteration state, so do not stream it.
            return Ok(());
        };
        if theta.len() < ne {
            return Err(Error::msg(
                "progress observer parameter vector is too short",
            ));
        }

        let xyz_flat = self.xyz_for(theta)?;
        let q = theta[..ne].to_vec();
        let should_continue = unsafe {
            (self.callback)(
                major_iteration,
                state.get_cost(),
                xyz_flat.as_ptr(),
                problem.topology.num_nodes,
                q.as_ptr(),
                ne,
            )
        };
        if should_continue == 0 {
            return Err(Error::msg(TheseusError::Cancelled.to_string()));
        }
        self.check_cancelled()
    }
}

fn should_report_major_iteration(major_iteration: usize, report_frequency: usize) -> bool {
    let frequency = report_frequency.max(1);
    major_iteration == 1 || major_iteration % frequency == 0
}

fn flatten_xyz(cache: &FdmCache, nn: usize) -> Vec<f64> {
    let nf = &cache.nf;
    (0..nn)
        .flat_map(|i| (0..3).map(move |d| nf[[i, d]]))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sparse::SparseColMatOwned;
    use crate::types::{AnchorInfo, Bounds, NetworkTopology, SolverOptions};
    use ariadne_lbfgsb::Stats;
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
    fn legacy_observer_does_not_desynchronize_evaluation_cache() {
        let problem = tiny_problem(Bounds {
            lower: vec![0.1],
            upper: vec![10.0],
        });
        let evaluation_cache = Rc::new(RefCell::new(FdmCache::new(&problem).unwrap()));
        let evaluation_snapshot = evaluation_cache.borrow().nf.clone();
        let last_eval = Rc::new(RefCell::new(Some((vec![2.0], 7.0, vec![3.0]))));
        let cancel = AtomicBool::new(false);
        let observer = MajorIterationProgressObserver {
            problem: &problem,
            cancel_flag: &cancel,
            evaluation_cache: evaluation_cache.clone(),
            observer_cache: RefCell::new(FdmCache::new(&problem).unwrap()),
            last_eval: last_eval.clone(),
            callback: record_progress,
            report_frequency: 1,
        };

        observer.xyz_for(&[3.0]).unwrap();

        assert_eq!(evaluation_cache.borrow().nf, evaluation_snapshot);
        assert_eq!(last_eval.borrow().as_ref().unwrap().0, vec![2.0]);
    }

    #[test]
    fn direct_box_wires_iteration_limits_above_one_thousand_exactly() {
        let mut problem = tiny_problem(Bounds {
            lower: vec![2.0],
            upper: vec![10.0],
        });
        problem.solver.max_iterations = 1_237;

        let options = direct_box_solver_options(&problem).unwrap();

        assert_eq!(options.max_iterations(), 1_237);
    }

    #[test]
    fn adapter_reports_only_configured_accepted_iterations_in_order() {
        REPORTED_ITERATIONS.lock().unwrap().clear();
        let problem = tiny_problem(Bounds {
            lower: vec![2.0],
            upper: vec![10.0],
        });
        let cancel = AtomicBool::new(false);
        let mut adapter =
            FdmSolverAdapter::new(&problem, &cancel, Some(record_progress), 3, Vec::new()).unwrap();
        let x = [2.0];
        let mut gradient = [0.0];
        let loss = adapter.value_and_gradient(&x, &mut gradient).unwrap();

        for iteration in 1..=7 {
            let mut stats = Stats::default();
            stats.iterations = iteration;
            adapter
                .accepted_iteration(Iteration {
                    x: &x,
                    value: loss,
                    projected_gradient_norm: 0.0,
                    stats,
                })
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
        let mut adapter =
            FdmSolverAdapter::new(&problem, &cancel, Some(cancel_progress), 1, Vec::new()).unwrap();
        let x = [2.0];
        let mut gradient = [0.0];
        let loss = adapter.value_and_gradient(&x, &mut gradient).unwrap();
        let mut stats = Stats::default();
        stats.iterations = 1;

        let error = adapter
            .accepted_iteration(Iteration {
                x: &x,
                value: loss,
                projected_gradient_norm: 0.0,
                stats,
            })
            .unwrap_err();

        assert!(matches!(error, TheseusError::Cancelled));
    }
}

impl<'a> CostFunction for FdmProblem<'a> {
    type Param = Vec<f64>;
    type Output = f64;

    fn cost(&self, theta: &Self::Param) -> Result<Self::Output, argmin::core::Error> {
        self.ensure_evaluated(theta)?;
        let cached = self.last_eval.borrow();
        Ok(cached.as_ref().unwrap().1)
    }
}

impl<'a> Gradient for FdmProblem<'a> {
    type Param = Vec<f64>;
    type Gradient = Vec<f64>;

    fn gradient(&self, theta: &Self::Param) -> Result<Self::Gradient, argmin::core::Error> {
        self.ensure_evaluated(theta)?;
        let cached = self.last_eval.borrow();
        Ok(cached.as_ref().unwrap().2.clone())
    }
}

// ─────────────────────────────────────────────────────────────
//  Parameter packing / unpacking
// ─────────────────────────────────────────────────────────────

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

fn lbfgs_tolerances(problem: &Problem) -> (f64, f64) {
    (
        problem.solver.absolute_tolerance,
        problem.solver.relative_tolerance,
    )
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
        let support_lambda = if support.saturation_lambda.is_finite()
            && support.saturation_lambda > 0.0
        {
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

struct FdmSolverAdapter<'a> {
    problem: &'a Problem,
    cancel_flag: &'a AtomicBool,
    progress_cb: Option<ProgressCallback>,
    report_frequency: usize,
    anchor_scales: Vec<f64>,
    evaluation_lower: Vec<f64>,
    evaluation_upper: Vec<f64>,
    lower_indices: Vec<usize>,
    upper_indices: Vec<usize>,
    cache: FdmCache,
    physical_parameters: Vec<f64>,
    physical_gradient: Vec<f64>,
    evaluated_x: Vec<f64>,
    has_evaluation: bool,
    xyz_flat: Vec<f64>,
    loss_trace: Vec<f64>,
}

impl<'a> FdmSolverAdapter<'a> {
    fn new(
        problem: &'a Problem,
        cancel_flag: &'a AtomicBool,
        progress_cb: Option<ProgressCallback>,
        report_frequency: usize,
        anchor_scales: Vec<f64>,
    ) -> Result<Self, TheseusError> {
        let n = problem.topology.num_edges + variable_supports::latent_dim(problem);
        Ok(Self {
            problem,
            cancel_flag,
            progress_cb,
            report_frequency: report_frequency.max(1),
            anchor_scales,
            evaluation_lower: vec![f64::NEG_INFINITY; n],
            evaluation_upper: vec![f64::INFINITY; n],
            lower_indices: Vec::new(),
            upper_indices: Vec::new(),
            cache: FdmCache::new(problem)?,
            physical_parameters: vec![0.0; n],
            physical_gradient: vec![0.0; n],
            evaluated_x: vec![0.0; n],
            has_evaluation: false,
            xyz_flat: vec![0.0; problem.topology.num_nodes * 3],
            loss_trace: Vec::new(),
        })
    }

    fn check_cancelled(&self) -> Result<(), TheseusError> {
        if self.cancel_flag.load(Ordering::Acquire) {
            Err(TheseusError::Cancelled)
        } else {
            Ok(())
        }
    }

    fn fill_physical_parameters(&mut self, x: &[f64]) {
        fill_physical_parameters(
            self.problem,
            x,
            &self.anchor_scales,
            &mut self.physical_parameters,
        );
    }

    fn sync_final_cache(&mut self, x: &[f64]) -> Result<(), TheseusError> {
        if self.has_evaluation && self.evaluated_x == x {
            return Ok(());
        }
        self.fill_physical_parameters(x);
        let ne = self.problem.topology.num_edges;
        let anchors = variable_supports::map_latents_to_positions(
            self.problem,
            &self.physical_parameters[ne..],
        )?;
        crate::fdm::solve_fdm_with_loads(
            &mut self.cache,
            &self.physical_parameters[..ne],
            self.problem,
            &anchors,
            1e-12,
        )?;
        self.evaluated_x.copy_from_slice(x);
        self.has_evaluation = true;
        Ok(())
    }

    fn finish(
        mut self,
        x: &[f64],
        report: ariadne_lbfgsb::Report,
        state: &mut OptimizationState,
    ) -> Result<SolverResult, TheseusError> {
        self.sync_final_cache(x)?;
        self.fill_physical_parameters(x);
        let ne = self.problem.topology.num_edges;
        let q = self.physical_parameters[..ne].to_vec();
        let latents = self.physical_parameters[ne..].to_vec();
        let anchors = variable_supports::map_latents_to_positions(self.problem, &latents)?;
        let termination_reason = direct_box_termination_text(&report);
        let converged = report.termination.converged();

        state.force_densities = q.clone();
        state.variable_anchor_positions = anchors.clone();
        state.variable_anchor_latents = latents;
        state.iterations = report.stats.iterations;
        state.loss_trace = self.loss_trace.clone();

        Ok(SolverResult {
            q,
            anchor_positions: anchors,
            xyz: self.cache.nf,
            member_lengths: self.cache.member_lengths,
            member_forces: self.cache.member_forces,
            reactions: self.cache.reactions,
            loss_trace: self.loss_trace,
            iterations: report.stats.iterations,
            converged,
            termination_reason,
            cross_section_areas: self.cache.cross_section_areas,
        })
    }
}

impl SolveAdapter for FdmSolverAdapter<'_> {
    type Error = TheseusError;

    fn value_and_gradient(&mut self, x: &[f64], gradient: &mut [f64]) -> Result<f64, Self::Error> {
        self.check_cancelled()?;
        if x.iter().any(|value| !value.is_finite()) {
            return Err(TheseusError::Solver(
                "DirectBoxBounds received non-finite solver parameters".into(),
            ));
        }

        self.fill_physical_parameters(x);
        let loss = value_and_gradient(
            &mut self.cache,
            self.problem,
            &self.physical_parameters,
            &mut self.physical_gradient,
            &self.evaluation_lower,
            &self.evaluation_upper,
            &self.lower_indices,
            &self.upper_indices,
        )?;
        fill_scaled_gradient(
            self.problem,
            &self.physical_gradient,
            &self.anchor_scales,
            gradient,
        );
        self.loss_trace.push(loss);
        self.evaluated_x.copy_from_slice(x);
        self.has_evaluation = true;
        self.check_cancelled()?;
        Ok(loss)
    }

    fn accepted_iteration(&mut self, iteration: Iteration<'_>) -> Result<Control, Self::Error> {
        self.check_cancelled()?;
        let Some(progress_callback) = self.progress_cb.filter(|_| {
            should_report_major_iteration(iteration.stats.iterations, self.report_frequency)
        }) else {
            return Ok(Control::Continue);
        };

        for (node, xyz) in self.xyz_flat.chunks_exact_mut(3).enumerate() {
            xyz.copy_from_slice(&[
                self.cache.nf[[node, 0]],
                self.cache.nf[[node, 1]],
                self.cache.nf[[node, 2]],
            ]);
        }
        let ne = self.problem.topology.num_edges;
        let should_continue = unsafe {
            progress_callback(
                iteration.stats.iterations,
                iteration.value,
                self.xyz_flat.as_ptr(),
                self.problem.topology.num_nodes,
                self.physical_parameters.as_ptr(),
                ne,
            )
        };
        if should_continue == 0 {
            return Err(TheseusError::Cancelled);
        }
        self.check_cancelled()?;
        Ok(Control::Continue)
    }
}

fn map_direct_box_solve_error(error: SolveError<TheseusError>) -> TheseusError {
    match error {
        SolveError::Objective(error) | SolveError::Callback(error) => error,
        other => TheseusError::Solver(other.to_string()),
    }
}

fn direct_box_termination_text(report: &ariadne_lbfgsb::Report) -> String {
    let reason = match report.termination {
        Termination::Converged(Convergence::ProjectedGradient) => {
            "converged: projected gradient tolerance reached"
        }
        Termination::Converged(Convergence::RelativeFunction) => {
            "converged: relative objective reduction tolerance reached"
        }
        Termination::Stopped(StopReason::MaximumIterations) => {
            "stopped: maximum iterations reached"
        }
        Termination::Stopped(StopReason::MaximumEvaluations) => {
            "stopped: maximum objective evaluations reached"
        }
        Termination::Stopped(StopReason::MaximumLineSearchEvaluations) => {
            "stopped: maximum line-search evaluations reached"
        }
        Termination::Stopped(StopReason::User) => "stopped: requested by callback",
        Termination::Failed(Failure::LineSearch) => "failed: line search could not find a step",
        Termination::Failed(Failure::Numerical) => "failed: numerical error",
        _ => "stopped: unrecognized solver termination",
    };
    format!(
        "{reason}; iterations={}; evaluations={}; projected_gradient={:.3e}",
        report.stats.iterations, report.stats.evaluations, report.projected_gradient_norm
    )
}

fn direct_box_solver_options(problem: &Problem) -> Result<Options, TheseusError> {
    Options::new()
        .with_history_size(10)
        .and_then(|options| options.with_backend(Backend::Faer))
        .and_then(|options| options.with_max_iterations(problem.solver.max_iterations.max(1)))
        .and_then(|options| {
            options.with_projected_gradient_tolerance(problem.solver.absolute_tolerance.max(0.0))
        })
        .and_then(|options| {
            options.with_relative_function_tolerance(problem.solver.relative_tolerance.max(0.0))
        })
        .map_err(|error| TheseusError::Solver(error.to_string()))
}

fn optimize_direct_box_bounds(
    problem: &Problem,
    state: &mut OptimizationState,
    progress_cb: Option<ProgressCallback>,
    report_freq: usize,
    cancel_flag: &AtomicBool,
) -> Result<SolverResult, TheseusError> {
    crate::objectives::validate_objectives(&problem.objectives)?;

    let ne = problem.topology.num_edges;
    let n_lat = variable_supports::latent_dim(problem);
    validate_direct_box_bounds(problem)?;
    let anchor_scales = anchor_optimizer_scales(problem);
    if anchor_scales.len() != n_lat {
        return Err(TheseusError::Shape(format!(
            "anchor optimizer scale length mismatch: got {}, expected {n_lat}",
            anchor_scales.len()
        )));
    }

    let mut x = pack_direct_box_scaled(problem, state, &anchor_scales);
    let (lower, upper) = direct_box_optimizer_bounds(problem, n_lat);
    let options = direct_box_solver_options(problem)?;
    let bounds = SolverBounds::new(&lower, &upper, ne + n_lat)
        .map_err(|error| TheseusError::Solver(error.to_string()))?;
    let mut adapter = FdmSolverAdapter::new(
        problem,
        cancel_flag,
        progress_cb,
        report_freq,
        anchor_scales,
    )?;
    let mut solver = Solver::new(options);
    let report = solver
        .minimize_with_adapter(&mut x, bounds, &mut adapter)
        .map_err(map_direct_box_solve_error)?;
    adapter.finish(&x, report, state)
}

// ─────────────────────────────────────────────────────────────
//  Top-level optimisation entry point
// ─────────────────────────────────────────────────────────────

/// Run L-BFGS optimisation on the FDM problem.
///
/// `progress_cb` / `report_freq` control an optional FFI callback invoked
/// every `report_freq` accepted L-BFGS iterations with the current node positions.
pub fn optimize(
    problem: &Problem,
    state: &mut OptimizationState,
    progress_cb: Option<ProgressCallback>,
    report_freq: usize,
    cancel_flag: &AtomicBool,
) -> Result<SolverResult, TheseusError> {
    if problem.solver.q_parameterization_mode == QParameterizationMode::DirectBoxBounds {
        return optimize_direct_box_bounds(problem, state, progress_cb, report_freq, cancel_flag);
    }

    crate::objectives::validate_objectives(&problem.objectives)?;

    let cache = FdmCache::new(problem)?;

    let (lb, ub) = parameter_bounds(problem);
    let lb_idx = finite_indices(&lb);
    let ub_idx = finite_indices(&ub);

    let init_param = pack_parameters(problem, state);

    let cache = Rc::new(RefCell::new(cache));
    let last_eval = Rc::new(RefCell::new(None));
    let report_frequency = if report_freq == 0 { 1 } else { report_freq };

    let fdm_problem = FdmProblem {
        problem,
        cache: cache.clone(),
        lb,
        ub,
        lb_idx,
        ub_idx,
        last_eval: last_eval.clone(),
        loss_trace: RefCell::new(Vec::new()),
        cancel_flag,
    };

    // Configure L-BFGS with user-specified tolerances
    let linesearch = MoreThuenteLineSearch::new();
    let (tolerance_grad, tolerance_cost) = lbfgs_tolerances(problem);
    let solver = LBFGS::new(linesearch, 10)
        .with_tolerance_grad(tolerance_grad)
        .map_err(|e| TheseusError::Solver(format!("tolerance_grad: {e}")))?
        .with_tolerance_cost(tolerance_cost)
        .map_err(|e| TheseusError::Solver(format!("tolerance_cost: {e}")))?;

    let mut executor = Executor::new(fdm_problem, solver).configure(|config| {
        config
            .param(init_param)
            .max_iters(problem.solver.max_iterations as u64)
            .target_cost(f64::NEG_INFINITY)
    });
    if let Some(callback) = progress_cb {
        executor = executor.add_observer(
            MajorIterationProgressObserver {
                problem: problem as *const Problem,
                cancel_flag,
                evaluation_cache: cache.clone(),
                observer_cache: RefCell::new(FdmCache::new(problem)?),
                last_eval,
                callback,
                report_frequency,
            },
            ObserverMode::Always,
        );
    }

    let result = executor.run()?;
    let loss_trace = result
        .problem
        .problem
        .as_ref()
        .map(|p| p.loss_trace.borrow().clone())
        .unwrap_or_default();

    // Extract solution
    let best_param = result
        .state()
        .get_best_param()
        .ok_or_else(|| TheseusError::Solver("L-BFGS returned no best parameters".into()))?;
    let (q, latents) = unpack_parameters(problem, best_param);
    let anchors = variable_supports::map_latents_to_positions(problem, &latents)?;

    // Reuse the evaluation cache from the final θ — no redundant forward solve.
    if let Some(ref p) = result.problem.problem {
        p.ensure_evaluated(best_param)?;
    }
    let final_cache = cache.borrow();

    let termination_status = result.state().get_termination_status();
    let converged = matches!(
        termination_status,
        TerminationStatus::Terminated(TerminationReason::SolverConverged)
    );

    let termination_reason = match termination_status {
        TerminationStatus::Terminated(reason) => format!("{reason}"),
        TerminationStatus::NotTerminated => "not terminated".to_string(),
    };

    state.force_densities = q.clone();
    state.variable_anchor_positions = anchors.clone();
    state.variable_anchor_latents = latents.clone();
    state.iterations = result.state().get_iter() as usize;
    state.loss_trace = loss_trace.clone();

    Ok(SolverResult {
        q,
        anchor_positions: anchors,
        xyz: final_cache.nf.clone(),
        member_lengths: final_cache.member_lengths.clone(),
        member_forces: final_cache.member_forces.clone(),
        reactions: final_cache.reactions.clone(),
        loss_trace,
        iterations: state.iterations,
        converged,
        termination_reason,
        cross_section_areas: final_cache.cross_section_areas.clone(),
    })
}

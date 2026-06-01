//! L-BFGS optimisation driver via the `argmin` crate.
//!
//! Wraps the hand-coded `value_and_gradient` into argmin's `CostFunction`
//! + `Gradient` traits, then runs L-BFGS with the user's solver options.
//!
//! Uses `Vec<f64>` as the argmin parameter type to avoid ndarray version
//! conflicts between our ndarray 0.16 and argmin-math's bundled ndarray.

use crate::ffi::ProgressCallback;
use crate::gradients::value_and_gradient;
use crate::q_parameterization;
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
use lbfgsb_rs_pure::{IterationControl, Status, LBFGSB};
use std::cell::RefCell;
use std::rc::Rc;
use std::sync::atomic::{AtomicBool, Ordering};

const DIRECT_BOX_SCALE_EPS: f64 = 1e-12;
const DIRECT_BOX_ERROR_COST: f64 = 1.0e300;
const DIRECT_BOX_MAX_BACKTRACK: usize = 32;

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
    last_eval: Rc<RefCell<Option<(Vec<f64>, f64, Vec<f64>)>>>,
    /// Loss value recorded at each unique evaluation.
    loss_trace: RefCell<Vec<f64>>,
    cancel_flag: &'a AtomicBool,
}

impl<'a> FdmProblem<'a> {
    fn check_cancelled(&self) -> Result<(), Error> {
        if self.cancel_flag.load(Ordering::Relaxed) {
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
    cache: Rc<RefCell<FdmCache>>,
    last_eval: Rc<RefCell<Option<(Vec<f64>, f64, Vec<f64>)>>>,
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

    fn xyz_for(&self, theta: &[f64]) -> Result<Vec<f64>, Error> {
        let problem = self.problem()?;
        let nn = problem.topology.num_nodes;
        let ne = problem.topology.num_edges;

        if let Some((cached_theta, _, _)) = self.last_eval.borrow().as_ref() {
            if cached_theta.as_slice() == theta {
                let cache = self.cache.borrow();
                return Ok(flatten_xyz(&cache, nn));
            }
        }

        let q_mode = problem.solver.q_parameterization_mode;
        let mut q = vec![0.0; ne];
        q_parameterization::map_q_slice(
            q_mode,
            &theta[..ne],
            &problem.bounds.lower,
            &problem.bounds.upper,
            &mut q,
        );
        let anchors = variable_supports::map_latents_to_positions(problem, &theta[ne..])
            .map_err(|e| Error::msg(e.to_string()))?;
        let mut cache = self.cache.borrow_mut();
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
        let q_mode = problem.solver.q_parameterization_mode;
        let mut q = vec![0.0; ne];
        q_parameterization::map_q_slice(
            q_mode,
            &theta[..ne],
            &problem.bounds.lower,
            &problem.bounds.upper,
            &mut q,
        );
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

        Ok(())
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
    use ndarray::Array2;

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
    fn direct_box_unit_coordinate_maps_to_physical_q() {
        let problem = tiny_problem(Bounds {
            lower: vec![2.0],
            upper: vec![10.0],
        });
        let q_spans = direct_box_spans(&problem).unwrap();
        let theta = direct_box_scaled_to_physical(&problem, &[0.25], &q_spans, &[]);

        assert_eq!(theta, vec![4.0]);
    }

    #[test]
    fn direct_box_requires_finite_two_sided_q_bounds() {
        let problem = tiny_problem(Bounds {
            lower: vec![0.1],
            upper: vec![f64::INFINITY],
        });

        assert!(direct_box_spans(&problem).is_err());
    }

    #[test]
    fn direct_box_gradient_scales_q_and_anchor_blocks() {
        let grad = direct_box_scaled_gradient(&[2.0, 3.0, -4.0], &[10.0], &[0.5, 2.0]);

        assert_eq!(grad, vec![20.0, 1.5, -8.0]);
    }

    #[test]
    fn projected_gradient_fallback_reduces_boxed_objective() {
        let mut x = vec![0.95];
        let lower = vec![0.0];
        let upper = vec![1.0];
        let mut objective = |x: &[f64]| {
            let residual = x[0] - 0.25;
            (residual * residual, vec![2.0 * residual])
        };

        let (loss, iterations, status) =
            projected_gradient_fallback(&mut x, &lower, &upper, 20, 1e-8, &mut objective);

        assert!(matches!(status, Status::Converged | Status::MaxIter));
        assert!(iterations > 0);
        assert!(loss < 1e-8);
        assert!((x[0] - 0.25).abs() < 1e-4);
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
        theta.push(q_parameterization::inverse_single(
            problem.solver.q_parameterization_mode,
            state.force_densities[i],
            problem.bounds.lower[i],
            problem.bounds.upper[i],
        ));
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
    let mut q = vec![0.0; ne];
    q_parameterization::map_q_slice(
        problem.solver.q_parameterization_mode,
        &theta[..ne],
        &problem.bounds.lower,
        &problem.bounds.upper,
        &mut q,
    );
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
    if problem.solver.q_parameterization_mode == QParameterizationMode::ImplicitBounded {
        // In implicit mode L-BFGS sees the latent z-gradient
        // grad_z = grad_q * dq/dz. Near box limits dq/dz can be small while
        // the physical q-space objective is still improving, so default
        // tolerances can report convergence much too early.
        return (
            (problem.solver.absolute_tolerance * 1e-3).max(f64::EPSILON),
            (problem.solver.relative_tolerance * 1e-3).max(f64::EPSILON),
        );
    }

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
        match &support.kind {
            VariableSupportKind::Sphere { radius } => {
                let scale = lambda * radius.abs().max(DIRECT_BOX_SCALE_EPS);
                scales.extend_from_slice(&[scale, scale, scale]);
            }
            VariableSupportKind::Roller {
                enabled,
                lower,
                upper,
            } => {
                for axis in 0..3 {
                    if enabled[axis] {
                        scales.push(
                            lambda * (upper[axis] - lower[axis]).abs().max(DIRECT_BOX_SCALE_EPS),
                        );
                    }
                }
            }
            VariableSupportKind::Rail { start, end } => {
                let dx = end[0] - start[0];
                let dy = end[1] - start[1];
                let dz = end[2] - start[2];
                let length = (dx * dx + dy * dy + dz * dz).sqrt();
                scales.push(lambda * length.max(DIRECT_BOX_SCALE_EPS));
            }
            VariableSupportKind::NurbsCurve { domain, .. } => {
                scales.push(lambda * (domain[1] - domain[0]).abs().max(DIRECT_BOX_SCALE_EPS));
            }
            VariableSupportKind::NurbsSurface {
                domain_u, domain_v, ..
            } => {
                scales.push(lambda * (domain_u[1] - domain_u[0]).abs().max(DIRECT_BOX_SCALE_EPS));
                scales.push(lambda * (domain_v[1] - domain_v[0]).abs().max(DIRECT_BOX_SCALE_EPS));
            }
        }
    }
    scales
}

fn direct_box_spans(problem: &Problem) -> Result<Vec<f64>, TheseusError> {
    let ne = problem.topology.num_edges;
    let mut spans = Vec::with_capacity(ne);
    for i in 0..ne {
        let lb = problem.bounds.lower[i];
        let ub = problem.bounds.upper[i];
        if !lb.is_finite() || !ub.is_finite() {
            return Err(TheseusError::Shape(
                "DirectBoxBounds requires finite lower and upper q bounds for every edge".into(),
            ));
        }
        let span = ub - lb;
        if !span.is_finite() || span <= DIRECT_BOX_SCALE_EPS {
            return Err(TheseusError::Shape(format!(
                "DirectBoxBounds requires increasing finite q bounds at edge {i}: [{lb}, {ub}]"
            )));
        }
        spans.push(span);
    }
    Ok(spans)
}

fn pack_direct_box_scaled(
    problem: &Problem,
    state: &OptimizationState,
    q_spans: &[f64],
    anchor_scales: &[f64],
) -> Vec<f64> {
    let ne = problem.topology.num_edges;
    let n_lat = variable_supports::latent_dim(problem);
    let mut x = Vec::with_capacity(ne + n_lat);
    for (i, span) in q_spans.iter().enumerate().take(ne) {
        let u = (state.force_densities[i] - problem.bounds.lower[i]) / span;
        x.push(u.clamp(0.0, 1.0));
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
    for i in 0..n_lat {
        let scale = anchor_scales[i].max(DIRECT_BOX_SCALE_EPS);
        x.push(latents.get(i).copied().unwrap_or(0.0) / scale);
    }
    x
}

fn direct_box_scaled_to_physical(
    problem: &Problem,
    x: &[f64],
    q_spans: &[f64],
    anchor_scales: &[f64],
) -> Vec<f64> {
    let ne = problem.topology.num_edges;
    let n_lat = variable_supports::latent_dim(problem);
    let mut theta = Vec::with_capacity(ne + n_lat);
    for i in 0..ne {
        let u = x[i].clamp(0.0, 1.0);
        theta.push(problem.bounds.lower[i] + q_spans[i] * u);
    }
    for i in 0..n_lat {
        theta.push(x[ne + i] * anchor_scales[i]);
    }
    theta
}

fn direct_box_scaled_gradient(
    grad_physical: &[f64],
    q_spans: &[f64],
    anchor_scales: &[f64],
) -> Vec<f64> {
    let ne = q_spans.len();
    let mut grad = vec![0.0; grad_physical.len()];
    for i in 0..ne {
        grad[i] = grad_physical[i] * q_spans[i];
    }
    for i in 0..anchor_scales.len() {
        grad[ne + i] = grad_physical[ne + i] * anchor_scales[i];
    }
    grad
}

fn direct_box_optimizer_bounds(problem: &Problem, n_lat: usize) -> (Vec<f64>, Vec<f64>) {
    let ne = problem.topology.num_edges;
    let mut lower = vec![0.0; ne];
    let mut upper = vec![1.0; ne];
    lower.extend(vec![f64::NEG_INFINITY; n_lat]);
    upper.extend(vec![f64::INFINITY; n_lat]);
    (lower, upper)
}

fn projected_gradient_norm(x: &[f64], grad: &[f64], lower: &[f64], upper: &[f64]) -> f64 {
    let mut norm = 0.0_f64;
    for i in 0..x.len() {
        let g = if x[i] <= lower[i] + DIRECT_BOX_SCALE_EPS && grad[i] > 0.0 {
            0.0
        } else if x[i] >= upper[i] - DIRECT_BOX_SCALE_EPS && grad[i] < 0.0 {
            0.0
        } else {
            grad[i]
        };
        norm = norm.max(g.abs());
    }
    norm
}

fn project_to_bounds(x: &mut [f64], lower: &[f64], upper: &[f64]) {
    for i in 0..x.len() {
        if lower[i].is_finite() && x[i] < lower[i] {
            x[i] = lower[i];
        }
        if upper[i].is_finite() && x[i] > upper[i] {
            x[i] = upper[i];
        }
    }
}

fn projected_gradient_fallback<F>(
    x: &mut Vec<f64>,
    lower: &[f64],
    upper: &[f64],
    max_iter: usize,
    pgtol: f64,
    f_and_grad: &mut F,
) -> (f64, usize, Status)
where
    F: FnMut(&[f64]) -> (f64, Vec<f64>),
{
    let (mut f, mut grad) = f_and_grad(x);
    if !f.is_finite() || grad.iter().any(|g| !g.is_finite()) {
        return (f, 0, Status::NumericalFailure);
    }

    for iter in 1..=max_iter {
        let pg_norm = projected_gradient_norm(x, &grad, lower, upper);
        if pg_norm <= pgtol {
            return (f, iter - 1, Status::Converged);
        }

        let mut direction = vec![0.0; x.len()];
        for i in 0..x.len() {
            if x[i] <= lower[i] + DIRECT_BOX_SCALE_EPS && grad[i] > 0.0 {
                continue;
            }
            if x[i] >= upper[i] - DIRECT_BOX_SCALE_EPS && grad[i] < 0.0 {
                continue;
            }
            direction[i] = -grad[i];
        }

        let max_dir = direction.iter().fold(0.0_f64, |acc, v| acc.max(v.abs()));
        if max_dir <= DIRECT_BOX_SCALE_EPS {
            return (f, iter - 1, Status::Converged);
        }
        if max_dir > 1.0 {
            for d in &mut direction {
                *d /= max_dir;
            }
        }

        let directional_derivative = grad
            .iter()
            .zip(direction.iter())
            .map(|(g, d)| g * d)
            .sum::<f64>();
        if directional_derivative >= 0.0 {
            return (f, iter - 1, Status::NumericalFailure);
        }

        let mut accepted = false;
        let mut step = 1.0;
        for _ in 0..DIRECT_BOX_MAX_BACKTRACK {
            let mut trial = x
                .iter()
                .zip(direction.iter())
                .map(|(xi, di)| xi + step * di)
                .collect::<Vec<_>>();
            project_to_bounds(&mut trial, lower, upper);
            if trial
                .iter()
                .zip(x.iter())
                .all(|(a, b)| (a - b).abs() <= DIRECT_BOX_SCALE_EPS)
            {
                step *= 0.5;
                continue;
            }

            let (trial_f, trial_grad) = f_and_grad(&trial);
            if trial_f.is_finite()
                && trial_grad.iter().all(|g| g.is_finite())
                && trial_f <= f + 1e-4 * step * directional_derivative
            {
                *x = trial;
                f = trial_f;
                grad = trial_grad;
                accepted = true;
                break;
            }
            step *= 0.5;
        }

        if !accepted {
            return (f, iter - 1, Status::LineSearchFailure);
        }
    }

    (f, max_iter, Status::MaxIter)
}

fn optimize_direct_box_bounds(
    problem: &Problem,
    state: &mut OptimizationState,
    progress_cb: Option<ProgressCallback>,
    report_freq: usize,
    cancel_flag: &AtomicBool,
) -> Result<SolverResult, TheseusError> {
    cancel_flag.store(false, Ordering::Relaxed);
    crate::objectives::validate_objectives(&problem.objectives)?;

    let ne = problem.topology.num_edges;
    let n_lat = variable_supports::latent_dim(problem);
    let q_spans = direct_box_spans(problem)?;
    let anchor_scales = anchor_optimizer_scales(problem);
    if anchor_scales.len() != n_lat {
        return Err(TheseusError::Shape(format!(
            "anchor optimizer scale length mismatch: got {}, expected {n_lat}",
            anchor_scales.len()
        )));
    }

    let mut x = pack_direct_box_scaled(problem, state, &q_spans, &anchor_scales);
    let (lower, upper) = direct_box_optimizer_bounds(problem, n_lat);
    let eval_lb = vec![f64::NEG_INFINITY; ne + n_lat];
    let eval_ub = vec![f64::INFINITY; ne + n_lat];
    let lb_idx = Vec::new();
    let ub_idx = Vec::new();
    let cache = Rc::new(RefCell::new(FdmCache::new(problem)?));
    let loss_trace = Rc::new(RefCell::new(Vec::new()));
    let pending_error: Rc<RefCell<Option<TheseusError>>> = Rc::new(RefCell::new(None));

    let q_spans_eval = q_spans.clone();
    let anchor_scales_eval = anchor_scales.clone();
    let cache_eval = cache.clone();
    let loss_trace_eval = loss_trace.clone();
    let pending_error_eval = pending_error.clone();
    let mut f_and_grad = |x_scaled: &[f64]| -> (f64, Vec<f64>) {
        if cancel_flag.load(Ordering::Relaxed) {
            *pending_error_eval.borrow_mut() = Some(TheseusError::Cancelled);
            return (DIRECT_BOX_ERROR_COST, vec![0.0; x_scaled.len()]);
        }
        if x_scaled.iter().any(|v| !v.is_finite()) {
            *pending_error_eval.borrow_mut() = Some(TheseusError::Solver(
                "L-BFGS-B produced NaN or Inf parameters".into(),
            ));
            return (DIRECT_BOX_ERROR_COST, vec![0.0; x_scaled.len()]);
        }

        let theta =
            direct_box_scaled_to_physical(problem, x_scaled, &q_spans_eval, &anchor_scales_eval);
        let mut grad_physical = vec![0.0; theta.len()];
        let val = {
            let mut fdm_cache = cache_eval.borrow_mut();
            value_and_gradient(
                &mut fdm_cache,
                problem,
                &theta,
                &mut grad_physical,
                &eval_lb,
                &eval_ub,
                &lb_idx,
                &ub_idx,
            )
        };

        match val {
            Ok(loss) if loss.is_finite() && grad_physical.iter().all(|g| g.is_finite()) => {
                loss_trace_eval.borrow_mut().push(loss);
                (
                    loss,
                    direct_box_scaled_gradient(&grad_physical, &q_spans_eval, &anchor_scales_eval),
                )
            }
            Ok(_) => (DIRECT_BOX_ERROR_COST, vec![0.0; x_scaled.len()]),
            Err(_) => (DIRECT_BOX_ERROR_COST, vec![0.0; x_scaled.len()]),
        }
    };

    let report_frequency = report_freq.max(1);
    let q_spans_cb = q_spans.clone();
    let anchor_scales_cb = anchor_scales.clone();
    let cache_cb = cache.clone();
    let pending_error_cb = pending_error.clone();
    let mut callback = |info: &lbfgsb_rs_pure::IterationInfo, x_scaled: &[f64]| {
        if cancel_flag.load(Ordering::Relaxed) {
            *pending_error_cb.borrow_mut() = Some(TheseusError::Cancelled);
            return IterationControl::StopCustom;
        }

        let major_iteration = info.iteration;
        if progress_cb.is_none()
            || !should_report_major_iteration(major_iteration, report_frequency)
        {
            return IterationControl::Continue;
        }

        let theta =
            direct_box_scaled_to_physical(problem, x_scaled, &q_spans_cb, &anchor_scales_cb);
        let q = theta[..ne].to_vec();
        let anchors = match variable_supports::map_latents_to_positions(problem, &theta[ne..]) {
            Ok(anchors) => anchors,
            Err(err) => {
                *pending_error_cb.borrow_mut() = Some(err);
                return IterationControl::StopCustom;
            }
        };
        let xyz_flat = {
            let mut fdm_cache = cache_cb.borrow_mut();
            if let Err(err) =
                crate::fdm::solve_fdm_with_loads(&mut fdm_cache, &q, problem, &anchors, 1e-12)
            {
                *pending_error_cb.borrow_mut() = Some(err);
                return IterationControl::StopCustom;
            }
            flatten_xyz(&fdm_cache, problem.topology.num_nodes)
        };

        let should_continue = unsafe {
            progress_cb.unwrap()(
                major_iteration,
                info.f,
                xyz_flat.as_ptr(),
                problem.topology.num_nodes,
                q.as_ptr(),
                ne,
            )
        };
        if should_continue == 0 {
            *pending_error_cb.borrow_mut() = Some(TheseusError::Cancelled);
            IterationControl::StopCustom
        } else {
            IterationControl::Continue
        }
    };

    let mut solver = LBFGSB::new(10)
        .with_max_iter(problem.solver.max_iterations)
        .with_pgtol(problem.solver.absolute_tolerance);
    let solution = solver
        .minimize_with_callback(&mut x, &lower, &upper, &mut f_and_grad, &mut callback)
        .map_err(|e| TheseusError::Solver(e.to_string()))?;

    if let Some(err) = pending_error.borrow_mut().take() {
        return Err(err);
    }

    let mut best_x = if solution.x.len() == ne + n_lat {
        solution.x.clone()
    } else {
        x.clone()
    };
    let mut iterations = solution.iterations;
    let mut status = solution.status;
    if matches!(status, Status::LineSearchFailure | Status::NumericalFailure) {
        let remaining_iterations = problem
            .solver
            .max_iterations
            .saturating_sub(iterations)
            .max(1);
        let (_fallback_loss, fallback_iterations, fallback_status) = projected_gradient_fallback(
            &mut best_x,
            &lower,
            &upper,
            remaining_iterations,
            problem.solver.absolute_tolerance,
            &mut f_and_grad,
        );
        iterations += fallback_iterations;
        if matches!(fallback_status, Status::Converged | Status::MaxIter) || fallback_iterations > 0
        {
            status = fallback_status;
        }
    }

    if let Some(err) = pending_error.borrow_mut().take() {
        return Err(err);
    }

    let theta = direct_box_scaled_to_physical(problem, &best_x, &q_spans, &anchor_scales);
    let q = theta[..ne].to_vec();
    let latents = theta[ne..].to_vec();
    let anchors = variable_supports::map_latents_to_positions(problem, &latents)?;
    let mut final_cache = FdmCache::new(problem)?;
    crate::fdm::solve_fdm_with_loads(&mut final_cache, &q, problem, &anchors, 1e-12)?;

    let converged = status == Status::Converged;
    let termination_reason = if solution.status != status {
        format!(
            "{:?}; projected-gradient fallback: {:?}",
            solution.status, status
        )
    } else {
        format!("{:?}", status)
    };
    let trace = loss_trace.borrow().clone();

    state.force_densities = q.clone();
    state.variable_anchor_positions = anchors.clone();
    state.variable_anchor_latents = latents.clone();
    state.iterations = iterations;
    state.loss_trace = trace.clone();

    Ok(SolverResult {
        q,
        anchor_positions: anchors,
        xyz: final_cache.nf.clone(),
        member_lengths: final_cache.member_lengths.clone(),
        member_forces: final_cache.member_forces.clone(),
        reactions: final_cache.reactions.clone(),
        loss_trace: trace,
        iterations: state.iterations,
        converged,
        termination_reason,
        cross_section_areas: final_cache.cross_section_areas.clone(),
    })
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

    cancel_flag.store(false, Ordering::Relaxed);
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
                cache: cache.clone(),
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

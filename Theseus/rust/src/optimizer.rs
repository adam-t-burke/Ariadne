//! L-BFGS optimisation driver via the `argmin` crate.
//!
//! Wraps the hand-coded `value_and_gradient` into argmin's `CostFunction`
//! + `Gradient` traits, then runs L-BFGS with the user's solver options.
//!
//! Uses `Vec<f64>` as the argmin parameter type to avoid ndarray version
//! conflicts between our ndarray 0.16 and argmin-math's bundled ndarray.

use crate::ffi::ProgressCallback;
use crate::gradients::value_and_gradient;
use crate::types::{FdmCache, Problem, SolverResult, OptimizationState, TheseusError};
use argmin::core::{CostFunction, Gradient, Executor, State, TerminationReason, TerminationStatus};
use argmin::solver::linesearch::MoreThuenteLineSearch;
use argmin::solver::quasinewton::LBFGS;
use ndarray::Array2;
use std::cell::RefCell;

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
    cache: RefCell<FdmCache>,
    lb: Vec<f64>,
    ub: Vec<f64>,
    lb_idx: Vec<usize>,
    ub_idx: Vec<usize>,
    /// Cached (θ, loss, gradient) from the last evaluation.
    last_eval: RefCell<Option<(Vec<f64>, f64, Vec<f64>)>>,
    /// Loss value recorded at each unique evaluation.
    loss_trace: RefCell<Vec<f64>>,
    /// Optional FFI callback for progress reporting.
    progress_callback: Option<ProgressCallback>,
    /// How often (in evaluations) to invoke the callback.
    report_frequency: usize,
}

impl<'a> FdmProblem<'a> {
    /// Ensure the cache contains results for `theta`.
    /// If θ matches the cached value, this is a no-op.
    /// Otherwise, runs the full forward + adjoint solve.
    fn ensure_evaluated(&self, theta: &[f64]) -> Result<(), argmin::core::Error> {
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
        ).map_err(|e| argmin::core::Error::msg(e.to_string()))?;

        // Guard against NaN/Inf in loss or gradient
        if !val.is_finite() || grad.iter().any(|g| !g.is_finite()) {
            return Err(argmin::core::Error::msg(
                "value_and_gradient produced NaN or Inf",
            ));
        }

        let eval_count = {
            let mut trace = self.loss_trace.borrow_mut();
            trace.push(val);
            trace.len()
        };

        if let Some(cb) = self.progress_callback {
            if eval_count == 1 || eval_count % self.report_frequency == 0 {
                let nn = self.problem.topology.num_nodes;
                let nf = &fdm_cache.nf;
                let xyz_flat: Vec<f64> = (0..nn)
                    .flat_map(|i| (0..3).map(move |d| nf[[i, d]]))
                    .collect();
                let should_continue = unsafe { cb(eval_count, val, xyz_flat.as_ptr(), nn) };
                if should_continue == 0 {
                    return Err(argmin::core::Error::msg("cancelled"));
                }
            }
        }

        *self.last_eval.borrow_mut() = Some((theta.to_vec(), val, grad));
        Ok(())
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

/// Pack q and anchor positions into a single θ vector.
pub fn pack_parameters(problem: &Problem, state: &OptimizationState) -> Vec<f64> {
    let ne = problem.topology.num_edges;
    let nvar = problem.anchors.variable_indices.len();
    let mut theta = Vec::with_capacity(ne + nvar * 3);
    theta.extend_from_slice(&state.force_densities);
    if nvar > 0 {
        for i in 0..nvar {
            theta.push(state.variable_anchor_positions[[i, 0]]);
            theta.push(state.variable_anchor_positions[[i, 1]]);
            theta.push(state.variable_anchor_positions[[i, 2]]);
        }
    }
    theta
}

/// Unpack θ into q and anchor positions.
pub fn unpack_parameters(problem: &Problem, theta: &[f64]) -> (Vec<f64>, Array2<f64>) {
    let ne = problem.topology.num_edges;
    let q = theta[..ne].to_vec();
    let nvar = problem.anchors.variable_indices.len();
    let anchors = if nvar > 0 {
        let mut a = Array2::zeros((nvar, 3));
        for i in 0..nvar {
            a[[i, 0]] = theta[ne + i * 3];
            a[[i, 1]] = theta[ne + i * 3 + 1];
            a[[i, 2]] = theta[ne + i * 3 + 2];
        }
        a
    } else {
        Array2::zeros((0, 3))
    };
    (q, anchors)
}

// ─────────────────────────────────────────────────────────────
//  Bound index precomputation
// ─────────────────────────────────────────────────────────────

fn parameter_bounds(problem: &Problem) -> (Vec<f64>, Vec<f64>) {
    let nvar = problem.anchors.variable_indices.len();
    let mut lb = problem.bounds.lower.clone();
    let mut ub = problem.bounds.upper.clone();
    if nvar > 0 {
        lb.extend(vec![f64::NEG_INFINITY; nvar * 3]);
        ub.extend(vec![f64::INFINITY; nvar * 3]);
    }
    (lb, ub)
}

fn finite_indices(v: &[f64]) -> Vec<usize> {
    v.iter().enumerate().filter(|(_, &x)| x.is_finite()).map(|(i, _)| i).collect()
}

// ─────────────────────────────────────────────────────────────
//  Top-level optimisation entry point
// ─────────────────────────────────────────────────────────────

/// Run L-BFGS optimisation on the FDM problem.
///
/// `progress_cb` / `report_freq` control an optional FFI callback invoked
/// every `report_freq` evaluations with the current node positions.
pub fn optimize(
    problem: &Problem,
    state: &mut OptimizationState,
    progress_cb: Option<ProgressCallback>,
    report_freq: usize,
) -> Result<SolverResult, TheseusError> {
    let cache = FdmCache::new(problem)?;

    let (lb, ub) = parameter_bounds(problem);
    let lb_idx = finite_indices(&lb);
    let ub_idx = finite_indices(&ub);

    let init_param = pack_parameters(problem, state);

    let fdm_problem = FdmProblem {
        problem,
        cache: RefCell::new(cache),
        lb,
        ub,
        lb_idx,
        ub_idx,
        last_eval: RefCell::new(None),
        loss_trace: RefCell::new(Vec::new()),
        progress_callback: progress_cb,
        report_frequency: if report_freq == 0 { 1 } else { report_freq },
    };

    // Configure L-BFGS with user-specified tolerances
    let linesearch = MoreThuenteLineSearch::new();
    let solver = LBFGS::new(linesearch, 10)
        .with_tolerance_grad(problem.solver.absolute_tolerance)
        .map_err(|e| TheseusError::Solver(format!("tolerance_grad: {e}")))?
        .with_tolerance_cost(problem.solver.relative_tolerance)
        .map_err(|e| TheseusError::Solver(format!("tolerance_cost: {e}")))?;

    let executor = Executor::new(fdm_problem, solver)
        .configure(|config| {
            config
                .param(init_param)
                .max_iters(problem.solver.max_iterations as u64)
                .target_cost(f64::NEG_INFINITY)
        });

    let result = executor.run()?;
    let loss_trace = result.problem.problem
        .as_ref()
        .map(|p| p.loss_trace.borrow().clone())
        .unwrap_or_default();

    // Extract solution
    let best_param = result.state().get_best_param()
        .ok_or_else(|| TheseusError::Solver("L-BFGS returned no best parameters".into()))?;
    let (q, anchors) = unpack_parameters(problem, best_param);

    // Final forward solve to get geometry
    let mut final_cache = FdmCache::new(problem)?;
    crate::fdm::solve_fdm(&mut final_cache, &q, problem, &anchors, 1e-12)?;
    crate::fdm::compute_geometry(&mut final_cache, problem);

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
    state.iterations = result.state().get_iter() as usize;
    state.loss_trace = loss_trace.clone();

    Ok(SolverResult {
        q,
        anchor_positions: anchors,
        xyz: final_cache.nf,
        member_lengths: final_cache.member_lengths,
        member_forces: final_cache.member_forces,
        reactions: final_cache.reactions,
        loss_trace,
        iterations: state.iterations,
        converged,
        termination_reason,
    })
}

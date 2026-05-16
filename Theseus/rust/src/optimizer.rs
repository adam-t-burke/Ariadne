//! L-BFGS optimisation driver via the `argmin` crate.
//!
//! Wraps the hand-coded `value_and_gradient` into argmin's `CostFunction`
//! + `Gradient` traits, then runs L-BFGS with the user's solver options.
//!
//! Uses `Vec<f64>` as the argmin parameter type to avoid ndarray version
//! conflicts between our ndarray 0.16 and argmin-math's bundled ndarray.

use crate::ffi::ProgressCallback;
use crate::gradients::value_and_gradient;
use crate::types::{FdmCache, OptimizationState, Problem, SolverResult, TheseusError};
use crate::variable_supports;
use argmin::core::observers::{Observe, ObserverMode};
use argmin::core::{
    CostFunction, Error, Executor, Gradient, IterState, State, TerminationReason,
    TerminationStatus, KV,
};
use argmin::solver::linesearch::MoreThuenteLineSearch;
use argmin::solver::quasinewton::LBFGS;
use std::cell::RefCell;
use std::rc::Rc;

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
        Ok(())
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

        let anchors = variable_supports::map_latents_to_positions(problem, &theta[ne..])
            .map_err(|e| Error::msg(e.to_string()))?;
        let mut cache = FdmCache::new(problem).map_err(|e| Error::msg(e.to_string()))?;
        crate::fdm::solve_fdm_with_loads(&mut cache, &theta[..ne], problem, &anchors, 1e-12)
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
        let q = &theta[..ne];
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
            return Err(Error::msg("cancelled"));
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
    use super::should_report_major_iteration;

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
    theta.extend_from_slice(&state.force_densities);
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
    let mut lb = problem.bounds.lower.clone();
    let mut ub = problem.bounds.upper.clone();
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
) -> Result<SolverResult, TheseusError> {
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
    };

    // Configure L-BFGS with user-specified tolerances
    let linesearch = MoreThuenteLineSearch::new();
    let solver = LBFGS::new(linesearch, 10)
        .with_tolerance_grad(problem.solver.absolute_tolerance)
        .map_err(|e| TheseusError::Solver(format!("tolerance_grad: {e}")))?
        .with_tolerance_cost(problem.solver.relative_tolerance)
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
                cache,
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

    // Final forward solve to get geometry (with self-weight/pressure if active)
    let mut final_cache = FdmCache::new(problem)?;
    crate::fdm::solve_fdm_with_loads(&mut final_cache, &q, problem, &anchors, 1e-12)?;

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

//! L-BFGS-B optimisation driver via the `lbfgsb-rs-pure` crate.
//!
//! Wraps the hand-coded `value_and_gradient` into a closure for the L-BFGS-B
//! solver, which directly supports box constraints on parameters.
//!
//! Key robustness features:
//! - Failed forward solves return a large finite penalty (not NaN) so the
//!   line search can backtrack instead of dying.
//! - Convergence is controlled entirely via the iteration callback, checking
//!   both projected-gradient norm AND relative function decrease.
//! - On LineSearchFailure / NumericalFailure the solver restarts from the
//!   best known point with fresh L-BFGS memory (up to MAX_RESTARTS times).

use crate::ffi::ProgressCallback;
use crate::gradients::value_and_gradient;
use crate::types::{FdmCache, Problem, SolverResult, OptimizationState, TheseusError};
use lbfgsb_rs_pure::{LBFGSB, IterationControl};
use ndarray::Array2;
use std::cell::RefCell;

const MAX_RESTARTS: usize = 3;
const MIN_ITERATIONS_BEFORE_CONVERGENCE: usize = 10;
const CONVERGENCE_WINDOW: usize = 5;

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
//  Bound precomputation
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

/// Clamp `x` into the feasible box `[lb, ub]`.
fn project_to_bounds(x: &mut [f64], lb: &[f64], ub: &[f64]) {
    for i in 0..x.len() {
        if x[i] < lb[i] { x[i] = lb[i]; }
        if x[i] > ub[i] { x[i] = ub[i]; }
    }
}

/// Apply a small deterministic perturbation to q-parameters that are strictly
/// interior to their bounds, nudging them toward the midpoint.  This helps
/// the solver escape a stale point after a restart.
fn perturb_interior(x: &mut [f64], lb: &[f64], ub: &[f64], ne: usize, strength: f64) {
    for i in 0..ne {
        let lo = lb[i];
        let hi = ub[i];
        if hi.is_infinite() || lo.is_infinite() { continue; }
        let range = hi - lo;
        if range <= 0.0 { continue; }
        let mid = (lo + hi) * 0.5;
        let nudge = (mid - x[i]) * strength;
        x[i] = (x[i] + nudge).max(lo).min(hi);
    }
}

// ─────────────────────────────────────────────────────────────
//  Top-level optimisation entry point
// ─────────────────────────────────────────────────────────────

/// Run L-BFGS-B optimisation on the FDM problem.
///
/// `progress_cb` / `report_freq` control an optional FFI callback invoked
/// every `report_freq` evaluations with the current node positions.
pub fn optimize(
    problem: &Problem,
    state: &mut OptimizationState,
    progress_cb: Option<ProgressCallback>,
    report_freq: usize,
) -> Result<SolverResult, TheseusError> {
    let report_freq = if report_freq == 0 { 1 } else { report_freq };
    let (lb, ub) = parameter_bounds(problem);
    let ne = problem.topology.num_edges;

    let mut x = pack_parameters(problem, state);
    project_to_bounds(&mut x, &lb, &ub);

    let mut global_best_x = x.clone();
    let mut global_best_f = f64::INFINITY;
    let mut all_traces: Vec<f64> = Vec::new();
    let mut total_iterations: usize = 0;
    let mut final_status = String::from("MaxIter");
    let mut was_cancelled = false;

    let abs_tol = problem.solver.absolute_tolerance;
    let rel_tol = problem.solver.relative_tolerance;
    let max_iter = problem.solver.max_iterations;

    for restart in 0..=MAX_RESTARTS {
        if restart > 0 {
            x = global_best_x.clone();
            let strength = 0.02 * (restart as f64);
            perturb_interior(&mut x, &lb, &ub, ne, strength);
        }

        let remaining_iter = max_iter.saturating_sub(total_iterations);
        if remaining_iter < 5 {
            break;
        }

        let cache = RefCell::new(FdmCache::new(problem)?);
        let loss_trace = RefCell::new(Vec::<f64>::new());
        let cancelled = RefCell::new(false);
        let best_x = RefCell::new(x.clone());
        let best_f = RefCell::new(f64::INFINITY);
        let last_valid_grad = RefCell::new(vec![0.0; x.len()]);
        let recent_f = RefCell::new(Vec::<f64>::with_capacity(CONVERGENCE_WINDOW + 1));

        // Disable the library's internal pgtol so we control convergence
        // entirely from the callback (Issue 2 fix).
        let mut solver = LBFGSB::new(10)
            .with_pgtol(0.0)
            .with_max_iter(remaining_iter);

        let solution_res = solver.minimize_with_callback(
            &mut x,
            &lb,
            &ub,
            // ── Objective closure ────────────────────────────
            &mut |theta: &[f64]| {
                let mut fdm_cache = cache.borrow_mut();
                let mut grad = vec![0.0; theta.len()];

                match value_and_gradient(&mut fdm_cache, problem, theta, &mut grad) {
                    Ok(val) => {
                        loss_trace.borrow_mut().push(val);
                        *last_valid_grad.borrow_mut() = grad.clone();
                        if val < *best_f.borrow() {
                            *best_f.borrow_mut() = val;
                            *best_x.borrow_mut() = theta.to_vec();
                        }
                        (val, grad)
                    }
                    Err(_) => {
                        // Issue 1 fix: return a large finite penalty so the
                        // line search sees "very bad point" and backtracks,
                        // instead of NaN which poisons dcstep permanently.
                        let penalty = best_f.borrow().abs().max(1.0) * 1e6;
                        let fallback_grad = last_valid_grad.borrow().clone();
                        (penalty, fallback_grad)
                    }
                }
            },
            // ── Iteration callback ───────────────────────────
            &mut |info, theta| {
                let eval_count = info.n_func_evals;
                let val = info.f;

                // Track rolling window for relative decrease check
                {
                    let mut rf = recent_f.borrow_mut();
                    rf.push(val);
                    if rf.len() > CONVERGENCE_WINDOW {
                        rf.remove(0);
                    }
                }

                // Issue 2 fix: callback-based convergence checking both
                // projected gradient AND relative function decrease.
                let iter_so_far = total_iterations + info.iteration;
                if iter_so_far >= MIN_ITERATIONS_BEFORE_CONVERGENCE
                    && info.proj_grad_norm <= abs_tol
                {
                    let rf = recent_f.borrow();
                    if rf.len() >= 2 {
                        let oldest = rf[0];
                        let newest = *rf.last().unwrap();
                        let denom = oldest.abs().max(newest.abs()).max(1.0);
                        let rel_change = (oldest - newest).abs() / denom;
                        if rel_change < rel_tol {
                            return IterationControl::StopConverged;
                        }
                    }
                }

                // Progress reporting via FFI callback
                if let Some(cb) = progress_cb {
                    if eval_count == 1 || eval_count % report_freq == 0 {
                        let nn = problem.topology.num_nodes;
                        let ne = problem.topology.num_edges;
                        let fdm_cache = cache.borrow();
                        let nf = &fdm_cache.nf;
                        let xyz_flat: Vec<f64> = (0..nn)
                            .flat_map(|i| (0..3).map(move |d| nf[[i, d]]))
                            .collect();
                        let q = &theta[..ne];

                        let should_continue = unsafe {
                            cb(eval_count, val, xyz_flat.as_ptr(), nn, q.as_ptr(), ne)
                        };

                        if should_continue == 0 {
                            *cancelled.borrow_mut() = true;
                            return IterationControl::StopCustom;
                        }
                    }
                }
                IterationControl::Continue
            },
        );

        if *cancelled.borrow() {
            was_cancelled = true;
        }

        // Harvest results from this run
        let local_best_x = best_x.into_inner();
        let local_best_f = *best_f.borrow();
        all_traces.append(&mut loss_trace.into_inner());

        if local_best_f < global_best_f {
            global_best_f = local_best_f;
            global_best_x = local_best_x;
        }

        match &solution_res {
            Ok(sol) => {
                total_iterations += sol.iterations;
                final_status = format!("{:?}", sol.status);
                if final_status.contains("Converged") || was_cancelled {
                    break;
                }
                // LineSearchFailure or MaxIter — try restart (Issue 3 fix)
            }
            Err(s) => {
                final_status = s.to_string();
                // Library returned Err — try restart
            }
        }

        if was_cancelled {
            break;
        }
    }

    if was_cancelled {
        return Err(TheseusError::Cancelled);
    }

    let (q, anchors) = unpack_parameters(problem, &global_best_x);

    // Final forward solve to get geometry at the best point
    let mut final_cache = FdmCache::new(problem)?;
    crate::fdm::solve_fdm(&mut final_cache, &q, problem, &anchors, 1e-12)?;
    crate::fdm::compute_geometry(&mut final_cache, problem);

    let converged = final_status.contains("Converged");

    state.force_densities = q.clone();
    state.variable_anchor_positions = anchors.clone();
    state.iterations = if total_iterations > 0 { total_iterations } else { all_traces.len() };
    state.loss_trace = all_traces.clone();

    Ok(SolverResult {
        q,
        anchor_positions: anchors,
        xyz: final_cache.nf,
        member_lengths: final_cache.member_lengths,
        member_forces: final_cache.member_forces,
        reactions: final_cache.reactions,
        loss_trace: all_traces,
        iterations: state.iterations,
        converged,
        termination_reason: final_status,
    })
}

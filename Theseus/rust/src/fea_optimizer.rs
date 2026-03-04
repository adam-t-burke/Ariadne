//! L-BFGS optimisation driver for FEA problems via `argmin`.
//!
//! Mirrors the FDM optimizer pattern: wraps `fea_value_and_gradient` into
//! argmin's CostFunction + Gradient traits, runs L-BFGS with user options.

use crate::fea_gradients::fea_value_and_gradient;
use crate::fea_solve::solve_fea;
use crate::fea_types::{
    FeaCache, FeaOptimizationState, FeaProblem, FeaSolverResult, FeaVariableMask,
};
use crate::types::TheseusError;
use argmin::core::{CostFunction, Executor, Gradient, State, TerminationReason, TerminationStatus};
use argmin::solver::linesearch::MoreThuenteLineSearch;
use argmin::solver::quasinewton::LBFGS;
use ndarray::Array2;
use std::cell::RefCell;

// ─────────────────────────────────────────────────────────────
//  Progress callback type
// ─────────────────────────────────────────────────────────────

pub type FeaProgressCallback = unsafe extern "C" fn(
    iteration: usize,
    loss: f64,
    xyz: *const f64,
    num_nodes: usize,
) -> u8;

// ─────────────────────────────────────────────────────────────
//  argmin problem wrapper
// ─────────────────────────────────────────────────────────────

struct FeaOptProblem<'a> {
    problem: &'a FeaProblem,
    cache: RefCell<FeaCache>,
    mask: &'a FeaVariableMask,
    lb: Vec<f64>,
    ub: Vec<f64>,
    lb_idx: Vec<usize>,
    ub_idx: Vec<usize>,
    last_eval: RefCell<Option<(Vec<f64>, f64, Vec<f64>)>>,
    loss_trace: RefCell<Vec<f64>>,
    progress_callback: Option<FeaProgressCallback>,
    report_frequency: usize,
}

impl<'a> FeaOptProblem<'a> {
    fn ensure_evaluated(&self, theta: &[f64]) -> Result<(), argmin::core::Error> {
        {
            let cached = self.last_eval.borrow();
            if let Some((ref t, _, _)) = *cached {
                if t == theta {
                    return Ok(());
                }
            }
        }

        if theta.iter().any(|v| !v.is_finite()) {
            return Err(argmin::core::Error::msg("theta contains NaN or Inf"));
        }

        let mut fea_cache = self.cache.borrow_mut();
        let mut grad = vec![0.0; theta.len()];
        let val = fea_value_and_gradient(
            &mut fea_cache,
            self.problem,
            theta,
            &mut grad,
            self.mask,
            &self.lb,
            &self.ub,
            &self.lb_idx,
            &self.ub_idx,
        )
        .map_err(|e| argmin::core::Error::msg(e.to_string()))?;

        if !val.is_finite() || grad.iter().any(|g| !g.is_finite()) {
            return Err(argmin::core::Error::msg(
                "fea_value_and_gradient produced NaN or Inf",
            ));
        }

        let eval_count = {
            let mut trace = self.loss_trace.borrow_mut();
            trace.push(val);
            trace.len()
        };

        if let Some(cb) = self.progress_callback {
            if eval_count == 1 || eval_count % self.report_frequency == 0 {
                // Reconstruct current positions from theta for preview
                let mut positions = self.problem.node_positions.clone();
                let mut offset = 0;
                if self.mask.node_positions {
                    for &ni in &self.mask.node_indices {
                        for d in 0..3 {
                            positions[[ni, d]] = theta[offset];
                            offset += 1;
                        }
                    }
                }
                let nn = self.problem.num_nodes;
                let pos_ref = &positions;
                let xyz_flat: Vec<f64> = (0..nn)
                    .flat_map(|i| (0..3).map(move |d| pos_ref[[i, d]]))
                    .collect();
                let should_continue =
                    unsafe { cb(eval_count, val, xyz_flat.as_ptr(), nn) };
                if should_continue == 0 {
                    return Err(argmin::core::Error::msg("cancelled"));
                }
            }
        }

        *self.last_eval.borrow_mut() = Some((theta.to_vec(), val, grad));
        Ok(())
    }
}

impl<'a> CostFunction for FeaOptProblem<'a> {
    type Param = Vec<f64>;
    type Output = f64;

    fn cost(&self, theta: &Self::Param) -> Result<Self::Output, argmin::core::Error> {
        self.ensure_evaluated(theta)?;
        let cached = self.last_eval.borrow();
        Ok(cached.as_ref().unwrap().1)
    }
}

impl<'a> Gradient for FeaOptProblem<'a> {
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

pub fn pack_fea_parameters(
    _problem: &FeaProblem,
    state: &FeaOptimizationState,
    mask: &FeaVariableMask,
) -> Vec<f64> {
    let mut theta = Vec::new();
    if mask.node_positions {
        for &ni in &mask.node_indices {
            for d in 0..3 {
                theta.push(state.node_positions[[ni, d]]);
            }
        }
    }
    if mask.cross_section_areas {
        for &ei in &mask.area_element_indices {
            theta.push(state.areas[ei]);
        }
    }
    if mask.support_positions {
        for &ni in &mask.support_node_indices {
            for d in 0..3 {
                theta.push(state.node_positions[[ni, d]]);
            }
        }
    }
    theta
}

pub fn unpack_fea_parameters(
    _problem: &FeaProblem,
    theta: &[f64],
    mask: &FeaVariableMask,
    base_positions: &Array2<f64>,
    base_areas: &[f64],
) -> (Array2<f64>, Vec<f64>) {
    let mut positions = base_positions.clone();
    let mut areas = base_areas.to_vec();
    let mut offset = 0;

    if mask.node_positions {
        for &ni in &mask.node_indices {
            for d in 0..3 {
                positions[[ni, d]] = theta[offset];
                offset += 1;
            }
        }
    }
    if mask.cross_section_areas {
        for &ei in &mask.area_element_indices {
            areas[ei] = theta[offset];
            offset += 1;
        }
    }
    if mask.support_positions {
        for &ni in &mask.support_node_indices {
            for d in 0..3 {
                positions[[ni, d]] = theta[offset];
                offset += 1;
            }
        }
    }
    (positions, areas)
}

fn parameter_bounds(mask: &FeaVariableMask) -> (Vec<f64>, Vec<f64>) {
    let mut lb = Vec::new();
    let mut ub = Vec::new();

    if mask.node_positions {
        let n = mask.node_indices.len() * 3;
        lb.extend(vec![f64::NEG_INFINITY; n]);
        ub.extend(vec![f64::INFINITY; n]);
    }
    if mask.cross_section_areas {
        let n = mask.area_element_indices.len();
        lb.extend(vec![1e-8; n]); // areas must be positive
        ub.extend(vec![f64::INFINITY; n]);
    }
    if mask.support_positions {
        let n = mask.support_node_indices.len() * 3;
        lb.extend(vec![f64::NEG_INFINITY; n]);
        ub.extend(vec![f64::INFINITY; n]);
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

pub fn fea_optimize(
    problem: &FeaProblem,
    state: &mut FeaOptimizationState,
    mask: &FeaVariableMask,
    progress_cb: Option<FeaProgressCallback>,
    report_freq: usize,
) -> Result<FeaSolverResult, TheseusError> {
    let cache = FeaCache::new(problem)?;

    let (lb, ub) = parameter_bounds(mask);
    let lb_idx = finite_indices(&lb);
    let ub_idx = finite_indices(&ub);

    let init_param = pack_fea_parameters(problem, state, mask);

    let fea_problem = FeaOptProblem {
        problem,
        cache: RefCell::new(cache),
        mask,
        lb,
        ub,
        lb_idx,
        ub_idx,
        last_eval: RefCell::new(None),
        loss_trace: RefCell::new(Vec::new()),
        progress_callback: progress_cb,
        report_frequency: if report_freq == 0 { 1 } else { report_freq },
    };

    let linesearch = MoreThuenteLineSearch::new();
    let solver = LBFGS::new(linesearch, 10)
        .with_tolerance_grad(problem.solver.absolute_tolerance)
        .map_err(|e| TheseusError::Solver(format!("tolerance_grad: {e}")))?
        .with_tolerance_cost(problem.solver.relative_tolerance)
        .map_err(|e| TheseusError::Solver(format!("tolerance_cost: {e}")))?;

    let executor = Executor::new(fea_problem, solver).configure(|config| {
        config
            .param(init_param)
            .max_iters(problem.solver.max_iterations as u64)
            .target_cost(f64::NEG_INFINITY)
    });

    let result = executor.run()?;
    let loss_trace = result
        .problem
        .problem
        .as_ref()
        .map(|p| p.loss_trace.borrow().clone())
        .unwrap_or_default();

    let best_param = result
        .state()
        .get_best_param()
        .ok_or_else(|| TheseusError::Solver("L-BFGS returned no best parameters".into()))?;

    let base_areas: Vec<f64> = problem.sections.iter().map(|s| s.area).collect();
    let (final_positions, final_areas) = unpack_fea_parameters(
        problem,
        best_param,
        mask,
        &problem.node_positions,
        &base_areas,
    );

    // Final forward solve
    let mut final_cache = FeaCache::new(problem)?;
    let fea_result = solve_fea(&mut final_cache, problem, &final_positions, &final_areas)?;

    let termination_status = result.state().get_termination_status();
    let converged = matches!(
        termination_status,
        TerminationStatus::Terminated(TerminationReason::SolverConverged)
    );
    let termination_reason = match termination_status {
        TerminationStatus::Terminated(reason) => format!("{reason}"),
        TerminationStatus::NotTerminated => "not terminated".to_string(),
    };

    state.node_positions = final_positions.clone();
    state.areas = final_areas.clone();
    state.iterations = result.state().get_iter() as usize;
    state.loss_trace = loss_trace.clone();

    Ok(FeaSolverResult {
        fea_result,
        node_positions: final_positions,
        areas: final_areas,
        loss_trace,
        iterations: state.iterations,
        converged,
        termination_reason,
    })
}

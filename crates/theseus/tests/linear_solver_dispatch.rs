//! The `FdmCache` linear-solver dispatch (WS-D): every forward, adjoint and
//! load-iteration solve goes through `Box<dyn LinearSystemSolver>`.
//!
//! * `Direct` must be bitwise identical to the pre-dispatch cache path —
//!   checked against `DirectSolver` driven stand-alone, exactly as the old
//!   `factor_and_solve` / `solve_adjoint` did.
//! * Statistics of every solve land in `FdmCache::linear_solver_totals`,
//!   `SolverResult::linear_solver_totals` and (per evaluation)
//!   `SolverResult::linear_solver_iterations`.
//! * Typed errors: a solve that does not converge is
//!   `IterativeSolverDidNotConverge` all the way up to the optimizer;
//!   cancellation inside a solve is `Cancelled`; bounds permitting `q ≤ 0`
//!   with an iterative kind are refused at `FdmCache::new`.
//! * The tolerance policy is validated at construction and the schedule
//!   behaves as §3 of the program plan specifies.
//! * Gradients match central finite differences at `Fixed(1e-12)` and CPU
//!   results are bitwise identical across thread counts, for every kind.
//!
//! Cases for `IterativeCpu` run once `LinearSolver::new` provides it and are
//! skipped with a printed reason until then (`support/linear_solver_kinds.rs`).

#[path = "support/grid.rs"]
#[allow(dead_code)]
mod grid;
#[path = "support/linear_solver_kinds.rs"]
#[allow(dead_code)]
mod kinds;

use kinds::{assert_kind_close, available, for_each_iterative_kind, for_each_kind};
use ndarray::Array2;
use std::sync::atomic::{AtomicBool, Ordering};
use theseus::linear_solver::{
    DirectSolver, LinearSolverKind, LinearSolverTotals, LinearSystemSolver, MemoryReport,
    SolveRequest, SolveStats, TolerancePolicy, ToleranceSchedule,
};
use theseus::types::*;

fn anchors() -> Array2<f64> {
    Array2::zeros((0, 3))
}

fn smooth_q(ne: usize, offset: f64) -> Vec<f64> {
    (0..ne)
        .map(|k| {
            let t = k as f64 / ne as f64;
            1.0 + 0.5 * (std::f64::consts::TAU * t + offset).sin()
        })
        .collect()
}

fn problem_for(kind: LinearSolverKind, n: usize) -> Problem {
    let mut problem = grid::make_recoverable_grid_problem(n);
    problem.solver.linear_solver = kind;
    problem
}

/// `value_and_gradient` on `cache` for `q` with the problem's soft bounds.
fn evaluate(cache: &mut FdmCache, problem: &Problem, q: &[f64], grad: &mut [f64]) -> f64 {
    let ne = problem.topology.num_edges;
    let idx: Vec<usize> = (0..ne).collect();
    theseus::gradients::value_and_gradient(
        cache,
        problem,
        q,
        grad,
        &problem.bounds.lower,
        &problem.bounds.upper,
        &idx,
        &idx,
    )
    .unwrap()
}

fn optimize(problem: &Problem, max_iterations: usize) -> SolverResult {
    let mut problem = problem.clone_shallow();
    problem.solver.q_parameterization_mode = QParameterizationMode::DirectBoxBounds;
    problem.solver.absolute_tolerance = 0.0;
    problem.solver.relative_tolerance = 0.0;
    problem.solver.max_iterations = max_iterations;
    let ne = problem.topology.num_edges;
    let mut state = OptimizationState::new(vec![1.0; ne], anchors());
    theseus::optimizer::optimize(&problem, &mut state, None, 1, &AtomicBool::new(false)).unwrap()
}

trait CloneShallow {
    fn clone_shallow(&self) -> Problem;
}

impl CloneShallow for Problem {
    /// A copy with the same topology, loads, bounds and options, and the
    /// same `TargetXYZ` objective (the fixtures only use that one).
    fn clone_shallow(&self) -> Problem {
        let mut problem = grid::make_recoverable_grid_problem(
            (self.topology.num_nodes as f64).sqrt().round() as usize,
        );
        problem.bounds = self.bounds.clone();
        problem.solver = self.solver.clone();
        problem
    }
}

// ─────────────────────────────────────────────────────────────
//  Direct: bitwise identical to the pre-dispatch path
// ─────────────────────────────────────────────────────────────

/// Drive `DirectSolver` stand-alone the way the pre-dispatch cache did
/// (assemble + perturb + refactor, three-column solve for `x`; the same
/// factorization for the adjoint `λ`) and compare with the cache bitwise.
fn assert_direct_matches_standalone(problem: &Problem, qs: &[Vec<f64>], perturbation: f64) {
    let ne = problem.topology.num_edges;
    let n = problem.topology.free_node_indices.len();
    let mut cache = FdmCache::new(problem).unwrap();
    assert_eq!(cache.linear_solver_kind, LinearSolverKind::Direct);
    assert!(cache.warm_x.is_empty() && cache.warm_lambda.is_empty());
    let mut solver = DirectSolver::from_bounds(&problem.topology, &problem.bounds).unwrap();
    solver.set_perturbation(perturbation);
    let mut x = vec![0.0; n * 3];
    let mut grad = vec![0.0; ne];
    let idx: Vec<usize> = (0..ne).collect();

    for (i, q) in qs.iter().enumerate() {
        // Forward solve through the dispatch.
        theseus::fdm::solve_fdm(&mut cache, q, problem, &anchors(), perturbation).unwrap();
        solver.update(q).unwrap();
        fn req(rhs: &[f64]) -> SolveRequest<'_> {
            SolveRequest {
                rhs,
                x0: None,
                tolerance: 0.0,
                max_iterations: 0,
                cancel: None,
            }
        }
        solver
            .solve(req(cache.rhs.as_slice().unwrap()), &mut x)
            .unwrap();
        assert_eq!(x, cache.x.as_slice().unwrap(), "forward x, solve {i}");
        assert_eq!(cache.strategy(), Some(solver.strategy()));
        assert_eq!(
            cache.a_matrix().unwrap().values,
            solver.a_matrix().values,
            "A(q) values, solve {i}"
        );

        // Full evaluation: the adjoint reuses the same factorization.
        let loss = theseus::gradients::value_and_gradient(
            &mut cache,
            problem,
            q,
            &mut grad,
            &problem.bounds.lower,
            &problem.bounds.upper,
            &idx,
            &idx,
        )
        .unwrap();
        assert!(loss.is_finite());
        // `value_and_gradient` always solves with the fixed 1e-12 perturbation.
        solver.set_perturbation(1e-12);
        solver.update(q).unwrap();
        solver
            .solve(req(cache.grad_x.as_slice().unwrap()), &mut x)
            .unwrap();
        assert_eq!(x, cache.lambda.as_slice().unwrap(), "adjoint λ, solve {i}");
        solver.set_perturbation(perturbation);
    }
}

#[test]
fn direct_dispatch_is_bitwise_identical_to_standalone_direct_solver() {
    let problem = grid::make_recoverable_grid_problem(12);
    let ne = problem.topology.num_edges;
    let qs = vec![vec![1.0; ne], smooth_q(ne, 0.0), smooth_q(ne, 1.3)];
    assert_direct_matches_standalone(&problem, &qs, 0.0);
    assert_direct_matches_standalone(&problem, &qs, 1e-12);

    // LDLᵀ strategy from mixed-sign bounds.
    let mut problem = grid::make_recoverable_grid_problem(9);
    let ne = problem.topology.num_edges;
    problem.bounds.lower = vec![-1.0; ne];
    let qs = vec![vec![1.0; ne], smooth_q(ne, 0.4)];
    assert_direct_matches_standalone(&problem, &qs, 1e-12);
}

#[test]
fn direct_optimization_is_reproducible_and_leaves_the_termination_string_alone() {
    let problem = problem_for(LinearSolverKind::Direct, 10);
    let a = optimize(&problem, 12);
    let b = optimize(&problem, 12);
    assert_eq!(a.loss_trace, b.loss_trace);
    assert_eq!(a.q, b.q);
    assert_eq!(a.xyz, b.xyz);
    assert_eq!(a.termination_reason, b.termination_reason);
    assert!(
        !a.termination_reason.contains("linear solver"),
        "Direct keeps the pre-dispatch termination string: {}",
        a.termination_reason
    );
    assert!(a
        .termination_reason
        .starts_with("stopped: maximum iterations reached; iterations=12; evaluations="));
    assert!(a.linear_solver_iterations.iter().all(|&i| i == 0));
    assert_eq!(a.linear_solver_iterations.len(), a.loss_trace.len());
}

// ─────────────────────────────────────────────────────────────
//  Totals recording
// ─────────────────────────────────────────────────────────────

#[test]
fn every_solve_is_recorded_in_the_cache_totals() {
    for_each_kind(|kind| {
        let problem = problem_for(kind, 8);
        let ne = problem.topology.num_edges;
        let mut cache = FdmCache::new(&problem).unwrap();
        assert_eq!(cache.linear_solver_kind, kind);
        assert_eq!(cache.linear_solver.kind(), kind);
        assert_eq!(cache.linear_solver_totals, LinearSolverTotals::new(kind));

        let q = smooth_q(ne, 0.2);
        theseus::fdm::solve_fdm(&mut cache, &q, &problem, &anchors(), 1e-12).unwrap();
        let t = cache.linear_solver_totals;
        assert_eq!(t.backend, kind);
        assert_eq!(t.solves, 1, "{kind}: one forward solve");
        assert!(t.converged_all);
        assert!(t.solve_ms_total >= 0.0 && t.setup_ms_total >= 0.0);
        if kind.is_iterative() {
            assert!(t.iterations_total >= 1, "{kind}: iterations are counted");
            assert_eq!(cache.has_warm_x, true);
            assert_eq!(cache.warm_x, cache.x.as_slice().unwrap());
        } else {
            assert_eq!(t.iterations_total, 0);
            assert_eq!(t.iterations_max, 0);
        }

        let mut grad = vec![0.0; ne];
        evaluate(&mut cache, &problem, &q, &mut grad);
        let t = cache.linear_solver_totals;
        assert_eq!(t.solves, 3, "{kind}: + forward + adjoint");
        if kind.is_iterative() {
            assert!(cache.has_warm_lambda);
            assert_eq!(cache.warm_lambda, cache.lambda.as_slice().unwrap());
        }

        // Geometry-dependent loads: the Picard loop and the Neumann adjoint
        // refinement run extra solves, all recorded.
        let mut problem = problem_for(kind, 8);
        problem.self_weight = Some(SelfWeightParams::Prescribed {
            linear_densities: vec![0.2; ne],
            gravity: [0.0, 0.0, -1.0],
            max_iters: 30,
            tolerance: 1e-10,
            relaxation: 1.0,
        });
        let mut cache = FdmCache::new(&problem).unwrap();
        evaluate(&mut cache, &problem, &q, &mut grad);
        assert!(
            cache.linear_solver_totals.solves > 2,
            "{kind}: load iteration solves are recorded ({})",
            cache.linear_solver_totals.solves
        );
        assert!(cache.linear_solver_totals.converged_all);
    });
}

#[test]
fn optimizer_reports_totals_and_per_evaluation_iterations() {
    for_each_kind(|kind| {
        let problem = problem_for(kind, 8);
        let result = optimize(&problem, 6);
        let totals = result.linear_solver_totals;
        assert_eq!(totals.backend, kind, "{kind}");
        assert!(totals.converged_all, "{kind}");
        let evaluations = result.loss_trace.len() as u64;
        assert!(
            totals.solves >= 2 * evaluations,
            "{kind}: at least forward + adjoint per evaluation ({} solves, {evaluations} evaluations)",
            totals.solves
        );
        assert_eq!(
            result.linear_solver_iterations.len(),
            result.loss_trace.len(),
            "{kind}: one entry per evaluation"
        );
        let per_evaluation: u64 = result
            .linear_solver_iterations
            .iter()
            .map(|&i| u64::from(i))
            .sum();
        assert!(
            per_evaluation <= totals.iterations_total,
            "{kind}: per-evaluation counts are a subset of the totals"
        );
        let summary = format!(
            "; linear solver: {kind}, {} solves, {} iterations",
            totals.solves, totals.iterations_total
        );
        if kind.is_iterative() {
            assert!(
                result.termination_reason.ends_with(&summary),
                "{kind}: {}",
                result.termination_reason
            );
            assert!(
                per_evaluation >= evaluations,
                "{kind}: every evaluation iterates"
            );
        } else {
            assert!(!result.termination_reason.contains("linear solver"));
            assert_eq!(per_evaluation, 0);
        }
    });
}

// ─────────────────────────────────────────────────────────────
//  Typed errors through the dispatch (scripted solver)
// ─────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum Script {
    /// Solve exactly but report `converged = false` as an iterative kind.
    NeverConverges,
    /// Raise the cancel flag and return `Cancelled` on the first solve.
    CancelsItself,
}

/// A `LinearSystemSolver` that delegates to `DirectSolver` but scripts its
/// outcome, standing in for an iterative solver to test the dispatch's
/// error handling before one exists.
struct Scripted {
    inner: DirectSolver,
    script: Script,
    solves: usize,
}

impl LinearSystemSolver for Scripted {
    fn update(&mut self, q: &[f64]) -> Result<(), TheseusError> {
        self.inner.update(q)
    }

    fn solve(&mut self, req: SolveRequest<'_>, x: &mut [f64]) -> Result<SolveStats, TheseusError> {
        self.solves += 1;
        match self.script {
            Script::NeverConverges => {
                let mut stats = self.inner.solve(req, x)?;
                stats.backend = LinearSolverKind::IterativeCpu;
                stats.iterations = [req.max_iterations; 3];
                stats.relative_residual = [1e-3, 2e-3, 3e-3];
                stats.converged = false;
                Ok(stats)
            }
            Script::CancelsItself => {
                if let Some(flag) = req.cancel {
                    flag.store(true, Ordering::Release);
                }
                Err(TheseusError::Cancelled)
            }
        }
    }

    fn kind(&self) -> LinearSolverKind {
        LinearSolverKind::IterativeCpu
    }

    fn memory_bytes(&self) -> MemoryReport {
        self.inner.memory_bytes()
    }
}

/// Replace the cache's solver by a scripted one posing as `IterativeCpu`.
fn script_cache(cache: &mut FdmCache, problem: &Problem, script: Script) {
    let n = problem.topology.free_node_indices.len();
    cache.linear_solver = Box::new(Scripted {
        inner: DirectSolver::from_bounds(&problem.topology, &problem.bounds).unwrap(),
        script,
        solves: 0,
    });
    cache.linear_solver_kind = LinearSolverKind::IterativeCpu;
    cache.linear_solver_totals = LinearSolverTotals::new(LinearSolverKind::IterativeCpu);
    cache.warm_x = vec![0.0; n * 3];
    cache.warm_lambda = vec![0.0; n * 3];
    cache.has_warm_x = false;
    cache.has_warm_lambda = false;
}

#[test]
fn unconverged_solve_is_a_hard_error_with_its_residual() {
    let problem = grid::make_recoverable_grid_problem(6);
    let ne = problem.topology.num_edges;
    let mut cache = FdmCache::new(&problem).unwrap();
    script_cache(&mut cache, &problem, Script::NeverConverges);
    cache.solve_max_iterations = 7;
    let q = smooth_q(ne, 0.0);

    let error = theseus::fdm::solve_fdm(&mut cache, &q, &problem, &anchors(), 1e-12).unwrap_err();
    match error {
        TheseusError::IterativeSolverDidNotConverge {
            iterations,
            relative_residual,
            kind,
        } => {
            assert_eq!(iterations, 7);
            assert_eq!(relative_residual, 3e-3, "largest column residual");
            assert_eq!(kind, LinearSolverKind::IterativeCpu);
        }
        other => panic!("expected IterativeSolverDidNotConverge, got {other}"),
    }
    // The failed solve is still accounted for.
    assert_eq!(cache.linear_solver_totals.solves, 1);
    assert!(!cache.linear_solver_totals.converged_all);
    assert_eq!(cache.linear_solver_totals.iterations_max, 7);

    // Through value_and_gradient (no silent continuation) …
    let idx: Vec<usize> = (0..ne).collect();
    let mut grad = vec![0.0; ne];
    let error = theseus::gradients::value_and_gradient(
        &mut cache,
        &problem,
        &q,
        &mut grad,
        &problem.bounds.lower,
        &problem.bounds.upper,
        &idx,
        &idx,
    )
    .unwrap_err();
    assert!(
        matches!(error, TheseusError::IterativeSolverDidNotConverge { .. }),
        "{error}"
    );

    // … and inside the load Newton iteration, which must not retry it as a
    // failed stage.
    let mut hydro = hydro_problem(5, 40.0);
    hydro.solver.linear_solver = LinearSolverKind::Direct;
    let mut cache = FdmCache::new(&hydro).unwrap();
    script_cache(&mut cache, &hydro, Script::NeverConverges);
    let q = vec![100.0; hydro.topology.num_edges];
    let error =
        theseus::fdm::solve_fdm_with_loads(&mut cache, &q, &hydro, &anchors(), 1e-12).unwrap_err();
    assert!(
        matches!(error, TheseusError::IterativeSolverDidNotConverge { .. }),
        "{error}"
    );
}

#[test]
fn cancellation_inside_a_solve_is_cancelled() {
    let problem = grid::make_recoverable_grid_problem(6);
    let ne = problem.topology.num_edges;
    let q = smooth_q(ne, 0.0);
    let idx: Vec<usize> = (0..ne).collect();
    let mut grad = vec![0.0; ne];

    let mut cache = FdmCache::new(&problem).unwrap();
    script_cache(&mut cache, &problem, Script::CancelsItself);
    let flag = AtomicBool::new(false);
    let error = theseus::gradients::value_and_gradient_cancellable(
        &mut cache,
        &problem,
        &q,
        &mut grad,
        &problem.bounds.lower,
        &problem.bounds.upper,
        &idx,
        &idx,
        Some(&flag),
    )
    .unwrap_err();
    assert!(matches!(error, TheseusError::Cancelled), "{error}");
    assert!(
        flag.load(Ordering::Acquire),
        "the solve saw the cancel flag"
    );

    // A pre-set flag stops every kind before or inside its first solve.
    for_each_kind(|kind| {
        let problem = problem_for(kind, 6);
        let mut cache = FdmCache::new(&problem).unwrap();
        let flag = AtomicBool::new(true);
        let error = theseus::fdm::solve_fdm_with_loads_cancellable(
            &mut cache,
            &q,
            &problem,
            &anchors(),
            1e-12,
            Some(&flag),
        )
        .unwrap_err();
        assert!(matches!(error, TheseusError::Cancelled), "{kind}: {error}");
    });
}

#[test]
fn iterative_max_iterations_one_does_not_converge() {
    for_each_iterative_kind(|kind| {
        let mut problem = problem_for(kind, 12);
        problem.solver.iterative.max_iterations = 1;
        problem.solver.iterative.tolerance = TolerancePolicy::Fixed(1e-12);
        let ne = problem.topology.num_edges;
        let mut cache = FdmCache::new(&problem).unwrap();
        let error =
            theseus::fdm::solve_fdm(&mut cache, &smooth_q(ne, 0.0), &problem, &anchors(), 1e-12)
                .unwrap_err();
        match error {
            TheseusError::IterativeSolverDidNotConverge {
                iterations,
                relative_residual,
                kind: reported,
            } => {
                assert_eq!(iterations, 1, "{kind}");
                assert!(relative_residual > 1e-12, "{kind}: carries the residual");
                assert_eq!(reported, kind);
            }
            other => panic!("{kind}: expected IterativeSolverDidNotConverge, got {other}"),
        }

        // The optimizer surfaces it unchanged.
        let mut problem = problem.clone_shallow();
        problem.solver.iterative.max_iterations = 1;
        problem.solver.iterative.tolerance = TolerancePolicy::Fixed(1e-12);
        problem.solver.q_parameterization_mode = QParameterizationMode::DirectBoxBounds;
        let mut state = OptimizationState::new(vec![1.0; ne], anchors());
        let error =
            theseus::optimizer::optimize(&problem, &mut state, None, 1, &AtomicBool::new(false))
                .unwrap_err();
        assert!(
            matches!(error, TheseusError::IterativeSolverDidNotConverge { .. }),
            "{kind}: {error}"
        );
        assert!(state.loss_trace.is_empty(), "{kind}: nothing published");
    });
}

// ─────────────────────────────────────────────────────────────
//  Construction-time checks
// ─────────────────────────────────────────────────────────────

#[test]
fn iterative_kinds_are_refused_for_bounds_permitting_non_positive_q() {
    for kind in [
        LinearSolverKind::IterativeCpu,
        LinearSolverKind::IterativeGpu,
    ] {
        for lower in [0.0, -1.0, f64::NEG_INFINITY] {
            let mut problem = problem_for(kind, 6);
            let ne = problem.topology.num_edges;
            problem.bounds.lower = vec![0.1; ne];
            problem.bounds.lower[ne / 2] = lower;
            // Authoritative at cache construction, before any factory call.
            match FdmCache::new(&problem) {
                Err(TheseusError::IterativeSolverUnsupported(msg)) => {
                    assert!(
                        msg.contains(&format!("edge {}", ne / 2)) && msg.contains("q > 0"),
                        "{kind} lower={lower}: {msg}"
                    );
                    let shown = TheseusError::IterativeSolverUnsupported(msg).to_string();
                    assert!(shown.contains("Direct"), "names the toggle: {shown}");
                }
                Err(other) => panic!("{kind} lower={lower}: unexpected error {other}"),
                Ok(_) => panic!("{kind} lower={lower}: must be refused"),
            }
            let mut state = OptimizationState::new(vec![1.0; ne], anchors());
            assert!(matches!(
                theseus::optimizer::optimize(
                    &problem,
                    &mut state,
                    None,
                    1,
                    &AtomicBool::new(false)
                ),
                Err(TheseusError::IterativeSolverUnsupported(_))
            ));
        }
    }
    // Direct accepts them (LDLᵀ).
    let mut problem = problem_for(LinearSolverKind::Direct, 6);
    problem.bounds.lower = vec![-1.0; problem.topology.num_edges];
    assert_eq!(
        FdmCache::new(&problem).unwrap().strategy(),
        Some(FactorizationStrategy::LDL)
    );
}

#[test]
fn tolerance_policy_is_validated_for_iterative_kinds_only() {
    for bad in [
        TolerancePolicy::Fixed(0.0),
        TolerancePolicy::Fixed(-1e-8),
        TolerancePolicy::Fixed(f64::NAN),
        TolerancePolicy::Adaptive {
            floor: 1e-6,
            ceiling: 1e-10,
            factor: 1e-2,
        },
        TolerancePolicy::Adaptive {
            floor: 1e-10,
            ceiling: 1e-6,
            factor: 0.0,
        },
    ] {
        assert!(bad.validate().is_err(), "{bad:?}");
        let mut problem = problem_for(LinearSolverKind::IterativeCpu, 6);
        problem.solver.iterative.tolerance = bad;
        assert!(
            matches!(FdmCache::new(&problem), Err(TheseusError::Shape(_))),
            "{bad:?} must be rejected at construction"
        );
        // Direct ignores the iterative options entirely.
        problem.solver.linear_solver = LinearSolverKind::Direct;
        FdmCache::new(&problem).unwrap();
    }
    for good in [TolerancePolicy::Fixed(1e-12), TolerancePolicy::default()] {
        good.validate().unwrap();
        let mut problem = problem_for(LinearSolverKind::Direct, 6);
        problem.solver.iterative.tolerance = good;
        let cache = FdmCache::new(&problem).unwrap();
        assert_eq!(cache.solve_tolerance, good.initial());
    }
}

// ─────────────────────────────────────────────────────────────
//  Tolerance schedule (public API)
// ─────────────────────────────────────────────────────────────

#[test]
fn tolerance_schedule_follows_the_projected_gradient_ratio() {
    let adaptive = TolerancePolicy::Adaptive {
        floor: 1e-10,
        ceiling: 1e-6,
        factor: 1e-2,
    };
    let mut schedule = ToleranceSchedule::new(adaptive);
    assert_eq!(schedule.current(), 1e-6, "first evaluation: ceiling");
    assert_eq!(schedule.reference(), None);
    // ‖g⁺_0‖ = 5: ratio 1 → 1e-2 → clamped to the ceiling.
    assert_eq!(schedule.observe(5.0), 1e-6);
    assert_eq!(schedule.reference(), Some(5.0));
    // ratio 1e-5 → 1e-7 inside the band.
    let tol = schedule.observe(5e-5);
    assert!((tol - 1e-7).abs() < 1e-21, "{tol}");
    assert_eq!(schedule.current(), tol);
    // ratio 1e-12 → below the floor.
    assert_eq!(schedule.observe(5e-12), 1e-10);
    // Growth (a rejected region) loosens again, never past the ceiling.
    assert_eq!(schedule.observe(500.0), 1e-6);
    // Degenerate norms fall back to the ceiling instead of panicking.
    assert_eq!(schedule.observe(f64::NAN), 1e-6);
    assert_eq!(adaptive.tolerance_for(1.0, 0.0), 1e-6);
    assert_eq!(adaptive.tolerance_for(f64::INFINITY, 1.0), 1e-6);
    schedule.reset();
    assert_eq!((schedule.reference(), schedule.current()), (None, 1e-6));

    let mut fixed = ToleranceSchedule::new(TolerancePolicy::Fixed(1e-12));
    assert_eq!(fixed.current(), 1e-12);
    assert_eq!(fixed.observe(3.0), 1e-12);
    assert_eq!(fixed.observe(0.0), 1e-12);
    assert_eq!(fixed.policy(), TolerancePolicy::Fixed(1e-12));
}

// ─────────────────────────────────────────────────────────────
//  Gradient correctness and determinism, per kind
// ─────────────────────────────────────────────────────────────

#[test]
fn gradient_matches_central_finite_differences_at_fixed_1e_12() {
    for_each_kind(|kind| {
        let mut problem = problem_for(kind, 6);
        problem.solver.iterative.tolerance = TolerancePolicy::Fixed(1e-12);
        problem.solver.max_iterations = 0;
        let ne = problem.topology.num_edges;
        let mut cache = FdmCache::new(&problem).unwrap();
        assert_eq!(cache.solve_tolerance, 1e-12);
        let q = smooth_q(ne, 0.9);
        let mut grad = vec![0.0; ne];
        evaluate(&mut cache, &problem, &q, &mut grad);

        let mut fd_cache = FdmCache::new(&problem).unwrap();
        let mut scratch = vec![0.0; ne];
        for k in (0..ne).step_by((ne / 12).max(1)) {
            let h = 1e-6 * q[k].abs().max(1.0);
            let mut plus = q.clone();
            plus[k] += h;
            let mut minus = q.clone();
            minus[k] -= h;
            let fp = evaluate(&mut fd_cache, &problem, &plus, &mut scratch);
            let fm = evaluate(&mut fd_cache, &problem, &minus, &mut scratch);
            let fd = (fp - fm) / (2.0 * h);
            let scale = grad[k].abs().max(fd.abs()).max(1e-3);
            assert!(
                (grad[k] - fd).abs() <= 1e-5 * scale,
                "{kind}: edge {k}: analytic {} vs finite difference {fd}",
                grad[k]
            );
        }
    });
}

#[test]
fn results_are_bitwise_identical_across_thread_counts() {
    for_each_kind(|kind| {
        let problem = problem_for(kind, 16);
        let ne = problem.topology.num_edges;
        let run = |threads: usize| {
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap()
                .install(|| {
                    let mut cache = FdmCache::new(&problem).unwrap();
                    let mut grad = vec![0.0; ne];
                    let loss = evaluate(&mut cache, &problem, &smooth_q(ne, 0.3), &mut grad);
                    let result = optimize(&problem, 6);
                    (loss, grad, cache.x, result)
                })
        };
        let (loss1, grad1, x1, r1) = run(1);
        let (loss4, grad4, x4, r4) = run(4);
        assert_eq!(loss1.to_bits(), loss4.to_bits(), "{kind}: loss");
        assert_eq!(grad1, grad4, "{kind}: gradient");
        assert_eq!(x1, x4, "{kind}: forward solution");
        assert_eq!(r1.loss_trace, r4.loss_trace, "{kind}: optimizer loss trace");
        assert_eq!(r1.q, r4.q, "{kind}: optimizer q");
        assert_eq!(r1.xyz, r4.xyz, "{kind}: optimizer xyz");
        assert_eq!(
            r1.linear_solver_iterations, r4.linear_solver_iterations,
            "{kind}: iteration counts"
        );
        assert_eq!(
            r1.linear_solver_totals.iterations_total, r4.linear_solver_totals.iterations_total,
            "{kind}"
        );
    });
}

#[test]
fn iterative_cpu_reaches_the_direct_loss_on_the_recoverable_grid() {
    if !available(LinearSolverKind::IterativeCpu) {
        return;
    }
    let direct = optimize(&problem_for(LinearSolverKind::Direct, 12), 30);
    let iterative = optimize(&problem_for(LinearSolverKind::IterativeCpu, 12), 30);
    let reference = *direct.loss_trace.last().unwrap();
    let reached = *iterative.loss_trace.last().unwrap();
    let initial = direct.loss_trace[0];
    // 1e-6 relative to the direct loss, with a floor relative to the initial
    // loss for the case where both are essentially at the (zero) optimum.
    let tolerance = 1e-6 * reference.abs().max(1e-6 * initial);
    assert!(
        (reached - reference).abs() <= tolerance,
        "IterativeCpu final loss {reached} vs Direct {reference} (initial {initial}); \
         iterations {} vs {}; {}",
        iterative.iterations,
        direct.iterations,
        iterative.termination_reason
    );
    assert!(iterative.linear_solver_totals.converged_all);
    assert!(iterative.linear_solver_totals.iterations_total > 0);
}

// ─────────────────────────────────────────────────────────────
//  Geometry-dependent loads through the dispatch
// ─────────────────────────────────────────────────────────────

/// Square grid with every boundary node fixed, faces clockwise from +Z, a
/// hydrostatic load and a target objective (the in-crate Newton fixture).
fn hydro_problem(size: usize, rho: f64) -> Problem {
    let mut problem = grid::make_grid_problem(size);
    let n = size;
    let mut fixed = Vec::new();
    let mut free = Vec::new();
    for row in 0..n {
        for col in 0..n {
            let node = row * n + col;
            if row == 0 || col == 0 || row + 1 == n || col + 1 == n {
                fixed.push(node);
            } else {
                free.push(node);
            }
        }
    }
    let incidence = problem.topology.incidence.clone();
    problem.topology.free_incidence = incidence.extract_columns(&free);
    problem.topology.fixed_incidence = incidence.extract_columns(&fixed);
    let fixed_positions = Array2::from_shape_fn((fixed.len(), 3), |(i, d)| {
        let node = fixed[i];
        match d {
            0 => (node % n) as f64,
            1 => (node / n) as f64,
            _ => 0.0,
        }
    });
    problem.topology.free_node_indices = free.clone();
    problem.topology.fixed_node_indices = fixed;
    problem.free_node_loads = Array2::zeros((free.len(), 3));
    problem.fixed_node_positions = fixed_positions.clone();
    problem.anchors = AnchorInfo::all_fixed(fixed_positions);
    let mut faces = Vec::new();
    for row in 0..(n - 1) {
        for col in 0..(n - 1) {
            let ll = row * n + col;
            faces.push(vec![ll, ll + n, ll + n + 1, ll + 1]);
        }
    }
    problem.pressure = Some(PressureParams::Hydrostatic {
        face_topology: FaceTopology { faces },
        rho_fluid: rho,
        g_magnitude: 1.0,
        z_datum: 1.0,
        up_direction: [0.0, 0.0, 1.0],
        max_iters: 40,
        tolerance: 1e-9,
        relaxation: 1.0,
    });
    let ne = problem.topology.num_edges;
    problem.bounds = Bounds {
        lower: vec![1.0; ne],
        upper: vec![1e6; ne],
    };
    problem.objectives = vec![Box::new(TargetXYZ {
        weight: 1.0,
        node_indices: free.clone(),
        target: Array2::from_shape_fn((free.len(), 3), |(i, d)| {
            let node = free[i];
            match d {
                0 => (node % n) as f64,
                1 => (node / n) as f64,
                _ => -0.1,
            }
        }),
        reduction: TargetGeometryReduction::Sse,
    })];
    problem
}

#[test]
fn load_newton_and_picard_iterations_run_through_the_dispatch() {
    let direct_hydro = {
        let problem = hydro_problem(5, 40.0);
        let q = vec![100.0; problem.topology.num_edges];
        let mut cache = FdmCache::new(&problem).unwrap();
        theseus::fdm::solve_fdm_with_loads(&mut cache, &q, &problem, &anchors(), 1e-12).unwrap();
        cache.x.as_slice().unwrap().to_vec()
    };
    for_each_kind(|kind| {
        // Hydrostatic follower load: Newton + GMRES with the solver as the
        // preconditioner (`precondition`).
        let mut problem = hydro_problem(5, 40.0);
        problem.solver.linear_solver = kind;
        problem.solver.iterative.tolerance = TolerancePolicy::Fixed(1e-12);
        let ne = problem.topology.num_edges;
        let q = vec![100.0; ne];
        let mut cache = FdmCache::new(&problem).unwrap();
        theseus::fdm::solve_fdm_with_loads(&mut cache, &q, &problem, &anchors(), 1e-12).unwrap();
        assert!(cache.x.column(2).iter().all(|z| *z < 0.0), "{kind}");
        assert_kind_close(kind, cache.x.as_slice().unwrap(), &direct_hydro, "hydro x");
        let newton_solves = cache.linear_solver_totals.solves;
        assert!(
            newton_solves > 1,
            "{kind}: preconditioner applications are recorded"
        );
        assert!(cache.linear_solver_totals.converged_all);
        let mut grad = vec![0.0; ne];
        let loss = evaluate(&mut cache, &problem, &q, &mut grad);
        assert!(
            loss.is_finite() && grad.iter().all(|g| g.is_finite()),
            "{kind}"
        );
        assert!(cache.linear_solver_totals.solves > newton_solves, "{kind}");

        // Self-weight Picard iteration + Neumann adjoint.
        let mut problem = problem_for(kind, 8);
        problem.solver.iterative.tolerance = TolerancePolicy::Fixed(1e-12);
        let ne = problem.topology.num_edges;
        problem.self_weight = Some(SelfWeightParams::Prescribed {
            linear_densities: vec![0.2; ne],
            gravity: [0.0, 0.0, -1.0],
            max_iters: 30,
            tolerance: 1e-10,
            relaxation: 1.0,
        });
        let mut cache = FdmCache::new(&problem).unwrap();
        let mut grad = vec![0.0; ne];
        let q = smooth_q(ne, 0.5);
        let loss = evaluate(&mut cache, &problem, &q, &mut grad);
        assert!(loss.is_finite(), "{kind}");
        // Central finite differences on a few edges.
        let mut fd_cache = FdmCache::new(&problem).unwrap();
        let mut scratch = vec![0.0; ne];
        for k in [0, ne / 3, ne - 1] {
            let h = 1e-6;
            let mut plus = q.clone();
            plus[k] += h;
            let mut minus = q.clone();
            minus[k] -= h;
            let fd = (evaluate(&mut fd_cache, &problem, &plus, &mut scratch)
                - evaluate(&mut fd_cache, &problem, &minus, &mut scratch))
                / (2.0 * h);
            let scale = grad[k].abs().max(fd.abs()).max(1e-3);
            assert!(
                (grad[k] - fd).abs() <= 1e-4 * scale,
                "{kind}: self-weight edge {k}: analytic {} vs finite difference {fd}",
                grad[k]
            );
        }
    });
}

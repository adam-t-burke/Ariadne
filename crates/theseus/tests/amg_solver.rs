//! `AmgSolver<CpuBackend>` behind `LinearSystemSolver` (program plan §4.1,
//! §4.4, §4.6): it must solve every regression fixture to the requested
//! tolerance and agree with `DirectSolver`, keep its iteration counts flat
//! with problem size, be a symmetric positive preconditioner, update
//! numerically on small `q` changes and re-setup on large ones, and report
//! the documented error paths.

#[path = "support/fixtures/mod.rs"]
mod fixtures;

use fixtures::grid;
use ndarray::Array2;
use std::sync::atomic::{AtomicBool, Ordering};
use theseus::amg::hierarchy::LevelRef;
use theseus::amg::AmgSolver;
use theseus::backend::CpuBackend;
use theseus::linear_solver::{
    CycleKind, DirectSolver, IterativeSolverOptions, LinearSolver, LinearSolverKind,
    LinearSystemSolver, Precision, SolveRequest, TolerancePolicy,
};
use theseus::types::*;

// ─────────────────────────────────────────────────────────────
//  Helpers
// ─────────────────────────────────────────────────────────────

struct Case {
    name: String,
    problem: Problem,
    q: Vec<f64>,
}

impl Case {
    fn new(name: &str, problem: Problem) -> Self {
        let ne = problem.topology.num_edges;
        let q = fixtures::smooth_q_star(ne);
        Self {
            name: name.to_string(),
            problem,
            q,
        }
    }

    fn n(&self) -> usize {
        self.problem.topology.free_node_indices.len()
    }

    /// `b = p + Σ q x_fixed` from the cache path for `q`.
    fn rhs(&self, q: &[f64]) -> Vec<f64> {
        let mut cache = FdmCache::new(&self.problem).unwrap();
        theseus::fdm::solve_fdm(&mut cache, q, &self.problem, &Array2::zeros((0, 3)), 0.0).unwrap();
        cache.rhs.as_slice().unwrap().to_vec()
    }

    fn direct_solution(&self, q: &[f64], rhs: &[f64], perturbation: f64) -> Vec<f64> {
        let mut direct =
            DirectSolver::from_bounds(&self.problem.topology, &self.problem.bounds).unwrap();
        direct.set_perturbation(perturbation);
        direct.update(q).unwrap();
        let mut x = vec![0.0; rhs.len()];
        direct.solve(request(rhs, 1e-10), &mut x).unwrap();
        x
    }

    fn amg(&self, options: &IterativeSolverOptions) -> AmgSolver<CpuBackend> {
        AmgSolver::cpu(&self.problem.topology, &self.problem.bounds, options).unwrap()
    }
}

fn request(rhs: &[f64], tolerance: f64) -> SolveRequest<'_> {
    SolveRequest {
        rhs,
        x0: None,
        tolerance,
        max_iterations: 200,
        cancel: None,
    }
}

/// The §3 configuration with a small coarsest level so that the test-sized
/// problems (a few hundred to a few thousand nodes) build real multilevel
/// hierarchies; the recommended 2000 would solve most of them directly on
/// level 0.
fn options() -> IterativeSolverOptions {
    IterativeSolverOptions {
        tolerance: TolerancePolicy::Fixed(1e-10),
        coarsest_size: 24,
        ..theseus::amg::recommended_options()
    }
}

fn max_abs(v: &[f64]) -> f64 {
    v.iter().fold(0.0f64, |m, x| m.max(x.abs()))
}

fn max_abs_diff(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b)
        .map(|(x, y)| (x - y).abs())
        .fold(0.0, f64::max)
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn pseudo_random(n: usize, seed: u64) -> Vec<f64> {
    let mut rng = fixtures::rng::Rng::new(seed);
    (0..n).map(|_| rng.next_signed()).collect()
}

/// Every fixture generator of `tests/support/fixtures`, at test sizes.
fn all_fixtures() -> Vec<Case> {
    vec![
        Case::new("grid 24", grid::make_recoverable_grid_problem(24)),
        Case::new(
            "irregular 18",
            fixtures::make_irregular_mesh_problem(18, fixtures::harness::IRREGULAR_SEED),
        ),
        Case::new("dome 10×24", fixtures::make_cable_dome_problem(10, 24)),
        Case::new("few-supports 20", fixtures::make_few_supports_problem(20)),
        Case::new(
            "disconnected 10×3",
            fixtures::make_disconnected_problem(10, 3),
        ),
        Case::new(
            "anisotropic 20 (100×)",
            fixtures::make_anisotropic_q_problem(20, 100.0),
        ),
    ]
}

fn dense_mul(a: &[Vec<f64>], b: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let (n, k, m) = (a.len(), b.len(), b[0].len());
    let mut c = vec![vec![0.0; m]; n];
    for i in 0..n {
        for (l, b_row) in b.iter().enumerate().take(k) {
            let ail = a[i][l];
            if ail != 0.0 {
                for (c_ij, b_lj) in c[i].iter_mut().zip(b_row) {
                    *c_ij += ail * b_lj;
                }
            }
        }
    }
    c
}

fn dense_transpose(a: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let (n, m) = (a.len(), a[0].len());
    (0..m).map(|j| (0..n).map(|i| a[i][j]).collect()).collect()
}

fn max_rel_entry_diff(a: &[Vec<f64>], b: &[Vec<f64>]) -> f64 {
    let scale = a.iter().flatten().fold(0.0f64, |m, v| m.max(v.abs()));
    a.iter()
        .flatten()
        .zip(b.iter().flatten())
        .map(|(x, y)| (x - y).abs() / scale)
        .fold(0.0, f64::max)
}

/// `max_l ‖A_{l+1} − P_lᵀ A_l P_l‖_max / ‖A_{l+1}‖_max` over the hierarchy.
fn worst_galerkin_error(solver: &AmgSolver<CpuBackend>) -> f64 {
    let sizes = solver.level_sizes();
    let mut worst = 0.0f64;
    for l in 0..sizes.len() - 1 {
        let a = if l == 0 {
            LevelRef::Graph(solver.level0()).to_dense()
        } else {
            solver.coarse_matrix(l).to_dense()
        };
        let p = solver.prolongator(l).to_dense();
        let ptap = dense_mul(&dense_transpose(&p), &dense_mul(&a, &p));
        worst = worst.max(max_rel_entry_diff(
            &ptap,
            &solver.coarse_matrix(l + 1).to_dense(),
        ));
    }
    worst
}

// ─────────────────────────────────────────────────────────────
//  Correctness against the direct solver
// ─────────────────────────────────────────────────────────────

#[test]
fn solves_every_fixture_to_tolerance_and_matches_direct() {
    for case in all_fixtures() {
        let n = case.n();
        let rhs = case.rhs(&case.q);
        let x_direct = case.direct_solution(&case.q, &rhs, 0.0);
        let mut solver: Box<dyn LinearSystemSolver> = LinearSolver::new(
            LinearSolverKind::IterativeCpu,
            &case.problem.topology,
            &case.problem.bounds,
            &options(),
        )
        .unwrap();
        assert_eq!(solver.kind(), LinearSolverKind::IterativeCpu);
        solver.update(&case.q).unwrap();
        let mut x = vec![0.0; n * 3];
        let stats = solver.solve(request(&rhs, 1e-10), &mut x).unwrap();
        assert!(stats.converged, "{}: {stats:?}", case.name);
        assert_eq!(stats.backend, LinearSolverKind::IterativeCpu);
        assert!(
            stats.setup_ms > 0.0 && stats.solve_ms > 0.0,
            "{}: {stats:?}",
            case.name
        );
        assert!(
            stats.relative_residual.iter().all(|&r| r <= 1.5e-10),
            "{}: residual {:?}",
            case.name,
            stats.relative_residual
        );
        assert!(
            stats.iterations.iter().all(|&i| (1..60).contains(&i)),
            "{}: iterations {:?}",
            case.name,
            stats.iterations
        );
        let err = max_abs_diff(&x, &x_direct) / max_abs(&x_direct);
        assert!(
            err < 1e-8,
            "{}: |x − x_direct| / |x| = {err:.2e}",
            case.name
        );

        // A second right-hand side (the adjoint) reuses the hierarchy.
        let rhs2 = pseudo_random(n * 3, 11);
        let x2_direct = case.direct_solution(&case.q, &rhs2, 0.0);
        let stats2 = solver.solve(request(&rhs2, 1e-10), &mut x).unwrap();
        assert!(
            stats2.converged && stats2.setup_ms == 0.0,
            "{}: {stats2:?}",
            case.name
        );
        let err2 = max_abs_diff(&x, &x2_direct) / max_abs(&x2_direct);
        assert!(err2 < 1e-8, "{}: random rhs error {err2:.2e}", case.name);

        let mem = solver.memory_bytes();
        assert!(mem.host_bytes > (n * 3 * 8 * 6) as u64 && mem.device_bytes == 0);
    }
}

#[test]
fn iteration_counts_are_flat_across_grid_sizes() {
    let mut per_size = Vec::new();
    for side in [64usize, 128, 224] {
        let case = Case::new(
            &format!("grid {side}"),
            grid::make_recoverable_grid_problem(side),
        );
        let rhs = case.rhs(&case.q);
        let mut solver = case.amg(&options());
        solver.update(&case.q).unwrap();
        let mut x = vec![0.0; case.n() * 3];
        let stats = solver.solve(request(&rhs, 1e-10), &mut x).unwrap();
        assert!(stats.converged);
        eprintln!(
            "grid {side}: levels {:?}, λ_max {:?}, iterations {:?}",
            solver.level_sizes(),
            solver.lambda_max(),
            stats.iterations
        );
        per_size.push(stats.iterations);
    }
    for k in 0..2 {
        let its: Vec<u32> = per_size.iter().map(|it| it[k]).collect();
        let (lo, hi) = (*its.iter().min().unwrap(), *its.iter().max().unwrap());
        assert!(
            hi as f64 <= 1.3 * lo as f64,
            "column {k}: iterations {its:?} spread more than +30%"
        );
    }
}

#[test]
fn coarse_operators_equal_ptap_through_the_solver() {
    let case = Case::new(
        "irregular 14",
        fixtures::make_irregular_mesh_problem(14, fixtures::harness::IRREGULAR_SEED),
    );
    let opts = IterativeSolverOptions {
        coarsest_size: 8,
        ..options()
    };
    let mut solver = case.amg(&opts);
    solver.update(&case.q).unwrap();
    assert!(
        solver.level_sizes().len() >= 3,
        "{:?}",
        solver.level_sizes()
    );
    let worst = worst_galerkin_error(&solver);
    assert!(worst < 1e-12, "coarse operator error {worst:.2e}");
    // Every prolongator is (I − ω D⁻¹ A) P₀ with the aggregate map's P₀.
    for l in 0..solver.level_sizes().len() - 1 {
        let p = solver.prolongator(l);
        let agg = solver.aggregate_of(l);
        assert_eq!(p.n, agg.len());
        assert_eq!(p.ncols, solver.level_sizes()[l + 1]);
        for (u, &agg_u) in agg.iter().enumerate() {
            let (cols, vals) = p.row(u);
            assert!(
                cols.contains(&agg_u),
                "row {u} of P_{l} misses its aggregate"
            );
            // Rows of (I − ω D⁻¹ A) P₀ sum to 1 − ω (D⁻¹ A 1)_u; with zero
            // anchors that is exactly 1 (constants are reproduced).
            let s: f64 = vals.iter().sum();
            assert!(s.is_finite() && s.abs() < 3.0);
        }
    }
}

#[test]
fn v_cycle_is_a_symmetric_positive_preconditioner() {
    let case = Case::new("dome 8×20", fixtures::make_cable_dome_problem(8, 20));
    let n3 = case.n() * 3;
    for degree in [1u8, 2, 3] {
        let opts = IterativeSolverOptions {
            smoother_degree: degree,
            coarsest_size: 20,
            ..options()
        };
        let mut solver = case.amg(&opts);
        assert!(matches!(
            solver.precondition(&vec![0.0; n3], &mut vec![0.0; n3]),
            Err(TheseusError::MissingFactorization)
        ));
        solver.update(&case.q).unwrap();
        let u = pseudo_random(n3, 1);
        let v = pseudo_random(n3, 2);
        let mut mu = vec![0.0; n3];
        let mut mv = vec![0.0; n3];
        solver.precondition(&u, &mut mu).unwrap();
        solver.precondition(&v, &mut mv).unwrap();
        let muv = dot(&mu, &v);
        let umv = dot(&u, &mv);
        let asym = (muv - umv).abs() / muv.abs().max(umv.abs());
        assert!(asym < 1e-10, "degree {degree}: asymmetry {asym:.2e}");
        assert!(dot(&mu, &u) > 0.0, "degree {degree}: ⟨Mu, u⟩ ≤ 0");
        // Deterministic and independent of the previous call.
        let mut mu2 = vec![0.0; n3];
        solver.precondition(&u, &mut mu2).unwrap();
        assert_eq!(mu, mu2);
    }
}

// ─────────────────────────────────────────────────────────────
//  update(q): numeric refill, drift and re-setup
// ─────────────────────────────────────────────────────────────

#[test]
fn numeric_update_matches_fresh_setup_and_keeps_galerkin_exact() {
    let case = Case::new("irregular 20", fixtures::make_irregular_mesh_problem(20, 3));
    let n = case.n();
    let opts = IterativeSolverOptions {
        coarsest_size: 30,
        ..options()
    };
    let mut solver = case.amg(&opts);
    solver.update(&case.q).unwrap();
    assert_eq!(solver.counters().setups, 1);
    let rhs = case.rhs(&case.q);
    let mut x = vec![0.0; n * 3];
    let cold = solver.solve(request(&rhs, 1e-10), &mut x).unwrap();

    // 2 % perturbation: numeric update only, frozen P, coarse operators
    // still exactly Pᵀ A(q') P.
    let noise = pseudo_random(case.q.len(), 5);
    let q2: Vec<f64> = case
        .q
        .iter()
        .zip(&noise)
        .map(|(q, r)| q * (1.0 + 0.02 * r))
        .collect();
    let p_before = solver.prolongator(0).clone();
    solver.update(&q2).unwrap();
    assert_eq!(solver.counters().setups, 1);
    assert_eq!(solver.counters().numeric_updates, 1);
    assert_eq!(solver.counters().lambda_reestimates, 0);
    assert_eq!(solver.prolongator(0), &p_before, "P must stay frozen");
    let worst = worst_galerkin_error(&solver);
    assert!(
        worst < 1e-12,
        "after update: coarse operator error {worst:.2e}"
    );

    let rhs2 = case.rhs(&q2);
    let x_direct = case.direct_solution(&q2, &rhs2, 0.0);
    // Warm start from the previous solution.
    let x0 = x.clone();
    let req = SolveRequest {
        x0: Some(&x0),
        ..request(&rhs2, 1e-10)
    };
    let warm = solver.solve(req, &mut x).unwrap();
    assert!(warm.converged);
    let err = max_abs_diff(&x, &x_direct) / max_abs(&x_direct);
    assert!(err < 1e-8, "updated solve vs direct {err:.2e}");
    assert!(
        warm.max_iterations() <= cold.max_iterations(),
        "warm {warm:?} vs cold {cold:?}"
    );

    // Same answer (to tolerance) as a solver set up from scratch on q2.
    let mut fresh = case.amg(&opts);
    fresh.update(&q2).unwrap();
    let mut xf = vec![0.0; n * 3];
    let fresh_stats = fresh.solve(request(&rhs2, 1e-10), &mut xf).unwrap();
    assert!(fresh_stats.converged);
    let diff = max_abs_diff(&x, &xf) / max_abs(&xf);
    assert!(diff < 1e-8, "updated vs fresh hierarchy {diff:.2e}");
    // The updated hierarchy is not slower than the fresh one by much.
    assert!(
        warm.max_iterations() <= fresh_stats.max_iterations() + 3,
        "warm {warm:?} vs fresh {fresh_stats:?}"
    );
}

#[test]
fn lambda_is_reestimated_after_moderate_drift() {
    let case = Case::new("grid 16", grid::make_recoverable_grid_problem(16));
    let mut solver = case.amg(&options());
    solver.update(&case.q).unwrap();
    let lambda0 = solver.lambda_max().to_vec();
    // 80 % change on one edge: below the re-setup threshold (2×), above the
    // λ threshold (0.5).
    let mut q2 = case.q.clone();
    q2[7] *= 1.8;
    solver.update(&q2).unwrap();
    assert_eq!(solver.counters().setups, 1);
    assert_eq!(solver.counters().lambda_reestimates, 1);
    assert_eq!(solver.lambda_max().len(), lambda0.len());
    // Small follow-up change relative to the new reference: no re-estimate.
    q2[8] *= 1.1;
    solver.update(&q2).unwrap();
    assert_eq!(solver.counters().lambda_reestimates, 1);
}

#[test]
fn resetup_triggers_on_a_3x_q_change() {
    let case = Case::new("grid 20", grid::make_recoverable_grid_problem(20));
    let n = case.n();
    let mut solver = case.amg(&options());
    solver.update(&case.q).unwrap();
    let sizes = solver.level_sizes();
    let q3: Vec<f64> = case.q.iter().map(|q| q * 2.9).collect();
    // Uniform scaling by 2.9: relative change 1.9 ≤ 2, still a numeric update.
    solver.update(&q3).unwrap();
    assert_eq!(solver.counters().setups, 1);
    // One edge at 3.5×: relative change 2.5 > 2 → re-setup.
    let mut q4 = case.q.clone();
    q4[3] *= 3.5;
    solver.update(&q4).unwrap();
    assert_eq!(solver.counters().setups, 2, "{:?}", solver.counters());
    assert_eq!(solver.level_sizes().len(), sizes.len());
    let rhs = case.rhs(&q4);
    let x_direct = case.direct_solution(&q4, &rhs, 0.0);
    let mut x = vec![0.0; n * 3];
    solver.solve(request(&rhs, 1e-10), &mut x).unwrap();
    assert!(max_abs_diff(&x, &x_direct) / max_abs(&x_direct) < 1e-8);
}

#[test]
fn resetup_triggers_after_two_slow_solves() {
    let case = Case::new("grid 20", grid::make_recoverable_grid_problem(20));
    let n = case.n();
    let mut solver = case.amg(&options());
    solver.update(&case.q).unwrap();
    let rhs = case.rhs(&case.q);
    let mut x = vec![0.0; n * 3];
    // Reference iteration count from a loose solve …
    let first = solver.solve(request(&rhs, 1e-1), &mut x).unwrap();
    assert_eq!(solver.setup_iterations(), Some(first.max_iterations()));
    // … then two tight solves need more than twice as many.
    solver.solve(request(&rhs, 1e-11), &mut x).unwrap();
    assert!(!solver.resetup_pending());
    solver.solve(request(&rhs, 1e-11), &mut x).unwrap();
    assert!(solver.resetup_pending());
    solver.update(&case.q).unwrap();
    assert_eq!(solver.counters().setups, 2);
    assert!(!solver.resetup_pending());
    assert_eq!(solver.setup_iterations(), None);
}

#[test]
fn perturbation_shifts_the_diagonal_like_the_direct_solver() {
    let case = Case::new("few-supports 12", fixtures::make_few_supports_problem(12));
    let n = case.n();
    let rhs = case.rhs(&case.q);
    let shift = 1e-3;
    let x_direct = case.direct_solution(&case.q, &rhs, shift);
    let x_unshifted = case.direct_solution(&case.q, &rhs, 0.0);
    assert!(
        max_abs_diff(&x_direct, &x_unshifted) > 1e-6,
        "shift should be visible"
    );
    let mut solver = case.amg(&options());
    solver.set_perturbation(shift);
    assert_eq!(solver.perturbation(), shift);
    solver.update(&case.q).unwrap();
    let mut x = vec![0.0; n * 3];
    solver.solve(request(&rhs, 1e-10), &mut x).unwrap();
    let err = max_abs_diff(&x, &x_direct) / max_abs(&x_direct);
    assert!(err < 1e-8, "shifted solve vs direct {err:.2e}");
}

// ─────────────────────────────────────────────────────────────
//  Options
// ─────────────────────────────────────────────────────────────

#[test]
fn zero_rhs_column_returns_zero_in_no_iterations() {
    let case = Case::new("grid 12", grid::make_recoverable_grid_problem(12));
    let n = case.n();
    let mut rhs = case.rhs(&case.q);
    for v in rhs.iter_mut().skip(1).step_by(3) {
        *v = 0.0;
    }
    let mut solver = case.amg(&options());
    solver.update(&case.q).unwrap();
    let mut x = vec![1.0; n * 3];
    let req = SolveRequest {
        x0: Some(&vec![1.0; n * 3]),
        ..request(&rhs, 1e-10)
    };
    let stats = solver.solve(req, &mut x).unwrap();
    assert_eq!(stats.iterations[1], 0);
    assert_eq!(stats.relative_residual[1], 0.0);
    assert!(x.iter().skip(1).step_by(3).all(|&v| v == 0.0));
    assert!(stats.iterations[0] > 0 && stats.iterations[2] > 0);
    // All-zero rhs: x = 0, no iterations, converged.
    let zero = vec![0.0; n * 3];
    let stats = solver.solve(request(&zero, 1e-10), &mut x).unwrap();
    assert_eq!(stats.iterations, [0; 3]);
    assert!(stats.converged && x.iter().all(|&v| v == 0.0));
}

#[test]
fn f32_preconditioner_with_f64_outer_loop_converges() {
    let case = Case::new("irregular 16", fixtures::make_irregular_mesh_problem(16, 9));
    let n = case.n();
    let rhs = case.rhs(&case.q);
    let x_direct = case.direct_solution(&case.q, &rhs, 0.0);
    let opts = IterativeSolverOptions {
        precondition_precision: Some(Precision::F32),
        ..options()
    };
    let mut solver = case.amg(&opts);
    assert_eq!(solver.precision(), Precision::F32);
    solver.update(&case.q).unwrap();
    let mut x = vec![0.0; n * 3];
    let stats = solver.solve(request(&rhs, 1e-10), &mut x).unwrap();
    assert!(stats.converged, "{stats:?}");
    let err = max_abs_diff(&x, &x_direct) / max_abs(&x_direct);
    assert!(err < 1e-8, "f32 preconditioner: {err:.2e}");
}

#[test]
fn k_cycle_option_solves_with_flexible_cg() {
    let case = Case::new("dome 8×16", fixtures::make_cable_dome_problem(8, 16));
    let n = case.n();
    let rhs = case.rhs(&case.q);
    let x_direct = case.direct_solution(&case.q, &rhs, 0.0);
    let opts = IterativeSolverOptions {
        cycle: CycleKind::K,
        coarsest_size: 20,
        ..options()
    };
    let mut solver = case.amg(&opts);
    solver.update(&case.q).unwrap();
    let mut x = vec![0.0; n * 3];
    let stats = solver.solve(request(&rhs, 1e-10), &mut x).unwrap();
    assert!(stats.converged, "{stats:?}");
    let err = max_abs_diff(&x, &x_direct) / max_abs(&x_direct);
    assert!(err < 1e-8, "K-cycle: {err:.2e}");
}

#[test]
fn single_level_hierarchy_solves_directly() {
    let case = Case::new("grid 10", grid::make_recoverable_grid_problem(10));
    let n = case.n();
    let opts = IterativeSolverOptions {
        coarsest_size: 10_000,
        ..options()
    };
    let mut solver = case.amg(&opts);
    solver.update(&case.q).unwrap();
    assert_eq!(solver.level_sizes(), vec![n]);
    let rhs = case.rhs(&case.q);
    let x_direct = case.direct_solution(&case.q, &rhs, 0.0);
    let mut x = vec![0.0; n * 3];
    let stats = solver.solve(request(&rhs, 1e-10), &mut x).unwrap();
    assert!(stats.iterations.iter().all(|&i| i <= 2), "{stats:?}");
    assert!(max_abs_diff(&x, &x_direct) / max_abs(&x_direct) < 1e-8);
    // Numeric update refreshes the explicit level-0 matrix.
    let q2: Vec<f64> = case.q.iter().map(|q| q * 1.3).collect();
    solver.update(&q2).unwrap();
    let rhs2 = case.rhs(&q2);
    let x_direct2 = case.direct_solution(&q2, &rhs2, 0.0);
    solver.solve(request(&rhs2, 1e-10), &mut x).unwrap();
    assert!(max_abs_diff(&x, &x_direct2) / max_abs(&x_direct2) < 1e-8);
}

// ─────────────────────────────────────────────────────────────
//  Error paths
// ─────────────────────────────────────────────────────────────

#[test]
fn bounds_permitting_nonpositive_q_are_rejected() {
    let mut problem = grid::make_grid_problem(8);
    let ne = problem.topology.num_edges;
    problem.bounds.lower[ne / 2] = 0.0;
    match LinearSolver::new(
        LinearSolverKind::IterativeCpu,
        &problem.topology,
        &problem.bounds,
        &IterativeSolverOptions::default(),
    ) {
        Err(TheseusError::IterativeSolverUnsupported(msg)) => {
            assert!(
                msg.contains("IterativeCpu") && msg.contains("Direct"),
                "{msg}"
            );
            assert!(msg.contains(&format!("edge {}", ne / 2)), "{msg}");
        }
        other => panic!(
            "expected IterativeSolverUnsupported, got {:?}",
            other.map(|_| ())
        ),
    }
    problem.bounds.lower[ne / 2] = -1.0;
    assert!(matches!(
        AmgSolver::cpu(
            &problem.topology,
            &problem.bounds,
            &IterativeSolverOptions::default()
        ),
        Err(TheseusError::IterativeSolverUnsupported(_))
    ));
    // The GPU kind is still unavailable (WS-H).
    match LinearSolver::new(
        LinearSolverKind::IterativeGpu,
        &problem.topology,
        &Bounds::default_for(ne),
        &IterativeSolverOptions::default(),
    ) {
        Err(TheseusError::IterativeSolverUnsupported(msg)) => assert!(msg.contains("GPU"), "{msg}"),
        other => panic!(
            "expected IterativeSolverUnsupported, got {:?}",
            other.map(|_| ())
        ),
    }
}

#[test]
fn cancellation_inside_pcg_is_reported() {
    let case = Case::new("grid 14", grid::make_recoverable_grid_problem(14));
    let n = case.n();
    let rhs = case.rhs(&case.q);
    let mut solver = case.amg(&options());
    solver.update(&case.q).unwrap();
    let mut x = vec![0.0; n * 3];
    let cancel = AtomicBool::new(true);
    let req = SolveRequest {
        cancel: Some(&cancel),
        ..request(&rhs, 1e-10)
    };
    assert!(matches!(
        solver.solve(req, &mut x),
        Err(TheseusError::Cancelled)
    ));
    cancel.store(false, Ordering::Release);
    let stats = solver.solve(req, &mut x).unwrap();
    assert!(stats.converged);
    // Cancel raised by another thread while the solve runs a long budget on
    // a tight tolerance: must come back as Cancelled, not hang.
    let cancel = AtomicBool::new(false);
    let mut long = request(&rhs, 1e-300);
    long.max_iterations = u32::MAX;
    long.cancel = Some(&cancel);
    std::thread::scope(|s| {
        s.spawn(|| {
            std::thread::sleep(std::time::Duration::from_millis(50));
            cancel.store(true, Ordering::Release);
        });
        assert!(matches!(
            solver.solve(long, &mut x),
            Err(TheseusError::Cancelled)
        ));
    });
}

#[test]
fn error_paths_name_the_problem() {
    let case = Case::new("grid 8", grid::make_recoverable_grid_problem(8));
    let n = case.n();
    let ne = case.q.len();
    let rhs = case.rhs(&case.q);
    let mut solver = case.amg(&options());
    let mut x = vec![0.0; n * 3];
    assert!(matches!(
        solver.solve(request(&rhs, 1e-10), &mut x),
        Err(TheseusError::MissingFactorization)
    ));
    assert!(matches!(
        solver.update(&vec![1.0; ne + 1]),
        Err(TheseusError::Shape(_))
    ));
    solver.update(&case.q).unwrap();
    assert!(matches!(
        solver.solve(request(&rhs[..3], 1e-10), &mut x),
        Err(TheseusError::Shape(_))
    ));
    let short = vec![0.0; 3];
    let req = SolveRequest {
        x0: Some(&short),
        ..request(&rhs, 1e-10)
    };
    assert!(matches!(
        solver.solve(req, &mut x),
        Err(TheseusError::Shape(_))
    ));
    // Budget of one iteration: not converged, error carries the residual.
    let mut req = request(&rhs, 1e-12);
    req.max_iterations = 1;
    match solver.solve(req, &mut x) {
        Err(TheseusError::IterativeSolverDidNotConverge {
            iterations,
            relative_residual,
            kind,
        }) => {
            assert_eq!(iterations, 1);
            assert!(relative_residual > 1e-12 && relative_residual.is_finite());
            assert_eq!(kind, LinearSolverKind::IterativeCpu);
        }
        other => panic!(
            "expected IterativeSolverDidNotConverge, got {:?}",
            other.map(|_| ())
        ),
    }
    // The partial iterate is still written.
    assert!(x.iter().any(|&v| v != 0.0));
    // Invalid options.
    let bad = IterativeSolverOptions {
        spectral_alpha: 1.0,
        ..options()
    };
    assert!(matches!(
        AmgSolver::cpu(&case.problem.topology, &case.problem.bounds, &bad),
        Err(TheseusError::Shape(_))
    ));
    let bad = IterativeSolverOptions {
        smoother_degree: 0,
        ..options()
    };
    assert!(matches!(
        AmgSolver::cpu(&case.problem.topology, &case.problem.bounds, &bad),
        Err(TheseusError::Shape(_))
    ));
}

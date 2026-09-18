//! `DirectSolver` behind `LinearSystemSolver` must reproduce the existing
//! `FdmCache` direct path exactly (same assembly, same factorization, same
//! specialised three-column triangular solve), and the `LinearSolver`
//! factory must refuse the iterative kinds until they exist.

#[path = "support/grid.rs"]
#[allow(dead_code)]
mod grid;

use ndarray::Array2;
use std::sync::atomic::{AtomicBool, Ordering};
use theseus::linear_solver::{
    DirectSolver, IterativeSolverOptions, LinearSolver, LinearSolverKind, LinearSystemSolver,
    SolveRequest,
};
use theseus::types::*;

fn smooth_q(ne: usize, offset: f64) -> Vec<f64> {
    (0..ne)
        .map(|k| {
            let t = k as f64 / ne as f64;
            1.0 + 0.5 * (std::f64::consts::TAU * t + offset).sin()
        })
        .collect()
}

fn request(rhs: &[f64]) -> SolveRequest<'_> {
    SolveRequest::new(rhs, &IterativeSolverOptions::default())
}

/// Run the cache path (`solve_fdm`) and the adapter for the same `q` and
/// `perturbation`; the free-node positions must agree bitwise.
fn assert_adapter_matches_cache(problem: &Problem, qs: &[Vec<f64>], perturbation: f64) {
    let anchors = Array2::zeros((0, 3));
    let mut cache = FdmCache::new(problem).unwrap();
    let mut solver = DirectSolver::from_bounds(&problem.topology, &problem.bounds).unwrap();
    solver.set_perturbation(perturbation);
    let n = problem.topology.free_node_indices.len();
    let mut x = vec![0.0; n * 3];

    for q in qs {
        theseus::fdm::solve_fdm(&mut cache, q, problem, &anchors, perturbation).unwrap();
        solver.update(q).unwrap();
        let stats = solver
            .solve(request(cache.rhs.as_slice().unwrap()), &mut x)
            .unwrap();

        assert_eq!(stats.backend, LinearSolverKind::Direct);
        assert!(stats.converged);
        assert_eq!(stats.iterations, [0; 3]);
        assert_eq!(Some(solver.strategy()), cache.strategy());
        assert_eq!(
            x,
            cache.x.as_slice().unwrap(),
            "adapter solution differs from the cache path"
        );

        // The adjoint system reuses the same factorization with another RHS.
        let rhs2: Vec<f64> = (0..n * 3).map(|i| ((i % 7) as f64) - 3.0).collect();
        let mut expected = vec![0.0; n * 3];
        let mut work = vec![0.0; n * 6];
        cache
            .factorization()
            .unwrap()
            .solve_slices::<3>(&rhs2, &mut expected, &mut work);
        solver.solve(request(&rhs2), &mut x).unwrap();
        assert_eq!(x, expected);
    }
}

#[test]
fn direct_solver_matches_cache_path_cholesky() {
    let problem = grid::make_grid_problem(12);
    let ne = problem.topology.num_edges;
    assert_eq!(
        FactorizationStrategy::from_bounds(&problem.bounds),
        FactorizationStrategy::Cholesky
    );
    // First update = fresh symbolic + numeric; later updates = refactor.
    let qs = vec![vec![1.0; ne], smooth_q(ne, 0.0), smooth_q(ne, 1.3)];
    assert_adapter_matches_cache(&problem, &qs, 0.0);
    assert_adapter_matches_cache(&problem, &qs, 1e-12);
}

#[test]
fn direct_solver_matches_cache_path_ldl() {
    let mut problem = grid::make_grid_problem(9);
    let ne = problem.topology.num_edges;
    problem.bounds = Bounds {
        lower: vec![-1.0; ne],
        upper: vec![10.0; ne],
    };
    assert_eq!(
        FactorizationStrategy::from_bounds(&problem.bounds),
        FactorizationStrategy::LDL
    );
    let qs = vec![vec![1.0; ne], smooth_q(ne, 0.4)];
    assert_adapter_matches_cache(&problem, &qs, 0.0);
}

#[test]
fn direct_solver_falls_back_to_ldl_like_the_cache() {
    // Cholesky strategy from the bounds, but a q that makes A indefinite
    // (negative definite here): both paths must switch to LDLᵀ and agree.
    let problem = grid::make_grid_problem(8);
    let ne = problem.topology.num_edges;
    let qs = vec![vec![-1.0; ne], vec![1.0; ne]];
    let anchors = Array2::zeros((0, 3));
    let mut cache = FdmCache::new(&problem).unwrap();
    let mut solver = DirectSolver::from_bounds(&problem.topology, &problem.bounds).unwrap();
    assert_eq!(solver.strategy(), FactorizationStrategy::Cholesky);
    let n = problem.topology.free_node_indices.len();
    let mut x = vec![0.0; n * 3];
    for q in &qs {
        theseus::fdm::solve_fdm(&mut cache, q, &problem, &anchors, 0.0).unwrap();
        solver.update(q).unwrap();
        solver
            .solve(request(cache.rhs.as_slice().unwrap()), &mut x)
            .unwrap();
        assert_eq!(cache.strategy(), Some(FactorizationStrategy::LDL));
        assert_eq!(solver.strategy(), FactorizationStrategy::LDL);
        assert_eq!(x, cache.x.as_slice().unwrap());
    }
}

#[test]
fn direct_solver_error_paths() {
    let problem = grid::make_grid_problem(5);
    let ne = problem.topology.num_edges;
    let n = problem.topology.free_node_indices.len();
    let mut solver = DirectSolver::from_bounds(&problem.topology, &problem.bounds).unwrap();
    let rhs = vec![1.0; n * 3];
    let mut x = vec![0.0; n * 3];

    assert!(matches!(
        solver.solve(request(&rhs), &mut x),
        Err(TheseusError::MissingFactorization)
    ));
    assert!(matches!(
        solver.update(&vec![1.0; ne + 1]),
        Err(TheseusError::Shape(_))
    ));
    solver.update(&vec![1.0; ne]).unwrap();
    assert!(matches!(
        solver.solve(request(&rhs[..3]), &mut x),
        Err(TheseusError::Shape(_))
    ));

    let cancel = AtomicBool::new(true);
    let mut req = request(&rhs);
    req.cancel = Some(&cancel);
    assert!(matches!(
        solver.solve(req, &mut x),
        Err(TheseusError::Cancelled)
    ));
    cancel.store(false, Ordering::Release);
    let stats = solver.solve(req, &mut x).unwrap();
    assert!(stats.solve_ms >= 0.0 && stats.setup_ms >= 0.0);

    let mem = solver.memory_bytes();
    assert!(mem.host_bytes > 0);
    assert_eq!(mem.device_bytes, 0);
    assert_eq!(mem.total_bytes(), mem.host_bytes);
}

#[test]
fn factory_builds_direct_and_refuses_iterative_kinds() {
    let problem = grid::make_grid_problem(6);
    let options = IterativeSolverOptions::default();

    let solver = LinearSolver::new(
        LinearSolverKind::Direct,
        &problem.topology,
        &problem.bounds,
        &options,
    )
    .unwrap();
    assert_eq!(solver.kind(), LinearSolverKind::Direct);

    for kind in [
        LinearSolverKind::IterativeCpu,
        LinearSolverKind::IterativeGpu,
    ] {
        match LinearSolver::new(kind, &problem.topology, &problem.bounds, &options) {
            Err(TheseusError::IterativeSolverUnsupported(msg)) => {
                assert!(msg.contains("not yet available"), "{msg}");
                assert!(msg.contains(&kind.to_string()), "{msg}");
            }
            Err(other) => panic!("unexpected error for {kind}: {other}"),
            Ok(_) => panic!("{kind} should not be available yet"),
        }
    }
}

#[test]
fn boxed_direct_solver_solves_through_the_trait() {
    let problem = grid::make_grid_problem(7);
    let ne = problem.topology.num_edges;
    let n = problem.topology.free_node_indices.len();
    let anchors = Array2::zeros((0, 3));
    let q = smooth_q(ne, 0.7);

    let mut cache = FdmCache::new(&problem).unwrap();
    theseus::fdm::solve_fdm(&mut cache, &q, &problem, &anchors, 0.0).unwrap();

    let mut solver: Box<dyn LinearSystemSolver> = LinearSolver::new(
        problem.solver.linear_solver,
        &problem.topology,
        &problem.bounds,
        &problem.solver.iterative,
    )
    .unwrap();
    solver.update(&q).unwrap();
    let mut x = vec![0.0; n * 3];
    solver
        .solve(request(cache.rhs.as_slice().unwrap()), &mut x)
        .unwrap();
    assert_eq!(x, cache.x.as_slice().unwrap());
}

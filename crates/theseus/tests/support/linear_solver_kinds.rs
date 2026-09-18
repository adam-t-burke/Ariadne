//! Parameterisation of integration tests over [`LinearSolverKind`].
//!
//! `Direct` always runs. `IterativeCpu` runs once `LinearSolver::new`
//! returns a solver for it; while the factory answers
//! `IterativeSolverUnsupported("… not yet available …")` the case is skipped
//! with a printed reason (any other construction error is a test failure).
//! `IterativeGpu` is opt-in through `THESEUS_TEST_GPU=1`, because it needs an
//! adapter the CI machine may not have.

use theseus::linear_solver::{IterativeSolverOptions, LinearSolver, LinearSolverKind};
use theseus::sparse::SparseColMatOwned;
use theseus::types::{Bounds, NetworkTopology, TheseusError};

/// The kinds every parameterised test is run for.
pub fn kinds() -> Vec<LinearSolverKind> {
    let mut kinds = vec![LinearSolverKind::Direct, LinearSolverKind::IterativeCpu];
    if std::env::var("THESEUS_TEST_GPU").is_ok_and(|v| v == "1") {
        kinds.push(LinearSolverKind::IterativeGpu);
    }
    kinds
}

/// Why a kind cannot run in this build, if it cannot: `Some(reason)` when
/// the factory reports the kind as not yet available, `None` when it builds.
/// Panics on any other factory error (that is a real failure, not a skip).
pub fn unavailable_reason(kind: LinearSolverKind) -> Option<String> {
    // The smallest well-formed topology (one edge, one free node) with a
    // positive lower bound, so only availability can fail.
    let incidence = SparseColMatOwned::from_coo(1, 2, &[0, 0], &[0, 1], &[-1.0, 1.0]).unwrap();
    let topology = NetworkTopology {
        free_incidence: incidence.extract_columns(&[0]),
        fixed_incidence: incidence.extract_columns(&[1]),
        incidence,
        num_edges: 1,
        num_nodes: 2,
        free_node_indices: vec![0],
        fixed_node_indices: vec![1],
    };
    let bounds = Bounds {
        lower: vec![1.0],
        upper: vec![10.0],
    };
    match LinearSolver::new(kind, &topology, &bounds, &IterativeSolverOptions::default()) {
        Ok(_) => None,
        Err(TheseusError::IterativeSolverUnsupported(msg)) if msg.contains("not yet available") => {
            Some(msg)
        }
        Err(TheseusError::GpuUnavailable(msg)) if kind == LinearSolverKind::IterativeGpu => {
            Some(msg)
        }
        Err(other) => panic!("LinearSolver::new({kind}) failed: {other}"),
    }
}

/// `true` when `kind` can be constructed in this build (prints the skip
/// reason otherwise).
pub fn available(kind: LinearSolverKind) -> bool {
    match unavailable_reason(kind) {
        None => true,
        Some(reason) => {
            println!("skipping {kind}: {reason}");
            false
        }
    }
}

/// Run `body` for every kind of [`kinds`] that this build provides.
pub fn for_each_kind(mut body: impl FnMut(LinearSolverKind)) {
    for kind in kinds() {
        if available(kind) {
            body(kind);
        }
    }
}

/// Run `body` for every *iterative* kind this build provides.
pub fn for_each_iterative_kind(mut body: impl FnMut(LinearSolverKind)) {
    for_each_kind(|kind| {
        if kind.is_iterative() {
            body(kind);
        }
    });
}

/// `assert_eq!` for `Direct` (its results are bitwise reproducible), a
/// relative/absolute tolerance for the iterative kinds (their solves stop
/// at a residual tolerance).
pub fn assert_kind_close(kind: LinearSolverKind, actual: &[f64], expected: &[f64], what: &str) {
    assert_eq!(actual.len(), expected.len(), "{kind}: {what} length");
    if !kind.is_iterative() {
        assert_eq!(actual, expected, "{kind}: {what} must be bitwise equal");
        return;
    }
    let scale = expected.iter().fold(0.0f64, |m, v| m.max(v.abs())).max(1.0);
    for (i, (a, e)) in actual.iter().zip(expected).enumerate() {
        assert!(
            (a - e).abs() <= 1e-8 * scale,
            "{kind}: {what}[{i}] = {a} differs from {e} (scale {scale})"
        );
    }
}

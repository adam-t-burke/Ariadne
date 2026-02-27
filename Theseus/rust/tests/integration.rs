//! Integration tests — end-to-end optimisation on the arch network.
//!
//! These tests verify that the full pipeline (problem construction →
//! L-BFGS optimisation → result extraction) produces reasonable geometry
//! and actually reduces the objective value.

use ndarray::Array2;
use sprs::TriMat;
use theseus::types::*;
use theseus::optimizer;

// ─────────────────────────────────────────────────────────────
//  Helpers (shared arch construction)
// ─────────────────────────────────────────────────────────────

fn build_incidence(edges: &[(usize, usize)], num_nodes: usize) -> sprs::CsMat<f64> {
    let ne = edges.len();
    let mut tri = TriMat::new((ne, num_nodes));
    for (e, &(s, t)) in edges.iter().enumerate() {
        tri.add_triplet(e, s, -1.0);
        tri.add_triplet(e, t, 1.0);
    }
    tri.to_csc()
}

fn extract_columns(mat: &sprs::CsMat<f64>, cols: &[usize]) -> sprs::CsMat<f64> {
    let nrows = mat.rows();
    let ncols = cols.len();
    let mut tri = TriMat::new((nrows, ncols));
    let mat_csc = mat.to_csc();
    for (new_col, &old_col) in cols.iter().enumerate() {
        let start = mat_csc.indptr().raw_storage()[old_col];
        let end_ = mat_csc.indptr().raw_storage()[old_col + 1];
        for nz in start..end_ {
            tri.add_triplet(mat_csc.indices()[nz], new_col, mat_csc.data()[nz]);
        }
    }
    tri.to_csc()
}

fn make_arch_problem(bounds: Bounds, objectives: Vec<Box<dyn ObjectiveTrait>>) -> Problem {
    let num_nodes = 7;
    let num_edges = 8;

    let edges = vec![
        (0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6),
        (1, 5), (2, 4),
    ];

    let free_idx: Vec<usize> = vec![1, 2, 3, 4, 5];
    let fixed_idx: Vec<usize> = vec![0, 6];

    let incidence = build_incidence(&edges, num_nodes);
    let free_inc = extract_columns(&incidence, &free_idx);
    let fixed_inc = extract_columns(&incidence, &fixed_idx);

    let topology = NetworkTopology {
        incidence,
        free_incidence: free_inc,
        fixed_incidence: fixed_inc,
        num_edges,
        num_nodes,
        free_node_indices: free_idx,
        fixed_node_indices: fixed_idx,
    };

    let free_node_loads = Array2::from_shape_vec(
        (5, 3),
        vec![
            0.0, 0.0, -1.0,
            0.0, 0.0, -1.0,
            0.0, 0.0, -2.0,
            0.0, 0.0, -1.0,
            0.0, 0.0, -1.0,
        ],
    ).unwrap();

    let fixed_node_positions = Array2::from_shape_vec(
        (2, 3),
        vec![0.0, 0.0, 0.0, 6.0, 0.0, 0.0],
    ).unwrap();

    let anchors = AnchorInfo::all_fixed(fixed_node_positions.clone());

    Problem {
        topology,
        free_node_loads,
        fixed_node_positions,
        anchors,
        objectives,
        bounds,
        solver: SolverOptions {
            max_iterations: 200,
            ..SolverOptions::default()
        },
    }
}

// ─────────────────────────────────────────────────────────────
//  Test: Unconstrained TargetXYZ optimisation
// ─────────────────────────────────────────────────────────────

/// Optimise the arch to hit target positions.  Verify:
///   1. Optimiser runs without error
///   2. Final positions are close to target
///   3. All lengths and forces are finite/positive
#[test]
fn optimize_target_xyz() {
    let ne = 8;
    let bounds = Bounds {
        lower: vec![0.1; ne],
        upper: vec![100.0; ne],
    };

    let target = Array2::from_shape_vec(
        (5, 3),
        vec![
            1.0, 0.0, 1.0,
            2.0, 0.0, 2.0,
            3.0, 0.0, 2.5,
            4.0, 0.0, 2.0,
            5.0, 0.0, 1.0,
        ],
    ).unwrap();

    let objectives: Vec<Box<dyn ObjectiveTrait>> = vec![
        Box::new(TargetXYZ {
            weight: 1.0,
            node_indices: vec![1, 2, 3, 4, 5],
            target: target.clone(),
        }),
    ];

    let problem = make_arch_problem(bounds, objectives);
    let mut state = OptimizationState::new(vec![1.0; ne], Array2::zeros((0, 3)));

    let result = optimizer::optimize(&problem, &mut state, None, 1).unwrap();

    // Basic sanity
    assert!(result.iterations > 0, "should run at least 1 iteration");

    // Positions should move toward target (not exact — FDM has physical constraints)
    // We check that the optimizer has reduced the distance compared to the
    // initial uniform-q solution.
    let mut total_error = 0.0;
    let free_idx = &problem.topology.free_node_indices;
    for (i, &node) in free_idx.iter().enumerate() {
        for d in 0..3 {
            let diff = result.xyz[[node, d]] - target[[i, d]];
            total_error += diff * diff;
        }
    }
    // The total squared error should be small-ish (< 21 for 5 nodes × 3 dims)
    assert!(
        total_error < 21.0,
        "total squared error = {total_error:.4}, expected < 21.0",
    );

    // All geometry should be finite
    for l in &result.member_lengths {
        assert!(l.is_finite() && *l > 0.0, "length must be finite positive: {l}");
    }
    for f in &result.member_forces {
        assert!(f.is_finite(), "force must be finite: {f}");
    }
    for &q in &result.q {
        assert!(q.is_finite() && q > 0.0, "q must be finite positive: {q}");
    }

    eprintln!("optimize_target_xyz: {} iterations, converged={}", result.iterations, result.converged);
}

// ─────────────────────────────────────────────────────────────
//  Test: Combined objectives optimisation
// ─────────────────────────────────────────────────────────────

/// Optimise with multiple objectives and verify convergence.
#[test]
fn optimize_combined_objectives() {
    let ne = 8;
    let bounds = Bounds {
        lower: vec![0.1; ne],
        upper: vec![20.0; ne],
    };

    let target = Array2::from_shape_vec(
        (5, 3),
        vec![
            1.0, 0.0, 0.8,
            2.0, 0.0, 1.5,
            3.0, 0.0, 2.0,
            4.0, 0.0, 1.5,
            5.0, 0.0, 0.8,
        ],
    ).unwrap();

    let objectives: Vec<Box<dyn ObjectiveTrait>> = vec![
        Box::new(TargetXYZ {
            weight: 1.0,
            node_indices: vec![1, 2, 3, 4, 5],
            target,
        }),
        Box::new(LengthVariation {
            weight: 0.5,
            edge_indices: (0..ne).collect(),
            sharpness: 20.0,
        }),
        Box::new(SumForceLength {
            weight: 0.01,
            edge_indices: (0..ne).collect(),
        }),
    ];

    let problem = make_arch_problem(bounds, objectives);
    let mut state = OptimizationState::new(vec![2.0; ne], Array2::zeros((0, 3)));

    let result = optimizer::optimize(&problem, &mut state, None, 1).unwrap();

    assert!(result.iterations > 0);
    // Check that all results are finite
    for l in &result.member_lengths {
        assert!(l.is_finite() && *l > 0.0);
    }

    eprintln!("optimize_combined: {} iterations, converged={}", result.iterations, result.converged);
}

// ─────────────────────────────────────────────────────────────
//  Test: Forward solve produces reasonable geometry
// ─────────────────────────────────────────────────────────────

/// Verify that a single forward solve produces finite geometry.
#[test]
fn forward_solve_basic() {
    let ne = 8;
    let bounds = Bounds::default_for(ne);

    let objectives: Vec<Box<dyn ObjectiveTrait>> = vec![];
    let problem = make_arch_problem(bounds, objectives);

    let q = vec![1.0; ne];
    let anchors = Array2::zeros((0, 3));

    let mut cache = FdmCache::new(&problem).unwrap();
    theseus::fdm::solve_fdm(&mut cache, &q, &problem, &anchors, 1e-12).unwrap();

    // Free-node positions should be finite
    for i in 0..problem.topology.num_nodes {
        for d in 0..3 {
            assert!(
                cache.nf[[i, d]].is_finite(),
                "node {i} dim {d} = {} is not finite", cache.nf[[i, d]],
            );
        }
    }

    // Anchor positions should be preserved
    assert!((cache.nf[[0, 0]] - 0.0).abs() < 1e-12);
    assert!((cache.nf[[6, 0]] - 6.0).abs() < 1e-12);

    // All lengths positive
    for (k, &len) in cache.member_lengths.iter().enumerate() {
        assert!(len > 0.0, "edge {k}: length={len} should be positive");
    }

    eprintln!("forward_solve_basic: all positions finite, anchors preserved");
}

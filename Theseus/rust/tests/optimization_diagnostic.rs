//! Diagnostic tests for the optimisation pipeline.
//!
//! Run with:   cargo test --release --test optimization_diagnostic -- --nocapture
//!
//! These tests exercise both Cholesky and LDL paths, vary barrier weights,
//! and print per-iteration loss traces so convergence behaviour is visible.

use ndarray::Array2;
use sprs::TriMat;
use theseus::types::*;

// ─────────────────────────────────────────────────────────────
//  Helpers
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

/// Build an n×n grid network.  4 corner nodes are fixed anchors.
fn make_grid_problem(n: usize, bounds: Bounds, objectives: Vec<Box<dyn ObjectiveTrait>>, solver: SolverOptions) -> Problem {
    let num_nodes = n * n;
    let mut edges = Vec::new();

    for row in 0..n {
        for col in 0..(n - 1) {
            edges.push((row * n + col, row * n + col + 1));
        }
    }
    for row in 0..(n - 1) {
        for col in 0..n {
            edges.push((row * n + col, (row + 1) * n + col));
        }
    }

    let num_edges = edges.len();
    let fixed_idx: Vec<usize> = vec![0, n - 1, n * (n - 1), n * n - 1];
    let free_idx: Vec<usize> = (0..num_nodes).filter(|i| !fixed_idx.contains(i)).collect();
    let nn_free = free_idx.len();

    let incidence = build_incidence(&edges, num_nodes);
    let free_inc = extract_columns(&incidence, &free_idx);
    let fixed_inc = extract_columns(&incidence, &fixed_idx);

    let topology = NetworkTopology {
        incidence,
        free_incidence: free_inc,
        fixed_incidence: fixed_inc,
        num_edges,
        num_nodes,
        free_node_indices: free_idx.clone(),
        fixed_node_indices: fixed_idx.clone(),
    };

    let mut loads_data = vec![0.0; nn_free * 3];
    for i in 0..nn_free {
        loads_data[i * 3 + 2] = -1.0;
    }
    let free_node_loads = Array2::from_shape_vec((nn_free, 3), loads_data).unwrap();

    let fixed_node_positions = Array2::from_shape_vec(
        (4, 3),
        vec![
            0.0, 0.0, 0.0,
            (n - 1) as f64, 0.0, 0.0,
            0.0, (n - 1) as f64, 0.0,
            (n - 1) as f64, (n - 1) as f64, 0.0,
        ],
    ).unwrap();

    let anchors = AnchorInfo::all_fixed(fixed_node_positions.clone());

    Problem {
        topology,
        free_node_loads,
        fixed_node_positions,
        anchors,
        objectives,
        bounds,
        solver,
    }
}

/// Build a TargetXYZ objective for the grid: each free node targets its grid
/// position with a slight z-sag.
fn make_target_xyz(free_idx: &[usize], n: usize, z_sag: f64) -> TargetXYZ {
    let nn_free = free_idx.len();
    let mut target_data = vec![0.0; nn_free * 3];
    for (i, &node) in free_idx.iter().enumerate() {
        let row = node / n;
        let col = node % n;
        target_data[i * 3] = col as f64;
        target_data[i * 3 + 1] = row as f64;
        target_data[i * 3 + 2] = z_sag;
    }
    TargetXYZ {
        weight: 1.0,
        node_indices: free_idx.to_vec(),
        target: Array2::from_shape_vec((nn_free, 3), target_data).unwrap(),
    }
}

fn print_loss_trace(label: &str, result: &SolverResult) {
    eprintln!("\n┌── {label}");
    eprintln!("│  iterations:  {}", result.iterations);
    eprintln!("│  converged:   {}", result.converged);
    eprintln!("│  termination: {}", result.termination_reason);
    eprintln!("│  final loss:  {:.6e}", result.loss_trace.last().copied().unwrap_or(f64::NAN));
    if result.loss_trace.len() > 1 {
        let initial = result.loss_trace[0];
        let final_ = *result.loss_trace.last().unwrap();
        let reduction = if initial > 0.0 { (initial - final_) / initial * 100.0 } else { 0.0 };
        eprintln!("│  loss reduction: {:.2}%  ({:.6e} → {:.6e})", reduction, initial, final_);
    }
    eprintln!("│  loss trace ({} evals):", result.loss_trace.len());
    for (i, &loss) in result.loss_trace.iter().enumerate() {
        if i < 10 || i == result.loss_trace.len() - 1 || i % 10 == 0 {
            eprintln!("│    [{:>4}] {:.8e}", i, loss);
        }
    }
    eprintln!("└──");
}

// ─────────────────────────────────────────────────────────────
//  Test: 10×10 grid with TargetXYZ — Cholesky path
// ─────────────────────────────────────────────────────────────

#[test]
fn diagnostic_grid_cholesky() {
    let n = 10;
    let num_edges = 2 * n * (n - 1);
    let fixed_idx: Vec<usize> = vec![0, n - 1, n * (n - 1), n * n - 1];
    let free_idx: Vec<usize> = (0..n * n).filter(|i| !fixed_idx.contains(i)).collect();

    let bounds = Bounds {
        lower: vec![0.1; num_edges],
        upper: vec![f64::INFINITY; num_edges],
    };

    let objectives: Vec<Box<dyn ObjectiveTrait>> = vec![
        Box::new(make_target_xyz(&free_idx, n, -0.2)),
    ];

    let solver_opts = SolverOptions {
        max_iterations: 200,
        ..SolverOptions::default()
    };

    assert_eq!(FactorizationStrategy::from_bounds(&bounds), FactorizationStrategy::Cholesky);

    let problem = make_grid_problem(n, bounds, objectives, solver_opts);
    let mut state = OptimizationState::new(vec![1.0; num_edges], Array2::zeros((0, 3)));

    let result = theseus::optimizer::optimize(&problem, &mut state, None, 1).unwrap();
    print_loss_trace("10×10 grid, Cholesky", &result);

    assert!(result.iterations > 3, "should run more than 3 iterations, got {}", result.iterations);
    assert!(result.loss_trace.len() >= 2, "should have at least 2 loss evaluations");

    let initial_loss = result.loss_trace[0];
    let min_loss = result.loss_trace.iter().cloned().fold(f64::INFINITY, f64::min);
    assert!(min_loss < initial_loss, "loss should decrease: {initial_loss:.6e} → {min_loss:.6e}");

    for &l in &result.member_lengths {
        assert!(l.is_finite() && l > 0.0, "length must be finite positive: {l}");
    }
    for &q in &result.q {
        assert!(q.is_finite() && q > 0.0, "q must be finite positive: {q}");
    }
}

// ─────────────────────────────────────────────────────────────
//  Test: 10×10 grid with TargetXYZ — LDL path (mixed bounds)
// ─────────────────────────────────────────────────────────────

#[test]
fn diagnostic_grid_ldl() {
    let n = 10;
    let num_edges = 2 * n * (n - 1);
    let fixed_idx: Vec<usize> = vec![0, n - 1, n * (n - 1), n * n - 1];
    let free_idx: Vec<usize> = (0..n * n).filter(|i| !fixed_idx.contains(i)).collect();

    let bounds = Bounds {
        lower: vec![-5.0; num_edges],
        upper: vec![20.0; num_edges],
    };

    let objectives: Vec<Box<dyn ObjectiveTrait>> = vec![
        Box::new(make_target_xyz(&free_idx, n, -0.2)),
    ];

    let solver_opts = SolverOptions {
        max_iterations: 200,
        ..SolverOptions::default()
    };

    assert_eq!(FactorizationStrategy::from_bounds(&bounds), FactorizationStrategy::LDL);

    let problem = make_grid_problem(n, bounds, objectives, solver_opts);
    let mut state = OptimizationState::new(vec![1.0; num_edges], Array2::zeros((0, 3)));

    let result = theseus::optimizer::optimize(&problem, &mut state, None, 1).unwrap();
    print_loss_trace("10×10 grid, LDL (mixed bounds)", &result);

    assert!(result.iterations >= 3, "should run at least 3 iterations, got {}", result.iterations);

    let initial_loss = result.loss_trace[0];
    let min_loss = result.loss_trace.iter().cloned().fold(f64::INFINITY, f64::min);
    assert!(min_loss < initial_loss, "loss should decrease: {initial_loss:.6e} → {min_loss:.6e}");

    for &l in &result.member_lengths {
        assert!(l.is_finite() && l > 0.0, "length must be finite positive: {l}");
    }
}

// ─────────────────────────────────────────────────────────────
//  Test: Cholesky fallback to LDL during optimisation
// ─────────────────────────────────────────────────────────────

/// Start with positive bounds (→ Cholesky) but very tight lower bounds
/// so the optimiser might push q close to zero.  Verify the fallback works.
#[test]
fn diagnostic_cholesky_fallback() {
    let n = 10;
    let num_edges = 2 * n * (n - 1);
    let fixed_idx: Vec<usize> = vec![0, n - 1, n * (n - 1), n * n - 1];
    let free_idx: Vec<usize> = (0..n * n).filter(|i| !fixed_idx.contains(i)).collect();

    let bounds = Bounds {
        lower: vec![1e-6; num_edges],
        upper: vec![f64::INFINITY; num_edges],
    };

    let objectives: Vec<Box<dyn ObjectiveTrait>> = vec![
        Box::new(make_target_xyz(&free_idx, n, -0.2)),
    ];

    let solver_opts = SolverOptions {
        max_iterations: 200,
        ..SolverOptions::default()
    };

    assert_eq!(FactorizationStrategy::from_bounds(&bounds), FactorizationStrategy::Cholesky);

    let problem = make_grid_problem(n, bounds, objectives, solver_opts);
    let mut state = OptimizationState::new(vec![1.0; num_edges], Array2::zeros((0, 3)));

    let result = theseus::optimizer::optimize(&problem, &mut state, None, 1);
    match result {
        Ok(result) => {
            print_loss_trace("Cholesky fallback test (lb=1e-6)", &result);
            assert!(result.iterations > 0, "should complete at least 1 iteration");
            for &l in &result.member_lengths {
                assert!(l.is_finite(), "length must be finite: {l}");
            }
        }
        Err(e) => {
            panic!("Optimisation should not fail with fallback: {e}");
        }
    }
}

// ─────────────────────────────────────────────────────────────
//  Test: Combined objectives (TargetXYZ + LengthVariation)
// ─────────────────────────────────────────────────────────────

#[test]
fn diagnostic_combined_objectives() {
    let n = 10;
    let num_edges = 2 * n * (n - 1);
    let fixed_idx: Vec<usize> = vec![0, n - 1, n * (n - 1), n * n - 1];
    let free_idx: Vec<usize> = (0..n * n).filter(|i| !fixed_idx.contains(i)).collect();

    let bounds = Bounds {
        lower: vec![0.1; num_edges],
        upper: vec![100.0; num_edges],
    };

    let objectives: Vec<Box<dyn ObjectiveTrait>> = vec![
        Box::new(make_target_xyz(&free_idx, n, -0.2)),
        Box::new(LengthVariation {
            weight: 0.1,
            edge_indices: (0..num_edges).collect(),
            sharpness: 10.0,
        }),
        Box::new(SumForceLength {
            weight: 0.001,
            edge_indices: (0..num_edges).collect(),
        }),
    ];

    let solver_opts = SolverOptions {
        max_iterations: 200,
        ..SolverOptions::default()
    };

    let problem = make_grid_problem(n, bounds, objectives, solver_opts);
    let mut state = OptimizationState::new(vec![1.0; num_edges], Array2::zeros((0, 3)));

    let result = theseus::optimizer::optimize(&problem, &mut state, None, 1).unwrap();
    print_loss_trace("10×10 combined (TargetXYZ + LengthVar + SumFL)", &result);

    assert!(result.iterations > 0, "combined should run at least 1 iteration, got {}", result.iterations);

    for &l in &result.member_lengths {
        assert!(l.is_finite() && l > 0.0, "length must be finite positive: {l}");
    }
    for &q in &result.q {
        assert!(q.is_finite() && q > 0.0, "q must be finite positive: {q}");
    }

    if result.loss_trace.len() >= 2 {
        let initial_loss = result.loss_trace[0];
        let min_loss = result.loss_trace.iter().cloned().fold(f64::INFINITY, f64::min);
        assert!(min_loss < initial_loss, "loss should decrease: {initial_loss:.6e} → {min_loss:.6e}");
    }
}

// ─────────────────────────────────────────────────────────────
//  Test: 7-node arch from dummy network geometry
// ─────────────────────────────────────────────────────────────

/// Mimics the arch network used in integration tests: 7 nodes, 8 edges,
/// optimise toward a target geometry within bounds.
#[test]
fn diagnostic_arch_network() {
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
        free_node_indices: free_idx.clone(),
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

    let bounds = Bounds {
        lower: vec![0.1; num_edges],
        upper: vec![100.0; num_edges],
    };

    let objectives: Vec<Box<dyn ObjectiveTrait>> = vec![
        Box::new(TargetXYZ {
            weight: 1.0,
            node_indices: free_idx.clone(),
            target,
        }),
    ];

    let problem = Problem {
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
    };

    let mut state = OptimizationState::new(vec![1.0; num_edges], Array2::zeros((0, 3)));

    let result = theseus::optimizer::optimize(&problem, &mut state, None, 1).unwrap();
    print_loss_trace("7-node arch, TargetXYZ", &result);

    assert!(result.iterations > 3, "arch should run >3 iters, got {}", result.iterations);
    assert!(result.loss_trace.len() >= 2);

    let initial_loss = result.loss_trace[0];
    let min_loss = result.loss_trace.iter().cloned().fold(f64::INFINITY, f64::min);
    assert!(min_loss < initial_loss * 0.5, "should reduce loss by at least 50%");

    eprintln!("\n  Final node positions (arch):");
    for (i, &node) in free_idx.iter().enumerate() {
        eprintln!(
            "    node {:>2}: ({:+.4}, {:+.4}, {:+.4})",
            node,
            result.xyz[[node, 0]],
            result.xyz[[node, 1]],
            result.xyz[[node, 2]],
        );
        let _ = i;
    }

    eprintln!("\n  Final force densities:");
    for (k, &q) in result.q.iter().enumerate() {
        eprintln!("    edge {:>2}: q = {:.6}", k, q);
    }
}

// ─────────────────────────────────────────────────────────────
//  Test: Out-of-bounds initialization
// ─────────────────────────────────────────────────────────────

#[test]
fn diagnostic_out_of_bounds_init() {
    let n = 5;
    let num_edges = 2 * n * (n - 1);
    let fixed_idx: Vec<usize> = vec![0, n - 1, n * (n - 1), n * n - 1];
    let free_idx: Vec<usize> = (0..n * n).filter(|i| !fixed_idx.contains(i)).collect();

    // Bounds: strictly positive q (0.1 to 100.0)
    let bounds = Bounds {
        lower: vec![0.1; num_edges],
        upper: vec![100.0; num_edges],
    };

    let objectives: Vec<Box<dyn ObjectiveTrait>> = vec![
        Box::new(make_target_xyz(&free_idx, n, -0.2)),
    ];

    let solver_opts = SolverOptions {
        max_iterations: 100,
        ..SolverOptions::default()
    };

    let problem = make_grid_problem(n, bounds, objectives, solver_opts);
    
    // Initialize with q = -10.0 (way out of bounds and would cause 
    // Cholesky failure if not clamped or fallen back to LDL).
    let mut state = OptimizationState::new(vec![-10.0; num_edges], Array2::zeros((0, 3)));

    let result = theseus::optimizer::optimize(&problem, &mut state, None, 1).unwrap();
    print_loss_trace("Out-of-bounds initialization (q_init=-10, bounds=[0.1, 100])", &result);

    assert!(result.iterations > 0);
    
    // Verify that the final q values are within the bounds [0.1, 100.0]
    for &q in &result.q {
        assert!(q >= 0.1 - 1e-9 && q <= 100.0 + 1e-9, "q value {} out of bounds", q);
    }
    
    // Verify that the state was updated with feasible values
    for &q in &state.force_densities {
        assert!(q >= 0.1 - 1e-9 && q <= 100.0 + 1e-9, "state q value {} out of bounds", q);
    }
}

//! Release-mode benchmarks for the Theseus solver.
//!
//! Run with:   cargo test --release --test bench_release -- --nocapture
//!
//! These are not criterion benchmarks (to avoid an extra dependency);
//! instead they time key operations using `std::time::Instant` and print
//! the results.

use ndarray::Array2;
use sprs::TriMat;
use std::time::Instant;
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

/// Build a larger grid network for benchmarking.
///
/// `n` × `n` grid of nodes, 4 corner nodes are fixed anchors.
/// Horizontal + vertical edges.
fn make_grid_problem(n: usize) -> Problem {
    let num_nodes = n * n;
    let mut edges = Vec::new();

    // Horizontal edges
    for row in 0..n {
        for col in 0..(n - 1) {
            edges.push((row * n + col, row * n + col + 1));
        }
    }
    // Vertical edges
    for row in 0..(n - 1) {
        for col in 0..n {
            edges.push((row * n + col, (row + 1) * n + col));
        }
    }

    let num_edges = edges.len();

    // Fixed: 4 corners
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
        free_node_indices: free_idx,
        fixed_node_indices: fixed_idx,
    };

    // Loads: uniform downward
    let mut loads_data = vec![0.0; nn_free * 3];
    for i in 0..nn_free {
        loads_data[i * 3 + 2] = -1.0; // z = -1
    }
    let free_node_loads = Array2::from_shape_vec((nn_free, 3), loads_data).unwrap();

    // Fixed positions at grid corners
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

    // Target: nodes at grid positions with slight sag in z
    let target_nodes: Vec<usize> = topology.free_node_indices.clone();
    let mut target_data = vec![0.0; nn_free * 3];
    for (i, &node) in target_nodes.iter().enumerate() {
        let row = node / n;
        let col = node % n;
        target_data[i * 3] = col as f64;
        target_data[i * 3 + 1] = row as f64;
        target_data[i * 3 + 2] = -0.2; // slight sag
    }
    let target = Array2::from_shape_vec((nn_free, 3), target_data).unwrap();

    let objectives: Vec<Box<dyn ObjectiveTrait>> = vec![
        Box::new(TargetXYZ {
            weight: 1.0,
            node_indices: target_nodes,
            target,
        }),
    ];

    let bounds = Bounds {
        lower: vec![0.1; num_edges],
        upper: vec![f64::INFINITY; num_edges],
    };

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
//  Scaling sweep helpers
// ─────────────────────────────────────────────────────────────

/// Grid sizes to test.  n×n grid → 2·n·(n-1) edges.
///   10  →     180 edges
///   32  →   1,984 edges
///  100  →  19,800 edges
///  158  →  49,612 edges   (~50k)
const GRID_SIZES: &[(usize, &str)] = &[
    (10,  "10×10"),
    (32,  "32×32"),
    (100, "100×100"),
    (158, "158×158"),
];

fn fmt_count(n: usize) -> String {
    if n >= 1_000_000 { format!("{:.2}M", n as f64 / 1e6) }
    else if n >= 1_000 { format!("{:.1}k", n as f64 / 1e3) }
    else { format!("{}", n) }
}

fn fmt_time(us: f64) -> String {
    if us >= 1_000_000.0 { format!("{:.2} s",  us / 1e6) }
    else if us >= 1_000.0 { format!("{:.2} ms", us / 1e3) }
    else { format!("{:.1} μs", us) }
}

// ─────────────────────────────────────────────────────────────
//  Benchmarks
// ─────────────────────────────────────────────────────────────

#[test]
fn bench_forward_solve_scaling() {
    eprintln!("\n┌─────────────────────────────────────────────────────────────────┐");
    eprintln!("│                 FORWARD SOLVE  (factor + triangular solve)      │");
    eprintln!("├──────────┬──────────┬───────────┬──────────────────────────── ───┤");
    eprintln!("│  grid    │  edges   │  per-solve│  total (iters)                │");
    eprintln!("├──────────┼──────────┼───────────┼───────────────────────────────┤");

    for &(n, label) in GRID_SIZES {
        let problem = make_grid_problem(n);
        let ne = problem.topology.num_edges;
        let q = vec![1.0; ne];
        let anchors = Array2::zeros((0, 3));

        // Warm-up
        let mut cache = FdmCache::new(&problem).unwrap();
        theseus::fdm::solve_fdm(&mut cache, &q, &problem, &anchors, 1e-12).unwrap();

        // Adaptive iteration count: keep total time ~1-3s per size
        let iters: usize = if ne < 1_000 { 5000 }
            else if ne < 20_000 { 500 }
            else if ne < 200_000 { 50 }
            else if ne < 600_000 { 10 }
            else { 3 };

        let start = Instant::now();
        for _ in 0..iters {
            cache.factorization = None; // force re-factor
            theseus::fdm::solve_fdm(&mut cache, &q, &problem, &anchors, 1e-12).unwrap();
        }
        let elapsed = start.elapsed();
        let per_us = elapsed.as_micros() as f64 / iters as f64;

        eprintln!(
            "│  {:<7} │ {:>8} │ {:>9} │  {:.2} ms  ({} iters){}│",
            label,
            fmt_count(ne),
            fmt_time(per_us),
            elapsed.as_secs_f64() * 1000.0,
            iters,
            " ".repeat(3usize.saturating_sub(format!("{}", iters).len())),
        );
    }
    eprintln!("└──────────┴──────────┴───────────┴───────────────────────────────┘\n");
}

#[test]
fn bench_value_and_gradient_scaling() {
    eprintln!("\n┌─────────────────────────────────────────────────────────────────┐");
    eprintln!("│             VALUE + GRADIENT  (forward + adjoint + explicit)    │");
    eprintln!("├──────────┬──────────┬───────────┬───────────────────────────────┤");
    eprintln!("│  grid    │  edges   │  per-eval │  total (iters)                │");
    eprintln!("├──────────┼──────────┼───────────┼───────────────────────────────┤");

    for &(n, label) in GRID_SIZES {
        let problem = make_grid_problem(n);
        let ne = problem.topology.num_edges;
        let theta = vec![1.0; ne];

        // Warm-up
        let mut cache = FdmCache::new(&problem).unwrap();
        let mut grad = vec![0.0; ne];
        theseus::gradients::value_and_gradient(
            &mut cache, &problem, &theta, &mut grad,
        ).unwrap();

        let iters: usize = if ne < 1_000 { 3000 }
            else if ne < 20_000 { 300 }
            else if ne < 200_000 { 30 }
            else if ne < 600_000 { 5 }
            else { 2 };

        let start = Instant::now();
        for _ in 0..iters {
            cache.factorization = None;
            grad.fill(0.0);
            theseus::gradients::value_and_gradient(
                &mut cache, &problem, &theta, &mut grad,
            ).unwrap();
        }
        let elapsed = start.elapsed();
        let per_us = elapsed.as_micros() as f64 / iters as f64;

        eprintln!(
            "│  {:<7} │ {:>8} │ {:>9} │  {:.2} ms  ({} iters){}│",
            label,
            fmt_count(ne),
            fmt_time(per_us),
            elapsed.as_secs_f64() * 1000.0,
            iters,
            " ".repeat(3usize.saturating_sub(format!("{}", iters).len())),
        );
    }
    eprintln!("└──────────┴──────────┴───────────┴───────────────────────────────┘\n");
}

#[test]
fn bench_full_optimize_scaling() {
    let opt_sizes: &[(usize, &str)] = &[
        (10,  "10×10"),
        (32,  "32×32"),
        (100, "100×100"),
        (158, "158×158"),
    ];

    eprintln!("\n┌─────────────────────────────────────────────────────────────────┐");
    eprintln!("│               FULL L-BFGS OPTIMIZE  (≤200 iters)               │");
    eprintln!("├──────────┬──────────┬───────────┬───────────────────────────────┤");
    eprintln!("│  grid    │  edges   │  per-run  │  total (runs)                 │");
    eprintln!("├──────────┼──────────┼───────────┼───────────────────────────────┤");

    for &(n, label) in opt_sizes {
        let problem = make_grid_problem(n);
        let ne = problem.topology.num_edges;

        let runs: usize = if ne < 1_000 { 20 }
            else if ne < 20_000 { 5 }
            else if ne < 200_000 { 2 }
            else { 1 };

        let start = Instant::now();
        for _ in 0..runs {
            let mut state = OptimizationState::new(vec![1.0; ne], Array2::zeros((0, 3)));
            let result = theseus::optimizer::optimize(&problem, &mut state, None, 1).unwrap();
            let _ = std::hint::black_box(result);
        }
        let elapsed = start.elapsed();
        let per_us = elapsed.as_micros() as f64 / runs as f64;

        eprintln!(
            "│  {:<7} │ {:>8} │ {:>9} │  {:.2} ms  ({} runs){}│",
            label,
            fmt_count(ne),
            fmt_time(per_us),
            elapsed.as_secs_f64() * 1000.0,
            runs,
            " ".repeat(4usize.saturating_sub(format!("{}", runs).len())),
        );
    }
    eprintln!("└──────────┴──────────┴───────────┴───────────────────────────────┘\n");
}

//! Scaling benchmark: full box-constrained L-BFGS-B solves on square grids
//! from ~10k to ~100k edges with a fixed iteration budget.
//!
//! ```text
//! cargo test --release -p theseus --test bench_scale -- --ignored --nocapture
//! THESEUS_SCALE_GRIDS=72,160,224,320 THESEUS_SCALE_ITERS=40 cargo test ...
//! ```
//!
//! Reported per grid: setup (symbolic analysis + first factorization), one
//! warm objective/gradient evaluation, and the full optimisation run with its
//! time per accepted iteration. Peak resident memory is read from
//! `/proc/self/status` where available.

#[path = "support/grid.rs"]
mod grid;

use ndarray::Array2;
use std::sync::atomic::AtomicBool;
use std::time::Instant;
use theseus::types::*;

fn env_list(name: &str, default: &[usize]) -> Vec<usize> {
    std::env::var(name)
        .ok()
        .map(|v| {
            v.split(',')
                .map(|s| s.trim().parse().expect("integer"))
                .collect()
        })
        .unwrap_or_else(|| default.to_vec())
}

fn peak_rss_mb() -> Option<f64> {
    let status = std::fs::read_to_string("/proc/self/status").ok()?;
    let line = status.lines().find(|l| l.starts_with("VmHWM:"))?;
    let kb: f64 = line.split_whitespace().nth(1)?.parse().ok()?;
    Some(kb / 1024.0)
}

#[test]
#[ignore]
fn bench_scale_grid_solves() {
    let grids = env_list("THESEUS_SCALE_GRIDS", &[72, 160, 224]);
    let iters = env_list("THESEUS_SCALE_ITERS", &[40])[0];
    let threads = std::env::var("RAYON_NUM_THREADS").unwrap_or_else(|_| "default".into());

    println!("rayon threads: {threads}; iteration budget: {iters}");
    println!(
        "{:>5} {:>7} {:>7} | {:>9} {:>9} | {:>5} {:>5} {:>9} {:>9} {:>11} {:>12} | {:>8}",
        "grid",
        "edges",
        "free",
        "setup ms",
        "eval ms",
        "iters",
        "evals",
        "total ms",
        "ms/iter",
        "non-eval ms",
        "final loss",
        "rss MB"
    );

    for n in grids {
        let mut problem = grid::make_recoverable_grid_problem(n);
        problem.solver.q_parameterization_mode = QParameterizationMode::DirectBoxBounds;
        problem.solver.absolute_tolerance = 0.0;
        problem.solver.relative_tolerance = 0.0;
        problem.solver.max_iterations = iters;
        let ne = problem.topology.num_edges;
        let nfree = problem.topology.free_node_indices.len();

        // Setup: cache construction, symbolic analysis and the first factorization.
        let anchors = Array2::zeros((0, 3));
        let q = vec![1.0; ne];
        let t = Instant::now();
        let mut cache = FdmCache::new(&problem).unwrap();
        theseus::fdm::solve_fdm(&mut cache, &q, &problem, &anchors, 1e-12).unwrap();
        let setup_ms = t.elapsed().as_secs_f64() * 1e3;

        // One warm fused evaluation (refactor + solve + adjoint).
        let lb = problem.bounds.lower.clone();
        let ub = problem.bounds.upper.clone();
        let lb_idx: Vec<usize> = (0..ne).collect();
        let ub_idx: Vec<usize> = (0..ne).collect();
        let mut grad = vec![0.0; ne];
        let reps = 3;
        let t = Instant::now();
        for _ in 0..reps {
            theseus::gradients::value_and_gradient(
                &mut cache, &problem, &q, &mut grad, &lb, &ub, &lb_idx, &ub_idx,
            )
            .unwrap();
        }
        let eval_ms = t.elapsed().as_secs_f64() * 1e3 / reps as f64;
        drop(cache);

        // Full solve.
        let cancel = AtomicBool::new(false);
        let mut state = OptimizationState::new(vec![1.0; ne], Array2::zeros((0, 3)));
        let t = Instant::now();
        let result = theseus::optimizer::optimize(&problem, &mut state, None, 1, &cancel).unwrap();
        let total_ms = t.elapsed().as_secs_f64() * 1e3;
        let final_loss = result.loss_trace.last().copied().unwrap_or(f64::NAN);
        // loss_trace holds one entry per fused evaluation, so the remainder is
        // setup plus optimizer bookkeeping (L-BFGS-B updates, line search,
        // parameter copies).
        let evals = result.loss_trace.len();
        let non_eval_ms = total_ms - evals as f64 * eval_ms;

        println!(
            "{:>5} {:>7} {:>7} | {:>9.1} {:>9.2} | {:>5} {:>5} {:>9.1} {:>9.2} {:>11.1} {:>12.6e} | {:>8}",
            n,
            ne,
            nfree,
            setup_ms,
            eval_ms,
            result.iterations,
            evals,
            total_ms,
            total_ms / result.iterations.max(1) as f64,
            non_eval_ms,
            final_loss,
            peak_rss_mb().map(|v| format!("{v:.0}")).unwrap_or_else(|| "n/a".into()),
        );
    }
}

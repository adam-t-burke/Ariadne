//! Scaling benchmark: full box-constrained L-BFGS-B solves on square grids
//! (or matched-size fixtures) from ~10k to ~1M edges with a fixed iteration
//! budget.
//!
//! ```text
//! cargo test --release -p theseus --test bench_scale -- --ignored --nocapture
//! THESEUS_SCALE_GRIDS=72,160,224,320 THESEUS_SCALE_ITERS=40 cargo test ...
//! ```
//!
//! Environment switches (§5.1 of `ITERATIVE_SOLVER_PROGRAM.md`):
//!
//! | variable | meaning | default |
//! |---|---|---|
//! | `THESEUS_FIXTURE` | `grid`, `irregular`, `dome`, `few-supports`, `anisotropic` | `grid` |
//! | `THESEUS_SCALE_GRIDS` | grid sides; other fixtures are built at the matching edge count | `72,160,224` |
//! | `THESEUS_SCALE_ITERS` | L-BFGS-B iteration budget (tolerances are zero) | `40` |
//! | `THESEUS_BENCH_REPS` | timed fused evaluations per size (median + IQR reported) | `5` |
//! | `THESEUS_LINEAR_SOLVER` | backend; only `direct` exists, anything else prints "not yet available" and skips | `direct` |
//! | `THESEUS_BENCH_JSON` | path of a JSON-lines file; one object per size is appended | unset |
//! | `THESEUS_MACHINE_ID` | machine id recorded in the JSON (else derived from OS/CPU/RAM) | unset |
//! | `THESEUS_ANISOTROPY` | `q*` ratio of the anisotropic fixture | `100` |
//! | `RAYON_NUM_THREADS` | thread count, recorded as `threads` | rayon default |
//!
//! Reported per size: setup (cache construction, symbolic analysis and first
//! factorization), one warm objective/gradient evaluation (median and IQR of
//! the repetitions, after one discarded warm-up), and the full optimisation
//! run with its time per accepted iteration. Peak resident memory is read from
//! `/proc/self/status` where available.

#[path = "support/fixtures/mod.rs"]
mod fixtures;

use fixtures::harness::{self, FixtureKind, Timing};
use ndarray::Array2;
use serde_json::json;
use std::sync::atomic::AtomicBool;
use std::time::Instant;
use theseus::types::*;

#[test]
#[ignore]
fn bench_scale_grid_solves() {
    let grids = harness::env_list("THESEUS_SCALE_GRIDS", &[72, 160, 224]);
    let iters = harness::env_usize("THESEUS_SCALE_ITERS", 40);
    let reps = harness::env_usize("THESEUS_BENCH_REPS", 5).max(1);
    let fixture_kind = FixtureKind::from_env();
    let backend = harness::linear_solver_from_env();
    let threads = harness::thread_count();
    let git_sha = harness::git_sha();
    let machine_id = harness::machine_id();

    println!(
        "fixture: {}; backend: {backend}; rayon threads: {threads}; iteration budget: {iters}; eval reps: {reps}",
        fixture_kind.name()
    );
    if !harness::is_direct(&backend) {
        println!("linear solver {backend:?}: not yet available; skipping");
        return;
    }
    println!(
        "{:>5} {:>7} {:>7} | {:>9} {:>9} {:>8} | {:>5} {:>5} {:>9} {:>9} {:>11} {:>12} | {:>8}",
        "grid",
        "edges",
        "free",
        "setup ms",
        "eval ms",
        "eval iqr",
        "iters",
        "evals",
        "total ms",
        "ms/iter",
        "non-eval ms",
        "final loss",
        "rss MB"
    );

    for n in grids {
        let t = Instant::now();
        let fixture = harness::build_fixture(fixture_kind, n);
        let build_ms = t.elapsed().as_secs_f64() * 1e3;
        let mut problem = fixture.problem;
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

        // Warm fused evaluations (refactor + solve + adjoint); the first is discarded.
        let lb = problem.bounds.lower.clone();
        let ub = problem.bounds.upper.clone();
        let lb_idx: Vec<usize> = (0..ne).collect();
        let ub_idx: Vec<usize> = (0..ne).collect();
        let mut grad = vec![0.0; ne];
        let mut samples = Vec::with_capacity(reps);
        for rep in 0..=reps {
            let t = Instant::now();
            theseus::gradients::value_and_gradient(
                &mut cache, &problem, &q, &mut grad, &lb, &ub, &lb_idx, &ub_idx,
            )
            .unwrap();
            if rep > 0 {
                samples.push(t.elapsed().as_secs_f64() * 1e3);
            }
        }
        let eval = Timing::from_samples(samples);
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
        let non_eval_ms = total_ms - evals as f64 * eval.median;
        let peak_rss = harness::peak_rss_bytes();

        println!(
            "{:>5} {:>7} {:>7} | {:>9.1} {:>9.2} {:>8.2} | {:>5} {:>5} {:>9.1} {:>9.2} {:>11.1} {:>12.6e} | {:>8}",
            n,
            ne,
            nfree,
            setup_ms,
            eval.median,
            eval.iqr(),
            result.iterations,
            evals,
            total_ms,
            total_ms / result.iterations.max(1) as f64,
            non_eval_ms,
            final_loss,
            peak_rss
                .map(|b| format!("{:.0}", b as f64 / (1024.0 * 1024.0)))
                .unwrap_or_else(|| "n/a".into()),
        );

        let mut parameters = fixture.parameters.clone();
        parameters.insert("iterations_budget".into(), json!(iters));
        parameters.insert("eval_reps".into(), json!(reps));
        parameters.insert("tolerances".into(), json!(0.0));
        parameters.insert("q_parameterization".into(), json!("direct-box-bounds"));
        parameters.insert("bounds".into(), json!([lb[0], ub[0]]));
        parameters.insert("q_start".into(), json!(1.0));
        let record = json!({
            "schema": 1,
            "harness": "bench_scale",
            "timestamp_utc": harness::timestamp_utc(),
            "fixture": fixture_kind.name(),
            "grid_side": n,
            "edges": ne,
            "free_nodes": nfree,
            "nodes": problem.topology.num_nodes,
            "fixed_nodes": problem.topology.fixed_node_indices.len(),
            "backend": backend,
            "adapter": serde_json::Value::Null,
            "threads": threads,
            "parameters": parameters,
            "build_ms": build_ms,
            "setup_ms": setup_ms,
            "eval_ms": eval.to_json(),
            "evaluations": evals,
            "iterations": result.iterations,
            "linear_solver_iterations": serde_json::Value::Null,
            "total_ms": total_ms,
            "ms_per_iteration": total_ms / result.iterations.max(1) as f64,
            "non_eval_ms": non_eval_ms,
            "peak_rss_bytes": peak_rss,
            "device_bytes": serde_json::Value::Null,
            "final_loss": final_loss,
            "termination": result.termination_reason,
            "git_sha": git_sha,
            "machine_id": machine_id,
        });
        harness::append_json_line(&record);
    }
}

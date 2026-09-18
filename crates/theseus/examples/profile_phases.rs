//! Per-phase wall-clock profile of one fused objective/gradient evaluation on
//! square grid networks (or matched-size fixtures), from a few thousand to
//! ~1M edges.
//!
//!   cargo run --release -p theseus --example profile_phases -- [grid sizes...]
//!
//! Default grid sizes: 72 (10k edges), 160 (51k), 224 (100k). Set
//! `RAYON_NUM_THREADS` to compare threading configurations. Environment
//! switches shared with `tests/bench_scale.rs`: `THESEUS_FIXTURE`
//! (`grid` | `irregular` | `dome` | `few-supports` | `anisotropic`; other
//! fixtures are built at the grid's edge count), `THESEUS_LINEAR_SOLVER`
//! (only `direct` is available; anything else prints "not yet available"),
//! `THESEUS_BENCH_REPS` (repetitions per size; default 5, 3 above 60k edges),
//! `THESEUS_BENCH_JSON` (append one JSON line per size with the per-phase
//! medians) and `THESEUS_MACHINE_ID`.

#[path = "../tests/support/fixtures/mod.rs"]
mod fixtures;

use fixtures::harness::{self, FixtureKind, Timing};
use ndarray::Array2;
use serde_json::json;
use std::time::Instant;
use theseus::types::*;

pub fn ms(start: Instant) -> f64 {
    start.elapsed().as_secs_f64() * 1e3
}

const PHASES: [&str; 9] = [
    "asm_A", "asm_rhs", "factor", "solve", "geom", "loss", "expl", "adjoint", "impl",
];

fn main() {
    let sizes: Vec<usize> = std::env::args()
        .skip(1)
        .map(|a| a.parse().expect("grid size"))
        .collect();
    let sizes = if sizes.is_empty() {
        vec![72, 160, 224]
    } else {
        sizes
    };
    let fixture_kind = FixtureKind::from_env();
    let backend = harness::linear_solver_from_env();
    let threads = harness::thread_count();
    let reps_override = std::env::var("THESEUS_BENCH_REPS")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .filter(|&r| r > 0);
    println!(
        "fixture: {}; backend: {backend}; rayon threads: {threads}",
        fixture_kind.name()
    );
    if !harness::is_direct(&backend) {
        println!("linear solver {backend:?}: not yet available; skipping");
        return;
    }
    let git_sha = harness::git_sha();
    let machine_id = harness::machine_id();
    println!(
        "{:>5} {:>7} {:>9} | {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} | {:>9}",
        "grid",
        "edges",
        "nodes",
        PHASES[0],
        PHASES[1],
        PHASES[2],
        PHASES[3],
        PHASES[4],
        PHASES[5],
        PHASES[6],
        PHASES[7],
        PHASES[8],
        "total"
    );
    for &n in &sizes {
        // The grid keeps the original flat-target fixture; the others use the
        // recoverable fixtures at the matching edge count.
        let (problem, parameters) = if fixture_kind == FixtureKind::Grid {
            let mut p = serde_json::Map::new();
            p.insert("grid_side".into(), json!(n));
            (fixtures::grid::make_grid_problem(n), p)
        } else {
            let f = harness::build_fixture(fixture_kind, n);
            (f.problem, f.parameters)
        };
        let ne = problem.topology.num_edges;
        let nn = problem.topology.num_nodes;
        let anchors = Array2::zeros((0, 3));
        let q = vec![1.0; ne];
        let mut cache = FdmCache::new(&problem).unwrap();
        // Warm up: symbolic analysis and first numeric factorization.
        theseus::fdm::solve_fdm(&mut cache, &q, &problem, &anchors, 1e-12).unwrap();

        let reps = reps_override.unwrap_or(if ne > 60_000 { 3 } else { 5 });
        let mut phase_samples: Vec<Vec<f64>> = vec![Vec::with_capacity(reps); PHASES.len()];
        let mut total_samples = Vec::with_capacity(reps);
        for _ in 0..reps {
            let mut acc = [0.0f64; 9];
            let t_all = Instant::now();
            cache.q.copy_from_slice(&q);
            let t = Instant::now();
            theseus::fdm::assemble_a(&mut cache);
            acc[0] += ms(t);
            theseus::fdm::update_fixed_positions(&mut cache, &problem, &anchors);
            let t = Instant::now();
            theseus::fdm::assemble_rhs(&mut cache, &problem);
            acc[1] += ms(t);
            // factor + solve are fused in factor_and_solve; time the refactor separately.
            let t = Instant::now();
            cache.direct_solver_mut().unwrap().factor().unwrap();
            acc[2] += ms(t);
            let t = Instant::now();
            theseus::fdm::solve_forward_system(&mut cache, None).unwrap();
            acc[3] += ms(t);
            for (i, &node) in problem.topology.free_node_indices.iter().enumerate() {
                for d in 0..3 {
                    cache.nf[[node, d]] = cache.x[[i, d]];
                }
            }
            let t = Instant::now();
            theseus::fdm::compute_geometry(&mut cache, &problem);
            acc[4] += ms(t);
            let t = Instant::now();
            let snap = GeometrySnapshot {
                xyz_full: &cache.nf,
                member_lengths: &cache.member_lengths,
                member_forces: &cache.member_forces,
                reactions: &cache.reactions,
            };
            let loss = theseus::objectives::total_loss(&problem.objectives, &snap);
            acc[5] += ms(t);
            std::hint::black_box(loss);
            let t = Instant::now();
            cache.grad_q.fill(0.0);
            cache.grad_nf.fill(0.0);
            theseus::gradients::accumulate_explicit_gradients(&mut cache, &problem);
            acc[6] += ms(t);
            let t = Instant::now();
            theseus::gradients::solve_adjoint(&mut cache).unwrap();
            acc[7] += ms(t);
            let t = Instant::now();
            theseus::gradients::accumulate_implicit_gradients(&mut cache, &problem);
            acc[8] += ms(t);
            total_samples.push(ms(t_all));
            for (samples, v) in phase_samples.iter_mut().zip(acc) {
                samples.push(v);
            }
        }
        let phases: Vec<Timing> = phase_samples
            .into_iter()
            .map(Timing::from_samples)
            .collect();
        let total = Timing::from_samples(total_samples);
        println!(
            "{:>5} {:>7} {:>9} | {:>8.2} {:>8.2} {:>8.2} {:>8.2} {:>8.2} {:>8.2} {:>8.2} {:>8.2} {:>8.2} | {:>9.2}",
            n,
            ne,
            nn,
            phases[0].median,
            phases[1].median,
            phases[2].median,
            phases[3].median,
            phases[4].median,
            phases[5].median,
            phases[6].median,
            phases[7].median,
            phases[8].median,
            total.median
        );

        let mut phase_json = serde_json::Map::new();
        for (name, timing) in PHASES.iter().zip(&phases) {
            phase_json.insert((*name).into(), timing.to_json());
        }
        let record = json!({
            "schema": 1,
            "harness": "profile_phases",
            "timestamp_utc": harness::timestamp_utc(),
            "fixture": fixture_kind.name(),
            "grid_side": n,
            "edges": ne,
            "free_nodes": problem.topology.free_node_indices.len(),
            "nodes": nn,
            "backend": backend,
            "adapter": serde_json::Value::Null,
            "threads": threads,
            "parameters": parameters,
            "reps": reps,
            "phases_ms": phase_json,
            "eval_ms": total.to_json(),
            "peak_rss_bytes": harness::peak_rss_bytes(),
            "git_sha": git_sha,
            "machine_id": machine_id,
        });
        harness::append_json_line(&record);
    }
}

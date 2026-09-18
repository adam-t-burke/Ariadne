//! WS-C measurements (program plan §7.3): `AmgSolver<CpuBackend>` setup,
//! numeric `update(q)`, cold and 2 %-warm solves against `DirectSolver`
//! update + solve on the same systems. Ignored by default; run with
//!
//! ```sh
//! env RAYON_NUM_THREADS=1 THESEUS_AMG_SIZES=224,448,708 THESEUS_AMG_FIXTURES=grid,irregular,dome \
//!   cargo test -p theseus --release --test amg_bench -- --ignored --nocapture
//! ```
//!
//! The configuration is [`theseus::amg::recommended_options`] (§3: V-cycle,
//! three passes, degree 2, α 10, coarsest 2000). Knobs: `THESEUS_AMG_SIZES`
//! (grid sides; other fixtures are matched by edge count),
//! `THESEUS_AMG_FIXTURES`, `THESEUS_AMG_COARSEST` (coarsest-level size,
//! default 2000), `THESEUS_AMG_REPS` (median of this many, default 3),
//! `THESEUS_AMG_TOL` (solve tolerance, default 1e-8 as in the Phase-0
//! tables), `THESEUS_AMG_DEGREE`, `THESEUS_AMG_PASSES`, `THESEUS_AMG_ALPHA`.
//! Prints one markdown row per fixture and the load average.

#[path = "support/fixtures/mod.rs"]
mod fixtures;

use fixtures::harness::{build_fixture, env_list, env_usize, thread_count, FixtureKind};
use ndarray::Array2;
use std::time::Instant;
use theseus::amg::AmgSolver;
use theseus::backend::CpuBackend;
use theseus::linear_solver::{
    DirectSolver, IterativeSolverOptions, LinearSystemSolver, SolveRequest, TolerancePolicy,
};
use theseus::types::*;

fn env_f64(name: &str, default: f64) -> f64 {
    std::env::var(name)
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(default)
}

fn median(mut v: Vec<f64>) -> f64 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v[v.len() / 2]
}

fn load_average() -> String {
    std::fs::read_to_string("/proc/loadavg")
        .map(|s| s.split_whitespace().take(3).collect::<Vec<_>>().join(" "))
        .unwrap_or_else(|_| "n/a".into())
}

fn rhs_for(problem: &Problem, q: &[f64]) -> Vec<f64> {
    let mut cache = FdmCache::new(problem).unwrap();
    theseus::fdm::solve_fdm(&mut cache, q, problem, &Array2::zeros((0, 3)), 0.0).unwrap();
    cache.rhs.as_slice().unwrap().to_vec()
}

fn perturbed(q: &[f64]) -> Vec<f64> {
    let mut rng = fixtures::rng::Rng::new(0xC0FFEE);
    q.iter()
        .map(|v| v * (1.0 + 0.02 * rng.next_signed()))
        .collect()
}

fn request<'a>(rhs: &'a [f64], x0: Option<&'a [f64]>, tol: f64) -> SolveRequest<'a> {
    SolveRequest {
        rhs,
        x0,
        tolerance: tol,
        max_iterations: 500,
        cancel: None,
    }
}

fn ms(start: Instant) -> f64 {
    start.elapsed().as_secs_f64() * 1e3
}

fn its(s: &[u32; 3]) -> String {
    format!("{}/{}/{}", s[0], s[1], s[2])
}

#[test]
#[ignore = "benchmark: run with --ignored --nocapture"]
fn amg_vs_direct() {
    let sizes = env_list("THESEUS_AMG_SIZES", &[224]);
    let kinds: Vec<FixtureKind> = std::env::var("THESEUS_AMG_FIXTURES")
        .ok()
        .map(|v| {
            v.split(',')
                .filter(|s| !s.trim().is_empty())
                .map(|s| FixtureKind::parse(s).unwrap_or_else(|| panic!("unknown fixture {s}")))
                .collect()
        })
        .unwrap_or_else(|| vec![FixtureKind::Grid]);
    let coarsest = env_usize("THESEUS_AMG_COARSEST", 2000);
    let reps = env_usize("THESEUS_AMG_REPS", 3).max(1);
    let tol = env_f64("THESEUS_AMG_TOL", 1e-8);
    let degree = env_usize("THESEUS_AMG_DEGREE", 2) as u8;
    let passes = env_usize("THESEUS_AMG_PASSES", 3) as u8;
    let alpha = env_f64("THESEUS_AMG_ALPHA", 10.0);
    let threads = thread_count();

    println!(
        "\nWS-C AMG vs direct — threads {threads}, coarsest {coarsest}, degree {degree}, passes {passes}, α {alpha}, tol {tol:e}, median of {reps}, load {}",
        load_average()
    );
    println!("| fixture | edges | free nodes | levels | direct update+solve ms | AMG setup ms | AMG update(q) ms | cold it x/y/z | cold ms | warm it x/y/z | warm ms | AMG host MiB |");
    println!("|---|---:|---:|---|---:|---:|---:|---|---:|---|---:|---:|");

    for &size in &sizes {
        for &kind in &kinds {
            let fx = build_fixture(kind, size);
            let problem = &fx.problem;
            let ne = problem.topology.num_edges;
            let n = problem.topology.free_node_indices.len();
            let q = fixtures::smooth_q_star(ne);
            let q2 = perturbed(&q);
            let rhs = rhs_for(problem, &q);
            let rhs2 = rhs_for(problem, &q2);

            // Direct: numeric refactorization + one 3-rhs solve.
            let mut direct = DirectSolver::from_bounds(&problem.topology, &problem.bounds).unwrap();
            direct.update(&q).unwrap();
            let mut x_direct = vec![0.0; n * 3];
            direct
                .solve(request(&rhs, None, tol), &mut x_direct)
                .unwrap();
            let mut direct_ms = Vec::new();
            for r in 0..reps {
                let qq = if r % 2 == 0 { &q } else { &q2 };
                let t = Instant::now();
                direct.update(qq).unwrap();
                direct
                    .solve(request(&rhs, None, tol), &mut x_direct)
                    .unwrap();
                direct_ms.push(ms(t));
            }

            let options = IterativeSolverOptions {
                tolerance: TolerancePolicy::Fixed(tol),
                coarsest_size: coarsest as u32,
                smoother_degree: degree,
                aggregation_passes: passes,
                spectral_alpha: alpha,
                ..theseus::amg::recommended_options()
            };
            let mut amg: AmgSolver<CpuBackend> =
                AmgSolver::cpu(&problem.topology, &problem.bounds, &options).unwrap();

            // Setup (median over fresh solvers would rebuild topology too;
            // measure the first update, which is the setup).
            let mut setup_ms = Vec::new();
            for _ in 0..reps {
                let mut fresh: AmgSolver<CpuBackend> =
                    AmgSolver::cpu(&problem.topology, &problem.bounds, &options).unwrap();
                let t = Instant::now();
                fresh.update(&q).unwrap();
                setup_ms.push(ms(t));
            }
            amg.update(&q).unwrap();
            let levels = amg.level_sizes();

            // Cold solve.
            let mut x = vec![0.0; n * 3];
            let mut cold_ms = Vec::new();
            let mut cold_stats = None;
            for _ in 0..reps {
                x.fill(0.0);
                let t = Instant::now();
                let s = amg.solve(request(&rhs, None, tol), &mut x).unwrap();
                cold_ms.push(ms(t));
                cold_stats = Some(s);
            }
            let cold_stats = cold_stats.unwrap();
            let x_cold = x.clone();
            let err = x_cold
                .iter()
                .zip(&x_direct)
                .map(|(a, b)| (a - b).abs())
                .fold(0.0, f64::max)
                / x_direct.iter().fold(0.0f64, |m, v| m.max(v.abs()));

            // Numeric update on the 2 % perturbed q (alternating so every
            // rep is a real change), then the warm solve from x_cold.
            let mut update_ms = Vec::new();
            let mut warm_ms = Vec::new();
            let mut warm_stats = None;
            let setups_before = amg.counters().setups;
            for r in 0..reps {
                let (qq, bb) = if r % 2 == 0 { (&q2, &rhs2) } else { (&q, &rhs) };
                let t = Instant::now();
                amg.update(qq).unwrap();
                update_ms.push(ms(t));
                x.copy_from_slice(&x_cold);
                let t = Instant::now();
                let s = amg.solve(request(bb, Some(&x_cold), tol), &mut x).unwrap();
                warm_ms.push(ms(t));
                warm_stats = Some(s);
            }
            assert_eq!(
                amg.counters().setups,
                setups_before,
                "update(q) must stay numeric"
            );
            let warm_stats = warm_stats.unwrap();
            let mem = amg.memory_bytes();

            println!(
                "| {} {} | {} | {} | {:?} | {:.1} | {:.1} | {:.1} | {} | {:.1} | {} | {:.1} | {:.0} |",
                kind.name(),
                size,
                ne,
                n,
                levels,
                median(direct_ms),
                median(setup_ms),
                median(update_ms),
                its(&cold_stats.iterations),
                median(cold_ms),
                its(&warm_stats.iterations),
                median(warm_ms),
                mem.host_bytes as f64 / (1u64 << 20) as f64,
            );
            eprintln!(
                "  {} {}: |x_amg − x_direct|/|x_direct| = {err:.2e}, λ_max {:?}, load {}",
                kind.name(),
                size,
                amg.lambda_max(),
                load_average()
            );
        }
    }
}

#[test]
#[ignore = "profile: run with --ignored --nocapture"]
fn update_breakdown() {
    use theseus::amg::hierarchy::{galerkin_numeric, LevelRef, ScratchPool};
    let sizes = env_list("THESEUS_AMG_SIZES", &[224]);
    let kinds: Vec<FixtureKind> = std::env::var("THESEUS_AMG_FIXTURES")
        .ok()
        .map(|v| {
            v.split(',')
                .map(|s| FixtureKind::parse(s).unwrap())
                .collect()
        })
        .unwrap_or_else(|| vec![FixtureKind::Grid]);
    for &size in &sizes {
        for &kind in &kinds {
            let fx = build_fixture(kind, size);
            let problem = &fx.problem;
            let q = fixtures::smooth_q_star(problem.topology.num_edges);
            let q2 = perturbed(&q);
            let mut amg: AmgSolver<CpuBackend> = AmgSolver::cpu(
                &problem.topology,
                &problem.bounds,
                &theseus::amg::recommended_options(),
            )
            .unwrap();
            amg.update(&q).unwrap();
            let t = Instant::now();
            amg.update(&q2).unwrap();
            let total = ms(t);
            let scratch = ScratchPool::new();
            let levels = amg.level_sizes().len();
            let mut parts = Vec::new();
            for l in 0..levels - 1 {
                let a = if l == 0 {
                    LevelRef::Graph(amg.level0())
                } else {
                    LevelRef::Csr(amg.coarse_matrix(l))
                };
                let p = amg.prolongator(l);
                let pt = p.transpose();
                let (mut ap, mut out) = theseus::amg::hierarchy::galerkin_pattern(&pt, a, p);
                galerkin_numeric(&pt, a, p, &mut ap, &mut out, &scratch);
                let t = Instant::now();
                for _ in 0..5 {
                    galerkin_numeric(&pt, a, p, &mut ap, &mut out, &scratch);
                }
                parts.push(format!(
                    "L{l}→{}: {:.2} ms (n {} nnz(P) {} nnz(AP) {} nnz(A_c) {} avg row {:.1})",
                    l + 1,
                    ms(t) / 5.0,
                    a.n(),
                    p.nnz(),
                    ap.nnz(),
                    out.nnz(),
                    out.nnz() as f64 / out.n as f64
                ));
            }
            println!(
                "{} {}: update {total:.1} ms; {}",
                kind.name(),
                size,
                parts.join("; ")
            );
        }
    }
}

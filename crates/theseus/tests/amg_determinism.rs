//! `AmgSolver<CpuBackend>` must produce bitwise identical hierarchies and
//! solutions regardless of the rayon thread count (program plan §4.1): the
//! aggregation is sequential, the Galerkin refill is row-parallel with a
//! fixed accumulation order per row, and every reduction goes through the
//! fixed-order chunked sums of `CpuBackend`.

#[path = "support/fixtures/mod.rs"]
mod fixtures;

use fixtures::grid;
use ndarray::Array2;
use theseus::amg::AmgSolver;
use theseus::backend::CpuBackend;
use theseus::linear_solver::{
    IterativeSolverOptions, LinearSystemSolver, SolveRequest, TolerancePolicy,
};
use theseus::types::*;

/// 160 × 160 grid (50,880 edges, 25,596 free nodes): every level-0 loop is
/// above the parallel threshold and the hierarchy has several levels.
const SIDE: usize = 160;

struct Run {
    lambda_max: Vec<f64>,
    coarse_values: Vec<Vec<f64>>,
    x_cold: Vec<f64>,
    x_warm: Vec<f64>,
    iterations: [[u32; 3]; 2],
}

fn fixture() -> (Problem, Vec<f64>, Vec<f64>) {
    let problem = grid::make_recoverable_grid_problem(SIDE);
    let q = fixtures::smooth_q_star(problem.topology.num_edges);
    let mut cache = FdmCache::new(&problem).unwrap();
    theseus::fdm::solve_fdm(&mut cache, &q, &problem, &Array2::zeros((0, 3)), 0.0).unwrap();
    let rhs = cache.rhs.as_slice().unwrap().to_vec();
    (problem, q, rhs)
}

fn run(problem: &Problem, q: &[f64], rhs: &[f64]) -> Run {
    let options = IterativeSolverOptions {
        tolerance: TolerancePolicy::Fixed(1e-10),
        coarsest_size: 200,
        ..theseus::amg::recommended_options()
    };
    let mut solver: AmgSolver<CpuBackend> =
        AmgSolver::cpu(&problem.topology, &problem.bounds, &options).unwrap();
    solver.update(q).unwrap();
    let n3 = rhs.len();
    let req = SolveRequest {
        rhs,
        x0: None,
        tolerance: 1e-10,
        max_iterations: 500,
        cancel: None,
    };
    let mut x_cold = vec![0.0; n3];
    let cold = solver.solve(req, &mut x_cold).unwrap();
    // 2 % perturbation → numeric update on the frozen pattern, warm start.
    let q2: Vec<f64> = q
        .iter()
        .enumerate()
        .map(|(e, v)| v * (1.0 + 0.02 * ((e % 7) as f64 / 7.0 - 0.5)))
        .collect();
    solver.update(&q2).unwrap();
    assert_eq!(solver.counters().setups, 1);
    let mut x_warm = x_cold.clone();
    let warm = solver
        .solve(
            SolveRequest {
                x0: Some(&x_cold),
                ..req
            },
            &mut x_warm,
        )
        .unwrap();
    let levels = solver.level_sizes().len();
    Run {
        lambda_max: solver.lambda_max().to_vec(),
        coarse_values: (1..levels)
            .map(|l| solver.coarse_matrix(l).values.clone())
            .collect(),
        x_cold,
        x_warm,
        iterations: [cold.iterations, warm.iterations],
    }
}

fn assert_bitwise(actual: &[f64], expected: &[f64], what: &str) {
    assert_eq!(actual.len(), expected.len(), "{what}: length");
    for (i, (a, e)) in actual.iter().zip(expected).enumerate() {
        assert_eq!(a.to_bits(), e.to_bits(), "{what}: entry {i} ({a} vs {e})");
    }
}

#[test]
fn hierarchy_and_solution_are_bitwise_identical_across_thread_counts() {
    let (problem, q, rhs) = fixture();
    let with_threads = |threads: usize| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap()
            .install(|| run(&problem, &q, &rhs))
    };
    let base = with_threads(1);
    assert!(
        base.coarse_values.len() >= 2,
        "expected a multilevel hierarchy"
    );
    for threads in [2, 4] {
        let other = with_threads(threads);
        let tag = format!("{threads} threads");
        assert_bitwise(
            &other.lambda_max,
            &base.lambda_max,
            &format!("λ_max, {tag}"),
        );
        for (l, (a, b)) in other
            .coarse_values
            .iter()
            .zip(&base.coarse_values)
            .enumerate()
        {
            assert_bitwise(a, b, &format!("A_{} values, {tag}", l + 1));
        }
        assert_eq!(other.iterations, base.iterations, "iterations, {tag}");
        assert_bitwise(&other.x_cold, &base.x_cold, &format!("cold x, {tag}"));
        assert_bitwise(&other.x_warm, &base.x_warm, &format!("warm x, {tag}"));
    }
}

/// FNV-1a over the bit patterns.
fn fingerprint(v: &[f64]) -> u64 {
    v.iter().fold(0xcbf2_9ce4_8422_2325u64, |h, x| {
        x.to_bits().to_le_bytes().iter().fold(h, |h, b| {
            (h ^ u64::from(*b)).wrapping_mul(0x0000_0100_0000_01b3)
        })
    })
}

/// Helper entry point for the subprocess test below: prints fingerprints of
/// the run under whatever `RAYON_NUM_THREADS` the process was started with.
#[test]
fn print_fingerprint_for_current_thread_count() {
    if std::env::var_os("THESEUS_AMG_FINGERPRINT").is_none() {
        return;
    }
    let (problem, q, rhs) = fixture();
    let r = run(&problem, &q, &rhs);
    println!(
        "FINGERPRINT threads={} cold={:016x} warm={:016x} iters={:?}",
        rayon::current_num_threads(),
        fingerprint(&r.x_cold),
        fingerprint(&r.x_warm),
        r.iterations
    );
}

/// The literal `RAYON_NUM_THREADS=1,2,4` check: re-run this test binary as a
/// subprocess with the variable set and compare the printed fingerprints.
#[test]
fn rayon_num_threads_env_gives_identical_fingerprints() {
    let exe = std::env::current_exe().unwrap();
    let mut lines = Vec::new();
    for threads in [1usize, 2, 4] {
        let out = std::process::Command::new(&exe)
            .args([
                "--exact",
                "print_fingerprint_for_current_thread_count",
                "--nocapture",
                "--test-threads=1",
            ])
            .env("RAYON_NUM_THREADS", threads.to_string())
            .env("THESEUS_AMG_FINGERPRINT", "1")
            .output()
            .unwrap();
        assert!(
            out.status.success(),
            "{}",
            String::from_utf8_lossy(&out.stderr)
        );
        let stdout = String::from_utf8_lossy(&out.stdout);
        // With `--nocapture` the harness may print the test name on the
        // same line, so locate the marker rather than the line start.
        let line = stdout
            .lines()
            .find_map(|l| l.find("FINGERPRINT").map(|i| l[i..].to_string()))
            .unwrap_or_else(|| panic!("no fingerprint in:\n{stdout}"));
        assert!(line.contains(&format!("threads={threads} ")), "{line}");
        lines.push(line);
    }
    let strip = |l: &str| l.splitn(3, ' ').nth(2).unwrap().to_string();
    let base = strip(&lines[0]);
    for l in &lines[1..] {
        assert_eq!(strip(l), base, "{lines:?}");
    }
}

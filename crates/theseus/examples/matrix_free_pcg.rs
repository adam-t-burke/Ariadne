//! Matrix-free conjugate gradients versus the direct sparse Cholesky solve.
//!
//! The FDM equilibrium system `A(q) x = b` with `A = Cnᵀ diag(q) Cn` can be
//! applied without ever assembling `A`: one pass over the edges computes
//! `Cn x`, scales by `q`, and scatters back. This example measures whether an
//! iterative solve built on that operator can beat re-factorising `A` for
//! every optimizer evaluation, on the same grid networks used by
//! `tests/bench_scale.rs`:
//!
//! * Jacobi-preconditioned CG, cold (zero) start and warm start from the
//!   solution of the previous force densities `q_prev`;
//! * CG preconditioned by the *existing* Cholesky factor of `A(q_prev)`
//!   ("frozen factor"), which converges in a handful of iterations but pays a
//!   triangular solve per iteration.
//!
//! ```text
//! RAYON_NUM_THREADS=1 cargo run --release -p theseus --example matrix_free_pcg -- 224
//! ```
//!
//! Findings on the 224×224 grid (100k edges, 50k free nodes) are recorded in
//! `BENCHMARKS.md`; in short, the grid Laplacian's condition number grows with
//! the node count, so Jacobi-CG needs ~1,000 iterations (>1 s) while the direct
//! refactor plus solve takes ~45 ms, and the frozen-factor variant is no better
//! than refactoring because each of its iterations costs a triangular solve.

#[path = "../tests/support/grid.rs"]
#[allow(dead_code)]
mod grid;

use ndarray::Array2;
use std::time::Instant;
use theseus::types::*;

/// Matrix-free `A = Cnᵀ diag(q) Cn` on the free nodes, three coordinates at once.
struct Operator<'a> {
    q: &'a [f64],
    starts: &'a [usize],
    ends: &'a [usize],
    node_to_free: &'a [Option<usize>],
    nfree: usize,
}

impl Operator<'_> {
    /// `y = A x`; `x` and `y` are `nfree × 3` row-major.
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        y.fill(0.0);
        for k in 0..self.q.len() {
            let s = self.node_to_free[self.starts[k]];
            let e = self.node_to_free[self.ends[k]];
            let qk = self.q[k];
            for d in 0..3 {
                let xs = s.map_or(0.0, |i| x[i * 3 + d]);
                let xe = e.map_or(0.0, |i| x[i * 3 + d]);
                let f = qk * (xe - xs);
                if let Some(i) = e {
                    y[i * 3 + d] += f;
                }
                if let Some(i) = s {
                    y[i * 3 + d] -= f;
                }
            }
        }
    }

    fn diagonal(&self) -> Vec<f64> {
        let mut diag = vec![0.0; self.nfree];
        for k in 0..self.q.len() {
            if let Some(i) = self.node_to_free[self.starts[k]] {
                diag[i] += self.q[k];
            }
            if let Some(i) = self.node_to_free[self.ends[k]] {
                diag[i] += self.q[k];
            }
        }
        diag
    }
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

/// Preconditioned CG on three independent right-hand sides sharing one
/// operator application per iteration. Returns `(iterations, relative residual)`.
fn pcg(
    op: &Operator,
    b: &[f64],
    x: &mut [f64],
    tol: f64,
    max_iterations: usize,
    precondition: &mut dyn FnMut(&[f64], &mut [f64]),
) -> (usize, f64) {
    let n = b.len();
    let mut r = vec![0.0; n];
    let mut z = vec![0.0; n];
    let mut p = vec![0.0; n];
    let mut ap = vec![0.0; n];
    op.apply(x, &mut ap);
    for i in 0..n {
        r[i] = b[i] - ap[i];
    }
    let b_norm = dot(b, b).sqrt().max(1e-300);
    precondition(&r, &mut z);
    p.copy_from_slice(&z);
    let mut rz = dot(&r, &z);
    let mut iterations = 0;
    let mut relative = dot(&r, &r).sqrt() / b_norm;
    while relative > tol && iterations < max_iterations {
        op.apply(&p, &mut ap);
        let alpha = rz / dot(&p, &ap);
        for i in 0..n {
            x[i] += alpha * p[i];
            r[i] -= alpha * ap[i];
        }
        precondition(&r, &mut z);
        let rz_next = dot(&r, &z);
        let beta = rz_next / rz;
        rz = rz_next;
        for i in 0..n {
            p[i] = z[i] + beta * p[i];
        }
        iterations += 1;
        relative = dot(&r, &r).sqrt() / b_norm;
    }
    (iterations, relative)
}

fn main() {
    let n: usize = std::env::args()
        .nth(1)
        .map(|s| s.parse().expect("grid size"))
        .unwrap_or(224);
    let problem = grid::make_grid_problem(n);
    let ne = problem.topology.num_edges;
    let nfree = problem.topology.free_node_indices.len();
    let anchors = Array2::zeros((0, 3));

    // Deterministic pseudo-random force densities: q_prev is the "previous
    // iterate", q is q_prev perturbed by a relative step `eps`.
    let mut seed = 12345u64;
    let mut noise = move || {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        (seed % 10_000) as f64 / 10_000.0 - 0.5
    };
    let q_prev: Vec<f64> = (0..ne).map(|_| 1.0 + 0.5 * noise()).collect();
    let mut prev = FdmCache::new(&problem).unwrap();
    theseus::fdm::solve_fdm(&mut prev, &q_prev, &problem, &anchors, 0.0).unwrap();
    let x_prev: Vec<f64> = prev.x.iter().copied().collect();

    println!("grid {n}: {ne} edges, {nfree} free nodes");
    for eps in [0.0, 0.01, 0.05, 0.2] {
        let q: Vec<f64> = q_prev.iter().map(|&v| v * (1.0 + eps * noise())).collect();

        // Direct reference (also provides the right-hand side).
        let mut direct = FdmCache::new(&problem).unwrap();
        theseus::fdm::solve_fdm(&mut direct, &q, &problem, &anchors, 0.0).unwrap();
        let t = Instant::now();
        theseus::fdm::solve_fdm(&mut direct, &q, &problem, &anchors, 0.0).unwrap();
        let direct_ms = t.elapsed().as_secs_f64() * 1e3;
        let x_ref: Vec<f64> = direct.x.iter().copied().collect();
        let b: Vec<f64> = direct.rhs.iter().copied().collect();

        let op = Operator {
            q: &q,
            starts: &prev.edge_starts,
            ends: &prev.edge_ends,
            node_to_free: &prev.node_to_free_idx,
            nfree,
        };
        let mut y = vec![0.0; nfree * 3];
        let t = Instant::now();
        for _ in 0..10 {
            op.apply(&x_prev, &mut y);
        }
        let matvec_ms = t.elapsed().as_secs_f64() * 1e3 / 10.0;
        println!(
            "--- relative q step {eps}: direct refactor+solve {direct_ms:.1} ms, matrix-free A·x (3 rhs) {matvec_ms:.2} ms"
        );

        let diag = op.diagonal();
        let mut jacobi = |r: &[f64], z: &mut [f64]| {
            for i in 0..nfree {
                for d in 0..3 {
                    z[i * 3 + d] = r[i * 3 + d] / diag[i];
                }
            }
        };
        let factor = prev.factorization().unwrap();
        let mut work = vec![0.0; nfree * 6];
        let mut frozen = |r: &[f64], z: &mut [f64]| factor.solve_slices::<3>(r, z, &mut work);

        let report = |label: &str,
                      tol: f64,
                      iterations: usize,
                      ms: f64,
                      relative: f64,
                      x: &[f64]| {
            let max_err = x
                .iter()
                .zip(&x_ref)
                .map(|(a, b)| (a - b).abs())
                .fold(0.0, f64::max);
            println!(
                "  {label:<22} tol={tol:.0e}: iterations={iterations:>5} time={ms:>8.1} ms rel_res={relative:.1e} max|x-x_direct|={max_err:.1e}"
            );
        };

        for tol in [1e-6, 1e-8, 1e-10] {
            let mut x = vec![0.0; nfree * 3];
            let t = Instant::now();
            let (it, rel) = pcg(&op, &b, &mut x, tol, 5000, &mut jacobi);
            report(
                "jacobi, cold start",
                tol,
                it,
                t.elapsed().as_secs_f64() * 1e3,
                rel,
                &x,
            );
        }
        for tol in [1e-6, 1e-8, 1e-10] {
            let mut x = x_prev.clone();
            let t = Instant::now();
            let (it, rel) = pcg(&op, &b, &mut x, tol, 5000, &mut jacobi);
            report(
                "jacobi, warm start",
                tol,
                it,
                t.elapsed().as_secs_f64() * 1e3,
                rel,
                &x,
            );
        }
        for tol in [1e-6, 1e-8, 1e-10] {
            let mut x = x_prev.clone();
            let t = Instant::now();
            let (it, rel) = pcg(&op, &b, &mut x, tol, 200, &mut frozen);
            report(
                "frozen factor, warm",
                tol,
                it,
                t.elapsed().as_secs_f64() * 1e3,
                rel,
                &x,
            );
        }
    }
}

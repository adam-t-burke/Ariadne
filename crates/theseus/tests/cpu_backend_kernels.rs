//! `backend::cpu`: every kernel against a scalar reference, in `f64` and
//! `f32`, and bitwise identical results under rayon pools of 1, 2 and 4
//! threads.

use theseus::backend::cpu::{CpuBackend, CpuBuf, CHUNK, PAR_MIN_LEN};
use theseus::backend::Backend;
use theseus::graph::{BlockVec, CsrAdjacency, LevelGraph};
use theseus::linear_solver::Precision;

struct Lcg(u64);

impl Lcg {
    fn next_u64(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0 >> 11
    }
    fn below(&mut self, n: usize) -> usize {
        (self.next_u64() % n as u64) as usize
    }
    fn unit(&mut self) -> f64 {
        (self.next_u64() % 1_000_000) as f64 / 1_000_000.0
    }
    fn signed(&mut self) -> f64 {
        2.0 * self.unit() - 1.0
    }
    fn vec(&mut self, len: usize) -> Vec<f64> {
        (0..len).map(|_| self.signed()).collect()
    }
}

/// Random connected weighted graph on `n` nodes with positive anchors on a
/// few nodes.
fn random_level(rng: &mut Lcg, n: usize) -> LevelGraph {
    let mut starts = Vec::new();
    let mut ends = Vec::new();
    for v in 1..n {
        starts.push(rng.below(v));
        ends.push(v);
    }
    for _ in 0..n {
        let u = rng.below(n);
        let v = rng.below(n);
        if u != v {
            starts.push(u);
            ends.push(v);
        }
    }
    let adjacency = CsrAdjacency::from_endpoints(n, &starts, &ends);
    let weight: Vec<f64> = (0..starts.len()).map(|_| 0.1 + 5.0 * rng.unit()).collect();
    let anchor: Vec<f64> = (0..n)
        .map(|_| {
            if rng.below(10) == 0 {
                1.0 + rng.unit()
            } else {
                0.0
            }
        })
        .collect();
    LevelGraph {
        n,
        adjacency,
        weight,
        anchor,
        aggregate_of: Vec::new(),
    }
}

fn random_aggregates(rng: &mut Lcg, n_fine: usize, n_coarse: usize) -> Vec<u32> {
    (0..n_fine).map(|_| rng.below(n_coarse) as u32).collect()
}

fn upload(be: &CpuBackend, v: &[f64], p: Precision) -> CpuBuf {
    let mut buf = be.alloc(v.len(), p);
    be.upload(v, &mut buf);
    buf
}

fn assert_close(actual: &[f64], expected: &[f64], tol: f64, what: &str) {
    assert_eq!(actual.len(), expected.len(), "{what}: length");
    let scale = expected
        .iter()
        .fold(0.0f64, |m, v| m.max(v.abs()))
        .max(1e-300);
    for (i, (a, e)) in actual.iter().zip(expected).enumerate() {
        assert!(
            (a - e).abs() <= tol * scale,
            "{what}: entry {i}: {a} vs {e} (tol {tol:e} × {scale:e})"
        );
    }
}

fn assert_bitwise(actual: &[f64], expected: &[f64], what: &str) {
    assert_eq!(actual.len(), expected.len(), "{what}: length");
    for (i, (a, e)) in actual.iter().zip(expected).enumerate() {
        assert_eq!(a.to_bits(), e.to_bits(), "{what}: entry {i}: {a} vs {e}");
    }
}

// ── Scalar references ────────────────────────────────────

fn ref_apply(level: &LevelGraph, x: &[f64]) -> Vec<f64> {
    let xb = BlockVec::<3>::from_vec(x.to_vec());
    let mut y = BlockVec::<3>::zeros(level.n);
    level.apply(&xb, &mut y);
    y.data
}

fn ref_chebyshev(
    level: &LevelGraph,
    alpha: f64,
    beta: f64,
    r: &[f64],
    d: &mut [f64],
    x: &mut [f64],
) {
    for u in 0..level.n {
        let inv = 1.0 / level.diagonal(u);
        for k in 0..3 {
            let i = u * 3 + k;
            let z = inv * r[i];
            d[i] = alpha * z + beta * d[i];
            x[i] += d[i];
        }
    }
}

fn ref_restrict(agg: &[u32], n_coarse: usize, fine: &[f64]) -> Vec<f64> {
    let mut coarse = vec![0.0; n_coarse * 3];
    for (i, &c) in agg.iter().enumerate() {
        for k in 0..3 {
            coarse[c as usize * 3 + k] += fine[i * 3 + k];
        }
    }
    coarse
}

fn ref_dot3(a: &[f64], b: &[f64]) -> [f64; 3] {
    let mut acc = [0.0; 3];
    for (ra, rb) in a.chunks_exact(3).zip(b.chunks_exact(3)) {
        for k in 0..3 {
            acc[k] += ra[k] * rb[k];
        }
    }
    acc
}

/// The documented reduction order: sequential per chunk of `CHUNK` rows,
/// partials added in chunk order.
fn chunked_dot3(a: &[f64], b: &[f64]) -> [f64; 3] {
    let mut total = [0.0; 3];
    for (ca, cb) in a.chunks(CHUNK * 3).zip(b.chunks(CHUNK * 3)) {
        let p = ref_dot3(ca, cb);
        for k in 0..3 {
            total[k] += p[k];
        }
    }
    total
}

// ── Tests ────────────────────────────────────────────────

/// Sizes straddling the sequential/parallel threshold and the chunk size.
fn sizes() -> Vec<usize> {
    vec![7, CHUNK - 1, CHUNK + 5, PAR_MIN_LEN / 3 + 100, 40_000]
}

#[test]
fn f64_kernels_match_scalar_reference_bitwise() {
    let be = CpuBackend::new();
    let mut rng = Lcg(3);
    for n in sizes() {
        let level = random_level(&mut rng, n);
        let lvl = be.upload_level(&level, Precision::F64);
        let x = rng.vec(n * 3);
        let b = rng.vec(n * 3);

        let xb = upload(&be, &x, Precision::F64);
        let bb = upload(&be, &b, Precision::F64);
        let mut y = be.alloc(n * 3, Precision::F64);
        be.apply_graph(&lvl, &xb, &mut y);
        let y_ref = ref_apply(&level, &x);
        assert_bitwise(y.as_f64().unwrap(), &y_ref, &format!("apply_graph n={n}"));

        let mut r = be.alloc(n * 3, Precision::F64);
        be.residual(&lvl, &xb, &bb, &mut r);
        let r_ref: Vec<f64> = b.iter().zip(&y_ref).map(|(b, y)| b - y).collect();
        assert_bitwise(r.as_f64().unwrap(), &r_ref, &format!("residual n={n}"));

        let (alpha, beta) = (0.37, -0.61);
        let d0 = rng.vec(n * 3);
        let mut d = upload(&be, &d0, Precision::F64);
        let mut xk = upload(&be, &x, Precision::F64);
        be.chebyshev_step(&lvl, alpha, beta, &r, &mut d, &mut xk);
        let mut d_ref = d0.clone();
        let mut x_ref = x.clone();
        ref_chebyshev(&level, alpha, beta, &r_ref, &mut d_ref, &mut x_ref);
        assert_bitwise(d.as_f64().unwrap(), &d_ref, &format!("chebyshev d n={n}"));
        assert_bitwise(xk.as_f64().unwrap(), &x_ref, &format!("chebyshev x n={n}"));

        let n_coarse = (n / 3).max(1);
        let agg = random_aggregates(&mut rng, n, n_coarse);
        let aggb = be.upload_aggregates(&agg, n_coarse);
        let mut coarse = be.alloc(n_coarse * 3, Precision::F64);
        be.restrict(&aggb, &xb, &mut coarse);
        let coarse_ref = ref_restrict(&agg, n_coarse, &x);
        assert_bitwise(
            coarse.as_f64().unwrap(),
            &coarse_ref,
            &format!("restrict n={n}"),
        );

        let mut fine = upload(&be, &b, Precision::F64);
        be.prolong_add(&aggb, &coarse, &mut fine);
        let fine_ref: Vec<f64> = (0..n * 3)
            .map(|i| b[i] + coarse_ref[agg[i / 3] as usize * 3 + i % 3])
            .collect();
        assert_bitwise(
            fine.as_f64().unwrap(),
            &fine_ref,
            &format!("prolong_add n={n}"),
        );

        let al = [1.5, -2.0, 0.25];
        let mut yy = upload(&be, &b, Precision::F64);
        be.axpy(al, &xb, &mut yy);
        let yy_ref: Vec<f64> = (0..n * 3).map(|i| b[i] + al[i % 3] * x[i]).collect();
        assert_bitwise(yy.as_f64().unwrap(), &yy_ref, &format!("axpy n={n}"));

        let mut sc = upload(&be, &x, Precision::F64);
        be.scale(al, &mut sc);
        let sc_ref: Vec<f64> = (0..n * 3).map(|i| al[i % 3] * x[i]).collect();
        assert_bitwise(sc.as_f64().unwrap(), &sc_ref, &format!("scale n={n}"));

        let dot = be.dot3(&xb, &bb);
        assert_eq!(dot, chunked_dot3(&x, &b), "dot3 chunk order n={n}");
        assert_close(&dot, &ref_dot3(&x, &b), 1e-12, &format!("dot3 n={n}"));
        let norm = be.norm3(&xb);
        let nref = ref_dot3(&x, &x).map(f64::sqrt);
        assert_close(&norm, &nref, 1e-12, &format!("norm3 n={n}"));

        // Buffer plumbing.
        let mut out = vec![0.0; n * 3];
        be.download(&xb, &mut out);
        assert_bitwise(&out, &x, "download");
        let mut copy = be.alloc(n * 3, Precision::F64);
        be.copy(&xb, &mut copy);
        assert_eq!(copy, xb);
        be.zero(&mut copy);
        assert!(copy.as_f64().unwrap().iter().all(|v| *v == 0.0));
        assert_eq!(be.len(&copy), n * 3);
    }
}

#[test]
fn f32_kernels_match_f64_reference() {
    let be = CpuBackend::new();
    let mut rng = Lcg(5);
    for n in [CHUNK + 5, 30_000] {
        let level = random_level(&mut rng, n);
        let lvl = be.upload_level(&level, Precision::F32);
        assert_eq!(lvl.precision(), Precision::F32);
        assert_eq!(lvl.n(), n);
        let x = rng.vec(n * 3);
        let b = rng.vec(n * 3);
        let xb = upload(&be, &x, Precision::F32);
        let bb = upload(&be, &b, Precision::F32);

        let mut y = be.alloc(n * 3, Precision::F32);
        be.apply_graph(&lvl, &xb, &mut y);
        assert_close(
            &y.to_vec_f64(),
            &ref_apply(&level, &x),
            2e-6,
            "f32 apply_graph",
        );

        let mut r = be.alloc(n * 3, Precision::F32);
        be.residual(&lvl, &xb, &bb, &mut r);
        let y_ref = ref_apply(&level, &x);
        let r_ref: Vec<f64> = b.iter().zip(&y_ref).map(|(b, y)| b - y).collect();
        assert_close(&r.to_vec_f64(), &r_ref, 2e-6, "f32 residual");

        let d0 = rng.vec(n * 3);
        let mut d = upload(&be, &d0, Precision::F32);
        let mut xk = upload(&be, &x, Precision::F32);
        be.chebyshev_step(&lvl, 0.8, 0.3, &r, &mut d, &mut xk);
        let mut d_ref = d0.clone();
        let mut x_ref = x.clone();
        ref_chebyshev(&level, 0.8, 0.3, &r_ref, &mut d_ref, &mut x_ref);
        assert_close(&d.to_vec_f64(), &d_ref, 2e-6, "f32 chebyshev d");
        assert_close(&xk.to_vec_f64(), &x_ref, 2e-6, "f32 chebyshev x");

        let n_coarse = n / 4;
        let agg = random_aggregates(&mut rng, n, n_coarse);
        let aggb = be.upload_aggregates(&agg, n_coarse);
        let mut coarse = be.alloc(n_coarse * 3, Precision::F32);
        be.restrict(&aggb, &xb, &mut coarse);
        assert_close(
            &coarse.to_vec_f64(),
            &ref_restrict(&agg, n_coarse, &x),
            2e-6,
            "f32 restrict",
        );

        let dot = be.dot3(&xb, &bb);
        // f32 inputs, f64 accumulation: only the rounding of the inputs.
        let xs: Vec<f64> = x.iter().map(|v| *v as f32 as f64).collect();
        let bs: Vec<f64> = b.iter().map(|v| *v as f32 as f64).collect();
        assert_close(&dot, &ref_dot3(&xs, &bs), 1e-12, "f32 dot3");

        // Weight refresh keeps the precision and the diagonal consistent.
        let mut lvl2 = lvl.clone();
        let w2: Vec<f64> = level.weight.iter().map(|w| 2.0 * w).collect();
        let a2: Vec<f64> = level.anchor.iter().map(|a| 2.0 * a).collect();
        be.update_level_weights(&mut lvl2, &w2, &a2);
        let mut y2 = be.alloc(n * 3, Precision::F32);
        be.apply_graph(&lvl2, &xb, &mut y2);
        let y2_ref: Vec<f64> = ref_apply(&level, &x).iter().map(|v| 2.0 * v).collect();
        assert_close(&y2.to_vec_f64(), &y2_ref, 2e-6, "f32 apply after update");
    }
}

#[test]
fn update_level_weights_matches_fresh_upload() {
    let be = CpuBackend::new();
    let mut rng = Lcg(9);
    let n = 20_000;
    let level = random_level(&mut rng, n);
    let mut level2 = level.clone();
    for w in &mut level2.weight {
        *w *= 0.5 + rng.unit();
    }
    for a in &mut level2.anchor {
        *a *= 0.5 + rng.unit();
    }
    let mut updated = be.upload_level(&level, Precision::F64);
    be.update_level_weights(&mut updated, &level2.weight, &level2.anchor);
    let fresh = be.upload_level(&level2, Precision::F64);
    assert_eq!(updated, fresh);
}

/// Every kernel output is bitwise identical for 1, 2 and 4 rayon threads.
#[test]
fn kernels_are_bitwise_identical_across_thread_counts() {
    let mut rng = Lcg(21);
    let n = 50_000;
    let level = random_level(&mut rng, n);
    let x = rng.vec(n * 3);
    let b = rng.vec(n * 3);
    let d0 = rng.vec(n * 3);
    let n_coarse = n / 3;
    let agg = random_aggregates(&mut rng, n, n_coarse);

    let run = |threads: usize| -> Vec<Vec<f64>> {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        pool.install(|| {
            let be = CpuBackend::new();
            let lvl = be.upload_level(&level, Precision::F64);
            let aggb = be.upload_aggregates(&agg, n_coarse);
            let xb = upload(&be, &x, Precision::F64);
            let bb = upload(&be, &b, Precision::F64);
            let mut y = be.alloc(n * 3, Precision::F64);
            be.apply_graph(&lvl, &xb, &mut y);
            let mut r = be.alloc(n * 3, Precision::F64);
            be.residual(&lvl, &xb, &bb, &mut r);
            let mut d = upload(&be, &d0, Precision::F64);
            let mut xk = upload(&be, &x, Precision::F64);
            be.chebyshev_step(&lvl, 0.9, -0.2, &r, &mut d, &mut xk);
            let mut coarse = be.alloc(n_coarse * 3, Precision::F64);
            be.restrict(&aggb, &xb, &mut coarse);
            let mut fine = upload(&be, &b, Precision::F64);
            be.prolong_add(&aggb, &coarse, &mut fine);
            let mut yy = upload(&be, &b, Precision::F64);
            be.axpy([0.1, 0.2, 0.3], &xb, &mut yy);
            let mut sc = upload(&be, &x, Precision::F64);
            be.scale([3.0, -1.0, 0.5], &mut sc);
            let dot = be.dot3(&xb, &bb);
            let norm = be.norm3(&r);
            vec![
                y.to_vec_f64(),
                r.to_vec_f64(),
                d.to_vec_f64(),
                xk.to_vec_f64(),
                coarse.to_vec_f64(),
                fine.to_vec_f64(),
                yy.to_vec_f64(),
                sc.to_vec_f64(),
                dot.to_vec(),
                norm.to_vec(),
            ]
        })
    };

    let base = run(1);
    for threads in [2, 4] {
        let other = run(threads);
        for (k, (a, b)) in base.iter().zip(&other).enumerate() {
            assert_bitwise(b, a, &format!("output {k}, {threads} threads"));
        }
    }
}

// ── CSR kernels (WS-C additions to `Backend`) ────────────

mod csr {
    use super::*;
    use theseus::amg::hierarchy::LevelMatrix;

    /// Random rectangular CSR (`nrows × ncols`, a few entries per row) and a
    /// random square SPD-like one with a full diagonal.
    fn random_rect(rng: &mut Lcg, nrows: usize, ncols: usize) -> LevelMatrix {
        let mut rows: Vec<Vec<(u32, f64)>> = (0..nrows)
            .map(|_| {
                (0..1 + rng.below(4))
                    .map(|_| (rng.below(ncols) as u32, rng.signed()))
                    .collect()
            })
            .collect();
        let m = LevelMatrix::from_rows(nrows, ncols, &mut rows);
        m.check();
        m
    }

    fn random_square(rng: &mut Lcg, n: usize) -> LevelMatrix {
        let mut rows: Vec<Vec<(u32, f64)>> = (0..n)
            .map(|u| {
                let mut row: Vec<(u32, f64)> = (0..rng.below(5))
                    .map(|_| (rng.below(n) as u32, -rng.unit()))
                    .collect();
                row.retain(|e| e.0 as usize != u);
                let off: f64 = row.iter().map(|e| e.1.abs()).sum();
                row.push((u as u32, off + 0.5 + rng.unit()));
                row
            })
            .collect();
        let m = LevelMatrix::from_rows(n, n, &mut rows);
        m.check();
        assert!(m.diag.iter().all(|&d| d > 0.0));
        m
    }

    fn ref_apply(m: &LevelMatrix, x: &[f64]) -> Vec<f64> {
        let mut y = vec![0.0; m.n * 3];
        m.apply(x, &mut y, 3);
        y
    }

    fn ref_chebyshev(
        m: &LevelMatrix,
        alpha: f64,
        beta: f64,
        r: &[f64],
        d: &mut [f64],
        x: &mut [f64],
    ) {
        for u in 0..m.n {
            let inv = 1.0 / m.diag[u];
            for k in 0..3 {
                let i = u * 3 + k;
                d[i] = alpha * (inv * r[i]) + beta * d[i];
                x[i] += d[i];
            }
        }
    }

    #[test]
    fn f64_csr_kernels_match_scalar_reference_bitwise() {
        let be = CpuBackend::new();
        let mut rng = Lcg(77);
        for &n in &[
            1usize,
            7,
            CHUNK - 1,
            CHUNK + 1,
            PAR_MIN_LEN + 5,
            3 * PAR_MIN_LEN + 11,
        ] {
            let nc = n.div_ceil(3).max(1);
            let p = random_rect(&mut rng, n, nc);
            let pt = p.transpose();
            let a = random_square(&mut rng, n);
            let xc = rng.vec(nc * 3);
            let x = rng.vec(n * 3);
            let b = rng.vec(n * 3);
            let d0 = rng.vec(n * 3);

            let pb = be.upload_csr(&p, Precision::F64);
            let ptb = be.upload_csr(&pt, Precision::F64);
            let ab = be.upload_csr(&a, Precision::F64);
            let xcb = upload(&be, &xc, Precision::F64);
            let xb = upload(&be, &x, Precision::F64);
            let bb = upload(&be, &b, Precision::F64);

            // y = P x_c  (prolongation)
            let mut y = be.alloc(n * 3, Precision::F64);
            be.apply_csr(&pb, &xcb, &mut y);
            assert_bitwise(
                &y.to_vec_f64(),
                &ref_apply(&p, &xc),
                &format!("apply_csr P, n={n}"),
            );
            // y += P x_c
            let mut y_add = upload(&be, &b, Precision::F64);
            be.apply_csr_add(&pb, &xcb, &mut y_add);
            let expected: Vec<f64> = b
                .iter()
                .zip(ref_apply(&p, &xc))
                .map(|(b, y)| b + y)
                .collect();
            assert_bitwise(
                &y_add.to_vec_f64(),
                &expected,
                &format!("apply_csr_add P, n={n}"),
            );
            // r_c = Pᵀ x  (restriction)
            let mut rc = be.alloc(nc * 3, Precision::F64);
            be.apply_csr(&ptb, &xb, &mut rc);
            assert_bitwise(
                &rc.to_vec_f64(),
                &ref_apply(&pt, &x),
                &format!("apply_csr Pᵀ, n={n}"),
            );
            // y = A x
            let mut ya = be.alloc(n * 3, Precision::F64);
            be.apply_csr(&ab, &xb, &mut ya);
            assert_bitwise(
                &ya.to_vec_f64(),
                &ref_apply(&a, &x),
                &format!("apply_csr A, n={n}"),
            );
            // r = b − A x
            let mut r = be.alloc(n * 3, Precision::F64);
            be.residual_csr(&ab, &xb, &bb, &mut r);
            let expected: Vec<f64> = b
                .iter()
                .zip(ref_apply(&a, &x))
                .map(|(b, y)| b - y)
                .collect();
            assert_bitwise(&r.to_vec_f64(), &expected, &format!("residual_csr, n={n}"));
            // Chebyshev step with diag(A)
            let mut d = upload(&be, &d0, Precision::F64);
            let mut xk = upload(&be, &x, Precision::F64);
            be.chebyshev_step_csr(&ab, 0.8, -0.3, &r, &mut d, &mut xk);
            let mut d_ref = d0.clone();
            let mut x_ref = x.clone();
            ref_chebyshev(&a, 0.8, -0.3, &expected, &mut d_ref, &mut x_ref);
            assert_bitwise(
                &d.to_vec_f64(),
                &d_ref,
                &format!("chebyshev_step_csr d, n={n}"),
            );
            assert_bitwise(
                &xk.to_vec_f64(),
                &x_ref,
                &format!("chebyshev_step_csr x, n={n}"),
            );
        }
    }

    #[test]
    fn f32_csr_kernels_match_f64_reference() {
        let be = CpuBackend::new();
        let mut rng = Lcg(78);
        let n = PAR_MIN_LEN + 3;
        let nc = n / 3;
        let p = random_rect(&mut rng, n, nc);
        let a = random_square(&mut rng, n);
        let xc = rng.vec(nc * 3);
        let x = rng.vec(n * 3);
        let b = rng.vec(n * 3);
        let pb = be.upload_csr(&p, Precision::F32);
        let ab = be.upload_csr(&a, Precision::F32);
        let xcb = upload(&be, &xc, Precision::F32);
        let xb = upload(&be, &x, Precision::F32);
        let bb = upload(&be, &b, Precision::F32);
        let mut y = be.alloc(n * 3, Precision::F32);
        be.apply_csr(&pb, &xcb, &mut y);
        assert_close(&y.to_vec_f64(), &ref_apply(&p, &xc), 1e-5, "apply_csr f32");
        let mut r = be.alloc(n * 3, Precision::F32);
        be.residual_csr(&ab, &xb, &bb, &mut r);
        let expected: Vec<f64> = b
            .iter()
            .zip(ref_apply(&a, &x))
            .map(|(b, y)| b - y)
            .collect();
        assert_close(&r.to_vec_f64(), &expected, 1e-5, "residual_csr f32");
        let mut d = be.alloc(n * 3, Precision::F32);
        let mut xk = upload(&be, &x, Precision::F32);
        be.chebyshev_step_csr(&ab, 0.8, 0.0, &r, &mut d, &mut xk);
        let mut d_ref = vec![0.0; n * 3];
        let mut x_ref = x.clone();
        ref_chebyshev(&a, 0.8, 0.0, &expected, &mut d_ref, &mut x_ref);
        assert_close(&xk.to_vec_f64(), &x_ref, 1e-5, "chebyshev_step_csr f32");
    }

    #[test]
    fn update_csr_values_matches_fresh_upload() {
        let be = CpuBackend::new();
        let mut rng = Lcg(79);
        let n = 5000;
        let a = random_square(&mut rng, n);
        let mut a2 = a.clone();
        for v in &mut a2.values {
            *v *= 0.5 + rng.unit();
        }
        a2.refresh_diag();
        for precision in [Precision::F64, Precision::F32] {
            let mut updated = be.upload_csr(&a, precision);
            be.update_csr_values(&a2, &mut updated);
            let fresh = be.upload_csr(&a2, precision);
            assert_eq!(updated, fresh, "{precision:?}");
        }
    }

    #[test]
    fn csr_kernels_are_bitwise_identical_across_thread_counts() {
        let mut rng = Lcg(80);
        let n = 50_000;
        let nc = n / 3;
        let p = random_rect(&mut rng, n, nc);
        let pt = p.transpose();
        let a = random_square(&mut rng, n);
        let xc = rng.vec(nc * 3);
        let x = rng.vec(n * 3);
        let b = rng.vec(n * 3);
        let run = |threads: usize| -> Vec<Vec<f64>> {
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap()
                .install(|| {
                    let be = CpuBackend::new();
                    let pb = be.upload_csr(&p, Precision::F64);
                    let ptb = be.upload_csr(&pt, Precision::F64);
                    let ab = be.upload_csr(&a, Precision::F64);
                    let xcb = upload(&be, &xc, Precision::F64);
                    let xb = upload(&be, &x, Precision::F64);
                    let bb = upload(&be, &b, Precision::F64);
                    let mut y = be.alloc(n * 3, Precision::F64);
                    be.apply_csr(&pb, &xcb, &mut y);
                    let mut rc = be.alloc(nc * 3, Precision::F64);
                    be.apply_csr(&ptb, &xb, &mut rc);
                    let mut r = be.alloc(n * 3, Precision::F64);
                    be.residual_csr(&ab, &xb, &bb, &mut r);
                    let mut d = upload(&be, &b, Precision::F64);
                    let mut xk = upload(&be, &x, Precision::F64);
                    be.chebyshev_step_csr(&ab, 0.7, -0.1, &r, &mut d, &mut xk);
                    vec![
                        y.to_vec_f64(),
                        rc.to_vec_f64(),
                        r.to_vec_f64(),
                        d.to_vec_f64(),
                        xk.to_vec_f64(),
                    ]
                })
        };
        let base = run(1);
        for threads in [2, 4] {
            let other = run(threads);
            for (k, (a, b)) in base.iter().zip(&other).enumerate() {
                assert_bitwise(b, a, &format!("csr output {k}, {threads} threads"));
            }
        }
    }
}

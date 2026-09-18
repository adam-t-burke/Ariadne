//! GPU kernel equivalence tests (program plan §4.2): every `GpuBackend`
//! kernel against a scalar reference on random inputs and random small
//! graphs, in `f32` (`1e-6` relative) and, where the adapter has
//! `SHADER_F64`, in `f64` (`1e-12` relative); plus an apply → residual →
//! norm chain on a 64×64 grid level built by hand from `CsrAdjacency`.
//!
//! Skipped (with a printed reason) when no adapter is usable. Software
//! adapters need `THESEUS_GPU_ALLOW_SOFTWARE=1`; on Linux CI that is
//! lavapipe (`mesa-vulkan-drivers`).
#![cfg(feature = "gpu")]

#[path = "support/fixtures/mod.rs"]
mod fixtures;

use theseus::backend::gpu::{members_csr, AdapterPolicy, GpuBackend, GpuContext};
use theseus::backend::{probe_gpu, Backend};
use theseus::graph::{BlockVec, CsrAdjacency, LevelGraph};
use theseus::linear_solver::Precision;

// ── Harness ────────────────────────────────────────────────────────────

fn backend() -> Option<GpuBackend> {
    match GpuBackend::from_env() {
        Ok(gpu) => Some(gpu),
        Err(e) => {
            eprintln!("skipping GPU test: {e}");
            None
        }
    }
}

fn precisions(gpu: &GpuBackend) -> Vec<Precision> {
    if gpu.shader_f64() {
        vec![Precision::F32, Precision::F64]
    } else {
        eprintln!("adapter has no SHADER_F64: f64 variants skipped");
        vec![Precision::F32]
    }
}

fn tol(precision: Precision) -> f64 {
    match precision {
        Precision::F32 => 1e-6,
        Precision::F64 => 1e-12,
    }
}

/// Round `v` to the storage precision so the reference sees the same inputs
/// the device does.
fn round(precision: Precision, v: f64) -> f64 {
    match precision {
        Precision::F32 => f64::from(v as f32),
        Precision::F64 => v,
    }
}

fn round_all(precision: Precision, v: &[f64]) -> Vec<f64> {
    v.iter().map(|&x| round(precision, x)).collect()
}

/// `max |a − b| ≤ tol · max(scale, tiny)`.
fn assert_close(what: &str, got: &[f64], want: &[f64], scale: f64, tol: f64) {
    assert_eq!(got.len(), want.len(), "{what}: length");
    let scale = scale.max(1e-300);
    let mut worst = 0.0f64;
    let mut at = 0;
    for (i, (g, w)) in got.iter().zip(want).enumerate() {
        let d = (g - w).abs();
        if d > worst {
            worst = d;
            at = i;
        }
    }
    assert!(
        worst <= tol * scale,
        "{what}: max |Δ| = {worst:.3e} at {at} (got {}, want {}), scale {scale:.3e}, tol {tol:.0e}",
        got[at],
        want[at]
    );
    // Visible with `--nocapture`: the relative error actually achieved.
    eprintln!(
        "{what}: max |Δ| / scale = {:.2e} (tol {tol:.0e})",
        worst / scale
    );
}

fn inf_norm(v: &[f64]) -> f64 {
    v.iter().fold(0.0, |m, x| m.max(x.abs()))
}

/// SplitMix64 — no dev-dependency, reproducible across platforms.
struct Rng(u64);

impl Rng {
    fn new(seed: u64) -> Self {
        Self(seed.wrapping_add(0x9E37_79B9_7F4A_7C15))
    }

    fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Uniform in `[0, 1)`.
    fn unit(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.unit()
    }

    fn below(&mut self, n: usize) -> usize {
        (self.next_u64() % n as u64) as usize
    }

    fn vec(&mut self, len: usize, lo: f64, hi: f64) -> Vec<f64> {
        (0..len).map(|_| self.range(lo, hi)).collect()
    }
}

// ── Graph construction ─────────────────────────────────────────────────

/// Node-centred CSR from an edge list `(start, end)` over `n` nodes.
fn csr_from_edges(n: usize, edges: &[(u32, u32)]) -> CsrAdjacency {
    let mut degree = vec![0u32; n + 1];
    for &(s, e) in edges {
        degree[s as usize + 1] += 1;
        degree[e as usize + 1] += 1;
    }
    for i in 0..n {
        degree[i + 1] += degree[i];
    }
    let offsets = degree;
    let mut next = offsets.clone();
    let total = offsets[n] as usize;
    let mut out_edges = vec![0u32; total];
    let mut sign = vec![0i8; total];
    let mut other = vec![0u32; total];
    for (e, &(s, t)) in edges.iter().enumerate() {
        let slot = next[s as usize] as usize;
        out_edges[slot] = e as u32;
        sign[slot] = -1;
        other[slot] = t;
        next[s as usize] += 1;
        let slot = next[t as usize] as usize;
        out_edges[slot] = e as u32;
        sign[slot] = 1;
        other[slot] = s;
        next[t as usize] += 1;
    }
    CsrAdjacency {
        offsets,
        edges: out_edges,
        sign,
        other,
    }
}

/// Random connected-ish graph: a random spanning path plus `extra` random
/// edges (multi-edges allowed), weights in `[0.1, 2)`, anchors in `[0, 1)`
/// on ~⅓ of the nodes.
fn random_level(rng: &mut Rng, n: usize, extra: usize) -> LevelGraph {
    let mut edges = Vec::with_capacity(n + extra);
    for i in 1..n {
        edges.push((rng.below(i) as u32, i as u32));
    }
    for _ in 0..extra {
        let u = rng.below(n);
        let mut v = rng.below(n);
        if u == v {
            v = (v + 1) % n;
        }
        edges.push((u as u32, v as u32));
    }
    let weight = rng.vec(edges.len(), 0.1, 2.0);
    let anchor = (0..n)
        .map(|_| {
            if rng.unit() < 0.34 {
                rng.range(0.0, 1.0)
            } else {
                0.0
            }
        })
        .collect();
    LevelGraph {
        n,
        adjacency: csr_from_edges(n, &edges),
        weight,
        anchor,
        aggregate_of: Vec::new(),
    }
}

/// `w × h` grid with unit weights (jittered) and anchors on the boundary
/// nodes, standing for edges to fixed nodes outside the grid.
fn grid_level(rng: &mut Rng, w: usize, h: usize) -> LevelGraph {
    let id = |i: usize, j: usize| (i * w + j) as u32;
    let mut edges = Vec::new();
    for i in 0..h {
        for j in 0..w {
            if j + 1 < w {
                edges.push((id(i, j), id(i, j + 1)));
            }
            if i + 1 < h {
                edges.push((id(i, j), id(i + 1, j)));
            }
        }
    }
    let weight = rng.vec(edges.len(), 0.8, 1.2);
    let mut anchor = vec![0.0; w * h];
    for i in 0..h {
        for j in 0..w {
            let mut a = 0.0;
            if i == 0 || i + 1 == h {
                a += rng.range(0.8, 1.2);
            }
            if j == 0 || j + 1 == w {
                a += rng.range(0.8, 1.2);
            }
            anchor[i * w + j] = a;
        }
    }
    LevelGraph {
        n: w * h,
        adjacency: csr_from_edges(w * h, &edges),
        weight,
        anchor,
        aggregate_of: Vec::new(),
    }
}

fn rounded_level(precision: Precision, level: &LevelGraph) -> LevelGraph {
    LevelGraph {
        weight: round_all(precision, &level.weight),
        anchor: round_all(precision, &level.anchor),
        ..level.clone()
    }
}

/// Random aggregate map with `n_coarse` non-empty aggregates.
fn random_aggregates(rng: &mut Rng, n_fine: usize, n_coarse: usize) -> Vec<u32> {
    assert!(n_coarse <= n_fine);
    let mut agg: Vec<u32> = (0..n_fine)
        .map(|i| {
            if i < n_coarse {
                i as u32
            } else {
                rng.below(n_coarse) as u32
            }
        })
        .collect();
    // Shuffle so aggregates are not index-sorted.
    for i in (1..n_fine).rev() {
        let j = rng.below(i + 1);
        agg.swap(i, j);
    }
    agg
}

// ── Scalar references ──────────────────────────────────────────────────

fn ref_apply(level: &LevelGraph, x: &[f64]) -> Vec<f64> {
    let xb = BlockVec::<3>::from_vec(x.to_vec());
    let mut y = BlockVec::<3>::zeros(level.n);
    level.apply(&xb, &mut y);
    y.data
}

fn ref_residual(level: &LevelGraph, x: &[f64], b: &[f64]) -> Vec<f64> {
    ref_apply(level, x)
        .iter()
        .zip(b)
        .map(|(ax, b)| b - ax)
        .collect()
}

fn ref_inv_diag(level: &LevelGraph) -> Vec<f64> {
    (0..level.n)
        .map(|u| {
            let d = level.diagonal(u);
            if d > 0.0 {
                1.0 / d
            } else {
                0.0
            }
        })
        .collect()
}

fn ref_chebyshev(
    level: &LevelGraph,
    alpha: f64,
    beta: f64,
    r: &[f64],
    d: &mut [f64],
    x: &mut [f64],
) {
    let inv = ref_inv_diag(level);
    for i in 0..r.len() {
        let dn = alpha * inv[i / 3] * r[i] + beta * d[i];
        d[i] = dn;
        x[i] += dn;
    }
}

fn ref_restrict(agg: &[u32], n_coarse: usize, fine: &[f64]) -> Vec<f64> {
    let mut coarse = vec![0.0; n_coarse * 3];
    for (i, &a) in agg.iter().enumerate() {
        for k in 0..3 {
            coarse[a as usize * 3 + k] += fine[i * 3 + k];
        }
    }
    coarse
}

fn ref_prolong_add(agg: &[u32], coarse: &[f64], fine: &mut [f64]) {
    for (i, &a) in agg.iter().enumerate() {
        for k in 0..3 {
            fine[i * 3 + k] += coarse[a as usize * 3 + k];
        }
    }
}

fn ref_dot3(a: &[f64], b: &[f64]) -> [f64; 3] {
    let mut s = [0.0; 3];
    for (i, (x, y)) in a.iter().zip(b).enumerate() {
        s[i % 3] += x * y;
    }
    s
}

fn abs_dot3(a: &[f64], b: &[f64]) -> [f64; 3] {
    let mut s = [0.0; 3];
    for (i, (x, y)) in a.iter().zip(b).enumerate() {
        s[i % 3] += (x * y).abs();
    }
    s
}

/// Coarse level of `fine` under `agg` (quotient graph): coarse edges in
/// ascending `(min(U,V), max(U,V))` order, each with the ascending list of
/// fine edges it represents; intra-aggregate fine edges vanish.
fn coarsen(fine: &LevelGraph, agg: &[u32], n_coarse: usize) -> (LevelGraph, Vec<u32>, Vec<u32>) {
    let mut pairs: Vec<((u32, u32), u32)> = Vec::new();
    for u in 0..fine.n {
        for (e, sign, v) in fine.adjacency.incident(u) {
            if sign < 0 {
                let (a, b) = (agg[u], agg[v as usize]);
                if a != b {
                    pairs.push(((a.min(b), a.max(b)), e));
                }
            }
        }
    }
    pairs.sort();
    let mut coarse_edges: Vec<(u32, u32)> = Vec::new();
    let mut offsets = vec![0u32];
    let mut fine_edges = Vec::new();
    let mut weight = Vec::new();
    for (key, e) in pairs {
        if coarse_edges.last() != Some(&key) {
            coarse_edges.push(key);
            offsets.push(fine_edges.len() as u32);
            weight.push(0.0);
        }
        fine_edges.push(e);
        *weight.last_mut().unwrap() += fine.weight[e as usize];
        *offsets.last_mut().unwrap() = fine_edges.len() as u32;
    }
    let mut anchor = vec![0.0; n_coarse];
    for (u, &a) in agg.iter().enumerate() {
        anchor[a as usize] += fine.anchor[u];
    }
    let level = LevelGraph {
        n: n_coarse,
        adjacency: csr_from_edges(n_coarse, &coarse_edges),
        weight,
        anchor,
        aggregate_of: Vec::new(),
    };
    (level, offsets, fine_edges)
}

// ── Tests ──────────────────────────────────────────────────────────────

#[test]
fn probe_reports_adapters() {
    let probe = probe_gpu();
    println!("gpu probe: {probe}");
    assert!(probe.built_with_gpu_feature);
    for a in &probe.adapters {
        println!(
            "  {} [{} {}] driver {:?} storage binding {} MiB buffer {} MiB f64 {} software {} selectable {}",
            a.name,
            a.backend,
            a.device_type,
            a.driver,
            a.max_storage_buffer_binding_size >> 20,
            a.max_buffer_size >> 20,
            a.shader_f64,
            a.software,
            a.selectable
        );
    }
    if let Some(chosen) = probe.chosen_adapter() {
        assert!(chosen.selectable);
    }
}

#[test]
fn software_adapters_need_opt_in() {
    let probe = probe_gpu();
    if probe.adapters.is_empty() || !probe.adapters.iter().all(|a| a.software) {
        eprintln!("skipping: needs a machine with only software adapters");
        return;
    }
    let strict = AdapterPolicy {
        allow_software: false,
        name_filter: None,
        ..AdapterPolicy::default()
    };
    let err = GpuContext::new(&strict).expect_err("software adapter rejected");
    let msg = err.to_string();
    assert!(msg.contains("software"), "{msg}");
    assert!(msg.contains("THESEUS_GPU_ALLOW_SOFTWARE"), "{msg}");
    let permissive = AdapterPolicy {
        allow_software: true,
        ..strict
    };
    GpuContext::new(&permissive).expect("software adapter accepted when allowed");
}

#[test]
fn adapter_name_filter_can_reject_everything() {
    let policy = AdapterPolicy {
        allow_software: true,
        name_filter: Some("no-such-adapter-name".to_string()),
        ..AdapterPolicy::default()
    };
    match GpuContext::new(&policy) {
        Err(theseus::TheseusError::GpuUnavailable(msg)) => {
            assert!(msg.contains("no-such-adapter-name"), "{msg}")
        }
        Err(e) => panic!("unexpected error {e}"),
        Ok(_) => panic!("adapter matched an impossible name"),
    }
}

#[test]
fn f64_requires_shader_f64() {
    let Some(gpu) = backend() else { return };
    let ok = gpu.try_alloc(6, Precision::F64).is_ok();
    assert_eq!(ok, gpu.shader_f64());
    assert_eq!(gpu.supports(Precision::F64), gpu.shader_f64());
    assert!(gpu.supports(Precision::F32));
}

#[test]
fn device_bytes_track_live_buffers() {
    let Some(gpu) = backend() else { return };
    let base = gpu.device_bytes();
    let a = gpu.alloc(3 * 100_000, Precision::F32);
    assert_eq!(gpu.device_bytes(), base + a.size_bytes());
    assert_eq!(a.size_bytes(), 3 * 100_000 * 4);
    let level = grid_level(&mut Rng::new(12), 32, 32);
    let dev = gpu.upload_level(&level, Precision::F32);
    assert!(gpu.device_bytes() > base + a.size_bytes());
    drop(dev);
    assert_eq!(gpu.device_bytes(), base + a.size_bytes());
    drop(a);
    assert_eq!(gpu.device_bytes(), base);
}

#[test]
fn upload_download_copy_zero() {
    let Some(gpu) = backend() else { return };
    let mut rng = Rng::new(1);
    for precision in precisions(&gpu) {
        for n in [0usize, 1, 7, 1000] {
            let src = rng.vec(n * 3, -5.0, 5.0);
            let mut a = gpu.alloc(n * 3, precision);
            assert_eq!(gpu.len(&a), n * 3);
            gpu.upload(&src, &mut a);
            let mut back = vec![0.0; n * 3];
            gpu.download(&a, &mut back);
            assert_eq!(back, round_all(precision, &src), "{precision:?} n={n}");

            let mut b = gpu.alloc(n * 3, precision);
            gpu.copy(&a, &mut b);
            gpu.download(&b, &mut back);
            assert_eq!(back, round_all(precision, &src));

            gpu.zero(&mut a);
            gpu.download(&a, &mut back);
            assert!(back.iter().all(|v| *v == 0.0));
        }
    }
}

#[test]
fn apply_graph_matches_reference() {
    let Some(gpu) = backend() else { return };
    let mut rng = Rng::new(2);
    for precision in precisions(&gpu) {
        for &(n, extra) in &[(1usize, 0usize), (2, 1), (17, 20), (200, 400), (3001, 6000)] {
            let level = rounded_level(precision, &random_level(&mut rng, n, extra));
            let x = round_all(precision, &rng.vec(n * 3, -1.0, 1.0));
            let want = ref_apply(&level, &x);

            let dev = gpu.upload_level(&level, precision);
            let mut xb = gpu.alloc(n * 3, precision);
            gpu.upload(&x, &mut xb);
            let mut yb = gpu.alloc(n * 3, precision);
            gpu.apply_graph(&dev, &xb, &mut yb);
            let mut got = vec![0.0; n * 3];
            gpu.download(&yb, &mut got);
            assert_close(
                &format!("apply_graph {precision:?} n={n}"),
                &got,
                &want,
                inf_norm(&want),
                tol(precision),
            );
        }
    }
}

#[test]
fn residual_matches_reference() {
    let Some(gpu) = backend() else { return };
    let mut rng = Rng::new(3);
    for precision in precisions(&gpu) {
        for &(n, extra) in &[(5usize, 3usize), (300, 500), (2000, 1000)] {
            let level = rounded_level(precision, &random_level(&mut rng, n, extra));
            let x = round_all(precision, &rng.vec(n * 3, -1.0, 1.0));
            let b = round_all(precision, &rng.vec(n * 3, -1.0, 1.0));
            let want = ref_residual(&level, &x, &b);

            let dev = gpu.upload_level(&level, precision);
            let mut xb = gpu.alloc(n * 3, precision);
            gpu.upload(&x, &mut xb);
            let mut bb = gpu.alloc(n * 3, precision);
            gpu.upload(&b, &mut bb);
            let mut rb = gpu.alloc(n * 3, precision);
            gpu.residual(&dev, &xb, &bb, &mut rb);
            let mut got = vec![0.0; n * 3];
            gpu.download(&rb, &mut got);
            let scale = inf_norm(&ref_apply(&level, &x)).max(inf_norm(&b));
            assert_close(
                &format!("residual {precision:?} n={n}"),
                &got,
                &want,
                scale,
                tol(precision),
            );
        }
    }
}

#[test]
fn inverse_diagonal_matches_reference() {
    let Some(gpu) = backend() else { return };
    let mut rng = Rng::new(4);
    for precision in precisions(&gpu) {
        let n = 500;
        let mut level = rounded_level(precision, &random_level(&mut rng, n, 800));
        // An isolated node without anchor has an empty diagonal → 0.
        level.n += 1;
        level
            .adjacency
            .offsets
            .push(*level.adjacency.offsets.last().unwrap());
        level.anchor.push(0.0);
        let want = ref_inv_diag(&level);
        assert_eq!(want[n], 0.0);

        let dev = gpu.upload_level(&level, precision);
        let mut got = vec![0.0; n + 1];
        gpu.download(dev.inv_diag(), &mut got);
        assert_close(
            &format!("inv_diag {precision:?}"),
            &got,
            &want,
            inf_norm(&want),
            tol(precision),
        );
    }
}

#[test]
fn chebyshev_step_matches_reference() {
    let Some(gpu) = backend() else { return };
    let mut rng = Rng::new(5);
    for precision in precisions(&gpu) {
        for &(n, extra) in &[(9usize, 4usize), (777, 900)] {
            let level = rounded_level(precision, &random_level(&mut rng, n, extra));
            let r = round_all(precision, &rng.vec(n * 3, -1.0, 1.0));
            let mut d = round_all(precision, &rng.vec(n * 3, -1.0, 1.0));
            let mut x = round_all(precision, &rng.vec(n * 3, -1.0, 1.0));
            let alpha = round(precision, 0.7321);
            let beta = round(precision, -0.4137);
            let scale = inf_norm(&x).max(inf_norm(&d)).max(inf_norm(&r));

            let dev = gpu.upload_level(&level, precision);
            let mut rb = gpu.alloc(n * 3, precision);
            gpu.upload(&r, &mut rb);
            let mut db = gpu.alloc(n * 3, precision);
            gpu.upload(&d, &mut db);
            let mut xb = gpu.alloc(n * 3, precision);
            gpu.upload(&x, &mut xb);
            // Two steps so the recurrence in d is exercised.
            gpu.chebyshev_step(&dev, alpha, beta, &rb, &mut db, &mut xb);
            gpu.chebyshev_step(&dev, beta, alpha, &rb, &mut db, &mut xb);
            ref_chebyshev(&level, alpha, beta, &r, &mut d, &mut x);
            ref_chebyshev(&level, beta, alpha, &r, &mut d, &mut x);

            let mut got = vec![0.0; n * 3];
            gpu.download(&db, &mut got);
            assert_close(
                &format!("chebyshev d {precision:?} n={n}"),
                &got,
                &d,
                scale,
                tol(precision),
            );
            gpu.download(&xb, &mut got);
            assert_close(
                &format!("chebyshev x {precision:?} n={n}"),
                &got,
                &x,
                scale,
                tol(precision),
            );
        }
    }
}

#[test]
fn restrict_and_prolong_match_reference() {
    let Some(gpu) = backend() else { return };
    let mut rng = Rng::new(6);
    for precision in precisions(&gpu) {
        for &(n_fine, n_coarse) in &[(1usize, 1usize), (10, 3), (1000, 260), (4097, 1024)] {
            let agg = random_aggregates(&mut rng, n_fine, n_coarse);
            let (offsets, members) = members_csr(&agg, n_coarse);
            assert_eq!(offsets.len(), n_coarse + 1);
            assert_eq!(members.len(), n_fine);
            for cu in 0..n_coarse {
                let m = &members[offsets[cu] as usize..offsets[cu + 1] as usize];
                assert!(m.windows(2).all(|w| w[0] < w[1]), "members sorted");
                assert!(m.iter().all(|&i| agg[i as usize] as usize == cu));
            }

            let fine = round_all(precision, &rng.vec(n_fine * 3, -1.0, 1.0));
            let want_coarse = ref_restrict(&agg, n_coarse, &fine);
            let dev = gpu.upload_aggregates(&agg, n_coarse);
            let mut fb = gpu.alloc(n_fine * 3, precision);
            gpu.upload(&fine, &mut fb);
            let mut cb = gpu.alloc(n_coarse * 3, precision);
            gpu.restrict(&dev, &fb, &mut cb);
            let mut got = vec![0.0; n_coarse * 3];
            gpu.download(&cb, &mut got);
            assert_close(
                &format!("restrict {precision:?} {n_fine}->{n_coarse}"),
                &got,
                &want_coarse,
                inf_norm(&fine) * (n_fine as f64 / n_coarse as f64).max(1.0),
                tol(precision),
            );

            let coarse = round_all(precision, &rng.vec(n_coarse * 3, -1.0, 1.0));
            let mut want_fine = fine.clone();
            ref_prolong_add(&agg, &coarse, &mut want_fine);
            gpu.upload(&coarse, &mut cb);
            gpu.prolong_add(&dev, &cb, &mut fb);
            let mut got = vec![0.0; n_fine * 3];
            gpu.download(&fb, &mut got);
            assert_close(
                &format!("prolong_add {precision:?} {n_coarse}->{n_fine}"),
                &got,
                &want_fine,
                2.0,
                tol(precision),
            );
        }
    }
}

#[test]
fn axpy_and_scale_match_reference() {
    let Some(gpu) = backend() else { return };
    let mut rng = Rng::new(7);
    for precision in precisions(&gpu) {
        for n in [1usize, 64, 1000, 70_001] {
            let x = round_all(precision, &rng.vec(n * 3, -1.0, 1.0));
            let y = round_all(precision, &rng.vec(n * 3, -1.0, 1.0));
            let alpha = [
                round(precision, 0.5),
                round(precision, -2.25),
                round(precision, 3.0),
            ];
            let gamma = [
                round(precision, 1.5),
                round(precision, 0.0),
                round(precision, -0.125),
            ];
            let mut want: Vec<f64> = y
                .iter()
                .zip(&x)
                .enumerate()
                .map(|(i, (y, x))| y + alpha[i % 3] * x)
                .collect();
            for (i, v) in want.iter_mut().enumerate() {
                *v *= gamma[i % 3];
            }

            let mut xb = gpu.alloc(n * 3, precision);
            gpu.upload(&x, &mut xb);
            let mut yb = gpu.alloc(n * 3, precision);
            gpu.upload(&y, &mut yb);
            gpu.axpy(alpha, &xb, &mut yb);
            gpu.scale(gamma, &mut yb);
            let mut got = vec![0.0; n * 3];
            gpu.download(&yb, &mut got);
            assert_close(
                &format!("axpy+scale {precision:?} n={n}"),
                &got,
                &want,
                inf_norm(&want).max(1.0),
                tol(precision),
            );
        }
    }
}

#[test]
fn dot3_and_norm3_match_reference() {
    let Some(gpu) = backend() else { return };
    let mut rng = Rng::new(8);
    for precision in precisions(&gpu) {
        // Sizes around workgroup and reduction-group boundaries, and one
        // beyond REDUCE_MAX_GROUPS × WG to exercise the grid-stride loop.
        for n in [
            0usize, 1, 63, 64, 65, 127, 128, 129, 255, 257, 4096, 5000, 300_000,
        ] {
            let a = round_all(precision, &rng.vec(n * 3, -1.0, 1.0));
            let b = round_all(precision, &rng.vec(n * 3, -1.0, 1.0));
            let mut ab = gpu.alloc(n * 3, precision);
            gpu.upload(&a, &mut ab);
            let mut bb = gpu.alloc(n * 3, precision);
            gpu.upload(&b, &mut bb);

            let want = ref_dot3(&a, &b);
            let bound = abs_dot3(&a, &b);
            let got = gpu.dot3(&ab, &bb);
            for k in 0..3 {
                assert_close(
                    &format!("dot3[{k}] {precision:?} n={n}"),
                    &[got[k]],
                    &[want[k]],
                    bound[k].max(1e-30),
                    tol(precision),
                );
            }
            // Deterministic: the same call gives the same bits.
            assert_eq!(
                gpu.dot3(&ab, &bb),
                got,
                "dot3 deterministic {precision:?} n={n}"
            );

            let want_n = ref_dot3(&a, &a).map(f64::sqrt);
            let got_n = gpu.norm3(&ab);
            for k in 0..3 {
                assert_close(
                    &format!("norm3[{k}] {precision:?} n={n}"),
                    &[got_n[k]],
                    &[want_n[k]],
                    want_n[k].max(1e-30),
                    tol(precision),
                );
            }
        }
    }
}

#[test]
fn coarse_weight_update_matches_reference() {
    let Some(gpu) = backend() else { return };
    let mut rng = Rng::new(9);
    for precision in precisions(&gpu) {
        for &(n_fine, extra, n_coarse) in &[(12usize, 10usize, 4usize), (2000, 3000, 500)] {
            let fine = rounded_level(precision, &random_level(&mut rng, n_fine, extra));
            let agg = random_aggregates(&mut rng, n_fine, n_coarse);
            let (coarse_ref, offsets, fine_edges) = coarsen(&fine, &agg, n_coarse);
            assert!(coarse_ref.num_edges() > 0);

            // Upload the coarse level with wrong weights/anchors, then let the
            // device recompute them from the fine level.
            let mut stale = coarse_ref.clone();
            stale.weight.iter_mut().for_each(|w| *w = -1.0);
            stale.anchor.iter_mut().for_each(|a| *a = -1.0);
            let fine_dev = gpu.upload_level(&fine, precision);
            let mut coarse_dev = gpu.upload_level(&stale, precision);
            let agg_dev = gpu.upload_aggregates(&agg, n_coarse);
            let map = gpu.upload_coarse_edge_map(&offsets, &fine_edges).unwrap();
            gpu.coarse_weight_update(&map, &agg_dev, &fine_dev, &mut coarse_dev);

            let mut got_w = vec![0.0; coarse_ref.num_edges()];
            gpu.download(coarse_dev.weight(), &mut got_w);
            assert_close(
                &format!("coarse weights {precision:?}"),
                &got_w,
                &coarse_ref.weight,
                inf_norm(&coarse_ref.weight),
                tol(precision),
            );
            let mut got_a = vec![0.0; n_coarse];
            gpu.download(coarse_dev.anchor(), &mut got_a);
            assert_close(
                &format!("coarse anchors {precision:?}"),
                &got_a,
                &coarse_ref.anchor,
                inf_norm(&coarse_ref.anchor).max(1e-30),
                tol(precision),
            );
            let mut got_d = vec![0.0; n_coarse];
            gpu.download(coarse_dev.inv_diag(), &mut got_d);
            let want_d = ref_inv_diag(&coarse_ref);
            assert_close(
                &format!("coarse inv_diag {precision:?}"),
                &got_d,
                &want_d,
                inf_norm(&want_d),
                tol(precision),
            );

            // The refreshed coarse operator equals Pᵀ A P applied to a vector.
            let xc = round_all(precision, &rng.vec(n_coarse * 3, -1.0, 1.0));
            let mut xf = vec![0.0; n_fine * 3];
            ref_prolong_add(&agg, &xc, &mut xf);
            let want = ref_restrict(&agg, n_coarse, &ref_apply(&fine, &xf));
            let mut xcb = gpu.alloc(n_coarse * 3, precision);
            gpu.upload(&xc, &mut xcb);
            let mut ycb = gpu.alloc(n_coarse * 3, precision);
            gpu.apply_graph(&coarse_dev, &xcb, &mut ycb);
            let mut got = vec![0.0; n_coarse * 3];
            gpu.download(&ycb, &mut got);
            assert_close(
                &format!("galerkin {precision:?}"),
                &got,
                &want,
                inf_norm(&want),
                // Pᵀ A P sums the same fine terms in a different order.
                tol(precision) * 4.0,
            );
        }
    }
}

#[test]
fn update_level_weights_refreshes_operator() {
    let Some(gpu) = backend() else { return };
    let mut rng = Rng::new(10);
    for precision in precisions(&gpu) {
        let n = 600;
        let level = rounded_level(precision, &random_level(&mut rng, n, 700));
        let mut dev = gpu.upload_level(&level, precision);
        let x = round_all(precision, &rng.vec(n * 3, -1.0, 1.0));
        let mut xb = gpu.alloc(n * 3, precision);
        gpu.upload(&x, &mut xb);
        let mut yb = gpu.alloc(n * 3, precision);
        gpu.apply_graph(&dev, &xb, &mut yb);

        // New q: weights scaled by a random factor, anchors doubled. The
        // pending apply must see the old weights, the next one the new.
        let mut fresh = level.clone();
        fresh.weight = round_all(
            precision,
            &fresh
                .weight
                .iter()
                .map(|w| w * rng.range(0.5, 1.5))
                .collect::<Vec<_>>(),
        );
        fresh.anchor = round_all(
            precision,
            &fresh.anchor.iter().map(|a| a * 2.0).collect::<Vec<_>>(),
        );
        gpu.update_level_weights(&mut dev, &fresh.weight, &fresh.anchor);
        let mut got_old = vec![0.0; n * 3];
        gpu.download(&yb, &mut got_old);
        let want_old = ref_apply(&level, &x);
        assert_close(
            &format!("apply before update {precision:?}"),
            &got_old,
            &want_old,
            inf_norm(&want_old),
            tol(precision),
        );

        gpu.apply_graph(&dev, &xb, &mut yb);
        let mut got_new = vec![0.0; n * 3];
        gpu.download(&yb, &mut got_new);
        let want_new = ref_apply(&fresh, &x);
        assert_close(
            &format!("apply after update {precision:?}"),
            &got_new,
            &want_new,
            inf_norm(&want_new),
            tol(precision),
        );
        let mut inv = vec![0.0; n];
        gpu.download(dev.inv_diag(), &mut inv);
        let want_inv = ref_inv_diag(&fresh);
        assert_close(
            &format!("inv_diag after update {precision:?}"),
            &inv,
            &want_inv,
            inf_norm(&want_inv),
            tol(precision),
        );
    }
}

#[test]
fn grid_64_apply_residual_norm_chain() {
    let Some(gpu) = backend() else { return };
    let mut rng = Rng::new(64);
    for precision in precisions(&gpu) {
        let level = rounded_level(precision, &grid_level(&mut rng, 64, 64));
        let n = level.n;
        assert_eq!(n, 4096);
        assert_eq!(level.num_edges(), 2 * 64 * 63);
        let x_true = round_all(precision, &rng.vec(n * 3, -1.0, 1.0));
        let b = round_all(precision, &ref_apply(&level, &x_true));
        let x0 = round_all(precision, &rng.vec(n * 3, -0.1, 0.1));

        let dev = gpu.upload_level(&level, precision);
        let mut xb = gpu.alloc(n * 3, precision);
        gpu.upload(&x_true, &mut xb);
        let mut bb = gpu.alloc(n * 3, precision);
        gpu.apply_graph(&dev, &xb, &mut bb); // b = A x_true on the device
        let mut x0b = gpu.alloc(n * 3, precision);
        gpu.upload(&x0, &mut x0b);
        let mut rb = gpu.alloc(n * 3, precision);
        let submits_before = gpu.stats().submits;
        gpu.residual(&dev, &x0b, &bb, &mut rb); // r = b − A x0
                                                // The upload of x0 flushed the apply; the residual stays recorded in
                                                // the pending encoder until a readback forces a submit.
        assert_eq!(gpu.pending_dispatches(), 1);
        assert_eq!(gpu.stats().submits, submits_before);
        let got_norm = gpu.norm3(&rb);
        assert!(gpu.stats().submits > submits_before);
        assert_eq!(gpu.pending_dispatches(), 0);

        let want_b = b.clone();
        let mut got_b = vec![0.0; n * 3];
        gpu.download(&bb, &mut got_b);
        assert_close(
            &format!("grid b {precision:?}"),
            &got_b,
            &want_b,
            inf_norm(&want_b),
            tol(precision),
        );
        let want_r = ref_residual(&level, &x0, &got_b);
        let mut got_r = vec![0.0; n * 3];
        gpu.download(&rb, &mut got_r);
        assert_close(
            &format!("grid r {precision:?}"),
            &got_r,
            &want_r,
            inf_norm(&want_b),
            tol(precision),
        );
        let want_norm = ref_dot3(&got_r, &got_r).map(f64::sqrt);
        for k in 0..3 {
            assert_close(
                &format!("grid ‖r‖[{k}] {precision:?}"),
                &[got_norm[k]],
                &[want_norm[k]],
                want_norm[k],
                tol(precision),
            );
        }
        let rel = got_norm
            .iter()
            .zip(ref_dot3(&b, &b).map(f64::sqrt))
            .map(|(r, b)| r / b)
            .collect::<Vec<_>>();
        println!("grid 64×64 {precision:?}: ‖b − A x0‖ / ‖b‖ = {rel:?}");
        assert!(rel.iter().all(|r| *r > 0.5 && *r < 2.0));
    }
}

#[test]
fn workgroup_sizes_agree() {
    let Some(first) = backend() else { return };
    let policy = AdapterPolicy::from_env();
    let mut rng = Rng::new(11);
    let level = rounded_level(Precision::F32, &random_level(&mut rng, 1500, 2500));
    let n = level.n;
    let x = round_all(Precision::F32, &rng.vec(n * 3, -1.0, 1.0));
    let want = ref_apply(&level, &x);
    let want_dot = ref_dot3(&x, &want);
    let bound = abs_dot3(&x, &want);
    drop(first);
    for wg in theseus::backend::gpu::WORKGROUP_SIZES {
        let ctx = GpuContext::new(&policy).expect("adapter");
        let gpu = match GpuBackend::with_workgroup_size(ctx, wg) {
            Ok(gpu) => gpu,
            Err(e) => {
                eprintln!("workgroup size {wg} unsupported: {e}");
                continue;
            }
        };
        assert_eq!(gpu.workgroup_size(), wg);
        let dev = gpu.upload_level(&level, Precision::F32);
        let mut xb = gpu.alloc(n * 3, Precision::F32);
        gpu.upload(&x, &mut xb);
        let mut yb = gpu.alloc(n * 3, Precision::F32);
        gpu.apply_graph(&dev, &xb, &mut yb);
        let mut got = vec![0.0; n * 3];
        gpu.download(&yb, &mut got);
        assert_close(
            &format!("apply wg={wg}"),
            &got,
            &want,
            inf_norm(&want),
            1e-6,
        );
        let got_dot = gpu.dot3(&xb, &yb);
        for k in 0..3 {
            assert_close(
                &format!("dot wg={wg} [{k}]"),
                &[got_dot[k]],
                &[want_dot[k]],
                bound[k],
                1e-6,
            );
        }
    }
}

// ── CSR kernels (WS-C) against `CpuBackend` ─────────────────────────────

mod csr {
    use super::*;
    use theseus::amg::hierarchy::LevelMatrix;
    use theseus::backend::cpu::CpuBackend;

    /// Random rectangular CSR with 1–4 entries per row.
    fn random_rect(rng: &mut Rng, nrows: usize, ncols: usize) -> LevelMatrix {
        let mut rows: Vec<Vec<(u32, f64)>> = (0..nrows)
            .map(|_| {
                (0..1 + rng.below(4))
                    .map(|_| (rng.below(ncols) as u32, rng.range(-1.0, 1.0)))
                    .collect()
            })
            .collect();
        LevelMatrix::from_rows(nrows, ncols, &mut rows)
    }

    /// Random square, diagonally dominant CSR (like a coarse AMG operator).
    fn random_square(rng: &mut Rng, n: usize) -> LevelMatrix {
        let mut rows: Vec<Vec<(u32, f64)>> = (0..n)
            .map(|u| {
                let mut row: Vec<(u32, f64)> = (0..rng.below(6))
                    .map(|_| (rng.below(n) as u32, -rng.range(0.1, 1.0)))
                    .collect();
                row.retain(|e| e.0 as usize != u);
                let off: f64 = row.iter().map(|e| e.1.abs()).sum();
                row.push((u as u32, off + rng.range(0.5, 1.5)));
                row
            })
            .collect();
        LevelMatrix::from_rows(n, n, &mut rows)
    }

    /// Round the stored values (and refresh the diagonal) so both backends
    /// see the same operator.
    fn rounded(precision: Precision, mut m: LevelMatrix) -> LevelMatrix {
        for v in &mut m.values {
            *v = round(precision, *v);
        }
        m.refresh_diag();
        m
    }

    fn cpu_upload(
        cpu: &CpuBackend,
        v: &[f64],
        precision: Precision,
    ) -> theseus::backend::cpu::CpuBuf {
        let mut b = cpu.alloc(v.len(), precision);
        cpu.upload(v, &mut b);
        b
    }

    #[test]
    fn apply_csr_matches_cpu_backend() {
        let Some(gpu) = backend() else { return };
        let cpu = CpuBackend::new();
        let mut rng = Rng::new(31);
        for precision in precisions(&gpu) {
            for &(n, nc) in &[(7usize, 3usize), (500, 170), (4099, 1400)] {
                let p = rounded(precision, random_rect(&mut rng, n, nc));
                let pt = p.transpose();
                let a = rounded(precision, random_square(&mut rng, n));
                let xc = round_all(precision, &rng.vec(nc * 3, -1.0, 1.0));
                let x = round_all(precision, &rng.vec(n * 3, -1.0, 1.0));
                let y0 = round_all(precision, &rng.vec(n * 3, -1.0, 1.0));

                let (gp, gpt, ga) = (
                    gpu.upload_csr(&p, precision),
                    gpu.upload_csr(&pt, precision),
                    gpu.upload_csr(&a, precision),
                );
                let (cp, cpt, ca) = (
                    cpu.upload_csr(&p, precision),
                    cpu.upload_csr(&pt, precision),
                    cpu.upload_csr(&a, precision),
                );
                assert_eq!(gp.nrows(), n);
                assert_eq!(gp.ncols(), nc);
                assert_eq!(gp.nnz(), p.nnz());
                assert_eq!(gp.precision(), precision);

                let mut gxc = gpu.alloc(nc * 3, precision);
                gpu.upload(&xc, &mut gxc);
                let mut gx = gpu.alloc(n * 3, precision);
                gpu.upload(&x, &mut gx);
                let cxc = cpu_upload(&cpu, &xc, precision);
                let cx = cpu_upload(&cpu, &x, precision);

                // Prolongation y = P x_c
                let mut gy = gpu.alloc(n * 3, precision);
                gpu.apply_csr(&gp, &gxc, &mut gy);
                let mut cy = cpu.alloc(n * 3, precision);
                cpu.apply_csr(&cp, &cxc, &mut cy);
                let mut got = vec![0.0; n * 3];
                gpu.download(&gy, &mut got);
                let want = cy.to_vec_f64();
                assert_close(
                    &format!("apply_csr P {precision:?} n={n}"),
                    &got,
                    &want,
                    inf_norm(&want),
                    tol(precision),
                );

                // y += P x_c
                let mut gy = gpu.alloc(n * 3, precision);
                gpu.upload(&y0, &mut gy);
                gpu.apply_csr_add(&gp, &gxc, &mut gy);
                let mut cy = cpu_upload(&cpu, &y0, precision);
                cpu.apply_csr_add(&cp, &cxc, &mut cy);
                gpu.download(&gy, &mut got);
                let want = cy.to_vec_f64();
                assert_close(
                    &format!("apply_csr_add P {precision:?} n={n}"),
                    &got,
                    &want,
                    inf_norm(&want),
                    tol(precision),
                );

                // Restriction r_c = Pᵀ x
                let mut grc = gpu.alloc(nc * 3, precision);
                gpu.apply_csr(&gpt, &gx, &mut grc);
                let mut crc = cpu.alloc(nc * 3, precision);
                cpu.apply_csr(&cpt, &cx, &mut crc);
                let mut got_c = vec![0.0; nc * 3];
                gpu.download(&grc, &mut got_c);
                let want = crc.to_vec_f64();
                assert_close(
                    &format!("apply_csr Pᵀ {precision:?} n={n}"),
                    &got_c,
                    &want,
                    inf_norm(&want),
                    tol(precision),
                );

                // Square operator y = A x
                let mut gya = gpu.alloc(n * 3, precision);
                gpu.apply_csr(&ga, &gx, &mut gya);
                let mut cya = cpu.alloc(n * 3, precision);
                cpu.apply_csr(&ca, &cx, &mut cya);
                gpu.download(&gya, &mut got);
                let want = cya.to_vec_f64();
                assert_close(
                    &format!("apply_csr A {precision:?} n={n}"),
                    &got,
                    &want,
                    inf_norm(&want),
                    tol(precision),
                );
            }
        }
    }

    #[test]
    fn residual_and_chebyshev_csr_match_cpu_backend() {
        let Some(gpu) = backend() else { return };
        let cpu = CpuBackend::new();
        let mut rng = Rng::new(32);
        for precision in precisions(&gpu) {
            for &n in &[11usize, 1234] {
                let a = rounded(precision, random_square(&mut rng, n));
                let x = round_all(precision, &rng.vec(n * 3, -1.0, 1.0));
                let b = round_all(precision, &rng.vec(n * 3, -1.0, 1.0));
                let d0 = round_all(precision, &rng.vec(n * 3, -1.0, 1.0));
                let alpha = round(precision, 0.7321);
                let beta = round(precision, -0.4137);

                let ga = gpu.upload_csr(&a, precision);
                let ca = cpu.upload_csr(&a, precision);
                let mut gx = gpu.alloc(n * 3, precision);
                gpu.upload(&x, &mut gx);
                let mut gb = gpu.alloc(n * 3, precision);
                gpu.upload(&b, &mut gb);
                let mut gr = gpu.alloc(n * 3, precision);
                gpu.residual_csr(&ga, &gx, &gb, &mut gr);
                let mut cx = cpu_upload(&cpu, &x, precision);
                let cb = cpu_upload(&cpu, &b, precision);
                let mut cr = cpu.alloc(n * 3, precision);
                cpu.residual_csr(&ca, &cx, &cb, &mut cr);
                let mut got = vec![0.0; n * 3];
                gpu.download(&gr, &mut got);
                let want = cr.to_vec_f64();
                assert_close(
                    &format!("residual_csr {precision:?} n={n}"),
                    &got,
                    &want,
                    inf_norm(&want).max(inf_norm(&b)),
                    tol(precision),
                );

                // Inverse diagonal as uploaded.
                let inv: Vec<f64> = a.diag.iter().map(|d| round(precision, 1.0 / d)).collect();
                let mut got_inv = vec![0.0; n];
                gpu.download(ga.inv_diag(), &mut got_inv);
                assert_close(
                    &format!("csr inv_diag {precision:?} n={n}"),
                    &got_inv,
                    &inv,
                    inf_norm(&inv),
                    tol(precision),
                );

                // Two Chebyshev steps against the CPU recurrence.
                let mut gd = gpu.alloc(n * 3, precision);
                gpu.upload(&d0, &mut gd);
                gpu.chebyshev_step_csr(&ga, alpha, beta, &gr, &mut gd, &mut gx);
                gpu.chebyshev_step_csr(&ga, beta, alpha, &gr, &mut gd, &mut gx);
                let mut cd = cpu_upload(&cpu, &d0, precision);
                cpu.chebyshev_step_csr(&ca, alpha, beta, &cr, &mut cd, &mut cx);
                cpu.chebyshev_step_csr(&ca, beta, alpha, &cr, &mut cd, &mut cx);
                let scale = inf_norm(&x).max(inf_norm(&d0)).max(inf_norm(&want));
                gpu.download(&gd, &mut got);
                assert_close(
                    &format!("chebyshev_step_csr d {precision:?} n={n}"),
                    &got,
                    &cd.to_vec_f64(),
                    scale,
                    tol(precision),
                );
                gpu.download(&gx, &mut got);
                assert_close(
                    &format!("chebyshev_step_csr x {precision:?} n={n}"),
                    &got,
                    &cx.to_vec_f64(),
                    scale,
                    tol(precision),
                );
            }
        }
    }

    #[test]
    fn update_csr_values_refreshes_operator() {
        let Some(gpu) = backend() else { return };
        let cpu = CpuBackend::new();
        let mut rng = Rng::new(33);
        for precision in precisions(&gpu) {
            let n = 900;
            let a = rounded(precision, random_square(&mut rng, n));
            let mut a2 = a.clone();
            for v in &mut a2.values {
                *v = round(precision, *v * rng.range(0.5, 1.5));
            }
            a2.refresh_diag();
            let x = round_all(precision, &rng.vec(n * 3, -1.0, 1.0));
            let b = round_all(precision, &rng.vec(n * 3, -1.0, 1.0));

            let mut ga = gpu.upload_csr(&a, precision);
            let mut gx = gpu.alloc(n * 3, precision);
            gpu.upload(&x, &mut gx);
            let mut gb = gpu.alloc(n * 3, precision);
            gpu.upload(&b, &mut gb);
            let mut gr = gpu.alloc(n * 3, precision);
            // Record a residual with the old values, then update: the pending
            // dispatch must see the old values, the next one the new ones.
            gpu.residual_csr(&ga, &gx, &gb, &mut gr);
            gpu.update_csr_values(&a2, &mut ga);
            let mut gr2 = gpu.alloc(n * 3, precision);
            gpu.residual_csr(&ga, &gx, &gb, &mut gr2);

            let ca = cpu.upload_csr(&a, precision);
            let ca2 = cpu.upload_csr(&a2, precision);
            let cx = cpu_upload(&cpu, &x, precision);
            let cb = cpu_upload(&cpu, &b, precision);
            let mut cr = cpu.alloc(n * 3, precision);
            cpu.residual_csr(&ca, &cx, &cb, &mut cr);
            let mut cr2 = cpu.alloc(n * 3, precision);
            cpu.residual_csr(&ca2, &cx, &cb, &mut cr2);

            let mut got = vec![0.0; n * 3];
            gpu.download(&gr, &mut got);
            let want = cr.to_vec_f64();
            assert_close(
                &format!("residual before update {precision:?}"),
                &got,
                &want,
                inf_norm(&want),
                tol(precision),
            );
            gpu.download(&gr2, &mut got);
            let want2 = cr2.to_vec_f64();
            assert_close(
                &format!("residual after update {precision:?}"),
                &got,
                &want2,
                inf_norm(&want2),
                tol(precision),
            );
            assert!(inf_norm(&want) > 0.0 && want != want2);

            // Chebyshev with the refreshed inverse diagonal.
            let mut gd = gpu.alloc(n * 3, precision);
            gpu.chebyshev_step_csr(&ga, 0.5, 0.0, &gr2, &mut gd, &mut gx);
            let mut cd = cpu.alloc(n * 3, precision);
            let mut cx2 = cpu_upload(&cpu, &x, precision);
            cpu.chebyshev_step_csr(&ca2, 0.5, 0.0, &cr2, &mut cd, &mut cx2);
            gpu.download(&gx, &mut got);
            let want = cx2.to_vec_f64();
            assert_close(
                &format!("chebyshev after update {precision:?}"),
                &got,
                &want,
                inf_norm(&want),
                tol(precision),
            );
        }
    }

    /// Whole hierarchy: `AmgSolver<GpuBackend>` (not wired into
    /// `LinearSolver::new` — WS-H) must reproduce the CPU solver's solution
    /// on a small grid, in every available precision.
    #[test]
    fn amg_solver_on_the_gpu_matches_the_cpu_solver() {
        use theseus::amg::AmgSolver;
        use theseus::linear_solver::{
            IterativeSolverOptions, LinearSolverKind, LinearSystemSolver, SolveRequest,
            TolerancePolicy,
        };

        let Some(gpu) = backend() else { return };
        let problem = fixtures::grid::make_recoverable_grid_problem(12);
        let topology = &problem.topology;
        let bounds = &problem.bounds;
        let q = fixtures::smooth_q_star(topology.num_edges);
        let nf = topology.free_node_indices.len();
        let mut rng = Rng::new(34);
        let rhs = rng.vec(nf * 3, -1.0, 1.0);

        for precision in precisions(&gpu) {
            let options = IterativeSolverOptions {
                tolerance: TolerancePolicy::Fixed(1e-10),
                coarsest_size: 8,
                precondition_precision: Some(precision),
                ..theseus::amg::recommended_options()
            };
            let mut cpu_solver = AmgSolver::cpu(topology, bounds, &options).unwrap();
            cpu_solver.update(&q).unwrap();
            let mut x_cpu = vec![0.0; nf * 3];
            let req = SolveRequest {
                rhs: &rhs,
                x0: None,
                tolerance: 1e-10,
                max_iterations: 200,
                cancel: None,
            };
            let stats_cpu = cpu_solver.solve(req, &mut x_cpu).unwrap();

            let gpu_ctx = GpuBackend::from_env().unwrap();
            let mut gpu_solver = AmgSolver::new(gpu_ctx, topology, bounds, &options).unwrap();
            assert_eq!(gpu_solver.kind(), LinearSolverKind::IterativeGpu);
            gpu_solver.update(&q).unwrap();
            assert_eq!(gpu_solver.level_sizes(), cpu_solver.level_sizes());
            let mut x_gpu = vec![0.0; nf * 3];
            let stats_gpu = gpu_solver.solve(req, &mut x_gpu).unwrap();
            assert!(stats_gpu.converged, "{stats_gpu:?}");
            assert_eq!(stats_gpu.backend, LinearSolverKind::IterativeGpu);
            assert_close(
                &format!("AMG solution gpu vs cpu {precision:?}"),
                &x_gpu,
                &x_cpu,
                inf_norm(&x_cpu),
                1e-8,
            );
            let diff = stats_gpu
                .iterations
                .iter()
                .zip(&stats_cpu.iterations)
                .map(|(a, b)| a.abs_diff(*b))
                .max()
                .unwrap();
            assert!(
                diff <= 2,
                "iterations gpu {:?} vs cpu {:?}",
                stats_gpu.iterations,
                stats_cpu.iterations
            );
            let mem = gpu_solver.memory_bytes();
            assert!(mem.device_bytes > 0 && mem.host_bytes > 0);
        }
    }
}

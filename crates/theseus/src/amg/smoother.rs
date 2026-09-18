//! Chebyshev smoother with Jacobi scaling and the `λ_max` power iteration
//! (program plan §3), written once over [`DeviceLevel`] so level 0 (graph
//! kernels) and the CSR levels share the code.
//!
//! The smoother applies the degree-`d` Chebyshev polynomial of `D⁻¹A` on the
//! interval `[λ_max/α, λ_max]` through the three-term recurrence of the
//! Phase-0 prototype. With `zero_initial` the first residual is `b` (one
//! operator application saved); the polynomial is the same either way, so
//! pre- and post-smoothing are adjoint and the V-cycle stays symmetric.

use crate::backend::Backend;

/// One level's operator on the device: the matrix-free graph (level 0) or
/// a CSR matrix (levels ≥ 1).
pub enum DeviceLevel<B: Backend> {
    /// Graph operator over `LevelGraph` weights and anchors.
    Graph(B::LevelGraphBuf),
    /// General sparse operator (`Pᵀ A P`).
    Csr(B::CsrBuf),
}

impl<B: Backend> DeviceLevel<B> {
    /// `y = A x`.
    pub fn apply(&self, backend: &B, x: &B::Buf, y: &mut B::Buf) {
        match self {
            Self::Graph(l) => backend.apply_graph(l, x, y),
            Self::Csr(m) => backend.apply_csr(m, x, y),
        }
    }

    /// `r = b − A x`.
    pub fn residual(&self, backend: &B, x: &B::Buf, b: &B::Buf, r: &mut B::Buf) {
        match self {
            Self::Graph(l) => backend.residual(l, x, b, r),
            Self::Csr(m) => backend.residual_csr(m, x, b, r),
        }
    }

    /// `d = alpha · D⁻¹ r + beta · d;  x += d`.
    pub fn chebyshev_step(
        &self,
        backend: &B,
        alpha: f64,
        beta: f64,
        r: &B::Buf,
        d: &mut B::Buf,
        x: &mut B::Buf,
    ) {
        match self {
            Self::Graph(l) => backend.chebyshev_step(l, alpha, beta, r, d, x),
            Self::Csr(m) => backend.chebyshev_step_csr(m, alpha, beta, r, d, x),
        }
    }
}

/// Chebyshev smoother parameters shared by all levels.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Chebyshev {
    /// Polynomial degree (`≥ 1`).
    pub degree: u32,
    /// Spectral ratio `α`: the interval is `[λ_max/α, λ_max]`.
    pub alpha: f64,
}

impl Chebyshev {
    /// Smoother of `degree` on `[λ_max/alpha, λ_max]`.
    pub fn new(degree: u32, alpha: f64) -> Self {
        assert!(degree >= 1, "Chebyshev degree must be ≥ 1");
        assert!(alpha > 1.0, "Chebyshev spectral ratio must be > 1");
        Self { degree, alpha }
    }

    /// Operator applications one sweep costs.
    pub fn applications(&self, zero_initial: bool) -> u32 {
        if zero_initial {
            self.degree - 1
        } else {
            self.degree
        }
    }

    /// One smoothing sweep on `A x ≈ b` at `level` with `lambda_max`; `r`
    /// and `d` are scratch of the level's size. With `zero_initial` the
    /// current `x` is ignored and treated as zero (and overwritten).
    #[allow(clippy::too_many_arguments)]
    pub fn smooth<B: Backend>(
        &self,
        backend: &B,
        level: &DeviceLevel<B>,
        lambda_max: f64,
        b: &B::Buf,
        x: &mut B::Buf,
        zero_initial: bool,
        r: &mut B::Buf,
        d: &mut B::Buf,
    ) {
        let lmax = lambda_max;
        let lmin = lmax / self.alpha;
        let theta = 0.5 * (lmax + lmin);
        let delta = 0.5 * (lmax - lmin);
        let sigma = theta / delta;
        let mut rho = 1.0 / sigma;
        if zero_initial {
            backend.zero(x);
            // r = b (the residual of x = 0); the D⁻¹ scaling is inside the step.
            backend.copy(b, r);
        } else {
            level.residual(backend, x, b, r);
        }
        // d = D⁻¹ r / θ;  x += d
        level.chebyshev_step(backend, 1.0 / theta, 0.0, r, d, x);
        for _ in 1..self.degree {
            level.residual(backend, x, b, r);
            let rho_new = 1.0 / (2.0 * sigma - rho);
            let c_d = rho_new * rho;
            let c_r = 2.0 * rho_new / delta;
            // d = c_r D⁻¹ r + c_d d;  x += d
            level.chebyshev_step(backend, c_r, c_d, r, d, x);
            rho = rho_new;
        }
    }
}

/// Deterministic xorshift generator for the power-iteration start vectors.
struct Xorshift(u64);

impl Xorshift {
    fn next_symmetric(&mut self) -> f64 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        let u = (x.wrapping_mul(0x2545_F491_4F6C_DD1D) >> 11) as f64 / (1u64 << 53) as f64;
        2.0 * u - 1.0
    }
}

/// Safety margin applied to the power-iteration estimate.
pub const LAMBDA_MARGIN: f64 = 1.1;

/// Relative `q` drift beyond which `λ_max` is re-estimated
/// (`max_e |Δq_e / q_e| > 0.5` since the last estimate, §3).
pub const LAMBDA_REESTIMATE_DRIFT: f64 = 0.5;

/// `λ_max(D⁻¹ A)` of `level` by `iterations` power iterations on three
/// independent deterministic start vectors (one per column), times
/// [`LAMBDA_MARGIN`]. `v`, `w`, `d`, `zero` are scratch of the level's size;
/// `host` is host scratch of the same length (`n * 3`).
///
/// The start vectors depend only on `seed` and `n`, and every kernel is
/// deterministic, so the estimate is reproducible on every thread count.
#[allow(clippy::too_many_arguments)]
pub fn estimate_lambda_max<B: Backend>(
    backend: &B,
    level: &DeviceLevel<B>,
    iterations: u32,
    seed: u64,
    host: &mut [f64],
    v: &mut B::Buf,
    w: &mut B::Buf,
    d: &mut B::Buf,
    zero: &mut B::Buf,
) -> f64 {
    let n3 = backend.len(v);
    assert_eq!(host.len(), n3, "estimate_lambda_max: host scratch length");
    if n3 == 0 {
        return 1.0;
    }
    let mut rng = Xorshift(seed.max(1) ^ 0x9E37_79B9_7F4A_7C15);
    for h in host.iter_mut() {
        *h = rng.next_symmetric();
    }
    backend.upload(host, v);
    let nv = backend.norm3(v);
    backend.scale(nv.map(|s| if s > 0.0 { 1.0 / s } else { 0.0 }), v);
    backend.zero(zero);
    let mut lambda = [1.0f64; 3];
    for _ in 0..iterations {
        // w = A v;  d = D⁻¹ w  (chebyshev_step with alpha = 1, beta = 0 onto x = 0).
        level.apply(backend, v, w);
        backend.zero(zero);
        level.chebyshev_step(backend, 1.0, 0.0, w, d, zero);
        lambda = backend.norm3(d);
        backend.copy(d, v);
        backend.scale(lambda.map(|l| if l > 0.0 { 1.0 / l } else { 0.0 }), v);
    }
    let lmax = lambda.iter().copied().fold(0.0, f64::max);
    if lmax > 0.0 && lmax.is_finite() {
        lmax * LAMBDA_MARGIN
    } else {
        // Degenerate level (zero operator): any positive value keeps the
        // recurrence finite.
        1.0
    }
}

/// Largest relative change `max_e |q_e − q_ref_e| / |q_ref_e|`; `∞` when a
/// reference entry is zero and the new one is not.
pub fn max_relative_change(q: &[f64], q_ref: &[f64]) -> f64 {
    debug_assert_eq!(q.len(), q_ref.len());
    let mut worst = 0.0f64;
    for (&a, &b) in q.iter().zip(q_ref) {
        if b != 0.0 {
            worst = worst.max(((a - b) / b).abs());
        } else if a != 0.0 {
            return f64::INFINITY;
        }
    }
    worst
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::amg::hierarchy::LevelMatrix;
    use crate::backend::CpuBackend;
    use crate::graph::{CsrAdjacency, LevelGraph};
    use crate::linear_solver::Precision;

    /// 1-D chain of `n` nodes, unit weights, both ends anchored.
    fn chain(n: usize) -> LevelGraph {
        let starts: Vec<usize> = (0..n - 1).collect();
        let ends: Vec<usize> = (1..n).collect();
        let mut anchor = vec![0.0; n];
        anchor[0] = 1.0;
        anchor[n - 1] = 1.0;
        LevelGraph {
            n,
            adjacency: CsrAdjacency::from_endpoints(n, &starts, &ends),
            weight: vec![1.0; n - 1],
            anchor,
            aggregate_of: Vec::new(),
        }
    }

    #[test]
    fn power_iteration_finds_the_top_of_the_jacobi_spectrum() {
        // D⁻¹A of the anchored chain is ½ (tridiagonal 2, −1) with unit
        // diagonal scaling: eigenvalues 1 − cos(kπ/(n+1)) → λ_max → 2.
        let n = 64;
        let g = chain(n);
        let be = CpuBackend::new();
        let level = DeviceLevel::<CpuBackend>::Graph(be.upload_level(&g, Precision::F64));
        let mut bufs: Vec<_> = (0..4).map(|_| be.alloc(n * 3, Precision::F64)).collect();
        let mut host = vec![0.0; n * 3];
        let [v, w, d, z] = bufs.as_mut_slice() else {
            unreachable!()
        };
        let lambda = estimate_lambda_max(&be, &level, 30, 7, &mut host, v, w, d, z);
        let exact = 1.0
            + (std::f64::consts::PI * n as f64 / (n as f64 + 1.0))
                .cos()
                .abs();
        assert!(lambda >= exact, "estimate {lambda} below λ_max {exact}");
        assert!(
            lambda <= 1.15 * exact,
            "estimate {lambda} too loose vs {exact}"
        );
        // Deterministic.
        let again = estimate_lambda_max(&be, &level, 30, 7, &mut host, v, w, d, z);
        assert_eq!(lambda, again);
        // Same answer through the CSR path.
        let m = LevelMatrix::from_graph(&g);
        let csr = DeviceLevel::<CpuBackend>::Csr(be.upload_csr(&m, Precision::F64));
        let lambda_csr = estimate_lambda_max(&be, &csr, 30, 7, &mut host, v, w, d, z);
        assert!((lambda_csr - lambda).abs() < 1e-12 * lambda);
    }

    #[test]
    fn chebyshev_damps_the_high_end_of_the_spectrum() {
        // Error propagation of one sweep from x = 0 on b = A e: the error
        // e − x = p(D⁻¹A) e. High-frequency modes (λ near 2) must be damped
        // strongly, the smooth mode (λ near 0) barely.
        let n = 64;
        let g = chain(n);
        let be = CpuBackend::new();
        let level = DeviceLevel::<CpuBackend>::Graph(be.upload_level(&g, Precision::F64));
        let alloc = || be.alloc(n * 3, Precision::F64);
        let (mut b, mut x, mut r, mut d) = (alloc(), alloc(), alloc(), alloc());
        let lambda = 2.0 * LAMBDA_MARGIN;
        let smoother = Chebyshev::new(2, 10.0);
        let mode = |k: usize| -> Vec<f64> {
            (0..n)
                .flat_map(|i| {
                    let v = (std::f64::consts::PI * k as f64 * (i as f64 + 1.0) / (n as f64 + 1.0))
                        .sin();
                    [v; 3]
                })
                .collect()
        };
        let mut damping = |k: usize| -> f64 {
            let e = mode(k);
            let mut ae = vec![0.0; n * 3];
            let eb = crate::graph::BlockVec::<3>::from_vec(e.clone());
            let mut aeb = crate::graph::BlockVec::<3>::zeros(n);
            g.apply(&eb, &mut aeb);
            ae.copy_from_slice(aeb.as_slice());
            be.upload(&ae, &mut b);
            smoother.smooth(&be, &level, lambda, &b, &mut x, true, &mut r, &mut d);
            let mut xh = vec![0.0; n * 3];
            be.download(&x, &mut xh);
            let err: f64 = e
                .iter()
                .zip(&xh)
                .map(|(a, b)| (a - b).powi(2))
                .sum::<f64>()
                .sqrt();
            let norm: f64 = e.iter().map(|v| v * v).sum::<f64>().sqrt();
            err / norm
        };
        let high = damping(n - 1);
        let mid = damping(n / 2);
        let low = damping(1);
        // Degree 2 on [λ/10, λ]: |p| ≤ 1/T₂(σ) ≈ 0.5 on the interval, and
        // the top mode sits near a root.
        assert!(high < 0.2, "highest mode damped by only {high}");
        assert!(mid < 0.5, "mid mode damped by only {mid}");
        assert!(
            low > 0.9,
            "smooth mode should pass almost unchanged, got {low}"
        );
        // Pre- (zero initial) and non-zero-initial sweeps apply the same
        // polynomial: smoothing from x = 0 explicitly must give the same x.
        let e = mode(5);
        let eb = crate::graph::BlockVec::<3>::from_vec(e);
        let mut aeb = crate::graph::BlockVec::<3>::zeros(n);
        g.apply(&eb, &mut aeb);
        be.upload(aeb.as_slice(), &mut b);
        smoother.smooth(&be, &level, lambda, &b, &mut x, true, &mut r, &mut d);
        let mut x1 = vec![0.0; n * 3];
        be.download(&x, &mut x1);
        be.zero(&mut x);
        smoother.smooth(&be, &level, lambda, &b, &mut x, false, &mut r, &mut d);
        let mut x2 = vec![0.0; n * 3];
        be.download(&x, &mut x2);
        for (a, b) in x1.iter().zip(&x2) {
            assert!((a - b).abs() < 1e-13, "{a} vs {b}");
        }
    }

    #[test]
    fn relative_change_handles_zero_reference() {
        assert_eq!(max_relative_change(&[1.0, 2.0], &[1.0, 1.0]), 1.0);
        assert_eq!(max_relative_change(&[0.0, 2.0], &[0.0, 1.0]), 1.0);
        assert_eq!(max_relative_change(&[1.0], &[0.0]), f64::INFINITY);
        assert_eq!(max_relative_change(&[], &[]), 0.0);
    }
}

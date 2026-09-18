//! Preconditioned conjugate gradients on the block of three right-hand
//! sides (program plan §3 "Outer"): every column carries its own scalars and
//! stops at `‖r_k‖ ≤ tol · ‖b_k‖` (`‖b_k‖` floored at `1e-300`; a column
//! with `b_k = 0` returns `x_k = 0` in zero iterations). Converged columns
//! are frozen (`α_k = 0`) while the others continue. With `flexible` the
//! recurrence is FCG(1) (`β` from `⟨z, Ap⟩`, `α` from `⟨p, r⟩`), which the
//! K-cycle needs because it is not a fixed operator.
//!
//! Everything runs on backend buffers; only per-column scalars come back to
//! the host, so the loop is allocation-free and its results are as
//! deterministic as the backend's reductions.

use super::smoother::DeviceLevel;
use crate::backend::Backend;
use crate::types::TheseusError;
use std::sync::atomic::{AtomicBool, Ordering};

/// Floor applied to `‖b‖` in the stopping test.
pub const RHS_NORM_FLOOR: f64 = 1e-300;

/// Vectors of the outer iteration (all `n * 3`, `f64`).
pub struct PcgBuffers<Buf> {
    pub r: Buf,
    pub z: Buf,
    pub p: Buf,
    pub ap: Buf,
}

impl<Buf> PcgBuffers<Buf> {
    /// Four zeroed `n * 3` vectors in `precision`.
    pub fn new<B: Backend<Buf = Buf>>(
        backend: &B,
        n: usize,
        precision: crate::linear_solver::Precision,
    ) -> Self {
        Self {
            r: backend.alloc(n * 3, precision),
            z: backend.alloc(n * 3, precision),
            p: backend.alloc(n * 3, precision),
            ap: backend.alloc(n * 3, precision),
        }
    }
}

/// What one PCG run did.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct PcgOutcome {
    /// Iterations until each column converged (or the budget ran out).
    pub iterations: [u32; 3],
    /// True relative residual `‖b − A x‖ / max(‖b‖, 1e-300)` per column,
    /// recomputed from scratch at the end.
    pub relative_residual: [f64; 3],
    /// Every column met `‖r‖ ≤ tol ‖b‖` (recurrence residual).
    pub converged: bool,
}

fn safe_div(a: f64, b: f64) -> f64 {
    if b == 0.0 {
        0.0
    } else {
        a / b
    }
}

fn div3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        safe_div(a[0], b[0]),
        safe_div(a[1], b[1]),
        safe_div(a[2], b[2]),
    ]
}

/// Solve `A x = b` on three columns from the initial guess in `x`, with
/// `precondition(r, z)` computing `z = M r`.
///
/// Returns [`TheseusError::Cancelled`] as soon as `cancel` is observed set
/// (checked once per iteration, before the preconditioner); `x` then holds
/// the last iterate. Non-convergence is *not* an error here: the caller
/// decides from `PcgOutcome::converged`.
#[allow(clippy::too_many_arguments)]
pub fn pcg<B: Backend>(
    backend: &B,
    a: &DeviceLevel<B>,
    b: &B::Buf,
    x: &mut B::Buf,
    work: &mut PcgBuffers<B::Buf>,
    tol: f64,
    max_iterations: u32,
    flexible: bool,
    cancel: Option<&AtomicBool>,
    mut precondition: impl FnMut(&B::Buf, &mut B::Buf),
) -> Result<PcgOutcome, TheseusError> {
    let PcgBuffers { r, z, p, ap } = work;
    let b_norm = backend.norm3(b);
    let mut converged = [false; 3];
    let mut iterations = [0u32; 3];
    // Columns with b = 0: x = 0, done.
    let mut keep = [1.0; 3];
    for k in 0..3 {
        if b_norm[k] == 0.0 {
            keep[k] = 0.0;
            converged[k] = true;
        }
    }
    if keep != [1.0; 3] {
        backend.scale(keep, x);
    }
    a.residual(backend, x, b, r);
    let floor = b_norm.map(|v| v.max(RHS_NORM_FLOOR));
    let mut rz_prev = [0.0; 3];
    let mut it = 0u32;
    loop {
        let rn = backend.norm3(r);
        for k in 0..3 {
            if !converged[k] && rn[k] <= tol * floor[k] {
                converged[k] = true;
                iterations[k] = it;
            }
        }
        if converged.iter().all(|&c| c) || it >= max_iterations {
            break;
        }
        if cancel.is_some_and(|flag| flag.load(Ordering::Acquire)) {
            return Err(TheseusError::Cancelled);
        }
        precondition(r, z);
        let rz = backend.dot3(r, z);
        if it == 0 {
            backend.copy(z, p);
        } else {
            let beta = if flexible {
                div3(backend.dot3(z, ap), backend.dot3(p, ap)).map(|v| -v)
            } else {
                div3(rz, rz_prev)
            };
            // p = z + beta p
            backend.scale(beta, p);
            backend.axpy([1.0; 3], z, p);
        }
        rz_prev = rz;
        a.apply(backend, p, ap);
        let pap = backend.dot3(p, ap);
        let num = if flexible { backend.dot3(p, r) } else { rz };
        let mut alpha = div3(num, pap);
        for k in 0..3 {
            if converged[k] {
                alpha[k] = 0.0;
            }
        }
        backend.axpy(alpha, p, x);
        backend.axpy(alpha.map(|v| -v), ap, r);
        it += 1;
    }
    for k in 0..3 {
        if !converged[k] {
            iterations[k] = it;
        }
    }
    // True residual, recomputed from scratch.
    a.residual(backend, x, b, r);
    let rn = backend.norm3(r);
    Ok(PcgOutcome {
        iterations,
        relative_residual: [
            safe_div(rn[0], floor[0]),
            safe_div(rn[1], floor[1]),
            safe_div(rn[2], floor[2]),
        ],
        converged: converged.iter().all(|&c| c),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::amg::hierarchy::LevelMatrix;
    use crate::backend::CpuBackend;
    use crate::graph::{CsrAdjacency, LevelGraph};
    use crate::linear_solver::Precision;

    fn chain(n: usize) -> LevelGraph {
        let starts: Vec<usize> = (0..n - 1).collect();
        let ends: Vec<usize> = (1..n).collect();
        let mut anchor = vec![0.0; n];
        anchor[0] = 1.0;
        LevelGraph {
            n,
            adjacency: CsrAdjacency::from_endpoints(n, &starts, &ends),
            weight: (0..n - 1).map(|k| 1.0 + (k % 3) as f64).collect(),
            anchor,
            aggregate_of: Vec::new(),
        }
    }

    #[test]
    fn jacobi_preconditioned_cg_solves_a_chain_per_column() {
        let n = 50;
        let g = chain(n);
        let be = CpuBackend::new();
        let level = DeviceLevel::<CpuBackend>::Graph(be.upload_level(&g, Precision::F64));
        let m = LevelMatrix::from_graph(&g);
        let mut b = be.alloc(n * 3, Precision::F64);
        let mut x = be.alloc(n * 3, Precision::F64);
        // Column 1 has b = 0 → x = 0; column 2 starts from a warm guess.
        let bh: Vec<f64> = (0..n * 3)
            .map(|i| {
                if i % 3 == 1 {
                    0.0
                } else {
                    ((i % 5) as f64) - 2.0
                }
            })
            .collect();
        be.upload(&bh, &mut b);
        let xh: Vec<f64> = (0..n * 3)
            .map(|i| if i % 3 == 0 { 0.0 } else { 0.3 })
            .collect();
        be.upload(&xh, &mut x);
        let mut work = PcgBuffers::new(&be, n, Precision::F64);
        let mut zero = be.alloc(n * 3, Precision::F64);
        let mut d = be.alloc(n * 3, Precision::F64);
        let out = pcg(
            &be,
            &level,
            &b,
            &mut x,
            &mut work,
            1e-12,
            500,
            false,
            None,
            |r, z| {
                be.zero(&mut zero);
                be.zero(&mut d);
                // z = D⁻¹ r via a Chebyshev step onto zero.
                level.chebyshev_step(&be, 1.0, 0.0, r, &mut d, &mut zero);
                be.copy(&d, z);
            },
        )
        .unwrap();
        assert!(out.converged, "{out:?}");
        assert_eq!(out.iterations[1], 0);
        assert!(out.iterations[0] > 0 && out.iterations[0] <= n as u32 + 2);
        assert!(out.relative_residual.iter().all(|&r| r <= 1e-11), "{out:?}");
        let mut xh = vec![0.0; n * 3];
        be.download(&x, &mut xh);
        let mut ax = vec![0.0; n * 3];
        m.apply(&xh, &mut ax, 3);
        for (i, (l, r)) in ax.iter().zip(&bh).enumerate() {
            assert!((l - r).abs() < 1e-9, "entry {i}: {l} vs {r}");
        }
        assert!(xh.iter().skip(1).step_by(3).all(|&v| v == 0.0));
    }

    #[test]
    fn budget_and_cancellation() {
        let n = 30;
        let g = chain(n);
        let be = CpuBackend::new();
        let level = DeviceLevel::<CpuBackend>::Graph(be.upload_level(&g, Precision::F64));
        let mut b = be.alloc(n * 3, Precision::F64);
        be.upload(&vec![1.0; n * 3], &mut b);
        let mut x = be.alloc(n * 3, Precision::F64);
        let mut work = PcgBuffers::new(&be, n, Precision::F64);
        let identity = |r: &_, z: &mut _| be.copy(r, z);
        let out = pcg(
            &be, &level, &b, &mut x, &mut work, 1e-12, 2, false, None, identity,
        )
        .unwrap();
        assert!(!out.converged);
        assert_eq!(out.iterations, [2; 3]);
        assert!(out
            .relative_residual
            .iter()
            .all(|&r| r > 1e-12 && r.is_finite()));

        let cancel = AtomicBool::new(true);
        let err = pcg(
            &be,
            &level,
            &b,
            &mut x,
            &mut work,
            1e-12,
            10,
            false,
            Some(&cancel),
            identity,
        );
        assert!(matches!(err, Err(TheseusError::Cancelled)));
    }
}

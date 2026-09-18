//! Multigrid cycles over the device hierarchy: the V-cycle (one pre- and one
//! post-smoothing sweep, the fixed SPD preconditioner of §3) and Notay's
//! K-cycle (`CycleKind::K`, kept for experiments; ported from the Phase-0
//! prototype).
//!
//! [`Cycle::apply`] computes `work[0].x ≈ A_0⁻¹ work[0].b`. Transfers use
//! the smoothed `P` / `Pᵀ` as CSR products (`apply_csr`), never the
//! piecewise-constant aggregate kernels.

use super::coarsest::CoarsestSolver;
use super::smoother::{Chebyshev, DeviceLevel};
use crate::backend::Backend;
use crate::linear_solver::CycleKind;

/// Per-level vectors of one cycle application (all `n_l * 3`). The K-cycle
/// needs the four extra vectors; they are empty (`len 0`) for the V-cycle.
pub struct LevelWork<Buf> {
    pub b: Buf,
    pub x: Buf,
    pub r: Buf,
    pub d: Buf,
    pub c1: Buf,
    pub v1: Buf,
    pub v2: Buf,
    pub rt: Buf,
    /// Coarsest level only: host copies of `b` and `x` (`n * 3` each) and
    /// the triangular-solve workspace (`6 n`).
    pub host_b: Vec<f64>,
    pub host_x: Vec<f64>,
    pub solve_work: Vec<f64>,
}

impl<Buf> LevelWork<Buf> {
    /// Buffers for a level of `n` nodes; `kcycle` adds the K-cycle scratch.
    pub fn new<B: Backend<Buf = Buf>>(
        backend: &B,
        n: usize,
        precision: crate::linear_solver::Precision,
        kcycle: bool,
    ) -> Self {
        let alloc = |len: usize| backend.alloc(len, precision);
        let extra = if kcycle { n * 3 } else { 0 };
        Self {
            b: alloc(n * 3),
            x: alloc(n * 3),
            r: alloc(n * 3),
            d: alloc(n * 3),
            c1: alloc(extra),
            v1: alloc(extra),
            v2: alloc(extra),
            rt: alloc(extra),
            host_b: Vec::new(),
            host_x: Vec::new(),
            solve_work: Vec::new(),
        }
    }

    /// Give the coarsest level its host-side solve scratch.
    pub fn with_coarsest_scratch(mut self, n: usize) -> Self {
        self.host_b = vec![0.0; n * 3];
        self.host_x = vec![0.0; n * 3];
        self.solve_work = vec![0.0; 6 * n];
        self
    }
}

/// Smoothed transfer of one level: `P` (`n_l × n_{l+1}`) and `Pᵀ`.
pub struct Transfer<B: Backend> {
    pub p: B::CsrBuf,
    pub pt: B::CsrBuf,
}

/// Borrowed view of the device hierarchy for one preconditioner
/// application.
pub struct Cycle<'a, B: Backend> {
    pub backend: &'a B,
    /// `levels[l]` is `A_l`; the last one is the coarsest.
    pub levels: &'a [DeviceLevel<B>],
    /// `transfers[l]` maps between levels `l` and `l + 1`.
    pub transfers: &'a [Transfer<B>],
    /// `λ_max(D⁻¹ A_l)` per level.
    pub lambda_max: &'a [f64],
    pub smoother: Chebyshev,
    pub kind: CycleKind,
    /// K-cycle: skip the second inner step when `‖r̃‖ ≤ ktol ‖r‖`.
    pub ktol: f64,
    pub coarsest: &'a CoarsestSolver,
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

impl<B: Backend> Cycle<'_, B> {
    /// Number of levels.
    pub fn num_levels(&self) -> usize {
        self.levels.len()
    }

    /// One cycle from level 0: `work[0].x ≈ A_0⁻¹ work[0].b`.
    pub fn apply(&self, work: &mut [LevelWork<B::Buf>]) {
        self.cycle(0, work);
    }

    fn coarsest_solve(&self, w: &mut LevelWork<B::Buf>) {
        let be = self.backend;
        be.download(&w.b, &mut w.host_b);
        self.coarsest
            .solve(&w.host_b, &mut w.host_x, &mut w.solve_work);
        be.upload(&w.host_x, &mut w.x);
    }

    /// One multigrid cycle at level `l` on `work[0]` (`work` starts at `l`).
    fn cycle(&self, l: usize, work: &mut [LevelWork<B::Buf>]) {
        let be = self.backend;
        let (cur, rest) = work.split_first_mut().expect("level work present");
        if l + 1 == self.levels.len() {
            self.coarsest_solve(cur);
            return;
        }
        let level = &self.levels[l];
        let lambda = self.lambda_max[l];
        self.smoother.smooth(
            be, level, lambda, &cur.b, &mut cur.x, true, &mut cur.r, &mut cur.d,
        );
        level.residual(be, &cur.x, &cur.b, &mut cur.r);
        let t = &self.transfers[l];
        be.apply_csr(&t.pt, &cur.r, &mut rest[0].b);
        match self.kind {
            CycleKind::V => self.cycle(l + 1, rest),
            CycleKind::K => self.kcycle_solve(l + 1, rest),
        }
        be.apply_csr_add(&t.p, &rest[0].x, &mut cur.x);
        self.smoother.smooth(
            be, level, lambda, &cur.b, &mut cur.x, false, &mut cur.r, &mut cur.d,
        );
    }

    /// Notay's K-cycle coarse solve at level `l`: at most two flexible-CG
    /// steps on `A_l x = b` preconditioned by `cycle(l)`, the second one
    /// skipped when the first reduces the residual below `ktol`.
    fn kcycle_solve(&self, l: usize, work: &mut [LevelWork<B::Buf>]) {
        let be = self.backend;
        if l + 1 == self.levels.len() {
            self.cycle(l, work);
            return;
        }
        let level = &self.levels[l];
        let b_norm = be.norm3(&work[0].b);
        self.cycle(l, work);
        let (rho1, alpha1, a);
        {
            let w = &mut work[0];
            be.copy(&w.x, &mut w.c1);
            level.apply(be, &w.c1, &mut w.v1);
            rho1 = be.dot3(&w.c1, &w.v1);
            alpha1 = be.dot3(&w.c1, &w.b);
            a = div3(alpha1, rho1);
            // rt = b − a v1
            be.copy(&w.b, &mut w.rt);
            be.axpy(a.map(|v| -v), &w.v1, &mut w.rt);
            let rt_norm = be.norm3(&w.rt);
            if (0..3).all(|k| rt_norm[k] <= self.ktol * b_norm[k]) {
                // x = a c1
                be.copy(&w.c1, &mut w.x);
                be.scale(a, &mut w.x);
                return;
            }
            be.copy(&w.rt, &mut w.b);
        }
        self.cycle(l, work);
        let w = &mut work[0];
        level.apply(be, &w.x, &mut w.v2);
        let gamma = be.dot3(&w.x, &w.v1);
        let beta = be.dot3(&w.x, &w.rt);
        let c2v2 = be.dot3(&w.x, &w.v2);
        let mut coef1 = [0.0; 3];
        let mut coef2 = [0.0; 3];
        for k in 0..3 {
            let rho2 = c2v2[k] - safe_div(gamma[k] * gamma[k], rho1[k]);
            coef1[k] = safe_div(alpha1[k], rho1[k]) - safe_div(gamma[k] * beta[k], rho1[k] * rho2);
            coef2[k] = safe_div(beta[k], rho2);
        }
        // x = coef1 c1 + coef2 x
        be.scale(coef2, &mut w.x);
        be.axpy(coef1, &w.c1, &mut w.x);
    }
}

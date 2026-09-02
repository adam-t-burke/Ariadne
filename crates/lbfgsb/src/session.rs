// SPDX-License-Identifier: BSD-3-Clause
//! Safe scalar translation of the compact L-BFGS-B 3.0 iteration blocks.
//!
//! Upstream: <http://users.iems.northwestern.edu/~nocedal/lbfgsb.html>.
//! The upstream distribution's verbatim license and attribution are preserved
//! in `UPSTREAM_LICENSE.txt` and `THIRD_PARTY_NOTICES.md`; no upstream
//! copyright holder is inferred from its blank license template.
#![allow(clippy::needless_range_loop, clippy::too_many_arguments)]

use crate::{
    kernel::{History, Kernel},
    Workspace,
};
#[cfg(feature = "benchmark-instrumentation")]
use std::time::Instant;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum FormkError {
    FirstFactorization,
    TriangularSolve,
    SecondFactorization,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum DirectionError {
    NoDescent,
    Factorization,
    Numerical,
}

#[derive(Debug)]
pub(crate) struct Session {
    n: usize,
    m: usize,
    head: usize,
    col: usize,
    updates: usize,
    nfree: usize,
    theta: f64,
    updated: bool,
    constrained: bool,
    accepted_updates: usize,
    skipped_updates: usize,
    history_resets: usize,
    kernel: Kernel,
    #[cfg(feature = "benchmark-instrumentation")]
    timings: SessionTimings,
}

#[cfg(feature = "benchmark-instrumentation")]
#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct SessionTimings {
    pub(crate) cauchy_nanoseconds: u64,
    pub(crate) freev_nanoseconds: u64,
    pub(crate) formk_nanoseconds: u64,
    pub(crate) cmprlb_nanoseconds: u64,
    pub(crate) subsm_nanoseconds: u64,
}

impl Session {
    pub(crate) fn new(n: usize, m: usize, kernel: Kernel) -> Self {
        Self {
            n,
            m,
            head: 0,
            col: 0,
            updates: 0,
            nfree: n,
            theta: 1.0,
            updated: false,
            constrained: false,
            accepted_updates: 0,
            skipped_updates: 0,
            history_resets: 0,
            kernel,
            #[cfg(feature = "benchmark-instrumentation")]
            timings: SessionTimings::default(),
        }
    }

    pub(crate) fn reset_history(&mut self) {
        self.history_resets += 1;
        self.head = 0;
        self.col = 0;
        self.updates = 0;
        self.theta = 1.0;
        self.updated = false;
    }

    pub(crate) fn has_history(&self) -> bool {
        self.col != 0
    }

    pub(crate) fn accepted_updates(&self) -> usize {
        self.accepted_updates
    }

    pub(crate) fn skipped_updates(&self) -> usize {
        self.skipped_updates
    }

    pub(crate) fn history_resets(&self) -> usize {
        self.history_resets
    }

    #[cfg(feature = "benchmark-instrumentation")]
    pub(crate) fn benchmark_timings(&self) -> SessionTimings {
        self.timings
    }

    pub(crate) fn direction(&mut self, x: &[f64], w: &mut Workspace) -> Result<(), DirectionError> {
        self.constrained = w
            .lower
            .iter()
            .zip(&w.upper)
            .any(|(&l, &u)| l.is_finite() || u.is_finite());
        #[cfg(feature = "benchmark-instrumentation")]
        let phase_start = Instant::now();
        let cauchy_ok = cauchy(self, x, w);
        #[cfg(feature = "benchmark-instrumentation")]
        {
            self.timings.cauchy_nanoseconds += elapsed_nanoseconds(phase_start);
        }
        if !cauchy_ok {
            if self.col > 0 {
                self.skipped_updates += 1;
            }
            self.reset_history();
            #[cfg(feature = "benchmark-instrumentation")]
            let phase_start = Instant::now();
            let cauchy_ok = cauchy(self, x, w);
            #[cfg(feature = "benchmark-instrumentation")]
            {
                self.timings.cauchy_nanoseconds += elapsed_nanoseconds(phase_start);
            }
            if !cauchy_ok {
                return Err(DirectionError::Factorization);
            }
        }
        #[cfg(feature = "benchmark-instrumentation")]
        let phase_start = Instant::now();
        let (nenter, ileave, work) = freev(self, w);
        #[cfg(feature = "benchmark-instrumentation")]
        {
            self.timings.freev_nanoseconds += elapsed_nanoseconds(phase_start);
        }
        if self.nfree != 0 && self.col != 0 {
            #[cfg(feature = "benchmark-instrumentation")]
            let phase_start = Instant::now();
            let formk_result = if work {
                formk(self, nenter, ileave, w)
            } else {
                Ok(())
            };
            #[cfg(feature = "benchmark-instrumentation")]
            {
                self.timings.formk_nanoseconds += elapsed_nanoseconds(phase_start);
            }
            if formk_result.is_err() {
                self.skipped_updates += 1;
                self.reset_history();
                return self.direction(x, w);
            }
            #[cfg(feature = "benchmark-instrumentation")]
            let phase_start = Instant::now();
            let cmprlb_ok = cmprlb(self, x, w);
            #[cfg(feature = "benchmark-instrumentation")]
            {
                self.timings.cmprlb_nanoseconds += elapsed_nanoseconds(phase_start);
            }
            #[cfg(feature = "benchmark-instrumentation")]
            let phase_start = Instant::now();
            let subsm_ok = cmprlb_ok && subsm(self, x, w);
            #[cfg(feature = "benchmark-instrumentation")]
            {
                self.timings.subsm_nanoseconds += elapsed_nanoseconds(phase_start);
            }
            if !cmprlb_ok || !subsm_ok {
                self.skipped_updates += 1;
                self.reset_history();
                return self.direction(x, w);
            }
        }
        for i in 0..self.n {
            w.direction[i] = w.cauchy[i] - x[i];
        }
        let directional_derivative = dot(&w.evaluation_gradient, &w.direction);
        if !directional_derivative.is_finite() {
            Err(DirectionError::Numerical)
        } else if directional_derivative < 0.0 {
            Ok(())
        } else {
            Err(DirectionError::NoDescent)
        }
    }

    pub(crate) fn update(
        &mut self,
        x: &[f64],
        old_directional_derivative: f64,
        w: &mut Workspace,
    ) -> bool {
        let mut dr = 0.0;
        let mut rr = 0.0;
        let mut gs = 0.0;
        for i in 0..self.n {
            w.displacement[i] = x[i] - w.old_x[i];
            w.product[i] = w.evaluation_gradient[i] - w.old_gradient[i];
            dr += w.displacement[i] * w.product[i];
            rr += w.product[i] * w.product[i];
            gs += w.old_gradient[i] * w.displacement[i];
        }
        if dr <= f64::EPSILON * (-old_directional_derivative).max(-gs) {
            self.updated = false;
            self.skipped_updates += 1;
            return false;
        }
        self.updates += 1;
        let slot = if self.updates <= self.m {
            self.col = self.updates;
            (self.head + self.updates - 1) % self.m
        } else {
            let slot = self.head;
            self.head = (self.head + 1) % self.m;
            slot
        };
        let base = slot * self.n;
        w.s[base..base + self.n].copy_from_slice(&w.displacement);
        w.y[base..base + self.n].copy_from_slice(&w.product);
        self.theta = rr / dr;
        matupd(self, w);
        self.updated = formt(self, w);
        if !self.updated {
            self.skipped_updates += 1;
            self.reset_history();
        } else {
            self.accepted_updates += 1;
        }
        self.updated
    }
}

fn at(a: &[f64], ld: usize, row: usize, col: usize) -> f64 {
    a[col * ld + row]
}

fn set(a: &mut [f64], ld: usize, row: usize, col: usize, value: f64) {
    a[col * ld + row] = value;
}

fn slot(s: &Session, logical: usize) -> usize {
    (s.head + logical) % s.m
}

fn vec_at(a: &[f64], n: usize, physical: usize, row: usize) -> f64 {
    a[physical * n + row]
}

fn matupd(s: &Session, w: &mut Workspace) {
    if s.updates > s.m {
        for j in 0..s.col - 1 {
            for i in j..s.col - 1 {
                let sy = at(&w.sy, s.m, i + 1, j + 1);
                set(&mut w.sy, s.m, i, j, sy);
                let ss = at(&w.ss, s.m, j + 1, i + 1);
                set(&mut w.ss, s.m, j, i, ss);
            }
        }
    }
    let newest = s.col - 1;
    let (sy_products, rest) = w.wa.split_at_mut(s.m);
    let ss_products = &mut rest[..s.m];
    s.kernel.newest_products(
        History {
            n: s.n,
            m: s.m,
            head: s.head,
            col: s.col,
            s: &w.s,
            y: &w.y,
        },
        &mut sy_products[..s.col],
        &mut ss_products[..s.col],
    );
    for j in 0..s.col {
        set(&mut w.sy, s.m, newest, j, sy_products[j]);
        set(&mut w.ss, s.m, j, newest, ss_products[j]);
    }
}

fn formt(s: &Session, w: &mut Workspace) -> bool {
    for j in 0..s.col {
        set(&mut w.wt, s.m, 0, j, s.theta * at(&w.ss, s.m, 0, j));
    }
    for i in 1..s.col {
        for j in i..s.col {
            let mut value = 0.0;
            for k in 0..i {
                value += at(&w.sy, s.m, i, k) * at(&w.sy, s.m, j, k) / at(&w.sy, s.m, k, k);
            }
            value += s.theta * at(&w.ss, s.m, i, j);
            set(&mut w.wt, s.m, i, j, value);
        }
    }
    cholesky_upper(&mut w.wt, s.m, 0, s.col)
}

fn bmv(s: &Session, sy: &[f64], wt: &[f64], v: &[f64], p: &mut [f64]) -> bool {
    if s.col == 0 {
        return true;
    }
    p[s.col] = v[s.col];
    for i in 1..s.col {
        let mut sum = 0.0;
        for k in 0..i {
            sum += at(sy, s.m, i, k) * v[k] / at(sy, s.m, k, k);
        }
        p[s.col + i] = v[s.col + i] + sum;
    }
    if !solve_upper(wt, s.m, 0, s.col, &mut p[s.col..2 * s.col], true) {
        return false;
    }
    for i in 0..s.col {
        p[i] = v[i] / at(sy, s.m, i, i).sqrt();
    }
    if !solve_upper(wt, s.m, 0, s.col, &mut p[s.col..2 * s.col], false) {
        return false;
    }
    for i in 0..s.col {
        p[i] = -p[i] / at(sy, s.m, i, i).sqrt();
    }
    for i in 0..s.col {
        let mut sum = 0.0;
        for k in i + 1..s.col {
            sum += at(sy, s.m, k, i) * p[s.col + k] / at(sy, s.m, i, i);
        }
        p[i] += sum;
    }
    true
}

fn cauchy(s: &Session, x: &[f64], w: &mut Workspace) -> bool {
    let c2 = 2 * s.col;
    let mut nbreak = 0;
    let mut nfree_tail = s.n;
    let mut ibkmin = 0;
    let mut bkmin = 0.0;
    let mut f1 = 0.0;
    let mut bounded = true;
    for i in 0..s.n {
        let neg = -w.evaluation_gradient[i];
        if w.lower[i] == w.upper[i] {
            w.status[i] = 3;
        } else if !(w.lower[i].is_finite() || w.upper[i].is_finite()) {
            w.status[i] = -1;
        } else {
            w.status[i] = 0;
            if x[i] <= w.lower[i] && neg <= 0.0 {
                w.status[i] = 1;
            } else if x[i] >= w.upper[i] && neg >= 0.0 {
                w.status[i] = 2;
            } else if neg == 0.0 {
                w.status[i] = -3;
            }
        }
        if w.status[i] != 0 && w.status[i] != -1 {
            w.direction[i] = 0.0;
            continue;
        }
        w.direction[i] = neg;
        f1 -= neg * neg;
        let bp = if neg < 0.0 && w.lower[i].is_finite() {
            Some((x[i] - w.lower[i]) / -neg)
        } else if neg > 0.0 && w.upper[i].is_finite() {
            Some((w.upper[i] - x[i]) / neg)
        } else {
            None
        };
        if let Some(t) = bp {
            w.order[nbreak] = i;
            w.breakpoints[nbreak] = t;
            if nbreak == 0 || t < bkmin {
                bkmin = t;
                ibkmin = nbreak;
            }
            nbreak += 1;
        } else {
            nfree_tail -= 1;
            w.order[nfree_tail] = i;
            if neg != 0.0 {
                bounded = false;
            }
        }
    }
    for j in 0..s.col {
        let p = slot(s, j);
        let y = &w.y[p * s.n..(p + 1) * s.n];
        let history_s = &w.s[p * s.n..(p + 1) * s.n];
        let mut y_dot = 0.0;
        let mut s_dot = 0.0;
        for i in 0..s.n {
            if w.status[i] == 0 || w.status[i] == -1 {
                let direction = w.direction[i];
                y_dot += y[i] * direction;
                s_dot += history_s[i] * direction;
            }
        }
        w.wa[j] = y_dot;
        w.wa[s.col + j] = s.theta * s_dot;
    }
    w.cauchy.copy_from_slice(x);
    w.wa[2 * s.m..2 * s.m + c2].fill(0.0);
    if f1 == 0.0 {
        return true;
    }
    let mut f2 = -s.theta * f1;
    if s.col > 0 {
        let (left, right) = w.wa.split_at_mut(6 * s.m);
        let p = &left[..c2];
        let v = &mut right[..c2];
        if !bmv(s, &w.sy, &w.wt, p, v) {
            return false;
        }
        f2 -= dot(p, v);
    }
    let f2_org = f2;
    let mut dtm = -f1 / f2;
    let mut tsum = 0.0;
    let mut tj = 0.0;
    let mut nleft = nbreak;
    let mut iter = 0;
    while nleft > 0 {
        let tj0 = tj;
        let ibp;
        if iter == 0 {
            tj = bkmin;
            ibp = w.order[ibkmin];
        } else {
            if iter == 1 && ibkmin != nbreak - 1 {
                w.breakpoints[ibkmin] = w.breakpoints[nbreak - 1];
                w.order[ibkmin] = w.order[nbreak - 1];
            }
            let pair = heap_pop_min(
                &mut w.breakpoints[..nleft],
                &mut w.order[..nleft],
                iter == 1,
            );
            tj = pair.0;
            ibp = pair.1;
        }
        let dt = tj - tj0;
        if dtm < dt {
            break;
        }
        tsum += dt;
        nleft -= 1;
        iter += 1;
        let dibp = w.direction[ibp];
        w.direction[ibp] = 0.0;
        let zibp = if dibp > 0.0 {
            w.cauchy[ibp] = w.upper[ibp];
            w.status[ibp] = 2;
            w.upper[ibp] - x[ibp]
        } else {
            w.cauchy[ibp] = w.lower[ibp];
            w.status[ibp] = 1;
            w.lower[ibp] - x[ibp]
        };
        if nleft == 0 && nbreak == s.n {
            dtm = dt;
            break;
        }
        let dibp2 = dibp * dibp;
        f1 = f1 + dt * f2 + dibp2 - s.theta * dibp * zibp;
        f2 -= s.theta * dibp2;
        if s.col > 0 {
            let (a, rest) = w.wa.split_at_mut(2 * s.m);
            let p = &mut a[..c2];
            let (cbuf, rest) = rest.split_at_mut(2 * s.m);
            let c = &mut cbuf[..c2];
            let (wbpbuf, vbuf) = rest.split_at_mut(2 * s.m);
            let wbp = &mut wbpbuf[..c2];
            let v = &mut vbuf[..c2];
            for j in 0..c2 {
                c[j] += dt * p[j];
            }
            for j in 0..s.col {
                let q = slot(s, j);
                wbp[j] = vec_at(&w.y, s.n, q, ibp);
                wbp[s.col + j] = s.theta * vec_at(&w.s, s.n, q, ibp);
            }
            if !bmv(s, &w.sy, &w.wt, wbp, v) {
                return false;
            }
            let wmc = dot(c, v);
            let wmp = dot(p, v);
            let wmw = dot(wbp, v);
            for j in 0..c2 {
                p[j] -= dibp * wbp[j];
            }
            f1 += dibp * wmc;
            f2 += 2.0 * dibp * wmp - dibp2 * wmw;
        }
        f2 = f2.max(f64::EPSILON * f2_org);
        dtm = if nleft == 0 && bounded { 0.0 } else { -f1 / f2 };
    }
    dtm = dtm.max(0.0);
    tsum += dtm;
    for i in 0..s.n {
        w.cauchy[i] += tsum * w.direction[i];
        w.displacement[i] = w.cauchy[i] - x[i];
    }
    for j in 0..c2 {
        w.wa[2 * s.m + j] += dtm * w.wa[j];
    }
    true
}

fn freev(s: &mut Session, w: &mut Workspace) -> (usize, usize, bool) {
    let mut nenter = 0;
    let mut ileave = s.n;
    if s.constrained {
        for i in 0..s.nfree {
            let k = w.index[i];
            if w.status[k] > 0 {
                ileave -= 1;
                w.changes[ileave] = k;
            }
        }
        for i in s.nfree..s.n {
            let k = w.index[i];
            if w.status[k] <= 0 {
                w.changes[nenter] = k;
                nenter += 1;
            }
        }
    }
    let work = ileave < s.n || nenter > 0 || s.updated;
    s.nfree = 0;
    let mut active = s.n;
    for i in 0..s.n {
        if w.status[i] <= 0 {
            w.index[s.nfree] = i;
            s.nfree += 1;
        } else {
            active -= 1;
            w.index[active] = i;
        }
    }
    (nenter, ileave, work)
}

fn formk(s: &Session, nenter: usize, ileave: usize, w: &mut Workspace) -> Result<(), FormkError> {
    let m2 = 2 * s.m;
    let upcl = if s.updated {
        if s.updates > s.m {
            for j in 0..s.m - 1 {
                let count = s.m - j - 1;
                for k in 0..count {
                    let yy = at(&w.wn1, m2, j + 1 + k, j + 1);
                    let ss = at(&w.wn1, m2, s.m + j + 1 + k, s.m + j + 1);
                    set(&mut w.wn1, m2, j + k, j, yy);
                    set(&mut w.wn1, m2, s.m + j + k, s.m + j, ss);
                }
                for k in 0..s.m - 1 {
                    let cross = at(&w.wn1, m2, s.m + 1 + k, j + 1);
                    set(&mut w.wn1, m2, s.m + k, j, cross);
                }
            }
        }
        let newest = slot(s, s.col - 1);
        for j in 0..s.col {
            let pj = slot(s, j);
            let mut yy = 0.0;
            let mut ss = 0.0;
            let mut active_sy = 0.0;
            for q in 0..s.nfree {
                let k = w.index[q];
                yy += vec_at(&w.y, s.n, newest, k) * vec_at(&w.y, s.n, pj, k);
            }
            for q in s.nfree..s.n {
                let k = w.index[q];
                ss += vec_at(&w.s, s.n, newest, k) * vec_at(&w.s, s.n, pj, k);
                active_sy += vec_at(&w.s, s.n, newest, k) * vec_at(&w.y, s.n, pj, k);
            }
            set(&mut w.wn1, m2, s.col - 1, j, yy);
            set(&mut w.wn1, m2, s.m + s.col - 1, s.m + j, ss);
            set(&mut w.wn1, m2, s.m + s.col - 1, j, active_sy);
        }
        for j in 0..s.col {
            let pj = slot(s, j);
            let mut free_sy = 0.0;
            for q in 0..s.nfree {
                let k = w.index[q];
                free_sy += vec_at(&w.s, s.n, pj, k) * vec_at(&w.y, s.n, newest, k);
            }
            set(&mut w.wn1, m2, s.m + j, s.col - 1, free_sy);
        }
        s.col - 1
    } else {
        s.col
    };
    for i in 0..upcl {
        let pi = slot(s, i);
        for j in 0..=i {
            let pj = slot(s, j);
            let mut enter_yy = 0.0;
            let mut enter_ss = 0.0;
            let mut leave_yy = 0.0;
            let mut leave_ss = 0.0;
            for q in 0..nenter {
                let k = w.changes[q];
                enter_yy += vec_at(&w.y, s.n, pi, k) * vec_at(&w.y, s.n, pj, k);
                enter_ss += vec_at(&w.s, s.n, pi, k) * vec_at(&w.s, s.n, pj, k);
            }
            for q in ileave..s.n {
                let k = w.changes[q];
                leave_yy += vec_at(&w.y, s.n, pi, k) * vec_at(&w.y, s.n, pj, k);
                leave_ss += vec_at(&w.s, s.n, pi, k) * vec_at(&w.s, s.n, pj, k);
            }
            let yy = at(&w.wn1, m2, i, j) + enter_yy - leave_yy;
            let ss = at(&w.wn1, m2, s.m + i, s.m + j) - enter_ss + leave_ss;
            set(&mut w.wn1, m2, i, j, yy);
            set(&mut w.wn1, m2, s.m + i, s.m + j, ss);
        }
    }
    for i in 0..upcl {
        let pi = slot(s, i);
        for j in 0..upcl {
            let pj = slot(s, j);
            let mut enter = 0.0;
            let mut leave = 0.0;
            for q in 0..nenter {
                let k = w.changes[q];
                enter += vec_at(&w.s, s.n, pi, k) * vec_at(&w.y, s.n, pj, k);
            }
            for q in ileave..s.n {
                let k = w.changes[q];
                leave += vec_at(&w.s, s.n, pi, k) * vec_at(&w.y, s.n, pj, k);
            }
            let old = at(&w.wn1, m2, s.m + i, j);
            let value = if i <= j {
                old + enter - leave
            } else {
                old - enter + leave
            };
            set(&mut w.wn1, m2, s.m + i, j, value);
        }
    }
    w.wn.fill(0.0);
    for i in 0..s.col {
        for j in 0..=i {
            set(&mut w.wn, m2, j, i, at(&w.wn1, m2, i, j) / s.theta);
            set(
                &mut w.wn,
                m2,
                s.col + j,
                s.col + i,
                at(&w.wn1, m2, s.m + i, s.m + j) * s.theta,
            );
        }
        for j in 0..i {
            set(&mut w.wn, m2, j, s.col + i, -at(&w.wn1, m2, s.m + i, j));
        }
        for j in i..s.col {
            set(&mut w.wn, m2, j, s.col + i, at(&w.wn1, m2, s.m + i, j));
        }
        let diagonal = at(&w.wn, m2, i, i) + at(&w.sy, s.m, i, i);
        set(&mut w.wn, m2, i, i, diagonal);
    }
    if !cholesky_upper(&mut w.wn, m2, 0, s.col) {
        return Err(FormkError::FirstFactorization);
    }
    for j in s.col..2 * s.col {
        if !solve_matrix_column(&mut w.wn, m2, 0, s.col, j, true) {
            return Err(FormkError::TriangularSolve);
        }
    }
    for i in s.col..2 * s.col {
        for j in i..2 * s.col {
            let mut value = at(&w.wn, m2, i, j);
            for k in 0..s.col {
                value += at(&w.wn, m2, k, i) * at(&w.wn, m2, k, j);
            }
            set(&mut w.wn, m2, i, j, value);
        }
    }
    if !cholesky_upper(&mut w.wn, m2, s.col, s.col) {
        return Err(FormkError::SecondFactorization);
    }
    Ok(())
}

fn cmprlb(s: &Session, x: &[f64], w: &mut Workspace) -> bool {
    for i in 0..s.nfree {
        let k = w.index[i];
        w.reduced_gradient[i] = -s.theta * (w.cauchy[k] - x[k]) - w.evaluation_gradient[k];
    }
    let c2 = 2 * s.col;
    if s.col == 0 {
        return true;
    }
    let (vbuf, rest) = w.wa.split_at_mut(2 * s.m);
    let v = &mut vbuf[..c2];
    let c = &rest[..c2];
    if !bmv(s, &w.sy, &w.wt, c, v) {
        return false;
    }
    for j in 0..s.col {
        let p = slot(s, j);
        let a1 = v[j];
        let a2 = s.theta * v[s.col + j];
        for i in 0..s.nfree {
            let k = w.index[i];
            w.reduced_gradient[i] += vec_at(&w.y, s.n, p, k) * a1 + vec_at(&w.s, s.n, p, k) * a2;
        }
    }
    true
}

fn subsm(s: &Session, xx: &[f64], w: &mut Workspace) -> bool {
    let c2 = 2 * s.col;
    for j in 0..s.col {
        let p = slot(s, j);
        let mut a = 0.0;
        let mut b = 0.0;
        for i in 0..s.nfree {
            let k = w.index[i];
            a += vec_at(&w.y, s.n, p, k) * w.reduced_gradient[i];
            b += vec_at(&w.s, s.n, p, k) * w.reduced_gradient[i];
        }
        w.wa[j] = a;
        w.wa[s.col + j] = s.theta * b;
    }
    if !solve_upper(&w.wn, 2 * s.m, 0, c2, &mut w.wa[..c2], true) {
        return false;
    }
    for i in 0..s.col {
        w.wa[i] = -w.wa[i];
    }
    if !solve_upper(&w.wn, 2 * s.m, 0, c2, &mut w.wa[..c2], false) {
        return false;
    }
    for i in 0..s.nfree {
        let k = w.index[i];
        let mut d = w.reduced_gradient[i];
        for j in 0..s.col {
            let p = slot(s, j);
            d += vec_at(&w.y, s.n, p, k) * w.wa[j] / s.theta
                + vec_at(&w.s, s.n, p, k) * w.wa[s.col + j];
        }
        w.direction[i] = d / s.theta;
    }
    w.row.copy_from_slice(&w.cauchy);
    let mut projected = false;
    for i in 0..s.nfree {
        let k = w.index[i];
        let candidate = w.cauchy[k] + w.direction[i];
        w.cauchy[k] = candidate.clamp(w.lower[k], w.upper[k]);
        projected |= w.cauchy[k] == w.lower[k] || w.cauchy[k] == w.upper[k];
    }
    if !projected
        || w.cauchy
            .iter()
            .zip(xx)
            .zip(&w.evaluation_gradient)
            .map(|((&z, &x), &g)| (z - x) * g)
            .sum::<f64>()
            <= 0.0
    {
        return true;
    }
    w.cauchy.copy_from_slice(&w.row);
    let mut alpha = 1.0;
    let mut ibd = None;
    for i in 0..s.nfree {
        let k = w.index[i];
        let d = w.direction[i];
        let candidate = if d < 0.0 && w.lower[k].is_finite() {
            (w.lower[k] - w.cauchy[k]) / d
        } else if d > 0.0 && w.upper[k].is_finite() {
            (w.upper[k] - w.cauchy[k]) / d
        } else {
            1.0
        };
        if candidate < alpha {
            alpha = candidate.max(0.0);
            ibd = Some(i);
        }
    }
    if let Some(i) = ibd {
        let k = w.index[i];
        w.cauchy[k] = if w.direction[i] > 0.0 {
            w.upper[k]
        } else {
            w.lower[k]
        };
        w.direction[i] = 0.0;
    }
    for i in 0..s.nfree {
        let k = w.index[i];
        w.cauchy[k] += alpha * w.direction[i];
    }
    true
}

fn cholesky_upper(a: &mut [f64], ld: usize, offset: usize, n: usize) -> bool {
    for j in 0..n {
        let mut sumsq = 0.0;
        for k in 0..j {
            let mut value = at(a, ld, offset + k, offset + j);
            for i in 0..k {
                value -= at(a, ld, offset + i, offset + k) * at(a, ld, offset + i, offset + j);
            }
            value /= at(a, ld, offset + k, offset + k);
            set(a, ld, offset + k, offset + j, value);
            sumsq += value * value;
        }
        let diagonal = at(a, ld, offset + j, offset + j) - sumsq;
        if diagonal <= 0.0 {
            return false;
        }
        set(a, ld, offset + j, offset + j, diagonal.sqrt());
    }
    true
}

fn solve_upper(
    a: &[f64],
    ld: usize,
    offset: usize,
    n: usize,
    b: &mut [f64],
    transpose: bool,
) -> bool {
    if transpose {
        for j in 0..n {
            let diagonal = at(a, ld, offset + j, offset + j);
            if diagonal == 0.0 {
                return false;
            }
            let mut value = b[j];
            for i in 0..j {
                value -= at(a, ld, offset + i, offset + j) * b[i];
            }
            b[j] = value / diagonal;
        }
    } else {
        for j in (0..n).rev() {
            let diagonal = at(a, ld, offset + j, offset + j);
            if diagonal == 0.0 {
                return false;
            }
            b[j] /= diagonal;
            let value = b[j];
            for i in 0..j {
                b[i] -= value * at(a, ld, offset + i, offset + j);
            }
        }
    }
    true
}

fn solve_matrix_column(
    a: &mut [f64],
    ld: usize,
    offset: usize,
    n: usize,
    column: usize,
    transpose: bool,
) -> bool {
    for j in 0..n {
        let diagonal = at(a, ld, offset + j, offset + j);
        if diagonal == 0.0 {
            return false;
        }
        if transpose {
            let mut value = at(a, ld, offset + j, column);
            for i in 0..j {
                value -= at(a, ld, offset + i, offset + j) * at(a, ld, offset + i, column);
            }
            set(a, ld, offset + j, column, value / diagonal);
        }
    }
    true
}

fn heap_pop_min(values: &mut [f64], indices: &mut [usize], initialize: bool) -> (f64, usize) {
    let n = values.len();
    if initialize {
        for k in 1..n {
            let value = values[k];
            let index = indices[k];
            let mut i = k;
            while i > 0 {
                let parent = (i - 1) / 2;
                if value >= values[parent] {
                    break;
                }
                values[i] = values[parent];
                indices[i] = indices[parent];
                i = parent;
            }
            values[i] = value;
            indices[i] = index;
        }
    }
    if n == 1 {
        return (values[0], indices[0]);
    }
    let result = (values[0], indices[0]);
    let value = values[n - 1];
    let index = indices[n - 1];
    let mut i = 0;
    loop {
        let mut child = 2 * i + 1;
        if child >= n - 1 {
            break;
        }
        if child + 1 < n - 1 && values[child + 1] < values[child] {
            child += 1;
        }
        if values[child] >= value {
            break;
        }
        values[i] = values[child];
        indices[i] = indices[child];
        i = child;
    }
    values[i] = value;
    indices[i] = index;
    values[n - 1] = result.0;
    indices[n - 1] = result.1;
    result
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(&a, &b)| a * b).sum()
}

#[cfg(feature = "benchmark-instrumentation")]
fn elapsed_nanoseconds(start: Instant) -> u64 {
    start.elapsed().as_nanos().min(u64::MAX as u128) as u64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn compact_cholesky_and_triangular_solves() {
        let mut a = vec![4.0, 0.0, 2.0, 3.0];
        assert!(cholesky_upper(&mut a, 2, 0, 2));
        let mut b = [1.0, 2.0];
        assert!(solve_upper(&a, 2, 0, 2, &mut b, true));
        assert!(solve_upper(&a, 2, 0, 2, &mut b, false));
        assert!((4.0 * b[0] + 2.0 * b[1] - 1.0).abs() < 1e-12);
        assert!((2.0 * b[0] + 3.0 * b[1] - 2.0).abs() < 1e-12);
    }

    #[test]
    fn curvature_skip_does_not_pollute_history() {
        let mut w = Workspace::new();
        w.prepare(&[-10.0], &[10.0], 2).unwrap();
        w.old_gradient[0] = -1.0;
        w.evaluation_gradient[0] = -1.0;
        let mut session = Session::new(1, 2, Kernel::Scalar);
        assert!(!session.update(&[1.0], -1.0, &mut w));
        assert!(!session.has_history());
    }

    #[test]
    fn curvature_skip_threshold_matches_fortran_inequality() {
        fn attempt(curvature: f64) -> Session {
            let mut w = Workspace::new();
            w.prepare(&[-10.0], &[10.0], 2).unwrap();
            w.old_x[0] = 0.0;
            w.old_gradient[0] = -1.0;
            w.evaluation_gradient[0] = -1.0 + curvature;
            let mut session = Session::new(1, 2, Kernel::Scalar);
            session.update(&[1.0], -1.0, &mut w);
            session
        }

        // mainlb skips when dr <= epsmch*(-gdold); equality is intentional.
        let skipped = attempt(f64::EPSILON);
        assert_eq!(skipped.skipped_updates(), 1);
        assert_eq!(skipped.accepted_updates(), 0);
        assert!(!skipped.has_history());

        let accepted = attempt(2.0 * f64::EPSILON);
        assert_eq!(accepted.skipped_updates(), 0);
        assert_eq!(accepted.accepted_updates(), 1);
        assert!(accepted.has_history());
    }

    #[test]
    fn tied_heap_breakpoints_follow_hpsolb_strict_order() {
        // hpsolb is not observable through NEW_X. This checks its defining
        // invariant directly: strict comparisons preserve insertion order for
        // equal minima while each pop leaves a valid minimum heap.
        let mut values = [1.0, 1.0, 1.0, 2.0];
        let mut indices = [10, 11, 12, 13];
        let first = heap_pop_min(&mut values, &mut indices, true);
        let second = heap_pop_min(&mut values[..3], &mut indices[..3], false);
        let third = heap_pop_min(&mut values[..2], &mut indices[..2], false);
        assert_eq!([first, second, third], [(1.0, 10), (1.0, 11), (1.0, 12)]);
    }

    #[test]
    fn subsm_positive_projection_uses_v3_backtracking_correction() {
        // The v3 correction decision is private to subsm and is not exported
        // by setulb. An artificial reduced problem isolates the invariant:
        // after a positive projected directional derivative, move only to the
        // first encountered bound instead of retaining the full projection.
        let mut w = Workspace::new();
        w.prepare(
            &[-1.0, f64::NEG_INFINITY],
            &[f64::INFINITY, f64::INFINITY],
            1,
        )
        .unwrap();
        let mut session = Session::new(2, 1, Kernel::Scalar);
        session.nfree = 2;
        w.index.copy_from_slice(&[0, 1]);
        w.cauchy.copy_from_slice(&[0.0, 0.0]);
        w.reduced_gradient.copy_from_slice(&[-2.0, 0.2]);
        w.evaluation_gradient.copy_from_slice(&[-2.0, 0.0]);

        assert!(subsm(&session, &[0.0, 0.0], &mut w));
        assert_eq!(w.cauchy, [-1.0, 0.1]);
    }

    #[test]
    fn compact_factorization_failure_resets_and_retries_without_history() {
        // setulb exposes only the eventual restarted iterate, not which
        // compact factorization failed. Corrupting the cached free-set block
        // deterministically exercises direction's one-reset recovery contract.
        let mut w = Workspace::new();
        w.prepare(&[0.0, -10.0], &[10.0, 10.0], 1).unwrap();
        w.evaluation_gradient.copy_from_slice(&[1.0, 1.0]);
        w.index.copy_from_slice(&[0, 1]);
        w.sy[0] = 1.0;
        w.wt[0] = 1.0;
        w.wn1[0] = -10.0;

        let mut session = Session::new(2, 1, Kernel::Scalar);
        session.col = 1;
        session.updates = 1;
        session.nfree = 2;

        assert_eq!(session.direction(&[0.0, 0.0], &mut w), Ok(()));
        assert_eq!(session.history_resets(), 1);
        assert_eq!(session.skipped_updates(), 1);
        assert!(!session.has_history());
        assert!(w.direction[1] < 0.0);
    }
}

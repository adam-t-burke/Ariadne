//! q-parameterization maps between optimizer variables and physical q.
//!
//! DirectSoftBounds mode uses identity (q = z).
//! ImplicitBounded mode maps z -> q with hard bounds and chain-rule derivative.

use crate::types::QParameterizationMode;
use std::f64::consts::FRAC_2_PI;

const TWO_SIDED_BETA: f64 = 1.0;
const ONE_SIDED_ALPHA: f64 = 1.0;
const INVERSE_CLAMP_EPS: f64 = 1e-9;

#[inline]
fn sigmoid(x: f64) -> f64 {
    if x >= 0.0 {
        let e = (-x).exp();
        1.0 / (1.0 + e)
    } else {
        let e = x.exp();
        e / (1.0 + e)
    }
}

#[inline]
fn softplus(x: f64) -> f64 {
    if x > 0.0 {
        x + (-x).exp().ln_1p()
    } else {
        x.exp().ln_1p()
    }
}

#[inline]
fn inverse_softplus(y: f64) -> f64 {
    // Inverse of softplus on y > 0: x = ln(exp(y) - 1)
    if y <= INVERSE_CLAMP_EPS {
        return (INVERSE_CLAMP_EPS.exp() - 1.0).ln();
    }
    if y > 50.0 {
        y
    } else {
        y.exp_m1().ln()
    }
}

#[inline]
fn map_bounded(z: f64, lb: f64, ub: f64) -> f64 {
    let lb_f = lb.is_finite();
    let ub_f = ub.is_finite();
    match (lb_f, ub_f) {
        (false, false) => z,
        (true, true) => {
            let center = 0.5 * (lb + ub);
            let half_span = 0.5 * (ub - lb);
            let s = FRAC_2_PI * (TWO_SIDED_BETA * z).atan();
            center + half_span * s
        }
        (true, false) => lb + softplus(ONE_SIDED_ALPHA * z) / ONE_SIDED_ALPHA,
        (false, true) => ub - softplus(-ONE_SIDED_ALPHA * z) / ONE_SIDED_ALPHA,
    }
}

#[inline]
fn dq_dz_bounded(z: f64, lb: f64, ub: f64) -> f64 {
    let lb_f = lb.is_finite();
    let ub_f = ub.is_finite();
    match (lb_f, ub_f) {
        (false, false) => 1.0,
        (true, true) => {
            let half_span = 0.5 * (ub - lb);
            let bz = TWO_SIDED_BETA * z;
            half_span * FRAC_2_PI * TWO_SIDED_BETA / (1.0 + bz * bz)
        }
        (true, false) => sigmoid(ONE_SIDED_ALPHA * z),
        (false, true) => sigmoid(-ONE_SIDED_ALPHA * z),
    }
}

#[inline]
fn inverse_bounded(q: f64, lb: f64, ub: f64) -> f64 {
    let lb_f = lb.is_finite();
    let ub_f = ub.is_finite();
    match (lb_f, ub_f) {
        (false, false) => q,
        (true, true) => {
            let center = 0.5 * (lb + ub);
            let half_span = 0.5 * (ub - lb);
            let mut s = (q - center) / half_span.max(f64::EPSILON);
            s = s.clamp(-1.0 + INVERSE_CLAMP_EPS, 1.0 - INVERSE_CLAMP_EPS);
            (s / FRAC_2_PI).tan() / TWO_SIDED_BETA
        }
        (true, false) => {
            let y = ONE_SIDED_ALPHA * (q - lb).max(INVERSE_CLAMP_EPS);
            inverse_softplus(y) / ONE_SIDED_ALPHA
        }
        (false, true) => {
            let y = ONE_SIDED_ALPHA * (ub - q).max(INVERSE_CLAMP_EPS);
            -inverse_softplus(y) / ONE_SIDED_ALPHA
        }
    }
}

#[inline]
pub fn map_single(mode: QParameterizationMode, z: f64, lb: f64, ub: f64) -> f64 {
    match mode {
        QParameterizationMode::DirectSoftBounds => z,
        QParameterizationMode::ImplicitBounded => map_bounded(z, lb, ub),
    }
}

#[inline]
pub fn dq_dz_single(mode: QParameterizationMode, z: f64, lb: f64, ub: f64) -> f64 {
    match mode {
        QParameterizationMode::DirectSoftBounds => 1.0,
        QParameterizationMode::ImplicitBounded => dq_dz_bounded(z, lb, ub),
    }
}

#[inline]
pub fn inverse_single(mode: QParameterizationMode, q: f64, lb: f64, ub: f64) -> f64 {
    match mode {
        QParameterizationMode::DirectSoftBounds => q,
        QParameterizationMode::ImplicitBounded => inverse_bounded(q, lb, ub),
    }
}

pub fn map_q_slice(
    mode: QParameterizationMode,
    z_q: &[f64],
    lb_q: &[f64],
    ub_q: &[f64],
    out_q: &mut [f64],
) {
    for i in 0..z_q.len() {
        out_q[i] = map_single(mode, z_q[i], lb_q[i], ub_q[i]);
    }
}


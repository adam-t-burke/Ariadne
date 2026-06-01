//! q-parameterization maps between optimizer variables and physical q.
//!
//! DirectSoftBounds mode uses identity (q = z).
//! ImplicitBounded mode maps z -> q with hard bounds and chain-rule derivative.

use crate::types::QParameterizationMode;

const TWO_SIDED_GAMMA: f64 = 1.0;
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
fn logit(p: f64) -> f64 {
    (p / (1.0 - p)).ln()
}

#[inline]
fn map_bounded(z: f64, lb: f64, ub: f64) -> f64 {
    let lb_f = lb.is_finite();
    let ub_f = ub.is_finite();
    match (lb_f, ub_f) {
        (false, false) => z,
        (true, true) => {
            let span = ub - lb;
            lb + span * sigmoid(TWO_SIDED_GAMMA * z)
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
            let span = ub - lb;
            let s = sigmoid(TWO_SIDED_GAMMA * z);
            span * TWO_SIDED_GAMMA * s * (1.0 - s)
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
            let span = ub - lb;
            if span.abs() <= f64::EPSILON {
                return 0.0;
            }
            let p = ((q - lb) / span).clamp(INVERSE_CLAMP_EPS, 1.0 - INVERSE_CLAMP_EPS);
            logit(p) / TWO_SIDED_GAMMA
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
        QParameterizationMode::DirectSoftBounds | QParameterizationMode::DirectBoxBounds => z,
        QParameterizationMode::ImplicitBounded => map_bounded(z, lb, ub),
    }
}

#[inline]
pub fn dq_dz_single(mode: QParameterizationMode, z: f64, lb: f64, ub: f64) -> f64 {
    match mode {
        QParameterizationMode::DirectSoftBounds | QParameterizationMode::DirectBoxBounds => 1.0,
        QParameterizationMode::ImplicitBounded => dq_dz_bounded(z, lb, ub),
    }
}

#[inline]
pub fn inverse_single(mode: QParameterizationMode, q: f64, lb: f64, ub: f64) -> f64 {
    match mode {
        QParameterizationMode::DirectSoftBounds | QParameterizationMode::DirectBoxBounds => q,
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

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1e-10;

    #[test]
    fn direct_soft_bounds_remains_identity() {
        let z = 42.0;
        let lb = -2.0;
        let ub = 3.0;

        assert_eq!(
            map_single(QParameterizationMode::DirectSoftBounds, z, lb, ub),
            z
        );
        assert_eq!(
            map_single(QParameterizationMode::DirectBoxBounds, z, lb, ub),
            z
        );
        assert_eq!(
            dq_dz_single(QParameterizationMode::DirectSoftBounds, z, lb, ub),
            1.0
        );
        assert_eq!(
            dq_dz_single(QParameterizationMode::DirectBoxBounds, z, lb, ub),
            1.0
        );
        assert_eq!(
            inverse_single(QParameterizationMode::DirectSoftBounds, z, lb, ub),
            z
        );
        assert_eq!(
            inverse_single(QParameterizationMode::DirectBoxBounds, z, lb, ub),
            z
        );
    }

    #[test]
    fn implicit_two_sided_logit_maps_center_and_stays_inside_box() {
        let lb = -2.0;
        let ub = 8.0;

        let center = map_single(QParameterizationMode::ImplicitBounded, 0.0, lb, ub);
        assert!((center - 3.0).abs() < TOL);

        for z in [-20.0, -2.0, 0.0, 2.0, 20.0] {
            let q = map_single(QParameterizationMode::ImplicitBounded, z, lb, ub);
            assert!(q > lb, "q={q} should be above lb={lb}");
            assert!(q < ub, "q={q} should be below ub={ub}");
        }
    }

    #[test]
    fn implicit_two_sided_round_trips_interior_values() {
        let lb = 0.1;
        let ub = 100.0;

        for q in [0.2, 1.0, 25.0, 50.05, 75.0, 99.0] {
            let z = inverse_single(QParameterizationMode::ImplicitBounded, q, lb, ub);
            let mapped = map_single(QParameterizationMode::ImplicitBounded, z, lb, ub);
            assert!(
                (mapped - q).abs() < 1e-9,
                "round trip failed: q={q}, z={z}, mapped={mapped}"
            );
        }
    }

    #[test]
    fn implicit_two_sided_derivative_matches_finite_difference() {
        let lb = -5.0;
        let ub = 7.0;
        let h = 1e-6;

        for z in [-4.0, -1.0, 0.0, 1.0, 4.0] {
            let analytic = dq_dz_single(QParameterizationMode::ImplicitBounded, z, lb, ub);
            let plus = map_single(QParameterizationMode::ImplicitBounded, z + h, lb, ub);
            let minus = map_single(QParameterizationMode::ImplicitBounded, z - h, lb, ub);
            let numeric = (plus - minus) / (2.0 * h);
            assert!(
                (analytic - numeric).abs() < 1e-8,
                "z={z}, analytic={analytic}, numeric={numeric}"
            );
        }
    }
}

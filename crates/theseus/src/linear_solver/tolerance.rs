//! Evaluation of [`TolerancePolicy`] and the per-evaluation tolerance
//! schedule of the iterative linear solvers (§3 "Tolerance policy" of the
//! program plan).
//!
//! * [`TolerancePolicy::Fixed`] uses the same relative residual for every
//!   solve.
//! * [`TolerancePolicy::Adaptive`] ties the accuracy of the linear solves to
//!   the progress of the outer optimisation:
//!   `tol_k = clamp(factor · ‖g⁺_k‖ / ‖g⁺_0‖, floor, ceiling)`, where `g⁺`
//!   is the projected gradient (box-bounded mode) or the gradient (soft
//!   bounds) of the **previous accepted iterate**. The first evaluation, and
//!   every evaluation before the reference norm `‖g⁺_0‖` is known, uses
//!   `ceiling`.
//!
//! [`ToleranceSchedule`] is the stateful form driven by `optimizer.rs`: the
//! optimizer feeds it the gradient norm of each accepted iterate and reads
//! the tolerance to store in `FdmCache::solve_tolerance` before the next
//! evaluations. `Direct` ignores the tolerance entirely.

use super::TolerancePolicy;
use crate::types::TheseusError;

impl TolerancePolicy {
    /// Check that the policy can produce a usable tolerance: `Fixed(tol)`
    /// needs `0 < tol < ∞`; `Adaptive` needs `0 < floor ≤ ceiling < ∞` and
    /// `0 < factor < ∞`.
    pub fn validate(&self) -> Result<(), TheseusError> {
        let positive_finite = |v: f64| v.is_finite() && v > 0.0;
        match *self {
            Self::Fixed(tol) if positive_finite(tol) => Ok(()),
            Self::Fixed(tol) => Err(TheseusError::Shape(format!(
                "iterative solver tolerance must be a positive finite number, got {tol}"
            ))),
            Self::Adaptive {
                floor,
                ceiling,
                factor,
            } => {
                if !(positive_finite(floor) && positive_finite(ceiling) && floor <= ceiling) {
                    return Err(TheseusError::Shape(format!(
                        "adaptive iterative solver tolerance needs 0 < floor <= ceiling (finite), \
                         got floor = {floor}, ceiling = {ceiling}"
                    )));
                }
                if !positive_finite(factor) {
                    return Err(TheseusError::Shape(format!(
                        "adaptive iterative solver tolerance factor must be a positive finite \
                         number, got {factor}"
                    )));
                }
                Ok(())
            }
        }
    }

    /// The tolerance for an evaluation whose previous accepted iterate had
    /// (projected) gradient norm `gradient_norm`, relative to the norm
    /// `reference_norm` at the initial point.
    ///
    /// `Fixed` returns its tolerance. `Adaptive` returns
    /// `clamp(factor · gradient_norm / reference_norm, floor, ceiling)`; when
    /// either norm is not finite, or the reference is zero (no information,
    /// or the optimisation is already stationary), it returns `ceiling`.
    pub fn tolerance_for(&self, gradient_norm: f64, reference_norm: f64) -> f64 {
        match *self {
            Self::Fixed(tol) => tol,
            Self::Adaptive {
                floor,
                ceiling,
                factor,
            } => {
                if !(reference_norm.is_finite() && reference_norm > 0.0 && gradient_norm.is_finite())
                {
                    return ceiling;
                }
                let tol = factor * (gradient_norm / reference_norm);
                if !tol.is_finite() {
                    return ceiling;
                }
                // `max` then `min` rather than `f64::clamp`, so an inverted
                // pair (rejected by `validate`, but the FFI may bypass it on
                // `Direct`) resolves to the ceiling instead of panicking.
                tol.max(floor).min(ceiling)
            }
        }
    }
}

/// The per-evaluation tolerance schedule of one optimisation run.
///
/// Starts at [`TolerancePolicy::initial`]. The first call to
/// [`ToleranceSchedule::observe`] fixes the reference norm `‖g⁺_0‖`; every
/// call recomputes the current tolerance from the norm it is given.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ToleranceSchedule {
    policy: TolerancePolicy,
    reference: Option<f64>,
    current: f64,
}

impl ToleranceSchedule {
    /// A schedule at its initial tolerance, with no reference norm yet.
    pub fn new(policy: TolerancePolicy) -> Self {
        Self {
            policy,
            reference: None,
            current: policy.initial(),
        }
    }

    /// The policy being evaluated.
    pub fn policy(&self) -> TolerancePolicy {
        self.policy
    }

    /// Tolerance for the next evaluations.
    pub fn current(&self) -> f64 {
        self.current
    }

    /// `‖g⁺_0‖` once the first accepted iterate has been observed.
    pub fn reference(&self) -> Option<f64> {
        self.reference
    }

    /// Record the (projected) gradient norm of the latest accepted iterate
    /// and return the tolerance for the evaluations that follow it. The
    /// first observation becomes the reference `‖g⁺_0‖` (a non-finite or
    /// non-positive first norm is ignored so a later one can take its place).
    pub fn observe(&mut self, gradient_norm: f64) -> f64 {
        if self.reference.is_none() && gradient_norm.is_finite() && gradient_norm > 0.0 {
            self.reference = Some(gradient_norm);
        }
        self.current = match self.reference {
            Some(reference) => self.policy.tolerance_for(gradient_norm, reference),
            None => self.policy.initial(),
        };
        self.current
    }

    /// Forget the reference norm and return to the initial tolerance.
    pub fn reset(&mut self) {
        self.reference = None;
        self.current = self.policy.initial();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const ADAPTIVE: TolerancePolicy = TolerancePolicy::Adaptive {
        floor: 1e-10,
        ceiling: 1e-6,
        factor: 1e-2,
    };

    #[test]
    fn fixed_policy_is_constant() {
        let p = TolerancePolicy::Fixed(1e-12);
        assert_eq!(p.initial(), 1e-12);
        assert_eq!(p.tolerance_for(1.0, 1.0), 1e-12);
        assert_eq!(p.tolerance_for(0.0, 1.0), 1e-12);
        assert_eq!(p.tolerance_for(f64::NAN, 0.0), 1e-12);
        let mut s = ToleranceSchedule::new(p);
        assert_eq!(s.observe(3.0), 1e-12);
        assert_eq!(s.observe(0.0), 1e-12);
        assert_eq!(s.current(), 1e-12);
    }

    #[test]
    fn adaptive_policy_clamps_the_gradient_ratio() {
        // ratio 1 → factor, capped at the ceiling.
        assert_eq!(ADAPTIVE.tolerance_for(1.0, 1.0), 1e-6);
        // ratio 1e-3 → 1e-5 → ceiling; ratio 1e-5 → 1e-7 (inside); 1e-9 → floor.
        assert_eq!(ADAPTIVE.tolerance_for(1e-3, 1.0), 1e-6);
        assert!((ADAPTIVE.tolerance_for(1e-5, 1.0) - 1e-7).abs() < 1e-22);
        assert_eq!(ADAPTIVE.tolerance_for(1e-9, 1.0), 1e-10);
        assert_eq!(ADAPTIVE.tolerance_for(0.0, 1.0), 1e-10);
        // Reference scales the ratio.
        assert!((ADAPTIVE.tolerance_for(2e-5, 2.0) - 1e-7).abs() < 1e-22);
        // No information → ceiling.
        assert_eq!(ADAPTIVE.tolerance_for(1.0, 0.0), 1e-6);
        assert_eq!(ADAPTIVE.tolerance_for(f64::NAN, 1.0), 1e-6);
        assert_eq!(ADAPTIVE.tolerance_for(1.0, f64::INFINITY), 1e-6);
        assert_eq!(ADAPTIVE.tolerance_for(f64::INFINITY, 1.0), 1e-6);
    }

    #[test]
    fn schedule_uses_the_first_observation_as_reference() {
        let mut s = ToleranceSchedule::new(ADAPTIVE);
        assert_eq!(s.current(), 1e-6, "first evaluation uses the ceiling");
        assert_eq!(s.reference(), None);
        assert_eq!(s.observe(10.0), 1e-6, "ratio 1 → factor 1e-2 → ceiling");
        assert_eq!(s.reference(), Some(10.0));
        assert!((s.observe(1e-4) - 1e-7).abs() < 1e-22);
        assert!((s.current() - 1e-7).abs() < 1e-22);
        assert_eq!(s.observe(1e-12), 1e-10);
        // The gradient may grow again (rejected region): tolerance loosens
        // but never beyond the ceiling.
        assert_eq!(s.observe(1e3), 1e-6);
        s.reset();
        assert_eq!(s.reference(), None);
        assert_eq!(s.current(), 1e-6);
    }

    #[test]
    fn schedule_skips_unusable_first_observations() {
        let mut s = ToleranceSchedule::new(ADAPTIVE);
        assert_eq!(s.observe(0.0), 1e-6);
        assert_eq!(s.reference(), None);
        assert_eq!(s.observe(f64::NAN), 1e-6);
        assert_eq!(s.reference(), None);
        s.observe(4.0);
        assert_eq!(s.reference(), Some(4.0));
    }

    #[test]
    fn validation_rejects_degenerate_policies() {
        assert!(TolerancePolicy::Fixed(1e-8).validate().is_ok());
        assert!(TolerancePolicy::default().validate().is_ok());
        for bad in [0.0, -1e-8, f64::NAN, f64::INFINITY] {
            assert!(TolerancePolicy::Fixed(bad).validate().is_err(), "{bad}");
        }
        let adaptive = |floor, ceiling, factor| TolerancePolicy::Adaptive {
            floor,
            ceiling,
            factor,
        };
        assert!(adaptive(1e-6, 1e-10, 1e-2).validate().is_err(), "inverted");
        assert!(adaptive(0.0, 1e-6, 1e-2).validate().is_err(), "zero floor");
        assert!(adaptive(1e-10, f64::INFINITY, 1e-2).validate().is_err());
        assert!(adaptive(1e-10, 1e-6, 0.0).validate().is_err(), "zero factor");
        assert!(adaptive(1e-10, 1e-6, f64::NAN).validate().is_err());
        assert!(adaptive(1e-8, 1e-8, 1.0).validate().is_ok(), "floor == ceiling");
        // An inverted pair still evaluates without panicking (ceiling wins).
        assert_eq!(adaptive(1e-6, 1e-10, 1e-2).tolerance_for(1.0, 1.0), 1e-10);
    }
}

//! Borrowed box constraints.

use crate::BoundsError;

/// Lower and upper bounds for each variable.
///
/// Infinite endpoints represent one-sided or unbounded variables. Bounds are
/// borrowed, so constructing this view does not allocate.
#[derive(Debug, Clone, Copy)]
pub struct Bounds<'a> {
    lower: &'a [f64],
    upper: &'a [f64],
}

impl<'a> Bounds<'a> {
    /// Creates bounds and checks their lengths and ordering.
    pub fn new(lower: &'a [f64], upper: &'a [f64], variables: usize) -> Result<Self, BoundsError> {
        let bounds = Self { lower, upper };
        bounds.validate(variables)?;
        Ok(bounds)
    }

    /// Returns the lower endpoints.
    pub fn lower(self) -> &'a [f64] {
        self.lower
    }

    /// Returns the upper endpoints.
    pub fn upper(self) -> &'a [f64] {
        self.upper
    }

    /// Returns the number of bounded variables represented by this view.
    pub fn len(self) -> usize {
        self.lower.len()
    }

    /// Returns whether this view contains no variables.
    pub fn is_empty(self) -> bool {
        self.lower.is_empty()
    }

    pub(crate) fn validate(self, variables: usize) -> Result<(), BoundsError> {
        if variables == 0 {
            return Err(BoundsError::Empty);
        }
        if self.lower.len() != variables || self.upper.len() != variables {
            return Err(BoundsError::LengthMismatch {
                variables,
                lower: self.lower.len(),
                upper: self.upper.len(),
            });
        }
        for (index, (&lower, &upper)) in self.lower.iter().zip(self.upper).enumerate() {
            if lower.is_nan() || upper.is_nan() {
                return Err(BoundsError::NotANumber { index });
            }
            if lower == upper && lower.is_infinite() {
                return Err(BoundsError::InfiniteFixed {
                    index,
                    value: lower,
                });
            }
            if lower == f64::INFINITY {
                return Err(BoundsError::PositiveInfiniteLower { index });
            }
            if upper == f64::NEG_INFINITY {
                return Err(BoundsError::NegativeInfiniteUpper { index });
            }
            if lower > upper {
                return Err(BoundsError::Inverted {
                    index,
                    lower,
                    upper,
                });
            }
        }
        Ok(())
    }
}

/// Computes the infinity norm of the gradient projected onto feasible
/// directions.
pub fn projected_gradient_norm(
    x: &[f64],
    gradient: &[f64],
    bounds: Bounds<'_>,
) -> Result<f64, BoundsError> {
    bounds.validate(x.len())?;
    if gradient.len() != x.len() {
        return Err(BoundsError::LengthMismatch {
            variables: x.len(),
            lower: gradient.len(),
            upper: bounds.len(),
        });
    }
    if let Some(index) = x.iter().position(|value| !value.is_finite()) {
        return Err(BoundsError::NonFiniteVariable { index });
    }
    if let Some(index) = gradient.iter().position(|value| !value.is_finite()) {
        return Err(BoundsError::NonFiniteGradient { index });
    }

    Ok(projected_gradient_and_active_count_unchecked(
        x,
        gradient,
        bounds.lower,
        bounds.upper,
    )
    .0)
}

pub(crate) fn projected_gradient_and_active_count_unchecked(
    x: &[f64],
    gradient: &[f64],
    lower: &[f64],
    upper: &[f64],
) -> (f64, usize) {
    x.iter()
        .zip(gradient)
        .zip(lower.iter().zip(upper))
        .fold(
            (0.0_f64, 0usize),
            |(norm, active), ((&x, &gradient), (&lower, &upper))| {
                let component = if gradient < 0.0 {
                    gradient.max(x - upper)
                } else {
                    gradient.min(x - lower)
                };
                let is_active = lower == upper
                    || (x <= lower && gradient >= 0.0)
                    || (x >= upper && gradient <= 0.0);
                (norm.max(component.abs()), active + usize::from(is_active))
            },
        )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn accepts_mixed_infinite_bounds() {
        let bounds = Bounds::new(&[f64::NEG_INFINITY, 0.0], &[2.0, f64::INFINITY], 2);
        assert!(bounds.is_ok());
    }

    #[test]
    fn rejects_inverted_bounds_with_index() {
        assert!(matches!(
            Bounds::new(&[2.0], &[1.0], 1),
            Err(BoundsError::Inverted { index: 0, .. })
        ));
    }

    #[test]
    fn projected_gradient_ignores_blocked_components() {
        let bounds = Bounds::new(&[0.0, -2.0], &[3.0, 2.0], 2).unwrap();
        assert_eq!(
            projected_gradient_norm(&[0.0, 1.0], &[4.0, -3.0], bounds).unwrap(),
            1.0
        );
    }

    #[test]
    fn rejects_invalid_infinite_endpoints() {
        assert!(matches!(
            Bounds::new(&[f64::INFINITY], &[f64::INFINITY], 1),
            Err(BoundsError::InfiniteFixed { .. })
        ));
        assert!(matches!(
            Bounds::new(&[f64::NEG_INFINITY], &[f64::NEG_INFINITY], 1),
            Err(BoundsError::InfiniteFixed { .. })
        ));
        assert!(matches!(
            Bounds::new(&[f64::INFINITY], &[1.0], 1),
            Err(BoundsError::PositiveInfiniteLower { index: 0 })
        ));
        assert!(matches!(
            Bounds::new(&[-1.0], &[f64::NEG_INFINITY], 1),
            Err(BoundsError::NegativeInfiniteUpper { index: 0 })
        ));
    }

    #[test]
    fn projected_gradient_rejects_nonfinite_inputs() {
        let bounds = Bounds::new(&[f64::NEG_INFINITY], &[f64::INFINITY], 1).unwrap();
        assert!(matches!(
            projected_gradient_norm(&[f64::INFINITY], &[1.0], bounds),
            Err(BoundsError::NonFiniteVariable { index: 0 })
        ));
        assert!(matches!(
            projected_gradient_norm(&[0.0], &[f64::NAN], bounds),
            Err(BoundsError::NonFiniteGradient { index: 0 })
        ));
    }
}

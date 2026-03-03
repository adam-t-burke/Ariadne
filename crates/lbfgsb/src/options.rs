//! Checked solver configuration.

use crate::OptionsError;

/// Dense-history arithmetic backend.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum Backend {
    /// Deterministic scalar arithmetic in a fixed accumulation order.
    Deterministic,
    /// Automatically choose a backend.
    ///
    /// Auto currently resolves to [`Backend::Deterministic`]. Future automatic
    /// selection changes will be based on representative cross-objective data.
    Auto,
    /// Faer matrix-vector products for dense history updates.
    ///
    /// Selecting this backend requires the `faer-backend` feature.
    Faer,
}

/// Validated L-BFGS-B configuration.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Options {
    history_size: usize,
    max_iterations: usize,
    max_evaluations: usize,
    max_line_search_evaluations: usize,
    relative_function_tolerance: f64,
    projected_gradient_tolerance: f64,
    backend: Backend,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            history_size: 10,
            max_iterations: 1_000,
            max_evaluations: usize::MAX,
            max_line_search_evaluations: 20,
            relative_function_tolerance: 1.0e7 * f64::EPSILON,
            projected_gradient_tolerance: 1.0e-5,
            backend: Backend::Auto,
        }
    }
}

impl Options {
    /// Largest supported correction-history size.
    ///
    /// L-BFGS methods conventionally use small histories; this generous bound
    /// prevents accidental quadratic workspace requests while retaining ample
    /// room for specialized workloads.
    pub const MAX_HISTORY_SIZE: usize = 1_024;

    /// Creates the default, valid configuration.
    pub fn new() -> Self {
        Self::default()
    }

    /// Sets the number of correction pairs retained by L-BFGS-B.
    pub fn with_history_size(mut self, value: usize) -> Result<Self, OptionsError> {
        if value == 0 {
            return Err(OptionsError::ZeroHistory);
        }
        if value > Self::MAX_HISTORY_SIZE {
            return Err(OptionsError::HistoryTooLarge {
                requested: value,
                maximum: Self::MAX_HISTORY_SIZE,
            });
        }
        self.history_size = value;
        Ok(self)
    }

    /// Selects the dense-history arithmetic backend.
    pub fn with_backend(mut self, value: Backend) -> Result<Self, OptionsError> {
        if matches!(value, Backend::Faer) && !cfg!(feature = "faer-backend") {
            return Err(OptionsError::UnsupportedBackend);
        }
        self.backend = value;
        Ok(self)
    }

    /// Sets the exact maximum number of accepted major iterations.
    pub fn with_max_iterations(mut self, value: usize) -> Result<Self, OptionsError> {
        if value == 0 {
            return Err(OptionsError::ZeroIterations);
        }
        self.max_iterations = value;
        Ok(self)
    }

    /// Sets the maximum number of objective/gradient evaluations.
    pub fn with_max_evaluations(mut self, value: usize) -> Result<Self, OptionsError> {
        if value == 0 {
            return Err(OptionsError::ZeroEvaluations);
        }
        self.max_evaluations = value;
        Ok(self)
    }

    /// Sets the maximum number of objective evaluations in one line search.
    pub fn with_max_line_search_evaluations(mut self, value: usize) -> Result<Self, OptionsError> {
        if value == 0 {
            return Err(OptionsError::ZeroLineSearchEvaluations);
        }
        self.max_line_search_evaluations = value;
        Ok(self)
    }

    /// Sets the relative accepted-iterate function-reduction tolerance.
    ///
    /// The test is `(f[k] - f[k+1]) / max(|f[k]|, |f[k+1]|, 1) <= value`.
    /// Set this to zero to disable convergence through this test.
    pub fn with_relative_function_tolerance(mut self, value: f64) -> Result<Self, OptionsError> {
        if !value.is_finite() || value < 0.0 {
            return Err(OptionsError::InvalidRelativeFunctionTolerance(value));
        }
        self.relative_function_tolerance = value;
        Ok(self)
    }

    /// Sets the infinity-norm tolerance for the projected gradient.
    ///
    /// Set this to zero to disable convergence through this test.
    pub fn with_projected_gradient_tolerance(mut self, value: f64) -> Result<Self, OptionsError> {
        if !value.is_finite() || value < 0.0 {
            return Err(OptionsError::InvalidProjectedGradientTolerance(value));
        }
        self.projected_gradient_tolerance = value;
        Ok(self)
    }

    /// Returns the number of correction pairs retained by the solver.
    pub fn history_size(self) -> usize {
        self.history_size
    }

    /// Returns the exact maximum number of accepted major iterations.
    pub fn max_iterations(self) -> usize {
        self.max_iterations
    }

    /// Returns the maximum number of objective/gradient evaluations.
    pub fn max_evaluations(self) -> usize {
        self.max_evaluations
    }

    /// Returns the per-line-search objective-evaluation limit.
    pub fn max_line_search_evaluations(self) -> usize {
        self.max_line_search_evaluations
    }

    /// Returns the relative accepted-iterate function-reduction tolerance.
    pub fn relative_function_tolerance(self) -> f64 {
        self.relative_function_tolerance
    }

    /// Returns the projected-gradient infinity-norm tolerance.
    pub fn projected_gradient_tolerance(self) -> f64 {
        self.projected_gradient_tolerance
    }

    /// Returns the configured dense-history backend policy.
    pub fn backend(self) -> Backend {
        self.backend
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_zero_limits() {
        assert_eq!(
            Options::new().with_max_evaluations(0),
            Err(OptionsError::ZeroEvaluations)
        );
        assert_eq!(
            Options::new().with_max_line_search_evaluations(0),
            Err(OptionsError::ZeroLineSearchEvaluations)
        );
    }

    #[test]
    fn uses_direct_tolerance_values() {
        let options = Options::new()
            .with_relative_function_tolerance(2.5e-9)
            .unwrap()
            .with_projected_gradient_tolerance(4.0e-7)
            .unwrap();
        assert_eq!(options.relative_function_tolerance(), 2.5e-9);
        assert_eq!(options.projected_gradient_tolerance(), 4.0e-7);
    }

    #[test]
    fn rejects_excessive_history() {
        assert!(matches!(
            Options::new().with_history_size(usize::MAX),
            Err(OptionsError::HistoryTooLarge { .. })
        ));
    }

    #[cfg(not(feature = "faer-backend"))]
    #[test]
    fn rejects_unavailable_backend() {
        assert_eq!(
            Options::new().with_backend(Backend::Faer),
            Err(OptionsError::UnsupportedBackend)
        );
    }
}

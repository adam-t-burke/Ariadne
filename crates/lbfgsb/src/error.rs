//! Error types returned while configuring and running a solve.

use std::fmt;

/// An invalid bound specification.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum BoundsError {
    /// A solve was requested with no variables.
    Empty,
    /// Variable, lower-bound, and upper-bound lengths differ.
    LengthMismatch {
        /// Number of optimization variables.
        variables: usize,
        /// Number of lower endpoints supplied.
        lower: usize,
        /// Number of upper endpoints supplied.
        upper: usize,
    },
    /// A bound endpoint is NaN.
    NotANumber {
        /// Index containing the invalid endpoint.
        index: usize,
    },
    /// A lower endpoint was positive infinity.
    PositiveInfiniteLower {
        /// Index containing the invalid endpoint.
        index: usize,
    },
    /// An upper endpoint was negative infinity.
    NegativeInfiniteUpper {
        /// Index containing the invalid endpoint.
        index: usize,
    },
    /// A fixed bound used an infinite coordinate.
    InfiniteFixed {
        /// Index containing the invalid fixed bound.
        index: usize,
        /// Infinite fixed value.
        value: f64,
    },
    /// A lower endpoint exceeds its upper endpoint.
    Inverted {
        /// Index containing the inverted interval.
        index: usize,
        /// Invalid lower endpoint.
        lower: f64,
        /// Invalid upper endpoint.
        upper: f64,
    },
    /// A variable coordinate was not finite.
    NonFiniteVariable {
        /// Index containing the invalid coordinate.
        index: usize,
    },
    /// A gradient component was not finite.
    NonFiniteGradient {
        /// Index containing the invalid component.
        index: usize,
    },
}

impl fmt::Display for BoundsError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("at least one variable is required"),
            Self::LengthMismatch {
                variables,
                lower,
                upper,
            } => write!(
                f,
                "bound lengths ({lower}, {upper}) do not match {variables} variables"
            ),
            Self::NotANumber { index } => write!(f, "bound {index} is NaN"),
            Self::PositiveInfiniteLower { index } => {
                write!(f, "lower bound {index} cannot be positive infinity")
            }
            Self::NegativeInfiniteUpper { index } => {
                write!(f, "upper bound {index} cannot be negative infinity")
            }
            Self::InfiniteFixed { index, value } => {
                write!(f, "fixed bound {index} cannot be infinite ({value})")
            }
            Self::Inverted {
                index,
                lower,
                upper,
            } => write!(
                f,
                "lower bound {lower} exceeds upper bound {upper} at index {index}"
            ),
            Self::NonFiniteVariable { index } => {
                write!(f, "variable {index} is not finite")
            }
            Self::NonFiniteGradient { index } => {
                write!(f, "gradient {index} is not finite")
            }
        }
    }
}

impl std::error::Error for BoundsError {}

/// An invalid solver option.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum OptionsError {
    /// History size was zero.
    ZeroHistory,
    /// History size exceeded the supported maximum.
    HistoryTooLarge {
        /// Requested history size.
        requested: usize,
        /// Largest supported history size.
        maximum: usize,
    },
    /// Maximum iterations was zero.
    ZeroIterations,
    /// Maximum evaluations was zero.
    ZeroEvaluations,
    /// Maximum line-search evaluations was zero.
    ZeroLineSearchEvaluations,
    /// Relative-function tolerance was negative or non-finite.
    InvalidRelativeFunctionTolerance(f64),
    /// Projected-gradient tolerance was negative or non-finite.
    InvalidProjectedGradientTolerance(f64),
    /// The selected backend is unavailable in this build.
    UnsupportedBackend,
}

impl fmt::Display for OptionsError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ZeroHistory => f.write_str("history size must be positive"),
            Self::HistoryTooLarge { requested, maximum } => write!(
                f,
                "history size {requested} exceeds the supported maximum of {maximum}"
            ),
            Self::ZeroIterations => f.write_str("maximum iterations must be positive"),
            Self::ZeroEvaluations => f.write_str("maximum evaluations must be positive"),
            Self::ZeroLineSearchEvaluations => {
                f.write_str("maximum line-search evaluations must be positive")
            }
            Self::InvalidRelativeFunctionTolerance(value) => write!(
                f,
                "relative function tolerance must be finite and non-negative, got {value}"
            ),
            Self::InvalidProjectedGradientTolerance(value) => write!(
                f,
                "projected gradient tolerance must be finite and non-negative, got {value}"
            ),
            Self::UnsupportedBackend => {
                f.write_str("the requested backend is unavailable in this build")
            }
        }
    }
}

impl std::error::Error for OptionsError {}

/// A workspace could not be sized safely for a solve.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum WorkspaceError {
    /// Required element counts overflowed `usize`.
    SizeOverflow {
        /// Number of variables requested.
        variables: usize,
        /// Number of correction pairs requested.
        history_size: usize,
    },
    /// Reserving a solver-owned buffer failed.
    AllocationFailed {
        /// Buffer that could not be reserved.
        buffer: &'static str,
        /// Required number of elements.
        elements: usize,
    },
}

impl fmt::Display for WorkspaceError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::SizeOverflow {
                variables,
                history_size,
            } => write!(
                f,
                "workspace size overflow for {variables} variables and history size {history_size}"
            ),
            Self::AllocationFailed { buffer, elements } => {
                write!(
                    f,
                    "could not reserve {elements} elements for workspace buffer {buffer}"
                )
            }
        }
    }
}

impl std::error::Error for WorkspaceError {}

/// A value/gradient evaluation violated the callback contract.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum EvaluationError {
    /// The objective returned NaN or infinity.
    NonFiniteValue,
    /// A gradient component was NaN or infinity.
    NonFiniteGradient {
        /// Index of the invalid gradient component.
        index: usize,
    },
}

impl fmt::Display for EvaluationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NonFiniteValue => f.write_str("objective returned a non-finite value"),
            Self::NonFiniteGradient { index } => {
                write!(
                    f,
                    "objective returned a non-finite gradient at index {index}"
                )
            }
        }
    }
}

impl std::error::Error for EvaluationError {}

/// An error that prevents the solver from producing a report.
#[derive(Debug)]
#[non_exhaustive]
pub enum SolveError<E> {
    /// Bound validation failed.
    Bounds(BoundsError),
    /// Option validation failed.
    Options(OptionsError),
    /// Reusable workspace sizing or allocation failed.
    Workspace(WorkspaceError),
    /// An initial coordinate was NaN or infinity.
    InvalidInitialValue {
        /// Index of the invalid initial coordinate.
        index: usize,
    },
    /// The objective violated the finite value/gradient contract.
    Evaluation(EvaluationError),
    /// The objective callback returned an application error.
    Objective(E),
    /// The accepted-iteration callback returned an application error.
    Callback(E),
    /// An invariant failed inside the numerical engine.
    Internal(String),
}

impl<E: fmt::Display> fmt::Display for SolveError<E> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Bounds(error) => error.fmt(f),
            Self::Options(error) => error.fmt(f),
            Self::Workspace(error) => error.fmt(f),
            Self::InvalidInitialValue { index } => {
                write!(f, "initial value at index {index} is not finite")
            }
            Self::Evaluation(error) => error.fmt(f),
            Self::Objective(error) => write!(f, "objective evaluation failed: {error}"),
            Self::Callback(error) => write!(f, "iteration callback failed: {error}"),
            Self::Internal(message) => write!(f, "solver failed internally: {message}"),
        }
    }
}

impl<E: std::error::Error + 'static> std::error::Error for SolveError<E> {}

impl<E> From<BoundsError> for SolveError<E> {
    fn from(value: BoundsError) -> Self {
        Self::Bounds(value)
    }
}

impl<E> From<OptionsError> for SolveError<E> {
    fn from(value: OptionsError) -> Self {
        Self::Options(value)
    }
}

impl<E> From<WorkspaceError> for SolveError<E> {
    fn from(value: WorkspaceError) -> Self {
        Self::Workspace(value)
    }
}

//! Standalone, callback-first L-BFGS-B interface.
//!
//! The public API uses `f64` slices throughout. The objective writes its
//! gradient into solver-provided storage, and [`Solver`] retains a [`Workspace`]
//! across solves.
//!
//! ```
//! use ariadne_lbfgsb::{Bounds, Options, Solver};
//! use std::convert::Infallible;
//!
//! let mut x = [3.0, -4.0];
//! let bounds = Bounds::new(&[-10.0; 2], &[10.0; 2], x.len())?;
//! let mut solver = Solver::new(
//!     Options::new().with_projected_gradient_tolerance(1e-8)?
//! );
//! let report = solver.minimize(&mut x, bounds, |x, gradient| {
//!     gradient[0] = 2.0 * (x[0] - 1.0);
//!     gradient[1] = 2.0 * (x[1] + 2.0);
//!     Ok::<_, Infallible>((x[0] - 1.0).powi(2) + (x[1] + 2.0).powi(2))
//! })?;
//! assert!(report.termination.converged());
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! The numerical engine is a safe scalar implementation of L-BFGS-B 3.0.
//!
//! # Stopping tests
//!
//! The projected-gradient test uses the infinity norm of the gradient
//! projected onto feasible directions. The relative-function test compares
//! consecutive accepted iterates using
//! `(f[k] - f[k+1]) / max(|f[k]|, |f[k+1]|, 1)`. Initial coordinates are
//! projected into the supplied box before the first objective evaluation.
#![warn(missing_docs)]

mod bounds;
mod core;
mod error;
mod kernel;
mod options;
mod report;
mod session;
mod solver;
mod workspace;

pub use bounds::{projected_gradient_norm, Bounds};
pub use error::{BoundsError, EvaluationError, OptionsError, SolveError, WorkspaceError};
pub use options::{Backend, Options};
#[cfg(feature = "benchmark-instrumentation")]
pub use report::BenchmarkTimings;
pub use report::{
    Control, Convergence, Failure, Iteration, Report, Stats, StopReason, Termination,
};
pub use solver::{SolveAdapter, Solver};
pub use workspace::Workspace;

//! **Theseus** — Force Density Method form-finding solver with hand-coded adjoints.
//!
//! This crate implements the complete FDM optimisation pipeline:
//!
//! 1. **Forward solve** (`fdm`): assemble A(q), factorise, triangular solve.
//! 2. **Objectives** (`objectives`): 13 loss functions on geometry / forces / reactions.
//! 3. **Gradients** (`gradients`): hand-coded adjoint + explicit derivatives.
//! 4. **Optimiser** (`optimizer`): L-BFGS-B via `lbfgsb-rs-pure`.
//! 5. **FFI** (`ffi`): C-compatible API for Grasshopper / C# P/Invoke.
//!
//! All public functions return `Result<_, TheseusError>` — the crate never
//! panics in normal operation.

pub mod types;
pub mod fdm;
pub mod objectives;
pub mod gradients;
pub mod optimizer;
pub mod ffi;

pub use types::TheseusError;
pub use types::ObjectiveTrait;

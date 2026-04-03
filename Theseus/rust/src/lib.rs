//! **Theseus** — Force Density Method form-finding solver with hand-coded adjoints.
//!
//! This crate implements the complete FDM optimisation pipeline:
//!
//! 1. **Forward solve** (`fdm`): assemble A(q), factorise, triangular solve.
//! 2. **Objectives** (`objectives`): 13 loss functions on geometry / forces / reactions.
//! 3. **Gradients** (`gradients`): hand-coded adjoint + explicit derivatives.
//! 4. **Optimiser** (`optimizer`): L-BFGS via `argmin`.
//! 5. **FFI** (`ffi`): C-compatible API for Grasshopper / C# P/Invoke.
//!
//! All public functions return `Result<_, TheseusError>` — the crate never
//! panics in normal operation.

pub mod sparse;
pub mod types;
pub mod fdm;
pub mod objectives;
pub mod gradients;
pub mod optimizer;
pub mod inverse;
pub mod ffi;

// FEA modules (Phase 1: bar/truss elements)
pub mod fea_types;
pub mod fea_assembly;
pub mod fea_solve;
pub mod fea_objectives;
pub mod fea_gradients;
pub mod fea_optimizer;
pub mod fea_ffi;

// Solid element modules (Phase 2: tet4 continuum elements)
pub mod solid_types;
pub mod solid_assembly;
pub mod solid_solve;
pub mod solid_ffi;
pub mod solid_spr;
pub mod solid_spr_ffi;

// Shell element modules (Phase 3: tri/quad shell elements)
pub mod shell_types;
pub mod shell_assembly;
pub mod shell_solve;
pub mod shell_ffi;
pub mod shell_spr;
pub mod shell_spr_ffi;

pub use types::TheseusError;
pub use types::ObjectiveTrait;
pub use fea_types::FeaObjectiveTrait;

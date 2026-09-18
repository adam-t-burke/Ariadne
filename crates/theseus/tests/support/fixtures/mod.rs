//! Regression and benchmark fixtures (§4.4 of `ITERATIVE_SOLVER_PROGRAM.md`).
//!
//! Every generator returns a [`network::Network`] (positions, edges, fixed
//! nodes) and has a `make_*_problem` companion returning a Theseus
//! [`theseus::types::Problem`] whose target is the equilibrium shape of a
//! smooth interior force-density field, mirroring
//! `grid::make_recoverable_grid_problem`, so an optimiser started at `q = 1`
//! has a reachable zero-loss minimum inside the box bounds.
//!
//! | generator | topology | supports |
//! |---|---|---|
//! | `grid::make_recoverable_grid_problem(n)` | `n × n` grid | 4 corners |
//! | `irregular::irregular_mesh(side, seed)` | jittered lattice, 5–7 nearest neighbours | boundary ring |
//! | `dome::cable_dome(rings, spokes)` | hub, rings, spokes, X-bracing | outer ring |
//! | `few_supports::few_supports(side)` | grid | 2 opposite corners |
//! | `disconnected::disconnected(side, k)` | `k` separate grids | 4 corners each |
//! | `anisotropic::anisotropic_q(side, ratio)` | grid, `q*` spanning `ratio` | 4 corners |
//!
//! Include from a test or example with
//! `#[path = "support/fixtures/mod.rs"] mod fixtures;` (the harness helpers in
//! [`harness`] need `serde_json`, a dev-dependency).

#![allow(dead_code, unused_imports)]

#[path = "../grid.rs"]
pub mod grid;

pub mod anisotropic;
pub mod disconnected;
pub mod dome;
pub mod few_supports;
pub mod harness;
pub mod irregular;
pub mod network;
pub mod rng;

pub use anisotropic::{anisotropic_q, anisotropic_q_star, make_anisotropic_q_problem};
pub use disconnected::{disconnected, make_disconnected_problem};
pub use dome::{cable_dome, make_cable_dome_problem};
pub use few_supports::{few_supports, few_supports_with, make_few_supports_problem};
pub use irregular::{irregular_mesh, make_irregular_mesh_problem};
pub use network::{grid_network, make_recoverable, smooth_q_star, Network};

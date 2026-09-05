//! Experimental 3D projective geometric algebra, `G(3,0,1)`.
//!
//! This crate is a sandbox for an Ariadne PGA branch. It is **not** linked
//! into Theseus. The FDM solve, adjoints, L-BFGS-B charts, and NURBS supports
//! stay Euclidean; PGA is a candidate language for *flats* and *incidence*:
//! points, lines, planes, meet/join, and Euclidean motors.
//!
//! # Why this algebra
//!
//! Plane-based PGA (`P(R*_{3,0,1})`) is the smallest geometric algebra that
//! treats Euclidean points, lines, and planes uniformly, including infinity.
//! That matches Ariadne's variable supports more closely than VGA (no points)
//! or CGA (spheres and circles, 32-D).
//!
//! # What it should not replace
//!
//! The force-density linear solve `A(q) x = b`, hydrostatic Newton loops, and
//! NURBS evaluation are not incidence problems. Those stay in `theseus` and
//! `nurbsbook`.
//!
//! # Dual conventions
//!
//! Meet is the outer product. Join is the Poincaré-dual regressive product,
//! not the metric dual. The degenerate `e0² = 0` metric makes that distinction
//! load-bearing.

pub mod ariadne;
pub mod multivector;
pub mod types;

pub use ariadne::{
    constraint_plane, planar_constraint_along_direction, plane_incidence, rail_incidence,
    rail_line, roller_flat, SupportFlatKind,
};
pub use multivector::MultiVector;
pub use types::{Line, Motor, Plane, Point};

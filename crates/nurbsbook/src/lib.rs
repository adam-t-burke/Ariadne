//! Core NURBS algorithms translated from Piegl & Tiller, The NURBS Book.
//!
//! This crate is intentionally standalone so it can be split out of Ariadne
//! later. The initial API focuses on the algorithms needed by Theseus variable
//! supports: basis functions, curve/surface evaluation, and first derivatives.

pub mod basis;
pub mod curves;
pub mod rational;
pub mod surfaces;
pub mod types;
pub mod utils;

pub use rational::{nurbs_curve_derivatives, nurbs_curve_point};
pub use rational::{nurbs_surface_derivatives, nurbs_surface_partials, nurbs_surface_point};
pub use types::{BSplineCurve, BSplineSurface, KnotVector, NurbsCurve, NurbsError, NurbsSurface};

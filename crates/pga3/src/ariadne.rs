//! Mapping from Ariadne / Theseus geometry onto PGA flats.
//!
//! This module is intentionally *not* wired into `theseus`. It exists to make
//! the correspondence testable and to keep the experimental branch honest about
//! what PGA can and cannot replace.

use crate::types::{Line, Plane, Point};

/// Theseus variable-support kinds, classified as PGA flats plus leftover charts.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum SupportFlatKind {
    /// Fully fixed node: a Euclidean point.
    Point,
    /// One free Euclidean direction: rail or single-axis roller.
    Line,
    /// Two free Euclidean directions: planar roller.
    Plane,
    /// Unbounded 3-space, spherical box, or NURBS. Not a PGA flat.
    NotAFlat,
}

/// Classify a Theseus roller from its enabled axes.
pub fn roller_flat(enabled: [bool; 3]) -> SupportFlatKind {
    match enabled.iter().filter(|&&b| b).count() {
        0 => SupportFlatKind::Point,
        1 => SupportFlatKind::Line,
        2 => SupportFlatKind::Plane,
        _ => SupportFlatKind::NotAFlat,
    }
}

/// Rail supporting line. The open-interval parameter `t ∈ (0, 1)` stays an
/// optimizer chart; PGA only sees the infinite line.
pub fn rail_line(start: [f64; 3], end: [f64; 3]) -> Line {
    Line::from_segment(start, end)
}

/// Plane used by `PlanarConstraintAlongDirection` (origin + in-plane axes).
pub fn constraint_plane(origin: [f64; 3], x_axis: [f64; 3], y_axis: [f64; 3]) -> Plane {
    Plane::from_origin_axes(origin, x_axis, y_axis)
}

/// Euclidean incidence residual of a node on a rail line.
pub fn rail_incidence(node: [f64; 3], start: [f64; 3], end: [f64; 3]) -> f64 {
    rail_line(start, end).distance(Point::from_xyz(node))
}

/// Euclidean incidence residual of a node on a constraint plane.
pub fn plane_incidence(
    node: [f64; 3],
    origin: [f64; 3],
    x_axis: [f64; 3],
    y_axis: [f64; 3],
) -> f64 {
    constraint_plane(origin, x_axis, y_axis).signed_distance(Point::from_xyz(node))
}

/// Theseus `PlanarConstraintAlongDirection` parameter
/// `t = n·(O−P) / (n·d)`, which is *not* PGA incidence unless `d` is the
/// plane normal (and then only up to scaling).
pub fn planar_constraint_along_direction(
    node: [f64; 3],
    origin: [f64; 3],
    x_axis: [f64; 3],
    y_axis: [f64; 3],
    direction: [f64; 3],
) -> Option<f64> {
    let n = [
        x_axis[1] * y_axis[2] - x_axis[2] * y_axis[1],
        x_axis[2] * y_axis[0] - x_axis[0] * y_axis[2],
        x_axis[0] * y_axis[1] - x_axis[1] * y_axis[0],
    ];
    let n_dot_d = n[0] * direction[0] + n[1] * direction[1] + n[2] * direction[2];
    if n_dot_d.abs() < 1e-12 {
        return None;
    }
    let n_dot_om_p =
        n[0] * (origin[0] - node[0]) + n[1] * (origin[1] - node[1]) + n[2] * (origin[2] - node[2]);
    Some(n_dot_om_p / n_dot_d)
}

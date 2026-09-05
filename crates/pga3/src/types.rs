//! Typed Euclidean flats in plane-based PGA.
//!
//! Embedding follows Klein / Gunn plane-based PGA:
//! - plane `ax + by + cz + d = 0` is `a e1 + b e2 + c e3 + d e0`
//! - Euclidean point `(x, y, z)` is `e123 + x e032 + y e013 + z e021`
//! - line is the corresponding simple bivector (Plücker direction + moment)

use crate::multivector::{
    MultiVector, E0, E01, E012, E013, E02, E023, E03, E1, E12, E123, E13, E2, E23, E3,
};

const EPS: f64 = 1e-14;

/// Homogeneous Euclidean point. Finite points have `w = 1`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Point {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub w: f64,
}

impl Point {
    pub fn euclidean(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z, w: 1.0 }
    }

    pub fn from_xyz(xyz: [f64; 3]) -> Self {
        Self::euclidean(xyz[0], xyz[1], xyz[2])
    }

    pub fn xyz(self) -> Option<[f64; 3]> {
        if self.w.abs() <= EPS {
            None
        } else {
            Some([self.x / self.w, self.y / self.w, self.z / self.w])
        }
    }

    pub fn to_mv(self) -> MultiVector {
        // Klein: w e123 + x e032 + y e013 + z e021
        // e032 = -e023, e021 = -e012
        let mut mv = MultiVector::zero();
        mv.c[E123] = self.w;
        mv.c[E023] = -self.x;
        mv.c[E013] = self.y;
        mv.c[E012] = -self.z;
        mv
    }

    pub fn from_mv(mv: MultiVector) -> Self {
        Self {
            x: -mv.c[E023],
            y: mv.c[E013],
            z: -mv.c[E012],
            w: mv.c[E123],
        }
    }

    pub fn join(self, other: Self) -> Line {
        Line::from_mv(self.to_mv().regressive(other.to_mv()))
    }

    pub fn join_line(self, line: Line) -> Plane {
        Plane::from_mv(self.to_mv().regressive(line.to_mv()))
    }
}

/// Plane `a x + b y + c z + d = 0`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Plane {
    pub a: f64,
    pub b: f64,
    pub c: f64,
    pub d: f64,
}

impl Plane {
    pub fn new(a: f64, b: f64, c: f64, d: f64) -> Self {
        Self { a, b, c, d }
    }

    pub fn from_point_normal(origin: [f64; 3], normal: [f64; 3]) -> Self {
        let d = -(normal[0] * origin[0] + normal[1] * origin[1] + normal[2] * origin[2]);
        Self::new(normal[0], normal[1], normal[2], d)
    }

    /// Theseus `PlanarConstraintAlongDirection` plane from origin + in-plane axes.
    pub fn from_origin_axes(origin: [f64; 3], x_axis: [f64; 3], y_axis: [f64; 3]) -> Self {
        let n = [
            x_axis[1] * y_axis[2] - x_axis[2] * y_axis[1],
            x_axis[2] * y_axis[0] - x_axis[0] * y_axis[2],
            x_axis[0] * y_axis[1] - x_axis[1] * y_axis[0],
        ];
        Self::from_point_normal(origin, n)
    }

    pub fn through(p: Point, q: Point, r: Point) -> Self {
        p.join(q).join_point(r)
    }

    pub fn normal(self) -> [f64; 3] {
        [self.a, self.b, self.c]
    }

    pub fn euclidean_norm(self) -> f64 {
        (self.a * self.a + self.b * self.b + self.c * self.c).sqrt()
    }

    pub fn to_mv(self) -> MultiVector {
        let mut mv = MultiVector::zero();
        mv.c[E0] = self.d;
        mv.c[E1] = self.a;
        mv.c[E2] = self.b;
        mv.c[E3] = self.c;
        mv
    }

    pub fn from_mv(mv: MultiVector) -> Self {
        Self {
            a: mv.c[E1],
            b: mv.c[E2],
            c: mv.c[E3],
            d: mv.c[E0],
        }
    }

    pub fn meet(self, other: Self) -> Line {
        Line::from_mv(self.to_mv().outer(other.to_mv()))
    }

    pub fn meet_line(self, line: Line) -> Point {
        Point::from_mv(self.to_mv().outer(line.to_mv()))
    }

    pub fn signed_distance(self, point: Point) -> f64 {
        let p = point.xyz().unwrap_or([point.x, point.y, point.z]);
        let num = self.a * p[0] + self.b * p[1] + self.c * p[2] + self.d;
        let den = self.euclidean_norm();
        if den <= EPS {
            num
        } else {
            num / den
        }
    }
}

/// Plücker line: `moment = p × direction` for any point `p` on the line.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Line {
    pub direction: [f64; 3],
    pub moment: [f64; 3],
}

impl Line {
    pub fn through(p: Point, q: Point) -> Self {
        p.join(q)
    }

    pub fn from_segment(start: [f64; 3], end: [f64; 3]) -> Self {
        Self::through(Point::from_xyz(start), Point::from_xyz(end))
    }

    pub fn join_point(self, point: Point) -> Plane {
        point.join_line(self)
    }

    pub fn to_mv(self) -> MultiVector {
        // Klein: vx e23 + vy e31 + vz e12 + mx e01 + my e02 + mz e03
        // e31 = -e13
        let mut mv = MultiVector::zero();
        mv.c[E23] = self.direction[0];
        mv.c[E13] = -self.direction[1];
        mv.c[E12] = self.direction[2];
        mv.c[E01] = self.moment[0];
        mv.c[E02] = self.moment[1];
        mv.c[E03] = self.moment[2];
        mv
    }

    pub fn from_mv(mv: MultiVector) -> Self {
        Self {
            direction: [mv.c[E23], -mv.c[E13], mv.c[E12]],
            moment: [mv.c[E01], mv.c[E02], mv.c[E03]],
        }
    }

    pub fn direction_norm(self) -> f64 {
        let d = self.direction;
        (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
    }

    pub fn distance(self, point: Point) -> f64 {
        let p = point.xyz().unwrap_or([point.x, point.y, point.z]);
        let qx = p[1] * self.direction[2] - p[2] * self.direction[1] - self.moment[0];
        let qy = p[2] * self.direction[0] - p[0] * self.direction[2] - self.moment[1];
        let qz = p[0] * self.direction[1] - p[1] * self.direction[0] - self.moment[2];
        let num = (qx * qx + qy * qy + qz * qz).sqrt();
        let den = self.direction_norm();
        if den <= EPS {
            num
        } else {
            num / den
        }
    }

    /// Parameter `t` of the closest point `start + t (end - start)` on a segment chart.
    pub fn segment_parameter(self, point: Point, start: [f64; 3], end: [f64; 3]) -> Option<f64> {
        let p = point.xyz()?;
        let d = [end[0] - start[0], end[1] - start[1], end[2] - start[2]];
        let len2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
        if len2 <= EPS {
            return None;
        }
        let rel = [p[0] - start[0], p[1] - start[1], p[2] - start[2]];
        Some((rel[0] * d[0] + rel[1] * d[1] + rel[2] * d[2]) / len2)
    }
}

/// Unit even versor (dual quaternion): translator and rotator.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Motor {
    mv: MultiVector,
}

impl Motor {
    pub fn identity() -> Self {
        Self {
            mv: MultiVector::scalar(1.0),
        }
    }

    /// Translator `1 - (1/2)(tx e01 + ty e02 + tz e03)`.
    pub fn translator(tx: f64, ty: f64, tz: f64) -> Self {
        let mut mv = MultiVector::scalar(1.0);
        mv.c[E01] = -0.5 * tx;
        mv.c[E02] = -0.5 * ty;
        mv.c[E03] = -0.5 * tz;
        Self { mv }
    }

    /// Rotor about a unit axis through the origin, angle in radians.
    pub fn rotator(axis: [f64; 3], angle: f64) -> Self {
        let half = 0.5 * angle;
        let s = half.sin();
        let mut mv = MultiVector::scalar(half.cos());
        // Bivector plane of rotation is dual to the axis: rx e23 + ry e31 + rz e12
        mv.c[E23] = -s * axis[0];
        mv.c[E13] = s * axis[1];
        mv.c[E12] = -s * axis[2];
        Self { mv }
    }

    pub fn apply_point(self, point: Point) -> Point {
        Point::from_mv(self.mv.sandwich(point.to_mv()))
    }

    pub fn apply_plane(self, plane: Plane) -> Plane {
        Plane::from_mv(self.mv.sandwich(plane.to_mv()))
    }

    pub fn apply_line(self, line: Line) -> Line {
        Line::from_mv(self.mv.sandwich(line.to_mv()))
    }
}

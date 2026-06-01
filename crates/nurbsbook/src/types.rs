use std::error::Error;
use std::fmt::{Display, Formatter};

#[derive(Debug, Clone, PartialEq)]
pub enum NurbsError {
    InvalidDegree,
    InvalidKnotVector(String),
    InvalidControlPoints(String),
    ZeroWeight,
    ParameterOutOfDomain,
}

impl Display for NurbsError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidDegree => write!(f, "degree must be at least one"),
            Self::InvalidKnotVector(msg) => write!(f, "invalid knot vector: {msg}"),
            Self::InvalidControlPoints(msg) => write!(f, "invalid control points: {msg}"),
            Self::ZeroWeight => write!(f, "rational denominator is zero"),
            Self::ParameterOutOfDomain => write!(f, "parameter is outside the spline domain"),
        }
    }
}

impl Error for NurbsError {}

#[derive(Debug, Clone, PartialEq)]
pub struct KnotVector {
    knots: Vec<f64>,
}

impl KnotVector {
    pub fn new(knots: Vec<f64>) -> Result<Self, NurbsError> {
        if knots.len() < 2 {
            return Err(NurbsError::InvalidKnotVector(
                "a knot vector needs at least two values".into(),
            ));
        }
        for (i, pair) in knots.windows(2).enumerate() {
            if !pair[0].is_finite() || !pair[1].is_finite() {
                return Err(NurbsError::InvalidKnotVector(
                    "knot values must be finite".into(),
                ));
            }
            if pair[1] < pair[0] {
                return Err(NurbsError::InvalidKnotVector(format!(
                    "values must be nondecreasing at index {i}"
                )));
            }
        }
        Ok(Self { knots })
    }

    pub fn unchecked(knots: Vec<f64>) -> Self {
        Self { knots }
    }

    pub fn as_slice(&self) -> &[f64] {
        &self.knots
    }

    pub fn len(&self) -> usize {
        self.knots.len()
    }

    pub fn is_empty(&self) -> bool {
        self.knots.is_empty()
    }

    pub fn domain(&self, degree: usize, num_control_points: usize) -> Result<[f64; 2], NurbsError> {
        if num_control_points < degree + 1 {
            return Err(NurbsError::InvalidControlPoints(
                "not enough control points for degree".into(),
            ));
        }
        let expected = num_control_points + degree + 1;
        if self.knots.len() != expected {
            return Err(NurbsError::InvalidKnotVector(format!(
                "got {}, expected {expected}",
                self.knots.len()
            )));
        }
        Ok([self.knots[degree], self.knots[num_control_points]])
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct BSplineCurve {
    pub degree: usize,
    pub knots: KnotVector,
    pub control_points: Vec<[f64; 3]>,
}

impl BSplineCurve {
    pub fn new(
        degree: usize,
        knots: Vec<f64>,
        control_points: Vec<[f64; 3]>,
    ) -> Result<Self, NurbsError> {
        validate_degree_and_points(degree, control_points.len())?;
        let knots = KnotVector::new(knots)?;
        knots.domain(degree, control_points.len())?;
        Ok(Self {
            degree,
            knots,
            control_points,
        })
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct NurbsCurve {
    pub degree: usize,
    pub knots: KnotVector,
    /// Homogeneous control points: `[x*w, y*w, z*w, w]`.
    pub control_points: Vec<[f64; 4]>,
}

impl NurbsCurve {
    pub fn new(
        degree: usize,
        knots: Vec<f64>,
        control_points: Vec<[f64; 4]>,
    ) -> Result<Self, NurbsError> {
        validate_degree_and_points(degree, control_points.len())?;
        for p in &control_points {
            if p[3] == 0.0 {
                return Err(NurbsError::InvalidControlPoints(
                    "weights must be nonzero".into(),
                ));
            }
        }
        let knots = KnotVector::new(knots)?;
        knots.domain(degree, control_points.len())?;
        Ok(Self {
            degree,
            knots,
            control_points,
        })
    }

    pub fn from_cartesian(
        degree: usize,
        knots: Vec<f64>,
        points: Vec<[f64; 3]>,
        weights: Vec<f64>,
    ) -> Result<Self, NurbsError> {
        if points.len() != weights.len() {
            return Err(NurbsError::InvalidControlPoints(
                "points and weights length mismatch".into(),
            ));
        }
        let control_points = points
            .into_iter()
            .zip(weights)
            .map(|(p, w)| [p[0] * w, p[1] * w, p[2] * w, w])
            .collect();
        Self::new(degree, knots, control_points)
    }

    pub fn domain(&self) -> [f64; 2] {
        self.knots
            .domain(self.degree, self.control_points.len())
            .expect("NurbsCurve invariants are checked at construction")
    }

    pub fn point(&self, u: f64) -> Result<[f64; 3], NurbsError> {
        crate::rational::nurbs_curve_point(self, u)
    }

    pub fn derivatives(&self, u: f64, d: usize) -> Result<Vec<[f64; 3]>, NurbsError> {
        crate::rational::nurbs_curve_derivatives(self, u, d)
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct BSplineSurface {
    pub degree_u: usize,
    pub degree_v: usize,
    pub knots_u: KnotVector,
    pub knots_v: KnotVector,
    pub num_u: usize,
    pub num_v: usize,
    pub control_points: Vec<[f64; 3]>,
}

impl BSplineSurface {
    pub fn new(
        degree_u: usize,
        degree_v: usize,
        knots_u: Vec<f64>,
        knots_v: Vec<f64>,
        num_u: usize,
        num_v: usize,
        control_points: Vec<[f64; 3]>,
    ) -> Result<Self, NurbsError> {
        validate_surface(degree_u, degree_v, num_u, num_v, control_points.len())?;
        let knots_u = KnotVector::new(knots_u)?;
        let knots_v = KnotVector::new(knots_v)?;
        knots_u.domain(degree_u, num_u)?;
        knots_v.domain(degree_v, num_v)?;
        Ok(Self {
            degree_u,
            degree_v,
            knots_u,
            knots_v,
            num_u,
            num_v,
            control_points,
        })
    }

    pub fn cp(&self, u: usize, v: usize) -> [f64; 3] {
        self.control_points[u * self.num_v + v]
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct NurbsSurface {
    pub degree_u: usize,
    pub degree_v: usize,
    pub knots_u: KnotVector,
    pub knots_v: KnotVector,
    pub num_u: usize,
    pub num_v: usize,
    /// Row-major homogeneous control points: `u * num_v + v`.
    pub control_points: Vec<[f64; 4]>,
}

impl NurbsSurface {
    pub fn new(
        degree_u: usize,
        degree_v: usize,
        knots_u: Vec<f64>,
        knots_v: Vec<f64>,
        num_u: usize,
        num_v: usize,
        control_points: Vec<[f64; 4]>,
    ) -> Result<Self, NurbsError> {
        validate_surface(degree_u, degree_v, num_u, num_v, control_points.len())?;
        for p in &control_points {
            if p[3] == 0.0 {
                return Err(NurbsError::InvalidControlPoints(
                    "weights must be nonzero".into(),
                ));
            }
        }
        let knots_u = KnotVector::new(knots_u)?;
        let knots_v = KnotVector::new(knots_v)?;
        knots_u.domain(degree_u, num_u)?;
        knots_v.domain(degree_v, num_v)?;
        Ok(Self {
            degree_u,
            degree_v,
            knots_u,
            knots_v,
            num_u,
            num_v,
            control_points,
        })
    }

    pub fn from_cartesian(
        degree_u: usize,
        degree_v: usize,
        knots_u: Vec<f64>,
        knots_v: Vec<f64>,
        num_u: usize,
        num_v: usize,
        points: Vec<[f64; 3]>,
        weights: Vec<f64>,
    ) -> Result<Self, NurbsError> {
        if points.len() != weights.len() {
            return Err(NurbsError::InvalidControlPoints(
                "points and weights length mismatch".into(),
            ));
        }
        let control_points = points
            .into_iter()
            .zip(weights)
            .map(|(p, w)| [p[0] * w, p[1] * w, p[2] * w, w])
            .collect();
        Self::new(
            degree_u,
            degree_v,
            knots_u,
            knots_v,
            num_u,
            num_v,
            control_points,
        )
    }

    pub fn domain_u(&self) -> [f64; 2] {
        self.knots_u
            .domain(self.degree_u, self.num_u)
            .expect("NurbsSurface invariants are checked at construction")
    }

    pub fn domain_v(&self) -> [f64; 2] {
        self.knots_v
            .domain(self.degree_v, self.num_v)
            .expect("NurbsSurface invariants are checked at construction")
    }

    pub fn cp(&self, u: usize, v: usize) -> [f64; 4] {
        self.control_points[u * self.num_v + v]
    }

    pub fn point(&self, u: f64, v: f64) -> Result<[f64; 3], NurbsError> {
        crate::rational::nurbs_surface_point(self, u, v)
    }

    pub fn partials(&self, u: f64, v: f64) -> Result<([f64; 3], [f64; 3]), NurbsError> {
        crate::rational::nurbs_surface_partials(self, u, v)
    }
}

fn validate_degree_and_points(degree: usize, num_points: usize) -> Result<(), NurbsError> {
    if degree == 0 {
        return Err(NurbsError::InvalidDegree);
    }
    if num_points < degree + 1 {
        return Err(NurbsError::InvalidControlPoints(
            "not enough control points for degree".into(),
        ));
    }
    Ok(())
}

fn validate_surface(
    degree_u: usize,
    degree_v: usize,
    num_u: usize,
    num_v: usize,
    num_points: usize,
) -> Result<(), NurbsError> {
    validate_degree_and_points(degree_u, num_u)?;
    validate_degree_and_points(degree_v, num_v)?;
    if num_points != num_u * num_v {
        return Err(NurbsError::InvalidControlPoints(format!(
            "surface point count got {num_points}, expected {}",
            num_u * num_v
        )));
    }
    Ok(())
}

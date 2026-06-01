use crate::basis::{basis_function_derivatives, basis_functions, find_span};
use crate::types::{BSplineCurve, NurbsError};
use crate::utils::{add_scaled3, add_scaled4};

/// Algorithm A3.1: compute a point on a B-spline curve.
pub fn curve_point(curve: &BSplineCurve, u: f64) -> Result<[f64; 3], NurbsError> {
    let span = find_span(curve.control_points.len(), curve.degree, u, &curve.knots)?;
    let n = basis_functions(span, u, curve.degree, &curve.knots);
    let mut c = [0.0; 3];
    for (j, basis) in n.iter().enumerate().take(curve.degree + 1) {
        let idx = span - curve.degree + j;
        add_scaled3(&mut c, *basis, curve.control_points[idx]);
    }
    Ok(c)
}

/// Algorithm A3.2: compute derivatives of a B-spline curve through order `d`.
pub fn curve_derivatives(
    curve: &BSplineCurve,
    u: f64,
    d: usize,
) -> Result<Vec<[f64; 3]>, NurbsError> {
    let du = d.min(curve.degree);
    let span = find_span(curve.control_points.len(), curve.degree, u, &curve.knots)?;
    let nders = basis_function_derivatives(span, u, curve.degree, du, &curve.knots);
    let mut ck = vec![[0.0; 3]; d + 1];
    for k in 0..=du {
        for j in 0..=curve.degree {
            let idx = span - curve.degree + j;
            add_scaled3(&mut ck[k], nders[k][j], curve.control_points[idx]);
        }
    }
    Ok(ck)
}

pub(crate) fn homogeneous_curve_derivatives(
    degree: usize,
    knots: &crate::types::KnotVector,
    control_points: &[[f64; 4]],
    u: f64,
    d: usize,
) -> Result<Vec<[f64; 4]>, NurbsError> {
    let du = d.min(degree);
    let span = find_span(control_points.len(), degree, u, knots)?;
    let nders = basis_function_derivatives(span, u, degree, du, knots);
    let mut ck = vec![[0.0; 4]; d + 1];
    for k in 0..=du {
        for j in 0..=degree {
            let idx = span - degree + j;
            add_scaled4(&mut ck[k], nders[k][j], control_points[idx]);
        }
    }
    Ok(ck)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cubic_bezier_midpoint_matches_julia_fixture() {
        let curve = BSplineCurve::new(
            3,
            vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
            vec![
                [0.0, 0.0, 0.0],
                [1.0, 1.0, 0.0],
                [2.0, 0.0, 0.0],
                [3.0, 1.0, 0.0],
            ],
        )
        .unwrap();
        let p = curve_point(&curve, 0.5).unwrap();
        assert!((p[0] - 1.5).abs() < 1e-12);
        assert!((p[1] - 0.5).abs() < 1e-12);
        assert!(p[2].abs() < 1e-12);
    }
}

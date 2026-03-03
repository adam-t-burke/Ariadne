use crate::curves::homogeneous_curve_derivatives;
use crate::surfaces::homogeneous_surface_derivatives;
use crate::types::{NurbsCurve, NurbsError, NurbsSurface};
use crate::utils::{binomial, scale3, sub3};

/// Algorithm A4.1: compute a point on a NURBS curve.
pub fn nurbs_curve_point(curve: &NurbsCurve, u: f64) -> Result<[f64; 3], NurbsError> {
    let cw =
        homogeneous_curve_derivatives(curve.degree, &curve.knots, &curve.control_points, u, 0)?;
    project_point(cw[0])
}

/// Algorithm A4.2: compute NURBS curve derivatives through order `d`.
pub fn nurbs_curve_derivatives(
    curve: &NurbsCurve,
    u: f64,
    d: usize,
) -> Result<Vec<[f64; 3]>, NurbsError> {
    let ckw =
        homogeneous_curve_derivatives(curve.degree, &curve.knots, &curve.control_points, u, d)?;
    let mut ck = vec![[0.0; 3]; d + 1];
    let w0 = ckw[0][3];
    if w0.abs() <= f64::EPSILON {
        return Err(NurbsError::ZeroWeight);
    }

    for k in 0..=d {
        let mut v = [ckw[k][0], ckw[k][1], ckw[k][2]];
        for i in 1..=k {
            let term = scale3(ck[k - i], binomial(k, i) * ckw[i][3]);
            v = sub3(v, term);
        }
        ck[k] = scale3(v, 1.0 / w0);
    }
    Ok(ck)
}

/// Algorithm A4.3: compute a point on a NURBS surface.
pub fn nurbs_surface_point(surface: &NurbsSurface, u: f64, v: f64) -> Result<[f64; 3], NurbsError> {
    let sw = homogeneous_surface_derivatives(
        surface.degree_u,
        surface.degree_v,
        &surface.knots_u,
        &surface.knots_v,
        surface.num_u,
        surface.num_v,
        &surface.control_points,
        u,
        v,
        0,
    )?;
    project_point(sw[0][0])
}

/// Algorithm A4.4: compute rational surface derivatives through total order `d`.
pub fn nurbs_surface_derivatives(
    surface: &NurbsSurface,
    u: f64,
    v: f64,
    d: usize,
) -> Result<Vec<Vec<[f64; 3]>>, NurbsError> {
    let sklw = homogeneous_surface_derivatives(
        surface.degree_u,
        surface.degree_v,
        &surface.knots_u,
        &surface.knots_v,
        surface.num_u,
        surface.num_v,
        &surface.control_points,
        u,
        v,
        d,
    )?;
    let mut skl = vec![vec![[0.0; 3]; d + 1]; d + 1];
    let w00 = sklw[0][0][3];
    if w00.abs() <= f64::EPSILON {
        return Err(NurbsError::ZeroWeight);
    }

    for k in 0..=d {
        for l in 0..=(d - k) {
            let mut v3 = [sklw[k][l][0], sklw[k][l][1], sklw[k][l][2]];

            for j in 1..=l {
                let term = scale3(skl[k][l - j], binomial(l, j) * sklw[0][j][3]);
                v3 = sub3(v3, term);
            }

            for i in 1..=k {
                let mut v2 = scale3(skl[k - i][l], sklw[i][0][3]);
                for j in 1..=l {
                    let term = scale3(skl[k - i][l - j], binomial(l, j) * sklw[i][j][3]);
                    v2[0] += term[0];
                    v2[1] += term[1];
                    v2[2] += term[2];
                }
                let term = scale3(v2, binomial(k, i));
                v3 = sub3(v3, term);
            }

            skl[k][l] = scale3(v3, 1.0 / w00);
        }
    }
    Ok(skl)
}

pub fn nurbs_surface_partials(
    surface: &NurbsSurface,
    u: f64,
    v: f64,
) -> Result<([f64; 3], [f64; 3]), NurbsError> {
    let ders = nurbs_surface_derivatives(surface, u, v, 1)?;
    Ok((ders[1][0], ders[0][1]))
}

fn project_point(pw: [f64; 4]) -> Result<[f64; 3], NurbsError> {
    if pw[3].abs() <= f64::EPSILON {
        return Err(NurbsError::ZeroWeight);
    }
    Ok([pw[0] / pw[3], pw[1] / pw[3], pw[2] / pw[3]])
}

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_close(a: [f64; 3], b: [f64; 3], tol: f64) {
        for i in 0..3 {
            assert!(
                (a[i] - b[i]).abs() <= tol,
                "component {i}: {} != {}",
                a[i],
                b[i]
            );
        }
    }

    #[test]
    fn unit_weight_curve_reduces_to_bspline_fixture() {
        let curve = NurbsCurve::from_cartesian(
            3,
            vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
            vec![
                [0.0, 0.0, 0.0],
                [1.0, 1.0, 0.0],
                [2.0, 0.0, 0.0],
                [3.0, 1.0, 0.0],
            ],
            vec![1.0; 4],
        )
        .unwrap();
        assert_close(
            nurbs_curve_point(&curve, 0.5).unwrap(),
            [1.5, 0.5, 0.0],
            1e-12,
        );
    }

    #[test]
    fn weighted_curve_pulls_toward_weighted_control_point() {
        let weighted = NurbsCurve::from_cartesian(
            2,
            vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            vec![[0.0, 0.0, 0.0], [1.0, 1.0, 0.0], [2.0, 0.0, 0.0]],
            vec![1.0, 2.0, 1.0],
        )
        .unwrap();
        let unit = NurbsCurve::from_cartesian(
            2,
            vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            vec![[0.0, 0.0, 0.0], [1.0, 1.0, 0.0], [2.0, 0.0, 0.0]],
            vec![1.0, 1.0, 1.0],
        )
        .unwrap();
        assert!(weighted.point(0.5).unwrap()[1] > unit.point(0.5).unwrap()[1]);
        assert_close(weighted.point(0.0).unwrap(), [0.0, 0.0, 0.0], 1e-12);
        assert_close(weighted.point(1.0).unwrap(), [2.0, 0.0, 0.0], 1e-12);
    }

    #[test]
    fn curve_derivative_matches_finite_difference() {
        let curve = NurbsCurve::from_cartesian(
            2,
            vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            vec![[0.0, 0.0, 0.0], [1.0, 1.0, 0.0], [2.0, 0.0, 0.0]],
            vec![1.0, 2.0, 1.0],
        )
        .unwrap();
        let u = 0.4;
        let h = 1e-7;
        let p0 = curve.point(u - h).unwrap();
        let p1 = curve.point(u + h).unwrap();
        let fd = [
            (p1[0] - p0[0]) / (2.0 * h),
            (p1[1] - p0[1]) / (2.0 * h),
            (p1[2] - p0[2]) / (2.0 * h),
        ];
        let d = curve.derivatives(u, 1).unwrap()[1];
        assert_close(fd, d, 1e-7);
    }

    #[test]
    fn bilinear_surface_fixture() {
        let surface = NurbsSurface::from_cartesian(
            1,
            1,
            vec![0.0, 0.0, 1.0, 1.0],
            vec![0.0, 0.0, 1.0, 1.0],
            2,
            2,
            vec![
                [0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [1.0, 0.0, 0.0],
                [1.0, 1.0, 0.0],
            ],
            vec![1.0; 4],
        )
        .unwrap();
        assert_close(surface.point(0.5, 0.5).unwrap(), [0.5, 0.5, 0.0], 1e-12);
        assert_close(surface.point(0.0, 0.0).unwrap(), [0.0, 0.0, 0.0], 1e-12);
        assert_close(surface.point(1.0, 1.0).unwrap(), [1.0, 1.0, 0.0], 1e-12);
    }

    #[test]
    fn surface_partials_match_finite_difference() {
        let surface = NurbsSurface::from_cartesian(
            1,
            1,
            vec![0.0, 0.0, 1.0, 1.0],
            vec![0.0, 0.0, 1.0, 1.0],
            2,
            2,
            vec![
                [0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [1.0, 0.0, 0.0],
                [1.0, 1.0, 1.0],
            ],
            vec![1.0; 4],
        )
        .unwrap();
        let u = 0.4;
        let v = 0.6;
        let h = 1e-7;
        let (du, dv) = surface.partials(u, v).unwrap();
        let pu0 = surface.point(u - h, v).unwrap();
        let pu1 = surface.point(u + h, v).unwrap();
        let pv0 = surface.point(u, v - h).unwrap();
        let pv1 = surface.point(u, v + h).unwrap();
        let fdu = [
            (pu1[0] - pu0[0]) / (2.0 * h),
            (pu1[1] - pu0[1]) / (2.0 * h),
            (pu1[2] - pu0[2]) / (2.0 * h),
        ];
        let fdv = [
            (pv1[0] - pv0[0]) / (2.0 * h),
            (pv1[1] - pv0[1]) / (2.0 * h),
            (pv1[2] - pv0[2]) / (2.0 * h),
        ];
        assert_close(fdu, du, 1e-7);
        assert_close(fdv, dv, 1e-7);
    }
}

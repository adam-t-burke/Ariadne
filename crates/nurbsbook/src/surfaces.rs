use crate::basis::{basis_function_derivatives, basis_functions, find_span};
use crate::types::{BSplineSurface, NurbsError};
use crate::utils::{add_scaled3, add_scaled4};

/// Algorithm A3.5: compute a point on a B-spline surface.
pub fn surface_point(surface: &BSplineSurface, u: f64, v: f64) -> Result<[f64; 3], NurbsError> {
    let uspan = find_span(surface.num_u, surface.degree_u, u, &surface.knots_u)?;
    let vspan = find_span(surface.num_v, surface.degree_v, v, &surface.knots_v)?;
    let nu = basis_functions(uspan, u, surface.degree_u, &surface.knots_u);
    let nv = basis_functions(vspan, v, surface.degree_v, &surface.knots_v);

    let mut s = [0.0; 3];
    for (l, nv_l) in nv.iter().enumerate().take(surface.degree_v + 1) {
        let vi = vspan - surface.degree_v + l;
        let mut temp = [0.0; 3];
        for (k, nu_k) in nu.iter().enumerate().take(surface.degree_u + 1) {
            let ui = uspan - surface.degree_u + k;
            add_scaled3(&mut temp, *nu_k, surface.cp(ui, vi));
        }
        add_scaled3(&mut s, *nv_l, temp);
    }
    Ok(s)
}

/// Algorithm A3.6: compute B-spline surface derivatives through total order `d`.
pub fn surface_derivatives(
    surface: &BSplineSurface,
    u: f64,
    v: f64,
    d: usize,
) -> Result<Vec<Vec<[f64; 3]>>, NurbsError> {
    let du = d.min(surface.degree_u);
    let dv = d.min(surface.degree_v);
    let uspan = find_span(surface.num_u, surface.degree_u, u, &surface.knots_u)?;
    let vspan = find_span(surface.num_v, surface.degree_v, v, &surface.knots_v)?;
    let nuders = basis_function_derivatives(uspan, u, surface.degree_u, du, &surface.knots_u);
    let nvders = basis_function_derivatives(vspan, v, surface.degree_v, dv, &surface.knots_v);

    let mut skl = vec![vec![[0.0; 3]; d + 1]; d + 1];
    for k in 0..=du {
        let dd = (d - k).min(dv);
        let mut temp = vec![[0.0; 3]; surface.degree_v + 1];
        for (s, temp_s) in temp.iter_mut().enumerate().take(surface.degree_v + 1) {
            let vi = vspan - surface.degree_v + s;
            for r in 0..=surface.degree_u {
                let ui = uspan - surface.degree_u + r;
                add_scaled3(temp_s, nuders[k][r], surface.cp(ui, vi));
            }
        }
        for l in 0..=dd {
            for (s, temp_s) in temp.iter().enumerate().take(surface.degree_v + 1) {
                add_scaled3(&mut skl[k][l], nvders[l][s], *temp_s);
            }
        }
    }
    Ok(skl)
}

pub(crate) fn homogeneous_surface_derivatives(
    degree_u: usize,
    degree_v: usize,
    knots_u: &crate::types::KnotVector,
    knots_v: &crate::types::KnotVector,
    num_u: usize,
    num_v: usize,
    control_points: &[[f64; 4]],
    u: f64,
    v: f64,
    d: usize,
) -> Result<Vec<Vec<[f64; 4]>>, NurbsError> {
    let du = d.min(degree_u);
    let dv = d.min(degree_v);
    let uspan = find_span(num_u, degree_u, u, knots_u)?;
    let vspan = find_span(num_v, degree_v, v, knots_v)?;
    let nuders = basis_function_derivatives(uspan, u, degree_u, du, knots_u);
    let nvders = basis_function_derivatives(vspan, v, degree_v, dv, knots_v);

    let mut skl = vec![vec![[0.0; 4]; d + 1]; d + 1];
    for k in 0..=du {
        let dd = (d - k).min(dv);
        let mut temp = vec![[0.0; 4]; degree_v + 1];
        for (s, temp_s) in temp.iter_mut().enumerate().take(degree_v + 1) {
            let vi = vspan - degree_v + s;
            for r in 0..=degree_u {
                let ui = uspan - degree_u + r;
                add_scaled4(temp_s, nuders[k][r], control_points[ui * num_v + vi]);
            }
        }
        for l in 0..=dd {
            for (s, temp_s) in temp.iter().enumerate().take(degree_v + 1) {
                add_scaled4(&mut skl[k][l], nvders[l][s], *temp_s);
            }
        }
    }
    Ok(skl)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bilinear_surface_midpoint() {
        let surface = BSplineSurface::new(
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
        )
        .unwrap();
        let p = surface_point(&surface, 0.5, 0.5).unwrap();
        assert!((p[0] - 0.5).abs() < 1e-12);
        assert!((p[1] - 0.5).abs() < 1e-12);
        assert!(p[2].abs() < 1e-12);
    }
}

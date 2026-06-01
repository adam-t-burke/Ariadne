//! Latent-variable maps for variable support constraints.
//!
//! We optimize unconstrained latent parameters and map them to admissible
//! anchor positions. This enforces support constraints implicitly (hard) while
//! keeping gradients smooth for L-BFGS.

use ndarray::Array2;

use crate::types::{AnchorInfo, Problem, TheseusError, VariableSupportKind};

const SIGMOID_EPS: f64 = 1e-9;
const SPHERE_NORM_EPS: f64 = 1e-6;
const RAIL_LENGTH_EPS: f64 = 1e-9;
const DOMAIN_LENGTH_EPS: f64 = 1e-9;
const ANCHOR_LAMBDA_EPS: f64 = 1e-12;

#[inline]
fn sigmoid(z: f64) -> f64 {
    1.0 / (1.0 + (-z).exp())
}

#[inline]
fn logit(p: f64) -> f64 {
    let pc = p.clamp(SIGMOID_EPS, 1.0 - SIGMOID_EPS);
    (pc / (1.0 - pc)).ln()
}

#[inline]
fn anchor_lambda(problem: &Problem) -> f64 {
    problem
        .solver
        .anchor_saturation_lambda
        .abs()
        .max(ANCHOR_LAMBDA_EPS)
}

#[inline]
fn scaled_sigmoid(latent: f64, scale: f64) -> f64 {
    sigmoid(latent / scale.max(ANCHOR_LAMBDA_EPS))
}

/// Total latent dimension for all variable supports.
///
/// Legacy fallback: if supports are not configured but variable indices exist,
/// treat anchor variables as direct XYZ parameters (3 per anchor).
pub fn latent_dim(problem: &Problem) -> usize {
    if problem.anchors.variable_supports.is_empty() {
        return problem.anchors.variable_indices.len() * 3;
    }
    problem
        .anchors
        .variable_supports
        .iter()
        .map(|s| s.kind.latent_dim())
        .sum()
}

fn rail_delta(start: &[f64; 3], end: &[f64; 3]) -> ([f64; 3], f64) {
    let delta = [end[0] - start[0], end[1] - start[1], end[2] - start[2]];
    let length = (delta[0] * delta[0] + delta[1] * delta[1] + delta[2] * delta[2]).sqrt();
    (delta, length)
}

fn domain_width(domain: &[f64; 2], label: &str) -> Result<f64, TheseusError> {
    let width = domain[1] - domain[0];
    if !width.is_finite() || width <= DOMAIN_LENGTH_EPS {
        return Err(TheseusError::Shape(format!(
            "{label} domain must be increasing and finite"
        )));
    }
    Ok(width)
}

fn latent_to_domain(
    latent: f64,
    domain: &[f64; 2],
    lambda: f64,
) -> Result<(f64, f64, f64), TheseusError> {
    let width = domain_width(domain, "NURBS support")?;
    let scale = lambda * width.max(DOMAIN_LENGTH_EPS);
    let a = scaled_sigmoid(latent, scale);
    Ok((domain[0] + width * a, a, scale))
}

/// Build a deterministic initial latent vector from support definitions.
pub fn initial_latents(
    anchors: &AnchorInfo,
    anchor_saturation_lambda: f64,
) -> Result<Vec<f64>, TheseusError> {
    let lambda = anchor_saturation_lambda.abs().max(ANCHOR_LAMBDA_EPS);
    let nvar = anchors.variable_indices.len();
    if nvar == 0 {
        return Ok(Vec::new());
    }
    if anchors.variable_supports.is_empty() {
        let mut lat = Vec::with_capacity(nvar * 3);
        for i in 0..nvar {
            lat.push(anchors.initial_variable_positions[[i, 0]]);
            lat.push(anchors.initial_variable_positions[[i, 1]]);
            lat.push(anchors.initial_variable_positions[[i, 2]]);
        }
        return Ok(lat);
    }
    if anchors.variable_supports.len() != nvar {
        return Err(TheseusError::Shape(format!(
            "variable_supports length ({}) must match variable_indices length ({nvar})",
            anchors.variable_supports.len()
        )));
    }

    let mut latents = Vec::with_capacity(
        anchors
            .variable_supports
            .iter()
            .map(|s| s.kind.latent_dim())
            .sum(),
    );
    for support in &anchors.variable_supports {
        match &support.kind {
            VariableSupportKind::Sphere { .. } => {
                // Zero latent maps to zero offset (at reference position).
                latents.extend_from_slice(&[0.0, 0.0, 0.0]);
            }
            VariableSupportKind::Roller {
                enabled,
                lower,
                upper,
            } => {
                for d in 0..3 {
                    if !enabled[d] {
                        continue;
                    }
                    let w = upper[d] - lower[d];
                    if w <= 0.0 {
                        return Err(TheseusError::Shape(format!(
                            "roller axis {d} has invalid bounds: [{}, {}]",
                            lower[d], upper[d]
                        )));
                    }
                    let alpha = (0.0 - lower[d]) / w;
                    latents.push(logit(alpha) * lambda * w.max(ANCHOR_LAMBDA_EPS));
                }
            }
            VariableSupportKind::Rail { start, end } => {
                let p = support.reference_position;
                let (delta, length) = rail_delta(start, end);
                if length <= RAIL_LENGTH_EPS {
                    return Err(TheseusError::Shape(
                        "rail support has degenerate segment (start=end)".into(),
                    ));
                }
                let px = p[0] - start[0];
                let py = p[1] - start[1];
                let pz = p[2] - start[2];
                let t = (px * delta[0] + py * delta[1] + pz * delta[2]) / (length * length);
                latents.push(logit(t) * lambda * length);
            }
            VariableSupportKind::NurbsCurve {
                domain, initial_t, ..
            } => {
                let width = domain_width(domain, "NURBS curve")?;
                let alpha = (*initial_t - domain[0]) / width;
                latents.push(logit(alpha) * lambda * width);
            }
            VariableSupportKind::NurbsSurface {
                domain_u,
                domain_v,
                initial_uv,
                ..
            } => {
                let wu = domain_width(domain_u, "NURBS surface U")?;
                let wv = domain_width(domain_v, "NURBS surface V")?;
                latents.push(logit((initial_uv[0] - domain_u[0]) / wu) * lambda * wu);
                latents.push(logit((initial_uv[1] - domain_v[0]) / wv) * lambda * wv);
            }
        }
    }
    Ok(latents)
}

/// Map latent support parameters to world anchor positions.
pub fn map_latents_to_positions(
    problem: &Problem,
    latents: &[f64],
) -> Result<Array2<f64>, TheseusError> {
    let nvar = problem.anchors.variable_indices.len();
    if nvar == 0 {
        return Ok(Array2::zeros((0, 3)));
    }

    if problem.anchors.variable_supports.is_empty() {
        if latents.len() != nvar * 3 {
            return Err(TheseusError::Shape(format!(
                "legacy anchor latent length mismatch: got {}, expected {}",
                latents.len(),
                nvar * 3
            )));
        }
        let mut out = Array2::<f64>::zeros((nvar, 3));
        for i in 0..nvar {
            out[[i, 0]] = latents[i * 3];
            out[[i, 1]] = latents[i * 3 + 1];
            out[[i, 2]] = latents[i * 3 + 2];
        }
        return Ok(out);
    }

    if problem.anchors.variable_supports.len() != nvar {
        return Err(TheseusError::Shape(format!(
            "variable_supports length ({}) must match variable_indices length ({nvar})",
            problem.anchors.variable_supports.len()
        )));
    }

    let expected = latent_dim(problem);
    if latents.len() != expected {
        return Err(TheseusError::Shape(format!(
            "latent length mismatch: got {}, expected {expected}",
            latents.len()
        )));
    }

    let mut out = Array2::<f64>::zeros((nvar, 3));
    let lambda = anchor_lambda(problem);
    let mut p = 0usize;
    for (i, support) in problem.anchors.variable_supports.iter().enumerate() {
        let x0 = support.reference_position;
        match &support.kind {
            VariableSupportKind::Sphere { radius } => {
                let ux = latents[p];
                let uy = latents[p + 1];
                let uz = latents[p + 2];
                p += 3;
                let n = (ux * ux + uy * uy + uz * uz + SPHERE_NORM_EPS * SPHERE_NORM_EPS).sqrt();
                let scale = lambda * radius.abs().max(SPHERE_NORM_EPS);
                let s = radius * (n / scale).tanh() / n;
                out[[i, 0]] = x0[0] + s * ux;
                out[[i, 1]] = x0[1] + s * uy;
                out[[i, 2]] = x0[2] + s * uz;
            }
            VariableSupportKind::Roller {
                enabled,
                lower,
                upper,
            } => {
                let mut d = [0.0f64; 3];
                for axis in 0..3 {
                    if enabled[axis] {
                        let span = upper[axis] - lower[axis];
                        let scale = lambda * span.abs().max(ANCHOR_LAMBDA_EPS);
                        let a = scaled_sigmoid(latents[p], scale);
                        p += 1;
                        d[axis] = lower[axis] + span * a;
                    }
                }
                out[[i, 0]] = x0[0] + d[0];
                out[[i, 1]] = x0[1] + d[1];
                out[[i, 2]] = x0[2] + d[2];
            }
            VariableSupportKind::Rail { start, end } => {
                let (delta, length) = rail_delta(start, end);
                let scale = lambda * length.max(RAIL_LENGTH_EPS);
                let t = scaled_sigmoid(latents[p], scale);
                p += 1;
                out[[i, 0]] = start[0] + t * delta[0];
                out[[i, 1]] = start[1] + t * delta[1];
                out[[i, 2]] = start[2] + t * delta[2];
            }
            VariableSupportKind::NurbsCurve { curve, domain, .. } => {
                let (t, _, _) = latent_to_domain(latents[p], domain, lambda)?;
                p += 1;
                let point = curve.point(t).map_err(|e| {
                    TheseusError::Shape(format!("NURBS curve evaluation failed: {e}"))
                })?;
                out[[i, 0]] = point[0];
                out[[i, 1]] = point[1];
                out[[i, 2]] = point[2];
            }
            VariableSupportKind::NurbsSurface {
                surface,
                domain_u,
                domain_v,
                ..
            } => {
                let (u, _, _) = latent_to_domain(latents[p], domain_u, lambda)?;
                let (v, _, _) = latent_to_domain(latents[p + 1], domain_v, lambda)?;
                p += 2;
                let point = surface.point(u, v).map_err(|e| {
                    TheseusError::Shape(format!("NURBS surface evaluation failed: {e}"))
                })?;
                out[[i, 0]] = point[0];
                out[[i, 1]] = point[1];
                out[[i, 2]] = point[2];
            }
        }
    }
    Ok(out)
}

/// Chain-rule projection from dJ/d(anchor_xyz) to dJ/d(latent).
pub fn accumulate_latent_gradient(
    problem: &Problem,
    latents: &[f64],
    grad_nf: &Array2<f64>,
    out_latent_grad: &mut [f64],
) -> Result<(), TheseusError> {
    let nvar = problem.anchors.variable_indices.len();
    if nvar == 0 {
        return Ok(());
    }

    if problem.anchors.variable_supports.is_empty() {
        if out_latent_grad.len() != nvar * 3 {
            return Err(TheseusError::Shape(format!(
                "legacy latent gradient length mismatch: got {}, expected {}",
                out_latent_grad.len(),
                nvar * 3
            )));
        }
        for i in 0..nvar {
            let node = problem.anchors.variable_indices[i];
            out_latent_grad[i * 3] += grad_nf[[node, 0]];
            out_latent_grad[i * 3 + 1] += grad_nf[[node, 1]];
            out_latent_grad[i * 3 + 2] += grad_nf[[node, 2]];
        }
        return Ok(());
    }

    if out_latent_grad.len() != latent_dim(problem) || latents.len() != latent_dim(problem) {
        return Err(TheseusError::Shape(
            "latent gradient/parameter length mismatch".into(),
        ));
    }

    let lambda = anchor_lambda(problem);
    let mut p = 0usize;
    for (i, support) in problem.anchors.variable_supports.iter().enumerate() {
        let node = problem.anchors.variable_indices[i];
        let dl_dx = [grad_nf[[node, 0]], grad_nf[[node, 1]], grad_nf[[node, 2]]];
        match &support.kind {
            VariableSupportKind::Sphere { radius } => {
                let ux = latents[p];
                let uy = latents[p + 1];
                let uz = latents[p + 2];
                let u = [ux, uy, uz];
                let n = (ux * ux + uy * uy + uz * uz + SPHERE_NORM_EPS * SPHERE_NORM_EPS).sqrt();
                let scale = lambda * radius.abs().max(SPHERE_NORM_EPS);
                let rho = n / scale;
                let tanh_rho = rho.tanh();
                let sech2 = 1.0 - tanh_rho * tanh_rho;
                let s = radius * tanh_rho / n;
                let gprime = (sech2 * n / scale - tanh_rho) / (n * n);
                let dot = dl_dx[0] * ux + dl_dx[1] * uy + dl_dx[2] * uz;
                for j in 0..3 {
                    let ds_duj = radius * gprime * u[j] / n;
                    out_latent_grad[p + j] += s * dl_dx[j] + dot * ds_duj;
                }
                p += 3;
            }
            VariableSupportKind::Roller {
                enabled,
                lower,
                upper,
            } => {
                for axis in 0..3 {
                    if !enabled[axis] {
                        continue;
                    }
                    let span = upper[axis] - lower[axis];
                    let scale = lambda * span.abs().max(ANCHOR_LAMBDA_EPS);
                    let a = scaled_sigmoid(latents[p], scale);
                    let ddelta = span * a * (1.0 - a) / scale;
                    out_latent_grad[p] += dl_dx[axis] * ddelta;
                    p += 1;
                }
            }
            VariableSupportKind::Rail { start, end } => {
                let (dx_dt, length) = rail_delta(start, end);
                let latent_scale = lambda * length.max(RAIL_LENGTH_EPS);
                let a = scaled_sigmoid(latents[p], latent_scale);
                let dt = a * (1.0 - a);
                let dldt = dl_dx[0] * dx_dt[0] + dl_dx[1] * dx_dt[1] + dl_dx[2] * dx_dt[2];
                out_latent_grad[p] += dldt * dt / latent_scale;
                p += 1;
            }
            VariableSupportKind::NurbsCurve { curve, domain, .. } => {
                let (t, a, scale) = latent_to_domain(latents[p], domain, lambda)?;
                let tangent = curve.derivatives(t, 1).map_err(|e| {
                    TheseusError::Shape(format!("NURBS curve derivative failed: {e}"))
                })?[1];
                let width = domain_width(domain, "NURBS curve")?;
                let dldt = dl_dx[0] * tangent[0] + dl_dx[1] * tangent[1] + dl_dx[2] * tangent[2];
                out_latent_grad[p] += dldt * width * a * (1.0 - a) / scale;
                p += 1;
            }
            VariableSupportKind::NurbsSurface {
                surface,
                domain_u,
                domain_v,
                ..
            } => {
                let (u, au, scale_u) = latent_to_domain(latents[p], domain_u, lambda)?;
                let (v, av, scale_v) = latent_to_domain(latents[p + 1], domain_v, lambda)?;
                let (du, dv) = surface.partials(u, v).map_err(|e| {
                    TheseusError::Shape(format!("NURBS surface derivative failed: {e}"))
                })?;
                let width_u = domain_width(domain_u, "NURBS surface U")?;
                let width_v = domain_width(domain_v, "NURBS surface V")?;
                let dl_du = dl_dx[0] * du[0] + dl_dx[1] * du[1] + dl_dx[2] * du[2];
                let dl_dv = dl_dx[0] * dv[0] + dl_dx[1] * dv[1] + dl_dx[2] * dv[2];
                out_latent_grad[p] += dl_du * width_u * au * (1.0 - au) / scale_u;
                out_latent_grad[p + 1] += dl_dv * width_v * av * (1.0 - av) / scale_v;
                p += 2;
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{initial_latents, map_latents_to_positions};
    use crate::sparse::SparseColMatOwned;
    use crate::types::{
        AnchorInfo, Bounds, NetworkTopology, Problem, SolverOptions, VariableSupport,
        VariableSupportKind,
    };
    use ndarray::Array2;

    fn tiny_problem_with_support(support: VariableSupport) -> Problem {
        // Minimal topology: one edge, two nodes, one fixed (node 1), one free (node 0)
        let incidence =
            SparseColMatOwned::from_coo(1, 2, &[0usize, 0usize], &[0usize, 1usize], &[-1.0, 1.0])
                .unwrap();
        let free_node_indices = vec![0usize];
        let fixed_node_indices = vec![1usize];
        let free_incidence = incidence.extract_columns(&free_node_indices);
        let fixed_incidence = incidence.extract_columns(&fixed_node_indices);
        let topology = NetworkTopology {
            incidence,
            free_incidence,
            fixed_incidence,
            num_edges: 1,
            num_nodes: 2,
            free_node_indices,
            fixed_node_indices,
        };

        Problem {
            topology,
            free_node_loads: Array2::zeros((1, 3)),
            fixed_node_positions: Array2::from_shape_vec(
                (1, 3),
                vec![
                    support.reference_position[0],
                    support.reference_position[1],
                    support.reference_position[2],
                ],
            )
            .unwrap(),
            anchors: AnchorInfo {
                variable_indices: vec![support.node_index],
                fixed_indices: vec![support.node_index],
                reference_positions: Array2::from_shape_vec(
                    (1, 3),
                    vec![
                        support.reference_position[0],
                        support.reference_position[1],
                        support.reference_position[2],
                    ],
                )
                .unwrap(),
                initial_variable_positions: Array2::from_shape_vec(
                    (1, 3),
                    vec![
                        support.reference_position[0],
                        support.reference_position[1],
                        support.reference_position[2],
                    ],
                )
                .unwrap(),
                variable_supports: vec![support],
            },
            objectives: Vec::new(),
            bounds: Bounds::default_for(1),
            solver: SolverOptions::default(),
            self_weight: None,
            pressure: None,
        }
    }

    #[test]
    fn sphere_map_stays_within_radius() {
        let p = tiny_problem_with_support(VariableSupport {
            node_index: 1,
            reference_position: [1.0, 2.0, 3.0],
            kind: VariableSupportKind::Sphere { radius: 0.5 },
        });
        let x = map_latents_to_positions(&p, &[10.0, -2.0, 5.0]).unwrap();
        let dx = x[[0, 0]] - 1.0;
        let dy = x[[0, 1]] - 2.0;
        let dz = x[[0, 2]] - 3.0;
        let dist = (dx * dx + dy * dy + dz * dz).sqrt();
        assert!(dist < 0.5 + 1e-12);
    }

    #[test]
    fn roller_map_respects_axis_bounds() {
        let p = tiny_problem_with_support(VariableSupport {
            node_index: 1,
            reference_position: [0.0, 0.0, 0.0],
            kind: VariableSupportKind::Roller {
                enabled: [true, false, true],
                lower: [-1.0, 0.0, -0.25],
                upper: [2.0, 0.0, 0.75],
            },
        });
        let x = map_latents_to_positions(&p, &[9.0, -9.0]).unwrap();
        assert!(x[[0, 0]] > -1.0 && x[[0, 0]] < 2.0);
        assert!(x[[0, 1]].abs() < 1e-12); // locked axis
        assert!(x[[0, 2]] > -0.25 && x[[0, 2]] < 0.75);
    }

    #[test]
    fn rail_map_stays_on_segment() {
        let a = [0.0, 0.0, 0.0];
        let b = [4.0, 0.0, 0.0];
        let p = tiny_problem_with_support(VariableSupport {
            node_index: 1,
            reference_position: [1.0, 0.0, 0.0],
            kind: VariableSupportKind::Rail { start: a, end: b },
        });
        let x = map_latents_to_positions(&p, &[1.0]).unwrap();
        assert!(x[[0, 0]] > 0.0 && x[[0, 0]] < 4.0);
        assert!(x[[0, 1]].abs() < 1e-12);
        assert!(x[[0, 2]].abs() < 1e-12);
    }

    #[test]
    fn nurbs_curve_map_stays_on_curve() {
        let curve = nurbsbook::NurbsCurve::from_cartesian(
            2,
            vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            vec![[0.0, 0.0, 0.0], [1.0, 1.0, 0.0], [2.0, 0.0, 0.0]],
            vec![1.0, 1.0, 1.0],
        )
        .unwrap();
        let p = tiny_problem_with_support(VariableSupport {
            node_index: 1,
            reference_position: [1.0, 0.5, 0.0],
            kind: VariableSupportKind::NurbsCurve {
                curve,
                domain: [0.0, 1.0],
                initial_t: 0.5,
            },
        });
        let x = map_latents_to_positions(&p, &[0.0]).unwrap();
        assert!((x[[0, 0]] - 1.0).abs() < 1e-12);
        assert!((x[[0, 1]] - 0.5).abs() < 1e-12);
        assert!(x[[0, 2]].abs() < 1e-12);
    }

    #[test]
    fn nurbs_surface_map_stays_on_surface() {
        let surface = nurbsbook::NurbsSurface::from_cartesian(
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
        let p = tiny_problem_with_support(VariableSupport {
            node_index: 1,
            reference_position: [0.5, 0.5, 0.0],
            kind: VariableSupportKind::NurbsSurface {
                surface,
                domain_u: [0.0, 1.0],
                domain_v: [0.0, 1.0],
                initial_uv: [0.5, 0.5],
            },
        });
        let x = map_latents_to_positions(&p, &[0.0, 0.0]).unwrap();
        assert!((x[[0, 0]] - 0.5).abs() < 1e-12);
        assert!((x[[0, 1]] - 0.5).abs() < 1e-12);
        assert!(x[[0, 2]].abs() < 1e-12);
    }

    #[test]
    fn anchor_lambda_initial_latents_round_trip_roller() {
        let mut p = tiny_problem_with_support(VariableSupport {
            node_index: 1,
            reference_position: [0.0, 0.0, 0.0],
            kind: VariableSupportKind::Roller {
                enabled: [true, false, true],
                lower: [-1.0, 0.0, -2.0],
                upper: [3.0, 0.0, 2.0],
            },
        });
        p.solver.anchor_saturation_lambda = 2.5;
        let latents = initial_latents(&p.anchors, p.solver.anchor_saturation_lambda).unwrap();
        let x = map_latents_to_positions(&p, &latents).unwrap();

        assert!(x[[0, 0]].abs() < 1e-12);
        assert!(x[[0, 1]].abs() < 1e-12);
        assert!(x[[0, 2]].abs() < 1e-12);
    }

    #[test]
    fn anchor_lambda_initial_latents_round_trip_rail() {
        let mut p = tiny_problem_with_support(VariableSupport {
            node_index: 1,
            reference_position: [1.0, 0.0, 0.0],
            kind: VariableSupportKind::Rail {
                start: [0.0, 0.0, 0.0],
                end: [4.0, 0.0, 0.0],
            },
        });
        p.solver.anchor_saturation_lambda = 0.5;
        let latents = initial_latents(&p.anchors, p.solver.anchor_saturation_lambda).unwrap();
        let x = map_latents_to_positions(&p, &latents).unwrap();

        assert!((x[[0, 0]] - 1.0).abs() < 1e-12);
        assert!(x[[0, 1]].abs() < 1e-12);
        assert!(x[[0, 2]].abs() < 1e-12);
    }

    #[test]
    fn anchor_lambda_initial_latents_round_trip_nurbs_curve() {
        let curve = nurbsbook::NurbsCurve::from_cartesian(
            2,
            vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            vec![[0.0, 0.0, 0.0], [1.0, 1.0, 0.0], [2.0, 0.0, 0.0]],
            vec![1.0, 1.0, 1.0],
        )
        .unwrap();
        let mut p = tiny_problem_with_support(VariableSupport {
            node_index: 1,
            reference_position: [1.0, 0.5, 0.0],
            kind: VariableSupportKind::NurbsCurve {
                curve,
                domain: [0.0, 1.0],
                initial_t: 0.5,
            },
        });
        p.solver.anchor_saturation_lambda = 3.0;
        let latents = initial_latents(&p.anchors, p.solver.anchor_saturation_lambda).unwrap();
        let x = map_latents_to_positions(&p, &latents).unwrap();

        assert!((x[[0, 0]] - 1.0).abs() < 1e-12);
        assert!((x[[0, 1]] - 0.5).abs() < 1e-12);
        assert!(x[[0, 2]].abs() < 1e-12);
    }

    #[test]
    fn anchor_lambda_initial_latents_round_trip_nurbs_surface() {
        let surface = nurbsbook::NurbsSurface::from_cartesian(
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
        let mut p = tiny_problem_with_support(VariableSupport {
            node_index: 1,
            reference_position: [0.25, 0.75, 0.0],
            kind: VariableSupportKind::NurbsSurface {
                surface,
                domain_u: [0.0, 1.0],
                domain_v: [0.0, 1.0],
                initial_uv: [0.25, 0.75],
            },
        });
        p.solver.anchor_saturation_lambda = 0.25;
        let latents = initial_latents(&p.anchors, p.solver.anchor_saturation_lambda).unwrap();
        let x = map_latents_to_positions(&p, &latents).unwrap();

        assert!((x[[0, 0]] - 0.25).abs() < 1e-12);
        assert!((x[[0, 1]] - 0.75).abs() < 1e-12);
        assert!(x[[0, 2]].abs() < 1e-12);
    }

    #[test]
    fn anchor_lambda_sphere_map_stays_within_radius() {
        let mut p = tiny_problem_with_support(VariableSupport {
            node_index: 1,
            reference_position: [1.0, 2.0, 3.0],
            kind: VariableSupportKind::Sphere { radius: 0.5 },
        });
        p.solver.anchor_saturation_lambda = 4.0;
        let x = map_latents_to_positions(&p, &[10.0, -2.0, 5.0]).unwrap();
        let dx = x[[0, 0]] - 1.0;
        let dy = x[[0, 1]] - 2.0;
        let dz = x[[0, 2]] - 3.0;
        let dist = (dx * dx + dy * dy + dz * dz).sqrt();

        assert!(dist < 0.5 + 1e-12);
    }
}

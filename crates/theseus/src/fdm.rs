//! Forward FDM solver: assemble A(q), build RHS, factorise, triangular solve.

use crate::backend::cpu::{for_each_chunk_mut, for_each_chunk_mut2, run_sequential, CHUNK};
use crate::linear_solver::direct;
use crate::types::{FdmCache, PressureParams, Problem, SelfWeightParams, TheseusError};
use ndarray::Array2;
use std::sync::atomic::{AtomicBool, Ordering};

fn check_cancelled(cancel: Option<&AtomicBool>) -> Result<(), TheseusError> {
    if cancel.is_some_and(|flag| flag.load(Ordering::Acquire)) {
        Err(TheseusError::Cancelled)
    } else {
        Ok(())
    }
}

// ─────────────────────────────────────────────────────────────
//  Fixed-node position assembly
// ─────────────────────────────────────────────────────────────

/// Write current fixed-node positions into `cache.nf`, overlaying reference
/// positions and (optionally) variable anchor positions.
pub fn update_fixed_positions(
    cache: &mut FdmCache,
    problem: &Problem,
    anchor_positions: &Array2<f64>,
) {
    let fixed = &problem.topology.fixed_node_indices;
    let nn_fixed = fixed.len();

    // 1. Copy reference positions for fixed nodes
    let ref_pos = &problem.anchors.reference_positions;
    if ref_pos.nrows() == nn_fixed {
        for (i, &node) in fixed.iter().enumerate() {
            for d in 0..3 {
                cache.nf[[node, d]] = ref_pos[[i, d]];
            }
        }
    } else if ref_pos.nrows() == problem.topology.num_nodes {
        for &node in fixed {
            for d in 0..3 {
                cache.nf[[node, d]] = ref_pos[[node, d]];
            }
        }
    }

    // 2. Overlay variable anchors
    for (i, &node) in problem.anchors.variable_indices.iter().enumerate() {
        for d in 0..3 {
            cache.nf[[node, d]] = anchor_positions[[i, d]];
        }
    }
}

// ─────────────────────────────────────────────────────────────
//  System matrix assembly  A(q)
// ─────────────────────────────────────────────────────────────

/// Zero-allocation in-place update of A's values from current q.
/// A = Cn^T diag(q) Cn  via the precomputed `q_to_nz` mapping.
///
/// Each nonzero is an independent gather over its contributing edges, summed
/// in map order exactly as before; the nonzeros are processed in fixed
/// [`CHUNK`]-sized chunks (parallel above [`PAR_MIN_LEN`]), so the values
/// are bitwise identical for any thread count.
///
/// [`CHUNK`]: crate::backend::cpu::CHUNK
/// [`PAR_MIN_LEN`]: crate::backend::cpu::PAR_MIN_LEN
pub fn assemble_a(cache: &mut FdmCache) {
    let map = &cache.q_to_nz;
    let q = &cache.q;
    for_each_chunk_mut(&mut cache.a_matrix.values, CHUNK, |start, chunk| {
        for (j, value) in chunk.iter_mut().enumerate() {
            let nz = start + j;
            let range = map.nz_offsets[nz]..map.nz_offsets[nz + 1];
            *value = map.edge[range.clone()]
                .iter()
                .zip(&map.coeff[range])
                .map(|(&k, &c)| q[k as usize] * c)
                .sum();
        }
    });
}

// ─────────────────────────────────────────────────────────────
//  RHS assembly:  b = Pn − Cn^T diag(q) Cf Nf_fixed
// ─────────────────────────────────────────────────────────────

/// Build the right-hand side for A x = b.
///
///   b = Pn − Cnᵀ diag(q) Cf Nf_fixed
///
/// Only edges with exactly one fixed endpoint have a nonzero entry in
/// `diag(q) Cf Nf_fixed`, so the product is accumulated directly from the
/// precomputed `boundary_edges` list: each such edge adds `q_k · x_fixed` to
/// its free endpoint.
///
/// The free rows are processed in fixed [`CHUNK`]-row chunks (parallel above
/// [`PAR_MIN_LEN`] scalars); each chunk copies its rows of `pn` and gathers
/// the boundary edges of its rows, which `boundary_edges` stores sorted by
/// free row then edge. Per row the additions happen in ascending edge order,
/// the same order as the sequential edge loop, so the result is bitwise
/// identical for any thread count.
///
/// [`CHUNK`]: crate::backend::cpu::CHUNK
/// [`PAR_MIN_LEN`]: crate::backend::cpu::PAR_MIN_LEN
pub fn assemble_rhs(cache: &mut FdmCache, problem: &Problem) {
    let _ = problem;
    let nf = cache
        .nf
        .as_slice()
        .expect("node positions must be contiguous row-major");
    let pn = cache
        .pn
        .as_slice()
        .expect("loads must be contiguous row-major");
    let rhs = cache
        .rhs
        .as_slice_mut()
        .expect("rhs must be contiguous row-major");
    let boundary = &cache.boundary_edges;
    let q = &cache.q;
    for_each_chunk_mut(rhs, CHUNK * 3, |start, chunk| {
        chunk.copy_from_slice(&pn[start..start + chunk.len()]);
        let row0 = start / 3;
        let row1 = row0 + chunk.len() / 3;
        let lo = boundary.partition_point(|&(_, free, _)| free < row0);
        let hi = boundary.partition_point(|&(_, free, _)| free < row1);
        for &(k, free, fixed_node) in &boundary[lo..hi] {
            let qk = q[k];
            let row = (free - row0) * 3;
            for d in 0..3 {
                chunk[row + d] += qk * nf[fixed_node * 3 + d];
            }
        }
    });
}

// ─────────────────────────────────────────────────────────────
//  Linear solve  (dense fallback — will be replaced with LDL)
// ─────────────────────────────────────────────────────────────

/// Factor A via sparse Cholesky or LDL^T and solve A x = rhs for all 3 columns.
///
/// On first call, performs a fresh factorization (symbolic + numeric).
/// On subsequent calls, reuses the symbolic structure via `Factorization::update()`
/// — only numeric values change.
///
/// If the preferred strategy (Cholesky) fails because the matrix is no longer
/// SPD (e.g. q values drifted negative during optimisation), we automatically
/// fall back to LDL and rebuild the factorization from scratch. The factor
/// step is the direct path shared with
/// [`crate::linear_solver::DirectSolver`] (`linear_solver::direct::refactor`).
pub fn factor_and_solve(cache: &mut FdmCache, perturbation: f64) -> Result<(), TheseusError> {
    // Add diagonal perturbation if requested
    if perturbation > 0.0 {
        cache.a_matrix.add_diagonal(perturbation);
    }

    // Try the current factorization, falling back to LDL on Cholesky failure.
    direct::refactor(
        &mut cache.factorization,
        &mut cache.strategy,
        &cache.a_matrix,
        &mut cache.factor_stack,
    )?;

    let fac = cache
        .factorization
        .as_ref()
        .ok_or(TheseusError::MissingFactorization)?;
    fac.solve_into(
        &cache.rhs,
        &mut cache.x,
        &mut cache.solve_workspace,
        &mut cache.solve_stack,
    )?;
    let n = cache.a_matrix.nrows;
    for d in 0..3 {
        for i in 0..n {
            if !cache.x[[i, d]].is_finite() {
                return Err(TheseusError::Solver(direct::NON_FINITE_SOLUTION_MSG.into()));
            }
        }
    }

    Ok(())
}

// ─────────────────────────────────────────────────────────────
//  Top-level forward solve
// ─────────────────────────────────────────────────────────────

/// Full forward FDM solve.  Updates `cache.x`, `cache.nf`,
/// `cache.member_lengths`, `cache.member_forces`, `cache.reactions`.
pub fn solve_fdm(
    cache: &mut FdmCache,
    q: &[f64],
    problem: &Problem,
    anchor_positions: &Array2<f64>,
    perturbation: f64,
) -> Result<(), TheseusError> {
    solve_fdm_cancellable(cache, q, problem, anchor_positions, perturbation, None)
}

fn solve_fdm_cancellable(
    cache: &mut FdmCache,
    q: &[f64],
    problem: &Problem,
    anchor_positions: &Array2<f64>,
    perturbation: f64,
    cancel: Option<&AtomicBool>,
) -> Result<(), TheseusError> {
    check_cancelled(cancel)?;

    // 0. Sync q
    cache.q.copy_from_slice(q);

    // 1. Assemble A
    assemble_a(cache);

    // 2. Update fixed positions in Nf
    update_fixed_positions(cache, problem, anchor_positions);

    // 3. Assemble RHS
    assemble_rhs(cache, problem);

    // 4. Factor A and solve A x = rhs
    factor_and_solve(cache, perturbation)?;
    // Sparse factorization/triangular solve is one non-interruptible boundary.
    check_cancelled(cancel)?;

    // 5. Write free-node positions back to Nf
    for (i, &node) in problem.topology.free_node_indices.iter().enumerate() {
        for d in 0..3 {
            cache.nf[[node, d]] = cache.x[[i, d]];
        }
    }

    // 6. Compute derived geometry
    compute_geometry(cache, problem);

    // 7. Ensure geometry is finite (catches NaN/Inf from positions or downstream)
    for (i, &len) in cache.member_lengths.iter().enumerate() {
        if !len.is_finite() {
            return Err(TheseusError::Solver(format!(
                "FDM geometry produced non-finite member length at edge {i} (node positions may be NaN/Inf; \
                 check network topology and initial q for singularity or ill-conditioning).",
            )));
        }
    }
    for (i, &f) in cache.member_forces.iter().enumerate() {
        if !f.is_finite() {
            return Err(TheseusError::Solver(format!(
                "FDM geometry produced non-finite member force at edge {i} (check network and initial q).",
            )));
        }
    }

    check_cancelled(cancel)?;
    Ok(())
}

/// Compute member lengths, forces, and reactions from current positions and q.
/// Uses max(0, …) before sqrt to avoid NaN from floating-point negative squared length.
///
/// Parallel form (more than [`PAR_MIN_LEN`] elements and a pool of more than
/// one thread): two chunked passes, one over the edges for lengths and
/// forces, one over the nodes gathering the reactions
/// `R_u = Σ_{e ∋ u} q_e (x_other − x_u)` from the CSR adjacency. Incident
/// edges are listed in ascending edge order, so each node's sum is
/// accumulated in exactly the order of the sequential edge-scatter loop.
///
/// Sequential form (small problems, or a single-thread pool): one fused
/// pass over the edges scattering the reactions to both endpoints, which
/// touches each edge once instead of three times. Both forms are bitwise
/// identical (`tests/graph_loops_determinism.rs`).
///
/// [`PAR_MIN_LEN`]: crate::backend::cpu::PAR_MIN_LEN
pub fn compute_geometry(cache: &mut FdmCache, problem: &Problem) {
    let _ = problem;
    let ne = cache.member_lengths.len();
    let nf = cache
        .nf
        .as_slice()
        .expect("node positions must be contiguous row-major");
    let q = &cache.q;
    let starts = &cache.edge_starts;
    let ends = &cache.edge_ends;
    let reactions = cache
        .reactions
        .as_slice_mut()
        .expect("reactions must be contiguous row-major");

    if run_sequential(reactions.len()) {
        reactions.fill(0.0);
        for i in 0..ne {
            let s = starts[i] * 3;
            let e = ends[i] * 3;
            let qi = q[i];
            let dx = nf[e] - nf[s];
            let dy = nf[e + 1] - nf[s + 1];
            let dz = nf[e + 2] - nf[s + 2];
            let len_sq = dx * dx + dy * dy + dz * dz;
            let len = len_sq.max(0.0).sqrt();
            cache.member_lengths[i] = len;
            cache.member_forces[i] = qi * len;
            let (rx, ry, rz) = (dx * qi, dy * qi, dz * qi);
            reactions[s] += rx;
            reactions[s + 1] += ry;
            reactions[s + 2] += rz;
            reactions[e] -= rx;
            reactions[e + 1] -= ry;
            reactions[e + 2] -= rz;
        }
        return;
    }

    for_each_chunk_mut2(
        &mut cache.member_lengths,
        CHUNK,
        &mut cache.member_forces,
        CHUNK,
        |start, lengths, forces| {
            for (j, (len_out, force_out)) in lengths.iter_mut().zip(forces.iter_mut()).enumerate() {
                let i = start + j;
                let s = starts[i] * 3;
                let e = ends[i] * 3;
                let dx = nf[e] - nf[s];
                let dy = nf[e + 1] - nf[s + 1];
                let dz = nf[e + 2] - nf[s + 2];
                let len_sq = dx * dx + dy * dy + dz * dz;
                let len = len_sq.max(0.0).sqrt();
                *len_out = len;
                *force_out = q[i] * len;
            }
        },
    );

    // Reaction gather. The scatter form adds `q (x_e − x_s)` at the start
    // node and subtracts it at the end node; both are `q (x_other − x_u)`
    // at node `u`, and since IEEE negation is exact the branch-free form is
    // bitwise identical to the scatter.
    let adjacency = &cache.adjacency;
    for_each_chunk_mut(reactions, CHUNK * 3, |start, chunk| {
        for (j, row) in chunk.chunks_exact_mut(3).enumerate() {
            let u = start / 3 + j;
            let xu: &[f64; 3] = nf[u * 3..u * 3 + 3].try_into().expect("3 coordinates");
            let range = adjacency.range(u);
            let mut acc = [0.0f64; 3];
            for (&k, &other) in adjacency.edges[range.clone()]
                .iter()
                .zip(&adjacency.other[range])
            {
                let qi = q[k as usize];
                let o = other as usize * 3;
                let xo: &[f64; 3] = nf[o..o + 3].try_into().expect("3 coordinates");
                acc[0] += (xo[0] - xu[0]) * qi;
                acc[1] += (xo[1] - xu[1]) * qi;
                acc[2] += (xo[2] - xu[2]) * qi;
            }
            row.copy_from_slice(&acc);
        }
    });
}

// ─────────────────────────────────────────────────────────────
//  Self-weight load computation
// ─────────────────────────────────────────────────────────────

/// Update `cache.sw_mu` from current forces (sizing mode only).
fn update_sizing_mu(cache: &mut FdmCache, rho: f64, sigma: f64) {
    let ne = cache.member_lengths.len();
    for k in 0..ne {
        let a_k = cache.member_forces[k].abs() / sigma;
        cache.sw_mu[k] = rho * a_k;
        cache.cross_section_areas[k] = a_k;
    }
}

/// Add lumped self-weight loads into `cache.pn` (which should already
/// contain base loads).  Does not reset `pn` -- caller handles that.
fn accumulate_self_weight_loads(cache: &mut FdmCache, gravity: &[f64; 3]) {
    let ne = cache.member_lengths.len();
    for k in 0..ne {
        let mu_k = cache.sw_mu[k];
        let l_k = cache.member_lengths[k];
        let w_half = 0.5 * mu_k * l_k;

        let s = cache.edge_starts[k];
        let e = cache.edge_ends[k];

        if let Some(sf) = cache.node_to_free_idx[s] {
            for d in 0..3 {
                cache.pn[[sf, d]] += w_half * gravity[d];
            }
        }
        if let Some(ef) = cache.node_to_free_idx[e] {
            for d in 0..3 {
                cache.pn[[ef, d]] += w_half * gravity[d];
            }
        }
    }
}

// ─────────────────────────────────────────────────────────────
//  Pressure load computation
// ─────────────────────────────────────────────────────────────

/// Compute Newell area-weighted normal for a face (un-normalized).
/// Returns `(nx, ny, nz)` where `||(nx,ny,nz)|| = face area`.
pub fn newell_normal(face: &[usize], nf: &Array2<f64>) -> [f64; 3] {
    let nv = face.len();
    let mut n = [0.0f64; 3];
    for i in 0..nv {
        let j = (i + 1) % nv;
        let vi = face[i];
        let vj = face[j];
        // Cross product contribution: (v_i × v_j)
        n[0] += (nf[[vi, 1]] - nf[[vj, 1]]) * (nf[[vi, 2]] + nf[[vj, 2]]);
        n[1] += (nf[[vi, 2]] - nf[[vj, 2]]) * (nf[[vi, 0]] + nf[[vj, 0]]);
        n[2] += (nf[[vi, 0]] - nf[[vj, 0]]) * (nf[[vi, 1]] + nf[[vj, 1]]);
    }
    // Newell gives 2× the area-weighted normal
    n[0] *= 0.5;
    n[1] *= 0.5;
    n[2] *= 0.5;
    n
}

/// Add pressure loads into `cache.pn` (which should already contain
/// base loads and any self-weight).  Does not reset `pn`.
/// Dispatches on the pressure mode (Normal / Hydrostatic / Directional).
fn accumulate_pressure_loads(cache: &mut FdmCache, pressure: &PressureParams) {
    let faces = &pressure.face_topology().faces;
    match pressure {
        PressureParams::Normal { pressures, .. } => {
            for (f_idx, face) in faces.iter().enumerate() {
                let p_f = pressures[f_idx];
                let n = newell_normal(face, &cache.nf);
                let nv = face.len() as f64;
                for &vi in face {
                    if let Some(fi) = cache.node_to_free_idx[vi] {
                        for d in 0..3 {
                            cache.pn[[fi, d]] += p_f * n[d] / nv;
                        }
                    }
                }
            }
        }
        PressureParams::Hydrostatic {
            rho_fluid,
            g_magnitude,
            z_datum,
            up_direction,
            ..
        } => {
            for face in faces {
                let nv = face.len() as f64;
                if nv < 3.0 {
                    continue;
                }

                // Face centroid projected onto the "up" direction
                let mut centroid_up = 0.0;
                for &vi in face {
                    for d in 0..3 {
                        centroid_up += cache.nf[[vi, d]] * up_direction[d];
                    }
                }
                centroid_up /= nv;

                let depth = z_datum - centroid_up;
                if depth <= 0.0 {
                    continue;
                }

                let p_f = rho_fluid * g_magnitude * depth;
                let n = newell_normal(face, &cache.nf);
                for &vi in face {
                    if let Some(fi) = cache.node_to_free_idx[vi] {
                        for d in 0..3 {
                            cache.pn[[fi, d]] += p_f * n[d] / nv;
                        }
                    }
                }
            }
        }
        PressureParams::Directional {
            pressures,
            direction,
            ..
        } => {
            for (f_idx, face) in faces.iter().enumerate() {
                let p_f = pressures[f_idx];
                let n = newell_normal(face, &cache.nf);
                let nv = face.len() as f64;

                // Projected area = n_f · d_hat
                let a_proj: f64 = (0..3).map(|d| n[d] * direction[d]).sum();
                if a_proj <= 0.0 {
                    continue;
                }

                let load_mag = p_f * a_proj / nv;
                for &vi in face {
                    if let Some(fi) = cache.node_to_free_idx[vi] {
                        for d in 0..3 {
                            cache.pn[[fi, d]] += load_mag * direction[d];
                        }
                    }
                }
            }
        }
    }
}

/// Apply the forward hydrostatic load Jacobian `(dp/dx) * direction`.
///
/// `direction` and `out` contain free-node vectors. Fixed vertices participate
/// in the current face geometry but have zero variation. `load_scale` is used
/// by continuation (`0 <= load_scale <= 1`).
pub fn apply_hydrostatic_load_jacobian(
    cache: &FdmCache,
    pressure: &PressureParams,
    direction: &Array2<f64>,
    load_scale: f64,
    out: &mut Array2<f64>,
) {
    out.fill(0.0);
    let PressureParams::Hydrostatic {
        face_topology,
        rho_fluid,
        g_magnitude,
        z_datum,
        up_direction,
        ..
    } = pressure
    else {
        return;
    };

    let rho_g = load_scale * rho_fluid * g_magnitude;
    for face in &face_topology.faces {
        let nv = face.len();
        if nv < 3 {
            continue;
        }
        let nv_f = nv as f64;

        let mut centroid_up = 0.0;
        for &vi in face {
            for d in 0..3 {
                centroid_up += cache.nf[[vi, d]] * up_direction[d];
            }
        }
        centroid_up /= nv_f;
        let depth = z_datum - centroid_up;
        if depth <= 0.0 {
            continue;
        }

        let normal = newell_normal(face, &cache.nf);
        let pressure_value = rho_g * depth;

        let mut delta_centroid_up = 0.0;
        let mut delta_normal = [0.0; 3];
        for j in 0..nv {
            let vj = face[j];
            let Some(jf) = cache.node_to_free_idx[vj] else {
                continue;
            };
            for d in 0..3 {
                delta_centroid_up += direction[[jf, d]] * up_direction[d] / nv_f;
            }

            let prev = face[(j + nv - 1) % nv];
            let next = face[(j + 1) % nv];
            let a = [
                cache.nf[[prev, 0]] - cache.nf[[next, 0]],
                cache.nf[[prev, 1]] - cache.nf[[next, 1]],
                cache.nf[[prev, 2]] - cache.nf[[next, 2]],
            ];
            let delta = [direction[[jf, 0]], direction[[jf, 1]], direction[[jf, 2]]];
            delta_normal[0] += 0.5 * (a[1] * delta[2] - a[2] * delta[1]);
            delta_normal[1] += 0.5 * (a[2] * delta[0] - a[0] * delta[2]);
            delta_normal[2] += 0.5 * (a[0] * delta[1] - a[1] * delta[0]);
        }

        let delta_pressure = -rho_g * delta_centroid_up;
        for &vi in face {
            if let Some(fi) = cache.node_to_free_idx[vi] {
                for d in 0..3 {
                    out[[fi, d]] +=
                        (delta_pressure * normal[d] + pressure_value * delta_normal[d]) / nv_f;
                }
            }
        }
    }
}

fn accumulate_hydrostatic_loads_scaled(
    cache: &mut FdmCache,
    pressure: &PressureParams,
    load_scale: f64,
) {
    let PressureParams::Hydrostatic {
        face_topology,
        rho_fluid,
        g_magnitude,
        z_datum,
        up_direction,
        ..
    } = pressure
    else {
        return;
    };

    for face in &face_topology.faces {
        let nv = face.len();
        if nv < 3 {
            continue;
        }
        let nv_f = nv as f64;
        let mut centroid_up = 0.0;
        for &vi in face {
            for d in 0..3 {
                centroid_up += cache.nf[[vi, d]] * up_direction[d];
            }
        }
        centroid_up /= nv_f;
        let depth = z_datum - centroid_up;
        if depth <= 0.0 {
            continue;
        }

        let pressure_value = load_scale * rho_fluid * g_magnitude * depth;
        let normal = newell_normal(face, &cache.nf);
        for &vi in face {
            if let Some(fi) = cache.node_to_free_idx[vi] {
                for d in 0..3 {
                    cache.pn[[fi, d]] += pressure_value * normal[d] / nv_f;
                }
            }
        }
    }
}

// ─────────────────────────────────────────────────────────────
//  Universal forward solve with geometry-dependent loads
// ─────────────────────────────────────────────────────────────

fn vector_norm(values: &[f64]) -> f64 {
    values.iter().map(|v| v * v).sum::<f64>().sqrt()
}

fn flatten_xyz(values: &Array2<f64>) -> Vec<f64> {
    let mut flat = vec![0.0; values.nrows() * 3];
    for i in 0..values.nrows() {
        for d in 0..3 {
            flat[i * 3 + d] = values[[i, d]];
        }
    }
    flat
}

fn unflatten_xyz(values: &[f64], out: &mut Array2<f64>) {
    for i in 0..out.nrows() {
        for d in 0..3 {
            out[[i, d]] = values[i * 3 + d];
        }
    }
}

fn sync_free_positions(cache: &mut FdmCache, problem: &Problem) {
    for (i, &node) in problem.topology.free_node_indices.iter().enumerate() {
        for d in 0..3 {
            cache.nf[[node, d]] = cache.x[[i, d]];
        }
    }
}

fn apply_a_xyz(cache: &FdmCache, input: &Array2<f64>, out: &mut Array2<f64>) {
    out.fill(0.0);
    for col in 0..cache.a_matrix.ncols {
        let start = cache.a_matrix.col_ptrs[col] as usize;
        let end = cache.a_matrix.col_ptrs[col + 1] as usize;
        for nz in start..end {
            let row = cache.a_matrix.row_indices[nz] as usize;
            let value = cache.a_matrix.values[nz];
            for d in 0..3 {
                out[[row, d]] += value * input[[col, d]];
            }
        }
    }
}

fn apply_a_preconditioner(cache: &mut FdmCache, input: &[f64]) -> Result<Vec<f64>, TheseusError> {
    let n = cache.x.nrows();
    let rhs = Array2::from_shape_fn((n, 3), |(i, d)| input[i * 3 + d]);
    let mut solution = Array2::zeros((n, 3));
    let factorization = cache
        .factorization
        .as_ref()
        .ok_or(TheseusError::MissingFactorization)?;
    factorization.solve_into(
        &rhs,
        &mut solution,
        &mut cache.solve_workspace,
        &mut cache.solve_stack,
    )?;
    Ok(flatten_xyz(&solution))
}

fn evaluate_hydrostatic_residual(
    cache: &mut FdmCache,
    problem: &Problem,
    pressure: &PressureParams,
    load_scale: f64,
) -> Result<(Vec<f64>, f64), TheseusError> {
    sync_free_positions(cache, problem);
    cache.pn.assign(&cache.pn_base);
    accumulate_hydrostatic_loads_scaled(cache, pressure, load_scale);
    if cache.pn.iter().any(|v| !v.is_finite()) || cache.x.iter().any(|v| !v.is_finite()) {
        return Err(TheseusError::Solver(
            "Hydrostatic follower solve diverged: non-finite load or position. \
             Increase q, reduce the load, or use smaller continuation steps."
                .into(),
        ));
    }

    assemble_rhs(cache, problem);
    let mut ax = Array2::zeros(cache.x.raw_dim());
    apply_a_xyz(cache, &cache.x, &mut ax);
    let mut residual = vec![0.0; cache.x.len()];
    let mut load_norm_sq = 0.0;
    for i in 0..cache.x.nrows() {
        for d in 0..3 {
            residual[i * 3 + d] = cache.rhs[[i, d]] - ax[[i, d]];
            load_norm_sq += cache.rhs[[i, d]] * cache.rhs[[i, d]];
        }
    }
    Ok((residual, load_norm_sq.sqrt()))
}

fn apply_hydrostatic_newton_operator(
    cache: &mut FdmCache,
    pressure: &PressureParams,
    load_scale: f64,
    input: &[f64],
) -> Result<Vec<f64>, TheseusError> {
    let n = cache.x.nrows();
    let direction = Array2::from_shape_fn((n, 3), |(i, d)| input[i * 3 + d]);
    let mut a_direction = Array2::zeros((n, 3));
    let mut load_direction = Array2::zeros((n, 3));
    apply_a_xyz(cache, &direction, &mut a_direction);
    apply_hydrostatic_load_jacobian(cache, pressure, &direction, load_scale, &mut load_direction);
    for i in 0..n {
        for d in 0..3 {
            a_direction[[i, d]] -= load_direction[[i, d]];
        }
    }
    apply_a_preconditioner(cache, &flatten_xyz(&a_direction))
}

/// Restarted left-preconditioned GMRES. The supplied operator is already
/// preconditioned, so this solves `operator(x) = rhs`.
fn gmres<F>(
    dimension: usize,
    rhs: &[f64],
    max_iterations: usize,
    tolerance: f64,
    cancel: Option<&AtomicBool>,
    mut operator: F,
) -> Result<Vec<f64>, TheseusError>
where
    F: FnMut(&[f64]) -> Result<Vec<f64>, TheseusError>,
{
    check_cancelled(cancel)?;
    if dimension == 0 {
        return Ok(Vec::new());
    }
    let restart = dimension.clamp(1, 30);
    let rhs_norm = vector_norm(rhs).max(1e-30);
    let mut x = vec![0.0; dimension];
    let mut iterations = 0;

    while iterations < max_iterations {
        check_cancelled(cancel)?;
        let ax = operator(&x)?;
        let residual: Vec<f64> = rhs.iter().zip(ax).map(|(b, a)| b - a).collect();
        let beta = vector_norm(&residual);
        if beta <= tolerance * rhs_norm {
            return Ok(x);
        }

        let cycle = restart.min(max_iterations - iterations);
        let mut basis = Vec::with_capacity(cycle + 1);
        basis.push(residual.iter().map(|v| v / beta).collect::<Vec<_>>());
        let mut h = vec![vec![0.0; cycle]; cycle + 1];
        let mut cosines = vec![0.0; cycle];
        let mut sines = vec![0.0; cycle];
        let mut g = vec![0.0; cycle + 1];
        g[0] = beta;
        let mut used = cycle;

        for j in 0..cycle {
            check_cancelled(cancel)?;
            let mut w = operator(&basis[j])?;
            for i in 0..=j {
                h[i][j] = basis[i].iter().zip(&w).map(|(a, b)| a * b).sum();
                for (wk, vik) in w.iter_mut().zip(&basis[i]) {
                    *wk -= h[i][j] * vik;
                }
            }
            h[j + 1][j] = vector_norm(&w);
            if h[j + 1][j] > 1e-14 {
                basis.push(w.iter().map(|v| v / h[j + 1][j]).collect());
            } else {
                basis.push(vec![0.0; dimension]);
            }

            for i in 0..j {
                let top = cosines[i] * h[i][j] + sines[i] * h[i + 1][j];
                let bottom = -sines[i] * h[i][j] + cosines[i] * h[i + 1][j];
                h[i][j] = top;
                h[i + 1][j] = bottom;
            }

            let hyp = h[j][j].hypot(h[j + 1][j]);
            if hyp > 0.0 {
                cosines[j] = h[j][j] / hyp;
                sines[j] = h[j + 1][j] / hyp;
            } else {
                cosines[j] = 1.0;
            }
            h[j][j] = cosines[j] * h[j][j] + sines[j] * h[j + 1][j];
            h[j + 1][j] = 0.0;
            g[j + 1] = -sines[j] * g[j];
            g[j] *= cosines[j];

            iterations += 1;
            if g[j + 1].abs() <= tolerance * rhs_norm || iterations >= max_iterations {
                used = j + 1;
                break;
            }
        }

        let mut y = vec![0.0; used];
        for i in (0..used).rev() {
            let tail: f64 = ((i + 1)..used).map(|j| h[i][j] * y[j]).sum();
            if h[i][i].abs() <= 1e-14 {
                return Err(TheseusError::Solver(
                    "Hydrostatic Newton GMRES encountered a singular Krylov system.".into(),
                ));
            }
            y[i] = (g[i] - tail) / h[i][i];
        }
        for j in 0..used {
            for (xk, vjk) in x.iter_mut().zip(&basis[j]) {
                *xk += y[j] * vjk;
            }
        }
    }

    let ax = operator(&x)?;
    let residual: Vec<f64> = rhs.iter().zip(ax).map(|(b, a)| b - a).collect();
    if vector_norm(&residual) <= tolerance * rhs_norm {
        Ok(x)
    } else {
        Err(TheseusError::Solver(format!(
            "Hydrostatic Newton GMRES did not converge within {max_iterations} iterations."
        )))
    }
}

fn solve_hydrostatic_newton(
    cache: &mut FdmCache,
    q: &[f64],
    problem: &Problem,
    anchor_positions: &Array2<f64>,
    pressure: &PressureParams,
    perturbation: f64,
    cancel: Option<&AtomicBool>,
) -> Result<(), TheseusError> {
    check_cancelled(cancel)?;
    let max_newton_iters = pressure.max_iters().max(1);
    let tolerance = pressure.tolerance().max(1e-12);

    // The zero-load FDM solution is the starting point and also assembles /
    // factorizes A(q), which is reused as the Newton preconditioner.
    cache.pn.assign(&cache.pn_base);
    solve_fdm_cancellable(cache, q, problem, anchor_positions, perturbation, cancel)?;

    let mut load_scale: f64 = 0.0;
    let mut step: f64 = 0.125;
    const MIN_STEP: f64 = 1.0 / 4096.0;

    while load_scale < 1.0 - 1e-12 {
        check_cancelled(cancel)?;
        let target_scale = (load_scale + step).min(1.0);
        let accepted_x = cache.x.clone();
        let mut stage_converged = false;
        let mut stage_error: Option<TheseusError> = None;

        for _ in 0..max_newton_iters {
            check_cancelled(cancel)?;
            let (residual, load_norm) =
                evaluate_hydrostatic_residual(cache, problem, pressure, target_scale)?;
            let residual_norm = vector_norm(&residual);
            if residual_norm <= tolerance * load_norm.max(1.0) {
                stage_converged = true;
                break;
            }

            let preconditioned_rhs = apply_a_preconditioner(cache, &residual)?;
            let dimension = residual.len();
            let gmres_result = gmres(
                dimension,
                &preconditioned_rhs,
                dimension.clamp(1, 60),
                (tolerance.sqrt()).clamp(1e-8, 1e-3),
                cancel,
                |v| apply_hydrostatic_newton_operator(cache, pressure, target_scale, v),
            );
            let delta = match gmres_result {
                Ok(delta) => delta,
                Err(error) => {
                    stage_error = Some(error);
                    break;
                }
            };

            let current_x = flatten_xyz(&cache.x);
            let mut alpha = 1.0;
            let mut line_search_accepted = false;
            while alpha >= 1.0 / 1024.0 {
                check_cancelled(cancel)?;
                let trial: Vec<f64> = current_x
                    .iter()
                    .zip(&delta)
                    .map(|(x, dx)| x + alpha * dx)
                    .collect();
                unflatten_xyz(&trial, &mut cache.x);
                match evaluate_hydrostatic_residual(cache, problem, pressure, target_scale) {
                    Ok((trial_residual, _))
                        if vector_norm(&trial_residual) <= (1.0 - 1e-4 * alpha) * residual_norm =>
                    {
                        line_search_accepted = true;
                        break;
                    }
                    _ => {
                        alpha *= 0.5;
                    }
                }
            }
            if !line_search_accepted {
                unflatten_xyz(&current_x, &mut cache.x);
                sync_free_positions(cache, problem);
                stage_error = Some(TheseusError::Solver(
                    "Hydrostatic Newton line search could not reduce the equilibrium residual."
                        .into(),
                ));
                break;
            }
        }

        if stage_converged {
            load_scale = target_scale;
            step = (step * 1.5).min(1.0 - load_scale);
            if step <= 0.0 {
                break;
            }
        } else {
            cache.x.assign(&accepted_x);
            sync_free_positions(cache, problem);
            step *= 0.5;
            if step < MIN_STEP {
                let detail = stage_error.map(|e| e.to_string()).unwrap_or_else(|| {
                    format!("Newton did not converge within {max_newton_iters} iterations.")
                });
                return Err(TheseusError::Solver(format!(
                    "Hydrostatic follower continuation failed near load factor \
                     {target_scale:.6}: {detail} Increase q or reduce the hydrostatic load."
                )));
            }
        }
    }

    check_cancelled(cancel)?;
    let (residual, load_norm) = evaluate_hydrostatic_residual(cache, problem, pressure, 1.0)?;
    let relative_residual = vector_norm(&residual) / load_norm.max(1.0);
    if !relative_residual.is_finite() || relative_residual > tolerance {
        return Err(TheseusError::Solver(format!(
            "Hydrostatic follower Newton solve ended with relative equilibrium residual \
             {relative_residual:.3e} (tolerance {tolerance:.3e})."
        )));
    }
    compute_geometry(cache, problem);
    Ok(())
}

/// Forward solve that handles self-weight and/or pressure loads via iteration.
///
/// When neither `problem.self_weight` nor `problem.pressure` is set, this
/// delegates directly to [`solve_fdm`] with zero overhead.
pub fn solve_fdm_with_loads(
    cache: &mut FdmCache,
    q: &[f64],
    problem: &Problem,
    anchor_positions: &Array2<f64>,
    perturbation: f64,
) -> Result<(), TheseusError> {
    solve_fdm_with_loads_cancellable(cache, q, problem, anchor_positions, perturbation, None)
}

/// Forward solve with cooperative cancellation between nonlinear iterations,
/// GMRES iterations, and immediately after each sparse factorization/solve.
///
/// A single sparse factorization or triangular solve cannot be interrupted.
pub fn solve_fdm_with_loads_cancellable(
    cache: &mut FdmCache,
    q: &[f64],
    problem: &Problem,
    anchor_positions: &Array2<f64>,
    perturbation: f64,
    cancel: Option<&AtomicBool>,
) -> Result<(), TheseusError> {
    check_cancelled(cancel)?;
    let has_sw = problem.self_weight.is_some();
    let has_pressure = problem.pressure.is_some();

    if !has_sw && !has_pressure {
        return solve_fdm_cancellable(cache, q, problem, anchor_positions, perturbation, cancel);
    }

    // Hydrostatic pressure couples all coordinate directions. Solve its true
    // follower equilibrium with Newton + load continuation rather than the
    // unstable lagged-load Picard map. Combined self-weight remains on the
    // legacy coupled Picard path until its Jacobian is available.
    if !has_sw {
        if let Some(pressure @ PressureParams::Hydrostatic { .. }) = &problem.pressure {
            return solve_hydrostatic_newton(
                cache,
                q,
                problem,
                anchor_positions,
                pressure,
                perturbation,
                cancel,
            );
        }
    }

    // Determine iteration parameters from whichever load type is active
    let max_iters = match (&problem.self_weight, &problem.pressure) {
        (Some(sw), Some(pr)) => sw.max_iters().max(pr.max_iters()),
        (Some(sw), None) => sw.max_iters(),
        (None, Some(pr)) => pr.max_iters(),
        (None, None) => unreachable!(),
    };
    let tolerance = match (&problem.self_weight, &problem.pressure) {
        (Some(sw), Some(pr)) => sw.tolerance().min(pr.tolerance()),
        (Some(sw), None) => sw.tolerance(),
        (None, Some(pr)) => pr.tolerance(),
        (None, None) => unreachable!(),
    };
    let relaxation = match &problem.self_weight {
        Some(sw) => sw.relaxation(),
        None => 1.0,
    };

    // Reset pn to base loads
    cache.pn.assign(&cache.pn_base);

    let nn_free = problem.topology.free_node_indices.len();
    let mut converged = false;

    for _iter in 0..max_iters {
        check_cancelled(cancel)?;
        // Inner linear FDM solve with current loads
        solve_fdm_cancellable(cache, q, problem, anchor_positions, perturbation, cancel)?;

        // Update sizing mu from current forces if applicable
        if let Some(SelfWeightParams::Sizing { rho, sigma, .. }) = &problem.self_weight {
            update_sizing_mu(cache, *rho, *sigma);
        }

        // Save current pn, then rebuild from base + geometry-dependent loads
        cache.pn_prev.assign(&cache.pn);
        cache.pn.assign(&cache.pn_base);

        if let Some(sw) = &problem.self_weight {
            accumulate_self_weight_loads(cache, sw.gravity());
        }
        if let Some(pr) = &problem.pressure {
            accumulate_pressure_loads(cache, pr);
        }

        // Convergence: relative change in load vector
        let mut delta_sq = 0.0;
        let mut pn_norm_sq = 0.0;
        for i in 0..nn_free {
            for d in 0..3 {
                let diff = cache.pn[[i, d]] - cache.pn_prev[[i, d]];
                delta_sq += diff * diff;
                pn_norm_sq += cache.pn[[i, d]] * cache.pn[[i, d]];
            }
        }
        if pn_norm_sq > 0.0 && delta_sq / pn_norm_sq < tolerance * tolerance {
            converged = true;
            break;
        }

        // Relaxation: blend new loads with old
        if relaxation < 1.0 {
            for i in 0..nn_free {
                for d in 0..3 {
                    cache.pn[[i, d]] =
                        relaxation * cache.pn[[i, d]] + (1.0 - relaxation) * cache.pn_prev[[i, d]];
                }
            }
        }
    }

    // Skip redundant final solve when the last iteration already converged.
    if !converged {
        solve_fdm_cancellable(cache, q, problem, anchor_positions, perturbation, cancel)?;
    }

    check_cancelled(cancel)?;
    Ok(())
}

#[cfg(test)]
mod hydrostatic_tests {
    use super::*;
    use crate::sparse::SparseColMatOwned;
    use crate::types::{AnchorInfo, Bounds, FaceTopology, NetworkTopology, SolverOptions};

    fn incidence(edges: &[(usize, usize)], num_nodes: usize) -> SparseColMatOwned {
        let mut rows = Vec::with_capacity(edges.len() * 2);
        let mut cols = Vec::with_capacity(edges.len() * 2);
        let mut values = Vec::with_capacity(edges.len() * 2);
        for (edge, &(start, end)) in edges.iter().enumerate() {
            rows.extend([edge, edge]);
            cols.extend([start, end]);
            values.extend([-1.0, 1.0]);
        }
        SparseColMatOwned::from_coo(edges.len(), num_nodes, &rows, &cols, &values).unwrap()
    }

    /// Square grid with every boundary node fixed. Faces are clockwise from
    /// +Z, so hydrostatic pressure acts toward -Z.
    fn anchored_grid_problem_size(size: usize, rho: f64, tolerance: f64) -> Problem {
        let mut edges = Vec::new();
        for row in 0..size {
            for col in 0..(size - 1) {
                edges.push((row * size + col, row * size + col + 1));
            }
        }
        for col in 0..size {
            for row in 0..(size - 1) {
                edges.push((row * size + col, (row + 1) * size + col));
            }
        }
        let num_nodes = size * size;
        let full = incidence(&edges, num_nodes);
        let mut free_indices = Vec::new();
        let mut fixed_indices = Vec::new();
        for row in 0..size {
            for col in 0..size {
                let node = row * size + col;
                if row == 0 || col == 0 || row + 1 == size || col + 1 == size {
                    fixed_indices.push(node);
                } else {
                    free_indices.push(node);
                }
            }
        }
        let free_incidence = full.extract_columns(&free_indices);
        let fixed_incidence = full.extract_columns(&fixed_indices);

        let fixed_positions = Array2::from_shape_fn((fixed_indices.len(), 3), |(i, d)| {
            let node = fixed_indices[i];
            match d {
                0 => (node % size) as f64,
                1 => (node / size) as f64,
                _ => 0.0,
            }
        });
        let mut faces = Vec::new();
        for row in 0..(size - 1) {
            for col in 0..(size - 1) {
                let lower_left = row * size + col;
                faces.push(vec![
                    lower_left,
                    lower_left + size,
                    lower_left + size + 1,
                    lower_left + 1,
                ]);
            }
        }
        let pressure = PressureParams::Hydrostatic {
            face_topology: FaceTopology { faces },
            rho_fluid: rho,
            g_magnitude: 1.0,
            z_datum: 1.0,
            up_direction: [0.0, 0.0, 1.0],
            max_iters: 40,
            tolerance,
            relaxation: 1.0,
        };

        Problem {
            topology: NetworkTopology {
                incidence: full,
                free_incidence,
                fixed_incidence,
                num_edges: edges.len(),
                num_nodes,
                free_node_indices: free_indices,
                fixed_node_indices: fixed_indices,
            },
            free_node_loads: Array2::zeros(((size - 2) * (size - 2), 3)),
            fixed_node_positions: fixed_positions.clone(),
            anchors: AnchorInfo::all_fixed(fixed_positions),
            objectives: Vec::new(),
            bounds: Bounds {
                lower: vec![1.0; edges.len()],
                upper: vec![1.0e6; edges.len()],
            },
            solver: SolverOptions::default(),
            self_weight: None,
            pressure: Some(pressure),
        }
    }

    fn anchored_grid_problem(rho: f64, tolerance: f64) -> Problem {
        anchored_grid_problem_size(3, rho, tolerance)
    }

    fn current_hydro_load(cache: &mut FdmCache, pressure: &PressureParams) -> Array2<f64> {
        cache.pn.fill(0.0);
        accumulate_hydrostatic_loads_scaled(cache, pressure, 1.0);
        cache.pn.clone()
    }

    #[test]
    fn hydrostatic_forward_jacobian_matches_finite_difference() {
        let problem = anchored_grid_problem(12.0, 1e-9);
        let pressure = problem.pressure.as_ref().unwrap();
        let mut cache = FdmCache::new(&problem).unwrap();
        update_fixed_positions(&mut cache, &problem, &Array2::zeros((0, 3)));
        cache.x[[0, 0]] = 1.05;
        cache.x[[0, 1]] = 0.93;
        cache.x[[0, 2]] = -0.2;
        sync_free_positions(&mut cache, &problem);

        let direction = Array2::from_shape_vec((1, 3), vec![0.2, -0.3, 0.4]).unwrap();
        let mut analytic = Array2::zeros((1, 3));
        apply_hydrostatic_load_jacobian(&cache, pressure, &direction, 1.0, &mut analytic);

        let base_x = cache.x.clone();
        let epsilon = 1e-6;
        for d in 0..3 {
            cache.x[[0, d]] = base_x[[0, d]] + epsilon * direction[[0, d]];
        }
        sync_free_positions(&mut cache, &problem);
        let plus = current_hydro_load(&mut cache, pressure);
        for d in 0..3 {
            cache.x[[0, d]] = base_x[[0, d]] - epsilon * direction[[0, d]];
        }
        sync_free_positions(&mut cache, &problem);
        let minus = current_hydro_load(&mut cache, pressure);

        for d in 0..3 {
            let finite_difference = (plus[[0, d]] - minus[[0, d]]) / (2.0 * epsilon);
            let error = (analytic[[0, d]] - finite_difference).abs();
            assert!(
                error <= 1e-7 * finite_difference.abs().max(1.0),
                "dimension {d}: analytic={}, finite_difference={}, error={error}",
                analytic[[0, d]],
                finite_difference,
            );
        }
    }

    #[test]
    fn hydrostatic_forward_jacobian_matches_multi_node_finite_difference() {
        let problem = anchored_grid_problem_size(5, 12.0, 1e-9);
        let pressure = problem.pressure.as_ref().unwrap();
        let mut cache = FdmCache::new(&problem).unwrap();
        update_fixed_positions(&mut cache, &problem, &Array2::zeros((0, 3)));
        for (i, &node) in problem.topology.free_node_indices.iter().enumerate() {
            cache.x[[i, 0]] = (node % 5) as f64 + 0.01 * i as f64;
            cache.x[[i, 1]] = (node / 5) as f64 - 0.015 * i as f64;
            cache.x[[i, 2]] = -0.05 * (i + 1) as f64;
        }
        sync_free_positions(&mut cache, &problem);

        let direction = Array2::from_shape_fn(cache.x.raw_dim(), |(i, d)| {
            0.03 * (1 + i * 3 + d) as f64 - 0.2
        });
        let mut analytic = Array2::zeros(cache.x.raw_dim());
        apply_hydrostatic_load_jacobian(&cache, pressure, &direction, 1.0, &mut analytic);

        let base_x = cache.x.clone();
        let epsilon = 1e-6;
        for i in 0..cache.x.nrows() {
            for d in 0..3 {
                cache.x[[i, d]] = base_x[[i, d]] + epsilon * direction[[i, d]];
            }
        }
        sync_free_positions(&mut cache, &problem);
        let plus = current_hydro_load(&mut cache, pressure);
        for i in 0..cache.x.nrows() {
            for d in 0..3 {
                cache.x[[i, d]] = base_x[[i, d]] - epsilon * direction[[i, d]];
            }
        }
        sync_free_positions(&mut cache, &problem);
        let minus = current_hydro_load(&mut cache, pressure);

        for i in 0..cache.x.nrows() {
            for d in 0..3 {
                let finite_difference = (plus[[i, d]] - minus[[i, d]]) / (2.0 * epsilon);
                let error = (analytic[[i, d]] - finite_difference).abs();
                assert!(
                    error <= 2e-7 * finite_difference.abs().max(1.0),
                    "node {i}, dimension {d}: analytic={}, finite_difference={}, error={error}",
                    analytic[[i, d]],
                    finite_difference,
                );
            }
        }
    }

    #[test]
    fn hydrostatic_jacobian_transpose_satisfies_dot_product_identity() {
        let problem = anchored_grid_problem(12.0, 1e-9);
        let pressure = problem.pressure.as_ref().unwrap();
        let mut cache = FdmCache::new(&problem).unwrap();
        update_fixed_positions(&mut cache, &problem, &Array2::zeros((0, 3)));
        cache.x[[0, 0]] = 1.05;
        cache.x[[0, 1]] = 0.93;
        cache.x[[0, 2]] = -0.2;
        sync_free_positions(&mut cache, &problem);

        let direction = Array2::from_shape_vec((1, 3), vec![0.2, -0.3, 0.4]).unwrap();
        let dual = Array2::from_shape_vec((1, 3), vec![-0.7, 0.5, 0.1]).unwrap();
        let mut j_direction = Array2::zeros((1, 3));
        let mut jt_dual = Array2::zeros((1, 3));
        apply_hydrostatic_load_jacobian(&cache, pressure, &direction, 1.0, &mut j_direction);
        crate::gradients::dpn_pressure_dx_transpose_matvec(&cache, pressure, &dual, &mut jt_dual);

        let lhs: f64 = j_direction
            .iter()
            .zip(dual.iter())
            .map(|(a, b)| a * b)
            .sum();
        let rhs: f64 = direction
            .iter()
            .zip(jt_dual.iter())
            .map(|(a, b)| a * b)
            .sum();
        assert!((lhs - rhs).abs() <= 1e-12, "lhs={lhs}, rhs={rhs}");
    }

    #[test]
    fn newton_hydrostatic_solve_reaches_equilibrium_for_fixed_q() {
        let problem = anchored_grid_problem(100.0, 1e-9);
        let q = vec![100.0; problem.topology.num_edges];
        let mut cache = FdmCache::new(&problem).unwrap();
        solve_fdm_with_loads(&mut cache, &q, &problem, &Array2::zeros((0, 3)), 1e-12).unwrap();

        let pressure = problem.pressure.as_ref().unwrap();
        let (residual, load_norm) =
            evaluate_hydrostatic_residual(&mut cache, &problem, pressure, 1.0).unwrap();
        assert!(cache.x.iter().all(|v| v.is_finite()));
        assert!(cache.x[[0, 2]] < 0.0);
        assert!(vector_norm(&residual) / load_norm.max(1.0) <= 1e-9);
    }

    #[test]
    fn newton_gmres_solves_multi_node_anchored_grid() {
        let problem = anchored_grid_problem_size(5, 40.0, 1e-8);
        let q = vec![100.0; problem.topology.num_edges];
        let mut cache = FdmCache::new(&problem).unwrap();
        solve_fdm_with_loads(&mut cache, &q, &problem, &Array2::zeros((0, 3)), 1e-12).unwrap();

        let pressure = problem.pressure.as_ref().unwrap();
        let (residual, load_norm) =
            evaluate_hydrostatic_residual(&mut cache, &problem, pressure, 1.0).unwrap();
        assert!(cache.x.iter().all(|v| v.is_finite()));
        assert!(cache.x.column(2).iter().all(|z| *z < 0.0));
        assert!(vector_norm(&residual) / load_norm.max(1.0) <= 1e-8);
    }

    #[test]
    fn continuation_failure_reports_hydrostatic_context() {
        let problem = anchored_grid_problem(2000.0, 1e-8);
        let q = vec![100.0; problem.topology.num_edges];
        let mut cache = FdmCache::new(&problem).unwrap();
        let error = solve_fdm_with_loads(&mut cache, &q, &problem, &Array2::zeros((0, 3)), 1e-12)
            .unwrap_err();
        assert!(error.to_string().contains("Hydrostatic follower"));
    }
}

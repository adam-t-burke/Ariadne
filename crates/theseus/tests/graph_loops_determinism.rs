//! The rewritten O(ne)/O(nn) loops (`fdm`, `gradients`, `objectives`):
//! bitwise identical results under rayon pools of 1, 2 and 4 threads, and
//! equal to straightforward sequential reference implementations (bitwise
//! where the summation order is unchanged, 1e-12 relative for the chunked
//! loss reductions).

#[path = "support/grid.rs"]
#[allow(dead_code)]
mod grid;

use ndarray::Array2;
use theseus::types::*;

/// Grid fixture with objectives exercising every rewritten loop: node lists
/// equal to the free-node list, a distinct but unsorted node list, a strided
/// node list, edge objectives and a reaction objective.
fn fixture(n: usize) -> (Problem, Vec<f64>) {
    let mut problem = grid::make_grid_problem(n);
    let ne = problem.topology.num_edges;
    let free = problem.topology.free_node_indices.clone();
    let mut rev = free.clone();
    rev.reverse();
    problem.objectives.push(Box::new(TargetXY {
        weight: 0.3,
        node_indices: rev.clone(),
        target: Array2::from_shape_fn((rev.len(), 3), |(i, d)| 0.01 * (i * 3 + d) as f64),
        reduction: TargetGeometryReduction::Mse,
    }));
    let strided: Vec<usize> = free.iter().copied().step_by(2).collect();
    problem.objectives.push(Box::new(TargetPlane {
        weight: 0.4,
        target: Array2::zeros((strided.len(), 3)),
        node_indices: strided,
        origin: [0.1, 0.2, 0.3],
        x_axis: [1.0, 0.0, 0.0],
        y_axis: [0.0, 0.0, 1.0],
        reduction: TargetGeometryReduction::Rmse,
    }));
    let sub: Vec<usize> = (0..ne).step_by(3).collect();
    problem.objectives.push(Box::new(TargetLength {
        weight: 0.7,
        target: vec![0.9; sub.len()],
        edge_indices: sub,
    }));
    problem.objectives.push(Box::new(SumForceLength {
        weight: 0.1,
        edge_indices: (0..ne).collect(),
    }));
    problem.objectives.push(Box::new(LengthVariation {
        weight: 0.2,
        edge_indices: (0..ne).collect(),
        sharpness: 3.0,
        use_normalized_variance: false,
        normalization_strategy: LengthVarianceNormalizationStrategy::SquaredMean,
    }));
    problem.objectives.push(Box::new(ReactionDirection {
        weight: 0.5,
        anchor_indices: problem.topology.fixed_node_indices.clone(),
        target_directions: Array2::from_shape_fn((4, 3), |(_, d)| f64::from(d == 2)),
    }));
    let q: Vec<f64> = (0..ne)
        .map(|k| 1.0 + 0.5 * ((k as f64) * 0.37).sin())
        .collect();
    (problem, q)
}

struct Outputs {
    loss: f64,
    grad: Vec<f64>,
    a_values: Vec<f64>,
    rhs: Vec<f64>,
    lengths: Vec<f64>,
    forces: Vec<f64>,
    reactions: Vec<f64>,
    grad_x: Vec<f64>,
    grad_nf: Vec<f64>,
}

fn evaluate(problem: &Problem, q: &[f64]) -> (FdmCache, Outputs) {
    let ne = q.len();
    let mut cache = FdmCache::new(problem).unwrap();
    let anchors = Array2::zeros((0, 3));
    theseus::fdm::solve_fdm(&mut cache, q, problem, &anchors, 1e-12).unwrap();
    let lb = vec![0.1; ne];
    let ub = vec![10.0; ne];
    let idx: Vec<usize> = (0..ne).collect();
    let mut grad = vec![0.0; ne];
    let loss = theseus::gradients::value_and_gradient(
        &mut cache, problem, q, &mut grad, &lb, &ub, &idx, &idx,
    )
    .unwrap();
    // `factor_and_solve` adds the 1e-12 diagonal perturbation in place;
    // re-assemble so the stored values are the pure `assemble_a` output.
    theseus::fdm::assemble_a(&mut cache);
    let out = Outputs {
        loss,
        grad,
        a_values: cache.a_matrix().unwrap().values.clone(),
        rhs: cache.rhs.as_slice().unwrap().to_vec(),
        lengths: cache.member_lengths.clone(),
        forces: cache.member_forces.clone(),
        reactions: cache.reactions.as_slice().unwrap().to_vec(),
        grad_x: cache.grad_x.as_slice().unwrap().to_vec(),
        grad_nf: cache.grad_nf.as_slice().unwrap().to_vec(),
    };
    (cache, out)
}

fn assert_bitwise(actual: &[f64], expected: &[f64], what: &str) {
    assert_eq!(actual.len(), expected.len(), "{what}: length");
    for (i, (a, e)) in actual.iter().zip(expected).enumerate() {
        assert_eq!(a.to_bits(), e.to_bits(), "{what}: entry {i}: {a} vs {e}");
    }
}

fn assert_rel(actual: f64, expected: f64, tol: f64, what: &str) {
    let scale = expected.abs().max(1e-300);
    assert!(
        (actual - expected).abs() <= tol * scale,
        "{what}: {actual} vs {expected} (rel {:e})",
        (actual - expected).abs() / scale
    );
}

/// 160 × 160 grid (50,880 edges, 25,596 free nodes): every loop is above
/// the parallel threshold.
#[test]
fn evaluation_is_bitwise_identical_across_thread_counts() {
    let (problem, q) = fixture(160);
    let run = |threads: usize| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap()
            .install(|| evaluate(&problem, &q).1)
    };
    let base = run(1);
    for threads in [2, 4] {
        let other = run(threads);
        let tag = format!("{threads} threads");
        assert_eq!(base.loss.to_bits(), other.loss.to_bits(), "loss, {tag}");
        assert_bitwise(&other.grad, &base.grad, &format!("gradient, {tag}"));
        assert_bitwise(&other.a_values, &base.a_values, &format!("A values, {tag}"));
        assert_bitwise(&other.rhs, &base.rhs, &format!("rhs, {tag}"));
        assert_bitwise(&other.lengths, &base.lengths, &format!("lengths, {tag}"));
        assert_bitwise(&other.forces, &base.forces, &format!("forces, {tag}"));
        assert_bitwise(
            &other.reactions,
            &base.reactions,
            &format!("reactions, {tag}"),
        );
        assert_bitwise(&other.grad_x, &base.grad_x, &format!("grad_x, {tag}"));
        assert_bitwise(&other.grad_nf, &base.grad_nf, &format!("grad_nf, {tag}"));
    }
}

/// Sequential edge-scatter references for the loops whose summation order
/// the gather form preserves exactly.
#[test]
fn loops_match_sequential_scatter_references_bitwise() {
    for n in [24usize, 160] {
        let (problem, q) = fixture(n);
        let (cache, out) = evaluate(&problem, &q);
        let ne = q.len();
        let nn = problem.topology.num_nodes;
        let nn_free = problem.topology.free_node_indices.len();
        let nf = cache.nf.as_slice().unwrap();
        let lambda = cache.lambda.as_slice().unwrap();
        let (starts, ends) = (&cache.edge_starts, &cache.edge_ends);

        // Geometry: one sequential pass, reactions scattered to endpoints.
        let mut lengths = vec![0.0; ne];
        let mut forces = vec![0.0; ne];
        let mut reactions = vec![0.0; nn * 3];
        for k in 0..ne {
            let (s, e) = (starts[k], ends[k]);
            let d = [
                nf[e * 3] - nf[s * 3],
                nf[e * 3 + 1] - nf[s * 3 + 1],
                nf[e * 3 + 2] - nf[s * 3 + 2],
            ];
            let len = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).max(0.0).sqrt();
            lengths[k] = len;
            forces[k] = q[k] * len;
            for c in 0..3 {
                reactions[s * 3 + c] += d[c] * q[k];
                reactions[e * 3 + c] -= d[c] * q[k];
            }
        }
        assert_bitwise(&out.lengths, &lengths, "lengths");
        assert_bitwise(&out.forces, &forces, "forces");
        assert_bitwise(&out.reactions, &reactions, "reactions");

        // Right-hand side: pn plus every boundary edge in edge order.
        let pn = cache.pn.as_slice().unwrap();
        let mut rhs = pn.to_vec();
        for k in 0..ne {
            let (s, e) = (starts[k], ends[k]);
            match (cache.node_to_free_idx[s], cache.node_to_free_idx[e]) {
                (Some(f), None) => {
                    for c in 0..3 {
                        rhs[f * 3 + c] += q[k] * nf[e * 3 + c];
                    }
                }
                (None, Some(f)) => {
                    for c in 0..3 {
                        rhs[f * 3 + c] += q[k] * nf[s * 3 + c];
                    }
                }
                _ => {}
            }
        }
        assert_bitwise(&out.rhs, &rhs, "rhs");

        // A values: sequential gather through the same map.
        let direct = cache.direct_solver().unwrap();
        let map = direct.q_to_nz();
        let a_values: Vec<f64> = (0..direct.a_matrix().values.len())
            .map(|nz| {
                let range = map.nz_offsets[nz]..map.nz_offsets[nz + 1];
                map.edge[range.clone()]
                    .iter()
                    .zip(&map.coeff[range])
                    .map(|(&k, &c)| q[k as usize] * c)
                    .sum()
            })
            .collect();
        assert_bitwise(&out.a_values, &a_values, "A values");

        // Implicit gradient: the sequential edge loop on top of the explicit
        // part. Recompute the explicit gradients alone to get the base.
        let mut cache2 = FdmCache::new(&problem).unwrap();
        cache2.q.copy_from_slice(&q);
        cache2.nf.assign(&cache.nf);
        cache2.x.assign(&cache.x);
        cache2.lambda.assign(&cache.lambda);
        cache2.member_lengths.copy_from_slice(&cache.member_lengths);
        cache2.member_forces.copy_from_slice(&cache.member_forces);
        cache2.reactions.assign(&cache.reactions);
        cache2.grad_q.fill(0.0);
        cache2.grad_nf.fill(0.0);
        theseus::gradients::accumulate_explicit_gradients(&mut cache2, &problem);
        assert_bitwise(cache2.grad_x.as_slice().unwrap(), &out.grad_x, "grad_x");
        let mut grad_q = cache2.grad_q.clone();
        let mut grad_nf = cache2.grad_nf.as_slice().unwrap().to_vec();
        for k in 0..ne {
            let (u, v) = (starts[k], ends[k]);
            let (u_free, v_free) = (cache.node_to_free_idx[u], cache.node_to_free_idx[v]);
            for c in 0..3 {
                let lam_u = u_free.map_or(0.0, |uf| lambda[uf * 3 + c]);
                let lam_v = v_free.map_or(0.0, |vf| lambda[vf * 3 + c]);
                let d_lam = lam_v - lam_u;
                let d_n = nf[v * 3 + c] - nf[u * 3 + c];
                grad_q[k] -= d_lam * d_n;
                let term = -q[k] * d_lam;
                if v_free.is_none() {
                    grad_nf[v * 3 + c] += term;
                }
                if u_free.is_none() {
                    grad_nf[u * 3 + c] -= term;
                }
            }
        }
        assert_bitwise(&out.grad_nf, &grad_nf, "grad_nf");
        // grad_q also carries the barrier term; compare the pre-barrier part
        // through the cache's own grad_q.
        assert_bitwise(&cache.grad_q, &grad_q, "grad_q (implicit + explicit)");

        // Explicit node-position gradient of the first objective (TargetXYZ on
        // every free node, weight 1, target (col, row, -0.2) as built by
        // `support/grid.rs`) against the sequential scatter.
        let mut cache3 = FdmCache::new(&problem).unwrap();
        cache3.nf.assign(&cache.nf);
        let mut only = grid::make_grid_problem(n);
        only.objectives.truncate(1);
        cache3.grad_q.fill(0.0);
        cache3.grad_nf.fill(0.0);
        theseus::gradients::accumulate_explicit_gradients(&mut cache3, &only);
        let free = &problem.topology.free_node_indices;
        let mut gx = vec![0.0; nn_free * 3];
        for (i, &node) in free.iter().enumerate() {
            let target = [(node % n) as f64, (node / n) as f64, -0.2];
            for c in 0..3 {
                gx[i * 3 + c] += 2.0 * 1.0 * (nf[node * 3 + c] - target[c]);
            }
        }
        assert_bitwise(cache3.grad_x.as_slice().unwrap(), &gx, "TargetXYZ grad_x");
    }
}

/// Chunked loss reductions against plain sequential sums.
#[test]
fn loss_reductions_match_sequential_sums() {
    let (problem, q) = fixture(160);
    let (cache, _) = evaluate(&problem, &q);
    let snap = GeometrySnapshot {
        xyz_full: &cache.nf,
        member_lengths: &cache.member_lengths,
        member_forces: &cache.member_forces,
        reactions: &cache.reactions,
    };
    let ne = q.len();
    let nf = cache.nf.as_slice().unwrap();
    let free = &problem.topology.free_node_indices;

    // Objectives: [0] TargetXYZ, [1] TargetXY, [2] TargetPlane, [3] TargetLength,
    // [4] SumForceLength, [5] LengthVariation, [6] ReactionDirection.
    // SumForceLength (weight 0.1 over all edges).
    let mut sfl = 0.0;
    for k in 0..ne {
        sfl += cache.member_lengths[k] * cache.member_forces[k].abs();
    }
    assert_rel(
        problem.objectives[4].loss(&snap),
        0.1 * sfl,
        1e-12,
        "SumForceLength",
    );

    // TargetLength (weight 0.7, every third edge, target 0.9).
    let mut tl = 0.0;
    for k in (0..ne).step_by(3) {
        let diff = cache.member_lengths[k] - 0.9;
        tl += diff * diff;
    }
    assert_rel(
        problem.objectives[3].loss(&snap),
        0.7 * tl,
        1e-12,
        "TargetLength",
    );

    // TargetXY, Mse, reversed free nodes.
    let mut sse = 0.0;
    let n = free.len();
    for (i, &node) in free.iter().rev().enumerate() {
        for d in 0..2 {
            let diff = nf[node * 3 + d] - 0.01 * (i * 3 + d) as f64;
            sse += diff * diff;
        }
    }
    assert_rel(
        problem.objectives[1].loss(&snap),
        0.3 * sse / n as f64,
        1e-12,
        "TargetXY",
    );

    // LengthVariation: smooth max − smooth min with β = 3.
    let beta = 3.0;
    let m = cache
        .member_lengths
        .iter()
        .fold(f64::NEG_INFINITY, |a, &b| a.max(b));
    let lo = cache
        .member_lengths
        .iter()
        .fold(f64::INFINITY, |a, &b| a.min(b));
    let smax = m + cache
        .member_lengths
        .iter()
        .map(|l| ((l - m) * beta).exp())
        .sum::<f64>()
        .ln()
        / beta;
    let smin = lo
        - cache
            .member_lengths
            .iter()
            .map(|l| ((lo - l) * beta).exp())
            .sum::<f64>()
            .ln()
            / beta;
    assert_rel(
        problem.objectives[5].loss(&snap),
        0.2 * (smax - smin),
        1e-12,
        "LengthVariation",
    );

    // total_loss is the sum of the objective losses in order.
    let total = theseus::objectives::total_loss(&problem.objectives, &snap);
    let seq = problem
        .objectives
        .iter()
        .map(|o| o.loss(&snap))
        .fold(0.0, |a, l| a + l);
    assert_eq!(total.to_bits(), seq.to_bits(), "total_loss order");
}

/// End to end: a full L-BFGS-B solve on the 160 × 160 grid (fixture with
/// every rewritten loop above the parallel threshold) produces bitwise the
/// same force densities, geometry and loss trace on pools of 1 and 4
/// threads.
#[test]
fn full_optimization_is_bitwise_identical_across_thread_counts() {
    use std::sync::atomic::AtomicBool;

    let (problem, _) = fixture(160);
    let mut problem = problem;
    problem.solver.q_parameterization_mode = QParameterizationMode::DirectBoxBounds;
    problem.solver.absolute_tolerance = 0.0;
    problem.solver.relative_tolerance = 0.0;
    problem.solver.max_iterations = 8;
    let ne = problem.topology.num_edges;

    let run = |threads: usize| -> SolverResult {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap()
            .install(|| {
                let cancel = AtomicBool::new(false);
                let mut state = OptimizationState::new(vec![1.0; ne], Array2::zeros((0, 3)));
                theseus::optimizer::optimize(&problem, &mut state, None, 1, &cancel).unwrap()
            })
    };
    let base = run(1);
    let other = run(4);
    assert!(base.loss_trace.len() > 1, "optimizer should take steps");
    assert_eq!(base.iterations, other.iterations, "iteration count");
    assert_bitwise(&other.loss_trace, &base.loss_trace, "loss trace");
    assert_bitwise(&other.q, &base.q, "final q");
    assert_bitwise(&other.member_lengths, &base.member_lengths, "final lengths");
    assert_bitwise(&other.member_forces, &base.member_forces, "final forces");
    assert_bitwise(
        other.xyz.as_slice().unwrap(),
        base.xyz.as_slice().unwrap(),
        "final xyz",
    );
    assert_bitwise(
        other.reactions.as_slice().unwrap(),
        base.reactions.as_slice().unwrap(),
        "final reactions",
    );
}

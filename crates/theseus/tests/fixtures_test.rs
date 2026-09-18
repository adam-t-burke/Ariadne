//! Unit tests of the regression/benchmark fixture generators in
//! `tests/support/fixtures/`: node and edge counts, symmetric simple
//! connectivity, at least one support per connected component, matched-size
//! helpers, determinism, and one successful FDM evaluation through the public
//! API on every `make_*` problem.

#[path = "support/fixtures/mod.rs"]
mod fixtures;

use fixtures::harness::{build_fixture, FixtureKind, Timing};
use fixtures::network::Network;
use fixtures::*;
use ndarray::Array2;
use theseus::types::*;

fn check_network(net: &Network) {
    assert!(
        net.is_symmetric_simple(),
        "edges must be unique with start < end"
    );
    for &(s, t) in &net.edges {
        assert!(t < net.num_nodes(), "edge endpoint out of range");
        assert_ne!(s, t);
    }
    assert!(!net.fixed.is_empty());
    assert!(
        net.fixed.windows(2).all(|w| w[0] < w[1]),
        "fixed sorted, unique"
    );
    assert!(*net.fixed.last().unwrap() < net.num_nodes());
    assert_eq!(
        net.unsupported_components(),
        0,
        "every connected component needs a fixed node"
    );
    assert!(net.degrees().iter().all(|&d| d > 0), "no isolated nodes");
}

/// One forward solve (`FdmCache::new` + `solve_fdm`) and one fused
/// objective/gradient evaluation at `q = 1`; the loss must be finite and the
/// solution at `q*` must have (near) zero loss. Box-bounds mode, as in the
/// benchmarks, so no soft-bound penalty is added to the target loss.
fn check_problem_solves(problem: &mut Problem, q_star: Option<&[f64]>) {
    problem.solver.q_parameterization_mode = QParameterizationMode::DirectBoxBounds;
    let problem = &*problem;
    let ne = problem.topology.num_edges;
    let anchors = Array2::zeros((0, 3));
    let q = vec![1.0; ne];
    let mut cache = FdmCache::new(problem).unwrap();
    theseus::fdm::solve_fdm(&mut cache, &q, problem, &anchors, 0.0).unwrap();
    assert!(cache.nf.iter().all(|v| v.is_finite()));

    let lb = problem.bounds.lower.clone();
    let ub = problem.bounds.upper.clone();
    let idx: Vec<usize> = (0..ne).collect();
    let mut grad = vec![0.0; ne];
    let loss = theseus::gradients::value_and_gradient(
        &mut cache, problem, &q, &mut grad, &lb, &ub, &idx, &idx,
    )
    .unwrap();
    assert!(
        loss.is_finite() && loss > 0.0,
        "loss at q = 1 should be positive: {loss}"
    );
    assert!(grad.iter().all(|g| g.is_finite()));

    if let Some(q_star) = q_star {
        // No soft-bound barrier terms: the geometric loss alone must vanish at q*.
        let loss_star = theseus::gradients::value_and_gradient(
            &mut cache,
            problem,
            q_star,
            &mut grad,
            &lb,
            &ub,
            &[],
            &[],
        )
        .unwrap();
        assert!(
            loss_star.abs() <= 1e-8 * loss.max(1.0),
            "loss at q* should vanish: {loss_star} vs {loss} at q = 1"
        );
    }
}

#[test]
fn grid_network_matches_grid_rs() {
    let net = grid_network(9);
    check_network(&net);
    assert_eq!(net.num_nodes(), 81);
    assert_eq!(net.num_edges(), grid::grid_edges(9));
    assert_eq!(net.fixed, vec![0, 8, 72, 80]);
    let problem = grid::make_grid_problem(9);
    assert_eq!(problem.topology.num_edges, net.num_edges());
    assert_eq!(problem.topology.fixed_node_indices, net.fixed);
}

#[test]
fn irregular_mesh_counts_and_connectivity() {
    let side = 24;
    let net = irregular_mesh(side, 7);
    check_network(&net);
    assert_eq!(net.num_nodes(), side * side);
    assert_eq!(net.fixed.len(), 4 * (side - 1), "boundary ring fixed");
    assert_eq!(net.components().len(), 1, "one connected component");
    // Every interior node asked for 5–7 neighbours; the union can only add.
    for (i, &d) in net.degrees().iter().enumerate() {
        if !net.fixed.contains(&i) {
            assert!(d >= irregular::MIN_K, "node {i} has degree {d}");
        }
    }
    let per_point = net.num_edges() as f64 / net.num_nodes() as f64;
    assert!(
        (per_point - irregular::EDGES_PER_POINT).abs() < 0.15,
        "edges per point {per_point} drifted from the calibration constant"
    );
    // Jitter bounded and boundary on the lattice.
    for (i, p) in net.positions.iter().enumerate() {
        let (c, r) = ((i % side) as f64, (i / side) as f64);
        assert!((p[0] - c).abs() <= irregular::JITTER + 1e-12);
        assert!((p[1] - r).abs() <= irregular::JITTER + 1e-12);
        if net.fixed.contains(&i) {
            assert_eq!((p[0], p[1]), (c, r));
        }
    }
}

#[test]
fn irregular_mesh_is_deterministic_and_seed_dependent() {
    let a = irregular_mesh(16, 1);
    let b = irregular_mesh(16, 1);
    let c = irregular_mesh(16, 2);
    assert_eq!(a.edges, b.edges);
    assert_eq!(a.positions, b.positions);
    assert_ne!(a.positions, c.positions);
}

#[test]
fn irregular_mesh_problem_solves() {
    let net = irregular_mesh(14, 3);
    let q_star = smooth_q_star(net.num_edges());
    let mut problem = make_irregular_mesh_problem(14, 3);
    assert_eq!(problem.topology.num_edges, net.num_edges());
    assert_eq!(problem.topology.free_node_indices.len(), net.num_free());
    check_problem_solves(&mut problem, Some(&q_star));
}

#[test]
fn cable_dome_counts_and_degrees() {
    let (rings, spokes) = (6, 12);
    let net = cable_dome(rings, spokes);
    check_network(&net);
    assert_eq!(net.num_nodes(), 1 + rings * spokes);
    assert_eq!(net.num_edges(), dome::dome_edges(rings, spokes));
    assert_eq!(net.fixed.len(), spokes, "outer ring fixed");
    assert_eq!(net.components().len(), 1);
    let deg = net.degrees();
    assert_eq!(deg[0], spokes, "hub has one edge per spoke");
    assert!((8..=12).contains(&deg[0]));
    for r in 1..=rings {
        for s in 0..spokes {
            let node = 1 + (r - 1) * spokes + s;
            let expected = if r == 1 {
                6
            } else if r == rings {
                5
            } else {
                8
            };
            assert_eq!(deg[node], expected, "ring {r} node degree");
        }
    }
}

#[test]
fn cable_dome_problem_solves() {
    let net = cable_dome(5, 10);
    let q_star = smooth_q_star(net.num_edges());
    let mut problem = make_cable_dome_problem(5, 10);
    assert_eq!(problem.topology.fixed_node_indices.len(), 10);
    check_problem_solves(&mut problem, Some(&q_star));
}

#[test]
fn few_supports_counts() {
    let side = 12;
    let net = few_supports(side);
    check_network(&net);
    assert_eq!(net.num_nodes(), side * side);
    assert_eq!(net.num_edges(), grid::grid_edges(side));
    assert_eq!(net.fixed, vec![0, side * side - 1]);
    for supports in 2..=4 {
        let n = few_supports_with(side, supports);
        check_network(&n);
        assert_eq!(n.fixed.len(), supports);
    }
}

#[test]
fn few_supports_problem_solves() {
    let net = few_supports(12);
    let q_star = smooth_q_star(net.num_edges());
    let mut problem = make_few_supports_problem(12);
    assert_eq!(problem.topology.fixed_node_indices.len(), 2);
    check_problem_solves(&mut problem, Some(&q_star));
}

#[test]
fn disconnected_components_each_supported() {
    let (side, k) = (7, 4);
    let net = disconnected(side, k);
    check_network(&net);
    assert_eq!(net.num_nodes(), k * side * side);
    assert_eq!(net.num_edges(), k * grid::grid_edges(side));
    assert_eq!(net.fixed.len(), 4 * k);
    let comps = net.components();
    assert_eq!(comps.len(), k);
    for c in &comps {
        assert_eq!(c.len(), side * side);
        assert_eq!(c.iter().filter(|u| net.fixed.contains(u)).count(), 4);
    }
    // Unsupported component detection works: drop the supports of one block.
    let mut broken = net.clone();
    broken.fixed.retain(|&f| f >= side * side);
    assert_eq!(broken.unsupported_components(), 1);
}

#[test]
fn disconnected_problem_solves() {
    let net = disconnected(7, 3);
    let q_star = smooth_q_star(net.num_edges());
    let mut problem = make_disconnected_problem(7, 3);
    assert_eq!(problem.topology.num_nodes, 3 * 49);
    check_problem_solves(&mut problem, Some(&q_star));
}

#[test]
fn anisotropic_q_spans_ratio_and_solves() {
    for &ratio in &[100.0, 1e4] {
        let side = 12;
        let q_star = anisotropic_q_star(side, ratio);
        assert_eq!(q_star.len(), grid::grid_edges(side));
        let (lo, hi) = q_star.iter().fold((f64::INFINITY, 0.0f64), |(lo, hi), &q| {
            (lo.min(q), hi.max(q))
        });
        // Smooth profile spans 7× on its own; the ramp adds the requested ratio.
        assert!(hi / lo >= ratio, "q* span {} < ratio {ratio}", hi / lo);
        assert!(hi / lo <= ratio * 7.5);
        let mut problem = anisotropic_q(side, ratio);
        assert!(problem
            .bounds
            .lower
            .iter()
            .zip(&q_star)
            .all(|(l, q)| l <= q));
        assert!(problem
            .bounds
            .upper
            .iter()
            .zip(&q_star)
            .all(|(u, q)| u >= q));
        assert!(problem.bounds.lower[0] <= 1.0 && problem.bounds.upper[0] >= 1.0);
        check_problem_solves(&mut problem, Some(&q_star));
    }
}

#[test]
fn matched_sizes_track_grid_edge_counts() {
    for &side in &[72usize, 160, 224] {
        let target = grid::grid_edges(side);
        let m = irregular::side_for_edges(target);
        let ne = irregular_mesh(m, 1).num_edges();
        let rel = (ne as f64 - target as f64).abs() / target as f64;
        assert!(
            rel < 0.1,
            "irregular side {m}: {ne} edges vs target {target}"
        );
        let (rings, spokes) = dome::size_for_edges(target);
        let ne = dome::dome_edges(rings, spokes);
        let rel = (ne as f64 - target as f64).abs() / target as f64;
        assert!(
            rel < 0.02,
            "dome {rings}×{spokes}: {ne} edges vs target {target}"
        );
    }
}

#[test]
fn build_fixture_by_kind() {
    for kind in FixtureKind::ALL {
        let mut f = build_fixture(kind, 10);
        assert_eq!(f.kind, kind);
        assert_eq!(FixtureKind::parse(kind.name()), Some(kind));
        assert!(f.parameters.contains_key("grid_side"));
        let ne = f.problem.topology.num_edges;
        assert!(ne > 0);
        check_problem_solves(&mut f.problem, None);
    }
    assert_eq!(
        FixtureKind::parse("Few_Supports"),
        Some(FixtureKind::FewSupports)
    );
    assert_eq!(FixtureKind::parse("nope"), None);
}

#[test]
fn large_fixtures_build_quickly() {
    // ~200k edges each; the neighbour search must be bucketed, not quadratic.
    let t = std::time::Instant::now();
    let net = irregular_mesh(240, 11);
    assert!(net.num_edges() > 150_000);
    check_network(&net);
    let (rings, spokes) = dome::size_for_edges(200_000);
    let net = cable_dome(rings, spokes);
    assert!(net.num_edges() > 150_000);
    assert_eq!(net.unsupported_components(), 0);
    assert!(
        t.elapsed().as_secs_f64() < 20.0,
        "fixture generation too slow: {:?}",
        t.elapsed()
    );
}

#[test]
fn timing_statistics() {
    let t = Timing::from_samples(vec![5.0, 1.0, 3.0, 2.0, 4.0]);
    assert_eq!(t.median, 3.0);
    assert_eq!(t.q25, 2.0);
    assert_eq!(t.q75, 4.0);
    assert_eq!(t.iqr(), 2.0);
    assert_eq!((t.min, t.max), (1.0, 5.0));
    let one = Timing::from_samples(vec![2.5]);
    assert_eq!(one.median, 2.5);
    assert_eq!(one.iqr(), 0.0);
}

//! Cable dome / spoke wheel: a central hub, `rings` concentric rings of
//! `spokes` nodes each, radial spokes, ring cables and X-bracing between
//! consecutive rings. The outer ring is fixed.
//!
//! Degrees: the hub has `spokes` incident edges (one per spoke of ring 1);
//! every interior ring node has 8 (two ring, two radial, four diagonals);
//! ring-1 nodes have 6 (two ring, hub, outward radial, two outward diagonals);
//! the fixed outer ring nodes have 5. Total edges: `spokes · (4·rings − 2)`.
//! With the matched sizing of [`size_for_edges`] the hub is a very
//! high-degree node (`spokes ≈ 2·rings`), which is the stress case a
//! spoke-wheel poses to aggregation; the unit test also covers a 12-spoke dome
//! whose hub has 12 incident edges.

use super::network::{smooth_q_star, Network};
use theseus::types::Problem;

/// Concentric-ring dome with a hub; `rings ≥ 2`, `spokes ≥ 3`.
pub fn cable_dome(rings: usize, spokes: usize) -> Network {
    assert!(rings >= 2, "dome needs rings >= 2");
    assert!(spokes >= 3, "dome needs spokes >= 3");
    let n = 1 + rings * spokes;
    let node = |r: usize, s: usize| 1 + (r - 1) * spokes + (s % spokes);
    let mut positions = Vec::with_capacity(n);
    positions.push([0.0, 0.0, 0.0]);
    let outer = rings as f64;
    for r in 1..=rings {
        let radius = r as f64;
        // Shallow dome: the reference rises towards the centre.
        let z = 0.25 * (outer * outer - radius * radius) / outer;
        for s in 0..spokes {
            let a = std::f64::consts::TAU * s as f64 / spokes as f64;
            positions.push([radius * a.cos(), radius * a.sin(), z]);
        }
    }
    let mut edges = Vec::with_capacity(spokes * (4 * rings - 2));
    for s in 0..spokes {
        edges.push((0, node(1, s)));
    }
    for r in 1..=rings {
        for s in 0..spokes {
            edges.push((node(r, s), node(r, s + 1)));
            if r < rings {
                edges.push((node(r, s), node(r + 1, s)));
                edges.push((node(r, s), node(r + 1, s + 1)));
                edges.push((node(r + 1, s), node(r, s + 1)));
            }
        }
    }
    let fixed = (0..spokes).map(|s| node(rings, s)).collect();
    Network {
        positions,
        edges,
        fixed,
    }
    .normalised()
}

/// Number of edges of [`cable_dome`].
pub fn dome_edges(rings: usize, spokes: usize) -> usize {
    spokes * (4 * rings - 2)
}

/// Recoverable optimisation problem on the dome.
pub fn make_cable_dome_problem(rings: usize, spokes: usize) -> Problem {
    let net = cable_dome(rings, spokes);
    let q_star = smooth_q_star(net.num_edges());
    net.into_recoverable_problem(&q_star)
}

/// `(rings, spokes)` with `spokes ≈ 2·rings` (outer arc spacing about π times
/// the ring spacing) whose edge count is close to `target_edges`.
pub fn size_for_edges(target_edges: usize) -> (usize, usize) {
    let rings = ((target_edges as f64 / 8.0).sqrt().round() as usize).max(2);
    let spokes = ((target_edges as f64 / (4 * rings - 2) as f64).round() as usize).max(3);
    (rings, spokes)
}

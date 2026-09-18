//! Irregular triangulated-like mesh: a `side × side` lattice of points jittered
//! by up to `JITTER` of the spacing, each connected to its 5–7 nearest
//! neighbours (count drawn per node from the seeded PRNG) and symmetrised, with
//! the outer ring of the lattice fixed.
//!
//! The neighbour search uses the lattice cells as buckets (a jittered point
//! never leaves its own cell), so construction is `O(n)` in the number of
//! points plus the final `O(ne log ne)` edge sort.

use super::network::{smooth_q_star, Network};
use super::rng::Rng;
use theseus::types::Problem;

/// Maximum displacement of a point from its lattice position, in units of the
/// lattice spacing (per coordinate).
pub const JITTER: f64 = 0.4;
/// Nearest-neighbour count is drawn uniformly from `MIN_K..=MAX_K` per node.
pub const MIN_K: usize = 5;
pub const MAX_K: usize = 7;
/// Half-width (in lattice cells) of the window searched for neighbours. Points
/// outside a radius-2 window are at least `3 − 2·JITTER = 2.2` away, which is
/// farther than the seventh neighbour on a jittered lattice.
const WINDOW: i64 = 2;

/// Jittered lattice points connected to their nearest neighbours; boundary ring
/// fixed. `side` is the number of lattice points per side (≥ 3).
pub fn irregular_mesh(side: usize, seed: u64) -> Network {
    assert!(side >= 3, "irregular mesh needs side >= 3");
    let mut rng = Rng::new(seed);
    let n = side * side;
    let mut positions = Vec::with_capacity(n);
    for i in 0..n {
        let (col, row) = ((i % side) as f64, (i / side) as f64);
        // Boundary points stay on the lattice so the fixed ring is a clean square.
        let on_ring =
            i % side == 0 || i % side == side - 1 || i / side == 0 || i / side == side - 1;
        let (dx, dy) = if on_ring {
            (0.0, 0.0)
        } else {
            (JITTER * rng.next_signed(), JITTER * rng.next_signed())
        };
        positions.push([col + dx, row + dy, 0.0]);
    }

    let mut edges = Vec::with_capacity(n * MAX_K);
    let mut candidates: Vec<(f64, usize)> =
        Vec::with_capacity(((2 * WINDOW + 1) * (2 * WINDOW + 1)) as usize);
    let s = side as i64;
    for i in 0..n {
        let k = rng.next_range(MIN_K, MAX_K);
        let (ci, ri) = ((i % side) as i64, (i / side) as i64);
        let p = positions[i];
        candidates.clear();
        for dr in -WINDOW..=WINDOW {
            let r = ri + dr;
            if r < 0 || r >= s {
                continue;
            }
            for dc in -WINDOW..=WINDOW {
                let c = ci + dc;
                if c < 0 || c >= s || (dr == 0 && dc == 0) {
                    continue;
                }
                let j = (r * s + c) as usize;
                let q = positions[j];
                let d2 = (p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2);
                candidates.push((d2, j));
            }
        }
        // Deterministic ordering: by distance, then by index.
        candidates.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
        for &(_, j) in candidates.iter().take(k) {
            edges.push((i.min(j), i.max(j)));
        }
    }

    let fixed = (0..n)
        .filter(|&i| i % side == 0 || i % side == side - 1 || i / side == 0 || i / side == side - 1)
        .collect();
    Network {
        positions,
        edges,
        fixed,
    }
    .normalised()
}

/// Recoverable optimisation problem on the irregular mesh (target = equilibrium
/// of the smooth interior `q*`).
pub fn make_irregular_mesh_problem(side: usize, seed: u64) -> Problem {
    let net = irregular_mesh(side, seed);
    let q_star = smooth_q_star(net.num_edges());
    net.into_recoverable_problem(&q_star)
}

/// Mean edges per lattice point of [`irregular_mesh`] (kNN union with `k ∈
/// 5..=7 on a 0.4-jittered lattice), measured over seeds and sizes 60–531
/// (3.36–3.39; the boundary ring lowers it slightly on small sides); used to
/// pick a side that matches a target edge count. The realised count is within
/// a few percent and is what the harness records.
pub const EDGES_PER_POINT: f64 = 3.37;

/// Lattice side whose irregular mesh has roughly `target_edges` edges.
pub fn side_for_edges(target_edges: usize) -> usize {
    ((target_edges as f64 / EDGES_PER_POINT).sqrt().round() as usize).max(3)
}

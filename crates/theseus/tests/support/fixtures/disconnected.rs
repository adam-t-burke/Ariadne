//! Several independent grids in one problem, each with its own four corner
//! supports, laid out side by side along `x`. The force-density matrix is
//! block diagonal; every block is non-singular.

use super::network::{grid_network, smooth_q_star, Network};
use theseus::types::Problem;

/// `components` copies of a `side × side` grid, each with four fixed corners.
pub fn disconnected(side: usize, components: usize) -> Network {
    assert!(components >= 1, "need at least one component");
    let unit = grid_network(side);
    let nn = unit.num_nodes();
    let gap = 1.5 * side as f64;
    let mut positions = Vec::with_capacity(nn * components);
    let mut edges = Vec::with_capacity(unit.num_edges() * components);
    let mut fixed = Vec::with_capacity(unit.fixed.len() * components);
    for c in 0..components {
        let offset = c * nn;
        let dx = c as f64 * gap;
        positions.extend(unit.positions.iter().map(|p| [p[0] + dx, p[1], p[2]]));
        edges.extend(unit.edges.iter().map(|&(s, t)| (s + offset, t + offset)));
        fixed.extend(unit.fixed.iter().map(|&f| f + offset));
    }
    Network {
        positions,
        edges,
        fixed,
    }
    .normalised()
}

/// Recoverable problem over all components at once.
pub fn make_disconnected_problem(side: usize, components: usize) -> Problem {
    let net = disconnected(side, components);
    let q_star = smooth_q_star(net.num_edges());
    net.into_recoverable_problem(&q_star)
}

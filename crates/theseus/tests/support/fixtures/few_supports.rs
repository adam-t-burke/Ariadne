//! Grid with only a few fixed corner nodes. With two supports the Laplacian
//! is anchored at two nodes only, so its smallest eigenvalue is tiny and the
//! system is near-singular — a conditioning stress case for iterative solvers.

use super::network::{grid_network, smooth_q_star, Network};
use theseus::types::Problem;

/// `side × side` grid fixed at two opposite corners.
pub fn few_supports(side: usize) -> Network {
    few_supports_with(side, 2)
}

/// `side × side` grid fixed at `supports ∈ 2..=4` corners: (0,0) and
/// (side−1, side−1) first, then (side−1, 0), then (0, side−1).
pub fn few_supports_with(side: usize, supports: usize) -> Network {
    assert!((2..=4).contains(&supports), "supports must be 2..=4");
    let mut net = grid_network(side);
    let corners = [0, side * side - 1, side - 1, side * (side - 1)];
    net.fixed = corners[..supports].to_vec();
    net.normalised()
}

/// Recoverable problem on the two-support grid.
pub fn make_few_supports_problem(side: usize) -> Problem {
    let net = few_supports(side);
    let q_star = smooth_q_star(net.num_edges());
    net.into_recoverable_problem(&q_star)
}

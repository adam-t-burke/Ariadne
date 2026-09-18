//! Grid whose recoverable force-density field `q*` spans a factor `ratio`
//! across the domain (geometric ramp in `x` on top of the smooth profile), so
//! the assembled Laplacian has strongly varying edge weights. The box bounds
//! are widened by `sqrt(ratio)` on both sides so `q*` and the unit start are
//! both feasible.

use super::network::{grid_network, make_recoverable, smooth_q_star};
use theseus::types::Problem;

/// `q*` for the `side × side` grid of [`grid_network`]: the smooth profile
/// multiplied by `ratio^(x_mid/(side-1) − 1/2)`, where `x_mid` is the edge
/// midpoint's `x` coordinate.
pub fn anisotropic_q_star(side: usize, ratio: f64) -> Vec<f64> {
    assert!(ratio >= 1.0, "ratio must be >= 1");
    let net = grid_network(side);
    let base = smooth_q_star(net.num_edges());
    let span = (side - 1).max(1) as f64;
    net.edges
        .iter()
        .zip(base)
        .map(|(&(s, t), b)| {
            let x_mid = 0.5 * (net.positions[s][0] + net.positions[t][0]);
            b * ratio.powf(x_mid / span - 0.5)
        })
        .collect()
}

/// Recoverable grid problem with the anisotropic `q*` and bounds
/// `[0.1/sqrt(ratio), 10·sqrt(ratio)]`.
pub fn anisotropic_q(side: usize, ratio: f64) -> Problem {
    let q_star = anisotropic_q_star(side, ratio);
    let mut problem = grid_network(side).into_problem();
    let half = ratio.sqrt();
    let ne = problem.topology.num_edges;
    problem.bounds.lower = vec![0.1 / half; ne];
    problem.bounds.upper = vec![10.0 * half; ne];
    make_recoverable(problem, &q_star)
}

/// Alias following the `make_*_problem` naming of the other fixtures.
pub fn make_anisotropic_q_problem(side: usize, ratio: f64) -> Problem {
    anisotropic_q(side, ratio)
}

//! Square cable-net grid fixtures shared by the scaling benchmarks and the
//! phase profiler: `n × n` nodes, `2n(n-1)` edges, four corner supports, unit
//! downward load on every free node and a target-position objective.

use ndarray::Array2;
use theseus::sparse::SparseColMatOwned;
use theseus::types::*;

/// Number of edges of an `n × n` grid.
pub fn grid_edges(n: usize) -> usize {
    2 * n * (n - 1)
}

/// Grid whose target geometry is the equilibrium shape of a smooth, strictly
/// interior force-density field, so the optimum (zero loss) is reachable
/// inside the box bounds and the optimiser exercises a real descent path.
pub fn make_recoverable_grid_problem(n: usize) -> Problem {
    let mut problem = make_grid_problem(n);
    let ne = problem.topology.num_edges;
    let q_star: Vec<f64> = (0..ne)
        .map(|k| {
            let t = k as f64 / ne as f64;
            let tau = std::f64::consts::TAU;
            1.0 + 0.5 * (tau * t).sin() + 0.25 * (2.0 * tau * t).cos()
        })
        .collect();
    let mut cache = FdmCache::new(&problem).unwrap();
    theseus::fdm::solve_fdm(&mut cache, &q_star, &problem, &Array2::zeros((0, 3)), 0.0).unwrap();
    let free = problem.topology.free_node_indices.clone();
    let mut target = Array2::zeros((free.len(), 3));
    for (i, &node) in free.iter().enumerate() {
        for d in 0..3 {
            target[[i, d]] = cache.nf[[node, d]];
        }
    }
    problem.objectives = vec![Box::new(TargetXYZ {
        weight: 1.0,
        node_indices: free,
        target,
        reduction: TargetGeometryReduction::Sse,
    })];
    problem
}

pub fn make_grid_problem(n: usize) -> Problem {
    let num_nodes = n * n;
    let mut edges = Vec::with_capacity(grid_edges(n));
    for row in 0..n {
        for col in 0..(n - 1) {
            edges.push((row * n + col, row * n + col + 1));
        }
    }
    for row in 0..(n - 1) {
        for col in 0..n {
            edges.push((row * n + col, (row + 1) * n + col));
        }
    }
    let ne = edges.len();
    let mut rows = Vec::with_capacity(ne * 2);
    let mut cols = Vec::with_capacity(ne * 2);
    let mut vals = Vec::with_capacity(ne * 2);
    for (e, &(s, t)) in edges.iter().enumerate() {
        rows.extend([e, e]);
        cols.extend([s, t]);
        vals.extend([-1.0, 1.0]);
    }
    let incidence = SparseColMatOwned::from_coo(ne, num_nodes, &rows, &cols, &vals).unwrap();
    let fixed_idx = vec![0, n - 1, n * (n - 1), n * n - 1];
    let free_idx: Vec<usize> = (0..num_nodes).filter(|i| !fixed_idx.contains(i)).collect();
    let nn_free = free_idx.len();
    let topology = NetworkTopology {
        free_incidence: incidence.extract_columns(&free_idx),
        fixed_incidence: incidence.extract_columns(&fixed_idx),
        incidence,
        num_edges: ne,
        num_nodes,
        free_node_indices: free_idx,
        fixed_node_indices: fixed_idx,
    };
    let mut loads = vec![0.0; nn_free * 3];
    for i in 0..nn_free {
        loads[i * 3 + 2] = -1.0;
    }
    let e = (n - 1) as f64;
    let fixed_node_positions = Array2::from_shape_vec(
        (4, 3),
        vec![0.0, 0.0, 0.0, e, 0.0, 0.0, 0.0, e, 0.0, e, e, 0.0],
    )
    .unwrap();
    let target_nodes = topology.free_node_indices.clone();
    let mut target = vec![0.0; nn_free * 3];
    for (i, &node) in target_nodes.iter().enumerate() {
        target[i * 3] = (node % n) as f64;
        target[i * 3 + 1] = (node / n) as f64;
        target[i * 3 + 2] = -0.2;
    }
    Problem {
        topology,
        free_node_loads: Array2::from_shape_vec((nn_free, 3), loads).unwrap(),
        anchors: AnchorInfo::all_fixed(fixed_node_positions.clone()),
        fixed_node_positions,
        objectives: vec![Box::new(TargetXYZ {
            weight: 1.0,
            node_indices: target_nodes,
            target: Array2::from_shape_vec((nn_free, 3), target).unwrap(),
            reduction: TargetGeometryReduction::Sse,
        })],
        bounds: Bounds {
            lower: vec![0.1; ne],
            upper: vec![10.0; ne],
        },
        solver: SolverOptions::default(),
        self_weight: None,
        pressure: None,
    }
}

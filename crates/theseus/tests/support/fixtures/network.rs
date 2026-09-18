//! Plain edge-list description of a cable network and its conversion into a
//! Theseus [`Problem`] with the same conventions as `support/grid.rs`: unit
//! downward load on every free node, box bounds `[0.1, 10]`, a `TargetXYZ`
//! objective over the free nodes.

use ndarray::Array2;
use theseus::sparse::SparseColMatOwned;
use theseus::types::*;

/// Node positions, undirected edges and the fixed (support) nodes.
#[derive(Debug, Clone)]
pub struct Network {
    /// Reference position of every node (fixed nodes use it as the support
    /// location, free nodes as the initial target).
    pub positions: Vec<[f64; 3]>,
    /// Undirected edges `(start, end)` with `start < end`, sorted, no duplicates.
    pub edges: Vec<(usize, usize)>,
    /// Sorted indices of the fixed nodes.
    pub fixed: Vec<usize>,
}

impl Network {
    /// Normalises the edge list (orders endpoints, sorts, removes duplicates and
    /// self-loops) and the fixed set (sorts, removes duplicates).
    pub fn normalised(mut self) -> Self {
        for e in &mut self.edges {
            if e.0 > e.1 {
                *e = (e.1, e.0);
            }
        }
        self.edges.retain(|e| e.0 != e.1);
        self.edges.sort_unstable();
        self.edges.dedup();
        self.fixed.sort_unstable();
        self.fixed.dedup();
        self
    }

    pub fn num_nodes(&self) -> usize {
        self.positions.len()
    }

    pub fn num_edges(&self) -> usize {
        self.edges.len()
    }

    pub fn num_free(&self) -> usize {
        self.num_nodes() - self.fixed.len()
    }

    /// Number of incident edges per node.
    pub fn degrees(&self) -> Vec<usize> {
        let mut deg = vec![0usize; self.num_nodes()];
        for &(s, t) in &self.edges {
            deg[s] += 1;
            deg[t] += 1;
        }
        deg
    }

    /// `true` when every edge appears exactly once with `start < end`, i.e. the
    /// connectivity is symmetric and free of duplicates.
    pub fn is_symmetric_simple(&self) -> bool {
        let mut sorted = self.edges.clone();
        sorted.sort_unstable();
        sorted.dedup();
        sorted.len() == self.edges.len() && self.edges.iter().all(|&(s, t)| s < t)
    }

    /// Connected components (by node index sets) via breadth-first search.
    pub fn components(&self) -> Vec<Vec<usize>> {
        let n = self.num_nodes();
        let mut offsets = vec![0usize; n + 1];
        for &(s, t) in &self.edges {
            offsets[s + 1] += 1;
            offsets[t + 1] += 1;
        }
        for i in 0..n {
            offsets[i + 1] += offsets[i];
        }
        let mut fill = offsets.clone();
        let mut adj = vec![0usize; offsets[n]];
        for &(s, t) in &self.edges {
            adj[fill[s]] = t;
            fill[s] += 1;
            adj[fill[t]] = s;
            fill[t] += 1;
        }
        let mut seen = vec![false; n];
        let mut comps = Vec::new();
        let mut queue = Vec::new();
        for root in 0..n {
            if seen[root] {
                continue;
            }
            seen[root] = true;
            queue.clear();
            queue.push(root);
            let mut head = 0;
            while head < queue.len() {
                let u = queue[head];
                head += 1;
                for &v in &adj[offsets[u]..offsets[u + 1]] {
                    if !seen[v] {
                        seen[v] = true;
                        queue.push(v);
                    }
                }
            }
            comps.push(queue.clone());
        }
        comps
    }

    /// Number of connected components without any fixed node (each makes the
    /// force-density system singular).
    pub fn unsupported_components(&self) -> usize {
        let mut is_fixed = vec![false; self.num_nodes()];
        for &f in &self.fixed {
            is_fixed[f] = true;
        }
        self.components()
            .iter()
            .filter(|c| !c.iter().any(|&u| is_fixed[u]))
            .count()
    }

    /// Builds the [`Problem`] with a flat target (`z = -0.2` below the
    /// reference positions), mirroring `grid::make_grid_problem`.
    pub fn into_problem(self) -> Problem {
        let nn = self.num_nodes();
        let ne = self.num_edges();
        let mut rows = Vec::with_capacity(ne * 2);
        let mut cols = Vec::with_capacity(ne * 2);
        let mut vals = Vec::with_capacity(ne * 2);
        for (e, &(s, t)) in self.edges.iter().enumerate() {
            rows.extend([e, e]);
            cols.extend([s, t]);
            vals.extend([-1.0, 1.0]);
        }
        let incidence = SparseColMatOwned::from_coo(ne, nn, &rows, &cols, &vals).unwrap();
        let mut is_fixed = vec![false; nn];
        for &f in &self.fixed {
            is_fixed[f] = true;
        }
        let fixed_idx = self.fixed.clone();
        let free_idx: Vec<usize> = (0..nn).filter(|&i| !is_fixed[i]).collect();
        let nn_free = free_idx.len();
        let topology = NetworkTopology {
            free_incidence: incidence.extract_columns(&free_idx),
            fixed_incidence: incidence.extract_columns(&fixed_idx),
            incidence,
            num_edges: ne,
            num_nodes: nn,
            free_node_indices: free_idx,
            fixed_node_indices: fixed_idx,
        };
        let mut loads = vec![0.0; nn_free * 3];
        for i in 0..nn_free {
            loads[i * 3 + 2] = -1.0;
        }
        let mut fixed_pos = Array2::zeros((self.fixed.len(), 3));
        for (i, &node) in self.fixed.iter().enumerate() {
            for d in 0..3 {
                fixed_pos[[i, d]] = self.positions[node][d];
            }
        }
        let target_nodes = topology.free_node_indices.clone();
        let mut target = Array2::zeros((nn_free, 3));
        for (i, &node) in target_nodes.iter().enumerate() {
            target[[i, 0]] = self.positions[node][0];
            target[[i, 1]] = self.positions[node][1];
            target[[i, 2]] = self.positions[node][2] - 0.2;
        }
        Problem {
            topology,
            free_node_loads: Array2::from_shape_vec((nn_free, 3), loads).unwrap(),
            anchors: AnchorInfo::all_fixed(fixed_pos.clone()),
            fixed_node_positions: fixed_pos,
            objectives: vec![Box::new(TargetXYZ {
                weight: 1.0,
                node_indices: target_nodes,
                target,
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

    /// Builds the problem and replaces its target by the equilibrium shape of
    /// `q_star`, so the optimum (zero loss) lies inside the bounds.
    pub fn into_recoverable_problem(self, q_star: &[f64]) -> Problem {
        make_recoverable(self.into_problem(), q_star)
    }
}

/// Smooth, strictly interior force-density field used by the recoverable
/// fixtures: the same `1 + 0.5 sin + 0.25 cos` profile over the edge index as
/// `grid::make_recoverable_grid_problem`, values in `[0.25, 1.75]`.
pub fn smooth_q_star(ne: usize) -> Vec<f64> {
    (0..ne)
        .map(|k| {
            let t = k as f64 / ne as f64;
            let tau = std::f64::consts::TAU;
            1.0 + 0.5 * (tau * t).sin() + 0.25 * (2.0 * tau * t).cos()
        })
        .collect()
}

/// Replaces the objectives of `problem` by a `TargetXYZ` over the free nodes
/// whose target is the FDM equilibrium under `q_star`.
pub fn make_recoverable(mut problem: Problem, q_star: &[f64]) -> Problem {
    assert_eq!(q_star.len(), problem.topology.num_edges);
    for (k, &q) in q_star.iter().enumerate() {
        assert!(
            problem.bounds.lower[k] <= q && q <= problem.bounds.upper[k],
            "q_star[{k}] = {q} outside the bounds"
        );
    }
    let mut cache = FdmCache::new(&problem).unwrap();
    theseus::fdm::solve_fdm(&mut cache, q_star, &problem, &Array2::zeros((0, 3)), 0.0).unwrap();
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

/// Regular `side × side` grid (spacing 1 in the xy-plane, z = 0) with the four
/// corners fixed; the same topology and node numbering as `grid.rs`.
pub fn grid_network(side: usize) -> Network {
    assert!(side >= 2, "grid needs side >= 2");
    let positions = (0..side * side)
        .map(|i| [(i % side) as f64, (i / side) as f64, 0.0])
        .collect();
    let mut edges = Vec::with_capacity(2 * side * (side - 1));
    for row in 0..side {
        for col in 0..(side - 1) {
            edges.push((row * side + col, row * side + col + 1));
        }
    }
    for row in 0..(side - 1) {
        for col in 0..side {
            edges.push((row * side + col, (row + 1) * side + col));
        }
    }
    Network {
        positions,
        edges,
        fixed: vec![0, side - 1, side * (side - 1), side * side - 1],
    }
    .normalised()
}

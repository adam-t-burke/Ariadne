//! `graph::build`: adjacency round-trips edge lists, and the level-0 graph
//! operator equals the assembled `A = Cnᵀ diag(q) Cn` on random problems.

use ndarray::Array2;
use theseus::graph::build::{edge_endpoints, free_node_map};
use theseus::graph::{BlockVec, CsrAdjacency, LevelGraph};
use theseus::sparse::SparseColMatOwned;
use theseus::types::*;

struct Lcg(u64);

impl Lcg {
    fn next_u64(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0 >> 11
    }
    fn below(&mut self, n: usize) -> usize {
        (self.next_u64() % n as u64) as usize
    }
    fn unit(&mut self) -> f64 {
        (self.next_u64() % 1_000_000) as f64 / 1_000_000.0
    }
}

/// Random connected multigraph on `nn` nodes with `ne ≥ nn − 1` edges (a
/// random spanning tree plus random extra edges, including some parallel
/// edges) and `n_fixed` random fixed nodes.
fn random_problem(rng: &mut Lcg, nn: usize, ne: usize, n_fixed: usize) -> Problem {
    let mut edges: Vec<(usize, usize)> = Vec::with_capacity(ne);
    for v in 1..nn {
        let u = rng.below(v);
        if rng.below(2) == 0 {
            edges.push((u, v));
        } else {
            edges.push((v, u));
        }
    }
    while edges.len() < ne {
        let u = rng.below(nn);
        let mut v = rng.below(nn);
        if u == v {
            v = (v + 1) % nn;
        }
        edges.push((u, v));
    }
    // Shuffle so edge indices are not correlated with node order.
    for i in (1..edges.len()).rev() {
        let j = rng.below(i + 1);
        edges.swap(i, j);
    }

    let mut rows = Vec::with_capacity(ne * 2);
    let mut cols = Vec::with_capacity(ne * 2);
    let mut vals = Vec::with_capacity(ne * 2);
    for (k, &(s, e)) in edges.iter().enumerate() {
        rows.extend([k, k]);
        cols.extend([s, e]);
        vals.extend([-1.0, 1.0]);
    }
    let incidence = SparseColMatOwned::from_coo(ne, nn, &rows, &cols, &vals).unwrap();

    let mut fixed: Vec<usize> = Vec::new();
    while fixed.len() < n_fixed {
        let node = rng.below(nn);
        if !fixed.contains(&node) {
            fixed.push(node);
        }
    }
    fixed.sort_unstable();
    let free: Vec<usize> = (0..nn).filter(|i| !fixed.contains(i)).collect();
    let nn_free = free.len();

    let topology = NetworkTopology {
        free_incidence: incidence.extract_columns(&free),
        fixed_incidence: incidence.extract_columns(&fixed),
        incidence,
        num_edges: ne,
        num_nodes: nn,
        free_node_indices: free,
        fixed_node_indices: fixed.clone(),
    };
    let fixed_positions = Array2::from_shape_fn((fixed.len(), 3), |(_, _)| 2.0 * rng.unit() - 1.0);
    Problem {
        topology,
        free_node_loads: Array2::from_shape_fn((nn_free, 3), |(_, _)| rng.unit() - 0.5),
        fixed_node_positions: fixed_positions.clone(),
        anchors: AnchorInfo::all_fixed(fixed_positions),
        objectives: Vec::new(),
        bounds: Bounds::default_for(ne),
        solver: SolverOptions::default(),
        self_weight: None,
        pressure: None,
    }
}

#[test]
fn adjacency_round_trips_random_edge_lists() {
    let mut rng = Lcg(7);
    for trial in 0..20 {
        let nn = 2 + rng.below(40);
        let ne = (nn - 1) + rng.below(3 * nn);
        let problem = random_problem(&mut rng, nn, ne, 1);
        let (starts, ends) = edge_endpoints(&problem.topology);
        let adj = CsrAdjacency::from_topology(&problem.topology);
        assert_eq!(adj.num_nodes(), nn, "trial {trial}");
        assert_eq!(adj.num_entries(), 2 * ne, "trial {trial}");
        let recovered = adj.to_endpoints();
        for k in 0..ne {
            assert_eq!(
                recovered[k],
                (starts[k], ends[k]),
                "trial {trial}, edge {k}"
            );
        }
        // Incident edges ascending, `other` consistent with the edge list.
        for u in 0..nn {
            let mut prev = None;
            for (e, s, other) in adj.incident(u) {
                assert!(prev.is_none_or(|p| p <= e), "node {u} not ascending");
                prev = Some(e);
                let (start, end) = recovered[e as usize];
                if s < 0 {
                    assert_eq!(start, u);
                    assert_eq!(other as usize, end);
                } else {
                    assert_eq!(end, u);
                    assert_eq!(other as usize, start);
                }
            }
        }
        assert_eq!(
            CsrAdjacency::from_endpoints(nn, &starts, &ends),
            adj,
            "from_endpoints and from_topology agree"
        );
    }
}

#[test]
fn level0_operator_matches_assembled_matrix() {
    let mut rng = Lcg(11);
    for trial in 0..25 {
        let nn = 3 + rng.below(30);
        let ne = (nn - 1) + rng.below(2 * nn);
        let n_fixed = 1 + rng.below((nn / 3).max(1));
        let problem = random_problem(&mut rng, nn, ne, n_fixed);
        let nn_free = problem.topology.free_node_indices.len();
        let q: Vec<f64> = (0..ne).map(|_| 0.1 + 10.0 * rng.unit()).collect();

        let mut cache = FdmCache::new(&problem).unwrap();
        cache.q.copy_from_slice(&q);
        theseus::fdm::assemble_a(&mut cache);
        let a = &cache.a_matrix;
        assert_eq!(a.nrows, nn_free);

        let level = LevelGraph::level0(&problem.topology, &free_node_map(&problem.topology), &q);
        assert_eq!(level.n, nn_free);
        assert_eq!(level.adjacency, {
            // The cache's all-node adjacency restricted to free–free edges
            // must be the level adjacency up to renumbering.
            let map = free_node_map(&problem.topology);
            let (starts, ends) = edge_endpoints(&problem.topology);
            let mut ls = Vec::new();
            let mut le = Vec::new();
            for k in 0..ne {
                if let (Some(s), Some(e)) = (map[starts[k]], map[ends[k]]) {
                    ls.push(s);
                    le.push(e);
                }
            }
            CsrAdjacency::from_endpoints(nn_free, &ls, &le)
        });

        let x = BlockVec::<3>::from_vec((0..nn_free * 3).map(|_| 2.0 * rng.unit() - 1.0).collect());
        let mut y = BlockVec::<3>::zeros(nn_free);
        level.apply(&x, &mut y);

        // Reference: y_ref = A x through the CSC matrix (both triangles stored).
        let mut y_ref = vec![0.0; nn_free * 3];
        for col in 0..a.ncols {
            for nz in a.col_ptrs[col] as usize..a.col_ptrs[col + 1] as usize {
                let row = a.row_indices[nz] as usize;
                let v = a.values[nz];
                for d in 0..3 {
                    y_ref[row * 3 + d] += v * x.data[col * 3 + d];
                }
            }
        }
        let scale = y_ref.iter().fold(0.0f64, |m, v| m.max(v.abs())).max(1e-300);
        for i in 0..nn_free * 3 {
            let err = (y.data[i] - y_ref[i]).abs();
            assert!(
                err <= 1e-12 * scale,
                "trial {trial}: entry {i}: level {} vs matrix {} (rel err {:e})",
                y.data[i],
                y_ref[i],
                err / scale
            );
        }

        // The cache's adjacency is the all-node adjacency.
        assert_eq!(
            cache.adjacency,
            CsrAdjacency::from_topology(&problem.topology)
        );
    }
}

//! Plain graph data shared by the iterative solvers, the parallel O(ne)
//! loops and the AMG hierarchy: node-centred CSR adjacency
//! ([`CsrAdjacency`]), one multigrid level as a weighted graph
//! ([`LevelGraph`]) and blocks of right-hand sides ([`BlockVec`]).
//!
//! These are data types only. Construction from `NetworkTopology`
//! (boundary edges, anchor weights, free-node renumbering) lives in
//! `graph/build.rs`; aggregation and coarse-level construction live in
//! `amg/`.
//!
//! Conventions: node and edge indices are `u32` (the 10M-edge target has
//! ~5M nodes and ~20M adjacency entries); an edge `e = (start, end)` has
//! sign `+1` at its `end` node and `−1` at its `start` node, matching the
//! `±1` entries of the incidence matrix (`−1` at start, `+1` at end).

pub mod build;

pub use build::Level0Map;

/// Block of `K` right-hand sides / solutions, row-major: `v[node * K + k]`.
#[derive(Debug, Clone, PartialEq)]
pub struct BlockVec<const K: usize> {
    /// `n * K` values.
    pub data: Vec<f64>,
    /// Number of nodes (rows).
    pub n: usize,
}

impl<const K: usize> BlockVec<K> {
    /// `n` rows of zeros.
    pub fn zeros(n: usize) -> Self {
        Self {
            data: vec![0.0; n * K],
            n,
        }
    }

    /// Wrap an existing row-major `n * K` vector.
    ///
    /// # Panics
    /// If `data.len()` is not a multiple of `K`.
    pub fn from_vec(data: Vec<f64>) -> Self {
        assert_eq!(data.len() % K, 0, "BlockVec: length not a multiple of K");
        let n = data.len() / K;
        Self { data, n }
    }

    /// Number of rows.
    pub fn len(&self) -> usize {
        self.n
    }

    /// `true` when there are no rows.
    pub fn is_empty(&self) -> bool {
        self.n == 0
    }

    /// The `K` values of row `i`.
    pub fn row(&self, i: usize) -> &[f64; K] {
        self.data[i * K..(i + 1) * K]
            .try_into()
            .expect("row slice has K elements")
    }

    /// Mutable `K` values of row `i`.
    pub fn row_mut(&mut self, i: usize) -> &mut [f64; K] {
        (&mut self.data[i * K..(i + 1) * K])
            .try_into()
            .expect("row slice has K elements")
    }

    /// Set every value to zero.
    pub fn fill_zero(&mut self) {
        self.data.fill(0.0);
    }

    /// Row-major view.
    pub fn as_slice(&self) -> &[f64] {
        &self.data
    }

    /// Mutable row-major view.
    pub fn as_mut_slice(&mut self) -> &mut [f64] {
        &mut self.data
    }
}

/// Node-centred adjacency of an undirected edge list.
///
/// For node `i`, incident entries are `offsets[i]..offsets[i+1]`; entry `j`
/// refers to edge `edges[j]`, with `sign[j] = +1` if `i` is the edge's end
/// node and `−1` if it is the start node, and `other[j]` the node at the far
/// end. Every edge appears exactly twice (once per endpoint). Built once per
/// topology; shared by the AMG levels, geometry, and gradient loops.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct CsrAdjacency {
    /// `num_nodes + 1` offsets into `edges` / `sign` / `other`.
    pub offsets: Vec<u32>,
    /// Incident edge index per entry.
    pub edges: Vec<u32>,
    /// `+1` (node is the edge end) or `−1` (node is the edge start).
    pub sign: Vec<i8>,
    /// Node at the other end of the edge.
    pub other: Vec<u32>,
}

impl CsrAdjacency {
    /// Number of nodes.
    pub fn num_nodes(&self) -> usize {
        self.offsets.len().saturating_sub(1)
    }

    /// Number of (node, edge) incidences, i.e. `2 × num_edges` for a graph
    /// with both endpoints present.
    pub fn num_entries(&self) -> usize {
        self.edges.len()
    }

    /// Number of incident edges of `node`.
    pub fn degree(&self, node: usize) -> usize {
        (self.offsets[node + 1] - self.offsets[node]) as usize
    }

    /// Entry range of `node` in `edges` / `sign` / `other`.
    pub fn range(&self, node: usize) -> std::ops::Range<usize> {
        self.offsets[node] as usize..self.offsets[node + 1] as usize
    }

    /// Incident edges of `node` as `(edge, sign, other_node)`.
    pub fn incident(&self, node: usize) -> impl Iterator<Item = (u32, i8, u32)> + '_ {
        let range = self.range(node);
        self.edges[range.clone()]
            .iter()
            .zip(&self.sign[range.clone()])
            .zip(&self.other[range])
            .map(|((&edge, &sign), &other)| (edge, sign, other))
    }

    /// Neighbouring nodes of `node` (with repetition for multi-edges).
    pub fn neighbours(&self, node: usize) -> &[u32] {
        &self.other[self.range(node)]
    }
}

/// One multigrid level as a weighted graph on free nodes with anchor weights.
///
/// The level operator is
/// `(A x)_u = anchor_u x_u + Σ_{e ∋ u} weight_e (x_u − x_other)`;
/// level 0 has `weight_e = q_e` over the free–free edges and `anchor_u = Σ q`
/// over the edges from `u` to fixed nodes. Coarser levels are quotient
/// graphs whose weights sum the fine edges they represent, so every level
/// uses the same edge kernel.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct LevelGraph {
    /// Number of nodes on this level.
    pub n: usize,
    /// Adjacency over the level's edges (both endpoints in `0..n`).
    pub adjacency: CsrAdjacency,
    /// Per edge: `Σ q` of the fine edges it represents.
    pub weight: Vec<f64>,
    /// Per node: `Σ q` of the fine edges to fixed nodes.
    pub anchor: Vec<f64>,
    /// Fine node → coarse node of the next level (empty on the coarsest).
    pub aggregate_of: Vec<u32>,
}

impl LevelGraph {
    /// Number of edges on this level.
    pub fn num_edges(&self) -> usize {
        self.weight.len()
    }

    /// `true` when no coarser level follows.
    pub fn is_coarsest(&self) -> bool {
        self.aggregate_of.is_empty()
    }

    /// Jacobi diagonal of node `u`: `anchor_u + Σ_{e ∋ u} weight_e`.
    pub fn diagonal(&self, u: usize) -> f64 {
        self.anchor[u]
            + self
                .adjacency
                .incident(u)
                .map(|(e, _, _)| self.weight[e as usize])
                .sum::<f64>()
    }

    /// Reference operator application `y = A x` on `K` columns, serial;
    /// used to validate backend kernels.
    pub fn apply<const K: usize>(&self, x: &BlockVec<K>, y: &mut BlockVec<K>) {
        debug_assert_eq!(x.n, self.n);
        debug_assert_eq!(y.n, self.n);
        for u in 0..self.n {
            let xu = *x.row(u);
            let mut acc = [0.0; K];
            for k in 0..K {
                acc[k] = self.anchor[u] * xu[k];
            }
            for (e, _, v) in self.adjacency.incident(u) {
                let w = self.weight[e as usize];
                let xv = x.row(v as usize);
                for k in 0..K {
                    acc[k] += w * (xu[k] - xv[k]);
                }
            }
            *y.row_mut(u) = acc;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Path 0 — 1 — 2 with edges e0 = (0, 1), e1 = (1, 2).
    fn path3() -> CsrAdjacency {
        CsrAdjacency {
            offsets: vec![0, 1, 3, 4],
            edges: vec![0, 0, 1, 1],
            sign: vec![-1, 1, -1, 1],
            other: vec![1, 0, 2, 1],
        }
    }

    #[test]
    fn adjacency_helpers() {
        let adj = path3();
        assert_eq!(adj.num_nodes(), 3);
        assert_eq!(adj.num_entries(), 4);
        assert_eq!(adj.degree(1), 2);
        assert_eq!(
            adj.incident(1).collect::<Vec<_>>(),
            vec![(0, 1, 0), (1, -1, 2)]
        );
        assert_eq!(adj.neighbours(2), &[1]);
        assert_eq!(CsrAdjacency::default().num_nodes(), 0);
    }

    #[test]
    fn level_graph_operator_matches_laplacian_plus_anchor() {
        let level = LevelGraph {
            n: 3,
            adjacency: path3(),
            weight: vec![2.0, 3.0],
            anchor: vec![1.0, 0.0, 4.0],
            aggregate_of: Vec::new(),
        };
        assert!(level.is_coarsest());
        assert_eq!(level.num_edges(), 2);
        assert_eq!(level.diagonal(0), 3.0);
        assert_eq!(level.diagonal(1), 5.0);
        assert_eq!(level.diagonal(2), 7.0);
        let x = BlockVec::<3>::from_vec(vec![1.0, 0.0, 1.0, 0.0, 1.0, 2.0, 0.0, 0.0, 3.0]);
        let mut y = BlockVec::<3>::zeros(3);
        level.apply(&x, &mut y);
        // A = [[3, -2, 0], [-2, 5, -3], [0, -3, 7]]
        assert_eq!(y.row(0), &[3.0, -2.0, -1.0]);
        assert_eq!(y.row(1), &[-2.0, 5.0, -1.0]);
        assert_eq!(y.row(2), &[0.0, -3.0, 15.0]);
    }

    #[test]
    fn block_vec_rows() {
        let mut v = BlockVec::<2>::zeros(2);
        assert_eq!(v.len(), 2);
        assert!(!v.is_empty());
        *v.row_mut(1) = [5.0, 6.0];
        assert_eq!(v.as_slice(), &[0.0, 0.0, 5.0, 6.0]);
        v.fill_zero();
        assert_eq!(v.as_slice(), &[0.0; 4]);
    }
}

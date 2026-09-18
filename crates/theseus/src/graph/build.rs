//! Construction of the shared graph data from a [`NetworkTopology`]:
//! node-centred [`CsrAdjacency`] over all nodes, and the level-0
//! [`LevelGraph`] on the free nodes whose `anchor` weights collect the edges
//! to fixed nodes.
//!
//! All builders are sequential counting sorts, so the incident edges of a
//! node are always listed in **ascending edge index** order. The parallel
//! node-gather loops in `fdm.rs` / `gradients.rs` rely on that order to
//! reproduce the sequential edge-scatter sums bit for bit.

use super::{CsrAdjacency, LevelGraph};
use crate::types::NetworkTopology;

/// Start and end node of every edge, read from the `±1` entries of the
/// full incidence matrix (`−1` at the start node, `+1` at the end node).
pub fn edge_endpoints(topo: &NetworkTopology) -> (Vec<usize>, Vec<usize>) {
    let ne = topo.num_edges;
    let mut starts = vec![0usize; ne];
    let mut ends = vec![0usize; ne];
    let inc = &topo.incidence;
    for col in 0..topo.num_nodes {
        let range = inc.col_ptrs[col] as usize..inc.col_ptrs[col + 1] as usize;
        for idx in range {
            let row = inc.row_indices[idx] as usize;
            if inc.values[idx] < 0.0 {
                starts[row] = col;
            } else {
                ends[row] = col;
            }
        }
    }
    (starts, ends)
}

/// Global node → free row (`None` for fixed nodes), as used by `FdmCache`.
pub fn free_node_map(topo: &NetworkTopology) -> Vec<Option<usize>> {
    let mut map = vec![None; topo.num_nodes];
    for (i, &node) in topo.free_node_indices.iter().enumerate() {
        map[node] = Some(i);
    }
    map
}

impl CsrAdjacency {
    /// Adjacency of `num_nodes` nodes from parallel start/end lists.
    ///
    /// Incident edges of every node are listed in ascending edge order; a
    /// self-loop (`start == end`) appears twice at its node, once per sign.
    ///
    /// # Panics
    /// If the lists differ in length or an endpoint is `≥ num_nodes`.
    pub fn from_endpoints(num_nodes: usize, starts: &[usize], ends: &[usize]) -> Self {
        assert_eq!(starts.len(), ends.len(), "edge lists differ in length");
        let ne = starts.len();
        let mut offsets = vec![0u32; num_nodes + 1];
        for k in 0..ne {
            assert!(
                starts[k] < num_nodes && ends[k] < num_nodes,
                "edge {k} out of range"
            );
            offsets[starts[k] + 1] += 1;
            offsets[ends[k] + 1] += 1;
        }
        for i in 0..num_nodes {
            offsets[i + 1] += offsets[i];
        }
        let total = offsets[num_nodes] as usize;
        let mut fill: Vec<u32> = offsets[..num_nodes].to_vec();
        let mut edges = vec![0u32; total];
        let mut sign = vec![0i8; total];
        let mut other = vec![0u32; total];
        for k in 0..ne {
            let (s, e) = (starts[k], ends[k]);
            let slot = fill[s] as usize;
            fill[s] += 1;
            edges[slot] = k as u32;
            sign[slot] = -1;
            other[slot] = e as u32;
            let slot = fill[e] as usize;
            fill[e] += 1;
            edges[slot] = k as u32;
            sign[slot] = 1;
            other[slot] = s as u32;
        }
        Self {
            offsets,
            edges,
            sign,
            other,
        }
    }

    /// Adjacency over **all** nodes of a topology (free and fixed), edge
    /// indices as in the incidence matrix rows.
    pub fn from_topology(topo: &NetworkTopology) -> Self {
        let (starts, ends) = edge_endpoints(topo);
        Self::from_endpoints(topo.num_nodes, &starts, &ends)
    }

    /// Reconstruct the `(start, end)` list the adjacency was built from.
    pub fn to_endpoints(&self) -> Vec<(usize, usize)> {
        let ne = self.num_entries() / 2;
        let mut edges = vec![(usize::MAX, usize::MAX); ne];
        for u in 0..self.num_nodes() {
            for (e, s, _) in self.incident(u) {
                if s < 0 {
                    edges[e as usize].0 = u;
                } else {
                    edges[e as usize].1 = u;
                }
            }
        }
        edges
    }
}

/// Fixed connectivity of the level-0 graph: which fine edges the level's
/// edges are, and which fine edges feed each free node's `anchor`.
///
/// Built once per topology; [`Level0Map::update_weights`] refreshes a
/// [`LevelGraph`]'s `weight` / `anchor` arrays from a new `q` without
/// touching the adjacency.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Level0Map {
    /// Level edge → fine edge (both endpoints free), ascending.
    pub fine_edge: Vec<u32>,
    /// `n_free + 1` offsets into `anchor_edges`.
    pub anchor_offsets: Vec<u32>,
    /// Fine edges with exactly one fixed endpoint, grouped by their free
    /// node, ascending within a node.
    pub anchor_edges: Vec<u32>,
}

impl Level0Map {
    /// Build the level-0 connectivity from the all-node adjacency and the
    /// free-node map. Returns the [`LevelGraph`] with zero weights and the
    /// map that fills them.
    ///
    /// Edges with both endpoints fixed are dropped (they influence neither
    /// `A` nor `b`).
    pub fn new(
        adjacency: &CsrAdjacency,
        node_to_free: &[Option<usize>],
        free_nodes: &[usize],
    ) -> (LevelGraph, Self) {
        let n = free_nodes.len();
        let edges = adjacency.to_endpoints();
        let mut fine_edge = Vec::new();
        let mut level_starts = Vec::new();
        let mut level_ends = Vec::new();
        for (k, &(s, e)) in edges.iter().enumerate() {
            if let (Some(sf), Some(ef)) = (node_to_free[s], node_to_free[e]) {
                fine_edge.push(k as u32);
                level_starts.push(sf);
                level_ends.push(ef);
            }
        }
        let level_adjacency = CsrAdjacency::from_endpoints(n, &level_starts, &level_ends);

        let mut anchor_offsets = vec![0u32; n + 1];
        let mut anchor_edges = Vec::new();
        for (f, &node) in free_nodes.iter().enumerate() {
            for (k, _, other) in adjacency.incident(node) {
                if node_to_free[other as usize].is_none() {
                    anchor_edges.push(k);
                }
            }
            anchor_offsets[f + 1] = anchor_edges.len() as u32;
        }

        let level = LevelGraph {
            n,
            adjacency: level_adjacency,
            weight: vec![0.0; fine_edge.len()],
            anchor: vec![0.0; n],
            aggregate_of: Vec::new(),
        };
        let map = Self {
            fine_edge,
            anchor_offsets,
            anchor_edges,
        };
        (level, map)
    }

    /// Fill `level.weight[e] = q[fine_edge[e]]` and
    /// `level.anchor[u] = Σ q` over `u`'s boundary edges (ascending edge
    /// order, sequential sum).
    pub fn update_weights(&self, q: &[f64], level: &mut LevelGraph) {
        debug_assert_eq!(level.weight.len(), self.fine_edge.len());
        debug_assert_eq!(level.anchor.len() + 1, self.anchor_offsets.len());
        for (w, &k) in level.weight.iter_mut().zip(&self.fine_edge) {
            *w = q[k as usize];
        }
        for (u, a) in level.anchor.iter_mut().enumerate() {
            let range = self.anchor_offsets[u] as usize..self.anchor_offsets[u + 1] as usize;
            let mut acc = 0.0;
            for &k in &self.anchor_edges[range] {
                acc += q[k as usize];
            }
            *a = acc;
        }
    }
}

impl LevelGraph {
    /// Level-0 graph of a topology for force densities `q`: nodes are the
    /// free nodes in `topo.free_node_indices` order, edges are the free–free
    /// edges with `weight = q`, and `anchor` sums `q` over each node's edges
    /// to fixed nodes. `node_to_free` is the map from [`free_node_map`] (or
    /// `FdmCache::node_to_free_idx`).
    pub fn level0(topo: &NetworkTopology, node_to_free: &[Option<usize>], q: &[f64]) -> Self {
        let adjacency = CsrAdjacency::from_topology(topo);
        let (mut level, map) = Level0Map::new(&adjacency, node_to_free, &topo.free_node_indices);
        map.update_weights(q, &mut level);
        level
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn from_endpoints_lists_edges_ascending_with_signs() {
        // 0 -e0-> 1 -e1-> 2, plus e2: 2 -> 0 and a self-loop e3 at node 1.
        let adj = CsrAdjacency::from_endpoints(3, &[0, 1, 2, 1], &[1, 2, 0, 1]);
        assert_eq!(adj.num_nodes(), 3);
        assert_eq!(adj.num_entries(), 8);
        assert_eq!(
            adj.incident(0).collect::<Vec<_>>(),
            vec![(0, -1, 1), (2, 1, 2)]
        );
        assert_eq!(
            adj.incident(1).collect::<Vec<_>>(),
            vec![(0, 1, 0), (1, -1, 2), (3, -1, 1), (3, 1, 1)]
        );
        assert_eq!(
            adj.incident(2).collect::<Vec<_>>(),
            vec![(1, 1, 1), (2, -1, 0)]
        );
        assert_eq!(adj.to_endpoints(), vec![(0, 1), (1, 2), (2, 0), (1, 1)]);
    }

    #[test]
    fn level0_map_anchor_and_weights() {
        // Path 0 - 1 - 2 - 3 with 0 and 3 fixed; free rows: 1 -> 0, 2 -> 1.
        let adj = CsrAdjacency::from_endpoints(4, &[0, 1, 2], &[1, 2, 3]);
        let node_to_free = vec![None, Some(0), Some(1), None];
        let (mut level, map) = Level0Map::new(&adj, &node_to_free, &[1, 2]);
        assert_eq!(map.fine_edge, vec![1]);
        assert_eq!(map.anchor_offsets, vec![0, 1, 2]);
        assert_eq!(map.anchor_edges, vec![0, 2]);
        map.update_weights(&[2.0, 3.0, 5.0], &mut level);
        assert_eq!(level.weight, vec![3.0]);
        assert_eq!(level.anchor, vec![2.0, 5.0]);
        assert_eq!(level.diagonal(0), 5.0);
        assert_eq!(level.diagonal(1), 8.0);
    }
}

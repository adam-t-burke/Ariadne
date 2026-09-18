//! Pairwise heavy-edge aggregation (program plan §3).
//!
//! One matching pass visits the nodes in decreasing degree order (ties by
//! lowest index) and pairs each unmatched node with its strongest unmatched
//! neighbour, strength `s_uv = w_uv / sqrt(d_u d_v)`; nodes left over join
//! the strongest neighbour's aggregate if it has fewer than three members,
//! otherwise they form singletons. `passes` passes are composed by running
//! the next pass on the quotient graph of the previous one, so three passes
//! give aggregates of about eight nodes.
//!
//! Everything here is sequential and index-ordered, so the aggregation of a
//! given graph is identical on every thread count and every platform.

use super::hierarchy::LevelMatrix;
use crate::graph::LevelGraph;

const NONE: u32 = u32::MAX;

/// Undirected weighted graph in node-centred CSR form, as the matching sees
/// it: entry `j` of node `u` points to `other[j]` with weight `weight[j]`,
/// and `diag[u]` is the row sum the strength is normalised by. Level 0 is
/// the force-density graph itself; coarser levels use `|a_uv|` of the
/// Galerkin operator, and intermediate matching passes use quotient graphs.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct StrengthGraph {
    pub n: usize,
    pub offsets: Vec<u32>,
    pub other: Vec<u32>,
    pub weight: Vec<f64>,
    pub diag: Vec<f64>,
}

impl StrengthGraph {
    /// The level-0 graph: entry weights are the edge weights, `diag_u =
    /// anchor_u + Σ w`. Self-loops (an edge from a node to itself) are kept
    /// as entries of weight zero so they never attract a match.
    pub fn from_level_graph(g: &LevelGraph) -> Self {
        let adj = &g.adjacency;
        let mut weight = Vec::with_capacity(adj.num_entries());
        let mut diag = vec![0.0; g.n];
        for u in 0..g.n {
            let mut d = g.anchor[u];
            for (e, _, v) in adj.incident(u) {
                let w = g.weight[e as usize];
                d += w;
                weight.push(if v as usize == u { 0.0 } else { w });
            }
            diag[u] = d;
        }
        Self {
            n: g.n,
            offsets: adj.offsets.clone(),
            other: adj.other.clone(),
            weight,
            diag,
        }
    }

    /// A general sparse level: off-diagonal entries with weight `|a_uv|`,
    /// `diag_u = a_uu`.
    pub fn from_matrix(m: &LevelMatrix) -> Self {
        assert_eq!(m.n, m.ncols, "strength graph needs a square matrix");
        let mut offsets = vec![0u32; m.n + 1];
        let mut other = Vec::with_capacity(m.nnz());
        let mut weight = Vec::with_capacity(m.nnz());
        for u in 0..m.n {
            let (cols, vals) = m.row(u);
            for (&c, &v) in cols.iter().zip(vals) {
                if c as usize != u {
                    other.push(c);
                    weight.push(v.abs());
                }
            }
            offsets[u + 1] = other.len() as u32;
        }
        Self {
            n: m.n,
            offsets,
            other,
            weight,
            diag: m.diag.clone(),
        }
    }

    /// Number of entries (twice the number of undirected edges).
    pub fn num_entries(&self) -> usize {
        self.other.len()
    }

    /// Number of incident entries of `u`.
    pub fn degree(&self, u: usize) -> usize {
        (self.offsets[u + 1] - self.offsets[u]) as usize
    }

    fn range(&self, u: usize) -> std::ops::Range<usize> {
        self.offsets[u] as usize..self.offsets[u + 1] as usize
    }

    /// Quotient graph under `aggregate_of` (values `< n_coarse`): inter-
    /// aggregate entries are summed per coarse pair, and the coarse diagonal
    /// is the anchor mass of the members (`Σ max(d_u − Σ_v w_uv, 0)`) plus
    /// the coarse entry weights, i.e. the diagonal of `P₀ᵀ A P₀` for a graph
    /// Laplacian.
    pub fn quotient(&self, aggregate_of: &[u32], n_coarse: usize) -> Self {
        debug_assert_eq!(aggregate_of.len(), self.n);
        let mut pairs: Vec<(u32, u32, f64)> = Vec::with_capacity(self.num_entries() / 2);
        let mut anchor = vec![0.0; n_coarse];
        for u in 0..self.n {
            let a = aggregate_of[u];
            let mut row_sum = 0.0;
            for j in self.range(u) {
                let v = self.other[j] as usize;
                let w = self.weight[j];
                row_sum += w;
                let b = aggregate_of[v];
                if a < b {
                    pairs.push((a, b, w));
                }
            }
            anchor[a as usize] += (self.diag[u] - row_sum).max(0.0);
        }
        pairs.sort_unstable_by(|x, y| (x.0, x.1).cmp(&(y.0, y.1)));
        // Merge duplicates, then emit both directions.
        let mut edges: Vec<(u32, u32, f64)> = Vec::with_capacity(pairs.len());
        for &(a, b, w) in &pairs {
            match edges.last_mut() {
                Some(last) if last.0 == a && last.1 == b => last.2 += w,
                _ => edges.push((a, b, w)),
            }
        }
        let mut offsets = vec![0u32; n_coarse + 1];
        for &(a, b, _) in &edges {
            offsets[a as usize + 1] += 1;
            offsets[b as usize + 1] += 1;
        }
        for c in 0..n_coarse {
            offsets[c + 1] += offsets[c];
        }
        let total = offsets[n_coarse] as usize;
        let mut fill: Vec<u32> = offsets[..n_coarse].to_vec();
        let mut other = vec![0u32; total];
        let mut weight = vec![0.0; total];
        let mut diag = anchor;
        for &(a, b, w) in &edges {
            let slot = fill[a as usize] as usize;
            fill[a as usize] += 1;
            other[slot] = b;
            weight[slot] = w;
            let slot = fill[b as usize] as usize;
            fill[b as usize] += 1;
            other[slot] = a;
            weight[slot] = w;
            diag[a as usize] += w;
            diag[b as usize] += w;
        }
        Self {
            n: n_coarse,
            offsets,
            other,
            weight,
            diag,
        }
    }
}

/// Strongest neighbour of `u` accepted by `accept`, by `w_uv / sqrt(d_u d_v)`
/// with the lowest index winning ties. Zero-weight entries never match.
fn strongest(g: &StrengthGraph, u: usize, accept: impl Fn(usize) -> bool) -> Option<usize> {
    let mut best: Option<usize> = None;
    let mut best_s = 0.0;
    let du = g.diag[u];
    for j in g.range(u) {
        let v = g.other[j] as usize;
        if v == u || !accept(v) {
            continue;
        }
        let denominator = (du * g.diag[v]).sqrt();
        let s = if denominator > 0.0 {
            g.weight[j] / denominator
        } else {
            0.0
        };
        if s > best_s || (s == best_s && s > 0.0 && best.is_some_and(|b| v < b)) {
            best = Some(v);
            best_s = s;
        }
    }
    best
}

/// One pass of pairwise matching. Returns `(aggregate_of, n_coarse)` with
/// aggregate ids numbered in the order they are created.
pub fn pairwise_match(g: &StrengthGraph) -> (Vec<u32>, usize) {
    let n = g.n;
    let mut order: Vec<u32> = (0..n as u32).collect();
    order.sort_by_key(|&u| (std::cmp::Reverse(g.degree(u as usize)), u));
    let mut agg = vec![NONE; n];
    let mut members: Vec<u8> = Vec::with_capacity(n / 2 + 1);

    for &u in &order {
        let u = u as usize;
        if agg[u] != NONE {
            continue;
        }
        if let Some(v) = strongest(g, u, |v| agg[v] == NONE) {
            let id = members.len() as u32;
            agg[u] = id;
            agg[v] = id;
            members.push(2);
        }
    }
    for &u in &order {
        let u = u as usize;
        if agg[u] != NONE {
            continue;
        }
        match strongest(g, u, |v| agg[v] != NONE && members[agg[v] as usize] < 3) {
            Some(v) => {
                agg[u] = agg[v];
                members[agg[v] as usize] += 1;
            }
            None => {
                agg[u] = members.len() as u32;
                members.push(1);
            }
        }
    }
    (agg, members.len())
}

/// Renumber aggregates in order of their lowest fine member, so the coarse
/// numbering inherits the fine one's locality. `pairwise_match` numbers
/// aggregates in visiting (degree-sorted) order, which on meshes with mixed
/// degrees scatters neighbouring aggregates over the whole index range and
/// makes every coarse-level gather a cache miss (the irregular mesh's
/// Galerkin refill ran 2× slower than the grid's for less work).
pub fn relabel_by_first_member(aggregate_of: &mut [u32], n_coarse: usize) {
    let mut new_id = vec![NONE; n_coarse];
    let mut next = 0u32;
    for a in aggregate_of.iter_mut() {
        let slot = &mut new_id[*a as usize];
        if *slot == NONE {
            *slot = next;
            next += 1;
        }
        *a = *slot;
    }
    debug_assert_eq!(next as usize, n_coarse, "every aggregate has a member");
}

/// `passes` matching passes composed into one fine → coarse map, numbered
/// by lowest fine member ([`relabel_by_first_member`]).
pub fn aggregate(g: &StrengthGraph, passes: usize) -> (Vec<u32>, usize) {
    let (mut agg, mut nc) = pairwise_match(g);
    for _ in 1..passes.max(1) {
        relabel_by_first_member(&mut agg, nc);
        let coarse = g.quotient(&agg, nc);
        let (agg2, nc2) = pairwise_match(&coarse);
        for a in agg.iter_mut() {
            *a = agg2[*a as usize];
        }
        nc = nc2;
    }
    relabel_by_first_member(&mut agg, nc);
    (agg, nc)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::CsrAdjacency;

    /// Path 0 — 1 — 2 — 3 with weights 1, 5, 1 and anchors at the ends.
    fn path4() -> LevelGraph {
        LevelGraph {
            n: 4,
            adjacency: CsrAdjacency::from_endpoints(4, &[0, 1, 2], &[1, 2, 3]),
            weight: vec![1.0, 5.0, 1.0],
            anchor: vec![1.0, 0.0, 0.0, 1.0],
            aggregate_of: Vec::new(),
        }
    }

    #[test]
    fn strength_graph_from_level_graph() {
        let g = StrengthGraph::from_level_graph(&path4());
        assert_eq!(g.n, 4);
        assert_eq!(g.diag, vec![2.0, 6.0, 6.0, 2.0]);
        assert_eq!(g.degree(1), 2);
        assert_eq!(g.num_entries(), 6);
    }

    #[test]
    fn matching_pairs_strongest_edge_first() {
        let g = StrengthGraph::from_level_graph(&path4());
        let (agg, nc) = pairwise_match(&g);
        // Nodes 1 and 2 (degree 2) are visited first and pair through the
        // weight-5 edge; 0 and 3 are left over and join their only
        // neighbour's pair (2 members → allowed), giving one aggregate of
        // 3 plus a singleton... unless the pair is already 3.
        assert_eq!(agg[1], agg[2]);
        assert_eq!(nc, 2);
        assert!(agg.iter().all(|&a| (a as usize) < nc));
        let sizes = {
            let mut s = vec![0; nc];
            for &a in &agg {
                s[a as usize] += 1;
            }
            s
        };
        assert_eq!(sizes.iter().sum::<usize>(), 4);
        assert!(sizes.iter().all(|&s| s <= 3));
    }

    #[test]
    fn passes_compose_and_quotient_preserves_mass() {
        // 8 × 8 grid, unit weights, anchors on the boundary.
        let side = 8;
        let mut starts = Vec::new();
        let mut ends = Vec::new();
        for i in 0..side {
            for j in 0..side {
                if j + 1 < side {
                    starts.push(i * side + j);
                    ends.push(i * side + j + 1);
                }
                if i + 1 < side {
                    starts.push(i * side + j);
                    ends.push((i + 1) * side + j);
                }
            }
        }
        let n = side * side;
        let g = LevelGraph {
            n,
            adjacency: CsrAdjacency::from_endpoints(n, &starts, &ends),
            weight: vec![1.0; starts.len()],
            anchor: (0..n)
                .map(|u| {
                    let (i, j) = (u / side, u % side);
                    if i == 0 || j == 0 || i + 1 == side || j + 1 == side {
                        1.0
                    } else {
                        0.0
                    }
                })
                .collect(),
            aggregate_of: Vec::new(),
        };
        let sg = StrengthGraph::from_level_graph(&g);
        let (agg1, nc1) = aggregate(&sg, 1);
        let (agg3, nc3) = aggregate(&sg, 3);
        assert!(nc1 <= n / 2 + 1 && nc1 >= n / 3);
        assert!(
            nc3 <= n / 6,
            "three passes should give aggregates of ~8, got {nc3}"
        );
        assert!(agg3.iter().all(|&a| (a as usize) < nc3));
        // Quotient keeps the total anchor mass and the total row-sum mass.
        let q = sg.quotient(&agg1, nc1);
        let mass = |g: &StrengthGraph| g.diag.iter().sum::<f64>();
        let entry_mass = |g: &StrengthGraph| g.weight.iter().sum::<f64>();
        let anchor_total: f64 = g.anchor.iter().sum();
        assert!((mass(&q) - entry_mass(&q) - anchor_total).abs() < 1e-12);
        // Determinism: same input, same output.
        assert_eq!(aggregate(&sg, 3), (agg3, nc3));
    }
}

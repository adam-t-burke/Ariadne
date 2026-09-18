//! Sparse level matrices and the building blocks of the smoothed-aggregation
//! hierarchy (program plan §3).
//!
//! * [`LevelMatrix`] — CSR with ascending column indices within a row (the
//!   convention the GPU kernels assume) and the diagonal kept separately for
//!   Jacobi scaling. Coarse operators `A_l`, `l ≥ 1`, and the transfers `P`,
//!   `Pᵀ` are all `LevelMatrix` values.
//! * [`LevelRef`] — one level operator seen as rows of `(column, value)`
//!   pairs, whether it is the matrix-free level-0 graph or a CSR level, so
//!   the prolongator and the triple product are written once.
//! * [`smoothed_prolongator`] — `P = (I − ω D⁻¹ A) P₀`.
//! * [`galerkin_pattern`] / [`galerkin_numeric`] — the symbolic pattern of
//!   `Pᵀ A P`, computed once, and the numeric-only refill on that frozen
//!   pattern, row-parallel over coarse rows with a fixed accumulation order.
//!
//! Everything symbolic is sequential and index-ordered; the numeric refill
//! sums every coarse entry in the same order on every thread count, so the
//! hierarchy is bitwise reproducible.

use crate::backend::cpu::{run_sequential, CHUNK};
use crate::graph::LevelGraph;
use crate::sparse::SparseColMatOwned;
use rayon::prelude::*;
use std::sync::Mutex;

/// Sparse matrix in CSR form with ascending column indices within each row.
///
/// `diag` holds `values` at `(u, u)` for square matrices (`0` where the
/// diagonal is structurally absent) and is empty for rectangular ones. The
/// diagonal entry is *also* stored in `values`, so `apply` needs no special
/// case.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct LevelMatrix {
    /// Rows.
    pub n: usize,
    /// Columns.
    pub ncols: usize,
    /// `n + 1` offsets into `col_idx` / `values`.
    pub row_ptr: Vec<u32>,
    /// Column of every entry, ascending within a row.
    pub col_idx: Vec<u32>,
    /// Value of every entry.
    pub values: Vec<f64>,
    /// `values[(u, u)]` per row for square matrices; empty otherwise.
    pub diag: Vec<f64>,
}

impl LevelMatrix {
    /// Number of stored entries.
    pub fn nnz(&self) -> usize {
        self.col_idx.len()
    }

    /// `true` for a square matrix.
    pub fn is_square(&self) -> bool {
        self.n == self.ncols
    }

    /// Entry range of row `u`.
    #[inline]
    pub fn row_range(&self, u: usize) -> std::ops::Range<usize> {
        self.row_ptr[u] as usize..self.row_ptr[u + 1] as usize
    }

    /// Columns and values of row `u`.
    #[inline]
    pub fn row(&self, u: usize) -> (&[u32], &[f64]) {
        let r = self.row_range(u);
        (&self.col_idx[r.clone()], &self.values[r])
    }

    /// Bytes of the index and value arrays.
    pub fn bytes(&self) -> usize {
        self.row_ptr.len() * 4 + self.col_idx.len() * 4 + (self.values.len() + self.diag.len()) * 8
    }

    /// Validate the CSR invariants (lengths, monotone `row_ptr`, ascending
    /// in-range columns, `diag` length). Panics on violation; meant for
    /// upload paths and tests.
    pub fn check(&self) {
        assert_eq!(
            self.row_ptr.len(),
            self.n + 1,
            "LevelMatrix: row_ptr length"
        );
        assert_eq!(self.row_ptr[0], 0, "LevelMatrix: row_ptr[0]");
        assert_eq!(
            *self.row_ptr.last().unwrap() as usize,
            self.col_idx.len(),
            "LevelMatrix: row_ptr end"
        );
        assert_eq!(
            self.values.len(),
            self.col_idx.len(),
            "LevelMatrix: values length"
        );
        assert!(
            self.diag.is_empty() || self.diag.len() == self.n,
            "LevelMatrix: diag length"
        );
        for u in 0..self.n {
            let r = self.row_range(u);
            assert!(
                r.start <= r.end,
                "LevelMatrix: row_ptr not monotone at row {u}"
            );
            let cols = &self.col_idx[r];
            for (k, &c) in cols.iter().enumerate() {
                assert!(
                    (c as usize) < self.ncols,
                    "LevelMatrix: column {c} ≥ {} in row {u}",
                    self.ncols
                );
                if k > 0 {
                    assert!(
                        cols[k - 1] < c,
                        "LevelMatrix: columns not ascending in row {u}"
                    );
                }
            }
        }
    }

    /// Recompute `diag` from `values` (square matrices only).
    pub fn refresh_diag(&mut self) {
        if !self.is_square() {
            self.diag.clear();
            return;
        }
        self.diag.resize(self.n, 0.0);
        for u in 0..self.n {
            let (cols, vals) = self.row(u);
            self.diag[u] = match cols.binary_search(&(u as u32)) {
                Ok(k) => vals[k],
                Err(_) => 0.0,
            };
        }
    }

    /// Build from per-row `(column, value)` lists: entries are sorted by
    /// column (stable, so duplicates keep their order) and duplicates are
    /// summed in that order.
    pub fn from_rows(nrows: usize, ncols: usize, rows: &mut [Vec<(u32, f64)>]) -> Self {
        assert_eq!(rows.len(), nrows);
        let mut row_ptr = vec![0u32; nrows + 1];
        let total: usize = rows.iter().map(Vec::len).sum();
        let mut col_idx = Vec::with_capacity(total);
        let mut values = Vec::with_capacity(total);
        for (u, row) in rows.iter_mut().enumerate() {
            row.sort_by_key(|e| e.0);
            let row_start = col_idx.len();
            for &(c, v) in row.iter() {
                assert!((c as usize) < ncols, "column out of range");
                if col_idx.len() > row_start && *col_idx.last().unwrap() == c {
                    *values.last_mut().unwrap() += v;
                } else {
                    col_idx.push(c);
                    values.push(v);
                }
            }
            row_ptr[u + 1] = col_idx.len() as u32;
        }
        let mut m = Self {
            n: nrows,
            ncols,
            row_ptr,
            col_idx,
            values,
            diag: Vec::new(),
        };
        m.refresh_diag();
        m
    }

    /// Explicit CSR of a level graph's operator
    /// `a_uu = anchor_u + Σ w`, `a_uv = −Σ w_uv` (multi-edges merged).
    pub fn from_graph(g: &LevelGraph) -> Self {
        let mut rows: Vec<Vec<(u32, f64)>> = (0..g.n)
            .map(|u| {
                let mut row = Vec::with_capacity(g.adjacency.degree(u) + 1);
                LevelRef::Graph(g).for_row(u, |v, a| row.push((v as u32, a)));
                row
            })
            .collect();
        Self::from_rows(g.n, g.n, &mut rows)
    }

    /// Numeric refresh from a level graph with the pattern of
    /// [`Self::from_graph`] (same topology, new weights).
    pub fn set_from_graph(&mut self, g: &LevelGraph) {
        assert_eq!(self.n, g.n, "set_from_graph: node count changed");
        let Self {
            col_idx, values, ..
        } = self;
        for u in 0..g.n {
            let r = self.row_ptr[u] as usize..self.row_ptr[u + 1] as usize;
            let cols = &col_idx[r.clone()];
            let vals = &mut values[r];
            vals.fill(0.0);
            LevelRef::Graph(g).for_row(u, |v, a| {
                let k = cols
                    .binary_search(&(v as u32))
                    .expect("set_from_graph: pattern changed");
                vals[k] += a;
            });
        }
        self.refresh_diag();
    }

    /// The transpose, also with ascending columns per row.
    pub fn transpose(&self) -> Self {
        let mut row_ptr = vec![0u32; self.ncols + 1];
        for &c in &self.col_idx {
            row_ptr[c as usize + 1] += 1;
        }
        for c in 0..self.ncols {
            row_ptr[c + 1] += row_ptr[c];
        }
        let mut fill: Vec<u32> = row_ptr[..self.ncols].to_vec();
        let mut col_idx = vec![0u32; self.nnz()];
        let mut values = vec![0.0; self.nnz()];
        for u in 0..self.n {
            let (cols, vals) = self.row(u);
            for (&c, &v) in cols.iter().zip(vals) {
                let slot = fill[c as usize] as usize;
                fill[c as usize] += 1;
                col_idx[slot] = u as u32;
                values[slot] = v;
            }
        }
        let mut t = Self {
            n: self.ncols,
            ncols: self.n,
            row_ptr,
            col_idx,
            values,
            diag: Vec::new(),
        };
        t.refresh_diag();
        t
    }

    /// Reference serial `y = M x` on `k` interleaved columns
    /// (`x.len() = ncols * k`, `y.len() = n * k`), ascending column order.
    pub fn apply(&self, x: &[f64], y: &mut [f64], k: usize) {
        assert_eq!(x.len(), self.ncols * k, "LevelMatrix::apply: x length");
        assert_eq!(y.len(), self.n * k, "LevelMatrix::apply: y length");
        for u in 0..self.n {
            let (cols, vals) = self.row(u);
            let out = &mut y[u * k..(u + 1) * k];
            out.fill(0.0);
            for (&c, &v) in cols.iter().zip(vals) {
                let xc = &x[c as usize * k..(c as usize + 1) * k];
                for (o, &xi) in out.iter_mut().zip(xc) {
                    *o += v * xi;
                }
            }
        }
    }

    /// Column-major copy for faer (square symmetric: the CSR arrays *are*
    /// the CSC arrays of the transpose, which equals the matrix).
    pub fn to_sparse(&self) -> SparseColMatOwned {
        assert!(self.is_square(), "to_sparse: square matrices only");
        SparseColMatOwned {
            nrows: self.n,
            ncols: self.ncols,
            col_ptrs: self.row_ptr.clone(),
            row_indices: self.col_idx.clone(),
            values: self.values.clone(),
        }
    }

    /// Copy the values into a matrix produced by [`Self::to_sparse`].
    pub fn write_sparse_values(&self, m: &mut SparseColMatOwned) {
        assert_eq!(
            m.values.len(),
            self.values.len(),
            "write_sparse_values: pattern changed"
        );
        m.values.copy_from_slice(&self.values);
    }

    /// Dense copy (tests only).
    pub fn to_dense(&self) -> Vec<Vec<f64>> {
        let mut a = vec![vec![0.0; self.ncols]; self.n];
        for u in 0..self.n {
            let (cols, vals) = self.row(u);
            for (&c, &v) in cols.iter().zip(vals) {
                a[u][c as usize] += v;
            }
        }
        a
    }
}

/// One level operator as rows of `(column, value)` pairs.
#[derive(Debug, Clone, Copy)]
pub enum LevelRef<'a> {
    /// The matrix-free level 0: `a_uu = anchor_u + Σ_{e ∋ u} w_e`,
    /// `a_uv = −w_e` per incident edge (entries of a multi-edge repeat).
    Graph(&'a LevelGraph),
    /// A CSR level.
    Csr(&'a LevelMatrix),
}

impl LevelRef<'_> {
    /// Rows (= columns).
    pub fn n(&self) -> usize {
        match self {
            Self::Graph(g) => g.n,
            Self::Csr(m) => m.n,
        }
    }

    /// Stored entries (graph: nodes + adjacency entries).
    pub fn nnz(&self) -> usize {
        match self {
            Self::Graph(g) => g.n + g.adjacency.num_entries(),
            Self::Csr(m) => m.nnz(),
        }
    }

    /// Jacobi diagonal of row `u`: the CSR diagonal, or for the graph
    /// `anchor_u + Σ_{e ∋ u} w_e` (the diagonal the backend kernels scale
    /// by; it equals `a_uu` unless the graph has self-loops, which carry no
    /// force in an FDM network).
    pub fn diag(&self, u: usize) -> f64 {
        match self {
            Self::Graph(g) => g.diagonal(u),
            Self::Csr(m) => m.diag[u],
        }
    }

    /// Call `f(v, a_uv)` for every stored entry of row `u`. For the graph the
    /// diagonal comes first, then the incident edges in adjacency order (a
    /// column may repeat); for CSR the entries come in ascending column order.
    #[inline]
    pub fn for_row(&self, u: usize, mut f: impl FnMut(usize, f64)) {
        match self {
            Self::Graph(g) => {
                f(u, g.diagonal(u));
                for (e, _, v) in g.adjacency.incident(u) {
                    f(v as usize, -g.weight[e as usize]);
                }
            }
            Self::Csr(m) => {
                let (cols, vals) = m.row(u);
                for (&c, &v) in cols.iter().zip(vals) {
                    f(c as usize, v);
                }
            }
        }
    }

    /// Dense copy (tests only).
    pub fn to_dense(&self) -> Vec<Vec<f64>> {
        match self {
            Self::Graph(g) => LevelMatrix::from_graph(g).to_dense(),
            Self::Csr(m) => m.to_dense(),
        }
    }
}

/// Smoothed-aggregation prolongator `P = (I − ω D⁻¹ A) P₀` as CSR
/// (`n × n_coarse`); `aggregate_of[u] < n_coarse` defines `P₀`.
pub fn smoothed_prolongator(
    a: LevelRef<'_>,
    aggregate_of: &[u32],
    n_coarse: usize,
    omega: f64,
) -> LevelMatrix {
    let n = a.n();
    assert_eq!(
        aggregate_of.len(),
        n,
        "smoothed_prolongator: aggregate map length"
    );
    let mut rows: Vec<Vec<(u32, f64)>> = Vec::with_capacity(n);
    for u in 0..n {
        let mut row: Vec<(u32, f64)> = Vec::with_capacity(8);
        row.push((aggregate_of[u], 1.0));
        let scale = -omega / a.diag(u);
        a.for_row(u, |v, a_uv| row.push((aggregate_of[v], scale * a_uv)));
        rows.push(row);
    }
    LevelMatrix::from_rows(n, n_coarse, &mut rows)
}

/// Symbolic pattern of `Pᵀ A P` (values and diagonal zero): row `I` holds
/// every `J` reachable as `Pᵀ[I, u] · A[u, v] · P[v, J]`, ascending.
pub fn galerkin_pattern(pt: &LevelMatrix, a: LevelRef<'_>, p: &LevelMatrix) -> LevelMatrix {
    let nc = pt.n;
    assert_eq!(pt.ncols, a.n(), "galerkin_pattern: Pᵀ columns ≠ A rows");
    assert_eq!(p.n, a.n(), "galerkin_pattern: P rows ≠ A columns");
    assert_eq!(p.ncols, nc, "galerkin_pattern: P columns ≠ Pᵀ rows");
    let mut mark = vec![u32::MAX; nc];
    let mut row_ptr = vec![0u32; nc + 1];
    let mut col_idx: Vec<u32> = Vec::with_capacity(a.nnz());
    for i in 0..nc {
        let start = col_idx.len();
        let (us, _) = pt.row(i);
        for &u in us {
            a.for_row(u as usize, |v, _| {
                let (js, _) = p.row(v);
                for &j in js {
                    if mark[j as usize] != i as u32 {
                        mark[j as usize] = i as u32;
                        col_idx.push(j);
                    }
                }
            });
        }
        col_idx[start..].sort_unstable();
        row_ptr[i + 1] = col_idx.len() as u32;
    }
    let nnz = col_idx.len();
    LevelMatrix {
        n: nc,
        ncols: nc,
        row_ptr,
        col_idx,
        values: vec![0.0; nnz],
        diag: vec![0.0; nc],
    }
}

/// Pool of dense accumulators for the row-parallel triple product: one per
/// concurrently running chunk, allocated on first use and reused after.
#[derive(Debug, Default)]
pub struct ScratchPool {
    pool: Mutex<Vec<Vec<f64>>>,
}

impl ScratchPool {
    /// Empty pool.
    pub fn new() -> Self {
        Self::default()
    }

    fn take(&self, len: usize) -> Vec<f64> {
        let mut pool = self.pool.lock().unwrap_or_else(|e| e.into_inner());
        match pool.pop() {
            Some(mut v) => {
                if v.len() < len {
                    v.resize(len, 0.0);
                }
                v
            }
            None => vec![0.0; len],
        }
    }

    fn give(&self, v: Vec<f64>) {
        let mut pool = self.pool.lock().unwrap_or_else(|e| e.into_inner());
        pool.push(v);
    }

    /// Bytes held.
    pub fn bytes(&self) -> usize {
        let pool = self.pool.lock().unwrap_or_else(|e| e.into_inner());
        pool.iter().map(|v| v.capacity() * 8).sum()
    }
}

/// Numeric `out = Pᵀ A P` on the frozen pattern of `out` (from
/// [`galerkin_pattern`]). Rows of `out` are processed in chunks of
/// [`CHUNK`] rows in parallel; within a row every contribution is summed in
/// the fixed order `u ∈ row(Pᵀ, I)`, `v ∈ row(A, u)`, `J ∈ row(P, v)`, so
/// the result is bitwise identical on every thread count. `out.diag` is
/// refreshed.
pub fn galerkin_numeric(
    pt: &LevelMatrix,
    a: LevelRef<'_>,
    p: &LevelMatrix,
    out: &mut LevelMatrix,
    scratch: &ScratchPool,
) {
    let nc = out.n;
    debug_assert_eq!(pt.n, nc);
    debug_assert_eq!(p.ncols, nc);
    debug_assert_eq!(out.diag.len(), nc);
    let row_ptr = &out.row_ptr;
    let col_idx = &out.col_idx;

    // Split `values` and `diag` into per-chunk slices aligned on row
    // boundaries (rows are contiguous in `values`).
    let n_chunks = nc.div_ceil(CHUNK);
    let mut pieces: Vec<(usize, &mut [f64], &mut [f64])> = Vec::with_capacity(n_chunks);
    {
        let mut values = out.values.as_mut_slice();
        let mut diag = out.diag.as_mut_slice();
        for c in 0..n_chunks {
            let r0 = c * CHUNK;
            let r1 = ((c + 1) * CHUNK).min(nc);
            let len = (row_ptr[r1] - row_ptr[r0]) as usize;
            let (vals, rest) = values.split_at_mut(len);
            values = rest;
            let (d, rest) = diag.split_at_mut(r1 - r0);
            diag = rest;
            pieces.push((r0, vals, d));
        }
    }

    let run_chunk = |r0: usize, vals: &mut [f64], diag: &mut [f64]| {
        let mut acc = scratch.take(nc);
        let base = row_ptr[r0] as usize;
        for (k, d) in diag.iter_mut().enumerate() {
            let i = r0 + k;
            let (us, pts) = pt.row(i);
            for (&u, &pt_iu) in us.iter().zip(pts) {
                a.for_row(u as usize, |v, a_uv| {
                    let c = pt_iu * a_uv;
                    let (js, ps) = p.row(v);
                    for (&j, &p_vj) in js.iter().zip(ps) {
                        acc[j as usize] += c * p_vj;
                    }
                });
            }
            let r = row_ptr[i] as usize - base..row_ptr[i + 1] as usize - base;
            let cols = &col_idx[row_ptr[i] as usize..row_ptr[i + 1] as usize];
            let mut di = 0.0;
            for (&j, out) in cols.iter().zip(&mut vals[r]) {
                let j = j as usize;
                *out = acc[j];
                if j == i {
                    di = acc[j];
                }
                acc[j] = 0.0;
            }
            *d = di;
        }
        scratch.give(acc);
    };

    if run_sequential(nc) {
        for (r0, vals, diag) in pieces {
            run_chunk(r0, vals, diag);
        }
    } else {
        pieces
            .into_par_iter()
            .for_each(|(r0, vals, diag)| run_chunk(r0, vals, diag));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::CsrAdjacency;

    /// Deterministic xorshift for test data.
    struct Rng(u64);

    impl Rng {
        fn next(&mut self) -> u64 {
            let mut x = self.0;
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            self.0 = x;
            x.wrapping_mul(0x2545_F491_4F6C_DD1D)
        }

        fn uniform(&mut self) -> f64 {
            (self.next() >> 11) as f64 / (1u64 << 53) as f64
        }
    }

    /// Random connected graph on `n` nodes with a chain plus `extra` random
    /// edges (multi-edges allowed, no self-loops), random positive weights
    /// and anchors on a few nodes.
    pub(crate) fn random_graph(n: usize, extra: usize, seed: u64) -> LevelGraph {
        let mut rng = Rng(seed.max(1));
        let mut starts = Vec::new();
        let mut ends = Vec::new();
        for i in 0..n - 1 {
            starts.push(i);
            ends.push(i + 1);
        }
        while starts.len() < n - 1 + extra {
            let (s, e) = (
                (rng.next() % n as u64) as usize,
                (rng.next() % n as u64) as usize,
            );
            if s != e {
                starts.push(s);
                ends.push(e);
            }
        }
        let ne = starts.len();
        let weight: Vec<f64> = (0..ne).map(|_| 0.1 + 2.0 * rng.uniform()).collect();
        let anchor: Vec<f64> = (0..n)
            .map(|u| if u % 7 == 3 { 0.5 + rng.uniform() } else { 0.0 })
            .collect();
        LevelGraph {
            n,
            adjacency: CsrAdjacency::from_endpoints(n, &starts, &ends),
            weight,
            anchor,
            aggregate_of: Vec::new(),
        }
    }

    fn dense_mul(a: &[Vec<f64>], b: &[Vec<f64>]) -> Vec<Vec<f64>> {
        let (n, k, m) = (a.len(), b.len(), b[0].len());
        let mut c = vec![vec![0.0; m]; n];
        for i in 0..n {
            for l in 0..k {
                let ail = a[i][l];
                if ail != 0.0 {
                    for j in 0..m {
                        c[i][j] += ail * b[l][j];
                    }
                }
            }
        }
        c
    }

    fn dense_transpose(a: &[Vec<f64>]) -> Vec<Vec<f64>> {
        let (n, m) = (a.len(), a[0].len());
        (0..m).map(|j| (0..n).map(|i| a[i][j]).collect()).collect()
    }

    fn max_rel_diff(a: &[Vec<f64>], b: &[Vec<f64>]) -> f64 {
        let scale = a.iter().flatten().fold(0.0f64, |m, v| m.max(v.abs()));
        a.iter()
            .flatten()
            .zip(b.iter().flatten())
            .map(|(x, y)| (x - y).abs() / scale)
            .fold(0.0, f64::max)
    }

    #[test]
    fn from_graph_matches_operator_and_merges_multi_edges() {
        // 0 -1- 1 -2- 2 plus a second edge 0-2 twice and a self-loop at 1.
        let g = LevelGraph {
            n: 3,
            adjacency: CsrAdjacency::from_endpoints(3, &[0, 1, 0, 2, 1], &[1, 2, 2, 0, 1]),
            weight: vec![1.0, 2.0, 3.0, 4.0, 5.0],
            anchor: vec![0.5, 0.0, 0.0],
            aggregate_of: Vec::new(),
        };
        let m = LevelMatrix::from_graph(&g);
        m.check();
        assert_eq!(m.n, 3);
        // Row 0: diag 0.5 + 1 + 3 + 4 = 8.5, col 1: −1, col 2: −7.
        assert_eq!(m.row(0), (&[0, 1, 2][..], &[8.5, -1.0, -7.0][..]));
        // Row 1: diag 1 + 2 + 5 + 5 (self-loop twice) − 5 − 5 = 3; col 0: −1, col 2: −2.
        assert_eq!(m.row(1), (&[0, 1, 2][..], &[-1.0, 3.0, -2.0][..]));
        assert_eq!(m.diag, vec![8.5, 3.0, 9.0]);
        // Same result as the gather operator on random x.
        let x: Vec<f64> = (0..9).map(|i| (i as f64 * 0.7).sin()).collect();
        let mut y = vec![0.0; 9];
        m.apply(&x, &mut y, 3);
        let xb = crate::graph::BlockVec::<3>::from_vec(x.clone());
        let mut yb = crate::graph::BlockVec::<3>::zeros(3);
        g.apply(&xb, &mut yb);
        for (a, b) in y.iter().zip(yb.as_slice()) {
            assert!((a - b).abs() < 1e-14, "{a} vs {b}");
        }
        // Numeric refresh reproduces from_graph on new weights.
        let mut g2 = g.clone();
        g2.weight.iter_mut().for_each(|w| *w *= 1.5);
        g2.anchor[2] = 2.0;
        let mut m2 = m.clone();
        m2.set_from_graph(&g2);
        assert_eq!(m2, LevelMatrix::from_graph(&g2));
        // Transpose of a symmetric matrix is itself.
        assert_eq!(m.transpose(), m);
        // to_sparse is the same pattern in CSC.
        let s = m.to_sparse();
        assert_eq!(s.col_ptrs, m.row_ptr);
        assert_eq!(s.row_indices, m.col_idx);
    }

    #[test]
    fn transpose_of_rectangular_matrix() {
        let mut rows = vec![
            vec![(2u32, 1.0), (0, 2.0)],
            vec![],
            vec![(1u32, 3.0), (1, 4.0)],
        ];
        let m = LevelMatrix::from_rows(3, 3, &mut rows);
        m.check();
        assert_eq!(m.row(0), (&[0, 2][..], &[2.0, 1.0][..]));
        assert_eq!(m.row(2), (&[1][..], &[7.0][..]));
        let t = m.transpose();
        t.check();
        assert_eq!(t.row(0), (&[0][..], &[2.0][..]));
        assert_eq!(t.row(1), (&[2][..], &[7.0][..]));
        assert_eq!(t.row(2), (&[0][..], &[1.0][..]));
        assert_eq!(dense_transpose(&m.to_dense()), t.to_dense());
    }

    #[test]
    fn galerkin_product_matches_dense_ptap() {
        use crate::amg::aggregate::{aggregate, StrengthGraph};
        for (n, extra, seed) in [(40usize, 30usize, 1u64), (120, 200, 2), (200, 100, 3)] {
            let g = random_graph(n, extra, seed);
            let a0 = LevelRef::Graph(&g);
            let sg = StrengthGraph::from_level_graph(&g);
            let (agg, nc) = aggregate(&sg, 2);
            let omega = 4.0 / (3.0 * 1.9);
            let p = smoothed_prolongator(a0, &agg, nc, omega);
            p.check();
            assert_eq!((p.n, p.ncols), (n, nc));
            // P dense = (I − ω D⁻¹ A) P₀.
            let a_dense = a0.to_dense();
            let mut p0 = vec![vec![0.0; nc]; n];
            for (u, &c) in agg.iter().enumerate() {
                p0[u][c as usize] = 1.0;
            }
            let mut s = a_dense.clone();
            for u in 0..n {
                let d = a_dense[u][u];
                for v in 0..n {
                    s[u][v] = -omega * a_dense[u][v] / d;
                }
                s[u][u] += 1.0;
            }
            let p_dense = dense_mul(&s, &p0);
            let pdiff = max_rel_diff(&p_dense, &p.to_dense());
            assert!(pdiff < 1e-13, "P mismatch {pdiff}");

            let pt = p.transpose();
            let mut a1 = galerkin_pattern(&pt, a0, &p);
            let scratch = ScratchPool::new();
            galerkin_numeric(&pt, a0, &p, &mut a1, &scratch);
            a1.check();
            let expected = dense_mul(&dense_transpose(&p_dense), &dense_mul(&a_dense, &p_dense));
            assert!(
                max_rel_diff(&expected, &a1.to_dense()) < 1e-12,
                "level 1 mismatch {}",
                max_rel_diff(&expected, &a1.to_dense())
            );
            for u in 0..nc {
                assert_eq!(a1.diag[u], a1.to_dense()[u][u]);
            }
            // Symmetric pattern and values.
            assert!(max_rel_diff(&a1.to_dense(), &dense_transpose(&a1.to_dense())) < 1e-12);

            // Second level from the CSR operator.
            let sg1 = StrengthGraph::from_matrix(&a1);
            let (agg1, nc1) = aggregate(&sg1, 2);
            let p1 = smoothed_prolongator(LevelRef::Csr(&a1), &agg1, nc1, omega);
            let pt1 = p1.transpose();
            let mut a2 = galerkin_pattern(&pt1, LevelRef::Csr(&a1), &p1);
            galerkin_numeric(&pt1, LevelRef::Csr(&a1), &p1, &mut a2, &scratch);
            let p1d = p1.to_dense();
            let expected2 = dense_mul(&dense_transpose(&p1d), &dense_mul(&a1.to_dense(), &p1d));
            assert!(max_rel_diff(&expected2, &a2.to_dense()) < 1e-12);

            // Numeric refill with new weights on the frozen pattern.
            let mut g2 = g.clone();
            for (k, w) in g2.weight.iter_mut().enumerate() {
                *w *= 1.0 + 0.02 * ((k % 5) as f64 - 2.0);
            }
            galerkin_numeric(&pt, LevelRef::Graph(&g2), &p, &mut a1, &scratch);
            let expected = dense_mul(
                &dense_transpose(&p_dense),
                &dense_mul(&LevelRef::Graph(&g2).to_dense(), &p_dense),
            );
            assert!(max_rel_diff(&expected, &a1.to_dense()) < 1e-12);
            // Refill is idempotent (scratch returned zeroed).
            let before = a1.clone();
            galerkin_numeric(&pt, LevelRef::Graph(&g2), &p, &mut a1, &scratch);
            assert_eq!(before, a1);
        }
    }

    #[test]
    fn galerkin_numeric_is_thread_count_independent_in_structure() {
        // More rows than one chunk so the parallel path splits; the result
        // must equal the sequential one (same code path per row).
        use crate::amg::aggregate::{aggregate, StrengthGraph};
        let n = 3 * CHUNK + 100;
        let g = random_graph(n, n, 9);
        let a0 = LevelRef::Graph(&g);
        let (agg, nc) = aggregate(&StrengthGraph::from_level_graph(&g), 1);
        let p = smoothed_prolongator(a0, &agg, nc, 0.7);
        let pt = p.transpose();
        let mut a1 = galerkin_pattern(&pt, a0, &p);
        let scratch = ScratchPool::new();
        galerkin_numeric(&pt, a0, &p, &mut a1, &scratch);
        // Reference: one row at a time with a fresh accumulator.
        let mut acc = vec![0.0; nc];
        for i in 0..nc {
            let (us, pts) = pt.row(i);
            for (&u, &pt_iu) in us.iter().zip(pts) {
                a0.for_row(u as usize, |v, a_uv| {
                    let c = pt_iu * a_uv;
                    let (js, ps) = p.row(v);
                    for (&j, &p_vj) in js.iter().zip(ps) {
                        acc[j as usize] += c * p_vj;
                    }
                });
            }
            let (cols, vals) = a1.row(i);
            for (&j, &v) in cols.iter().zip(vals) {
                assert_eq!(v, acc[j as usize], "row {i} col {j}");
                acc[j as usize] = 0.0;
            }
        }
    }
}

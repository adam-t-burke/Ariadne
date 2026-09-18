//! Phase 0 of the iterative-solver program: aggregation AMG prototype.
//!
//! A self-contained, single-threaded prototype of the matrix-free solver
//! described in `AMG_PLAN.md` and `ITERATIVE_SOLVER_PROGRAM.md` §3, used to
//! decide the algorithm choices before any production code is written:
//!
//! * the free-node weighted graph Laplacian `A(q)` in node-gather (CSR
//!   incidence) form, built from a `Problem`'s topology plus `q`;
//! * pairwise heavy-edge aggregation (1–3 passes per level) with the coarse
//!   operator as the quotient graph (`w_c = Σ q`, `anchor_U = Σ anchor_u`);
//! * a Chebyshev smoother (degree 1–3) with Jacobi scaling on the spectral
//!   interval `[λ_max/α, λ_max]`, `λ_max` from power iteration;
//! * V-cycle and K-cycle (Notay: up to two flexible-CG steps at each coarse
//!   level preconditioned by the next-coarser cycle), coarsest level solved
//!   with faer's sparse Cholesky through `theseus::types::Factorization`;
//! * PCG (V-cycle) and flexible CG FCG(1) (K-cycle) on the block of three
//!   right-hand sides with per-column convergence;
//! * optional smoothed aggregation (`P = (I − ω D⁻¹A) P₀`), whose coarse
//!   operators are general sparse matrices, as the documented fallback.
//!
//! Fixtures: the grid from `tests/support/grid.rs`, an irregular
//! triangulated mesh (jittered lattice, symmetric k-nearest-neighbour
//! connectivity, boundary ring fixed) and a cable dome / spoke wheel
//! (concentric rings, radial spokes, diagonal bracing, ring doubling, outer
//! ring fixed), each with a smooth non-constant `q` field.
//!
//! ```text
//! export RUSTUP_TOOLCHAIN=1.89.0
//! RAYON_NUM_THREADS=1 cargo run -p theseus --release --example amg_prototype -- \
//!     --fixture grid --size 224 --cycle v,k --degree 1,2,3 --passes 1,2 --warm 2 --reps 3
//! cargo run -p theseus --release --example amg_prototype -- --self-test
//! ```
//!
//! Results and the recommendation are recorded in `BENCHMARKS.md`
//! ("Phase 0: aggregation AMG prototype").

#[path = "../tests/support/grid.rs"]
#[allow(dead_code)]
mod grid;

use ndarray::Array2;
use std::cell::Cell;
use std::time::Instant;
use theseus::sparse::SparseColMatOwned;
use theseus::types::*;

const NONE: u32 = u32::MAX;

// ─────────────────────────────────────────────────────────────
//  Small helpers: deterministic RNG and n×3 block-vector kernels
// ─────────────────────────────────────────────────────────────

struct Rng(u64);

impl Rng {
    fn new(seed: u64) -> Self {
        Self(seed.max(1))
    }

    fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        x.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }

    /// Uniform in `[0, 1)`.
    fn uniform(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    /// Uniform in `[-1, 1)`.
    fn symmetric(&mut self) -> f64 {
        2.0 * self.uniform() - 1.0
    }

    fn below(&mut self, n: usize) -> usize {
        (self.next_u64() % n as u64) as usize
    }
}

fn dot3(a: &[f64], b: &[f64]) -> [f64; 3] {
    let mut s = [0.0; 3];
    for (x, y) in a.chunks_exact(3).zip(b.chunks_exact(3)) {
        s[0] += x[0] * y[0];
        s[1] += x[1] * y[1];
        s[2] += x[2] * y[2];
    }
    s
}

fn norm3(a: &[f64]) -> [f64; 3] {
    dot3(a, a).map(f64::sqrt)
}

fn dot1(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn safe_div(a: f64, b: f64) -> f64 {
    if b == 0.0 {
        0.0
    } else {
        a / b
    }
}

fn div3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        safe_div(a[0], b[0]),
        safe_div(a[1], b[1]),
        safe_div(a[2], b[2]),
    ]
}

/// `y += alpha_d x` per column.
fn axpy3(y: &mut [f64], alpha: [f64; 3], x: &[f64]) {
    for (yy, xx) in y.chunks_exact_mut(3).zip(x.chunks_exact(3)) {
        yy[0] += alpha[0] * xx[0];
        yy[1] += alpha[1] * xx[1];
        yy[2] += alpha[2] * xx[2];
    }
}

fn restrict(aggregate_of: &[u32], fine: &[f64], coarse: &mut [f64]) {
    coarse.fill(0.0);
    for (u, &c) in aggregate_of.iter().enumerate() {
        let c = c as usize * 3;
        coarse[c] += fine[u * 3];
        coarse[c + 1] += fine[u * 3 + 1];
        coarse[c + 2] += fine[u * 3 + 2];
    }
}

fn prolong_add(aggregate_of: &[u32], coarse: &[f64], fine: &mut [f64]) {
    for (u, &c) in aggregate_of.iter().enumerate() {
        let c = c as usize * 3;
        fine[u * 3] += coarse[c];
        fine[u * 3 + 1] += coarse[c + 1];
        fine[u * 3 + 2] += coarse[c + 2];
    }
}

fn median(v: &mut [f64]) -> f64 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v[v.len() / 2]
}

fn ms(t: Instant) -> f64 {
    t.elapsed().as_secs_f64() * 1e3
}

// ─────────────────────────────────────────────────────────────
//  Level operators
// ─────────────────────────────────────────────────────────────

/// One multigrid level as a weighted graph on free nodes with anchor
/// weights: `(A x)_u = anchor_u x_u + Σ_{(u,v)} w (x_u − x_v)`.
///
/// Edges are stored twice: as an edge list (for coarsening and the coarse
/// weight update) and as a CSR node-gather adjacency (for the operator).
struct LevelGraph {
    n: usize,
    eu: Vec<u32>,
    ev: Vec<u32>,
    weight: Vec<f64>,
    anchor: Vec<f64>,
    offsets: Vec<u32>,
    other: Vec<u32>,
    edge: Vec<u32>,
    /// `anchor_u + Σ w` per node; refreshed with the weights.
    diag: Vec<f64>,
    /// Edge weight per adjacency entry; refreshed with the weights.
    adj_w: Vec<f64>,
}

impl LevelGraph {
    fn new(n: usize, eu: Vec<u32>, ev: Vec<u32>, weight: Vec<f64>, anchor: Vec<f64>) -> Self {
        let ne = eu.len();
        let mut offsets = vec![0u32; n + 1];
        for j in 0..ne {
            offsets[eu[j] as usize + 1] += 1;
            offsets[ev[j] as usize + 1] += 1;
        }
        for i in 0..n {
            offsets[i + 1] += offsets[i];
        }
        let mut fill = offsets.clone();
        let mut other = vec![0u32; 2 * ne];
        let mut edge = vec![0u32; 2 * ne];
        for j in 0..ne {
            let (u, v) = (eu[j] as usize, ev[j] as usize);
            let slot = fill[u] as usize;
            other[slot] = v as u32;
            edge[slot] = j as u32;
            fill[u] += 1;
            let slot = fill[v] as usize;
            other[slot] = u as u32;
            edge[slot] = j as u32;
            fill[v] += 1;
        }
        let mut g = Self {
            n,
            eu,
            ev,
            weight,
            anchor,
            offsets,
            other,
            edge,
            diag: vec![0.0; n],
            adj_w: vec![0.0; 2 * ne],
        };
        g.refresh();
        g
    }

    fn ne(&self) -> usize {
        self.eu.len()
    }

    fn degree(&self, u: usize) -> usize {
        (self.offsets[u + 1] - self.offsets[u]) as usize
    }

    /// Recompute `diag` and the per-entry weights from `weight`/`anchor`.
    fn refresh(&mut self) {
        for u in 0..self.n {
            let mut d = self.anchor[u];
            for j in self.offsets[u] as usize..self.offsets[u + 1] as usize {
                let w = self.weight[self.edge[j] as usize];
                self.adj_w[j] = w;
                d += w;
            }
            self.diag[u] = d;
        }
    }

    /// `y = A x` for `K` interleaved columns, node-gather form.
    fn apply<const K: usize>(&self, x: &[f64], y: &mut [f64]) {
        for u in 0..self.n {
            let du = self.diag[u];
            let mut acc = [0.0; K];
            for k in 0..K {
                acc[k] = du * x[u * K + k];
            }
            for j in self.offsets[u] as usize..self.offsets[u + 1] as usize {
                let w = self.adj_w[j];
                let o = self.other[j] as usize * K;
                for k in 0..K {
                    acc[k] -= w * x[o + k];
                }
            }
            y[u * K..u * K + K].copy_from_slice(&acc);
        }
    }

    /// Dense copy (tests only).
    fn dense(&self) -> Vec<Vec<f64>> {
        let mut a = vec![vec![0.0; self.n]; self.n];
        for u in 0..self.n {
            a[u][u] = self.diag[u];
            for j in self.offsets[u] as usize..self.offsets[u + 1] as usize {
                a[u][self.other[j] as usize] -= self.adj_w[j];
            }
        }
        a
    }

    /// Full symmetric CSC matrix (coarsest-level factorization).
    fn to_sparse(&self) -> SparseColMatOwned {
        let mut triplets = Vec::with_capacity(self.n + 2 * self.ne());
        for u in 0..self.n {
            triplets.push((u as u32, u as u32, self.diag[u]));
        }
        for j in 0..self.ne() {
            let w = self.weight[j];
            triplets.push((self.eu[j], self.ev[j], -w));
            triplets.push((self.ev[j], self.eu[j], -w));
        }
        SparseColMatOwned::from_triplets(self.n, self.n, &triplets).unwrap()
    }
}

/// General sparse symmetric level (smoothed aggregation only): CSR with the
/// diagonal stored separately, `(A x)_u = diag_u x_u + Σ_j val_j x_{col_j}`.
struct SparseLevel {
    n: usize,
    offsets: Vec<u32>,
    cols: Vec<u32>,
    vals: Vec<f64>,
    diag: Vec<f64>,
}

impl SparseLevel {
    fn nnz(&self) -> usize {
        self.cols.len() + self.n
    }

    fn apply<const K: usize>(&self, x: &[f64], y: &mut [f64]) {
        for u in 0..self.n {
            let du = self.diag[u];
            let mut acc = [0.0; K];
            for k in 0..K {
                acc[k] = du * x[u * K + k];
            }
            for j in self.offsets[u] as usize..self.offsets[u + 1] as usize {
                let w = self.vals[j];
                let o = self.cols[j] as usize * K;
                for k in 0..K {
                    acc[k] += w * x[o + k];
                }
            }
            y[u * K..u * K + K].copy_from_slice(&acc);
        }
    }

    fn from_triplets(n: usize, mut triplets: Vec<(u32, u32, f64)>) -> Self {
        triplets.sort_unstable_by(|a, b| (a.0, a.1).cmp(&(b.0, b.1)));
        let mut offsets = vec![0u32; n + 1];
        let mut cols = Vec::with_capacity(triplets.len());
        let mut vals = Vec::with_capacity(triplets.len());
        let mut diag = vec![0.0; n];
        let mut i = 0;
        while i < triplets.len() {
            let (r, c, _) = triplets[i];
            let mut sum = 0.0;
            while i < triplets.len() && triplets[i].0 == r && triplets[i].1 == c {
                sum += triplets[i].2;
                i += 1;
            }
            if r == c {
                diag[r as usize] = sum;
            } else if sum != 0.0 {
                cols.push(c);
                vals.push(sum);
                offsets[r as usize + 1] += 1;
            }
        }
        for u in 0..n {
            offsets[u + 1] += offsets[u];
        }
        Self {
            n,
            offsets,
            cols,
            vals,
            diag,
        }
    }

    fn to_sparse(&self) -> SparseColMatOwned {
        let mut triplets = Vec::with_capacity(self.nnz());
        for u in 0..self.n {
            triplets.push((u as u32, u as u32, self.diag[u]));
            for j in self.offsets[u] as usize..self.offsets[u + 1] as usize {
                triplets.push((u as u32, self.cols[j], self.vals[j]));
            }
        }
        SparseColMatOwned::from_triplets(self.n, self.n, &triplets).unwrap()
    }

    fn dense(&self) -> Vec<Vec<f64>> {
        let mut a = vec![vec![0.0; self.n]; self.n];
        for u in 0..self.n {
            a[u][u] = self.diag[u];
            for j in self.offsets[u] as usize..self.offsets[u + 1] as usize {
                a[u][self.cols[j] as usize] += self.vals[j];
            }
        }
        a
    }
}

/// Prolongation with a general sparse pattern (smoothed aggregation):
/// row `u` of `P` has entries `cols[offsets[u]..offsets[u+1]]`.
struct SparseProlongation {
    nfine: usize,
    offsets: Vec<u32>,
    cols: Vec<u32>,
    vals: Vec<f64>,
}

impl SparseProlongation {
    /// `coarse = Pᵀ fine`.
    fn restrict(&self, fine: &[f64], coarse: &mut [f64]) {
        coarse.fill(0.0);
        for u in 0..self.nfine {
            let f = &fine[u * 3..u * 3 + 3];
            for j in self.offsets[u] as usize..self.offsets[u + 1] as usize {
                let c = self.cols[j] as usize * 3;
                let p = self.vals[j];
                coarse[c] += p * f[0];
                coarse[c + 1] += p * f[1];
                coarse[c + 2] += p * f[2];
            }
        }
    }

    /// `fine += P coarse`.
    fn prolong_add(&self, coarse: &[f64], fine: &mut [f64]) {
        for u in 0..self.nfine {
            let mut acc = [0.0; 3];
            for j in self.offsets[u] as usize..self.offsets[u + 1] as usize {
                let c = self.cols[j] as usize * 3;
                let p = self.vals[j];
                acc[0] += p * coarse[c];
                acc[1] += p * coarse[c + 1];
                acc[2] += p * coarse[c + 2];
            }
            fine[u * 3] += acc[0];
            fine[u * 3 + 1] += acc[1];
            fine[u * 3 + 2] += acc[2];
        }
    }
}

enum Operator {
    Graph(LevelGraph),
    Sparse(SparseLevel),
}

impl Operator {
    fn n(&self) -> usize {
        match self {
            Self::Graph(g) => g.n,
            Self::Sparse(s) => s.n,
        }
    }

    /// Work of one application, in edge-equivalents.
    fn cost(&self) -> f64 {
        match self {
            Self::Graph(g) => g.ne() as f64,
            Self::Sparse(s) => 0.5 * s.nnz() as f64,
        }
    }

    fn diag(&self) -> &[f64] {
        match self {
            Self::Graph(g) => &g.diag,
            Self::Sparse(s) => &s.diag,
        }
    }

    fn apply<const K: usize>(&self, x: &[f64], y: &mut [f64]) {
        match self {
            Self::Graph(g) => g.apply::<K>(x, y),
            Self::Sparse(s) => s.apply::<K>(x, y),
        }
    }

    fn to_sparse(&self) -> SparseColMatOwned {
        match self {
            Self::Graph(g) => g.to_sparse(),
            Self::Sparse(s) => s.to_sparse(),
        }
    }

    fn dense(&self) -> Vec<Vec<f64>> {
        match self {
            Self::Graph(g) => g.dense(),
            Self::Sparse(s) => s.dense(),
        }
    }
}

enum Transfer {
    /// Piecewise-constant prolongation: fine node → aggregate.
    Aggregate {
        aggregate_of: Vec<u32>,
        /// Fine edges feeding each coarse edge (gather form for `update`).
        fine_offsets: Vec<u32>,
        fine_list: Vec<u32>,
    },
    Smoothed {
        aggregate_of: Vec<u32>,
        p: SparseProlongation,
    },
}

impl Transfer {
    fn aggregate_of(&self) -> &[u32] {
        match self {
            Self::Aggregate { aggregate_of, .. } | Self::Smoothed { aggregate_of, .. } => {
                aggregate_of
            }
        }
    }

    fn restrict(&self, fine: &[f64], coarse: &mut [f64]) {
        match self {
            Self::Aggregate { aggregate_of, .. } => restrict(aggregate_of, fine, coarse),
            Self::Smoothed { p, .. } => p.restrict(fine, coarse),
        }
    }

    fn prolong_add(&self, coarse: &[f64], fine: &mut [f64]) {
        match self {
            Self::Aggregate { aggregate_of, .. } => prolong_add(aggregate_of, coarse, fine),
            Self::Smoothed { p, .. } => p.prolong_add(coarse, fine),
        }
    }

    /// Dense `P` (tests only).
    fn dense(&self, nfine: usize, ncoarse: usize) -> Vec<Vec<f64>> {
        let mut p = vec![vec![0.0; ncoarse]; nfine];
        match self {
            Self::Aggregate { aggregate_of, .. } => {
                for (u, &c) in aggregate_of.iter().enumerate() {
                    p[u][c as usize] = 1.0;
                }
            }
            Self::Smoothed { p: sp, .. } => {
                for u in 0..nfine {
                    for j in sp.offsets[u] as usize..sp.offsets[u + 1] as usize {
                        p[u][sp.cols[j] as usize] += sp.vals[j];
                    }
                }
            }
        }
        p
    }
}

// ─────────────────────────────────────────────────────────────
//  Aggregation
// ─────────────────────────────────────────────────────────────

/// One pass of pairwise heavy-edge matching (program document §3):
/// strength `s_uv = w_uv / sqrt(d_u d_v)`, greedy in decreasing node degree
/// with lowest-index tie-break; unmatched nodes join the strongest
/// neighbour's aggregate if it has fewer than three members, otherwise they
/// form singletons. Returns `(aggregate_of, n_coarse)`.
fn pairwise_match(g: &LevelGraph) -> (Vec<u32>, usize) {
    let n = g.n;
    let mut order: Vec<u32> = (0..n as u32).collect();
    order.sort_by_key(|&u| (std::cmp::Reverse(g.degree(u as usize)), u));
    let mut agg = vec![NONE; n];
    let mut members: Vec<u8> = Vec::with_capacity(n / 2 + 1);

    let strongest = |u: usize, accept: &dyn Fn(usize) -> bool| -> Option<usize> {
        let mut best: Option<usize> = None;
        let mut best_s = 0.0;
        for j in g.offsets[u] as usize..g.offsets[u + 1] as usize {
            let v = g.other[j] as usize;
            if v == u || !accept(v) {
                continue;
            }
            let s = g.adj_w[j] / (g.diag[u] * g.diag[v]).sqrt();
            if s > best_s || (s == best_s && best.is_some_and(|b| v < b)) {
                best = Some(v);
                best_s = s;
            }
        }
        best
    };

    for &u in &order {
        let u = u as usize;
        if agg[u] != NONE {
            continue;
        }
        if let Some(v) = strongest(u, &|v| agg[v] == NONE) {
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
        match strongest(u, &|v| agg[v] != NONE && members[agg[v] as usize] < 3) {
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

/// Quotient graph of `g` under `aggregate_of`: inter-aggregate edges summed,
/// anchors summed. Also returns the fine-edge gather lists per coarse edge.
fn coarsen(g: &LevelGraph, aggregate_of: &[u32], nc: usize) -> (LevelGraph, Vec<u32>, Vec<u32>) {
    let mut pairs: Vec<(u64, u32)> = Vec::with_capacity(g.ne());
    for j in 0..g.ne() {
        let a = aggregate_of[g.eu[j] as usize];
        let b = aggregate_of[g.ev[j] as usize];
        if a != b {
            let (lo, hi) = if a < b { (a, b) } else { (b, a) };
            pairs.push((((lo as u64) << 32) | hi as u64, j as u32));
        }
    }
    pairs.sort_unstable();
    let mut eu = Vec::new();
    let mut ev = Vec::new();
    let mut weight = Vec::new();
    let mut fine_offsets = vec![0u32];
    let mut fine_list = Vec::with_capacity(pairs.len());
    let mut last = u64::MAX;
    for &(key, j) in &pairs {
        if key != last {
            eu.push((key >> 32) as u32);
            ev.push((key & 0xffff_ffff) as u32);
            weight.push(0.0);
            fine_offsets.push(fine_list.len() as u32);
            last = key;
        }
        *weight.last_mut().unwrap() += g.weight[j as usize];
        fine_list.push(j);
        *fine_offsets.last_mut().unwrap() = fine_list.len() as u32;
    }
    let mut anchor = vec![0.0; nc];
    for (u, &c) in aggregate_of.iter().enumerate() {
        anchor[c as usize] += g.anchor[u];
    }
    (
        LevelGraph::new(nc, eu, ev, weight, anchor),
        fine_offsets,
        fine_list,
    )
}

/// `passes` matching passes composed into one fine → coarse map.
fn aggregate(g: &LevelGraph, passes: usize) -> (Vec<u32>, usize) {
    let (mut agg, mut nc) = pairwise_match(g);
    for _ in 1..passes {
        let (coarse, _, _) = coarsen(g, &agg, nc);
        let (agg2, nc2) = pairwise_match(&coarse);
        for a in agg.iter_mut() {
            *a = agg2[*a as usize];
        }
        nc = nc2;
    }
    (agg, nc)
}

/// Smoothed-aggregation prolongation `P = (I − ω D⁻¹ A) P₀` and the Galerkin
/// coarse operator `Pᵀ A P` as a general sparse matrix.
fn smoothed_transfer(
    op: &Operator,
    aggregate_of: &[u32],
    nc: usize,
    omega: f64,
) -> (SparseProlongation, SparseLevel) {
    let n = op.n();
    let diag = op.diag();
    // Row u of P: P₀ contributes 1 at agg[u]; −ω D⁻¹ A contributes
    // −ω a_uv / d_u at agg[v] for every v (including u).
    let mut rows: Vec<Vec<(u32, f64)>> = vec![Vec::new(); n];
    let mut push = |u: usize, c: u32, v: f64| {
        let row = &mut rows[u];
        match row.iter_mut().find(|(cc, _)| *cc == c) {
            Some(entry) => entry.1 += v,
            None => row.push((c, v)),
        }
    };
    match op {
        Operator::Graph(g) => {
            for u in 0..n {
                push(u, aggregate_of[u], 1.0 - omega);
                let scale = omega / diag[u];
                for j in g.offsets[u] as usize..g.offsets[u + 1] as usize {
                    let v = g.other[j] as usize;
                    push(u, aggregate_of[v], scale * g.adj_w[j]);
                }
            }
        }
        Operator::Sparse(s) => {
            for u in 0..n {
                push(u, aggregate_of[u], 1.0 - omega);
                let scale = -omega / diag[u];
                for j in s.offsets[u] as usize..s.offsets[u + 1] as usize {
                    let v = s.cols[j] as usize;
                    push(u, aggregate_of[v], scale * s.vals[j]);
                }
            }
        }
    }
    let mut offsets = vec![0u32; n + 1];
    let mut cols = Vec::new();
    let mut vals = Vec::new();
    for u in 0..n {
        rows[u].sort_by_key(|e| e.0);
        for &(c, v) in &rows[u] {
            cols.push(c);
            vals.push(v);
        }
        offsets[u + 1] = cols.len() as u32;
    }
    let p = SparseProlongation {
        nfine: n,
        offsets,
        cols,
        vals,
    };
    // Coarse operator: Pᵀ A P via (A P) rows then Pᵀ.
    let ap: Vec<Vec<(u32, f64)>> = (0..n)
        .map(|u| {
            let mut row: Vec<(u32, f64)> = Vec::new();
            let mut add = |c: u32, v: f64| match row.iter_mut().find(|(cc, _)| *cc == c) {
                Some(e) => e.1 += v,
                None => row.push((c, v)),
            };
            let mut visit = |v: usize, a_uv: f64| {
                for j in p.offsets[v] as usize..p.offsets[v + 1] as usize {
                    add(p.cols[j], a_uv * p.vals[j]);
                }
            };
            visit(u, diag[u]);
            match op {
                Operator::Graph(g) => {
                    for j in g.offsets[u] as usize..g.offsets[u + 1] as usize {
                        visit(g.other[j] as usize, -g.adj_w[j]);
                    }
                }
                Operator::Sparse(s) => {
                    for j in s.offsets[u] as usize..s.offsets[u + 1] as usize {
                        visit(s.cols[j] as usize, s.vals[j]);
                    }
                }
            }
            row
        })
        .collect();
    let mut triplets = Vec::new();
    for u in 0..n {
        for j in p.offsets[u] as usize..p.offsets[u + 1] as usize {
            let r = p.cols[j];
            let pv = p.vals[j];
            for &(c, v) in &ap[u] {
                triplets.push((r, c, pv * v));
            }
        }
    }
    (p, SparseLevel::from_triplets(nc, triplets))
}

// ─────────────────────────────────────────────────────────────
//  Hierarchy, smoother, cycles, outer Krylov loops
// ─────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Cycle {
    V,
    K,
}

#[derive(Clone, Debug)]
struct Params {
    cycle: Cycle,
    degree: usize,
    passes: usize,
    coarsest: usize,
    alpha: f64,
    /// K-cycle: skip the second inner step when `‖r̃‖ ≤ ktol ‖r‖`.
    ktol: f64,
    smoothed: bool,
    /// Over-correction factor applied to the coarse-grid correction.
    over: f64,
    power_iterations: usize,
}

impl Params {
    fn label(&self) -> String {
        format!(
            "{}{} d{} p{} a{:.0}{}",
            match self.cycle {
                Cycle::V => "V",
                Cycle::K => "K",
            },
            if self.cycle == Cycle::K {
                format!("(t{})", self.ktol)
            } else {
                String::new()
            },
            self.degree,
            self.passes,
            self.alpha,
            if self.smoothed { " SA" } else { "" }
        ) + &if self.over != 1.0 {
            format!(" o{}", self.over)
        } else {
            String::new()
        }
    }
}

/// Free-node view of a `Problem`'s topology.
struct Topo {
    starts: Vec<usize>,
    ends: Vec<usize>,
    node_to_free: Vec<Option<usize>>,
    nfree: usize,
}

impl Topo {
    fn from_problem(problem: &Problem) -> Self {
        let topo = &problem.topology;
        let nn = topo.num_nodes;
        let ne = topo.num_edges;
        let mut starts = vec![0usize; ne];
        let mut ends = vec![0usize; ne];
        let inc = &topo.incidence;
        for col in 0..nn {
            for idx in inc.col_ptrs[col] as usize..inc.col_ptrs[col + 1] as usize {
                let row = inc.row_indices[idx] as usize;
                if inc.values[idx] < 0.0 {
                    starts[row] = col;
                } else {
                    ends[row] = col;
                }
            }
        }
        let mut node_to_free = vec![None; nn];
        for (i, &node) in topo.free_node_indices.iter().enumerate() {
            node_to_free[node] = Some(i);
        }
        Self {
            starts,
            ends,
            node_to_free,
            nfree: topo.free_node_indices.len(),
        }
    }
}

struct LevelWork {
    b: Vec<f64>,
    x: Vec<f64>,
    r: Vec<f64>,
    d: Vec<f64>,
    c1: Vec<f64>,
    v1: Vec<f64>,
    v2: Vec<f64>,
    rt: Vec<f64>,
    solve_work: Vec<f64>,
}

impl LevelWork {
    fn new(n: usize) -> Self {
        let z = || vec![0.0; n * 3];
        Self {
            b: z(),
            x: z(),
            r: z(),
            d: z(),
            c1: z(),
            v1: z(),
            v2: z(),
            rt: z(),
            solve_work: Vec::new(),
        }
    }
}

struct Coarsest {
    factorization: Factorization,
    stack: dyn_stack::GlobalPodBuffer,
}

struct Amg {
    params: Params,
    levels: Vec<Operator>,
    transfers: Vec<Transfer>,
    lambda_max: Vec<f64>,
    /// Level-0 graph edge → problem edge.
    graph_edge_of: Vec<u32>,
    /// `(free node, problem edge)` for edges with one fixed end.
    anchor_edges: Vec<(u32, u32)>,
    coarsest: Option<Coarsest>,
    work: Vec<LevelWork>,
    /// Operator applications in level-0 edge-equivalents (work counter).
    ops: Cell<f64>,
}

#[derive(Clone, Copy, Debug, Default)]
struct SolveStats {
    iterations: [usize; 3],
    relative: [f64; 3],
    ops: f64,
}

impl Amg {
    fn setup(topo: &Topo, q: &[f64], params: Params) -> Self {
        let nfree = topo.nfree;
        let mut eu = Vec::new();
        let mut ev = Vec::new();
        let mut weight = Vec::new();
        let mut graph_edge_of = Vec::new();
        let mut anchor = vec![0.0; nfree];
        let mut anchor_edges = Vec::new();
        for e in 0..q.len() {
            match (
                topo.node_to_free[topo.starts[e]],
                topo.node_to_free[topo.ends[e]],
            ) {
                (Some(u), Some(v)) => {
                    eu.push(u as u32);
                    ev.push(v as u32);
                    weight.push(q[e]);
                    graph_edge_of.push(e as u32);
                }
                (Some(u), None) | (None, Some(u)) => {
                    anchor[u] += q[e];
                    anchor_edges.push((u as u32, e as u32));
                }
                (None, None) => {}
            }
        }
        let fine = LevelGraph::new(nfree, eu, ev, weight, anchor);
        let mut levels = vec![Operator::Graph(fine)];
        let mut transfers = Vec::new();
        let mut lambda_max = Vec::new();
        let mut rng = Rng::new(0x5eed);

        while levels.last().unwrap().n() >= params.coarsest && levels.len() < 60 {
            let op = levels.last().unwrap();
            let n = op.n();
            // Aggregation uses the graph structure; for smoothed aggregation
            // on general sparse levels the strength graph is |a_uv|.
            let strength_graph = match op {
                Operator::Graph(_) => None,
                Operator::Sparse(s) => Some(sparse_to_strength_graph(s)),
            };
            let g = match op {
                Operator::Graph(g) => g,
                Operator::Sparse(_) => strength_graph.as_ref().unwrap(),
            };
            let (agg, nc) = aggregate(g, params.passes);
            if nc as f64 > 0.7 * n as f64 {
                break;
            }
            let lmax = estimate_lambda_max(op, params.power_iterations, &mut rng);
            lambda_max.push(lmax);
            if params.smoothed {
                let omega = 4.0 / (3.0 * lmax);
                let (p, coarse) = smoothed_transfer(op, &agg, nc, omega);
                transfers.push(Transfer::Smoothed {
                    aggregate_of: agg,
                    p,
                });
                levels.push(Operator::Sparse(coarse));
            } else {
                let (coarse, fine_offsets, fine_list) = coarsen(g, &agg, nc);
                transfers.push(Transfer::Aggregate {
                    aggregate_of: agg,
                    fine_offsets,
                    fine_list,
                });
                levels.push(Operator::Graph(coarse));
            }
        }
        let work = levels
            .iter()
            .map(|op| LevelWork::new(op.n()))
            .collect::<Vec<_>>();
        let mut amg = Self {
            params,
            levels,
            transfers,
            lambda_max,
            graph_edge_of,
            anchor_edges,
            coarsest: None,
            work,
            ops: Cell::new(0.0),
        };
        let nc = amg.levels.last().unwrap().n();
        amg.work.last_mut().unwrap().solve_work = vec![0.0; 6 * nc];
        amg.factor_coarsest();
        amg
    }

    fn factor_coarsest(&mut self) {
        let matrix = self.levels.last().unwrap().to_sparse();
        match &mut self.coarsest {
            Some(c) => c.factorization.update(&matrix, &mut c.stack).unwrap(),
            None => {
                let mut stack = dyn_stack::GlobalPodBuffer::new(dyn_stack::StackReq::empty());
                let factorization =
                    Factorization::new(&matrix, FactorizationStrategy::Cholesky, &mut stack)
                        .unwrap();
                self.coarsest = Some(Coarsest {
                    factorization,
                    stack,
                });
            }
        }
    }

    /// New force densities on the same topology: level weights, anchors,
    /// coarsest refactorization. The Chebyshev bounds are kept (they are
    /// upper bounds with a 10% margin; a large `q` drift would re-estimate).
    fn update(&mut self, q: &[f64]) {
        if self.params.smoothed {
            // General sparse levels have no cheap weight update; rebuild the
            // hierarchy with the same aggregation (prototype only).
            self.update_level0(q);
            for l in 0..self.transfers.len() {
                let omega = 4.0 / (3.0 * self.lambda_max[l]);
                let nc = self.levels[l + 1].n();
                let agg = self.transfers[l].aggregate_of().to_vec();
                let (p, coarse) = smoothed_transfer(&self.levels[l], &agg, nc, omega);
                self.transfers[l] = Transfer::Smoothed {
                    aggregate_of: agg,
                    p,
                };
                self.levels[l + 1] = Operator::Sparse(coarse);
            }
            self.factor_coarsest();
            return;
        }
        self.update_level0(q);
        for l in 0..self.transfers.len() {
            let (fine, coarse) = self.levels.split_at_mut(l + 1);
            let (Operator::Graph(fine), Operator::Graph(coarse)) = (&fine[l], &mut coarse[0])
            else {
                unreachable!()
            };
            let Transfer::Aggregate {
                aggregate_of,
                fine_offsets,
                fine_list,
            } = &self.transfers[l]
            else {
                unreachable!()
            };
            for j in 0..coarse.ne() {
                let mut w = 0.0;
                for &f in &fine_list[fine_offsets[j] as usize..fine_offsets[j + 1] as usize] {
                    w += fine.weight[f as usize];
                }
                coarse.weight[j] = w;
            }
            coarse.anchor.fill(0.0);
            for (u, &c) in aggregate_of.iter().enumerate() {
                coarse.anchor[c as usize] += fine.anchor[u];
            }
            coarse.refresh();
        }
        self.factor_coarsest();
    }

    fn update_level0(&mut self, q: &[f64]) {
        let Operator::Graph(g) = &mut self.levels[0] else {
            unreachable!()
        };
        for (j, &e) in self.graph_edge_of.iter().enumerate() {
            g.weight[j] = q[e as usize];
        }
        g.anchor.fill(0.0);
        for &(u, e) in &self.anchor_edges {
            g.anchor[u as usize] += q[e as usize];
        }
        g.refresh();
    }

    fn sizes(&self) -> Vec<usize> {
        self.levels.iter().map(|op| op.n()).collect()
    }

    /// Bytes held by the hierarchy above level 0 (coarse operators and
    /// transfers), for the memory budget.
    fn coarse_bytes(&self) -> usize {
        let mut total = 0;
        for op in &self.levels[1..] {
            total += match op {
                Operator::Graph(g) => {
                    g.eu.len() * 8 + g.weight.len() * 8 + g.anchor.len() * 8 + g.other.len() * 16
                }
                Operator::Sparse(s) => s.cols.len() * 12 + s.diag.len() * 8,
            };
        }
        for t in &self.transfers {
            total += match t {
                Transfer::Aggregate {
                    aggregate_of,
                    fine_list,
                    ..
                } => 4 * (aggregate_of.len() + 2 * fine_list.len()),
                Transfer::Smoothed { p, .. } => 12 * p.cols.len(),
            };
        }
        total
    }

    fn apply(&self, l: usize, x: &[f64], y: &mut [f64]) {
        let op = &self.levels[l];
        op.apply::<3>(x, y);
        self.ops.set(self.ops.get() + op.cost());
    }

    fn residual(&self, l: usize, b: &[f64], x: &[f64], r: &mut [f64]) {
        self.apply(l, x, r);
        for (ri, bi) in r.iter_mut().zip(b) {
            *ri = bi - *ri;
        }
    }

    /// Chebyshev polynomial smoother of fixed degree on `D⁻¹A` over
    /// `[λ_max/α, λ_max]` (three-term recurrence). With `zero_initial` the
    /// first residual is `b` and one operator application is saved; the
    /// polynomial is the same, so pre- and post-smoothing are adjoint.
    fn smooth(
        &self,
        l: usize,
        b: &[f64],
        x: &mut [f64],
        zero_initial: bool,
        r: &mut [f64],
        d: &mut [f64],
    ) {
        let diag = self.levels[l].diag();
        let lmax = self.lambda_max[l];
        let lmin = lmax / self.params.alpha;
        let theta = 0.5 * (lmax + lmin);
        let delta = 0.5 * (lmax - lmin);
        let sigma = theta / delta;
        let mut rho = 1.0 / sigma;
        if zero_initial {
            x.fill(0.0);
            for u in 0..diag.len() {
                let inv = 1.0 / diag[u];
                for k in 0..3 {
                    r[u * 3 + k] = b[u * 3 + k] * inv;
                }
            }
        } else {
            self.apply(l, x, r);
            for u in 0..diag.len() {
                let inv = 1.0 / diag[u];
                for k in 0..3 {
                    r[u * 3 + k] = (b[u * 3 + k] - r[u * 3 + k]) * inv;
                }
            }
        }
        let inv_theta = 1.0 / theta;
        for (di, ri) in d.iter_mut().zip(r.iter()) {
            *di = ri * inv_theta;
        }
        for (xi, di) in x.iter_mut().zip(d.iter()) {
            *xi += di;
        }
        for _ in 1..self.params.degree {
            self.apply(l, x, r);
            for u in 0..diag.len() {
                let inv = 1.0 / diag[u];
                for k in 0..3 {
                    r[u * 3 + k] = (b[u * 3 + k] - r[u * 3 + k]) * inv;
                }
            }
            let rho_new = 1.0 / (2.0 * sigma - rho);
            let c_d = rho_new * rho;
            let c_r = 2.0 * rho_new / delta;
            for ((di, ri), xi) in d.iter_mut().zip(r.iter()).zip(x.iter_mut()) {
                *di = c_d * *di + c_r * ri;
                *xi += *di;
            }
            rho = rho_new;
        }
    }

    fn coarsest_solve(&self, b: &[f64], x: &mut [f64], work: &mut [f64]) {
        let c = self.coarsest.as_ref().unwrap();
        c.factorization.solve_slices::<3>(b, x, work);
    }

    /// One multigrid cycle at level `l`: `work[0].x ≈ A_l⁻¹ work[0].b`.
    fn cycle(&self, l: usize, work: &mut [LevelWork]) {
        let (cur, rest) = work.split_first_mut().unwrap();
        if l + 1 == self.levels.len() {
            self.coarsest_solve(&cur.b, &mut cur.x, &mut cur.solve_work);
            return;
        }
        let mut x = std::mem::take(&mut cur.x);
        self.smooth(l, &cur.b, &mut x, true, &mut cur.r, &mut cur.d);
        self.residual(l, &cur.b, &x, &mut cur.r);
        let t = &self.transfers[l];
        t.restrict(&cur.r, &mut rest[0].b);
        match self.params.cycle {
            Cycle::V => self.cycle(l + 1, rest),
            Cycle::K => self.kcycle_solve(l + 1, rest),
        }
        if self.params.over != 1.0 {
            for v in rest[0].x.iter_mut() {
                *v *= self.params.over;
            }
        }
        t.prolong_add(&rest[0].x, &mut x);
        self.smooth(l, &cur.b, &mut x, false, &mut cur.r, &mut cur.d);
        cur.x = x;
    }

    /// Notay's K-cycle coarse solve at level `l`: at most two flexible-CG
    /// steps on `A_l x = b` preconditioned by `cycle(l)`, the second one
    /// skipped when the first reduces the residual below `ktol`.
    fn kcycle_solve(&self, l: usize, work: &mut [LevelWork]) {
        if l + 1 == self.levels.len() {
            self.cycle(l, work);
            return;
        }
        let b_norm = norm3(&work[0].b);
        self.cycle(l, work);
        let (rho1, alpha1, a);
        {
            let w = &mut work[0];
            w.c1.copy_from_slice(&w.x);
            self.apply(l, &w.c1, &mut w.v1);
            rho1 = dot3(&w.c1, &w.v1);
            alpha1 = dot3(&w.c1, &w.b);
            a = div3(alpha1, rho1);
            for ((rt, bi), v) in
                w.rt.chunks_exact_mut(3)
                    .zip(w.b.chunks_exact(3))
                    .zip(w.v1.chunks_exact(3))
            {
                for k in 0..3 {
                    rt[k] = bi[k] - a[k] * v[k];
                }
            }
            let rt_norm = norm3(&w.rt);
            if (0..3).all(|k| rt_norm[k] <= self.params.ktol * b_norm[k]) {
                for (xi, ci) in w.x.chunks_exact_mut(3).zip(w.c1.chunks_exact(3)) {
                    for k in 0..3 {
                        xi[k] = a[k] * ci[k];
                    }
                }
                return;
            }
            w.b.copy_from_slice(&w.rt);
        }
        self.cycle(l, work);
        let w = &mut work[0];
        self.apply(l, &w.x, &mut w.v2);
        let gamma = dot3(&w.x, &w.v1);
        let beta = dot3(&w.x, &w.rt);
        let c2v2 = dot3(&w.x, &w.v2);
        let mut coef1 = [0.0; 3];
        let mut coef2 = [0.0; 3];
        for k in 0..3 {
            let rho2 = c2v2[k] - safe_div(gamma[k] * gamma[k], rho1[k]);
            coef1[k] = safe_div(alpha1[k], rho1[k]) - safe_div(gamma[k] * beta[k], rho1[k] * rho2);
            coef2[k] = safe_div(beta[k], rho2);
        }
        for (xi, ci) in w.x.chunks_exact_mut(3).zip(w.c1.chunks_exact(3)) {
            for k in 0..3 {
                xi[k] = coef1[k] * ci[k] + coef2[k] * xi[k];
            }
        }
    }

    /// `z = M r` with one cycle at level 0.
    fn precondition(&self, r: &[f64], z: &mut [f64], work: &mut [LevelWork]) {
        work[0].b.copy_from_slice(r);
        self.cycle(0, work);
        z.copy_from_slice(&work[0].x);
    }

    /// PCG (V-cycle) or FCG(1) (K-cycle) on three right-hand sides with
    /// per-column convergence to `‖r‖ ≤ tol ‖b‖`. `x` holds the initial
    /// guess on entry.
    fn solve(&mut self, b: &[f64], x: &mut [f64], tol: f64, max_iterations: usize) -> SolveStats {
        let flexible = self.params.cycle == Cycle::K;
        let n3 = b.len();
        let mut work = std::mem::take(&mut self.work);
        let ops0 = self.ops.get();
        let mut r = vec![0.0; n3];
        let mut z = vec![0.0; n3];
        let mut p = vec![0.0; n3];
        let mut ap = vec![0.0; n3];
        self.residual(0, b, x, &mut r);
        let b_norm = norm3(b);
        let mut converged = [false; 3];
        let mut iterations = [0usize; 3];
        for k in 0..3 {
            if b_norm[k] == 0.0 {
                for xi in x.chunks_exact_mut(3) {
                    xi[k] = 0.0;
                }
                for ri in r.chunks_exact_mut(3) {
                    ri[k] = 0.0;
                }
                converged[k] = true;
            }
        }
        let floor = b_norm.map(|v| v.max(1e-300));
        let mut rz_prev = [0.0; 3];
        let mut it = 0;
        loop {
            let rn = norm3(&r);
            for k in 0..3 {
                if !converged[k] && rn[k] <= tol * floor[k] {
                    converged[k] = true;
                    iterations[k] = it;
                }
            }
            if converged.iter().all(|&c| c) || it >= max_iterations {
                break;
            }
            self.precondition(&r, &mut z, &mut work);
            let rz = dot3(&r, &z);
            if it == 0 {
                p.copy_from_slice(&z);
            } else if flexible {
                let beta = div3(dot3(&z, &ap), dot3(&p, &ap));
                for (pi, zi) in p.chunks_exact_mut(3).zip(z.chunks_exact(3)) {
                    for k in 0..3 {
                        pi[k] = zi[k] - beta[k] * pi[k];
                    }
                }
            } else {
                let beta = div3(rz, rz_prev);
                for (pi, zi) in p.chunks_exact_mut(3).zip(z.chunks_exact(3)) {
                    for k in 0..3 {
                        pi[k] = zi[k] + beta[k] * pi[k];
                    }
                }
            }
            rz_prev = rz;
            self.apply(0, &p, &mut ap);
            let pap = dot3(&p, &ap);
            let num = if flexible { dot3(&p, &r) } else { rz };
            let mut alpha = div3(num, pap);
            for k in 0..3 {
                if converged[k] {
                    alpha[k] = 0.0;
                }
            }
            axpy3(x, alpha, &p);
            axpy3(&mut r, alpha.map(|a| -a), &ap);
            it += 1;
        }
        for k in 0..3 {
            if !converged[k] {
                iterations[k] = it;
            }
        }
        // True residual, recomputed from scratch.
        self.residual(0, b, x, &mut r);
        let rn = norm3(&r);
        let ops = (self.ops.get() - ops0) / self.levels[0].cost();
        self.work = work;
        SolveStats {
            iterations,
            relative: [
                safe_div(rn[0], floor[0]),
                safe_div(rn[1], floor[1]),
                safe_div(rn[2], floor[2]),
            ],
            ops,
        }
    }
}

/// Strength graph of a general sparse level: `w_uv = |a_uv|`, anchor so
/// that `diag` is preserved (aggregation only; never applied).
fn sparse_to_strength_graph(s: &SparseLevel) -> LevelGraph {
    let mut eu = Vec::new();
    let mut ev = Vec::new();
    let mut weight = Vec::new();
    let mut offdiag = vec![0.0; s.n];
    for u in 0..s.n {
        for j in s.offsets[u] as usize..s.offsets[u + 1] as usize {
            let v = s.cols[j] as usize;
            let w = s.vals[j].abs();
            offdiag[u] += w;
            if u < v {
                eu.push(u as u32);
                ev.push(v as u32);
                weight.push(w);
            }
        }
    }
    let anchor: Vec<f64> = (0..s.n)
        .map(|u| (s.diag[u] - offdiag[u]).max(0.0))
        .collect();
    LevelGraph::new(s.n, eu, ev, weight, anchor)
}

/// Power iteration for `λ_max(D⁻¹A)` with a 10% safety margin.
fn estimate_lambda_max(op: &Operator, iterations: usize, rng: &mut Rng) -> f64 {
    let n = op.n();
    let diag = op.diag();
    let mut v: Vec<f64> = (0..n).map(|_| rng.symmetric()).collect();
    let mut w = vec![0.0; n];
    let mut lambda = 1.0;
    let nv = dot1(&v, &v).sqrt();
    for vi in v.iter_mut() {
        *vi /= nv;
    }
    for _ in 0..iterations {
        op.apply::<1>(&v, &mut w);
        for u in 0..n {
            w[u] /= diag[u];
        }
        lambda = dot1(&w, &w).sqrt();
        for u in 0..n {
            v[u] = w[u] / lambda;
        }
    }
    lambda * 1.1
}

// ─────────────────────────────────────────────────────────────
//  Fixtures
// ─────────────────────────────────────────────────────────────

struct Fixture {
    name: String,
    problem: Problem,
    /// Planar coordinates of every node (for the smooth `q` field).
    xy: Vec<[f64; 2]>,
}

/// `Problem` with unit downward loads on free nodes, all supports fixed,
/// box bounds `[0.1, 10]` (Cholesky path) and no objectives.
fn build_problem(
    num_nodes: usize,
    edges: &[(usize, usize)],
    fixed: &[usize],
    positions: &[[f64; 3]],
) -> Problem {
    let ne = edges.len();
    let mut rows = Vec::with_capacity(2 * ne);
    let mut cols = Vec::with_capacity(2 * ne);
    let mut vals = Vec::with_capacity(2 * ne);
    for (e, &(s, t)) in edges.iter().enumerate() {
        rows.extend([e, e]);
        cols.extend([s, t]);
        vals.extend([-1.0, 1.0]);
    }
    let incidence = SparseColMatOwned::from_coo(ne, num_nodes, &rows, &cols, &vals).unwrap();
    let mut fixed_idx = fixed.to_vec();
    fixed_idx.sort_unstable();
    fixed_idx.dedup();
    let mut is_fixed = vec![false; num_nodes];
    for &f in &fixed_idx {
        is_fixed[f] = true;
    }
    let free_idx: Vec<usize> = (0..num_nodes).filter(|&i| !is_fixed[i]).collect();
    let nn_free = free_idx.len();
    let topology = NetworkTopology {
        free_incidence: incidence.extract_columns(&free_idx),
        fixed_incidence: incidence.extract_columns(&fixed_idx),
        incidence,
        num_edges: ne,
        num_nodes,
        free_node_indices: free_idx,
        fixed_node_indices: fixed_idx.clone(),
    };
    let mut loads = vec![0.0; nn_free * 3];
    for i in 0..nn_free {
        loads[i * 3 + 2] = -1.0;
    }
    let mut fixed_pos = Array2::zeros((fixed_idx.len(), 3));
    for (i, &f) in fixed_idx.iter().enumerate() {
        for d in 0..3 {
            fixed_pos[[i, d]] = positions[f][d];
        }
    }
    Problem {
        topology,
        free_node_loads: Array2::from_shape_vec((nn_free, 3), loads).unwrap(),
        anchors: AnchorInfo::all_fixed(fixed_pos.clone()),
        fixed_node_positions: fixed_pos,
        objectives: Vec::new(),
        bounds: Bounds {
            lower: vec![0.1; ne],
            upper: vec![10.0; ne],
        },
        solver: SolverOptions::default(),
        self_weight: None,
        pressure: None,
    }
}

fn make_grid_fixture(n: usize) -> Fixture {
    let problem = grid::make_grid_problem(n);
    let xy = (0..n * n)
        .map(|i| [(i % n) as f64, (i / n) as f64])
        .collect();
    Fixture {
        name: format!("grid {n}"),
        problem,
        xy,
    }
}

/// Irregular triangulated mesh: an `m × m` lattice jittered by up to 0.4 of
/// the spacing, each point connected to its 5–7 nearest neighbours, edge set
/// symmetrised, boundary ring fixed.
fn make_irregular_fixture(m: usize) -> Fixture {
    let mut rng = Rng::new(0x1a2b3c);
    let nn = m * m;
    let mut xy = Vec::with_capacity(nn);
    for j in 0..m {
        for i in 0..m {
            xy.push([
                i as f64 + 0.4 * rng.symmetric(),
                j as f64 + 0.4 * rng.symmetric(),
            ]);
        }
    }
    // Bucket points by lattice cell; the k nearest of any point lie within
    // two cells (density one point per unit area, jitter ≤ 0.4).
    let cell = |p: [f64; 2]| -> (i64, i64) { (p[0].round() as i64, p[1].round() as i64) };
    let mut buckets: Vec<Vec<u32>> = vec![Vec::new(); nn];
    let bucket_index = |c: (i64, i64)| -> Option<usize> {
        if c.0 < 0 || c.1 < 0 || c.0 >= m as i64 || c.1 >= m as i64 {
            None
        } else {
            Some(c.1 as usize * m + c.0 as usize)
        }
    };
    for (i, &p) in xy.iter().enumerate() {
        buckets[bucket_index(cell(p)).unwrap()].push(i as u32);
    }
    let mut pairs: Vec<(u32, u32)> = Vec::with_capacity(nn * 4);
    let mut cand: Vec<(f64, u32)> = Vec::new();
    for (u, &p) in xy.iter().enumerate() {
        let k = 5 + (rng.next_u64() % 3) as usize;
        let c = cell(p);
        cand.clear();
        for dj in -2..=2 {
            for di in -2..=2 {
                if let Some(bi) = bucket_index((c.0 + di, c.1 + dj)) {
                    for &v in &buckets[bi] {
                        if v as usize != u {
                            let q = xy[v as usize];
                            let d2 = (q[0] - p[0]).powi(2) + (q[1] - p[1]).powi(2);
                            cand.push((d2, v));
                        }
                    }
                }
            }
        }
        cand.sort_by(|a, b| a.partial_cmp(b).unwrap());
        for &(_, v) in cand.iter().take(k) {
            let (a, b) = (u as u32, v);
            pairs.push(if a < b { (a, b) } else { (b, a) });
        }
    }
    pairs.sort_unstable();
    pairs.dedup();
    let edges: Vec<(usize, usize)> = pairs
        .iter()
        .map(|&(a, b)| (a as usize, b as usize))
        .collect();
    let fixed: Vec<usize> = (0..nn)
        .filter(|&i| {
            let (x, y) = (i % m, i / m);
            x == 0 || y == 0 || x == m - 1 || y == m - 1
        })
        .collect();
    let positions: Vec<[f64; 3]> = xy.iter().map(|p| [p[0], p[1], 0.0]).collect();
    Fixture {
        name: format!("irregular {m}"),
        problem: build_problem(nn, &edges, &fixed, &positions),
        xy,
    }
}

/// Cable dome / spoke wheel: a hub, `rings` concentric rings whose node
/// count doubles whenever the circumferential spacing exceeds 1.5× the
/// radial spacing, radial spokes, diagonal bracing, outer ring fixed. Hub
/// degree 12; nodes on a doubling ring have degree 8.
fn make_dome_fixture(rings: usize) -> Fixture {
    let mut counts = vec![0usize; rings + 1];
    counts[1] = 12;
    for r in 1..rings {
        let spacing = std::f64::consts::TAU * (r + 1) as f64 / counts[r] as f64;
        counts[r + 1] = if spacing > 1.5 {
            2 * counts[r]
        } else {
            counts[r]
        };
    }
    let mut first = vec![0usize; rings + 2];
    first[1] = 1;
    for r in 1..=rings {
        first[r + 1] = first[r] + counts[r];
    }
    let nn = first[rings + 1];
    let mut xy = vec![[0.0, 0.0]; nn];
    let node = |r: usize, s: usize| first[r] + s.rem_euclid(counts[r]);
    for r in 1..=rings {
        for s in 0..counts[r] {
            let angle = std::f64::consts::TAU * s as f64 / counts[r] as f64;
            xy[node(r, s)] = [r as f64 * angle.cos(), r as f64 * angle.sin()];
        }
    }
    let mut edges = Vec::with_capacity(3 * nn);
    for s in 0..counts[1] {
        edges.push((0, node(1, s)));
    }
    for r in 1..=rings {
        for s in 0..counts[r] {
            edges.push((node(r, s), node(r, s + 1)));
        }
        if r < rings {
            for s in 0..counts[r] {
                let u = node(r, s);
                if counts[r + 1] == counts[r] {
                    edges.push((u, node(r + 1, s)));
                    edges.push((u, node(r + 1, s + 1)));
                } else {
                    edges.push((u, node(r + 1, 2 * s)));
                    edges.push((u, node(r + 1, 2 * s + 1)));
                    edges.push((u, node(r + 1, 2 * s + counts[r + 1] - 1)));
                    edges.push((u, node(r + 1, 2 * s + 2)));
                }
            }
        }
    }
    let fixed: Vec<usize> = (0..counts[rings]).map(|s| node(rings, s)).collect();
    let h = 0.25 * rings as f64;
    let positions: Vec<[f64; 3]> = xy
        .iter()
        .map(|p| {
            let rr = (p[0] * p[0] + p[1] * p[1]).sqrt() / rings as f64;
            [p[0], p[1], h * (1.0 - rr * rr)]
        })
        .collect();
    Fixture {
        name: format!("dome {rings}"),
        problem: build_problem(nn, &edges, &fixed, &positions),
        xy,
    }
}

/// Fixture whose edge count matches the `size × size` grid.
fn make_fixture(kind: &str, size: usize) -> Fixture {
    let target = grid::grid_edges(size) as f64;
    match kind {
        "grid" => make_grid_fixture(size),
        "irregular" => {
            let mut m = ((target / 3.7).sqrt() as usize).max(6);
            let f = make_irregular_fixture(m);
            let ratio = target / f.problem.topology.num_edges as f64;
            m = ((m as f64 * ratio.sqrt()).round() as usize).max(6);
            make_irregular_fixture(m)
        }
        "dome" => {
            let mut r = ((target / 8.0).sqrt() as usize).max(3);
            let f = make_dome_fixture(r);
            let ratio = target / f.problem.topology.num_edges as f64;
            r = ((r as f64 * ratio.sqrt()).round() as usize).max(3);
            make_dome_fixture(r)
        }
        other => panic!("unknown fixture {other}"),
    }
}

/// Smooth force-density field `q ∈ [1, ratio]` varying over the domain.
fn q_field(fx: &Fixture, topo: &Topo, ratio: f64) -> Vec<f64> {
    let (mut lo, mut hi) = ([f64::MAX; 2], [f64::MIN; 2]);
    for p in &fx.xy {
        for d in 0..2 {
            lo[d] = lo[d].min(p[d]);
            hi[d] = hi[d].max(p[d]);
        }
    }
    let span = [(hi[0] - lo[0]).max(1e-12), (hi[1] - lo[1]).max(1e-12)];
    let ln = ratio.ln();
    (0..topo.starts.len())
        .map(|e| {
            let a = fx.xy[topo.starts[e]];
            let b = fx.xy[topo.ends[e]];
            let x = (0.5 * (a[0] + b[0]) - lo[0]) / span[0];
            let y = (0.5 * (a[1] + b[1]) - lo[1]) / span[1];
            let f =
                0.5 + 0.5 * (std::f64::consts::TAU * x).sin() * (std::f64::consts::TAU * y).cos();
            (ln * f).exp()
        })
        .collect()
}

/// `b = p + Σ_{edges to fixed nodes} q_e x_fixed`.
fn rhs(fx: &Fixture, topo: &Topo, q: &[f64]) -> Vec<f64> {
    let problem = &fx.problem;
    let mut b: Vec<f64> = problem.free_node_loads.iter().copied().collect();
    let mut fixed_row = vec![usize::MAX; problem.topology.num_nodes];
    for (i, &f) in problem.topology.fixed_node_indices.iter().enumerate() {
        fixed_row[f] = i;
    }
    for e in 0..q.len() {
        let (s, t) = (topo.starts[e], topo.ends[e]);
        let (free, fixed) = match (topo.node_to_free[s], topo.node_to_free[t]) {
            (Some(u), None) => (u, t),
            (None, Some(u)) => (u, s),
            _ => continue,
        };
        let row = fixed_row[fixed];
        for d in 0..3 {
            b[free * 3 + d] += q[e] * problem.fixed_node_positions[[row, d]];
        }
    }
    b
}

fn max_abs_diff(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b)
        .map(|(x, y)| (x - y).abs())
        .fold(0.0, f64::max)
}

// ─────────────────────────────────────────────────────────────
//  Self-test
// ─────────────────────────────────────────────────────────────

fn dense_ptap(a: &[Vec<f64>], p: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let n = a.len();
    let nc = p[0].len();
    let mut ap = vec![vec![0.0; nc]; n];
    for i in 0..n {
        for k in 0..n {
            let aik = a[i][k];
            if aik != 0.0 {
                for j in 0..nc {
                    ap[i][j] += aik * p[k][j];
                }
            }
        }
    }
    let mut out = vec![vec![0.0; nc]; nc];
    for i in 0..nc {
        for k in 0..n {
            let pki = p[k][i];
            if pki != 0.0 {
                for j in 0..nc {
                    out[i][j] += pki * ap[k][j];
                }
            }
        }
    }
    out
}

/// Random connected 200-node graph with a few fixed nodes, as a `Topo`.
fn random_topo(rng: &mut Rng, nfree: usize, nfixed: usize) -> (Topo, Vec<f64>) {
    let nn = nfree + nfixed;
    let mut edges: Vec<(usize, usize)> = Vec::new();
    for i in 0..nn - 1 {
        edges.push((i, i + 1));
    }
    for u in 0..nn {
        for _ in 0..3 {
            let v = rng.below(nn);
            if v != u {
                edges.push((u, v));
            }
        }
    }
    let fixed: Vec<usize> = (0..nfixed).map(|k| (k * 37 + 11) % nn).collect();
    let mut node_to_free = vec![None; nn];
    let mut next = 0;
    for i in 0..nn {
        if !fixed.contains(&i) {
            node_to_free[i] = Some(next);
            next += 1;
        }
    }
    let q: Vec<f64> = (0..edges.len())
        .map(|_| 0.1 * (rng.uniform() * 100.0f64.ln()).exp())
        .collect();
    (
        Topo {
            starts: edges.iter().map(|e| e.0).collect(),
            ends: edges.iter().map(|e| e.1).collect(),
            node_to_free,
            nfree: next,
        },
        q,
    )
}

fn self_test() {
    let mut rng = Rng::new(42);
    let mut failures = 0;
    let mut check = |name: &str, ok: bool, detail: String| {
        println!("{} {name}: {detail}", if ok { "ok  " } else { "FAIL" });
        if !ok {
            failures += 1;
        }
    };
    for smoothed in [false, true] {
        let (topo, q) = random_topo(&mut rng, 200, 7);
        let params = Params {
            cycle: Cycle::V,
            degree: 2,
            passes: 2,
            coarsest: 12,
            alpha: 30.0,
            ktol: 0.25,
            smoothed,
            over: 1.0,
            power_iterations: 10,
        };
        let mut amg = Amg::setup(&topo, &q, params.clone());
        let tag = if smoothed { "smoothed" } else { "plain" };
        println!("[{tag}] levels {:?}", amg.sizes());
        // 1. Coarse operator equals Pᵀ A P on every level.
        let mut worst = 0.0;
        for l in 0..amg.transfers.len() {
            let a = amg.levels[l].dense();
            let p = amg.transfers[l].dense(amg.levels[l].n(), amg.levels[l + 1].n());
            let ptap = dense_ptap(&a, &p);
            let ac = amg.levels[l + 1].dense();
            let scale = ac.iter().flatten().fold(0.0f64, |m, v| m.max(v.abs()));
            for i in 0..ac.len() {
                for j in 0..ac.len() {
                    worst = f64::max(worst, (ptap[i][j] - ac[i][j]).abs() / scale);
                }
            }
        }
        check(
            &format!("[{tag}] coarse operator = PᵀAP"),
            worst < 1e-12,
            format!("max relative entry error {worst:.2e}"),
        );
        // 2. V-cycle preconditioner symmetric and positive definite.
        let n3 = topo.nfree * 3;
        for degree in [1, 2, 3] {
            amg.params.degree = degree;
            let u: Vec<f64> = (0..n3).map(|_| rng.symmetric()).collect();
            let v: Vec<f64> = (0..n3).map(|_| rng.symmetric()).collect();
            let mut mu = vec![0.0; n3];
            let mut mv = vec![0.0; n3];
            let mut work = std::mem::take(&mut amg.work);
            amg.precondition(&u, &mut mu, &mut work);
            amg.precondition(&v, &mut mv, &mut work);
            amg.work = work;
            let muv = dot1(&mu, &v);
            let umv = dot1(&u, &mv);
            let muu = dot1(&mu, &u);
            let asym = (muv - umv).abs() / muv.abs().max(umv.abs());
            check(
                &format!("[{tag}] V-cycle symmetric, degree {degree}"),
                asym < 1e-10,
                format!("|<Mu,v> − <u,Mv>| / |<Mu,v>| = {asym:.2e}"),
            );
            check(
                &format!("[{tag}] V-cycle positive, degree {degree}"),
                muu > 0.0,
                format!("<Mu,u> = {muu:.3e}"),
            );
        }
        amg.params.degree = 2;
        // 3. PCG (V) and FCG (K) solve to the direct solution.
        let b: Vec<f64> = (0..n3).map(|_| rng.symmetric()).collect();
        let matrix = amg.levels[0].to_sparse();
        let mut stack = dyn_stack::GlobalPodBuffer::new(dyn_stack::StackReq::empty());
        let fac = Factorization::new(&matrix, FactorizationStrategy::Cholesky, &mut stack).unwrap();
        let mut x_ref = vec![0.0; n3];
        let mut w = vec![0.0; 2 * n3];
        fac.solve_slices::<3>(&b, &mut x_ref, &mut w);
        let x_scale = x_ref.iter().fold(0.0f64, |m, v| m.max(v.abs()));
        for cycle in [Cycle::V, Cycle::K] {
            amg.params.cycle = cycle;
            let mut x = vec![0.0; n3];
            let stats = amg.solve(&b, &mut x, 1e-10, 200);
            let err = max_abs_diff(&x, &x_ref) / x_scale;
            check(
                &format!("[{tag}] {:?}-cycle solve matches direct", cycle),
                err < 1e-8 && stats.iterations.iter().all(|&i| i < 200),
                format!(
                    "iterations {:?}, relative residual {:.1e}, max |x − x_direct| / max |x| = {err:.1e}",
                    stats.iterations,
                    stats.relative.iter().cloned().fold(0.0, f64::max)
                ),
            );
        }
        // 4. update(q') keeps every coarse operator equal to Pᵀ A(q') P.
        let q2: Vec<f64> = q
            .iter()
            .map(|v| v * (1.0 + 0.3 * rng.symmetric()))
            .collect();
        amg.update(&q2);
        let mut worst = 0.0;
        for l in 0..amg.transfers.len() {
            let a = amg.levels[l].dense();
            let p = amg.transfers[l].dense(amg.levels[l].n(), amg.levels[l + 1].n());
            let ptap = dense_ptap(&a, &p);
            let ac = amg.levels[l + 1].dense();
            let scale = ac.iter().flatten().fold(0.0f64, |m, v| m.max(v.abs()));
            for i in 0..ac.len() {
                for j in 0..ac.len() {
                    worst = f64::max(worst, (ptap[i][j] - ac[i][j]).abs() / scale);
                }
            }
        }
        let fresh = Amg::setup(&topo, &q2, params.clone());
        let a0 = amg.levels[0].dense();
        let a0_fresh = fresh.levels[0].dense();
        let scale = a0.iter().flatten().fold(0.0f64, |m, v| m.max(v.abs()));
        let mut diff = 0.0f64;
        for i in 0..a0.len() {
            for j in 0..a0.len() {
                diff = diff.max((a0[i][j] - a0_fresh[i][j]).abs() / scale);
            }
        }
        check(
            &format!("[{tag}] update(q'): level 0 = A(q'), coarse levels = PᵀAP"),
            worst < 1e-12 && diff < 1e-12,
            format!("max relative entry error {worst:.2e} (coarse), {diff:.2e} (fine)"),
        );
    }
    if failures > 0 {
        eprintln!("{failures} self-test failure(s)");
        std::process::exit(1);
    }
    println!("self-test passed");
}

// ─────────────────────────────────────────────────────────────
//  Command line and measurement loop
// ─────────────────────────────────────────────────────────────

struct Cli {
    fixture: String,
    size: usize,
    cycles: Vec<Cycle>,
    degrees: Vec<usize>,
    passes: Vec<usize>,
    coarsest: Vec<usize>,
    alphas: Vec<f64>,
    ktol: Vec<f64>,
    smoothed: Vec<bool>,
    over: Vec<f64>,
    warm: f64,
    reps: usize,
    qratio: f64,
    tol: f64,
    max_iterations: usize,
    direct: bool,
    self_test: bool,
}

fn parse_list<T: std::str::FromStr>(s: &str) -> Vec<T>
where
    T::Err: std::fmt::Debug,
{
    s.split(',')
        .map(|v| v.trim().parse().expect("list value"))
        .collect()
}

fn parse_cli() -> Cli {
    let mut cli = Cli {
        fixture: "grid".into(),
        size: 64,
        cycles: vec![Cycle::K],
        degrees: vec![2],
        passes: vec![2],
        coarsest: vec![2000],
        alphas: vec![30.0],
        ktol: vec![0.25],
        smoothed: vec![false],
        over: vec![1.0],
        warm: 2.0,
        reps: 1,
        qratio: 10.0,
        tol: 1e-8,
        max_iterations: 500,
        direct: true,
        self_test: false,
    };
    let args: Vec<String> = std::env::args().skip(1).collect();
    let mut i = 0;
    while i < args.len() {
        let flag = args[i].as_str();
        let mut value = || {
            i += 1;
            args.get(i)
                .cloned()
                .unwrap_or_else(|| panic!("{flag} needs a value"))
        };
        match flag {
            "--fixture" => cli.fixture = value(),
            "--size" => cli.size = value().parse().unwrap(),
            "--cycle" => {
                cli.cycles = value()
                    .split(',')
                    .map(|c| match c.trim() {
                        "v" | "V" => Cycle::V,
                        "k" | "K" => Cycle::K,
                        other => panic!("cycle {other}"),
                    })
                    .collect()
            }
            "--degree" => cli.degrees = parse_list(&value()),
            "--passes" => cli.passes = parse_list(&value()),
            "--coarsest" => cli.coarsest = parse_list(&value()),
            "--alpha" => cli.alphas = parse_list(&value()),
            "--ktol" => cli.ktol = parse_list(&value()),
            "--smoothed" => {
                cli.smoothed = value()
                    .split(',')
                    .map(|v| match v.trim() {
                        "0" | "false" | "plain" => false,
                        "1" | "true" | "sa" => true,
                        other => panic!("smoothed {other}"),
                    })
                    .collect()
            }
            "--over" => cli.over = parse_list(&value()),
            "--warm" => cli.warm = value().parse().unwrap(),
            "--reps" => cli.reps = value().parse().unwrap(),
            "--qratio" => cli.qratio = value().parse().unwrap(),
            "--tol" => cli.tol = value().parse().unwrap(),
            "--maxit" => cli.max_iterations = value().parse().unwrap(),
            "--no-direct" => cli.direct = false,
            "--self-test" => cli.self_test = true,
            other => panic!("unknown argument {other}"),
        }
        i += 1;
    }
    cli
}

fn fmt_iters(it: [usize; 3]) -> String {
    format!("{}/{}/{}", it[0], it[1], it[2])
}

fn fmt_sizes(sizes: &[usize]) -> String {
    sizes
        .iter()
        .map(|&n| {
            if n >= 10_000 {
                format!("{}k", (n as f64 / 1000.0).round() as usize)
            } else {
                n.to_string()
            }
        })
        .collect::<Vec<_>>()
        .join("/")
}

fn main() {
    let cli = parse_cli();
    if cli.self_test {
        self_test();
        return;
    }
    let t = Instant::now();
    let fx = make_fixture(&cli.fixture, cli.size);
    let topo = Topo::from_problem(&fx.problem);
    let ne = fx.problem.topology.num_edges;
    let nfree = topo.nfree;
    let q = q_field(&fx, &topo, cli.qratio);
    let b = rhs(&fx, &topo, &q);
    let mut rng = Rng::new(777);
    let q_warm: Vec<f64> = q
        .iter()
        .map(|v| v * (1.0 + 0.01 * cli.warm * rng.symmetric()))
        .collect();
    let b_warm = rhs(&fx, &topo, &q_warm);
    let mut degree_hist = vec![0usize; nfree];
    for e in 0..ne {
        for node in [topo.starts[e], topo.ends[e]] {
            if let Some(u) = topo.node_to_free[node] {
                degree_hist[u] += 1;
            }
        }
    }
    let max_deg = degree_hist.iter().copied().max().unwrap_or(0);
    let mean_deg = degree_hist.iter().sum::<usize>() as f64 / nfree.max(1) as f64;
    println!(
        "{}: {ne} edges, {nfree} free nodes, {} fixed, degree mean {mean_deg:.2} max {max_deg}, q ∈ [{:.2}, {:.2}] (built in {:.0} ms)",
        fx.name,
        fx.problem.topology.fixed_node_indices.len(),
        q.iter().cloned().fold(f64::MAX, f64::min),
        q.iter().cloned().fold(f64::MIN, f64::max),
        ms(t)
    );

    // Direct reference: FdmCache::new + solve_fdm (symbolic + first factor),
    // then time refactor + one 3-RHS solve.
    let mut x_ref: Option<Vec<f64>> = None;
    let mut x_ref_warm: Option<Vec<f64>> = None;
    let mut direct_ms = f64::NAN;
    if cli.direct {
        let anchors = Array2::zeros((0, 3));
        let t = Instant::now();
        let mut cache = FdmCache::new(&fx.problem).unwrap();
        theseus::fdm::solve_fdm(&mut cache, &q, &fx.problem, &anchors, 0.0).unwrap();
        let first_ms = ms(t);
        let cache_rhs: Vec<f64> = cache.rhs.iter().copied().collect();
        let rhs_diff = max_abs_diff(&cache_rhs, &b);
        let mut times = Vec::new();
        for _ in 0..cli.reps.max(1) {
            let t = Instant::now();
            theseus::fdm::solve_fdm(&mut cache, &q, &fx.problem, &anchors, 0.0).unwrap();
            times.push(ms(t));
        }
        direct_ms = median(&mut times);
        x_ref = Some(cache.x.iter().copied().collect());
        theseus::fdm::solve_fdm(&mut cache, &q_warm, &fx.problem, &anchors, 0.0).unwrap();
        x_ref_warm = Some(cache.x.iter().copied().collect());
        println!(
            "direct: setup + symbolic + first factor {first_ms:.0} ms, refactor + solve (3 rhs) {direct_ms:.1} ms (median of {}), |rhs − rhs_direct| = {rhs_diff:.1e}",
            cli.reps.max(1)
        );
    }

    for &smoothed in &cli.smoothed {
        for &over in &cli.over {
            for &passes in &cli.passes {
                for &coarsest in &cli.coarsest {
                    for &alpha in &cli.alphas {
                        for &cycle in &cli.cycles {
                            for &degree in &cli.degrees {
                                let ktols: Vec<f64> = if cycle == Cycle::K {
                                    cli.ktol.clone()
                                } else {
                                    vec![0.0]
                                };
                                for &ktol in &ktols {
                                    let params = Params {
                                        cycle,
                                        degree,
                                        passes,
                                        coarsest,
                                        alpha,
                                        ktol,
                                        smoothed,
                                        over,
                                        power_iterations: 10,
                                    };
                                    run_one(
                                        &cli,
                                        &fx,
                                        &topo,
                                        &q,
                                        &b,
                                        &q_warm,
                                        &b_warm,
                                        x_ref.as_deref(),
                                        x_ref_warm.as_deref(),
                                        direct_ms,
                                        params,
                                    );
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn run_one(
    cli: &Cli,
    fx: &Fixture,
    topo: &Topo,
    q: &[f64],
    b: &[f64],
    q_warm: &[f64],
    b_warm: &[f64],
    x_ref: Option<&[f64]>,
    x_ref_warm: Option<&[f64]>,
    direct_ms: f64,
    params: Params,
) {
    let ne = fx.problem.topology.num_edges;
    let nfree = topo.nfree;
    let n3 = nfree * 3;
    let t = Instant::now();
    let mut amg = Amg::setup(topo, q, params.clone());
    let setup_ms = ms(t);
    let sizes = amg.sizes();
    let coarse_edges: usize = amg.levels[1..].iter().map(|op| op.cost() as usize).sum();

    let mut x = vec![0.0; n3];
    let apply_ms = {
        let mut y = vec![0.0; n3];
        let t = Instant::now();
        for _ in 0..5 {
            amg.levels[0].apply::<3>(b, &mut y);
        }
        ms(t) / 5.0
    };
    let mut cold_times = Vec::new();
    let mut cold = SolveStats::default();
    for _ in 0..cli.reps.max(1) {
        x.fill(0.0);
        let t = Instant::now();
        cold = amg.solve(b, &mut x, cli.tol, cli.max_iterations);
        cold_times.push(ms(t));
    }
    let cold_ms = median(&mut cold_times);
    let x_scale = x.iter().fold(0.0f64, |m, v| m.max(v.abs())).max(1e-300);
    let cold_err = x_ref.map(|r| max_abs_diff(&x, r) / x_scale);
    let x_cold = x.clone();

    let mut warm_times = Vec::new();
    let mut update_times = Vec::new();
    let mut warm = SolveStats::default();
    for _ in 0..cli.reps.max(1) {
        amg.update(q);
        let t = Instant::now();
        amg.update(q_warm);
        update_times.push(ms(t));
        x.copy_from_slice(&x_cold);
        let t = Instant::now();
        warm = amg.solve(b_warm, &mut x, cli.tol, cli.max_iterations);
        warm_times.push(ms(t));
    }
    let warm_ms = median(&mut warm_times);
    let update_ms = median(&mut update_times);
    let warm_err = x_ref_warm.map(|r| max_abs_diff(&x, r) / x_scale);
    let rel = cold
        .relative
        .iter()
        .chain(warm.relative.iter())
        .cloned()
        .fold(0.0, f64::max);

    let err_str = |e: Option<f64>| e.map_or("n/a".to_string(), |e| format!("{e:.1e}"));
    println!(
        "  {:<16} levels {} [{}] coarse edges {:.2}×ne, {:.0} MB | setup {setup_ms:.0} ms, A·x {apply_ms:.2} ms | cold it {} {cold_ms:.1} ms ({:.1} ops/it) | warm {}% it {} {warm_ms:.1} ms + update {update_ms:.1} ms | rel_res {rel:.1e} | err vs direct {} / {} | direct {direct_ms:.1} ms",
        params.label(),
        sizes.len(),
        fmt_sizes(&sizes),
        coarse_edges as f64 / ne as f64,
        amg.coarse_bytes() as f64 / 1e6,
        fmt_iters(cold.iterations),
        cold.ops / cold.iterations.iter().copied().max().unwrap_or(1).max(1) as f64,
        cli.warm,
        fmt_iters(warm.iterations),
        err_str(cold_err),
        err_str(warm_err),
    );
    println!(
        "RESULT\t{}\t{}\t{}\t{}\t{}\t{:.1}\t{}\t{}\t{}\t{:.1}\t{}\t{}\t{}\t{:.1}\t{:.1}\t{:.1e}\t{:.1}\t{:.1}\t{:.2}",
        fx.name,
        ne,
        nfree,
        params.label(),
        sizes.len(),
        setup_ms,
        cold.iterations[0],
        cold.iterations[1],
        cold.iterations[2],
        cold_ms,
        warm.iterations[0],
        warm.iterations[1],
        warm.iterations[2],
        warm_ms,
        update_ms,
        rel,
        direct_ms,
        cold.ops / cold.iterations.iter().copied().max().unwrap_or(1).max(1) as f64,
        coarse_edges as f64 / ne as f64,
    );
}

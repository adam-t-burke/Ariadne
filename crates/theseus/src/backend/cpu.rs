//! CPU implementation of [`Backend`]: rayon node-gather kernels over host
//! vectors, and the chunked, fixed-order reductions the rest of the crate
//! uses to stay bitwise reproducible across thread counts.
//!
//! # Determinism
//!
//! Every parallel loop here works on **fixed-size chunks** ([`CHUNK`] rows)
//! whose boundaries depend only on the data length, never on the number of
//! threads. Inside a chunk the work is a plain sequential loop; a
//! reduction produces one partial per chunk and the partials are combined
//! **sequentially in chunk order**. Below [`PAR_MIN_LEN`] elements the same
//! chunk loop runs on the calling thread, so small problems pay no rayon
//! overhead and still compute exactly what the parallel path would.
//!
//! Scatter-shaped operations (`coarse[agg[i]] += fine[i]`, reactions,
//! gradient accumulation) are expressed as gathers over a CSR map whose
//! entries are in ascending source order, which reproduces the sequential
//! scatter's per-target summation order exactly.
//!
//! # Precision
//!
//! [`CpuBuf`] and [`CpuLevelGraph`] exist in `f64` and `f32` variants
//! selected by [`Precision`]. Arithmetic runs in the buffer precision;
//! `dot3` / `norm3` accumulate in `f64` regardless. A kernel is given
//! buffers of one precision only — mixing them is a programming error and
//! panics with a message naming the kernel.

use crate::amg::hierarchy::LevelMatrix;
use crate::backend::{Backend, BackendHandle};
use crate::graph::{CsrAdjacency, LevelGraph};
use crate::linear_solver::Precision;
use rayon::prelude::*;
use std::ops::{Add, AddAssign, Div, Mul, Range, Sub};
use std::sync::Mutex;

/// Rows per chunk of every parallel loop and every fixed-order reduction.
pub const CHUNK: usize = 4096;

/// Below this many elements a loop runs sequentially (same chunking).
pub const PAR_MIN_LEN: usize = 16_384;

// ─────────────────────────────────────────────────────────────
//  Chunked loops and fixed-order reductions
// ─────────────────────────────────────────────────────────────

/// Whether a loop over `len` elements should run on the calling thread:
/// below [`PAR_MIN_LEN`] elements, or when the current rayon pool has a
/// single thread (there is nothing to gain and the scheduling round trip
/// costs more than the work). The chunked computation is the same either
/// way, so this never changes a result.
#[inline]
pub fn run_sequential(len: usize) -> bool {
    len < PAR_MIN_LEN || rayon::current_num_threads() == 1
}

/// Run `f(start, chunk)` over consecutive `chunk`-element pieces of `data`
/// (the last one may be shorter); `start` is the index of `chunk[0]` in
/// `data`. Parallel when `data.len() >= PAR_MIN_LEN` and the pool has more
/// than one thread.
pub fn for_each_chunk_mut<T: Send>(
    data: &mut [T],
    chunk: usize,
    f: impl Fn(usize, &mut [T]) + Sync,
) {
    if run_sequential(data.len()) {
        for (c, piece) in data.chunks_mut(chunk).enumerate() {
            f(c * chunk, piece);
        }
    } else {
        data.par_chunks_mut(chunk)
            .enumerate()
            .for_each(|(c, piece)| f(c * chunk, piece));
    }
}

/// Two-output variant of [`for_each_chunk_mut`]: `a` and `b` are split into
/// aligned chunks of `chunk_a` and `chunk_b` elements (so `a.len() /
/// chunk_a == b.len() / chunk_b` up to rounding). `start` is the chunk
/// index times `chunk_a`.
///
/// # Panics
/// If the two slices do not yield the same number of chunks.
pub fn for_each_chunk_mut2<A: Send, B: Send>(
    a: &mut [A],
    chunk_a: usize,
    b: &mut [B],
    chunk_b: usize,
    f: impl Fn(usize, &mut [A], &mut [B]) + Sync,
) {
    let n_chunks = a.len().div_ceil(chunk_a);
    assert_eq!(
        n_chunks,
        b.len().div_ceil(chunk_b),
        "for_each_chunk_mut2: slices yield different chunk counts"
    );
    if run_sequential(a.len()) {
        for (c, (pa, pb)) in a.chunks_mut(chunk_a).zip(b.chunks_mut(chunk_b)).enumerate() {
            f(c * chunk_a, pa, pb);
        }
    } else {
        a.par_chunks_mut(chunk_a)
            .zip(b.par_chunks_mut(chunk_b))
            .enumerate()
            .for_each(|(c, (pa, pb))| f(c * chunk_a, pa, pb));
    }
}

/// Fixed-order chunked sum: `f(range)` returns the **sequential** sum over
/// one chunk of `0..len` ([`CHUNK`] indices each); the partials are added in
/// chunk order. The result does not depend on the thread count, and for
/// `len <= CHUNK` it is exactly the plain sequential sum.
pub fn deterministic_sum(len: usize, f: impl Fn(Range<usize>) -> f64 + Sync) -> f64 {
    let chunk_range = |c: usize| c * CHUNK..((c + 1) * CHUNK).min(len);
    let n_chunks = len.div_ceil(CHUNK);
    if run_sequential(len) {
        let mut total = 0.0;
        for c in 0..n_chunks {
            total += f(chunk_range(c));
        }
        total
    } else {
        let partials: Vec<f64> = (0..n_chunks)
            .into_par_iter()
            .map(|c| f(chunk_range(c)))
            .collect();
        partials.iter().fold(0.0, |acc, p| acc + p)
    }
}

/// Element-wise sum over `values` in fixed chunk order (see
/// [`deterministic_sum`]).
pub fn deterministic_slice_sum(values: &[f64]) -> f64 {
    deterministic_sum(values.len(), |r| values[r].iter().fold(0.0, |a, v| a + v))
}

// ─────────────────────────────────────────────────────────────
//  Scalars and buffers
// ─────────────────────────────────────────────────────────────

/// Element type of a CPU buffer (`f64` or `f32`).
pub trait CpuScalar:
    Copy
    + Send
    + Sync
    + Default
    + PartialEq
    + std::fmt::Debug
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + AddAssign
    + 'static
{
    /// The [`Precision`] tag of this scalar.
    const PRECISION: Precision;
    /// Convert from `f64` (rounding for `f32`).
    fn from_f64(v: f64) -> Self;
    /// Widen to `f64` (exact).
    fn to_f64(self) -> f64;
}

impl CpuScalar for f64 {
    const PRECISION: Precision = Precision::F64;
    #[inline]
    fn from_f64(v: f64) -> Self {
        v
    }
    #[inline]
    fn to_f64(self) -> f64 {
        self
    }
}

impl CpuScalar for f32 {
    const PRECISION: Precision = Precision::F32;
    #[inline]
    fn from_f64(v: f64) -> Self {
        v as f32
    }
    #[inline]
    fn to_f64(self) -> f64 {
        self as f64
    }
}

/// Host vector of `n * 3` scalars in one precision.
#[derive(Debug, Clone, PartialEq)]
pub enum CpuBuf {
    F64(Vec<f64>),
    F32(Vec<f32>),
}

impl CpuBuf {
    /// Precision of the stored scalars.
    pub fn precision(&self) -> Precision {
        match self {
            Self::F64(_) => Precision::F64,
            Self::F32(_) => Precision::F32,
        }
    }

    /// Number of scalars.
    pub fn len(&self) -> usize {
        match self {
            Self::F64(v) => v.len(),
            Self::F32(v) => v.len(),
        }
    }

    /// `true` when there are no scalars.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The `f64` storage, if this is an `F64` buffer.
    pub fn as_f64(&self) -> Option<&[f64]> {
        match self {
            Self::F64(v) => Some(v),
            Self::F32(_) => None,
        }
    }

    /// The `f32` storage, if this is an `F32` buffer.
    pub fn as_f32(&self) -> Option<&[f32]> {
        match self {
            Self::F32(v) => Some(v),
            Self::F64(_) => None,
        }
    }

    /// Contents widened to `f64`.
    pub fn to_vec_f64(&self) -> Vec<f64> {
        match self {
            Self::F64(v) => v.clone(),
            Self::F32(v) => v.iter().map(|&x| x as f64).collect(),
        }
    }
}

/// One uploaded level in precision `T`: adjacency, per-edge weights,
/// per-node anchors and the reciprocal Jacobi diagonal.
#[derive(Debug, Clone, PartialEq)]
pub struct CpuLevelData<T> {
    pub adjacency: CsrAdjacency,
    pub weight: Vec<T>,
    pub anchor: Vec<T>,
    /// `1 / (anchor_u + Σ_{e ∋ u} weight_e)`, computed in `f64` then rounded.
    pub inv_diag: Vec<T>,
}

impl<T: CpuScalar> CpuLevelData<T> {
    fn new(level: &LevelGraph) -> Self {
        let mut data = Self {
            adjacency: level.adjacency.clone(),
            weight: vec![T::default(); level.weight.len()],
            anchor: vec![T::default(); level.anchor.len()],
            inv_diag: vec![T::default(); level.n],
        };
        data.set_weights(&level.weight, &level.anchor);
        data
    }

    fn set_weights(&mut self, weight: &[f64], anchor: &[f64]) {
        assert_eq!(
            weight.len(),
            self.weight.len(),
            "level weight length changed"
        );
        assert_eq!(
            anchor.len(),
            self.anchor.len(),
            "level anchor length changed"
        );
        for (dst, &w) in self.weight.iter_mut().zip(weight) {
            *dst = T::from_f64(w);
        }
        for (dst, &a) in self.anchor.iter_mut().zip(anchor) {
            *dst = T::from_f64(a);
        }
        let adjacency = &self.adjacency;
        for_each_chunk_mut(&mut self.inv_diag, CHUNK, |start, chunk| {
            for (j, out) in chunk.iter_mut().enumerate() {
                let u = start + j;
                // Same order as `LevelGraph::diagonal`: anchor + Σ weights.
                let mut sum = 0.0;
                for (e, _, _) in adjacency.incident(u) {
                    sum += weight[e as usize];
                }
                *out = T::from_f64(1.0 / (anchor[u] + sum));
            }
        });
    }

    /// Number of nodes.
    pub fn n(&self) -> usize {
        self.inv_diag.len()
    }
}

/// Uploaded [`LevelGraph`] in the precision chosen at upload.
#[derive(Debug, Clone, PartialEq)]
pub enum CpuLevelGraph {
    F64(CpuLevelData<f64>),
    F32(CpuLevelData<f32>),
}

impl CpuLevelGraph {
    /// Precision of the stored weights.
    pub fn precision(&self) -> Precision {
        match self {
            Self::F64(_) => Precision::F64,
            Self::F32(_) => Precision::F32,
        }
    }

    /// Number of nodes.
    pub fn n(&self) -> usize {
        match self {
            Self::F64(l) => l.n(),
            Self::F32(l) => l.n(),
        }
    }
}

/// One uploaded CSR matrix in precision `T` (see [`LevelMatrix`]): pattern
/// shared with the host copy, values converted, reciprocal diagonal for the
/// Jacobi-scaled smoother (empty for rectangular transfers).
#[derive(Debug, Clone, PartialEq)]
pub struct CpuCsrData<T> {
    pub nrows: usize,
    pub ncols: usize,
    pub row_ptr: Vec<u32>,
    pub col_idx: Vec<u32>,
    pub values: Vec<T>,
    /// `1 / diag[u]` (0 where `diag[u] == 0`), computed in `f64` then rounded.
    pub inv_diag: Vec<T>,
}

impl<T: CpuScalar> CpuCsrData<T> {
    fn new(m: &LevelMatrix) -> Self {
        m.check();
        let mut data = Self {
            nrows: m.n,
            ncols: m.ncols,
            row_ptr: m.row_ptr.clone(),
            col_idx: m.col_idx.clone(),
            values: vec![T::default(); m.values.len()],
            inv_diag: vec![T::default(); m.diag.len()],
        };
        data.set_values(m);
        data
    }

    fn set_values(&mut self, m: &LevelMatrix) {
        assert_eq!(m.values.len(), self.values.len(), "csr value count changed");
        assert_eq!(
            m.diag.len(),
            self.inv_diag.len(),
            "csr diagonal length changed"
        );
        assert_eq!(m.n, self.nrows, "csr row count changed");
        for_each_chunk_mut(&mut self.values, CHUNK * 4, |start, chunk| {
            for (dst, &v) in chunk.iter_mut().zip(&m.values[start..]) {
                *dst = T::from_f64(v);
            }
        });
        for_each_chunk_mut(&mut self.inv_diag, CHUNK, |start, chunk| {
            for (dst, &d) in chunk.iter_mut().zip(&m.diag[start..]) {
                *dst = T::from_f64(if d != 0.0 { 1.0 / d } else { 0.0 });
            }
        });
    }
}

/// Uploaded [`LevelMatrix`] in the precision chosen at upload.
#[derive(Debug, Clone, PartialEq)]
pub enum CpuCsr {
    F64(CpuCsrData<f64>),
    F32(CpuCsrData<f32>),
}

impl CpuCsr {
    /// Precision of the stored values.
    pub fn precision(&self) -> Precision {
        match self {
            Self::F64(_) => Precision::F64,
            Self::F32(_) => Precision::F32,
        }
    }

    /// Rows of the matrix.
    pub fn nrows(&self) -> usize {
        match self {
            Self::F64(m) => m.nrows,
            Self::F32(m) => m.nrows,
        }
    }

    /// Columns of the matrix.
    pub fn ncols(&self) -> usize {
        match self {
            Self::F64(m) => m.ncols,
            Self::F32(m) => m.ncols,
        }
    }
}

/// Fine → coarse map of one level, with the inverse CSR (coarse → its fine
/// nodes, ascending) so restriction is a gather.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CpuAggregates {
    pub aggregate_of: Vec<u32>,
    pub n_coarse: usize,
    /// `n_coarse + 1` offsets into `coarse_fine`.
    pub coarse_offsets: Vec<u32>,
    /// Fine nodes of each coarse node, ascending within a coarse node.
    pub coarse_fine: Vec<u32>,
}

impl CpuAggregates {
    fn new(aggregate_of: &[u32], n_coarse: usize) -> Self {
        let mut coarse_offsets = vec![0u32; n_coarse + 1];
        for &c in aggregate_of {
            assert!(
                (c as usize) < n_coarse,
                "aggregate index {c} >= n_coarse {n_coarse}"
            );
            coarse_offsets[c as usize + 1] += 1;
        }
        for c in 0..n_coarse {
            coarse_offsets[c + 1] += coarse_offsets[c];
        }
        let mut fill: Vec<u32> = coarse_offsets[..n_coarse].to_vec();
        let mut coarse_fine = vec![0u32; aggregate_of.len()];
        for (i, &c) in aggregate_of.iter().enumerate() {
            let slot = fill[c as usize] as usize;
            fill[c as usize] += 1;
            coarse_fine[slot] = i as u32;
        }
        Self {
            aggregate_of: aggregate_of.to_vec(),
            n_coarse,
            coarse_offsets,
            coarse_fine,
        }
    }
}

// ─────────────────────────────────────────────────────────────
//  Generic kernels (one precision)
// ─────────────────────────────────────────────────────────────

/// `acc = anchor_u x_u + Σ_{e ∋ u} w_e (x_u − x_other)` for node `u`, in the
/// same operation order as [`LevelGraph::apply`].
#[inline]
fn gather_row<T: CpuScalar>(level: &CpuLevelData<T>, x: &[T], u: usize) -> [T; 3] {
    let xu = [x[u * 3], x[u * 3 + 1], x[u * 3 + 2]];
    let anchor = level.anchor[u];
    let mut acc = [anchor * xu[0], anchor * xu[1], anchor * xu[2]];
    for (e, _, v) in level.adjacency.incident(u) {
        let w = level.weight[e as usize];
        let v = v as usize * 3;
        acc[0] += w * (xu[0] - x[v]);
        acc[1] += w * (xu[1] - x[v + 1]);
        acc[2] += w * (xu[2] - x[v + 2]);
    }
    acc
}

fn apply_graph_t<T: CpuScalar>(level: &CpuLevelData<T>, x: &[T], y: &mut [T]) {
    let n = level.n();
    assert_eq!(x.len(), n * 3, "apply_graph: x length");
    assert_eq!(y.len(), n * 3, "apply_graph: y length");
    for_each_chunk_mut(y, CHUNK * 3, |start, chunk| {
        for (j, row) in chunk.chunks_exact_mut(3).enumerate() {
            let acc = gather_row(level, x, start / 3 + j);
            row.copy_from_slice(&acc);
        }
    });
}

fn residual_t<T: CpuScalar>(level: &CpuLevelData<T>, x: &[T], b: &[T], r: &mut [T]) {
    let n = level.n();
    assert_eq!(x.len(), n * 3, "residual: x length");
    assert_eq!(b.len(), n * 3, "residual: b length");
    assert_eq!(r.len(), n * 3, "residual: r length");
    for_each_chunk_mut(r, CHUNK * 3, |start, chunk| {
        for (j, row) in chunk.chunks_exact_mut(3).enumerate() {
            let i = start + j * 3;
            let acc = gather_row(level, x, i / 3);
            row[0] = b[i] - acc[0];
            row[1] = b[i + 1] - acc[1];
            row[2] = b[i + 2] - acc[2];
        }
    });
}

fn chebyshev_step_t<T: CpuScalar>(
    level: &CpuLevelData<T>,
    alpha: f64,
    beta: f64,
    r: &[T],
    d: &mut [T],
    x: &mut [T],
) {
    chebyshev_step_inv(&level.inv_diag, alpha, beta, r, d, x);
}

/// `Σ_j values[j] · x[col[j]]` over row `u`, ascending column order.
#[inline]
fn csr_row<T: CpuScalar>(m: &CpuCsrData<T>, x: &[T], u: usize) -> [T; 3] {
    let mut acc = [T::default(); 3];
    let range = m.row_ptr[u] as usize..m.row_ptr[u + 1] as usize;
    for (&c, &v) in m.col_idx[range.clone()].iter().zip(&m.values[range]) {
        let c = c as usize * 3;
        acc[0] += v * x[c];
        acc[1] += v * x[c + 1];
        acc[2] += v * x[c + 2];
    }
    acc
}

fn apply_csr_t<T: CpuScalar>(m: &CpuCsrData<T>, x: &[T], y: &mut [T], accumulate: bool) {
    assert_eq!(x.len(), m.ncols * 3, "apply_csr: x length");
    assert_eq!(y.len(), m.nrows * 3, "apply_csr: y length");
    for_each_chunk_mut(y, CHUNK * 3, |start, chunk| {
        for (j, row) in chunk.chunks_exact_mut(3).enumerate() {
            let acc = csr_row(m, x, start / 3 + j);
            if accumulate {
                row[0] += acc[0];
                row[1] += acc[1];
                row[2] += acc[2];
            } else {
                row.copy_from_slice(&acc);
            }
        }
    });
}

fn residual_csr_t<T: CpuScalar>(m: &CpuCsrData<T>, x: &[T], b: &[T], r: &mut [T]) {
    assert_eq!(m.nrows, m.ncols, "residual_csr: matrix must be square");
    let n = m.nrows;
    assert_eq!(x.len(), n * 3, "residual_csr: x length");
    assert_eq!(b.len(), n * 3, "residual_csr: b length");
    assert_eq!(r.len(), n * 3, "residual_csr: r length");
    for_each_chunk_mut(r, CHUNK * 3, |start, chunk| {
        for (j, row) in chunk.chunks_exact_mut(3).enumerate() {
            let i = start + j * 3;
            let acc = csr_row(m, x, i / 3);
            row[0] = b[i] - acc[0];
            row[1] = b[i + 1] - acc[1];
            row[2] = b[i + 2] - acc[2];
        }
    });
}

/// Shared body of `chebyshev_step` / `chebyshev_step_csr` over any
/// reciprocal diagonal.
fn chebyshev_step_inv<T: CpuScalar>(
    inv_diag: &[T],
    alpha: f64,
    beta: f64,
    r: &[T],
    d: &mut [T],
    x: &mut [T],
) {
    let n = inv_diag.len();
    assert_eq!(r.len(), n * 3, "chebyshev_step: r length");
    assert_eq!(d.len(), n * 3, "chebyshev_step: d length");
    assert_eq!(x.len(), n * 3, "chebyshev_step: x length");
    let alpha = T::from_f64(alpha);
    let beta = T::from_f64(beta);
    for_each_chunk_mut2(d, CHUNK * 3, x, CHUNK * 3, |start, dc, xc| {
        for (j, (dv, xv)) in dc.iter_mut().zip(xc.iter_mut()).enumerate() {
            let i = start + j;
            let z = inv_diag[i / 3] * r[i];
            *dv = alpha * z + beta * *dv;
            *xv += *dv;
        }
    });
}

fn restrict_t<T: CpuScalar>(agg: &CpuAggregates, fine: &[T], coarse: &mut [T]) {
    assert_eq!(
        fine.len(),
        agg.aggregate_of.len() * 3,
        "restrict: fine length"
    );
    assert_eq!(coarse.len(), agg.n_coarse * 3, "restrict: coarse length");
    for_each_chunk_mut(coarse, CHUNK * 3, |start, chunk| {
        for (j, row) in chunk.chunks_exact_mut(3).enumerate() {
            let c = start / 3 + j;
            let range = agg.coarse_offsets[c] as usize..agg.coarse_offsets[c + 1] as usize;
            let mut acc = [T::default(); 3];
            for &i in &agg.coarse_fine[range] {
                let i = i as usize * 3;
                acc[0] += fine[i];
                acc[1] += fine[i + 1];
                acc[2] += fine[i + 2];
            }
            row.copy_from_slice(&acc);
        }
    });
}

fn prolong_add_t<T: CpuScalar>(agg: &CpuAggregates, coarse: &[T], fine: &mut [T]) {
    assert_eq!(
        fine.len(),
        agg.aggregate_of.len() * 3,
        "prolong_add: fine length"
    );
    assert_eq!(coarse.len(), agg.n_coarse * 3, "prolong_add: coarse length");
    for_each_chunk_mut(fine, CHUNK * 3, |start, chunk| {
        for (j, row) in chunk.chunks_exact_mut(3).enumerate() {
            let c = agg.aggregate_of[start / 3 + j] as usize * 3;
            row[0] += coarse[c];
            row[1] += coarse[c + 1];
            row[2] += coarse[c + 2];
        }
    });
}

fn axpy_t<T: CpuScalar>(alpha: [f64; 3], x: &[T], y: &mut [T]) {
    assert_eq!(x.len(), y.len(), "axpy: lengths");
    let a = alpha.map(T::from_f64);
    for_each_chunk_mut(y, CHUNK * 3, |start, chunk| {
        for (j, row) in chunk.chunks_exact_mut(3).enumerate() {
            let i = start + j * 3;
            row[0] += a[0] * x[i];
            row[1] += a[1] * x[i + 1];
            row[2] += a[2] * x[i + 2];
        }
    });
}

fn scale_t<T: CpuScalar>(alpha: [f64; 3], x: &mut [T]) {
    let a = alpha.map(T::from_f64);
    for_each_chunk_mut(x, CHUNK * 3, |_, chunk| {
        for row in chunk.chunks_exact_mut(3) {
            row[0] = a[0] * row[0];
            row[1] = a[1] * row[1];
            row[2] = a[2] * row[2];
        }
    });
}

/// Per-column dot product of one chunk of rows, accumulated in `f64`.
#[inline]
fn dot_chunk<T: CpuScalar>(a: &[T], b: &[T]) -> [f64; 3] {
    let mut acc = [0.0f64; 3];
    for (ra, rb) in a.chunks_exact(3).zip(b.chunks_exact(3)) {
        acc[0] += ra[0].to_f64() * rb[0].to_f64();
        acc[1] += ra[1].to_f64() * rb[1].to_f64();
        acc[2] += ra[2].to_f64() * rb[2].to_f64();
    }
    acc
}

fn dot3_t<T: CpuScalar>(partials: &Mutex<Vec<[f64; 3]>>, a: &[T], b: &[T]) -> [f64; 3] {
    assert_eq!(a.len(), b.len(), "dot3: lengths");
    assert_eq!(a.len() % 3, 0, "dot3: length not a multiple of 3");
    let chunk = CHUNK * 3;
    let n_chunks = a.len().div_ceil(chunk);
    let mut total = [0.0f64; 3];
    if run_sequential(a.len()) {
        for (ca, cb) in a.chunks(chunk).zip(b.chunks(chunk)) {
            let p = dot_chunk(ca, cb);
            total[0] += p[0];
            total[1] += p[1];
            total[2] += p[2];
        }
    } else {
        let mut partials = partials.lock().unwrap_or_else(|e| e.into_inner());
        partials.clear();
        partials.resize(n_chunks, [0.0; 3]);
        partials
            .par_iter_mut()
            .zip(a.par_chunks(chunk).zip(b.par_chunks(chunk)))
            .for_each(|(out, (ca, cb))| *out = dot_chunk(ca, cb));
        for p in partials.iter() {
            total[0] += p[0];
            total[1] += p[1];
            total[2] += p[2];
        }
    }
    total
}

// ─────────────────────────────────────────────────────────────
//  The backend
// ─────────────────────────────────────────────────────────────

/// Host backend: rayon kernels over [`CpuBuf`] with fixed-order reductions.
///
/// `Send + Sync`; the only state is a grow-only scratch of per-chunk
/// partials so `dot3` / `norm3` do not allocate once warmed up.
#[derive(Debug, Default)]
pub struct CpuBackend {
    partials: Mutex<Vec<[f64; 3]>>,
}

impl CpuBackend {
    /// A backend with an empty partials scratch.
    pub fn new() -> Self {
        Self::default()
    }

    /// Reserve the reduction scratch for vectors of up to `len` scalars.
    pub fn reserve(&self, len: usize) {
        let n_chunks = len.div_ceil(CHUNK * 3);
        let mut partials = self.partials.lock().unwrap_or_else(|e| e.into_inner());
        let len = partials.len();
        if partials.capacity() < n_chunks {
            partials.reserve(n_chunks - len);
        }
    }
}

#[cold]
#[track_caller]
fn precision_mismatch(kernel: &str) -> ! {
    panic!("CpuBackend::{kernel}: buffers and level must share one precision")
}

impl Backend for CpuBackend {
    type Buf = CpuBuf;
    type LevelGraphBuf = CpuLevelGraph;
    type AggBuf = CpuAggregates;
    type CsrBuf = CpuCsr;

    fn handle(&self) -> BackendHandle {
        BackendHandle::Cpu
    }

    fn alloc(&self, len: usize, precision: Precision) -> CpuBuf {
        match precision {
            Precision::F64 => CpuBuf::F64(vec![0.0; len]),
            Precision::F32 => CpuBuf::F32(vec![0.0; len]),
        }
    }

    fn upload(&self, src: &[f64], dst: &mut CpuBuf) {
        assert_eq!(src.len(), dst.len(), "upload: lengths");
        match dst {
            CpuBuf::F64(v) => v.copy_from_slice(src),
            CpuBuf::F32(v) => {
                for (d, &s) in v.iter_mut().zip(src) {
                    *d = s as f32;
                }
            }
        }
    }

    fn download(&self, src: &CpuBuf, dst: &mut [f64]) {
        assert_eq!(src.len(), dst.len(), "download: lengths");
        match src {
            CpuBuf::F64(v) => dst.copy_from_slice(v),
            CpuBuf::F32(v) => {
                for (d, &s) in dst.iter_mut().zip(v) {
                    *d = s as f64;
                }
            }
        }
    }

    fn len(&self, buf: &CpuBuf) -> usize {
        buf.len()
    }

    fn copy(&self, src: &CpuBuf, dst: &mut CpuBuf) {
        match (src, dst) {
            (CpuBuf::F64(s), CpuBuf::F64(d)) => d.copy_from_slice(s),
            (CpuBuf::F32(s), CpuBuf::F32(d)) => d.copy_from_slice(s),
            _ => precision_mismatch("copy"),
        }
    }

    fn zero(&self, buf: &mut CpuBuf) {
        match buf {
            CpuBuf::F64(v) => v.fill(0.0),
            CpuBuf::F32(v) => v.fill(0.0),
        }
    }

    fn upload_level(&self, level: &LevelGraph, precision: Precision) -> CpuLevelGraph {
        match precision {
            Precision::F64 => CpuLevelGraph::F64(CpuLevelData::new(level)),
            Precision::F32 => CpuLevelGraph::F32(CpuLevelData::new(level)),
        }
    }

    fn update_level_weights(&self, level: &mut CpuLevelGraph, weight: &[f64], anchor: &[f64]) {
        match level {
            CpuLevelGraph::F64(l) => l.set_weights(weight, anchor),
            CpuLevelGraph::F32(l) => l.set_weights(weight, anchor),
        }
    }

    fn upload_aggregates(&self, aggregate_of: &[u32], n_coarse: usize) -> CpuAggregates {
        CpuAggregates::new(aggregate_of, n_coarse)
    }

    fn upload_csr(&self, m: &LevelMatrix, precision: Precision) -> CpuCsr {
        match precision {
            Precision::F64 => CpuCsr::F64(CpuCsrData::new(m)),
            Precision::F32 => CpuCsr::F32(CpuCsrData::new(m)),
        }
    }

    fn update_csr_values(&self, m: &LevelMatrix, dst: &mut CpuCsr) {
        match dst {
            CpuCsr::F64(d) => d.set_values(m),
            CpuCsr::F32(d) => d.set_values(m),
        }
    }

    fn apply_csr(&self, m: &CpuCsr, x: &CpuBuf, y: &mut CpuBuf) {
        match (m, x, y) {
            (CpuCsr::F64(m), CpuBuf::F64(x), CpuBuf::F64(y)) => apply_csr_t(m, x, y, false),
            (CpuCsr::F32(m), CpuBuf::F32(x), CpuBuf::F32(y)) => apply_csr_t(m, x, y, false),
            _ => precision_mismatch("apply_csr"),
        }
    }

    fn apply_csr_add(&self, m: &CpuCsr, x: &CpuBuf, y: &mut CpuBuf) {
        match (m, x, y) {
            (CpuCsr::F64(m), CpuBuf::F64(x), CpuBuf::F64(y)) => apply_csr_t(m, x, y, true),
            (CpuCsr::F32(m), CpuBuf::F32(x), CpuBuf::F32(y)) => apply_csr_t(m, x, y, true),
            _ => precision_mismatch("apply_csr_add"),
        }
    }

    fn residual_csr(&self, m: &CpuCsr, x: &CpuBuf, b: &CpuBuf, r: &mut CpuBuf) {
        match (m, x, b, r) {
            (CpuCsr::F64(m), CpuBuf::F64(x), CpuBuf::F64(b), CpuBuf::F64(r)) => {
                residual_csr_t(m, x, b, r)
            }
            (CpuCsr::F32(m), CpuBuf::F32(x), CpuBuf::F32(b), CpuBuf::F32(r)) => {
                residual_csr_t(m, x, b, r)
            }
            _ => precision_mismatch("residual_csr"),
        }
    }

    fn chebyshev_step_csr(
        &self,
        m: &CpuCsr,
        alpha: f64,
        beta: f64,
        r: &CpuBuf,
        d: &mut CpuBuf,
        x: &mut CpuBuf,
    ) {
        match (m, r, d, x) {
            (CpuCsr::F64(m), CpuBuf::F64(r), CpuBuf::F64(d), CpuBuf::F64(x)) => {
                chebyshev_step_inv(&m.inv_diag, alpha, beta, r, d, x)
            }
            (CpuCsr::F32(m), CpuBuf::F32(r), CpuBuf::F32(d), CpuBuf::F32(x)) => {
                chebyshev_step_inv(&m.inv_diag, alpha, beta, r, d, x)
            }
            _ => precision_mismatch("chebyshev_step_csr"),
        }
    }

    fn apply_graph(&self, level: &CpuLevelGraph, x: &CpuBuf, y: &mut CpuBuf) {
        match (level, x, y) {
            (CpuLevelGraph::F64(l), CpuBuf::F64(x), CpuBuf::F64(y)) => apply_graph_t(l, x, y),
            (CpuLevelGraph::F32(l), CpuBuf::F32(x), CpuBuf::F32(y)) => apply_graph_t(l, x, y),
            _ => precision_mismatch("apply_graph"),
        }
    }

    fn residual(&self, level: &CpuLevelGraph, x: &CpuBuf, b: &CpuBuf, r: &mut CpuBuf) {
        match (level, x, b, r) {
            (CpuLevelGraph::F64(l), CpuBuf::F64(x), CpuBuf::F64(b), CpuBuf::F64(r)) => {
                residual_t(l, x, b, r)
            }
            (CpuLevelGraph::F32(l), CpuBuf::F32(x), CpuBuf::F32(b), CpuBuf::F32(r)) => {
                residual_t(l, x, b, r)
            }
            _ => precision_mismatch("residual"),
        }
    }

    fn chebyshev_step(
        &self,
        level: &CpuLevelGraph,
        alpha: f64,
        beta: f64,
        r: &CpuBuf,
        d: &mut CpuBuf,
        x: &mut CpuBuf,
    ) {
        match (level, r, d, x) {
            (CpuLevelGraph::F64(l), CpuBuf::F64(r), CpuBuf::F64(d), CpuBuf::F64(x)) => {
                chebyshev_step_t(l, alpha, beta, r, d, x)
            }
            (CpuLevelGraph::F32(l), CpuBuf::F32(r), CpuBuf::F32(d), CpuBuf::F32(x)) => {
                chebyshev_step_t(l, alpha, beta, r, d, x)
            }
            _ => precision_mismatch("chebyshev_step"),
        }
    }

    fn restrict(&self, agg: &CpuAggregates, fine: &CpuBuf, coarse: &mut CpuBuf) {
        match (fine, coarse) {
            (CpuBuf::F64(f), CpuBuf::F64(c)) => restrict_t(agg, f, c),
            (CpuBuf::F32(f), CpuBuf::F32(c)) => restrict_t(agg, f, c),
            _ => precision_mismatch("restrict"),
        }
    }

    fn prolong_add(&self, agg: &CpuAggregates, coarse: &CpuBuf, fine: &mut CpuBuf) {
        match (coarse, fine) {
            (CpuBuf::F64(c), CpuBuf::F64(f)) => prolong_add_t(agg, c, f),
            (CpuBuf::F32(c), CpuBuf::F32(f)) => prolong_add_t(agg, c, f),
            _ => precision_mismatch("prolong_add"),
        }
    }

    fn axpy(&self, alpha: [f64; 3], x: &CpuBuf, y: &mut CpuBuf) {
        match (x, y) {
            (CpuBuf::F64(x), CpuBuf::F64(y)) => axpy_t(alpha, x, y),
            (CpuBuf::F32(x), CpuBuf::F32(y)) => axpy_t(alpha, x, y),
            _ => precision_mismatch("axpy"),
        }
    }

    fn scale(&self, alpha: [f64; 3], x: &mut CpuBuf) {
        match x {
            CpuBuf::F64(x) => scale_t(alpha, x),
            CpuBuf::F32(x) => scale_t(alpha, x),
        }
    }

    fn dot3(&self, a: &CpuBuf, b: &CpuBuf) -> [f64; 3] {
        match (a, b) {
            (CpuBuf::F64(a), CpuBuf::F64(b)) => dot3_t(&self.partials, a, b),
            (CpuBuf::F32(a), CpuBuf::F32(b)) => dot3_t(&self.partials, a, b),
            _ => precision_mismatch("dot3"),
        }
    }

    fn norm3(&self, a: &CpuBuf) -> [f64; 3] {
        self.dot3(a, a).map(f64::sqrt)
    }

    fn sync(&self) {}
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn deterministic_sum_is_plain_sum_for_one_chunk() {
        let values: Vec<f64> = (0..1000).map(|i| (i as f64).sin()).collect();
        let plain = values.iter().fold(0.0, |a, v| a + v);
        assert_eq!(deterministic_slice_sum(&values), plain);
    }

    #[test]
    fn deterministic_sum_matches_chunk_partials() {
        let n = 3 * CHUNK + 17;
        let values: Vec<f64> = (0..n).map(|i| 1.0 / (i as f64 + 1.0)).collect();
        let mut expected = 0.0;
        for c in values.chunks(CHUNK) {
            expected += c.iter().fold(0.0, |a, v| a + v);
        }
        assert_eq!(deterministic_slice_sum(&values), expected);
        // Parallel path (above PAR_MIN_LEN) has the same chunk structure.
        let n = 5 * CHUNK + 3;
        let values: Vec<f64> = (0..n).map(|i| 1.0 / (i as f64 + 1.0)).collect();
        let mut expected = 0.0;
        for c in values.chunks(CHUNK) {
            expected += c.iter().fold(0.0, |a, v| a + v);
        }
        assert!(n >= PAR_MIN_LEN);
        assert_eq!(deterministic_slice_sum(&values), expected);
    }

    #[test]
    fn for_each_chunk_reports_starts() {
        let mut data = vec![0usize; 10_000];
        for_each_chunk_mut(&mut data, 3000, |start, chunk| {
            for (j, v) in chunk.iter_mut().enumerate() {
                *v = start + j;
            }
        });
        assert!(data.iter().enumerate().all(|(i, &v)| i == v));
        let mut data = vec![0usize; 20_000];
        for_each_chunk_mut(&mut data, 4096, |start, chunk| {
            for (j, v) in chunk.iter_mut().enumerate() {
                *v = start + j;
            }
        });
        assert!(data.iter().enumerate().all(|(i, &v)| i == v));
    }

    #[test]
    fn aggregates_inverse_map_is_ascending() {
        let agg = CpuAggregates::new(&[1, 0, 1, 2, 0], 3);
        assert_eq!(agg.coarse_offsets, vec![0, 2, 4, 5]);
        assert_eq!(agg.coarse_fine, vec![1, 4, 0, 2, 3]);
    }

    #[test]
    #[should_panic(expected = "share one precision")]
    fn mixed_precision_panics() {
        let be = CpuBackend::new();
        let a = be.alloc(6, Precision::F64);
        let mut b = be.alloc(6, Precision::F32);
        be.copy(&a, &mut b);
    }
}

//! Multi-right-hand-side triangular solves on `faer-sparse` Cholesky factors.
//!
//! `faer`'s own `solve_in_place_with_conj` is generic over the entity type,
//! works on column-major right-hand sides and calls dense BLAS-3 style kernels
//! per supernode. For the equilibrium systems solved here — three right-hand
//! sides (x, y, z) and factors whose supernodes are mostly a handful of columns
//! wide — that overhead is comparable to the arithmetic itself. The routines in
//! this module operate directly on `faer`'s public factor layout with the `K`
//! right-hand sides of a row stored contiguously (`x[row * K + j]`), so a factor
//! column is streamed once and the inner loops vectorise over `K`.
//!
//! Layout facts relied upon (faer-sparse 0.17):
//! * simplicial: `L` (or unit `L` with `D` on the diagonal for LDLᵀ) is stored
//!   in CSC with the diagonal entry first in every column;
//! * supernodal: supernode `s` owns columns `begin[s]..begin[s+1]`, its values
//!   are a column-major `(pattern.len() + size) × size` block whose top `size`
//!   rows are the (lower) diagonal block, and `pattern` lists the rows of the
//!   off-diagonal block; for LDLᵀ the diagonal block is unit lower and holds
//!   `D⁻¹` on the diagonal.
//!
//! The solve on the permuted system is exposed for testing; the `solve` entry
//! point applies the fill-reducing permutation on the way in and out.

use faer_sparse::cholesky::{SymbolicCholesky, SymbolicCholeskyRaw};

/// Which factorization the numeric values represent.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(crate) enum Kind {
    Llt,
    Ldlt,
}

/// Solve `A X = B` for `K` right-hand sides using a `faer` Cholesky factor.
///
/// `rhs` and `out` are row-major `n × K` slices. `work` must hold at least
/// `2 * n * K` values and is used for the permuted solution and supernodal
/// gather/scatter staging.
pub(crate) fn solve<const K: usize>(
    symbolic: &SymbolicCholesky<u32>,
    values: &[f64],
    kind: Kind,
    rhs: &[f64],
    out: &mut [f64],
    work: &mut [f64],
) {
    let n = symbolic.nrows();
    debug_assert_eq!(rhs.len(), n * K);
    debug_assert_eq!(out.len(), n * K);
    assert!(work.len() >= 2 * n * K, "solve workspace too small");

    let (x, tmp) = work.split_at_mut(n * K);
    let tmp = &mut tmp[..n * K];
    let (fwd, inv) = symbolic.perm().into_arrays();

    for (i, &src) in fwd.iter().enumerate() {
        let src = src as usize;
        x[i * K..(i + 1) * K].copy_from_slice(&rhs[src * K..(src + 1) * K]);
    }

    solve_permuted_in_place::<K>(symbolic, values, kind, x, tmp);

    for (i, &src) in inv.iter().enumerate() {
        let src = src as usize;
        out[i * K..(i + 1) * K].copy_from_slice(&x[src * K..(src + 1) * K]);
    }
}

/// Solve `(P A Pᵀ) Y = X` in place on the already-permuted right-hand side.
pub(crate) fn solve_permuted_in_place<const K: usize>(
    symbolic: &SymbolicCholesky<u32>,
    values: &[f64],
    kind: Kind,
    x: &mut [f64],
    tmp: &mut [f64],
) {
    match symbolic.raw() {
        SymbolicCholeskyRaw::Simplicial(sym) => {
            let col_ptrs = sym.col_ptrs();
            let row_indices = sym.row_indices();
            simplicial::<K>(col_ptrs, row_indices, values, kind, x);
        }
        SymbolicCholeskyRaw::Supernodal(sym) => {
            let layout = Supernodal {
                begin: sym.supernode_begin(),
                end: sym.supernode_end(),
                col_ptrs_for_row_indices: sym.col_ptrs_for_row_indices(),
                col_ptrs_for_values: sym.col_ptrs_for_values(),
                row_indices: sym.row_indices(),
                n_supernodes: sym.n_supernodes(),
            };
            supernodal::<K>(&layout, values, kind, x, tmp);
        }
    }
}

// ─────────────────────────────────────────────────────────────
//  Simplicial factors (CSC, diagonal first in each column)
// ─────────────────────────────────────────────────────────────

fn simplicial<const K: usize>(
    col_ptrs: &[u32],
    row_indices: &[u32],
    values: &[f64],
    kind: Kind,
    x: &mut [f64],
) {
    let n = col_ptrs.len() - 1;

    // Forward: L y = b (LLᵀ) or L y = b with unit L (LDLᵀ).
    for j in 0..n {
        let start = col_ptrs[j] as usize;
        let end = col_ptrs[j + 1] as usize;
        let mut xj = [0.0; K];
        xj.copy_from_slice(&x[j * K..(j + 1) * K]);
        if kind == Kind::Llt {
            let inv_d = 1.0 / values[start];
            for v in xj.iter_mut() {
                *v *= inv_d;
            }
            x[j * K..(j + 1) * K].copy_from_slice(&xj);
        }
        for p in start + 1..end {
            let row = row_indices[p] as usize;
            let l = values[p];
            let xr = &mut x[row * K..(row + 1) * K];
            for c in 0..K {
                xr[c] -= l * xj[c];
            }
        }
    }

    // Diagonal scaling for LDLᵀ (D stored as the first entry of each column).
    if kind == Kind::Ldlt {
        for j in 0..n {
            let inv_d = 1.0 / values[col_ptrs[j] as usize];
            for v in &mut x[j * K..(j + 1) * K] {
                *v *= inv_d;
            }
        }
    }

    // Backward: Lᵀ x = y.
    for j in (0..n).rev() {
        let start = col_ptrs[j] as usize;
        let end = col_ptrs[j + 1] as usize;
        let mut acc = [0.0; K];
        for p in start + 1..end {
            let row = row_indices[p] as usize;
            let l = values[p];
            let xr = &x[row * K..(row + 1) * K];
            for c in 0..K {
                acc[c] += l * xr[c];
            }
        }
        let xj = &mut x[j * K..(j + 1) * K];
        if kind == Kind::Llt {
            let inv_d = 1.0 / values[start];
            for c in 0..K {
                xj[c] = (xj[c] - acc[c]) * inv_d;
            }
        } else {
            for c in 0..K {
                xj[c] -= acc[c];
            }
        }
    }
}

// ─────────────────────────────────────────────────────────────
//  Supernodal factors
// ─────────────────────────────────────────────────────────────

struct Supernodal<'a> {
    begin: &'a [u32],
    end: &'a [u32],
    col_ptrs_for_row_indices: &'a [u32],
    col_ptrs_for_values: &'a [u32],
    row_indices: &'a [u32],
    n_supernodes: usize,
}

impl Supernodal<'_> {
    /// `(first column, column count, off-diagonal row pattern, values block)`.
    #[inline]
    fn node<'v>(&self, s: usize, values: &'v [f64]) -> (usize, usize, &[u32], &'v [f64]) {
        let start = self.begin[s] as usize;
        let size = self.end[s] as usize - start;
        let pattern = &self.row_indices[self.col_ptrs_for_row_indices[s] as usize
            ..self.col_ptrs_for_row_indices[s + 1] as usize];
        let vals =
            &values[self.col_ptrs_for_values[s] as usize..self.col_ptrs_for_values[s + 1] as usize];
        (start, size, pattern, vals)
    }
}

fn supernodal<const K: usize>(
    layout: &Supernodal<'_>,
    values: &[f64],
    kind: Kind,
    x: &mut [f64],
    tmp: &mut [f64],
) {
    // Forward substitution, supernodes in postorder (children before parents).
    for s in 0..layout.n_supernodes {
        let (start, size, pattern, vals) = layout.node(s, values);
        let nrows = size + pattern.len();
        let top = &mut x[start * K..(start + size) * K];

        // Dense lower-triangular solve on the diagonal block.
        for c in 0..size {
            let col = &vals[c * nrows..(c + 1) * nrows];
            let mut xc = [0.0; K];
            xc.copy_from_slice(&top[c * K..(c + 1) * K]);
            if kind == Kind::Llt {
                let inv_d = 1.0 / col[c];
                for v in xc.iter_mut() {
                    *v *= inv_d;
                }
                top[c * K..(c + 1) * K].copy_from_slice(&xc);
            }
            for r in c + 1..size {
                let l = col[r];
                let xr = &mut top[r * K..(r + 1) * K];
                for k in 0..K {
                    xr[k] -= l * xc[k];
                }
            }
        }

        if pattern.is_empty() {
            continue;
        }

        // tmp = L_bot · x_top, accumulated column by column (contiguous reads).
        let stage = &mut tmp[..pattern.len() * K];
        stage.fill(0.0);
        for c in 0..size {
            let col = &vals[c * nrows + size..(c + 1) * nrows];
            let mut xc = [0.0; K];
            xc.copy_from_slice(&top[c * K..(c + 1) * K]);
            for (idx, &l) in col.iter().enumerate() {
                let t = &mut stage[idx * K..(idx + 1) * K];
                for k in 0..K {
                    t[k] += l * xc[k];
                }
            }
        }
        for (idx, &row) in pattern.iter().enumerate() {
            let row = row as usize;
            let xr = &mut x[row * K..(row + 1) * K];
            let t = &stage[idx * K..(idx + 1) * K];
            for k in 0..K {
                xr[k] -= t[k];
            }
        }
    }

    // LDLᵀ: the diagonal block stores D⁻¹ on its diagonal.
    if kind == Kind::Ldlt {
        for s in 0..layout.n_supernodes {
            let (start, size, pattern, vals) = layout.node(s, values);
            let nrows = size + pattern.len();
            for c in 0..size {
                let inv_d = vals[c * nrows + c];
                for v in &mut x[(start + c) * K..(start + c + 1) * K] {
                    *v *= inv_d;
                }
            }
        }
    }

    // Backward substitution, parents before children.
    for s in (0..layout.n_supernodes).rev() {
        let (start, size, pattern, vals) = layout.node(s, values);
        let nrows = size + pattern.len();

        // Gather the (already final) pattern rows once.
        let stage = &mut tmp[..pattern.len() * K];
        for (idx, &row) in pattern.iter().enumerate() {
            let row = row as usize;
            stage[idx * K..(idx + 1) * K].copy_from_slice(&x[row * K..(row + 1) * K]);
        }

        let top = &mut x[start * K..(start + size) * K];
        for c in (0..size).rev() {
            let col = &vals[c * nrows..(c + 1) * nrows];
            let mut acc = [0.0; K];
            // Off-diagonal block: L_botᵀ · x_pattern.
            for (idx, &l) in col[size..].iter().enumerate() {
                let t = &stage[idx * K..(idx + 1) * K];
                for k in 0..K {
                    acc[k] += l * t[k];
                }
            }
            // Strictly-lower part of the diagonal block, rows already final.
            for r in c + 1..size {
                let l = col[r];
                let xr = &top[r * K..(r + 1) * K];
                for k in 0..K {
                    acc[k] += l * xr[k];
                }
            }
            let xc = &mut top[c * K..(c + 1) * K];
            if kind == Kind::Llt {
                let inv_d = 1.0 / col[c];
                for k in 0..K {
                    xc[k] = (xc[k] - acc[k]) * inv_d;
                }
            } else {
                for k in 0..K {
                    xc[k] -= acc[k];
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use dyn_stack::{GlobalPodBuffer, PodStack};
    use faer_core::{Conj, Mat, Parallelism, Side};
    use faer_sparse::cholesky::{
        factorize_symbolic_cholesky, CholeskySymbolicParams, LdltRef, LdltRegularization, LltRef,
        LltRegularization,
    };
    use faer_sparse::SupernodalThreshold;

    /// Weighted grid Laplacian plus a diagonal shift: SPD with a 2-D fill pattern.
    fn grid_matrix(n: usize) -> crate::sparse::SparseColMatOwned {
        let mut rows = Vec::new();
        let mut cols = Vec::new();
        let mut vals = Vec::new();
        let mut seed = 0x9e37_79b9u64;
        let mut rnd = || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            0.5 + (seed % 1000) as f64 / 1000.0
        };
        let idx = |r: usize, c: usize| r * n + c;
        let dim = n * n;
        let mut diag = vec![0.1; dim];
        let mut push = |a: usize, b: usize, w: f64, diag: &mut Vec<f64>| {
            rows.push(a);
            cols.push(b);
            vals.push(-w);
            rows.push(b);
            cols.push(a);
            vals.push(-w);
            diag[a] += w;
            diag[b] += w;
        };
        for r in 0..n {
            for c in 0..n {
                if c + 1 < n {
                    push(idx(r, c), idx(r, c + 1), rnd(), &mut diag);
                }
                if r + 1 < n {
                    push(idx(r, c), idx(r + 1, c), rnd(), &mut diag);
                }
            }
        }
        for (i, d) in diag.iter().enumerate() {
            rows.push(i);
            cols.push(i);
            vals.push(*d);
        }
        crate::sparse::SparseColMatOwned::from_coo(dim, dim, &rows, &cols, &vals).unwrap()
    }

    fn check<const K: usize>(n: usize, threshold: SupernodalThreshold, kind: Kind) {
        let a = grid_matrix(n);
        let a_ref = a.as_faer_ref();
        let dim = a.nrows;
        let params = CholeskySymbolicParams {
            supernodal_flop_ratio_threshold: threshold,
            ..Default::default()
        };
        let symbolic = factorize_symbolic_cholesky(a_ref.symbolic(), Side::Upper, params).unwrap();
        let mut values = vec![0.0; symbolic.len_values()];
        match kind {
            Kind::Llt => {
                let mut stack = GlobalPodBuffer::new(
                    symbolic
                        .factorize_numeric_llt_req::<f64>(Parallelism::None)
                        .unwrap(),
                );
                symbolic
                    .factorize_numeric_llt(
                        values.as_mut_slice(),
                        a_ref,
                        Side::Upper,
                        LltRegularization::default(),
                        Parallelism::None,
                        PodStack::new(&mut stack),
                    )
                    .unwrap();
            }
            Kind::Ldlt => {
                let mut stack = GlobalPodBuffer::new(
                    symbolic
                        .factorize_numeric_ldlt_req::<f64>(false, Parallelism::None)
                        .unwrap(),
                );
                symbolic.factorize_numeric_ldlt(
                    values.as_mut_slice(),
                    a_ref,
                    Side::Upper,
                    LdltRegularization::default(),
                    Parallelism::None,
                    PodStack::new(&mut stack),
                );
            }
        }

        let rhs: Vec<f64> = (0..dim * K)
            .map(|i| ((i * 7919) % 101) as f64 / 50.0 - 1.0)
            .collect();

        // Reference: faer's own solve on a column-major copy.
        let mut reference = Mat::<f64>::from_fn(dim, K, |i, j| rhs[i * K + j]);
        let mut stack = GlobalPodBuffer::new(symbolic.solve_in_place_req::<f64>(K).unwrap());
        match kind {
            Kind::Llt => LltRef::new(&symbolic, values.as_slice()).solve_in_place_with_conj(
                Conj::No,
                reference.as_mut(),
                Parallelism::None,
                PodStack::new(&mut stack),
            ),
            Kind::Ldlt => LdltRef::new(&symbolic, values.as_slice()).solve_in_place_with_conj(
                Conj::No,
                reference.as_mut(),
                Parallelism::None,
                PodStack::new(&mut stack),
            ),
        }

        let mut out = vec![0.0; dim * K];
        let mut work = vec![0.0; 2 * dim * K];
        solve::<K>(&symbolic, &values, kind, &rhs, &mut out, &mut work);

        let scale = reference.norm_max().max(1.0);
        for i in 0..dim {
            for j in 0..K {
                let diff = (out[i * K + j] - reference.read(i, j)).abs();
                assert!(
                    diff <= 1e-10 * scale,
                    "{kind:?} {threshold:?} n={n} K={K}: mismatch at ({i},{j}): {} vs {}",
                    out[i * K + j],
                    reference.read(i, j)
                );
            }
        }
    }

    #[test]
    fn matches_faer_simplicial() {
        for kind in [Kind::Llt, Kind::Ldlt] {
            check::<3>(9, SupernodalThreshold::FORCE_SIMPLICIAL, kind);
            check::<1>(23, SupernodalThreshold::FORCE_SIMPLICIAL, kind);
            check::<3>(40, SupernodalThreshold::FORCE_SIMPLICIAL, kind);
        }
    }

    #[test]
    fn matches_faer_supernodal() {
        for kind in [Kind::Llt, Kind::Ldlt] {
            check::<3>(9, SupernodalThreshold::FORCE_SUPERNODAL, kind);
            check::<1>(23, SupernodalThreshold::FORCE_SUPERNODAL, kind);
            check::<3>(40, SupernodalThreshold::FORCE_SUPERNODAL, kind);
        }
    }

    #[test]
    fn matches_faer_auto_selection() {
        for kind in [Kind::Llt, Kind::Ldlt] {
            check::<3>(60, SupernodalThreshold::AUTO, kind);
        }
    }
}

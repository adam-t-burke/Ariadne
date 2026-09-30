//! Faer-native sparse CSC matrix and utilities.
//!
//! Replaces sprs for zero-copy integration with faer-sparse factorization.

use faer_core::sparse::{SparseColMatRef, SymbolicSparseColMatRef};
use std::cmp::Ordering;

/// Owned sparse matrix in CSC (Compressed Sparse Column) format.
/// Compatible with faer's SparseColMatRef.
#[derive(Debug, Clone)]
pub struct SparseColMatOwned {
    pub nrows: usize,
    pub ncols: usize,
    /// Column pointers: length ncols + 1, col_ptrs[j]..col_ptrs[j+1] gives indices into row_indices/values for column j.
    pub col_ptrs: Vec<u32>,
    /// Row index for each nonzero (same order as values).
    pub row_indices: Vec<u32>,
    /// Nonzero values.
    pub values: Vec<f64>,
}

impl SparseColMatOwned {
    /// Build CSC from COO (row, col, val) arrays.
    /// Sorts by (col, row) and merges duplicates by summing values.
    pub fn from_coo(
        nrows: usize,
        ncols: usize,
        rows: &[usize],
        cols: &[usize],
        vals: &[f64],
    ) -> Result<Self, String> {
        if rows.len() != cols.len() || rows.len() != vals.len() {
            return Err("COO arrays must have same length".into());
        }
        if nrows > u32::MAX as usize || ncols > u32::MAX as usize {
            return Err("Matrix dimensions exceed u32::MAX".into());
        }

        let triplets: Vec<(u32, u32, f64)> = rows
            .iter()
            .zip(cols.iter())
            .zip(vals.iter())
            .map(|((&r, &c), &v)| {
                if r >= nrows || c >= ncols {
                    panic!("COO index out of bounds");
                }
                (r as u32, c as u32, v)
            })
            .collect();

        Self::from_triplets(nrows, ncols, &triplets)
    }

    /// Build CSC from triplets (row, col, val). Sorts and merges duplicates.
    pub fn from_triplets(
        nrows: usize,
        ncols: usize,
        triplets: &[(u32, u32, f64)],
    ) -> Result<Self, String> {
        if nrows > u32::MAX as usize || ncols > u32::MAX as usize {
            return Err("Matrix dimensions exceed u32::MAX".into());
        }

        let mut triplets = triplets.to_vec();
        triplets.sort_by(|a, b| match a.1.cmp(&b.1) {
            Ordering::Equal => a.0.cmp(&b.0),
            o => o,
        });

        let mut merged: Vec<(u32, u32, f64)> = Vec::with_capacity(triplets.len());
        for (r, c, v) in triplets {
            if c >= ncols as u32 || r >= nrows as u32 {
                return Err("Triplet index out of bounds".into());
            }
            if let Some((lr, lc, lv)) = merged.last_mut() {
                if *lr == r && *lc == c {
                    *lv += v;
                    continue;
                }
            }
            merged.push((r, c, v));
        }

        let mut col_ptrs = vec![0u32; ncols + 1];
        let mut row_indices = Vec::with_capacity(merged.len());
        let mut values = Vec::with_capacity(merged.len());

        let mut current_col = 0u32;
        for (r, c, v) in merged {
            if c > current_col {
                let len = row_indices.len() as u32;
                for j in (current_col + 1)..=c {
                    col_ptrs[j as usize] = len;
                }
                current_col = c;
            }
            row_indices.push(r);
            values.push(v);
        }
        let len = row_indices.len() as u32;
        for j in (current_col as usize + 1)..=ncols {
            col_ptrs[j] = len;
        }

        Ok(Self {
            nrows,
            ncols,
            col_ptrs,
            row_indices,
            values,
        })
    }

    /// Create from triplets without merging duplicates (for building M in inverse).
    /// Triplets must be sorted by (col, row).
    pub fn from_triplets_sorted(
        nrows: usize,
        ncols: usize,
        triplets: &[(u32, u32, f64)],
    ) -> Result<Self, String> {
        if nrows > u32::MAX as usize || ncols > u32::MAX as usize {
            return Err("Matrix dimensions exceed u32::MAX".into());
        }

        let mut col_ptrs = vec![0u32; ncols + 1];
        let mut row_indices = Vec::with_capacity(triplets.len());
        let mut values = Vec::with_capacity(triplets.len());

        let mut current_col = 0u32;
        for &(r, c, v) in triplets {
            if c >= ncols as u32 || r >= nrows as u32 {
                return Err("Triplet index out of bounds".into());
            }
            if c > current_col {
                let len = row_indices.len() as u32;
                for j in (current_col + 1)..=c {
                    col_ptrs[j as usize] = len;
                }
                current_col = c;
            }
            row_indices.push(r);
            values.push(v);
        }
        let len = row_indices.len() as u32;
        for j in (current_col as usize + 1)..=ncols {
            col_ptrs[j] = len;
        }

        Ok(Self {
            nrows,
            ncols,
            col_ptrs,
            row_indices,
            values,
        })
    }

    /// Number of nonzeros.
    pub fn nnz(&self) -> usize {
        self.values.len()
    }

    /// Create a faer SparseColMatRef for use with faer-sparse APIs.
    pub fn as_faer_ref(&self) -> SparseColMatRef<'_, u32, f64> {
        let symbolic = SymbolicSparseColMatRef::new_checked(
            self.nrows,
            self.ncols,
            &self.col_ptrs,
            None,
            &self.row_indices,
        );
        SparseColMatRef::<_, f64>::new(symbolic, &*self.values)
    }

    /// Transpose: returns CSC of A^T (row-major of A = column-major of A^T).
    pub fn transpose(&self) -> Self {
        let m = self.nrows;
        let n = self.ncols;
        let mut triplets = Vec::with_capacity(self.nnz());
        for j in 0..n {
            let start = self.col_ptrs[j] as usize;
            let end = self.col_ptrs[j + 1] as usize;
            for idx in start..end {
                let i = self.row_indices[idx] as usize;
                triplets.push((j as u32, i as u32, self.values[idx]));
            }
        }
        Self::from_triplets(n, m, &triplets).expect("transpose from valid matrix")
    }

    /// Sparse matrix-vector product: y = A * x.
    pub fn matvec(&self, x: &[f64]) -> Vec<f64> {
        assert_eq!(x.len(), self.ncols);
        let mut y = vec![0.0; self.nrows];
        for j in 0..self.ncols {
            let start = self.col_ptrs[j] as usize;
            let end = self.col_ptrs[j + 1] as usize;
            let xj = x[j];
            for idx in start..end {
                let i = self.row_indices[idx] as usize;
                y[i] += self.values[idx] * xj;
            }
        }
        y
    }

    /// Sparse * sparse: C = A * B. A is m×k, B is k×n, C is m×n.
    pub fn sparse_times_sparse(a: &Self, b: &Self) -> Result<Self, String> {
        if a.ncols != b.nrows {
            return Err("Dimension mismatch for sparse * sparse".into());
        }
        let m = a.nrows;
        let n = b.ncols;

        let mut triplets = Vec::new();
        for j in 0..n {
            for b_idx in (b.col_ptrs[j] as usize)..(b.col_ptrs[j + 1] as usize) {
                let i_b = b.row_indices[b_idx] as usize;
                let val_b = b.values[b_idx];
                for a_idx in (a.col_ptrs[i_b] as usize)..(a.col_ptrs[i_b + 1] as usize) {
                    let i_a = a.row_indices[a_idx] as usize;
                    let val_a = a.values[a_idx];
                    triplets.push((i_a as u32, j as u32, val_a * val_b));
                }
            }
        }
        Self::from_triplets(m, n, &triplets)
    }

    /// Add value to diagonal: A += d * I. Modifies in place.
    pub fn add_diagonal(&mut self, d: f64) {
        let n = self.nrows.min(self.ncols);
        let mut missing = Vec::new();
        for j in 0..n {
            let start = self.col_ptrs[j] as usize;
            let end = self.col_ptrs[j + 1] as usize;
            match (start..end).find(|&idx| self.row_indices[idx] == j as u32) {
                Some(idx) => self.values[idx] += d,
                None => missing.push(j),
            }
        }
        if missing.is_empty() {
            return;
        }
        // Structurally empty diagonals (e.g. the Gram column of an edge whose
        // both ends are fixed) must receive the shift too, or a factorisation
        // meets an exact zero pivot.
        let mut triplets: Vec<(u32, u32, f64)> = Vec::with_capacity(self.nnz() + missing.len());
        for j in 0..self.ncols {
            for idx in self.col_ptrs[j] as usize..self.col_ptrs[j + 1] as usize {
                triplets.push((self.row_indices[idx], j as u32, self.values[idx]));
            }
        }
        for j in missing {
            triplets.push((j as u32, j as u32, d));
        }
        *self = Self::from_triplets(self.nrows, self.ncols, &triplets).expect("add_diagonal");
    }

    /// Extract selected columns into a new matrix.
    pub fn extract_columns(&self, col_indices: &[usize]) -> Self {
        let ncols = col_indices.len();
        let mut triplets = Vec::new();
        for (out_j, &in_j) in col_indices.iter().enumerate() {
            if in_j >= self.ncols {
                continue;
            }
            let start = self.col_ptrs[in_j] as usize;
            let end = self.col_ptrs[in_j + 1] as usize;
            for idx in start..end {
                let row = self.row_indices[idx];
                let val = self.values[idx];
                triplets.push((row, out_j as u32, val));
            }
        }
        Self::from_triplets(self.nrows, ncols, &triplets).expect("extract_columns")
    }

    /// Row of `nz` in column `col`, if that entry exists.
    pub fn nz_at(&self, row: usize, col: usize) -> Option<usize> {
        let start = self.col_ptrs[col] as usize;
        let end = self.col_ptrs[col + 1] as usize;
        (start..end).find(|&idx| self.row_indices[idx] as usize == row)
    }

    pub fn same_pattern(&self, other: &Self) -> bool {
        self.nrows == other.nrows
            && self.ncols == other.ncols
            && self.col_ptrs == other.col_ptrs
            && self.row_indices == other.row_indices
    }
}

/// `Aᵀ` whose pattern is fixed. Later calls overwrite values from `A`.
pub struct FixedTranspose {
    pub matrix: SparseColMatOwned,
    /// `at[src_nz]` is the slot in `matrix.values`.
    at: Vec<usize>,
}

impl FixedTranspose {
    pub fn new(src: &SparseColMatOwned) -> Self {
        let matrix = src.transpose();
        let mut at = vec![0; src.nnz()];
        for col in 0..src.ncols {
            for idx in src.col_ptrs[col] as usize..src.col_ptrs[col + 1] as usize {
                let row = src.row_indices[idx] as usize;
                at[idx] = matrix
                    .nz_at(col, row)
                    .expect("transpose contains every source nonzero");
            }
        }
        Self { matrix, at }
    }

    pub fn refresh(&mut self, src_values: &[f64]) {
        for (src, &dst) in self.at.iter().enumerate() {
            self.matrix.values[dst] = src_values[src];
        }
    }
}

/// Sparse product whose pattern was captured by [`SparseColMatOwned::sparse_times_sparse`].
///
/// Numeric refreshes use a Gustavson column accumulator. The first matrix is
/// the triplet product, so a caller can check the refresh against it.
pub struct FixedProduct {
    pub result: SparseColMatOwned,
    acc: Vec<f64>,
}

impl FixedProduct {
    pub fn new(a: &SparseColMatOwned, b: &SparseColMatOwned) -> Result<Self, String> {
        let result = SparseColMatOwned::sparse_times_sparse(a, b)?;
        let acc = vec![0.0; result.nrows];
        Ok(Self { result, acc })
    }

    pub fn refresh(&mut self, a: &SparseColMatOwned, b: &SparseColMatOwned) {
        for j in 0..b.ncols {
            let start = self.result.col_ptrs[j] as usize;
            let end = self.result.col_ptrs[j + 1] as usize;
            for idx in start..end {
                self.acc[self.result.row_indices[idx] as usize] = 0.0;
            }
            for b_idx in b.col_ptrs[j] as usize..b.col_ptrs[j + 1] as usize {
                let row_b = b.row_indices[b_idx] as usize;
                let val_b = b.values[b_idx];
                for a_idx in a.col_ptrs[row_b] as usize..a.col_ptrs[row_b + 1] as usize {
                    self.acc[a.row_indices[a_idx] as usize] += a.values[a_idx] * val_b;
                }
            }
            for idx in start..end {
                let row = self.result.row_indices[idx] as usize;
                self.result.values[idx] = self.acc[row];
            }
        }
    }
}

/// Sum of two matrices with a fixed union pattern.
pub struct FixedSum {
    pub result: SparseColMatOwned,
    left_at: Vec<usize>,
    right_at: Vec<usize>,
}

impl FixedSum {
    pub fn new(left: &SparseColMatOwned, right: &SparseColMatOwned) -> Result<Self, String> {
        let mut triplets = Vec::with_capacity(left.nnz() + right.nnz());
        for mat in [left, right] {
            for col in 0..mat.ncols {
                for nz in mat.col_ptrs[col] as usize..mat.col_ptrs[col + 1] as usize {
                    triplets.push((mat.row_indices[nz], col as u32, mat.values[nz]));
                }
            }
        }
        let result = SparseColMatOwned::from_triplets(left.nrows, left.ncols, &triplets)?;
        Ok(Self {
            left_at: locate_nonzeros(left, &result),
            right_at: locate_nonzeros(right, &result),
            result,
        })
    }

    pub fn refresh(&mut self, left: &SparseColMatOwned, right: &SparseColMatOwned) {
        self.result.values.fill(0.0);
        for (src, &dst) in self.left_at.iter().enumerate() {
            self.result.values[dst] += left.values[src];
        }
        for (src, &dst) in self.right_at.iter().enumerate() {
            self.result.values[dst] += right.values[src];
        }
    }
}

fn locate_nonzeros(src: &SparseColMatOwned, dst: &SparseColMatOwned) -> Vec<usize> {
    let mut at = Vec::with_capacity(src.nnz());
    for col in 0..src.ncols {
        for idx in src.col_ptrs[col] as usize..src.col_ptrs[col + 1] as usize {
            at.push(
                dst.nz_at(src.row_indices[idx] as usize, col)
                    .expect("sum pattern covers both terms"),
            );
        }
    }
    at
}

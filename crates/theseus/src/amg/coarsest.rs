//! Coarsest-level direct solve: faer sparse LLᵀ through
//! [`types::Factorization`](crate::types::Factorization) on the coarsest
//! [`LevelMatrix`], refactored numerically on every `update(q)` (the
//! symbolic analysis is reused) and applied to the three columns with the
//! specialised triangular solve of `factor_solve`.

use super::hierarchy::LevelMatrix;
use crate::sparse::SparseColMatOwned;
use crate::types::{Factorization, FactorizationStrategy, TheseusError};

/// Factorization of the coarsest operator and its host-side copies.
pub struct CoarsestSolver {
    matrix: SparseColMatOwned,
    factorization: Factorization,
    stack: dyn_stack::GlobalPodBuffer,
    n: usize,
}

impl std::fmt::Debug for CoarsestSolver {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("CoarsestSolver")
            .field("n", &self.n)
            .field("nnz", &self.matrix.nnz())
            .field("factorization", &self.factorization)
            .finish()
    }
}

impl CoarsestSolver {
    /// Symbolic + numeric Cholesky of `a` (square SPD).
    pub fn new(a: &LevelMatrix) -> Result<Self, TheseusError> {
        assert!(a.is_square(), "coarsest matrix must be square");
        let matrix = a.to_sparse();
        let mut stack = dyn_stack::GlobalPodBuffer::new(dyn_stack::StackReq::empty());
        let factorization =
            Factorization::new(&matrix, FactorizationStrategy::Cholesky, &mut stack)?;
        Ok(Self {
            matrix,
            factorization,
            stack,
            n: a.n,
        })
    }

    /// Numeric refactorization with new values on the same pattern.
    pub fn update(&mut self, a: &LevelMatrix) -> Result<(), TheseusError> {
        assert_eq!(a.n, self.n, "coarsest matrix size changed");
        a.write_sparse_values(&mut self.matrix);
        self.factorization.update(&self.matrix, &mut self.stack)
    }

    /// Nodes on the coarsest level.
    pub fn n(&self) -> usize {
        self.n
    }

    /// Length the `work` slice of [`Self::solve`] must have.
    pub fn work_len(&self) -> usize {
        6 * self.n
    }

    /// `x = A⁻¹ b` on three interleaved columns (`b`, `x`: `n * 3`; `work`:
    /// `6 n`). Allocation-free.
    pub fn solve(&self, b: &[f64], x: &mut [f64], work: &mut [f64]) {
        debug_assert_eq!(b.len(), self.n * 3);
        debug_assert_eq!(x.len(), self.n * 3);
        debug_assert!(work.len() >= self.work_len());
        self.factorization.solve_slices::<3>(b, x, work);
    }

    /// Host bytes held (matrix copy, factor values, faer stack).
    pub fn bytes(&self) -> usize {
        self.matrix.col_ptrs.len() * 4
            + self.matrix.row_indices.len() * 4
            + self.matrix.values.len() * 8
            + self.factorization.len_values() * 8
            + self.stack.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::{CsrAdjacency, LevelGraph};

    #[test]
    fn coarsest_solve_and_refactor_match_the_operator() {
        // Small anchored grid as an explicit CSR matrix.
        let side = 6;
        let n = side * side;
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
        let mut g = LevelGraph {
            n,
            adjacency: CsrAdjacency::from_endpoints(n, &starts, &ends),
            weight: (0..starts.len())
                .map(|k| 1.0 + 0.1 * (k % 4) as f64)
                .collect(),
            anchor: (0..n)
                .map(|u| if u % side == 0 { 1.0 } else { 0.0 })
                .collect(),
            aggregate_of: Vec::new(),
        };
        let mut a = LevelMatrix::from_graph(&g);
        let mut solver = CoarsestSolver::new(&a).unwrap();
        assert_eq!(solver.n(), n);
        let b: Vec<f64> = (0..n * 3).map(|i| ((i * 7 % 11) as f64) - 5.0).collect();
        let mut x = vec![0.0; n * 3];
        let mut work = vec![0.0; solver.work_len()];
        solver.solve(&b, &mut x, &mut work);
        let mut ax = vec![0.0; n * 3];
        a.apply(&x, &mut ax, 3);
        for (l, r) in ax.iter().zip(&b) {
            assert!((l - r).abs() < 1e-10, "{l} vs {r}");
        }
        // Refactor with new weights.
        g.weight.iter_mut().for_each(|w| *w *= 3.0);
        g.anchor[0] = 5.0;
        a.set_from_graph(&g);
        solver.update(&a).unwrap();
        solver.solve(&b, &mut x, &mut work);
        a.apply(&x, &mut ax, 3);
        for (l, r) in ax.iter().zip(&b) {
            assert!((l - r).abs() < 1e-10, "{l} vs {r}");
        }
        assert!(solver.bytes() > 0);
    }

    #[test]
    fn singular_coarsest_is_an_error() {
        // Two nodes, one edge, no anchors: A is singular.
        let g = LevelGraph {
            n: 2,
            adjacency: CsrAdjacency::from_endpoints(2, &[0], &[1]),
            weight: vec![1.0],
            anchor: vec![0.0, 0.0],
            aggregate_of: Vec::new(),
        };
        let a = LevelMatrix::from_graph(&g);
        assert!(CoarsestSolver::new(&a).is_err());
    }
}

//! Reusable solve storage.

use crate::WorkspaceError;

/// Reusable allocation owned by a [`Solver`](crate::Solver).
///
/// A workspace grows when needed and retains its buffers for later solves.
#[derive(Debug, Default)]
pub struct Workspace {
    pub(crate) lower: Vec<f64>,
    pub(crate) upper: Vec<f64>,
    pub(crate) evaluation_gradient: Vec<f64>,
    pub(crate) final_gradient: Vec<f64>,
    pub(crate) old_x: Vec<f64>,
    pub(crate) old_gradient: Vec<f64>,
    pub(crate) direction: Vec<f64>,
    pub(crate) cauchy: Vec<f64>,
    pub(crate) displacement: Vec<f64>,
    pub(crate) product: Vec<f64>,
    pub(crate) row: Vec<f64>,
    pub(crate) reduced_gradient: Vec<f64>,
    pub(crate) trial_gradient: Vec<f64>,
    pub(crate) breakpoints: Vec<f64>,
    pub(crate) order: Vec<usize>,
    pub(crate) status: Vec<i8>,
    pub(crate) s: Vec<f64>,
    pub(crate) y: Vec<f64>,
    pub(crate) sy: Vec<f64>,
    pub(crate) ss: Vec<f64>,
    pub(crate) wt: Vec<f64>,
    pub(crate) wn: Vec<f64>,
    pub(crate) wn1: Vec<f64>,
    pub(crate) wa: Vec<f64>,
    pub(crate) index: Vec<usize>,
    pub(crate) changes: Vec<usize>,
}

impl Workspace {
    /// Creates an empty workspace that grows on its first solve.
    pub fn new() -> Self {
        Self::default()
    }

    /// Creates a workspace with per-variable buffers reserved in advance.
    ///
    /// History-dependent buffers are sized when the workspace is first used.
    pub fn with_capacity(variables: usize) -> Result<Self, WorkspaceError> {
        let mut workspace = Self::new();
        reserve(&mut workspace.lower, variables, "lower")?;
        reserve(&mut workspace.upper, variables, "upper")?;
        reserve(
            &mut workspace.evaluation_gradient,
            variables,
            "evaluation_gradient",
        )?;
        reserve(&mut workspace.final_gradient, variables, "final_gradient")?;
        reserve(&mut workspace.old_x, variables, "old_x")?;
        reserve(&mut workspace.old_gradient, variables, "old_gradient")?;
        reserve(&mut workspace.direction, variables, "direction")?;
        reserve(&mut workspace.cauchy, variables, "cauchy")?;
        reserve(&mut workspace.displacement, variables, "displacement")?;
        reserve(&mut workspace.product, variables, "product")?;
        reserve(&mut workspace.row, variables, "row")?;
        reserve(
            &mut workspace.reduced_gradient,
            variables,
            "reduced_gradient",
        )?;
        reserve(&mut workspace.trial_gradient, variables, "trial_gradient")?;
        reserve(&mut workspace.breakpoints, variables, "breakpoints")?;
        reserve(&mut workspace.order, variables, "order")?;
        reserve(&mut workspace.status, variables, "status")?;
        reserve(&mut workspace.index, variables, "index")?;
        reserve(&mut workspace.changes, variables, "changes")?;
        Ok(workspace)
    }

    /// Number of variables for which the workspace currently holds a result.
    pub fn dimension(&self) -> usize {
        self.final_gradient.len()
    }

    /// Gradient at the final point of the most recent successful solve.
    pub fn gradient(&self) -> &[f64] {
        &self.final_gradient
    }

    pub(crate) fn prepare(
        &mut self,
        lower: &[f64],
        upper: &[f64],
        history: usize,
    ) -> Result<(), WorkspaceError> {
        let n = lower.len();
        let n_history = n.checked_mul(history).ok_or(WorkspaceError::SizeOverflow {
            variables: n,
            history_size: history,
        })?;
        let history_squared = history
            .checked_mul(history)
            .ok_or(WorkspaceError::SizeOverflow {
                variables: n,
                history_size: history,
            })?;
        let four_history_squared =
            history_squared
                .checked_mul(4)
                .ok_or(WorkspaceError::SizeOverflow {
                    variables: n,
                    history_size: history,
                })?;
        let eight_history = history.checked_mul(8).ok_or(WorkspaceError::SizeOverflow {
            variables: n,
            history_size: history,
        })?;

        reserve(&mut self.lower, n, "lower")?;
        reserve(&mut self.upper, n, "upper")?;
        reserve(&mut self.evaluation_gradient, n, "evaluation_gradient")?;
        reserve(&mut self.trial_gradient, n, "trial_gradient")?;
        reserve(&mut self.old_x, n, "old_x")?;
        reserve(&mut self.old_gradient, n, "old_gradient")?;
        reserve(&mut self.direction, n, "direction")?;
        reserve(&mut self.cauchy, n, "cauchy")?;
        reserve(&mut self.displacement, n, "displacement")?;
        reserve(&mut self.product, n, "product")?;
        reserve(&mut self.row, n, "row")?;
        reserve(&mut self.reduced_gradient, n, "reduced_gradient")?;
        reserve(&mut self.breakpoints, n, "breakpoints")?;
        reserve(&mut self.order, n, "order")?;
        reserve(&mut self.status, n, "status")?;
        reserve(&mut self.s, n_history, "s")?;
        reserve(&mut self.y, n_history, "y")?;
        reserve(&mut self.sy, history_squared, "sy")?;
        reserve(&mut self.ss, history_squared, "ss")?;
        reserve(&mut self.wt, history_squared, "wt")?;
        reserve(&mut self.wn, four_history_squared, "wn")?;
        reserve(&mut self.wn1, four_history_squared, "wn1")?;
        reserve(&mut self.wa, eight_history, "wa")?;
        reserve(&mut self.index, n, "index")?;
        reserve(&mut self.changes, n, "changes")?;

        self.lower.clear();
        self.lower.extend_from_slice(lower);
        self.upper.clear();
        self.upper.extend_from_slice(upper);
        self.evaluation_gradient.resize(n, 0.0);
        self.trial_gradient.resize(n, 0.0);
        self.old_x.resize(n, 0.0);
        self.old_gradient.resize(n, 0.0);
        self.direction.resize(n, 0.0);
        self.cauchy.resize(n, 0.0);
        self.displacement.resize(n, 0.0);
        self.product.resize(n, 0.0);
        self.row.resize(n, 0.0);
        self.reduced_gradient.resize(n, 0.0);
        self.breakpoints.resize(n, 0.0);
        self.order.resize(n, 0);
        self.status.resize(n, 0);
        self.s.resize(n_history, 0.0);
        self.y.resize(n_history, 0.0);
        self.sy.resize(history_squared, 0.0);
        self.ss.resize(history_squared, 0.0);
        self.wt.resize(history_squared, 0.0);
        self.wn.resize(four_history_squared, 0.0);
        self.wn1.resize(four_history_squared, 0.0);
        self.wa.resize(eight_history, 0.0);
        self.index.resize(n, 0);
        self.changes.resize(n, 0);
        Ok(())
    }

    pub(crate) fn commit_gradient(&mut self) -> Result<(), WorkspaceError> {
        reserve(
            &mut self.final_gradient,
            self.evaluation_gradient.len(),
            "final_gradient",
        )?;
        self.final_gradient.clear();
        self.final_gradient
            .extend_from_slice(&self.evaluation_gradient);
        Ok(())
    }
}

fn reserve<T>(
    buffer: &mut Vec<T>,
    elements: usize,
    name: &'static str,
) -> Result<(), WorkspaceError> {
    let additional = elements.saturating_sub(buffer.len());
    buffer
        .try_reserve_exact(additional)
        .map_err(|_| WorkspaceError::AllocationFailed {
            buffer: name,
            elements,
        })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn size_overflow_is_typed() {
        let mut workspace = Workspace::new();
        assert!(matches!(
            workspace.prepare(&[0.0], &[1.0], usize::MAX),
            Err(WorkspaceError::SizeOverflow { .. })
        ));
    }

    #[test]
    fn capacity_failure_is_typed() {
        assert!(matches!(
            Workspace::with_capacity(usize::MAX),
            Err(WorkspaceError::AllocationFailed { .. })
        ));
    }
}

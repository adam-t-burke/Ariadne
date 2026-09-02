//! Private dense-history kernels.  Solver control flow is backend-independent.

use crate::Backend;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Kernel {
    Scalar,
    #[cfg(feature = "faer-backend")]
    Faer,
}

pub(crate) struct History<'a> {
    pub(crate) n: usize,
    pub(crate) m: usize,
    pub(crate) head: usize,
    pub(crate) col: usize,
    pub(crate) s: &'a [f64],
    pub(crate) y: &'a [f64],
}

impl Kernel {
    pub(crate) fn selected(backend: Backend) -> Self {
        match backend {
            Backend::Deterministic | Backend::Auto => Self::Scalar,
            #[cfg(feature = "faer-backend")]
            Backend::Faer => Self::Faer,
            #[cfg(not(feature = "faer-backend"))]
            Backend::Faer => unreachable!("unsupported backend rejected by Options"),
        }
    }

    pub(crate) fn public(self) -> Backend {
        match self {
            Self::Scalar => Backend::Deterministic,
            #[cfg(feature = "faer-backend")]
            Self::Faer => Backend::Faer,
        }
    }

    /// Computes dot products of the newest S vector against all logical Y and
    /// S history vectors. History is split at the physical ring boundary, so
    /// neither backend repacks it.
    pub(crate) fn newest_products(
        self,
        history: History<'_>,
        sy_out: &mut [f64],
        ss_out: &mut [f64],
    ) {
        let History {
            n,
            m,
            head,
            col,
            s,
            y,
        } = history;
        debug_assert!(col > 0);
        let newest = (head + col - 1) % m;
        let newest_s = &s[newest * n..(newest + 1) * n];
        match self {
            Self::Scalar => {
                for logical in 0..col {
                    let physical = (head + logical) % m;
                    sy_out[logical] = ordered_dot(newest_s, &y[physical * n..(physical + 1) * n]);
                    ss_out[logical] = ordered_dot(&s[physical * n..(physical + 1) * n], newest_s);
                }
            }
            #[cfg(feature = "faer-backend")]
            Self::Faer => {
                let first = col.min(m - head);
                faer_products(
                    n,
                    &s[head * n..(head + first) * n],
                    &y[head * n..(head + first) * n],
                    newest_s,
                    &mut sy_out[..first],
                    &mut ss_out[..first],
                );
                if first < col {
                    let second = col - first;
                    faer_products(
                        n,
                        &s[..second * n],
                        &y[..second * n],
                        newest_s,
                        &mut sy_out[first..col],
                        &mut ss_out[first..col],
                    );
                }
            }
        }
    }
}

#[inline]
fn ordered_dot(a: &[f64], b: &[f64]) -> f64 {
    let mut sum = 0.0;
    for i in 0..a.len() {
        sum += a[i] * b[i];
    }
    sum
}

#[cfg(feature = "faer-backend")]
fn faer_products(
    n: usize,
    s: &[f64],
    y: &[f64],
    newest_s: &[f64],
    sy_out: &mut [f64],
    ss_out: &mut [f64],
) {
    use faer_core::{mat, mul::matmul, Parallelism};

    let cols = sy_out.len();
    let s_mat = mat::from_column_major_slice::<f64>(s, n, cols);
    let y_mat = mat::from_column_major_slice::<f64>(y, n, cols);
    let vector = mat::from_column_major_slice::<f64>(newest_s, n, 1);
    let sy = mat::from_column_major_slice_mut::<f64>(sy_out, cols, 1);
    let ss = mat::from_column_major_slice_mut::<f64>(ss_out, cols, 1);
    matmul(sy, y_mat.transpose(), vector, None, 1.0, Parallelism::None);
    matmul(ss, s_mat.transpose(), vector, None, 1.0, Parallelism::None);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn scalar_products_follow_wrapped_logical_order() {
        let n = 3;
        let m = 4;
        let s = [
            1.0, 2.0, 3.0, // physical 0, newest
            99.0, 99.0, 99.0, // unused
            2.0, 0.0, 1.0, // logical 0
            0.0, 3.0, 1.0, // logical 1
        ];
        let y = [
            2.0, 1.0, 0.0, // physical 0
            99.0, 99.0, 99.0, // unused
            1.0, 1.0, 1.0, // logical 0
            0.0, 1.0, 2.0, // logical 1
        ];
        let mut sy = [0.0; 3];
        let mut ss = [0.0; 3];
        Kernel::Scalar.newest_products(
            History {
                n,
                m,
                head: 2,
                col: 3,
                s: &s,
                y: &y,
            },
            &mut sy,
            &mut ss,
        );
        assert_eq!(sy, [6.0, 8.0, 4.0]);
        assert_eq!(ss, [5.0, 9.0, 14.0]);
    }

    #[cfg(feature = "faer-backend")]
    #[test]
    fn faer_matches_scalar_across_ring_wrap() {
        let n = 513;
        let m = 10;
        let s: Vec<_> = (0..n * m).map(|i| (i % 37) as f64 - 18.0).collect();
        let y: Vec<_> = (0..n * m).map(|i| (i % 29) as f64 - 14.0).collect();
        let mut scalar_sy = [0.0; 10];
        let mut scalar_ss = [0.0; 10];
        let mut faer_sy = [0.0; 10];
        let mut faer_ss = [0.0; 10];
        let history = || History {
            n,
            m,
            head: 7,
            col: 10,
            s: &s,
            y: &y,
        };
        Kernel::Scalar.newest_products(history(), &mut scalar_sy, &mut scalar_ss);
        Kernel::Faer.newest_products(history(), &mut faer_sy, &mut faer_ss);
        for (actual, expected) in faer_sy
            .iter()
            .chain(&faer_ss)
            .zip(scalar_sy.iter().chain(&scalar_ss))
        {
            assert!((actual - expected).abs() <= 2.0e-12 * expected.abs().max(1.0));
        }
    }
}

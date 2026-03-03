use crate::types::{KnotVector, NurbsError};

/// Algorithm A2.1: determine the knot span index containing `u`.
///
/// `num_control_points` is `n + 1` in the book's notation.
pub fn find_span(
    num_control_points: usize,
    degree: usize,
    u: f64,
    knots: &KnotVector,
) -> Result<usize, NurbsError> {
    let u_knots = knots.as_slice();
    knots.domain(degree, num_control_points)?;

    let n = num_control_points - 1;
    if u >= u_knots[n + 1] {
        return Ok(n);
    }
    if u <= u_knots[degree] {
        return Ok(degree);
    }

    let mut low = degree;
    let mut high = n + 1;
    let mut mid = (low + high) / 2;
    while u < u_knots[mid] || u >= u_knots[mid + 1] {
        if u < u_knots[mid] {
            high = mid;
        } else {
            low = mid;
        }
        mid = (low + high) / 2;
    }
    Ok(mid)
}

/// Algorithm A2.2: evaluate the nonzero basis functions for a span.
pub fn basis_functions(span: usize, u: f64, degree: usize, knots: &KnotVector) -> Vec<f64> {
    let u_knots = knots.as_slice();
    let mut n = vec![0.0; degree + 1];
    let mut left = vec![0.0; degree + 1];
    let mut right = vec![0.0; degree + 1];
    n[0] = 1.0;

    for j in 1..=degree {
        left[j] = u - u_knots[span + 1 - j];
        right[j] = u_knots[span + j] - u;
        let mut saved = 0.0;
        for r in 0..j {
            let denom = right[r + 1] + left[j - r];
            let temp = if denom.abs() > f64::EPSILON {
                n[r] / denom
            } else {
                0.0
            };
            n[r] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        n[j] = saved;
    }
    n
}

/// Algorithm A2.3: derivatives of the nonzero basis functions through order `derivative_order`.
pub fn basis_function_derivatives(
    span: usize,
    u: f64,
    degree: usize,
    derivative_order: usize,
    knots: &KnotVector,
) -> Vec<Vec<f64>> {
    let u_knots = knots.as_slice();
    let n = derivative_order.min(degree);
    let mut ndu = vec![vec![0.0; degree + 1]; degree + 1];
    let mut left = vec![0.0; degree + 1];
    let mut right = vec![0.0; degree + 1];
    ndu[0][0] = 1.0;

    for j in 1..=degree {
        left[j] = u - u_knots[span + 1 - j];
        right[j] = u_knots[span + j] - u;
        let mut saved = 0.0;
        for r in 0..j {
            ndu[j][r] = right[r + 1] + left[j - r];
            let temp = if ndu[j][r].abs() > f64::EPSILON {
                ndu[r][j - 1] / ndu[j][r]
            } else {
                0.0
            };
            ndu[r][j] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        ndu[j][j] = saved;
    }

    let mut ders = vec![vec![0.0; degree + 1]; derivative_order + 1];
    for j in 0..=degree {
        ders[0][j] = ndu[j][degree];
    }

    let mut a = vec![vec![0.0; degree + 1]; 2];
    for r in 0..=degree {
        let mut s1 = 0usize;
        let mut s2 = 1usize;
        a[0][0] = 1.0;

        for k in 1..=n {
            let mut d = 0.0;
            let rk = r as isize - k as isize;
            let pk = degree as isize - k as isize;

            if r >= k {
                a[s2][0] = a[s1][0] / ndu[(pk + 1) as usize][rk as usize];
                d = a[s2][0] * ndu[rk as usize][pk as usize];
            }

            let j1 = if rk >= -1 { 1 } else { (-rk) as usize };
            let j2 = if (r as isize - 1) <= pk {
                k - 1
            } else {
                degree - r
            };

            for j in j1..=j2 {
                a[s2][j] =
                    (a[s1][j] - a[s1][j - 1]) / ndu[(pk + 1) as usize][(rk + j as isize) as usize];
                d += a[s2][j] * ndu[(rk + j as isize) as usize][pk as usize];
            }

            if r <= degree - k {
                a[s2][k] = -a[s1][k - 1] / ndu[(pk + 1) as usize][r];
                d += a[s2][k] * ndu[r][pk as usize];
            }

            ders[k][r] = d;
            std::mem::swap(&mut s1, &mut s2);
        }
    }

    let mut factor = degree as f64;
    for k in 1..=n {
        for j in 0..=degree {
            ders[k][j] *= factor;
        }
        factor *= (degree - k) as f64;
    }
    ders
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cubic_bezier_basis_sums_to_one() {
        let knots = KnotVector::new(vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]).unwrap();
        let span = find_span(4, 3, 0.5, &knots).unwrap();
        let n = basis_functions(span, 0.5, 3, &knots);
        let sum: f64 = n.iter().sum();
        assert!((sum - 1.0).abs() < 1e-12);
        assert!((n[0] - 0.125).abs() < 1e-12);
        assert!((n[1] - 0.375).abs() < 1e-12);
        assert!((n[2] - 0.375).abs() < 1e-12);
        assert!((n[3] - 0.125).abs() < 1e-12);
    }
}

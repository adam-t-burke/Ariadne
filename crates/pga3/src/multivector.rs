//! G(3,0,1) multivectors with a bitwise geometric product.
//!
//! Basis vectors are `e0, e1, e2, e3` with metric `(0, 1, 1, 1)`. Coefficients
//! are stored in bitmask order: bit `i` set means factor `e_i`, and factors are
//! always written in increasing index (so `e13` not `e31`).

use std::ops::{Add, Mul, Neg, Sub};

pub const BASIS: usize = 16;

/// Bitmask blades: `e0=1`, `e1=2`, `e2=4`, `e3=8`.
pub const E0: usize = 1;
pub const E1: usize = 2;
pub const E2: usize = 4;
pub const E3: usize = 8;
pub const E01: usize = 3;
pub const E02: usize = 5;
pub const E03: usize = 9;
pub const E12: usize = 6;
pub const E13: usize = 10;
pub const E23: usize = 12;
pub const E012: usize = 7;
pub const E013: usize = 11;
pub const E023: usize = 13;
pub const E123: usize = 14;
pub const E0123: usize = 15;

/// 16-coefficient projective multivector.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MultiVector {
    pub c: [f64; BASIS],
}

impl Default for MultiVector {
    fn default() -> Self {
        Self::zero()
    }
}

impl MultiVector {
    pub const fn zero() -> Self {
        Self { c: [0.0; BASIS] }
    }

    pub const fn scalar(s: f64) -> Self {
        let mut c = [0.0; BASIS];
        c[0] = s;
        Self { c }
    }

    pub fn from_blade(blade: usize, coeff: f64) -> Self {
        let mut mv = Self::zero();
        mv.c[blade] = coeff;
        mv
    }

    pub fn grade(&self, k: u32) -> Self {
        let mut out = Self::zero();
        for (i, &coeff) in self.c.iter().enumerate() {
            if (i as u8).count_ones() == k {
                out.c[i] = coeff;
            }
        }
        out
    }

    pub fn reverse(self) -> Self {
        let mut out = Self::zero();
        for (i, &coeff) in self.c.iter().enumerate() {
            out.c[i] = coeff * reverse_sign(i as u8) as f64;
        }
        out
    }

    pub fn complement(self) -> Self {
        let mut out = Self::zero();
        for (i, &coeff) in self.c.iter().enumerate() {
            if coeff == 0.0 {
                continue;
            }
            let (dual, sign) = COMPLEMENT[i];
            out.c[dual as usize] += sign as f64 * coeff;
        }
        out
    }

    pub fn geometric(self, other: Self) -> Self {
        product(self, other, &GP)
    }

    pub fn outer(self, other: Self) -> Self {
        product(self, other, &OP)
    }

    pub fn left_contract(self, other: Self) -> Self {
        product(self, other, &LC)
    }

    /// Poincaré-dual regressive product `A ∨ B = (A* ∧ B*)*`.
    pub fn regressive(self, other: Self) -> Self {
        self.complement().outer(other.complement()).complement()
    }

    pub fn sandwich(self, x: Self) -> Self {
        self.geometric(x).geometric(self.reverse())
    }

    pub fn norm_sq(self) -> f64 {
        self.c.iter().map(|v| v * v).sum()
    }

    pub fn approx_eq(self, other: Self, tol: f64) -> bool {
        self.c
            .iter()
            .zip(other.c.iter())
            .all(|(a, b)| (a - b).abs() <= tol)
    }
}

impl Add for MultiVector {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        let mut c = [0.0; BASIS];
        for i in 0..BASIS {
            c[i] = self.c[i] + rhs.c[i];
        }
        Self { c }
    }
}

impl Sub for MultiVector {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        let mut c = [0.0; BASIS];
        for i in 0..BASIS {
            c[i] = self.c[i] - rhs.c[i];
        }
        Self { c }
    }
}

impl Neg for MultiVector {
    type Output = Self;
    fn neg(self) -> Self {
        let mut c = [0.0; BASIS];
        for i in 0..BASIS {
            c[i] = -self.c[i];
        }
        Self { c }
    }
}

impl Mul<f64> for MultiVector {
    type Output = Self;
    fn mul(self, rhs: f64) -> Self {
        let mut c = [0.0; BASIS];
        for i in 0..BASIS {
            c[i] = self.c[i] * rhs;
        }
        Self { c }
    }
}

fn product(a: MultiVector, b: MultiVector, table: &[[(u8, i8); BASIS]; BASIS]) -> MultiVector {
    let mut out = MultiVector::zero();
    for i in 0..BASIS {
        if a.c[i] == 0.0 {
            continue;
        }
        for j in 0..BASIS {
            if b.c[j] == 0.0 {
                continue;
            }
            let (blade, sign) = table[i][j];
            if sign == 0 {
                continue;
            }
            out.c[blade as usize] += sign as f64 * a.c[i] * b.c[j];
        }
    }
    out
}

const fn gp_basis(a: u8, b: u8) -> (u8, i8) {
    let mut sign: i8 = 1;
    let result = a ^ b;
    let mut i = 0;
    while i < 4 {
        let bit = 1u8 << i;
        if a & bit != 0 {
            let lower = bit - 1;
            if (b & lower).count_ones() % 2 == 1 {
                sign = -sign;
            }
            if b & bit != 0 {
                if i == 0 {
                    return (0, 0);
                }
            }
        }
        i += 1;
    }
    (result, sign)
}

const fn op_basis(a: u8, b: u8) -> (u8, i8) {
    if a & b != 0 {
        (0, 0)
    } else {
        gp_basis(a, b)
    }
}

const fn lc_basis(a: u8, b: u8) -> (u8, i8) {
    let (result, sign) = gp_basis(a, b);
    if sign == 0 {
        return (0, 0);
    }
    let ga = a.count_ones();
    let gb = b.count_ones();
    let gr = result.count_ones();
    if gb >= ga && gr == gb - ga {
        (result, sign)
    } else {
        (0, 0)
    }
}

const fn reverse_sign(bits: u8) -> i8 {
    let k = bits.count_ones();
    if (k * (k.wrapping_sub(1)) / 2) % 2 == 1 {
        -1
    } else {
        1
    }
}

const fn complement_basis(a: u8) -> (u8, i8) {
    let dual = (!a) & 0b1111;
    let (_, sign) = op_basis(a, dual);
    (dual, sign)
}

const fn make_table(kind: u8) -> [[(u8, i8); BASIS]; BASIS] {
    let mut table = [[(0u8, 0i8); BASIS]; BASIS];
    let mut i = 0;
    while i < BASIS {
        let mut j = 0;
        while j < BASIS {
            table[i][j] = match kind {
                0 => gp_basis(i as u8, j as u8),
                1 => op_basis(i as u8, j as u8),
                _ => lc_basis(i as u8, j as u8),
            };
            j += 1;
        }
        i += 1;
    }
    table
}

const GP: [[(u8, i8); BASIS]; BASIS] = make_table(0);
const OP: [[(u8, i8); BASIS]; BASIS] = make_table(1);
const LC: [[(u8, i8); BASIS]; BASIS] = make_table(2);

const COMPLEMENT: [(u8, i8); BASIS] = {
    let mut table = [(0u8, 0i8); BASIS];
    let mut i = 0;
    while i < BASIS {
        table[i] = complement_basis(i as u8);
        i += 1;
    }
    table
};

#[cfg(test)]
mod tests {
    use super::*;

    fn e(blade: usize) -> MultiVector {
        MultiVector::from_blade(blade, 1.0)
    }

    #[test]
    fn metric_squares() {
        assert!(e(E0).geometric(e(E0)).approx_eq(MultiVector::zero(), 0.0));
        assert!(e(E1)
            .geometric(e(E1))
            .approx_eq(MultiVector::scalar(1.0), 1e-15));
        assert!(e(E2)
            .geometric(e(E2))
            .approx_eq(MultiVector::scalar(1.0), 1e-15));
        assert!(e(E3)
            .geometric(e(E3))
            .approx_eq(MultiVector::scalar(1.0), 1e-15));
    }

    #[test]
    fn distinct_vectors_anticommute() {
        let e01 = e(E0).geometric(e(E1));
        let e10 = e(E1).geometric(e(E0));
        assert!(e01.approx_eq(-e10, 1e-15));
        assert!(
            (e(E1).geometric(e(E2)) + e(E2).geometric(e(E1))).approx_eq(MultiVector::zero(), 1e-15)
        );
    }

    #[test]
    fn bivector_square() {
        let e12 = e(E1).geometric(e(E2));
        assert!(e12
            .geometric(e12)
            .approx_eq(MultiVector::scalar(-1.0), 1e-15));
    }

    #[test]
    fn geometric_product_associative() {
        let a = e(E1) + e(E01) * 0.5 + MultiVector::from_blade(E123, 2.0);
        let b = e(E2) + e(E0) * -1.0 + e(E23);
        let c = e(E3) + MultiVector::scalar(0.25) + e(E0123);
        let left = a.geometric(b).geometric(c);
        let right = a.geometric(b.geometric(c));
        assert!(left.approx_eq(right, 1e-12));
    }

    #[test]
    fn reverse_involution() {
        let x = e(E12) + e(E123) + e(E01) + MultiVector::scalar(3.0);
        assert!(x.reverse().reverse().approx_eq(x, 0.0));
    }

    #[test]
    fn complement_recovers_pseudoscalar() {
        for i in 0..BASIS {
            let a = e(i);
            let star = a.complement();
            let meet = a.outer(star);
            assert!(meet.approx_eq(e(E0123), 1e-15), "blade {i}: {:?}", meet.c);
        }
    }
}

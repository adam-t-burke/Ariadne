/// Return `n choose k` as `f64`.
pub fn binomial(n: usize, k: usize) -> f64 {
    if k > n {
        return 0.0;
    }
    let k = k.min(n - k);
    let mut result = 1.0;
    for i in 1..=k {
        result *= (n - k + i) as f64;
        result /= i as f64;
    }
    result
}

pub fn add_scaled3(acc: &mut [f64; 3], scale: f64, p: [f64; 3]) {
    acc[0] += scale * p[0];
    acc[1] += scale * p[1];
    acc[2] += scale * p[2];
}

pub fn add_scaled4(acc: &mut [f64; 4], scale: f64, p: [f64; 4]) {
    acc[0] += scale * p[0];
    acc[1] += scale * p[1];
    acc[2] += scale * p[2];
    acc[3] += scale * p[3];
}

pub fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

pub fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

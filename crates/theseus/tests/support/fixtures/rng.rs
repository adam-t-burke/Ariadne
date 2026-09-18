//! Small deterministic PRNG (SplitMix64 seeding, xorshift64* stream) so the
//! fixtures are reproducible without pulling a random-number crate into the
//! dev-dependencies.

#[derive(Debug, Clone)]
pub struct Rng {
    state: u64,
}

impl Rng {
    /// Seeds the generator; any seed (including 0) produces a usable stream.
    pub fn new(seed: u64) -> Self {
        // SplitMix64 finaliser on the seed so nearby seeds give unrelated streams.
        let mut z = seed.wrapping_add(0x9E37_79B9_7F4A_7C15);
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^= z >> 31;
        Self { state: z | 1 }
    }

    pub fn next_u64(&mut self) -> u64 {
        // xorshift64*
        let mut x = self.state;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.state = x;
        x.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }

    /// Uniform in `[0, 1)` with 53 bits of resolution.
    pub fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
    }

    /// Uniform in `[-1, 1)`.
    pub fn next_signed(&mut self) -> f64 {
        2.0 * self.next_f64() - 1.0
    }

    /// Uniform integer in `lo..=hi`.
    pub fn next_range(&mut self, lo: usize, hi: usize) -> usize {
        debug_assert!(lo <= hi);
        let span = (hi - lo + 1) as u64;
        lo + (self.next_u64() % span) as usize
    }
}

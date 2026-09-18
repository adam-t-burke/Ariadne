//! Compute backends for the iterative solvers: the [`Backend`] kernel set
//! over an abstract buffer type, and the [`BackendHandle`] markers.
//!
//! The AMG preconditioner and the outer flexible-CG iteration are written
//! once against [`Backend`]; the CPU implementation (`backend/cpu.rs`, rayon
//! kernels with fixed-order reductions) and the GPU implementation
//! (`backend/gpu/`, wgpu + WGSL) are supplied by later workstreams. Every
//! kernel operates on **blocks of three columns** (`x`, `y`, `z` right-hand
//! sides) stored row-major, `v[node * 3 + k]`; per-column scalars are passed
//! as `[f64; 3]`.
//!
//! Buffers are allocated once per topology by the solver and reused; no
//! kernel allocates. Scalars that leave the device (`dot3`, `norm3`) are
//! returned as `f64` regardless of the buffer precision.

use crate::graph::LevelGraph;
use crate::linear_solver::Precision;

#[cfg(feature = "gpu")]
pub mod gpu;

/// One adapter seen by [`probe_gpu`], with the limits the solver cares about.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct GpuAdapterReport {
    /// Adapter (device) name as reported by the driver.
    pub name: String,
    /// wgpu backend: `vulkan`, `dx12`, `metal`, `gl`, ….
    pub backend: String,
    /// `discrete`, `integrated`, `virtual`, `cpu` or `other`.
    pub device_type: String,
    /// Driver name and version string.
    pub driver: String,
    /// Largest single storage-buffer binding, in bytes.
    pub max_storage_buffer_binding_size: u64,
    /// Largest buffer allocation, in bytes.
    pub max_buffer_size: u64,
    /// Whether `wgpu::Features::SHADER_F64` is available.
    pub shader_f64: bool,
    /// `true` for software rasterisers (`DeviceType::Cpu`: lavapipe,
    /// SwiftShader, WARP).
    pub software: bool,
    /// `true` when the selection policy accepts this adapter.
    pub selectable: bool,
}

/// Result of enumerating GPU adapters (§2.5 of the program plan): every
/// adapter seen, the one the policy would pick, and a human-readable message
/// suitable for `GpuUnavailable` errors and the FFI probe.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct GpuProbe {
    /// `false` when the crate was compiled without the `gpu` feature.
    pub built_with_gpu_feature: bool,
    /// Adapters in the order the selection policy ranks them.
    pub adapters: Vec<GpuAdapterReport>,
    /// Index into `adapters` of the selected adapter, if any is usable.
    pub chosen: Option<usize>,
    /// Summary of the outcome (chosen adapter, or why none was usable).
    pub message: String,
}

impl GpuProbe {
    /// `true` when a usable adapter was found.
    pub fn is_available(&self) -> bool {
        self.chosen.is_some()
    }

    /// The selected adapter's report.
    pub fn chosen_adapter(&self) -> Option<&GpuAdapterReport> {
        self.chosen.and_then(|i| self.adapters.get(i))
    }
}

impl std::fmt::Display for GpuProbe {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.message)
    }
}

/// Enumerate GPU adapters and report which one the solver would use.
///
/// Honours `THESEUS_GPU_ALLOW_SOFTWARE=1` (accept `DeviceType::Cpu`
/// adapters), `WGPU_BACKEND` and `WGPU_ADAPTER_NAME`. Never fails: without
/// the `gpu` feature, or when no adapter is usable, the report says so in
/// `message` and `chosen` is `None`.
pub fn probe_gpu() -> GpuProbe {
    #[cfg(feature = "gpu")]
    {
        gpu::probe(&gpu::AdapterPolicy::from_env())
    }
    #[cfg(not(feature = "gpu"))]
    {
        GpuProbe {
            built_with_gpu_feature: false,
            adapters: Vec::new(),
            chosen: None,
            message: "theseus was not built with the `gpu` feature".to_string(),
        }
    }
}

/// Marker for the compute backend an iterative solver runs on.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum BackendHandle {
    /// Host memory, rayon kernels.
    Cpu,
    /// wgpu device (DX12 / Vulkan / Metal).
    Gpu,
}

impl std::fmt::Display for BackendHandle {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::Cpu => "cpu",
            Self::Gpu => "gpu",
        })
    }
}

/// The kernel set an iterative solver needs, over backend-owned buffers.
///
/// Vector buffers (`Buf`) hold `n * 3` scalars in row-major block layout.
/// Level operators and aggregate maps are uploaded once per topology
/// (`upload_level`, `upload_aggregates`); per-`q` weight changes go through
/// `update_level_weights` so no allocation happens inside a solve.
pub trait Backend {
    /// A vector of `n * 3` scalars (host `Vec` or device buffer) in the
    /// precision it was allocated with.
    type Buf;
    /// Backend-resident form of a [`LevelGraph`]: adjacency, per-edge weights
    /// and per-node anchor weights (plus any derived data such as the Jacobi
    /// diagonal).
    type LevelGraphBuf;
    /// Backend-resident fine-node → coarse-node map of one level.
    type AggBuf;

    /// Which backend this is.
    fn handle(&self) -> BackendHandle;

    // ── Buffers ────────────────────────────────────────────

    /// Allocate a zeroed vector buffer of `len` scalars in `precision`.
    fn alloc(&self, len: usize, precision: Precision) -> Self::Buf;
    /// Copy `src` (host `f64`) into `dst`, converting to the buffer precision.
    fn upload(&self, src: &[f64], dst: &mut Self::Buf);
    /// Copy `src` into `dst` (host `f64`), converting from the buffer precision.
    fn download(&self, src: &Self::Buf, dst: &mut [f64]);
    /// Number of scalars in `buf`.
    fn len(&self, buf: &Self::Buf) -> usize;
    /// `dst = src` (same length and precision).
    fn copy(&self, src: &Self::Buf, dst: &mut Self::Buf);
    /// `buf = 0`.
    fn zero(&self, buf: &mut Self::Buf);

    // ── Level operators ────────────────────────────────────

    /// Upload one level's graph (adjacency, weights, anchors) in `precision`.
    fn upload_level(&self, level: &LevelGraph, precision: Precision) -> Self::LevelGraphBuf;
    /// Refresh the per-edge `weight` and per-node `anchor` arrays of an
    /// uploaded level after `q` changed (lengths as at upload time).
    fn update_level_weights(&self, level: &mut Self::LevelGraphBuf, weight: &[f64], anchor: &[f64]);
    /// Upload the fine → coarse map of a level (`aggregate_of[i] < n_coarse`).
    fn upload_aggregates(&self, aggregate_of: &[u32], n_coarse: usize) -> Self::AggBuf;

    // ── Kernels ────────────────────────────────────────────

    /// `y = A_l x` in gather form:
    /// `(A x)_u = anchor_u x_u + Σ_{e ∋ u} w_e (x_u − x_other)`.
    fn apply_graph(&self, level: &Self::LevelGraphBuf, x: &Self::Buf, y: &mut Self::Buf);
    /// `r = b − A_l x`.
    fn residual(
        &self,
        level: &Self::LevelGraphBuf,
        x: &Self::Buf,
        b: &Self::Buf,
        r: &mut Self::Buf,
    );
    /// One Chebyshev / Jacobi smoothing update on all three columns:
    /// `d = alpha · D_l⁻¹ r + beta · d;  x += d`, where `D_l` is the level's
    /// diagonal (`anchor_u + Σ_{e ∋ u} w_e`) held in `level`.
    fn chebyshev_step(
        &self,
        level: &Self::LevelGraphBuf,
        alpha: f64,
        beta: f64,
        r: &Self::Buf,
        d: &mut Self::Buf,
        x: &mut Self::Buf,
    );
    /// `coarse = 0; coarse[agg[i]] += fine[i]` (piecewise-constant restriction).
    fn restrict(&self, agg: &Self::AggBuf, fine: &Self::Buf, coarse: &mut Self::Buf);
    /// `fine[i] += coarse[agg[i]]` (piecewise-constant prolongation, accumulated).
    fn prolong_add(&self, agg: &Self::AggBuf, coarse: &Self::Buf, fine: &mut Self::Buf);
    /// `y[:, k] += alpha[k] · x[:, k]` for each column `k`.
    fn axpy(&self, alpha: [f64; 3], x: &Self::Buf, y: &mut Self::Buf);
    /// `x[:, k] *= alpha[k]` for each column `k`.
    fn scale(&self, alpha: [f64; 3], x: &mut Self::Buf);
    /// Per-column dot products `Σ_i a[i, k] · b[i, k]`, accumulated in `f64`
    /// in a fixed, thread-count-independent order.
    fn dot3(&self, a: &Self::Buf, b: &Self::Buf) -> [f64; 3];
    /// Per-column Euclidean norms, accumulated in `f64` in a fixed order.
    fn norm3(&self, a: &Self::Buf) -> [f64; 3];
    /// Wait for all queued work (no-op on the CPU).
    fn sync(&self);
}

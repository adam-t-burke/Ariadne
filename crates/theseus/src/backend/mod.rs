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

pub mod cpu;

use crate::graph::LevelGraph;
use crate::linear_solver::Precision;

pub use cpu::CpuBackend;

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

    /// A probe that found nothing usable, with `reason` as the message.
    pub fn unavailable(reason: impl Into<String>) -> Self {
        Self {
            built_with_gpu_feature: cfg!(feature = "gpu"),
            adapters: Vec::new(),
            chosen: None,
            message: reason.into(),
        }
    }

    /// Serialise for the FFI (`theseus_gpu_probe`). Fixed shape consumed by
    /// the C# side:
    ///
    /// ```json
    /// {"available": bool, "adapters": [{"name", "backend", "device_type",
    ///   "shader_f64", "max_storage_buffer_binding_size", "driver",
    ///   "max_buffer_size", "software", "selectable"}],
    ///  "chosen": index|null, "reason": string|null,
    ///  "built_with_gpu_feature": bool}
    /// ```
    ///
    /// `reason` carries `message` when no adapter is usable and is `null`
    /// otherwise; readers must ignore keys they do not know.
    pub fn to_json(&self) -> String {
        let mut out = String::with_capacity(160 + 240 * self.adapters.len());
        out.push_str("{\"available\":");
        out.push_str(if self.is_available() { "true" } else { "false" });
        out.push_str(",\"adapters\":[");
        for (i, a) in self.adapters.iter().enumerate() {
            if i > 0 {
                out.push(',');
            }
            out.push_str("{\"name\":");
            push_json_string(&mut out, &a.name);
            out.push_str(",\"backend\":");
            push_json_string(&mut out, &a.backend);
            out.push_str(",\"device_type\":");
            push_json_string(&mut out, &a.device_type);
            out.push_str(",\"shader_f64\":");
            out.push_str(if a.shader_f64 { "true" } else { "false" });
            out.push_str(",\"max_storage_buffer_binding_size\":");
            out.push_str(&a.max_storage_buffer_binding_size.to_string());
            out.push_str(",\"driver\":");
            push_json_string(&mut out, &a.driver);
            out.push_str(",\"max_buffer_size\":");
            out.push_str(&a.max_buffer_size.to_string());
            out.push_str(",\"software\":");
            out.push_str(if a.software { "true" } else { "false" });
            out.push_str(",\"selectable\":");
            out.push_str(if a.selectable { "true" } else { "false" });
            out.push('}');
        }
        out.push_str("],\"chosen\":");
        match self.chosen {
            Some(i) => out.push_str(&i.to_string()),
            None => out.push_str("null"),
        }
        out.push_str(",\"reason\":");
        if self.is_available() {
            out.push_str("null");
        } else {
            push_json_string(&mut out, &self.message);
        }
        out.push_str(",\"built_with_gpu_feature\":");
        out.push_str(if self.built_with_gpu_feature {
            "true"
        } else {
            "false"
        });
        out.push('}');
        out
    }
}

/// Append `s` as a JSON string literal (RFC 8259 escaping).
fn push_json_string(out: &mut String, s: &str) {
    out.push('"');
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if (c as u32) < 0x20 => {
                out.push_str(&format!("\\u{:04x}", c as u32));
            }
            c => out.push(c),
        }
    }
    out.push('"');
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn unavailable_probe_serialises_with_reason() {
        let probe = GpuProbe::unavailable("no adapter");
        assert!(!probe.is_available());
        assert_eq!(probe.chosen_adapter(), None);
        assert_eq!(probe.built_with_gpu_feature, cfg!(feature = "gpu"));
        let json = probe.to_json();
        let parsed: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(parsed["available"], false);
        assert_eq!(parsed["adapters"].as_array().unwrap().len(), 0);
        assert!(parsed["chosen"].is_null());
        assert_eq!(parsed["reason"], "no adapter");
    }

    #[cfg(not(feature = "gpu"))]
    #[test]
    fn probe_without_feature_reports_not_built() {
        let probe = probe_gpu();
        assert!(!probe.built_with_gpu_feature);
        assert!(!probe.is_available());
        assert!(probe.message.contains("gpu"));
        let parsed: serde_json::Value = serde_json::from_str(&probe.to_json()).unwrap();
        assert_eq!(parsed["available"], false);
        assert_eq!(parsed["built_with_gpu_feature"], false);
    }

    #[test]
    fn json_serialises_adapters_and_escapes_strings() {
        let probe = GpuProbe {
            built_with_gpu_feature: true,
            adapters: vec![
                GpuAdapterReport {
                    name: "Fake \"GPU\"\n".into(),
                    backend: "vulkan".into(),
                    device_type: "discrete".into(),
                    driver: "test 1.0".into(),
                    max_storage_buffer_binding_size: 2_147_483_648,
                    max_buffer_size: 4_294_967_296,
                    shader_f64: true,
                    software: false,
                    selectable: true,
                },
                GpuAdapterReport {
                    name: "llvmpipe".into(),
                    backend: "vulkan".into(),
                    device_type: "cpu".into(),
                    driver: "llvmpipe".into(),
                    max_storage_buffer_binding_size: 134_217_728,
                    max_buffer_size: 2_147_483_648,
                    shader_f64: false,
                    software: true,
                    selectable: false,
                },
            ],
            chosen: Some(0),
            message: "using Fake GPU".into(),
        };
        assert_eq!(
            probe.chosen_adapter().map(|a| a.name.as_str()),
            Some("Fake \"GPU\"\n")
        );
        let json = probe.to_json();
        assert!(json.starts_with(
            r#"{"available":true,"adapters":[{"name":"Fake \"GPU\"\n","backend":"vulkan","device_type":"discrete","shader_f64":true,"max_storage_buffer_binding_size":2147483648,"#
        ));
        let parsed: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(parsed["available"], true);
        assert_eq!(parsed["adapters"].as_array().unwrap().len(), 2);
        assert_eq!(parsed["adapters"][0]["name"], "Fake \"GPU\"\n");
        assert_eq!(parsed["adapters"][1]["software"], true);
        assert_eq!(parsed["adapters"][1]["selectable"], false);
        assert_eq!(parsed["chosen"], 0);
        assert!(parsed["reason"].is_null());
        assert_eq!(parsed["built_with_gpu_feature"], true);
    }
}

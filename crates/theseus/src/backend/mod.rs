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

// ─────────────────────────────────────────────────────────────
//  GPU adapter probe
// ─────────────────────────────────────────────────────────────

/// One GPU adapter seen by the probe.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GpuAdapterInfo {
    /// Driver-reported adapter name.
    pub name: String,
    /// Graphics API the adapter is reached through (`"vulkan"`, `"dx12"`,
    /// `"metal"`, ...).
    pub backend: String,
    /// `"discrete"`, `"integrated"`, `"cpu"` (software) or `"other"`.
    pub device_type: String,
    /// Whether the adapter exposes 64-bit floats in shaders (`SHADER_F64`).
    pub shader_f64: bool,
    /// `max_storage_buffer_binding_size` limit, in bytes.
    pub max_storage_buffer_binding_size: u64,
}

/// Result of enumerating GPU adapters for the `IterativeGpu` backend.
///
/// Serialised to JSON by [`GpuProbe::to_json`] for the FFI
/// (`theseus_gpu_probe`). The document has a fixed shape:
///
/// ```json
/// {"available": false, "adapters": [], "chosen": null, "reason": "..."}
/// ```
///
/// `adapters` lists every adapter seen (fields of [`GpuAdapterInfo`] in
/// snake_case), `chosen` is the index into `adapters` of the one the solver
/// would use (or `null`), and `reason` explains an unavailable result (or is
/// `null`).
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct GpuProbe {
    /// A usable adapter exists.
    pub available: bool,
    /// All adapters enumerated, in enumeration order.
    pub adapters: Vec<GpuAdapterInfo>,
    /// Index into `adapters` of the adapter the solver would use.
    pub chosen: Option<usize>,
    /// Why no adapter is usable, when `available` is `false`.
    pub reason: Option<String>,
}

impl GpuProbe {
    /// The probe result of a build without a GPU backend.
    pub fn unavailable(reason: impl Into<String>) -> Self {
        Self {
            available: false,
            adapters: Vec::new(),
            chosen: None,
            reason: Some(reason.into()),
        }
    }

    /// The adapter the solver would use, if any.
    pub fn chosen_adapter(&self) -> Option<&GpuAdapterInfo> {
        self.chosen.and_then(|i| self.adapters.get(i))
    }

    /// Serialise to the JSON document described on [`GpuProbe`].
    pub fn to_json(&self) -> String {
        let mut out = String::with_capacity(128 + 160 * self.adapters.len());
        out.push_str("{\"available\":");
        out.push_str(if self.available { "true" } else { "false" });
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
            out.push('}');
        }
        out.push_str("],\"chosen\":");
        match self.chosen {
            Some(i) => out.push_str(&i.to_string()),
            None => out.push_str("null"),
        }
        out.push_str(",\"reason\":");
        match &self.reason {
            Some(r) => push_json_string(&mut out, r),
            None => out.push_str("null"),
        }
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

/// Enumerate GPU adapters and pick the one `IterativeGpu` would run on.
///
/// This build has no GPU backend (the wgpu layer is added by the GPU
/// workstream), so the probe always reports `available: false` with a
/// reason. The GPU workstream replaces the body with the real wgpu
/// enumeration (prefer discrete, then integrated; reject software adapters
/// unless `THESEUS_GPU_ALLOW_SOFTWARE=1`); the signature and the JSON shape
/// stay as they are.
pub fn probe_gpu() -> GpuProbe {
    GpuProbe::unavailable("GPU backend not yet built")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn probe_without_backend_is_unavailable_with_reason() {
        let probe = probe_gpu();
        assert!(!probe.available);
        assert!(probe.adapters.is_empty());
        assert_eq!(probe.chosen, None);
        assert_eq!(probe.chosen_adapter(), None);
        assert_eq!(probe.reason.as_deref(), Some("GPU backend not yet built"));
        assert_eq!(
            probe.to_json(),
            r#"{"available":false,"adapters":[],"chosen":null,"reason":"GPU backend not yet built"}"#
        );
    }

    #[test]
    fn json_serialises_adapters_and_escapes_strings() {
        let probe = GpuProbe {
            available: true,
            adapters: vec![
                GpuAdapterInfo {
                    name: "Fake \"GPU\"\n".into(),
                    backend: "vulkan".into(),
                    device_type: "discrete".into(),
                    shader_f64: true,
                    max_storage_buffer_binding_size: 2_147_483_648,
                },
                GpuAdapterInfo {
                    name: "llvmpipe".into(),
                    backend: "vulkan".into(),
                    device_type: "cpu".into(),
                    shader_f64: false,
                    max_storage_buffer_binding_size: 134_217_728,
                },
            ],
            chosen: Some(0),
            reason: None,
        };
        assert_eq!(
            probe.chosen_adapter().map(|a| a.name.as_str()),
            Some("Fake \"GPU\"\n")
        );
        let json = probe.to_json();
        assert_eq!(
            json,
            concat!(
                r#"{"available":true,"adapters":["#,
                r#"{"name":"Fake \"GPU\"\n","backend":"vulkan","device_type":"discrete","shader_f64":true,"max_storage_buffer_binding_size":2147483648},"#,
                r#"{"name":"llvmpipe","backend":"vulkan","device_type":"cpu","shader_f64":false,"max_storage_buffer_binding_size":134217728}"#,
                r#"],"chosen":0,"reason":null}"#
            )
        );
        // The dev-dependency parser confirms the hand-written document is valid JSON.
        let parsed: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(parsed["available"], true);
        assert_eq!(parsed["adapters"].as_array().unwrap().len(), 2);
        assert_eq!(parsed["adapters"][0]["name"], "Fake \"GPU\"\n");
        assert_eq!(parsed["chosen"], 0);
        assert!(parsed["reason"].is_null());
    }
}

# `backend::gpu` — wgpu compute backend (WS-G)

The wgpu implementation of the `Backend` kernel set from
`ITERATIVE_SOLVER_PROGRAM.md` §2.3, designed per §2.5. Built only with the
`gpu` cargo feature (`cargo build -p theseus --release --features gpu`);
the default build and the C# native artefacts do not link wgpu.

| file | contents |
|---|---|
| `adapter.rs` | adapter enumeration and ranking (`AdapterPolicy`), env overrides, `GpuProbe` report, `GpuContext` (device + queue + limits + `shader_f64`) |
| `buffers.rs` | `BufferPool` (storage buffers, OOM error scopes, byte accounting), `GpuBuf` / `IndexBuf`, staging ring for readback, uniform parameter ring |
| `kernels.wgsl` | every kernel, written once against a `Real` alias and a `WG` constant |
| `pipelines.rs` | source specialisation (`f32` / `f64`, workgroup size), bind-group layouts per kernel family, pipelines |
| `backend.rs` | `GpuBackend: Backend`, command batching, host-side reduction finish, `coarse_weight_update` |

`crate::backend::probe_gpu()` is available in every build; without the
feature it reports "not built with the `gpu` feature".

## Dependency pins

`wgpu = "=30.0.1"` (MSRV 1.87, matches the program's 1.87 / toolchain 1.89).
Its transitive `ordered-float` is held at 5.4.0 in `Cargo.lock` because
5.5.0 requires rustc 1.90; `cargo update` on a 1.89 toolchain will refuse
the bump, on newer toolchains it is harmless. wgpu 27 (MSRV 1.88) and 28
(MSRV 1.92) were skipped for the same reason; 29 and 30 declare 1.87.
`pollster = "=1.0.1"` blocks on wgpu futures; `bytemuck` casts `u32`/`f64`
slices. wgpu features enabled: `std`, `parking_lot`, `wgsl`, `vulkan`,
`dx12`, `metal` (no GL, no WebGPU).

## Adapter selection

1. `WGPU_BACKEND` (wgpu's own variable) restricts the backends enumerated;
   `WGPU_ADAPTER_NAME` is a case-insensitive substring filter, with the same
   semantics as `wgpu::util::initialize_adapter_from_env`.
2. Adapters are ranked by `linear_solver::AdapterPreference`: `Discrete`
   (discrete → integrated → virtual → other → cpu), `Integrated`, or `Any`
   (enumeration order). Ties keep enumeration order, so the choice is
   reproducible.
3. `DeviceType::Cpu` adapters (lavapipe, SwiftShader, WARP) are rejected
   unless `THESEUS_GPU_ALLOW_SOFTWARE=1`; the `GpuUnavailable` message says
   so and lists every adapter seen.
4. The device is requested with the **adapter's own limits**, not wgpu's
   defaults: the default `max_storage_buffer_binding_size` is 128 MiB, below
   a 10M-node `f32` vector (120 MB) plus headroom, and the default
   `max_buffer_size` is 256 MiB. `Features::SHADER_F64` is requested when
   the adapter offers it and recorded in `GpuContext::shader_f64` /
   `GpuAdapterReport::shader_f64`.

## Buffers and batching

* Vectors are flat `array<Real>` storage buffers with `v[node * 3 + k]`; no
  `vec3` padding. A runtime-sized array binding must hold ≥ 1 element, so a
  zero-length buffer is still allocated with one element.
* `upload_level` uploads `offsets`, `edges`, `other` (`u32`), `weight`,
  `anchor` (`Real`) and computes `inv_diag` on the device; `sign` is not
  needed by any kernel. `update_level_weights` rewrites `weight`/`anchor`
  in place and recomputes `inv_diag`.
* `upload_aggregates` uploads `aggregate_of` and the inverse CSR (coarse →
  ascending fine members, counting sort on the host) so `restrict` is a
  deterministic gather with no atomics.
* Every kernel records into a pending `CommandEncoder`; `dot3`, `norm3`,
  `download` and `sync` submit it. `upload` and `update_level_weights` also
  submit first, because `queue.write_buffer` is applied at the start of the
  *next* submit and would otherwise overtake pending reads. The ring of
  2048 uniform slots (one per dispatch, dynamic offsets) is reset at each
  submit; `GpuBackend::pending_dispatches()` and `stats()` expose the
  batching for tests and benchmarks.
* Reductions: each thread accumulates a grid-strided subset of nodes, the
  workgroup reduces in shared memory with a fixed tree, thread 0 writes 3
  partials; at most 1024 workgroups; the host sums the partials in `f64` in
  workgroup order. Same input, same bits.
* Readback: two mappable staging buffers (double buffering), grown only when
  a larger readback than before is requested; `BufferPool::reserve_staging`
  preallocates for a known vector size. Each readback is a submit + wait.
* Not yet pooled (WS-H): bind groups are created per dispatch (host-side
  objects, no device allocation), and `queue.write_buffer` stages the
  48/96-byte parameter block internally.

## Kernels

| entry point | threads | notes |
|---|---|---|
| `apply_graph` | per node | anchor term, then incident edges in CSR order — same summation order as `LevelGraph::apply` |
| `residual` | per node | `r = b − A x`, fused |
| `inv_diag` | per node | `1 / (anchor + Σ w)`, `0` for an empty diagonal |
| `chebyshev_step` | per scalar | `d = α D⁻¹ r + β d; x += d` |
| `restrict_sum` | per coarse node | `restrict` is a WGSL reserved word; `cols` = 1 or 3 |
| `prolong_add` | per fine node | |
| `axpy`, `scale` | per scalar | per-column scalars from the uniform `vec4<Real>` |
| `dot_partial`, `norm_partial` | per node + workgroup reduction | |
| `coarse_weight_update` | per coarse edge | plus `restrict_sum(cols = 1)` for anchors and `inv_diag`, in `GpuBackend::coarse_weight_update` |
| `apply_csr` | per row | `y = M x` (`flag = 0`) or `y += M x` (`flag = 1`) over a `LevelMatrix` (CSR, ascending columns), entries in stored order — same order as the CPU `csr_row`; used for the coarse AMG operators and the smoothed `P` / `Pᵀ` (WS-C) |
| `residual_csr` | per row | `r = b − M x` for a square CSR level |

`chebyshev_step_csr` reuses the `chebyshev_step` entry point with the CSR
level's inverse diagonal bound at binding 9; that diagonal is computed on
the host (`1 / diag` in `f64`, rounded to the buffer precision, like
`CpuCsrData`) at `upload_csr` / `update_csr_values`.

All kernels use a grid-stride loop, so sizes above
`max_compute_workgroups_per_dimension × WG` (65535 × 128 ≈ 8.4M) still
dispatch in one call. Binding numbers are unique across the module (0 =
uniform params, 1–26 storage), so each kernel family's bind-group layout is
a strict subset of group 0 and every entry point compiles from one module.

## `f64` status (WGSL / naga)

* naga's WGSL front end accepts the `f64` type keyword without any `enable`
  directive (there is no `enable f64;` extension in WGSL or naga 30); wgpu
  validates that the device has `Features::SHADER_F64`. The `f64` kernel set
  is therefore generated by replacing `alias Real = f32;` with
  `alias Real = f64;` and compiled only when `shader_f64` is true.
* Verified under lavapipe (Mesa 25.2.8, `shaderFloat64 = true`): all
  kernels pass at `1e-12` relative in `f64` — `vec3<f64>`, `vec4<f64>` in a
  uniform struct (32-byte alignment, see `Params::encode`), `f64` workgroup
  arrays and `Real(…)` constructors through the alias all lower correctly
  to SPIR-V.
* wgpu exposes `SHADER_F64` only on Vulkan (`wgpu-hal/src/vulkan`); DX12 and
  Metal never report it, so on Windows-DX12 and macOS the outer loop runs on
  the host (§2.5). `Precision::F64` buffers on such a device fail with
  `IterativeSolverUnsupported` (`try_alloc`) or panic (`Backend::alloc`).
* `f64` literals: only `Real(0)` / `Real(1)` constructors and abstract-float
  expressions are used in the shader, so no `f32`-suffixed literal is ever
  narrowed.

## Adapter limits encountered

| adapter | backend | type | `max_storage_buffer_binding_size` | `max_buffer_size` | `SHADER_F64` | notes |
|---|---|---|---|---|---|---|
| llvmpipe (LLVM 20.1.2, 256 bits), Mesa 25.2.8 | vulkan | cpu | 128 MiB | 2047 MiB (`0xFFFFFFFF`) | yes | `maxComputeWorkGroupInvocations` 1024, shared memory 32 KiB, uniform range 64 KiB; the 10M-node grid's `edges`/`other` arrays (2·ne·4 B = 160 MB) exceed one binding, so `upload_level` reports `IterativeSolverUnsupported` |

Consequences:

* Any array above `max_storage_buffer_binding_size` is rejected up front
  with an `IterativeSolverUnsupported` message naming the array and the
  limit; real allocation failure is `GpuOutOfMemory` (wgpu error scope).
  Splitting large arrays over several bindings (§2.5 "if the adapter limit
  is below 2 GB") is left to WS-H; on discrete Vulkan adapters
  `maxStorageBufferRange` is 4 GiB, on Metal it equals the maximum buffer
  length, on DX12 wgpu reports 4 GiB. Older Intel Vulkan drivers report
  128 MiB — to verify on the §5.4 laptop-class machine.
* At 10M nodes in `f32` one vector is 120 MB and the CSR `edges`+`other`
  arrays 320 MB; in `f64` a vector is 240 MB. Budget per level at 10M edges
  (5M nodes, f32): ≈ 60 MB per vector, 160 MB CSR, 40 MB weights — device
  memory for a 7-level hierarchy with ~8 vectors on the fine level is
  ≈ 1.2 GB, inside the `AMG_PLAN.md` envelope.

Windows (WARP through DX12) and macOS (Metal) limits are not yet recorded:
the CI matrix (`gpu-software` job) prints them in the
`probe_reports_adapters` test with `--nocapture`, and WS-H's validation
machines fill in the real-GPU rows.

## Running

```sh
# Linux, software adapter
sudo apt-get install -y mesa-vulkan-drivers libvulkan1 vulkan-tools
export VK_ICD_FILENAMES=$(ls /usr/share/vulkan/icd.d/lvp_icd*.json | head -n1)
export THESEUS_GPU_ALLOW_SOFTWARE=1
cargo test -p theseus --release --features gpu --test backend_gpu_equivalence -- --nocapture
cargo run  -p theseus --release --features gpu --example gpu_kernels_bench

# Windows, WARP:  WGPU_BACKEND=dx12  WGPU_ADAPTER_NAME="Microsoft Basic Render Driver"
# macOS:          WGPU_BACKEND=metal  (the runner's GPU)
```

Benchmark knobs: `THESEUS_GPU_BENCH_NODES=1000000,10000000`,
`THESEUS_GPU_BENCH_REPS=10`, `THESEUS_GPU_BENCH_WG=64,128,256`,
`THESEUS_GPU_BENCH_PRECISION=f32,f64`.

## Deviations and open questions for WS-H / the integrator

* `Backend::alloc`, `upload_level`, `upload_aggregates` return values, not
  `Result`; the GPU versions panic on allocation failure with the
  `GpuOutOfMemory` / `IterativeSolverUnsupported` message. `try_alloc`,
  `try_upload_level`, `try_upload_aggregates` are the fallible forms the
  solver should use to size its working set. A trait-level `Result` would
  remove the panic path (§0.2 "loud failure").
* `coarse_weight_update` (per coarse edge over its fine-edge list) has no
  slot in the `Backend` trait, whose `update_level_weights` takes host
  slices; it is a `GpuBackend` method taking a `GpuCoarseEdgeMap` uploaded
  from a host CSR. WS-C's `hierarchy.rs` decides the coarse edge order; the
  test builds the quotient graph with coarse edges in ascending
  `(min(U,V), max(U,V))` order.
* No trait method reads back a level's arrays; `GpuLevel::{weight, anchor,
  inv_diag}` expose the buffers for `download`.
* Bind-group caching and a mapped parameter ring (instead of
  `queue.write_buffer`) are the two obvious per-dispatch host costs to
  remove once WS-H measures a full K-cycle.

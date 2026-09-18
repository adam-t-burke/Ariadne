//! wgpu compute backend (program plan §2.5): adapter probing and selection,
//! device buffers, WGSL kernels and the [`GpuBackend`] implementation of
//! [`Backend`](crate::backend::Backend).
//!
//! Built only with the `gpu` cargo feature. Compute-only (no surface); the
//! WGSL kernels run in `f32` on every adapter and in `f64` where the adapter
//! offers `Features::SHADER_F64` (Vulkan on discrete GPUs and lavapipe;
//! never Metal). See `README.md` in this directory for adapter limits and
//! the WGSL/naga `f64` status.
//!
//! Typical use:
//!
//! ```no_run
//! use theseus::backend::gpu::{AdapterPolicy, GpuBackend, GpuContext};
//! use theseus::backend::Backend;
//! use theseus::linear_solver::Precision;
//!
//! let ctx = GpuContext::new(&AdapterPolicy::from_env())?;
//! let gpu = GpuBackend::new(ctx)?;
//! let mut x = gpu.alloc(3 * 1000, Precision::F32);
//! gpu.upload(&vec![1.0; 3 * 1000], &mut x);
//! assert_eq!(gpu.norm3(&x)[0], (1000f64).sqrt());
//! # Ok::<(), theseus::TheseusError>(())
//! ```

mod adapter;
mod backend;
mod buffers;
mod pipelines;

pub use adapter::{enumerate, instance, probe, AdapterPolicy, GpuContext, ALLOW_SOFTWARE_ENV};
pub use backend::{members_csr, GpuAggregates, GpuBackend, GpuCoarseEdgeMap, GpuLevel, GpuStats};
pub use buffers::{BufferPool, GpuBuf, IndexBuf};
pub use pipelines::{shader_source, KernelSet, KERNELS_WGSL, WORKGROUP_SIZES};

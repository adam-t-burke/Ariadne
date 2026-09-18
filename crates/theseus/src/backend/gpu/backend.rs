//! [`GpuBackend`]: the [`Backend`] kernel set on a wgpu device.
//!
//! Every kernel call records one dispatch into a pending command encoder
//! (the "encoder scope"); nothing is submitted until a result must leave the
//! device (`dot3`, `norm3`, `download`, `sync`), a host write must be ordered
//! after pending reads (`upload`, `update_level_weights`), or the per-dispatch
//! parameter ring fills up. Reductions finish on the host in `f64` in fixed
//! workgroup order, so results are deterministic for a given size and
//! workgroup size.

use std::cell::RefCell;
use std::sync::mpsc;

use super::adapter::GpuContext;
use super::buffers::{BufferPool, GpuBuf, IndexBuf, Params};
use super::pipelines::{Kernel, KernelSet};
use crate::backend::{Backend, BackendHandle, GpuProbe};
use crate::graph::LevelGraph;
use crate::linear_solver::Precision;
use crate::types::TheseusError;

/// Maximum workgroups a reduction dispatches; the host sums this many
/// partials per column.
const REDUCE_MAX_GROUPS: u32 = 1024;

/// Dispatches recorded before the pending encoder is submitted anyway.
const MAX_BATCHED_DISPATCHES: u32 = 4096;

/// One uploaded [`LevelGraph`]: CSR adjacency (`sign` is not needed by the
/// kernels), per-edge weights, per-node anchors and the Jacobi inverse
/// diagonal computed on the device.
#[derive(Debug)]
pub struct GpuLevel {
    n: usize,
    num_edges: usize,
    precision: Precision,
    offsets: IndexBuf,
    edges: IndexBuf,
    other: IndexBuf,
    weight: GpuBuf,
    anchor: GpuBuf,
    inv_diag: GpuBuf,
}

impl GpuLevel {
    /// Nodes on the level.
    pub fn n(&self) -> usize {
        self.n
    }

    /// Edges on the level.
    pub fn num_edges(&self) -> usize {
        self.num_edges
    }

    /// Storage precision of weights, anchors and the inverse diagonal.
    pub fn precision(&self) -> Precision {
        self.precision
    }

    /// Per-edge weights (readable through `GpuBackend::download`).
    pub fn weight(&self) -> &GpuBuf {
        &self.weight
    }

    /// Per-node anchor weights.
    pub fn anchor(&self) -> &GpuBuf {
        &self.anchor
    }

    /// Per-node `1 / (anchor_u + Σ w_e)` (0 where the diagonal is 0).
    pub fn inv_diag(&self) -> &GpuBuf {
        &self.inv_diag
    }
}

/// One uploaded aggregate map: fine → coarse for `prolong_add`, and the
/// inverse CSR (coarse → sorted fine members) so `restrict` is a gather.
#[derive(Debug)]
pub struct GpuAggregates {
    n_fine: usize,
    n_coarse: usize,
    agg: IndexBuf,
    coarse_offsets: IndexBuf,
    members: IndexBuf,
}

impl GpuAggregates {
    /// Fine nodes.
    pub fn n_fine(&self) -> usize {
        self.n_fine
    }

    /// Coarse nodes.
    pub fn n_coarse(&self) -> usize {
        self.n_coarse
    }
}

/// Coarse edge → fine edges CSR for [`GpuBackend::coarse_weight_update`].
#[derive(Debug)]
pub struct GpuCoarseEdgeMap {
    n_coarse_edges: usize,
    offsets: IndexBuf,
    fine_edges: IndexBuf,
}

impl GpuCoarseEdgeMap {
    /// Coarse edges.
    pub fn n_coarse_edges(&self) -> usize {
        self.n_coarse_edges
    }
}

/// Counters for benchmarks and batching tests.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct GpuStats {
    /// Compute dispatches recorded.
    pub dispatches: u64,
    /// Queue submissions.
    pub submits: u64,
    /// Host readbacks (map + wait).
    pub readbacks: u64,
    /// Bytes written host → device (vector uploads, not parameters).
    pub bytes_uploaded: u64,
    /// Bytes read device → host.
    pub bytes_downloaded: u64,
}

struct Encoder {
    encoder: wgpu::CommandEncoder,
    dispatches: u32,
}

/// The wgpu implementation of [`Backend`].
///
/// `Precision::F32` is always available; `Precision::F64` only when the
/// adapter reports `SHADER_F64` (`supports`). Buffers, levels and aggregate
/// maps are created through the `try_*` constructors, which report
/// allocation failures as [`TheseusError::GpuOutOfMemory`]; the infallible
/// [`Backend`] methods panic with the same message, so solvers should size
/// their working set through `try_*` first.
pub struct GpuBackend {
    ctx: GpuContext,
    pool: BufferPool,
    f32_kernels: KernelSet,
    f64_kernels: Option<KernelSet>,
    workgroup_size: u32,
    max_groups: u32,
    pending: RefCell<Option<Encoder>>,
    partials_f32: GpuBuf,
    partials_f64: Option<GpuBuf>,
    dummy: wgpu::Buffer,
    partials_host: RefCell<Vec<f64>>,
    stats: RefCell<GpuStats>,
}

impl std::fmt::Debug for GpuBackend {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("GpuBackend")
            .field("adapter", &self.ctx.report.name)
            .field("backend", &self.ctx.report.backend)
            .field("shader_f64", &self.ctx.shader_f64)
            .field("workgroup_size", &self.workgroup_size)
            .field("device_bytes", &self.pool.device_bytes())
            .finish()
    }
}

impl GpuBackend {
    /// Default workgroup size (Stage-D tuning sweeps 64 / 128 / 256).
    pub const DEFAULT_WORKGROUP_SIZE: u32 = 128;

    /// Build the backend on an opened device with the default workgroup size.
    pub fn new(ctx: GpuContext) -> Result<Self, TheseusError> {
        Self::with_workgroup_size(ctx, Self::DEFAULT_WORKGROUP_SIZE)
    }

    /// Build the backend with kernels compiled for `workgroup_size` (a power
    /// of two ≤ `max_compute_workgroup_size_x`).
    pub fn with_workgroup_size(ctx: GpuContext, workgroup_size: u32) -> Result<Self, TheseusError> {
        if !workgroup_size.is_power_of_two()
            || workgroup_size > ctx.limits.max_compute_workgroup_size_x
            || workgroup_size > ctx.limits.max_compute_invocations_per_workgroup
        {
            return Err(TheseusError::GpuUnavailable(format!(
                "workgroup size {workgroup_size} not supported by {} (max {})",
                ctx.report.name,
                ctx.limits
                    .max_compute_workgroup_size_x
                    .min(ctx.limits.max_compute_invocations_per_workgroup)
            )));
        }
        let device = &ctx.device;
        let pool = BufferPool::new(device, &ctx.queue, &ctx.limits);
        let f32_kernels = KernelSet::new(device, Precision::F32, workgroup_size);
        let f64_kernels = ctx
            .shader_f64
            .then(|| KernelSet::new(device, Precision::F64, workgroup_size));
        let partials_len = REDUCE_MAX_GROUPS as usize * 3;
        let partials_f32 = pool.create_vec(partials_len, Precision::F32, "theseus partials f32")?;
        let partials_f64 = if ctx.shader_f64 {
            Some(pool.create_vec(partials_len, Precision::F64, "theseus partials f64")?)
        } else {
            None
        };
        let dummy = pool.create_storage(16, "theseus dummy binding")?;
        let max_groups = ctx.limits.max_compute_workgroups_per_dimension.max(1);
        Ok(Self {
            ctx,
            pool,
            f32_kernels,
            f64_kernels,
            workgroup_size,
            max_groups,
            pending: RefCell::new(None),
            partials_f32,
            partials_f64,
            dummy,
            partials_host: RefCell::new(vec![0.0; partials_len]),
            stats: RefCell::new(GpuStats::default()),
        })
    }

    /// Open the adapter selected by the environment and build the backend.
    pub fn from_env() -> Result<Self, TheseusError> {
        Self::new(GpuContext::from_env()?)
    }

    /// The device context.
    pub fn context(&self) -> &GpuContext {
        &self.ctx
    }

    /// Probe report the device was created from.
    pub fn probe(&self) -> &GpuProbe {
        &self.ctx.probe
    }

    /// `true` when `SHADER_F64` kernels are compiled.
    pub fn shader_f64(&self) -> bool {
        self.f64_kernels.is_some()
    }

    /// Whether buffers and kernels in `precision` are available.
    pub fn supports(&self, precision: Precision) -> bool {
        match precision {
            Precision::F32 => true,
            Precision::F64 => self.shader_f64(),
        }
    }

    /// Workgroup size the kernels were compiled with.
    pub fn workgroup_size(&self) -> u32 {
        self.workgroup_size
    }

    /// Device bytes allocated so far (buffers, levels, aggregates, staging,
    /// parameter ring).
    pub fn device_bytes(&self) -> u64 {
        self.pool.device_bytes()
    }

    /// Buffer pool (for `reserve_staging` and limits).
    pub fn pool(&self) -> &BufferPool {
        &self.pool
    }

    /// Counters since creation.
    pub fn stats(&self) -> GpuStats {
        *self.stats.borrow()
    }

    /// Dispatches recorded but not yet submitted.
    pub fn pending_dispatches(&self) -> u32 {
        self.pending.borrow().as_ref().map_or(0, |e| e.dispatches)
    }

    fn kernels(&self, precision: Precision) -> &KernelSet {
        match precision {
            Precision::F32 => &self.f32_kernels,
            Precision::F64 => self.f64_kernels.as_ref().unwrap_or_else(|| {
                panic!(
                    "gpu: Precision::F64 requested but {} has no SHADER_F64",
                    self.ctx.report.name
                )
            }),
        }
    }

    fn check_precision(&self, precision: Precision) -> Result<(), TheseusError> {
        if self.supports(precision) {
            Ok(())
        } else {
            Err(TheseusError::IterativeSolverUnsupported(format!(
                "Precision::F64 on the GPU needs SHADER_F64, which {} ({}) does not offer",
                self.ctx.report.name, self.ctx.report.backend
            )))
        }
    }

    fn groups_for(&self, n: usize) -> u32 {
        let per = self.workgroup_size as usize;
        let groups = n.div_ceil(per).max(1);
        u32::try_from(groups)
            .unwrap_or(u32::MAX)
            .min(self.max_groups)
    }

    fn reduce_groups(&self, n: usize) -> u32 {
        self.groups_for(n).min(REDUCE_MAX_GROUPS)
    }

    // ── Fallible constructors ────────────────────────────────

    /// Zeroed vector of `len` scalars in `precision`.
    pub fn try_alloc(&self, len: usize, precision: Precision) -> Result<GpuBuf, TheseusError> {
        self.check_precision(precision)?;
        self.pool.create_vec(len, precision, "theseus vector")
    }

    /// Upload one level (adjacency, weights, anchors) and compute its inverse
    /// diagonal on the device.
    pub fn try_upload_level(
        &self,
        level: &LevelGraph,
        precision: Precision,
    ) -> Result<GpuLevel, TheseusError> {
        self.check_precision(precision)?;
        let adj = &level.adjacency;
        assert_eq!(
            adj.num_nodes(),
            level.n,
            "gpu: adjacency/node count mismatch"
        );
        assert_eq!(level.anchor.len(), level.n, "gpu: anchor length mismatch");
        let offsets = self
            .pool
            .create_u32(&adj.offsets, "theseus level offsets")?;
        let edges = self.pool.create_u32(&adj.edges, "theseus level edges")?;
        let other = self.pool.create_u32(&adj.other, "theseus level other")?;
        let weight = self
            .pool
            .create_vec_from(&level.weight, precision, "theseus level weight")?;
        let anchor = self
            .pool
            .create_vec_from(&level.anchor, precision, "theseus level anchor")?;
        let inv_diag = self
            .pool
            .create_vec(level.n, precision, "theseus level inv_diag")?;
        let uploaded = GpuLevel {
            n: level.n,
            num_edges: level.weight.len(),
            precision,
            offsets,
            edges,
            other,
            weight,
            anchor,
            inv_diag,
        };
        self.stats.borrow_mut().bytes_uploaded += uploaded.weight.data_bytes()
            + uploaded.anchor.data_bytes()
            + (adj.offsets.len() + adj.edges.len() + adj.other.len()) as u64 * 4;
        self.compute_inv_diag(&uploaded);
        Ok(uploaded)
    }

    /// Upload a fine → coarse map together with its inverse CSR.
    pub fn try_upload_aggregates(
        &self,
        aggregate_of: &[u32],
        n_coarse: usize,
    ) -> Result<GpuAggregates, TheseusError> {
        let (coarse_offsets, members) = members_csr(aggregate_of, n_coarse);
        Ok(GpuAggregates {
            n_fine: aggregate_of.len(),
            n_coarse,
            agg: self.pool.create_u32(aggregate_of, "theseus aggregate_of")?,
            coarse_offsets: self
                .pool
                .create_u32(&coarse_offsets, "theseus coarse offsets")?,
            members: self.pool.create_u32(&members, "theseus coarse members")?,
        })
    }

    /// Upload the coarse edge → fine edges CSR (`offsets.len() = coarse edges
    /// + 1`, `fine_edges` sorted within each list for determinism).
    pub fn upload_coarse_edge_map(
        &self,
        offsets: &[u32],
        fine_edges: &[u32],
    ) -> Result<GpuCoarseEdgeMap, TheseusError> {
        assert!(
            !offsets.is_empty(),
            "gpu: coarse edge offsets need n + 1 entries"
        );
        assert_eq!(
            *offsets.last().unwrap() as usize,
            fine_edges.len(),
            "gpu: coarse edge CSR is inconsistent"
        );
        Ok(GpuCoarseEdgeMap {
            n_coarse_edges: offsets.len() - 1,
            offsets: self
                .pool
                .create_u32(offsets, "theseus coarse edge offsets")?,
            fine_edges: self
                .pool
                .create_u32(fine_edges, "theseus coarse edge fine edges")?,
        })
    }

    // ── Extra kernels (not in the `Backend` trait) ────────────

    /// Recompute the coarse level's weights, anchors and inverse diagonal on
    /// the device from the fine level: `w_E = Σ w_e` over the fine edges of
    /// each coarse edge, `anchor_U = Σ anchor_u` over the members of each
    /// aggregate.
    pub fn coarse_weight_update(
        &self,
        map: &GpuCoarseEdgeMap,
        agg: &GpuAggregates,
        fine: &GpuLevel,
        coarse: &mut GpuLevel,
    ) {
        assert_eq!(
            fine.precision, coarse.precision,
            "gpu: level precision mismatch"
        );
        assert_eq!(
            map.n_coarse_edges, coarse.num_edges,
            "gpu: coarse edge count mismatch"
        );
        assert_eq!(
            agg.n_fine, fine.n,
            "gpu: aggregate map / fine level mismatch"
        );
        assert_eq!(
            agg.n_coarse, coarse.n,
            "gpu: aggregate map / coarse level mismatch"
        );
        let precision = fine.precision;
        if coarse.num_edges > 0 {
            self.dispatch(
                Kernel::CoarseWeightUpdate,
                precision,
                Params::count(coarse.num_edges),
                &[
                    &map.offsets.buffer,
                    &map.fine_edges.buffer,
                    &fine.weight.buffer,
                    &coarse.weight.buffer,
                ],
                self.groups_for(coarse.num_edges),
            );
        }
        if coarse.n > 0 {
            self.dispatch(
                Kernel::Restrict,
                precision,
                Params::count(coarse.n).cols(1),
                &[
                    &agg.coarse_offsets.buffer,
                    &agg.members.buffer,
                    &agg.agg.buffer,
                    &fine.anchor.buffer,
                    &coarse.anchor.buffer,
                ],
                self.groups_for(coarse.n),
            );
        }
        self.compute_inv_diag(coarse);
    }

    fn compute_inv_diag(&self, level: &GpuLevel) {
        if level.n == 0 {
            return;
        }
        self.dispatch(
            Kernel::InvDiag,
            level.precision,
            Params::count(level.n),
            &[
                &level.offsets.buffer,
                &level.edges.buffer,
                &level.other.buffer,
                &level.weight.buffer,
                &level.anchor.buffer,
                &self.dummy,
                &self.dummy,
                &level.inv_diag.buffer,
            ],
            self.groups_for(level.n),
        );
    }

    // ── Command batching ─────────────────────────────────────

    /// Submit the pending encoder (if any) without waiting.
    pub fn flush(&self) -> Option<wgpu::SubmissionIndex> {
        let pending = self.pending.borrow_mut().take()?;
        let index = self.ctx.queue.submit([pending.encoder.finish()]);
        // Every parameter slot written so far belongs to this submission;
        // later writes are applied at the start of the next one.
        self.pool.params.reset();
        self.stats.borrow_mut().submits += 1;
        Some(index)
    }

    fn wait(&self, index: Option<wgpu::SubmissionIndex>) {
        let poll = match index {
            Some(index) => wgpu::PollType::Wait {
                submission_index: Some(index),
                timeout: None,
            },
            None => wgpu::PollType::wait_indefinitely(),
        };
        self.ctx
            .device
            .poll(poll)
            .expect("gpu: device poll failed (device lost?)");
    }

    fn with_encoder<R>(&self, f: impl FnOnce(&mut wgpu::CommandEncoder) -> R) -> R {
        let mut pending = self.pending.borrow_mut();
        let enc = pending.get_or_insert_with(|| Encoder {
            encoder: self
                .ctx
                .device
                .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                    label: Some("theseus batch"),
                }),
            dispatches: 0,
        });
        f(&mut enc.encoder)
    }

    /// Record one dispatch of `kernel`; `resources` follow the binding order
    /// of the kernel's family.
    fn dispatch(
        &self,
        kernel: Kernel,
        precision: Precision,
        params: Params,
        resources: &[&wgpu::Buffer],
        groups: u32,
    ) {
        if self.pool.params.is_full() || self.pending_dispatches() >= MAX_BATCHED_DISPATCHES {
            self.flush();
        }
        let offset = self.pool.params.push(&self.ctx.queue, &params, precision);
        let set = self.kernels(precision);
        let family = kernel.family();
        debug_assert_eq!(resources.len(), family.bindings().len());
        let mut entries = Vec::with_capacity(resources.len() + 1);
        entries.push(wgpu::BindGroupEntry {
            binding: 0,
            resource: self.pool.params.binding(precision),
        });
        for (&(binding, _), buffer) in family.bindings().iter().zip(resources) {
            entries.push(wgpu::BindGroupEntry {
                binding,
                resource: buffer.as_entire_binding(),
            });
        }
        let bind_group = self
            .ctx
            .device
            .create_bind_group(&wgpu::BindGroupDescriptor {
                label: Some(kernel.entry_point()),
                layout: set.layout(family),
                entries: &entries,
            });
        let mut pending = self.pending.borrow_mut();
        let enc = pending.get_or_insert_with(|| Encoder {
            encoder: self
                .ctx
                .device
                .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                    label: Some("theseus batch"),
                }),
            dispatches: 0,
        });
        {
            let mut pass = enc
                .encoder
                .begin_compute_pass(&wgpu::ComputePassDescriptor {
                    label: Some(kernel.entry_point()),
                    timestamp_writes: None,
                });
            pass.set_pipeline(set.pipeline(kernel));
            pass.set_bind_group(0, &bind_group, &[offset]);
            pass.dispatch_workgroups(groups, 1, 1);
        }
        enc.dispatches += 1;
        self.stats.borrow_mut().dispatches += 1;
    }

    /// Copy `bytes` of `src` into a staging buffer, submit everything pending,
    /// wait, and decode into `dst`.
    fn read_back(&self, src: &wgpu::Buffer, bytes: u64, precision: Precision, dst: &mut [f64]) {
        if bytes == 0 {
            return;
        }
        let staging = self.pool.staging(bytes);
        staging.with_buffer(|staging| {
            self.with_encoder(|enc| enc.copy_buffer_to_buffer(src, 0, staging, 0, bytes));
        });
        let index = self.flush();
        staging.with_buffer(|staging| {
            let slice = staging.slice(..bytes);
            let (tx, rx) = mpsc::channel();
            slice.map_async(wgpu::MapMode::Read, move |r| {
                let _ = tx.send(r);
            });
            self.wait(index);
            rx.recv()
                .expect("gpu: map callback dropped")
                .expect("gpu: buffer map failed");
            {
                let view = slice.get_mapped_range().expect("gpu: mapped range");
                BufferPool::decode(precision, &view, dst);
            }
            staging.unmap();
        });
        let mut stats = self.stats.borrow_mut();
        stats.readbacks += 1;
        stats.bytes_downloaded += bytes;
    }

    fn reduce(&self, kernel: Kernel, a: &GpuBuf, b: &GpuBuf) -> [f64; 3] {
        assert_eq!(a.precision, b.precision, "gpu: precision mismatch");
        assert_eq!(a.len, b.len, "gpu: length mismatch");
        assert_eq!(a.len % 3, 0, "gpu: vectors hold blocks of 3");
        let n = a.len / 3;
        if n == 0 {
            return [0.0; 3];
        }
        let precision = a.precision;
        let partials = match precision {
            Precision::F32 => &self.partials_f32,
            Precision::F64 => self.partials_f64.as_ref().expect("f64 partials"),
        };
        let groups = self.reduce_groups(n);
        self.dispatch(
            kernel,
            precision,
            Params::count(n),
            &[&a.buffer, &b.buffer, &partials.buffer],
            groups,
        );
        let count = groups as usize * 3;
        let mut host = self.partials_host.borrow_mut();
        self.read_back(
            &partials.buffer,
            (count * precision.size_of()) as u64,
            precision,
            &mut host[..count],
        );
        let mut out = [0.0f64; 3];
        for g in 0..groups as usize {
            for (k, o) in out.iter_mut().enumerate() {
                *o += host[3 * g + k];
            }
        }
        out
    }
}

/// Inverse CSR of an aggregate map: `coarse_offsets[U]..coarse_offsets[U+1]`
/// lists the fine members of `U` in ascending order (counting sort).
pub fn members_csr(aggregate_of: &[u32], n_coarse: usize) -> (Vec<u32>, Vec<u32>) {
    let mut counts = vec![0u32; n_coarse + 1];
    for &a in aggregate_of {
        assert!((a as usize) < n_coarse, "gpu: aggregate index out of range");
        counts[a as usize + 1] += 1;
    }
    for u in 0..n_coarse {
        counts[u + 1] += counts[u];
    }
    let mut next = counts.clone();
    let mut members = vec![0u32; aggregate_of.len()];
    for (i, &a) in aggregate_of.iter().enumerate() {
        let slot = &mut next[a as usize];
        members[*slot as usize] = i as u32;
        *slot += 1;
    }
    (counts, members)
}

impl Backend for GpuBackend {
    type Buf = GpuBuf;
    type LevelGraphBuf = GpuLevel;
    type AggBuf = GpuAggregates;

    fn handle(&self) -> BackendHandle {
        BackendHandle::Gpu
    }

    fn alloc(&self, len: usize, precision: Precision) -> GpuBuf {
        self.try_alloc(len, precision)
            .unwrap_or_else(|e| panic!("gpu: alloc of {len} scalars failed: {e}"))
    }

    fn upload(&self, src: &[f64], dst: &mut GpuBuf) {
        // Pending commands may still read `dst`; write_buffer would otherwise
        // be applied before them at the next submit.
        self.flush();
        self.pool.write_vec(dst, src);
        self.stats.borrow_mut().bytes_uploaded += dst.data_bytes();
    }

    fn download(&self, src: &GpuBuf, dst: &mut [f64]) {
        assert_eq!(src.len, dst.len(), "gpu: download length mismatch");
        self.read_back(&src.buffer, src.data_bytes(), src.precision, dst);
    }

    fn len(&self, buf: &GpuBuf) -> usize {
        buf.len
    }

    fn copy(&self, src: &GpuBuf, dst: &mut GpuBuf) {
        assert_eq!(src.len, dst.len, "gpu: copy length mismatch");
        assert_eq!(src.precision, dst.precision, "gpu: copy precision mismatch");
        let bytes = src.data_bytes();
        if bytes == 0 {
            return;
        }
        self.with_encoder(|enc| enc.copy_buffer_to_buffer(&src.buffer, 0, &dst.buffer, 0, bytes));
    }

    fn zero(&self, buf: &mut GpuBuf) {
        let bytes = buf.data_bytes();
        if bytes == 0 {
            return;
        }
        self.with_encoder(|enc| enc.clear_buffer(&buf.buffer, 0, Some(bytes)));
    }

    fn upload_level(&self, level: &LevelGraph, precision: Precision) -> GpuLevel {
        self.try_upload_level(level, precision)
            .unwrap_or_else(|e| panic!("gpu: level upload failed: {e}"))
    }

    fn update_level_weights(&self, level: &mut GpuLevel, weight: &[f64], anchor: &[f64]) {
        assert_eq!(weight.len(), level.num_edges, "gpu: weight length mismatch");
        assert_eq!(anchor.len(), level.n, "gpu: anchor length mismatch");
        self.flush();
        self.pool.write_vec(&level.weight, weight);
        self.pool.write_vec(&level.anchor, anchor);
        self.stats.borrow_mut().bytes_uploaded +=
            level.weight.data_bytes() + level.anchor.data_bytes();
        self.compute_inv_diag(level);
    }

    fn upload_aggregates(&self, aggregate_of: &[u32], n_coarse: usize) -> GpuAggregates {
        self.try_upload_aggregates(aggregate_of, n_coarse)
            .unwrap_or_else(|e| panic!("gpu: aggregate upload failed: {e}"))
    }

    fn apply_graph(&self, level: &GpuLevel, x: &GpuBuf, y: &mut GpuBuf) {
        assert_eq!(x.len, level.n * 3, "gpu: x length mismatch");
        assert_eq!(y.len, level.n * 3, "gpu: y length mismatch");
        assert!(
            x.precision == level.precision && y.precision == level.precision,
            "gpu: precision mismatch"
        );
        if level.n == 0 {
            return;
        }
        self.dispatch(
            Kernel::ApplyGraph,
            level.precision,
            Params::count(level.n),
            &[
                &level.offsets.buffer,
                &level.edges.buffer,
                &level.other.buffer,
                &level.weight.buffer,
                &level.anchor.buffer,
                &x.buffer,
                &self.dummy,
                &y.buffer,
            ],
            self.groups_for(level.n),
        );
    }

    fn residual(&self, level: &GpuLevel, x: &GpuBuf, b: &GpuBuf, r: &mut GpuBuf) {
        assert_eq!(x.len, level.n * 3, "gpu: x length mismatch");
        assert_eq!(b.len, level.n * 3, "gpu: b length mismatch");
        assert_eq!(r.len, level.n * 3, "gpu: r length mismatch");
        assert!(
            x.precision == level.precision
                && b.precision == level.precision
                && r.precision == level.precision,
            "gpu: precision mismatch"
        );
        if level.n == 0 {
            return;
        }
        self.dispatch(
            Kernel::Residual,
            level.precision,
            Params::count(level.n),
            &[
                &level.offsets.buffer,
                &level.edges.buffer,
                &level.other.buffer,
                &level.weight.buffer,
                &level.anchor.buffer,
                &x.buffer,
                &b.buffer,
                &r.buffer,
            ],
            self.groups_for(level.n),
        );
    }

    fn chebyshev_step(
        &self,
        level: &GpuLevel,
        alpha: f64,
        beta: f64,
        r: &GpuBuf,
        d: &mut GpuBuf,
        x: &mut GpuBuf,
    ) {
        let len = level.n * 3;
        assert!(
            r.len == len && d.len == len && x.len == len,
            "gpu: chebyshev length mismatch"
        );
        assert!(
            r.precision == level.precision
                && d.precision == level.precision
                && x.precision == level.precision,
            "gpu: precision mismatch"
        );
        if len == 0 {
            return;
        }
        self.dispatch(
            Kernel::ChebyshevStep,
            level.precision,
            Params::count(len)
                .alpha([alpha, 0.0, 0.0])
                .beta([beta, 0.0, 0.0]),
            &[&level.inv_diag.buffer, &r.buffer, &d.buffer, &x.buffer],
            self.groups_for(len),
        );
    }

    fn restrict(&self, agg: &GpuAggregates, fine: &GpuBuf, coarse: &mut GpuBuf) {
        assert_eq!(fine.len, agg.n_fine * 3, "gpu: fine length mismatch");
        assert_eq!(coarse.len, agg.n_coarse * 3, "gpu: coarse length mismatch");
        assert_eq!(fine.precision, coarse.precision, "gpu: precision mismatch");
        if agg.n_coarse == 0 {
            return;
        }
        self.dispatch(
            Kernel::Restrict,
            fine.precision,
            Params::count(agg.n_coarse).cols(3),
            &[
                &agg.coarse_offsets.buffer,
                &agg.members.buffer,
                &agg.agg.buffer,
                &fine.buffer,
                &coarse.buffer,
            ],
            self.groups_for(agg.n_coarse),
        );
    }

    fn prolong_add(&self, agg: &GpuAggregates, coarse: &GpuBuf, fine: &mut GpuBuf) {
        assert_eq!(fine.len, agg.n_fine * 3, "gpu: fine length mismatch");
        assert_eq!(coarse.len, agg.n_coarse * 3, "gpu: coarse length mismatch");
        assert_eq!(fine.precision, coarse.precision, "gpu: precision mismatch");
        if agg.n_fine == 0 {
            return;
        }
        self.dispatch(
            Kernel::ProlongAdd,
            fine.precision,
            Params::count(agg.n_fine).cols(3),
            &[
                &agg.coarse_offsets.buffer,
                &agg.members.buffer,
                &agg.agg.buffer,
                &coarse.buffer,
                &fine.buffer,
            ],
            self.groups_for(agg.n_fine),
        );
    }

    fn axpy(&self, alpha: [f64; 3], x: &GpuBuf, y: &mut GpuBuf) {
        assert_eq!(x.len, y.len, "gpu: axpy length mismatch");
        assert_eq!(x.precision, y.precision, "gpu: precision mismatch");
        if x.len == 0 {
            return;
        }
        self.dispatch(
            Kernel::Axpy,
            x.precision,
            Params::count(x.len).alpha(alpha),
            &[&x.buffer, &y.buffer],
            self.groups_for(x.len),
        );
    }

    fn scale(&self, alpha: [f64; 3], x: &mut GpuBuf) {
        if x.len == 0 {
            return;
        }
        self.dispatch(
            Kernel::Scale,
            x.precision,
            Params::count(x.len).alpha(alpha),
            &[&self.dummy, &x.buffer],
            self.groups_for(x.len),
        );
    }

    fn dot3(&self, a: &GpuBuf, b: &GpuBuf) -> [f64; 3] {
        self.reduce(Kernel::DotPartial, a, b)
    }

    fn norm3(&self, a: &GpuBuf) -> [f64; 3] {
        let s = self.reduce(Kernel::NormPartial, a, a);
        [s[0].sqrt(), s[1].sqrt(), s[2].sqrt()]
    }

    fn sync(&self) {
        let index = self.flush();
        self.wait(index);
    }
}

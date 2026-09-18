//! Device buffers: [`GpuBuf`] (one vector), the [`BufferPool`] that creates
//! storage buffers and accounts for device bytes, a persistent staging ring
//! for readback, and the uniform ring the per-dispatch [`Params`] go through.
//!
//! All device allocation happens here; kernels never allocate. The staging
//! ring grows only when a readback larger than any previous one is requested
//! (`BufferPool::reserve_staging` preallocates for a known vector size).
//!
//! Ordering rule: `queue.write_buffer` is applied at the start of the *next*
//! submit, so the backend submits any pending command encoder before writing
//! to a buffer that pending commands may read (see `GpuBackend::flush`).

use std::cell::{Cell, RefCell};
use std::num::NonZeroU64;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Arc;

use crate::linear_solver::Precision;
use crate::types::TheseusError;

/// Storage-buffer usage for vectors and graph arrays.
const STORAGE_USAGE: wgpu::BufferUsages = wgpu::BufferUsages::STORAGE
    .union(wgpu::BufferUsages::COPY_SRC)
    .union(wgpu::BufferUsages::COPY_DST);

/// Round a byte count up to `wgpu::COPY_BUFFER_ALIGNMENT` and at least 4
/// bytes (zero-sized bindings are invalid).
pub(crate) fn buffer_bytes(bytes: u64) -> u64 {
    let a = wgpu::COPY_BUFFER_ALIGNMENT;
    bytes.max(4).div_ceil(a) * a
}

/// Shared counter of live device bytes; decremented when a tracked buffer
/// drops so `BufferPool::device_bytes` is the current working set.
#[derive(Debug, Default)]
struct ByteCounter(AtomicU64);

impl ByteCounter {
    fn add(&self, bytes: u64) {
        self.0.fetch_add(bytes, Ordering::Relaxed);
    }

    fn sub(&self, bytes: u64) {
        self.0.fetch_sub(bytes, Ordering::Relaxed);
    }

    fn get(&self) -> u64 {
        self.0.load(Ordering::Relaxed)
    }
}

/// A storage buffer whose size is accounted in its pool for as long as it
/// lives.
#[derive(Debug)]
pub(crate) struct Tracked {
    pub(crate) buffer: wgpu::Buffer,
    counter: Arc<ByteCounter>,
}

impl Drop for Tracked {
    fn drop(&mut self) {
        self.counter.sub(self.buffer.size());
    }
}

impl std::ops::Deref for Tracked {
    type Target = wgpu::Buffer;

    fn deref(&self) -> &wgpu::Buffer {
        &self.buffer
    }
}

/// A vector buffer of `len` scalars in `precision`.
#[derive(Debug)]
pub struct GpuBuf {
    pub(crate) buffer: Tracked,
    pub(crate) len: usize,
    pub(crate) precision: Precision,
}

impl GpuBuf {
    /// Number of scalars.
    pub fn len(&self) -> usize {
        self.len
    }

    /// `true` when the buffer holds no scalars.
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Storage precision.
    pub fn precision(&self) -> Precision {
        self.precision
    }

    /// Allocated size in bytes.
    pub fn size_bytes(&self) -> u64 {
        self.buffer.size()
    }

    /// Bytes holding the `len` scalars (≤ `size_bytes`).
    pub(crate) fn data_bytes(&self) -> u64 {
        (self.len * self.precision.size_of()) as u64
    }

    /// The underlying wgpu buffer.
    pub fn raw(&self) -> &wgpu::Buffer {
        &self.buffer
    }
}

/// A `u32` index array (CSR offsets, edge ids, aggregate maps).
#[derive(Debug)]
pub struct IndexBuf {
    pub(crate) buffer: Tracked,
    pub(crate) len: usize,
}

impl IndexBuf {
    /// Number of `u32` entries.
    pub fn len(&self) -> usize {
        self.len
    }

    /// `true` when there are no entries.
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// The underlying wgpu buffer.
    pub fn raw(&self) -> &wgpu::Buffer {
        &self.buffer
    }
}

/// Uniform parameters of one dispatch, encoded on the host for the module's
/// `Real` type (see `kernels.wgsl`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Params {
    pub n: u32,
    pub cols: u32,
    pub flag: u32,
    pub alpha: [f64; 3],
    pub beta: [f64; 3],
}

impl Params {
    pub fn count(n: usize) -> Self {
        Self {
            n: u32::try_from(n).expect("gpu: element count exceeds u32"),
            cols: 3,
            flag: 0,
            alpha: [0.0; 3],
            beta: [0.0; 3],
        }
    }

    pub fn cols(mut self, cols: u32) -> Self {
        self.cols = cols;
        self
    }

    pub fn alpha(mut self, alpha: [f64; 3]) -> Self {
        self.alpha = alpha;
        self
    }

    pub fn beta(mut self, beta: [f64; 3]) -> Self {
        self.beta = beta;
        self
    }

    /// Size of the `Params` struct in the shader for `precision` (WGSL
    /// uniform layout: `vec4<f32>` is 16-byte aligned, `vec4<f64>` 32-byte).
    pub const fn encoded_size(precision: Precision) -> u64 {
        match precision {
            Precision::F32 => 48,
            Precision::F64 => 96,
        }
    }

    /// Encode for `precision`; returns the bytes and their length.
    pub fn encode(&self, precision: Precision) -> ([u8; 96], usize) {
        let mut out = [0u8; 96];
        out[0..4].copy_from_slice(&self.n.to_le_bytes());
        out[4..8].copy_from_slice(&self.cols.to_le_bytes());
        out[8..12].copy_from_slice(&self.flag.to_le_bytes());
        match precision {
            Precision::F32 => {
                for (k, a) in self.alpha.iter().enumerate() {
                    out[16 + 4 * k..20 + 4 * k].copy_from_slice(&(*a as f32).to_le_bytes());
                }
                for (k, b) in self.beta.iter().enumerate() {
                    out[32 + 4 * k..36 + 4 * k].copy_from_slice(&(*b as f32).to_le_bytes());
                }
                (out, 48)
            }
            Precision::F64 => {
                for (k, a) in self.alpha.iter().enumerate() {
                    out[32 + 8 * k..40 + 8 * k].copy_from_slice(&a.to_le_bytes());
                }
                for (k, b) in self.beta.iter().enumerate() {
                    out[64 + 8 * k..72 + 8 * k].copy_from_slice(&b.to_le_bytes());
                }
                (out, 96)
            }
        }
    }
}

/// Ring of uniform slots (one per dispatch) written with `queue.write_buffer`
/// and bound with a dynamic offset. When the ring is exhausted the backend
/// submits the pending work and resets it.
pub(crate) struct ParamsRing {
    buffer: wgpu::Buffer,
    slot_bytes: u32,
    slots: u32,
    next: Cell<u32>,
}

impl ParamsRing {
    /// Dispatches recordable between two submits.
    pub const SLOTS: u32 = 2048;

    fn new(device: &wgpu::Device, min_alignment: u32) -> Self {
        let align = min_alignment.max(1);
        let slot_bytes = 96u32.next_multiple_of(align);
        let buffer = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("theseus params ring"),
            size: u64::from(slot_bytes) * u64::from(Self::SLOTS),
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        Self {
            buffer,
            slot_bytes,
            slots: Self::SLOTS,
            next: Cell::new(0),
        }
    }

    pub fn size_bytes(&self) -> u64 {
        self.buffer.size()
    }

    pub fn is_full(&self) -> bool {
        self.next.get() >= self.slots
    }

    pub fn reset(&self) {
        self.next.set(0);
    }

    /// Write `params` into the next slot and return its dynamic offset.
    /// Callers check [`Self::is_full`] first.
    pub fn push(&self, queue: &wgpu::Queue, params: &Params, precision: Precision) -> u32 {
        let slot = self.next.get();
        assert!(
            slot < self.slots,
            "gpu: params ring exhausted without a flush"
        );
        self.next.set(slot + 1);
        let offset = slot * self.slot_bytes;
        let (bytes, len) = params.encode(precision);
        queue.write_buffer(&self.buffer, u64::from(offset), &bytes[..len]);
        offset
    }

    pub fn binding(&self, precision: Precision) -> wgpu::BindingResource<'_> {
        wgpu::BindingResource::Buffer(wgpu::BufferBinding {
            buffer: &self.buffer,
            offset: 0,
            size: NonZeroU64::new(Params::encoded_size(precision)),
        })
    }
}

/// Creates device buffers, converts host data to the buffer precision, and
/// accounts for device bytes.
pub struct BufferPool {
    device: wgpu::Device,
    queue: wgpu::Queue,
    device_bytes: Arc<ByteCounter>,
    /// Two mappable readback buffers (double buffering), grown on demand.
    staging: RefCell<[Option<wgpu::Buffer>; 2]>,
    staging_next: Cell<usize>,
    pub(crate) params: ParamsRing,
    max_buffer_size: u64,
    max_binding_size: u64,
    /// Host scratch for precision conversion on upload.
    scratch: RefCell<Vec<u8>>,
}

impl std::fmt::Debug for BufferPool {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("BufferPool")
            .field("device_bytes", &self.device_bytes.get())
            .finish()
    }
}

impl BufferPool {
    pub(crate) fn new(device: &wgpu::Device, queue: &wgpu::Queue, limits: &wgpu::Limits) -> Self {
        let params = ParamsRing::new(device, limits.min_uniform_buffer_offset_alignment);
        let pool = Self {
            device: device.clone(),
            queue: queue.clone(),
            device_bytes: Arc::new(ByteCounter::default()),
            staging: RefCell::new([None, None]),
            staging_next: Cell::new(0),
            params,
            max_buffer_size: limits.max_buffer_size,
            max_binding_size: limits.max_storage_buffer_binding_size,
            scratch: RefCell::new(Vec::new()),
        };
        pool.account(pool.params.size_bytes());
        pool
    }

    fn account(&self, bytes: u64) {
        self.device_bytes.add(bytes);
    }

    /// Device bytes currently held through this pool (vectors, index
    /// arrays, staging buffers, parameter ring); tracked buffers subtract
    /// themselves when dropped.
    pub fn device_bytes(&self) -> u64 {
        self.device_bytes.get()
    }

    /// Largest storage buffer a single binding may cover on this device.
    pub fn max_binding_size(&self) -> u64 {
        self.max_binding_size
    }

    /// Create a storage buffer of `bytes` (rounded up). Out-of-memory is
    /// caught through a wgpu error scope and reported as
    /// [`TheseusError::GpuOutOfMemory`] instead of being fatal.
    pub(crate) fn create_storage(&self, bytes: u64, label: &str) -> Result<Tracked, TheseusError> {
        let size = buffer_bytes(bytes);
        if size > self.max_binding_size {
            // Not a memory shortage: one array is larger than a single
            // storage binding on this adapter (lavapipe: 128 MiB). Splitting
            // arrays over several bindings is a WS-H item (§2.5).
            return Err(TheseusError::IterativeSolverUnsupported(format!(
                "GPU array {label:?} of {size} bytes exceeds the adapter's \
                 max_storage_buffer_binding_size of {} bytes ({} MiB); this network is too \
                 large for a single binding on this adapter",
                self.max_binding_size,
                self.max_binding_size >> 20
            )));
        }
        if size > self.max_buffer_size {
            return Err(TheseusError::GpuOutOfMemory {
                requested: size,
                available: self.max_buffer_size,
            });
        }
        let scope = self.device.push_error_scope(wgpu::ErrorFilter::OutOfMemory);
        let buffer = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some(label),
            size,
            usage: STORAGE_USAGE,
            mapped_at_creation: false,
        });
        if pollster::block_on(scope.pop()).is_some() {
            drop(buffer);
            return Err(TheseusError::GpuOutOfMemory {
                requested: size,
                available: self.max_buffer_size.saturating_sub(self.device_bytes.get()),
            });
        }
        self.account(size);
        Ok(Tracked {
            buffer,
            counter: Arc::clone(&self.device_bytes),
        })
    }

    /// Zeroed vector of `len` scalars in `precision`.
    pub(crate) fn create_vec(
        &self,
        len: usize,
        precision: Precision,
        label: &str,
    ) -> Result<GpuBuf, TheseusError> {
        // A runtime-sized `array<Real>` binding must hold at least one element.
        let buffer = self.create_storage((len.max(1) * precision.size_of()) as u64, label)?;
        Ok(GpuBuf {
            buffer,
            len,
            precision,
        })
    }

    /// Vector in `precision` initialised from `src` (converted on the host).
    pub(crate) fn create_vec_from(
        &self,
        src: &[f64],
        precision: Precision,
        label: &str,
    ) -> Result<GpuBuf, TheseusError> {
        let buf = self.create_vec(src.len(), precision, label)?;
        self.write_vec(&buf, src);
        Ok(buf)
    }

    /// `u32` index array initialised from `src`.
    pub(crate) fn create_u32(&self, src: &[u32], label: &str) -> Result<IndexBuf, TheseusError> {
        let buffer = self.create_storage((src.len() * 4) as u64, label)?;
        if !src.is_empty() {
            self.queue
                .write_buffer(&buffer, 0, bytemuck::cast_slice(src));
        }
        Ok(IndexBuf {
            buffer,
            len: src.len(),
        })
    }

    /// Convert `src` to the buffer precision and enqueue the write. The
    /// caller must have flushed any pending commands that read `dst`.
    pub(crate) fn write_vec(&self, dst: &GpuBuf, src: &[f64]) {
        assert_eq!(src.len(), dst.len, "gpu: upload length mismatch");
        if src.is_empty() {
            return;
        }
        let mut scratch = self.scratch.borrow_mut();
        match dst.precision {
            Precision::F32 => {
                scratch.clear();
                scratch.reserve(src.len() * 4);
                for v in src {
                    scratch.extend_from_slice(&(*v as f32).to_le_bytes());
                }
                self.queue.write_buffer(&dst.buffer, 0, &scratch);
            }
            Precision::F64 => {
                self.queue
                    .write_buffer(&dst.buffer, 0, bytemuck::cast_slice(src));
            }
        }
    }

    /// Decode `bytes` (device layout of `precision`) into `dst` as `f64`.
    pub(crate) fn decode(precision: Precision, bytes: &[u8], dst: &mut [f64]) {
        match precision {
            Precision::F32 => {
                for (d, c) in dst.iter_mut().zip(bytes.chunks_exact(4)) {
                    *d = f64::from(f32::from_le_bytes([c[0], c[1], c[2], c[3]]));
                }
            }
            Precision::F64 => {
                for (d, c) in dst.iter_mut().zip(bytes.chunks_exact(8)) {
                    *d = f64::from_le_bytes([c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7]]);
                }
            }
        }
    }

    /// Acquire a mappable staging buffer of at least `bytes`, alternating
    /// between the two slots. Grows a slot only when it is too small.
    pub(crate) fn staging(&self, bytes: u64) -> StagingHandle<'_> {
        let size = buffer_bytes(bytes);
        let idx = self.staging_next.get();
        self.staging_next.set((idx + 1) % 2);
        let mut slots = self.staging.borrow_mut();
        let needs_new = slots[idx].as_ref().is_none_or(|b| b.size() < size);
        if needs_new {
            let buffer = self.device.create_buffer(&wgpu::BufferDescriptor {
                label: Some("theseus staging"),
                size,
                usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
                mapped_at_creation: false,
            });
            if let Some(old) = &slots[idx] {
                self.device_bytes.sub(old.size());
            }
            self.account(size);
            slots[idx] = Some(buffer);
        }
        StagingHandle { pool: self, idx }
    }

    /// Preallocate both staging slots for readbacks of `bytes`.
    pub fn reserve_staging(&self, bytes: u64) {
        let _ = self.staging(bytes);
        let _ = self.staging(bytes);
    }
}

/// Borrow of one staging slot.
pub(crate) struct StagingHandle<'a> {
    pool: &'a BufferPool,
    idx: usize,
}

impl StagingHandle<'_> {
    pub fn with_buffer<R>(&self, f: impl FnOnce(&wgpu::Buffer) -> R) -> R {
        let slots = self.pool.staging.borrow();
        f(slots[self.idx].as_ref().expect("staging slot allocated"))
    }
}

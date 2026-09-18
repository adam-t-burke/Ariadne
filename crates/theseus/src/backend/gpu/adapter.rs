//! Adapter enumeration, selection policy and device creation (§2.5).
//!
//! Selection: enumerate every adapter of the backends `WGPU_BACKEND` allows,
//! rank them by the caller's [`AdapterPreference`] (discrete → integrated →
//! virtual → other → cpu, or integrated first), drop software rasterisers
//! (`DeviceType::Cpu`) unless `THESEUS_GPU_ALLOW_SOFTWARE=1`, and honour a
//! `WGPU_ADAPTER_NAME` substring filter with the same case-insensitive
//! semantics as `wgpu::util::initialize_adapter_from_env`. The device is
//! requested with the adapter's own limits (the wgpu defaults cap storage
//! bindings at 128 MiB, below a 10M-node `f32` vector) and with
//! `Features::SHADER_F64` when the adapter offers it.

use std::fmt::Write as _;

use crate::backend::{GpuAdapterReport, GpuProbe};
use crate::linear_solver::AdapterPreference;
use crate::types::TheseusError;

/// Environment variable that admits `DeviceType::Cpu` adapters (lavapipe,
/// SwiftShader, WARP) for CI correctness runs.
pub const ALLOW_SOFTWARE_ENV: &str = "THESEUS_GPU_ALLOW_SOFTWARE";

/// How adapters are filtered and ranked.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AdapterPolicy {
    /// Device-class ranking.
    pub preference: AdapterPreference,
    /// Accept software rasterisers.
    pub allow_software: bool,
    /// Case-insensitive substring the adapter name must contain
    /// (`WGPU_ADAPTER_NAME`).
    pub name_filter: Option<String>,
}

impl Default for AdapterPolicy {
    fn default() -> Self {
        Self {
            preference: AdapterPreference::Discrete,
            allow_software: false,
            name_filter: None,
        }
    }
}

impl AdapterPolicy {
    /// Default policy with the `THESEUS_GPU_ALLOW_SOFTWARE` and
    /// `WGPU_ADAPTER_NAME` overrides applied.
    pub fn from_env() -> Self {
        Self::default().with_env()
    }

    /// Policy for a given preference, with the environment overrides applied.
    pub fn for_preference(preference: AdapterPreference) -> Self {
        Self {
            preference,
            ..Self::default()
        }
        .with_env()
    }

    fn with_env(mut self) -> Self {
        if let Ok(v) = std::env::var(ALLOW_SOFTWARE_ENV) {
            let v = v.trim().to_ascii_lowercase();
            self.allow_software = matches!(v.as_str(), "1" | "true" | "yes" | "on");
        }
        if let Ok(name) = std::env::var("WGPU_ADAPTER_NAME") {
            if !name.is_empty() {
                self.name_filter = Some(name);
            }
        }
        self
    }

    /// Rank of a device type under this policy (lower is better).
    fn rank(&self, device_type: wgpu::DeviceType) -> u8 {
        use wgpu::DeviceType::*;
        match self.preference {
            AdapterPreference::Discrete => match device_type {
                DiscreteGpu => 0,
                IntegratedGpu => 1,
                VirtualGpu => 2,
                Other => 3,
                Cpu => 4,
            },
            AdapterPreference::Integrated => match device_type {
                IntegratedGpu => 0,
                DiscreteGpu => 1,
                VirtualGpu => 2,
                Other => 3,
                Cpu => 4,
            },
            AdapterPreference::Any => 0,
        }
    }

    fn accepts(&self, info: &wgpu::AdapterInfo) -> bool {
        if info.device_type == wgpu::DeviceType::Cpu && !self.allow_software {
            return false;
        }
        match &self.name_filter {
            Some(filter) => info.name.to_lowercase().contains(&filter.to_lowercase()),
            None => true,
        }
    }
}

fn device_type_name(t: wgpu::DeviceType) -> &'static str {
    match t {
        wgpu::DeviceType::DiscreteGpu => "discrete",
        wgpu::DeviceType::IntegratedGpu => "integrated",
        wgpu::DeviceType::VirtualGpu => "virtual",
        wgpu::DeviceType::Cpu => "cpu",
        wgpu::DeviceType::Other => "other",
    }
}

fn report_for(adapter: &wgpu::Adapter, policy: &AdapterPolicy) -> GpuAdapterReport {
    let info = adapter.get_info();
    let limits = adapter.limits();
    let driver = if info.driver_info.is_empty() {
        info.driver.clone()
    } else {
        format!("{} {}", info.driver, info.driver_info)
    };
    GpuAdapterReport {
        name: info.name.clone(),
        backend: info.backend.to_str().to_string(),
        device_type: device_type_name(info.device_type).to_string(),
        driver,
        max_storage_buffer_binding_size: limits.max_storage_buffer_binding_size,
        max_buffer_size: limits.max_buffer_size,
        shader_f64: adapter.features().contains(wgpu::Features::SHADER_F64),
        software: info.device_type == wgpu::DeviceType::Cpu,
        selectable: policy.accepts(&info),
    }
}

/// wgpu instance honouring `WGPU_BACKEND` (and the other `WGPU_*` instance
/// variables). Compute-only: no surface is ever created.
pub fn instance() -> wgpu::Instance {
    wgpu::Instance::new(wgpu::InstanceDescriptor::new_without_display_handle_from_env())
}

/// Adapters ranked by `policy`, each paired with its report, and the index
/// of the first selectable one.
pub fn enumerate(
    instance: &wgpu::Instance,
    policy: &AdapterPolicy,
) -> (Vec<(wgpu::Adapter, GpuAdapterReport)>, Option<usize>) {
    let adapters = pollster::block_on(instance.enumerate_adapters(wgpu::Backends::all()));
    let mut ranked: Vec<(u8, usize, wgpu::Adapter)> = adapters
        .into_iter()
        .enumerate()
        .map(|(i, a)| (policy.rank(a.get_info().device_type), i, a))
        .collect();
    // Stable by (rank, enumeration order) so the choice is reproducible.
    ranked.sort_by_key(|(rank, i, _)| (*rank, *i));
    let list: Vec<(wgpu::Adapter, GpuAdapterReport)> = ranked
        .into_iter()
        .map(|(_, _, a)| {
            let report = report_for(&a, policy);
            (a, report)
        })
        .collect();
    let chosen = list.iter().position(|(_, r)| r.selectable);
    (list, chosen)
}

fn summary(adapters: &[GpuAdapterReport], chosen: Option<usize>, policy: &AdapterPolicy) -> String {
    let mut s = String::new();
    match chosen {
        Some(i) => {
            let a = &adapters[i];
            let _ = write!(
                s,
                "using {} ({}, {}, {}; storage binding ≤ {} MiB, buffer ≤ {} MiB, f64 {})",
                a.name,
                a.backend,
                a.device_type,
                a.driver,
                a.max_storage_buffer_binding_size >> 20,
                a.max_buffer_size >> 20,
                if a.shader_f64 { "yes" } else { "no" }
            );
        }
        None if adapters.is_empty() => {
            s.push_str("no GPU adapter found");
            if let Ok(b) = std::env::var("WGPU_BACKEND") {
                let _ = write!(s, " (WGPU_BACKEND={b})");
            }
        }
        None => {
            s.push_str("no usable GPU adapter");
            if let Some(f) = &policy.name_filter {
                let _ = write!(s, " matching WGPU_ADAPTER_NAME={f:?}");
            }
            if adapters.iter().all(|a| a.software) && !policy.allow_software {
                let _ = write!(
                    s,
                    "; only software adapters were found (set {ALLOW_SOFTWARE_ENV}=1 to allow them)"
                );
            }
        }
    }
    if !adapters.is_empty() {
        s.push_str("; adapters seen: ");
        for (i, a) in adapters.iter().enumerate() {
            if i > 0 {
                s.push_str(", ");
            }
            let _ = write!(s, "{} [{} {}]", a.name, a.backend, a.device_type);
            if !a.selectable {
                s.push_str(" (rejected)");
            }
        }
    }
    s
}

/// Enumerate adapters under `policy` and build the report.
pub fn probe(policy: &AdapterPolicy) -> GpuProbe {
    let instance = instance();
    let (list, chosen) = enumerate(&instance, policy);
    let adapters: Vec<GpuAdapterReport> = list.into_iter().map(|(_, r)| r).collect();
    let message = summary(&adapters, chosen, policy);
    GpuProbe {
        built_with_gpu_feature: true,
        adapters,
        chosen,
        message,
    }
}

/// An opened device: the wgpu handles plus the facts the backend needs
/// (limits, `f64` support, the probe report it was created from).
pub struct GpuContext {
    pub(crate) instance: wgpu::Instance,
    pub(crate) adapter: wgpu::Adapter,
    pub(crate) device: wgpu::Device,
    pub(crate) queue: wgpu::Queue,
    /// Report of the adapter that was opened.
    pub report: GpuAdapterReport,
    /// Full probe (all adapters seen) at creation time.
    pub probe: GpuProbe,
    /// Whether the device was created with `SHADER_F64`.
    pub shader_f64: bool,
    /// Limits the device was created with.
    pub limits: wgpu::Limits,
}

impl std::fmt::Debug for GpuContext {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("GpuContext")
            .field("adapter", &self.report.name)
            .field("backend", &self.report.backend)
            .field("shader_f64", &self.shader_f64)
            .finish()
    }
}

impl GpuContext {
    /// Open the adapter `policy` selects.
    ///
    /// Errors with [`TheseusError::GpuUnavailable`] (message = probe summary)
    /// when no adapter is usable or the device request fails.
    pub fn new(policy: &AdapterPolicy) -> Result<Self, TheseusError> {
        let instance = instance();
        let (mut list, chosen) = enumerate(&instance, policy);
        let adapters: Vec<GpuAdapterReport> = list.iter().map(|(_, r)| r.clone()).collect();
        let message = summary(&adapters, chosen, policy);
        let probe = GpuProbe {
            built_with_gpu_feature: true,
            adapters,
            chosen,
            message: message.clone(),
        };
        let Some(i) = chosen else {
            return Err(TheseusError::GpuUnavailable(message));
        };
        let (adapter, report) = list.swap_remove(i);

        let shader_f64 = report.shader_f64;
        let mut required_features = wgpu::Features::empty();
        if shader_f64 {
            required_features |= wgpu::Features::SHADER_F64;
        }
        let limits = adapter.limits();
        let desc = wgpu::DeviceDescriptor {
            label: Some("theseus"),
            required_features,
            required_limits: limits.clone(),
            experimental_features: wgpu::ExperimentalFeatures::disabled(),
            memory_hints: wgpu::MemoryHints::Performance,
            trace: wgpu::Trace::Off,
        };
        let (device, queue) = pollster::block_on(adapter.request_device(&desc)).map_err(|e| {
            TheseusError::GpuUnavailable(format!(
                "device request failed on {} ({}): {e}",
                report.name, report.backend
            ))
        })?;
        Ok(Self {
            instance,
            adapter,
            device,
            queue,
            report,
            probe,
            shader_f64,
            limits,
        })
    }

    /// Open the adapter selected by the default policy and the environment.
    pub fn from_env() -> Result<Self, TheseusError> {
        Self::new(&AdapterPolicy::from_env())
    }

    /// The wgpu device.
    pub fn device(&self) -> &wgpu::Device {
        &self.device
    }

    /// The wgpu queue.
    pub fn queue(&self) -> &wgpu::Queue {
        &self.queue
    }

    /// The adapter the device was created from.
    pub fn adapter(&self) -> &wgpu::Adapter {
        &self.adapter
    }

    /// The instance (kept alive for the device's lifetime).
    pub fn instance(&self) -> &wgpu::Instance {
        &self.instance
    }
}

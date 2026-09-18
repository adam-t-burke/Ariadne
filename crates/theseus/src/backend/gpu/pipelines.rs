//! Shader compilation and compute pipelines.
//!
//! One [`KernelSet`] is built per (precision, workgroup size): the WGSL source
//! in `kernels.wgsl` is specialised by string substitution of the `Real`
//! alias (`f32` / `f64`) and the `WG` constant, then compiled once. The `f64`
//! set is created only on devices with `Features::SHADER_F64` (naga accepts
//! the `f64` type without an `enable` directive; wgpu validates the feature).

use crate::linear_solver::Precision;

/// WGSL source of every kernel (before specialisation).
pub const KERNELS_WGSL: &str = include_str!("kernels.wgsl");

const REAL_LINE: &str = "alias Real = f32;";
const WG_LINE: &str = "const WG: u32 = 64u;";

/// Workgroup sizes the Stage-D tuning sweeps over.
pub const WORKGROUP_SIZES: [u32; 3] = [64, 128, 256];

/// Specialise the kernel source for `precision` and workgroup size `wg`.
pub fn shader_source(precision: Precision, wg: u32) -> String {
    debug_assert!(KERNELS_WGSL.contains(REAL_LINE) && KERNELS_WGSL.contains(WG_LINE));
    let real = match precision {
        Precision::F32 => "alias Real = f32;",
        Precision::F64 => "alias Real = f64;",
    };
    KERNELS_WGSL.replacen(REAL_LINE, real, 1).replacen(
        WG_LINE,
        &format!("const WG: u32 = {wg}u;"),
        1,
    )
}

/// Bind-group families: which storage bindings of group 0 an entry point
/// uses (binding numbers match `kernels.wgsl`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub(crate) enum Family {
    /// offsets, edges, other, weight, anchor, x, b, out (1..=8).
    Graph,
    /// inv_diag, r, d, x (9..=12).
    Chebyshev,
    /// coarse_offsets, members, agg, src, dst (13..=17).
    Transfer,
    /// x, y (18, 19).
    Vector,
    /// a, b, partials (20..=22).
    Reduce,
    /// offsets, fine_edges, fine_weight, coarse_weight (23..=26).
    CoarseWeights,
    /// row_ptr, col_idx, values, x, b, out (27..=32).
    Csr,
}

impl Family {
    /// `(binding, read_only)` storage entries of the family.
    pub(crate) fn bindings(self) -> &'static [(u32, bool)] {
        match self {
            Family::Graph => &[
                (1, true),
                (2, true),
                (3, true),
                (4, true),
                (5, true),
                (6, true),
                (7, true),
                (8, false),
            ],
            Family::Chebyshev => &[(9, true), (10, true), (11, false), (12, false)],
            Family::Transfer => &[(13, true), (14, true), (15, true), (16, true), (17, false)],
            Family::Vector => &[(18, true), (19, false)],
            Family::Reduce => &[(20, true), (21, true), (22, false)],
            Family::CoarseWeights => &[(23, true), (24, true), (25, true), (26, false)],
            Family::Csr => &[
                (27, true),
                (28, true),
                (29, true),
                (30, true),
                (31, true),
                (32, false),
            ],
        }
    }

    const ALL: [Family; 7] = [
        Family::Graph,
        Family::Chebyshev,
        Family::Transfer,
        Family::Vector,
        Family::Reduce,
        Family::CoarseWeights,
        Family::Csr,
    ];

    fn index(self) -> usize {
        Self::ALL
            .iter()
            .position(|f| *f == self)
            .expect("family listed")
    }
}

/// Entry points of `kernels.wgsl`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub(crate) enum Kernel {
    ApplyGraph,
    Residual,
    InvDiag,
    ChebyshevStep,
    Restrict,
    ProlongAdd,
    Axpy,
    Scale,
    DotPartial,
    NormPartial,
    CoarseWeightUpdate,
    ApplyCsr,
    ResidualCsr,
}

impl Kernel {
    const ALL: [Kernel; 13] = [
        Kernel::ApplyGraph,
        Kernel::Residual,
        Kernel::InvDiag,
        Kernel::ChebyshevStep,
        Kernel::Restrict,
        Kernel::ProlongAdd,
        Kernel::Axpy,
        Kernel::Scale,
        Kernel::DotPartial,
        Kernel::NormPartial,
        Kernel::CoarseWeightUpdate,
        Kernel::ApplyCsr,
        Kernel::ResidualCsr,
    ];

    /// WGSL entry-point name.
    pub(crate) fn entry_point(self) -> &'static str {
        match self {
            Kernel::ApplyGraph => "apply_graph",
            Kernel::Residual => "residual",
            Kernel::InvDiag => "inv_diag",
            Kernel::ChebyshevStep => "chebyshev_step",
            Kernel::Restrict => "restrict_sum",
            Kernel::ProlongAdd => "prolong_add",
            Kernel::Axpy => "axpy",
            Kernel::Scale => "scale",
            Kernel::DotPartial => "dot_partial",
            Kernel::NormPartial => "norm_partial",
            Kernel::CoarseWeightUpdate => "coarse_weight_update",
            Kernel::ApplyCsr => "apply_csr",
            Kernel::ResidualCsr => "residual_csr",
        }
    }

    pub(crate) fn family(self) -> Family {
        match self {
            Kernel::ApplyGraph | Kernel::Residual | Kernel::InvDiag => Family::Graph,
            Kernel::ChebyshevStep => Family::Chebyshev,
            Kernel::Restrict | Kernel::ProlongAdd => Family::Transfer,
            Kernel::Axpy | Kernel::Scale => Family::Vector,
            Kernel::DotPartial | Kernel::NormPartial => Family::Reduce,
            Kernel::CoarseWeightUpdate => Family::CoarseWeights,
            Kernel::ApplyCsr | Kernel::ResidualCsr => Family::Csr,
        }
    }

    fn index(self) -> usize {
        Self::ALL
            .iter()
            .position(|k| *k == self)
            .expect("kernel listed")
    }
}

/// Compiled pipelines of one precision at one workgroup size.
pub struct KernelSet {
    precision: Precision,
    workgroup_size: u32,
    layouts: Vec<wgpu::BindGroupLayout>,
    pipelines: Vec<wgpu::ComputePipeline>,
}

impl std::fmt::Debug for KernelSet {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("KernelSet")
            .field("precision", &self.precision)
            .field("workgroup_size", &self.workgroup_size)
            .finish()
    }
}

impl KernelSet {
    /// Compile the kernels for `precision` at workgroup size `wg`.
    pub(crate) fn new(device: &wgpu::Device, precision: Precision, wg: u32) -> Self {
        let source = shader_source(precision, wg);
        let label = format!("theseus kernels {precision:?} wg{wg}");
        let module = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some(&label),
            source: wgpu::ShaderSource::Wgsl(source.into()),
        });

        let layouts: Vec<wgpu::BindGroupLayout> = Family::ALL
            .iter()
            .map(|family| {
                let mut entries = vec![wgpu::BindGroupLayoutEntry {
                    binding: 0,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Uniform,
                        has_dynamic_offset: true,
                        min_binding_size: None,
                    },
                    count: None,
                }];
                entries.extend(family.bindings().iter().map(|&(binding, read_only)| {
                    wgpu::BindGroupLayoutEntry {
                        binding,
                        visibility: wgpu::ShaderStages::COMPUTE,
                        ty: wgpu::BindingType::Buffer {
                            ty: wgpu::BufferBindingType::Storage { read_only },
                            has_dynamic_offset: false,
                            min_binding_size: None,
                        },
                        count: None,
                    }
                }));
                device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
                    label: Some(&format!("theseus {family:?}")),
                    entries: &entries,
                })
            })
            .collect();

        let pipelines = Kernel::ALL
            .iter()
            .map(|kernel| {
                let layout = &layouts[kernel.family().index()];
                let pipeline_layout =
                    device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
                        label: Some(kernel.entry_point()),
                        bind_group_layouts: &[Some(layout)],
                        immediate_size: 0,
                    });
                device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
                    label: Some(kernel.entry_point()),
                    layout: Some(&pipeline_layout),
                    module: &module,
                    entry_point: Some(kernel.entry_point()),
                    compilation_options: wgpu::PipelineCompilationOptions::default(),
                    cache: None,
                })
            })
            .collect();

        Self {
            precision,
            workgroup_size: wg,
            layouts,
            pipelines,
        }
    }

    /// Precision of the `Real` type these pipelines were compiled with.
    pub fn precision(&self) -> Precision {
        self.precision
    }

    /// Workgroup size these pipelines were compiled with.
    pub fn workgroup_size(&self) -> u32 {
        self.workgroup_size
    }

    pub(crate) fn layout(&self, family: Family) -> &wgpu::BindGroupLayout {
        &self.layouts[family.index()]
    }

    pub(crate) fn pipeline(&self, kernel: Kernel) -> &wgpu::ComputePipeline {
        &self.pipelines[kernel.index()]
    }
}

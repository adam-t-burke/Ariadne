// Theseus iterative-solver kernels (program plan §2.5).
//
// Vectors are blocks of three columns stored row-major, `v[node * 3 + k]`,
// as flat `array<Real>` storage buffers (no vec3 padding). Every kernel is
// one thread per node / element with a grid-stride loop so that sizes above
// `max_compute_workgroups_per_dimension * WG` still dispatch in one call.
//
// `Real` and `WG` are substituted by the host at pipeline creation
// (`pipelines::shader_source`): `Real` is `f32` or, on adapters with
// `SHADER_F64`, `f64`; `WG` is the workgroup size (64 / 128 / 256).
// Binding numbers are unique across the module so that each entry point's
// bind-group layout is a strict subset of group 0.

alias Real = f32;
const WG: u32 = 64u;

struct Params {
    // Primary element count of the dispatch (nodes, scalars or edges).
    n: u32,
    // Columns per node for the transfer kernels (1 or 3).
    cols: u32,
    // Kernel-specific flag.
    flag: u32,
    _pad: u32,
    // Per-column scalars (x, y, z; w unused). `alpha.x` is the Chebyshev
    // alpha, `beta.x` the Chebyshev beta.
    alpha: vec4<Real>,
    beta: vec4<Real>,
}

@group(0) @binding(0) var<uniform> params: Params;

// ── Level operator (apply_graph / residual / inv_diag) ───────────────
@group(0) @binding(1) var<storage, read> g_offsets: array<u32>;
@group(0) @binding(2) var<storage, read> g_edges: array<u32>;
@group(0) @binding(3) var<storage, read> g_other: array<u32>;
@group(0) @binding(4) var<storage, read> g_weight: array<Real>;
@group(0) @binding(5) var<storage, read> g_anchor: array<Real>;
@group(0) @binding(6) var<storage, read> g_x: array<Real>;
@group(0) @binding(7) var<storage, read> g_b: array<Real>;
@group(0) @binding(8) var<storage, read_write> g_out: array<Real>;

// ── Chebyshev step ───────────────────────────────────────────────────
@group(0) @binding(9) var<storage, read> c_inv_diag: array<Real>;
@group(0) @binding(10) var<storage, read> c_r: array<Real>;
@group(0) @binding(11) var<storage, read_write> c_d: array<Real>;
@group(0) @binding(12) var<storage, read_write> c_x: array<Real>;

// ── Transfers ────────────────────────────────────────────────────────
@group(0) @binding(13) var<storage, read> t_coarse_offsets: array<u32>;
@group(0) @binding(14) var<storage, read> t_members: array<u32>;
@group(0) @binding(15) var<storage, read> t_agg: array<u32>;
@group(0) @binding(16) var<storage, read> t_src: array<Real>;
@group(0) @binding(17) var<storage, read_write> t_dst: array<Real>;

// ── Vector updates ───────────────────────────────────────────────────
@group(0) @binding(18) var<storage, read> v_x: array<Real>;
@group(0) @binding(19) var<storage, read_write> v_y: array<Real>;

// ── Reductions ───────────────────────────────────────────────────────
@group(0) @binding(20) var<storage, read> r_a: array<Real>;
@group(0) @binding(21) var<storage, read> r_b: array<Real>;
@group(0) @binding(22) var<storage, read_write> r_partials: array<Real>;

// ── Coarse weight update ─────────────────────────────────────────────
@group(0) @binding(23) var<storage, read> w_offsets: array<u32>;
@group(0) @binding(24) var<storage, read> w_fine_edges: array<u32>;
@group(0) @binding(25) var<storage, read> w_fine_weight: array<Real>;
@group(0) @binding(26) var<storage, read_write> w_coarse_weight: array<Real>;

// ── CSR level matrices (apply_csr / residual_csr) ────────────────────
@group(0) @binding(27) var<storage, read> m_row_ptr: array<u32>;
@group(0) @binding(28) var<storage, read> m_col_idx: array<u32>;
@group(0) @binding(29) var<storage, read> m_values: array<Real>;
@group(0) @binding(30) var<storage, read> m_x: array<Real>;
@group(0) @binding(31) var<storage, read> m_b: array<Real>;
@group(0) @binding(32) var<storage, read_write> m_out: array<Real>;

// (A x)_u for the three columns, summed in the same order as the scalar
// reference `LevelGraph::apply`: anchor term first, then incident edges in
// CSR order.
fn apply_row(u: u32) -> vec3<Real> {
    let xu = vec3<Real>(g_x[3u * u], g_x[3u * u + 1u], g_x[3u * u + 2u]);
    var acc = g_anchor[u] * xu;
    let begin = g_offsets[u];
    let end = g_offsets[u + 1u];
    for (var j = begin; j < end; j++) {
        let w = g_weight[g_edges[j]];
        let v = g_other[j];
        let xv = vec3<Real>(g_x[3u * v], g_x[3u * v + 1u], g_x[3u * v + 2u]);
        acc += w * (xu - xv);
    }
    return acc;
}

// y = A x
@compute @workgroup_size(WG)
fn apply_graph(
    @builtin(global_invocation_id) gid: vec3<u32>,
    @builtin(num_workgroups) nwg: vec3<u32>,
) {
    let stride = nwg.x * WG;
    for (var u = gid.x; u < params.n; u += stride) {
        let y = apply_row(u);
        g_out[3u * u] = y.x;
        g_out[3u * u + 1u] = y.y;
        g_out[3u * u + 2u] = y.z;
    }
}

// r = b − A x
@compute @workgroup_size(WG)
fn residual(
    @builtin(global_invocation_id) gid: vec3<u32>,
    @builtin(num_workgroups) nwg: vec3<u32>,
) {
    let stride = nwg.x * WG;
    for (var u = gid.x; u < params.n; u += stride) {
        let y = apply_row(u);
        g_out[3u * u] = g_b[3u * u] - y.x;
        g_out[3u * u + 1u] = g_b[3u * u + 1u] - y.y;
        g_out[3u * u + 2u] = g_b[3u * u + 2u] - y.z;
    }
}

// inv_diag_u = 1 / (anchor_u + Σ_{e ∋ u} w_e); 0 for an empty diagonal.
// Uses the operator bindings (offsets, edges, weight, anchor) and writes
// `g_out` (bound to the level's inv_diag buffer).
@compute @workgroup_size(WG)
fn inv_diag(
    @builtin(global_invocation_id) gid: vec3<u32>,
    @builtin(num_workgroups) nwg: vec3<u32>,
) {
    let stride = nwg.x * WG;
    for (var u = gid.x; u < params.n; u += stride) {
        var d = g_anchor[u];
        let begin = g_offsets[u];
        let end = g_offsets[u + 1u];
        for (var j = begin; j < end; j++) {
            d += g_weight[g_edges[j]];
        }
        if (d > Real(0)) {
            g_out[u] = Real(1) / d;
        } else {
            g_out[u] = Real(0);
        }
    }
}

// d = alpha · D⁻¹ r + beta · d;  x += d      (params.n = 3 × nodes)
@compute @workgroup_size(WG)
fn chebyshev_step(
    @builtin(global_invocation_id) gid: vec3<u32>,
    @builtin(num_workgroups) nwg: vec3<u32>,
) {
    let stride = nwg.x * WG;
    let alpha = params.alpha.x;
    let beta = params.beta.x;
    for (var i = gid.x; i < params.n; i += stride) {
        let node = i / 3u;
        let dn = alpha * c_inv_diag[node] * c_r[i] + beta * c_d[i];
        c_d[i] = dn;
        c_x[i] = c_x[i] + dn;
    }
}

// coarse[U] = Σ_{i ∈ members(U)} fine[i], one thread per coarse node, members
// in ascending fine order (deterministic, no atomics). params.n = n_coarse,
// params.cols = columns per node.
@compute @workgroup_size(WG)
fn restrict_sum(
    @builtin(global_invocation_id) gid: vec3<u32>,
    @builtin(num_workgroups) nwg: vec3<u32>,
) {
    let stride = nwg.x * WG;
    let cols = params.cols;
    for (var cu = gid.x; cu < params.n; cu += stride) {
        let begin = t_coarse_offsets[cu];
        let end = t_coarse_offsets[cu + 1u];
        for (var k = 0u; k < cols; k++) {
            var acc = Real(0);
            for (var j = begin; j < end; j++) {
                acc += t_src[t_members[j] * cols + k];
            }
            t_dst[cu * cols + k] = acc;
        }
    }
}

// fine[i] += coarse[agg[i]], one thread per fine node. params.n = n_fine.
@compute @workgroup_size(WG)
fn prolong_add(
    @builtin(global_invocation_id) gid: vec3<u32>,
    @builtin(num_workgroups) nwg: vec3<u32>,
) {
    let stride = nwg.x * WG;
    let cols = params.cols;
    for (var i = gid.x; i < params.n; i += stride) {
        let cu = t_agg[i];
        for (var k = 0u; k < cols; k++) {
            t_dst[i * cols + k] = t_dst[i * cols + k] + t_src[cu * cols + k];
        }
    }
}

// y[i] += alpha[i mod 3] · x[i]      (params.n = 3 × nodes)
@compute @workgroup_size(WG)
fn axpy(
    @builtin(global_invocation_id) gid: vec3<u32>,
    @builtin(num_workgroups) nwg: vec3<u32>,
) {
    let stride = nwg.x * WG;
    for (var i = gid.x; i < params.n; i += stride) {
        let k = i % 3u;
        v_y[i] = v_y[i] + params.alpha[k] * v_x[i];
    }
}

// y[i] *= alpha[i mod 3]              (params.n = 3 × nodes)
@compute @workgroup_size(WG)
fn scale(
    @builtin(global_invocation_id) gid: vec3<u32>,
    @builtin(num_workgroups) nwg: vec3<u32>,
) {
    let stride = nwg.x * WG;
    for (var i = gid.x; i < params.n; i += stride) {
        let k = i % 3u;
        v_y[i] = v_y[i] * params.alpha[k];
    }
}

// Workgroup reduction of per-column products. Each thread accumulates a
// fixed, grid-strided subset of nodes; the workgroup then reduces in shared
// memory with a fixed tree, and thread 0 writes the three partial sums to
// `r_partials[3 · workgroup + k]`. The host finishes the sum in f64 in
// workgroup order, so the result is deterministic for a given (n, WG,
// workgroup count). params.n = nodes.
var<workgroup> red: array<vec3<Real>, WG>;

fn reduce_write(
    acc: vec3<Real>,
    lid: u32,
    wid: u32,
) {
    red[lid] = acc;
    workgroupBarrier();
    for (var s = WG / 2u; s > 0u; s = s >> 1u) {
        if (lid < s) {
            red[lid] = red[lid] + red[lid + s];
        }
        workgroupBarrier();
    }
    if (lid == 0u) {
        r_partials[3u * wid] = red[0].x;
        r_partials[3u * wid + 1u] = red[0].y;
        r_partials[3u * wid + 2u] = red[0].z;
    }
}

@compute @workgroup_size(WG)
fn dot_partial(
    @builtin(global_invocation_id) gid: vec3<u32>,
    @builtin(local_invocation_id) lid: vec3<u32>,
    @builtin(workgroup_id) wid: vec3<u32>,
    @builtin(num_workgroups) nwg: vec3<u32>,
) {
    let stride = nwg.x * WG;
    var acc = vec3<Real>(Real(0), Real(0), Real(0));
    for (var u = gid.x; u < params.n; u += stride) {
        let a = vec3<Real>(r_a[3u * u], r_a[3u * u + 1u], r_a[3u * u + 2u]);
        let b = vec3<Real>(r_b[3u * u], r_b[3u * u + 1u], r_b[3u * u + 2u]);
        acc += a * b;
    }
    reduce_write(acc, lid.x, wid.x);
}

@compute @workgroup_size(WG)
fn norm_partial(
    @builtin(global_invocation_id) gid: vec3<u32>,
    @builtin(local_invocation_id) lid: vec3<u32>,
    @builtin(workgroup_id) wid: vec3<u32>,
    @builtin(num_workgroups) nwg: vec3<u32>,
) {
    let stride = nwg.x * WG;
    var acc = vec3<Real>(Real(0), Real(0), Real(0));
    for (var u = gid.x; u < params.n; u += stride) {
        let a = vec3<Real>(r_a[3u * u], r_a[3u * u + 1u], r_a[3u * u + 2u]);
        acc += a * a;
    }
    reduce_write(acc, lid.x, wid.x);
}

// (M x)_u for the three columns of a CSR matrix, entries in stored
// (ascending column) order — the same order as the CPU `csr_row`.
fn csr_row(u: u32) -> vec3<Real> {
    var acc = vec3<Real>(Real(0), Real(0), Real(0));
    let begin = m_row_ptr[u];
    let end = m_row_ptr[u + 1u];
    for (var j = begin; j < end; j++) {
        let c = 3u * m_col_idx[j];
        let v = m_values[j];
        acc += v * vec3<Real>(m_x[c], m_x[c + 1u], m_x[c + 2u]);
    }
    return acc;
}

// y = M x (params.flag = 0) or y += M x (params.flag = 1), one thread per
// row of M (params.n = rows). Used for the coarse operators and for the
// smoothed prolongation / restriction (`P`, `Pᵀ` as CSR).
@compute @workgroup_size(WG)
fn apply_csr(
    @builtin(global_invocation_id) gid: vec3<u32>,
    @builtin(num_workgroups) nwg: vec3<u32>,
) {
    let stride = nwg.x * WG;
    let accumulate = params.flag != 0u;
    for (var u = gid.x; u < params.n; u += stride) {
        let y = csr_row(u);
        if (accumulate) {
            m_out[3u * u] = m_out[3u * u] + y.x;
            m_out[3u * u + 1u] = m_out[3u * u + 1u] + y.y;
            m_out[3u * u + 2u] = m_out[3u * u + 2u] + y.z;
        } else {
            m_out[3u * u] = y.x;
            m_out[3u * u + 1u] = y.y;
            m_out[3u * u + 2u] = y.z;
        }
    }
}

// r = b − M x for a square CSR matrix (params.n = rows).
@compute @workgroup_size(WG)
fn residual_csr(
    @builtin(global_invocation_id) gid: vec3<u32>,
    @builtin(num_workgroups) nwg: vec3<u32>,
) {
    let stride = nwg.x * WG;
    for (var u = gid.x; u < params.n; u += stride) {
        let y = csr_row(u);
        m_out[3u * u] = m_b[3u * u] - y.x;
        m_out[3u * u + 1u] = m_b[3u * u + 1u] - y.y;
        m_out[3u * u + 2u] = m_b[3u * u + 2u] - y.z;
    }
}

// coarse_weight[E] = Σ_{e ∈ fine_edges(E)} fine_weight[e], one thread per
// coarse edge; fine-edge lists in ascending order. params.n = coarse edges.
@compute @workgroup_size(WG)
fn coarse_weight_update(
    @builtin(global_invocation_id) gid: vec3<u32>,
    @builtin(num_workgroups) nwg: vec3<u32>,
) {
    let stride = nwg.x * WG;
    for (var ce = gid.x; ce < params.n; ce += stride) {
        let begin = w_offsets[ce];
        let end = w_offsets[ce + 1u];
        var acc = Real(0);
        for (var j = begin; j < end; j++) {
            acc += w_fine_weight[w_fine_edges[j]];
        }
        w_coarse_weight[ce] = acc;
    }
}

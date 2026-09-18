//! GPU kernel microbenchmarks (WS-G): `apply_graph`, `dot3` and readback
//! latency at 1M and 10M free nodes on a synthetic grid graph, with achieved
//! GB/s against a bandwidth estimate from the bytes each kernel must move.
//!
//! ```text
//! cargo run -p theseus --release --features gpu --example gpu_kernels_bench
//! THESEUS_GPU_BENCH_NODES=1000000,10000000 THESEUS_GPU_BENCH_REPS=20 ...
//! THESEUS_GPU_BENCH_WG=64,128,256          # workgroup sizes to sweep
//! THESEUS_GPU_BENCH_PRECISION=f32,f64      # f64 only where SHADER_F64
//! ```
//!
//! Under a software adapter (lavapipe, `THESEUS_GPU_ALLOW_SOFTWARE=1`) the
//! numbers say nothing about GPU speed; the harness exists so WS-H can run
//! it unchanged on the validation machines. Output is one JSON line per
//! measurement after the table, for `benchmarks/results/`.
//!
//! Requires the `gpu` feature (`required-features` in `Cargo.toml`).

use std::time::Instant;

use theseus::backend::gpu::{GpuBackend, GpuContext};
use theseus::backend::{probe_gpu, Backend};
use theseus::graph::{CsrAdjacency, LevelGraph};
use theseus::linear_solver::Precision;

fn env_list<T: std::str::FromStr + Clone>(name: &str, default: &[T]) -> Vec<T> {
    match std::env::var(name) {
        Ok(v) => v.split(',').filter_map(|s| s.trim().parse().ok()).collect(),
        Err(_) => default.to_vec(),
    }
}

/// `w × w` grid with unit weights and boundary anchors, like the test
/// fixture but sized from the requested node count.
fn grid_level(nodes: usize) -> LevelGraph {
    let w = (nodes as f64).sqrt().ceil() as usize;
    let n = w * w;
    let id = |i: usize, j: usize| (i * w + j) as u32;
    let mut degree = vec![0u32; n + 1];
    let mut edges: Vec<(u32, u32)> = Vec::with_capacity(2 * n);
    for i in 0..w {
        for j in 0..w {
            if j + 1 < w {
                edges.push((id(i, j), id(i, j + 1)));
            }
            if i + 1 < w {
                edges.push((id(i, j), id(i + 1, j)));
            }
        }
    }
    for &(s, e) in &edges {
        degree[s as usize + 1] += 1;
        degree[e as usize + 1] += 1;
    }
    for i in 0..n {
        degree[i + 1] += degree[i];
    }
    let offsets = degree;
    let mut next = offsets.clone();
    let total = offsets[n] as usize;
    let mut adj_edges = vec![0u32; total];
    let mut sign = vec![0i8; total];
    let mut other = vec![0u32; total];
    for (e, &(s, t)) in edges.iter().enumerate() {
        let slot = next[s as usize] as usize;
        adj_edges[slot] = e as u32;
        sign[slot] = -1;
        other[slot] = t;
        next[s as usize] += 1;
        let slot = next[t as usize] as usize;
        adj_edges[slot] = e as u32;
        sign[slot] = 1;
        other[slot] = s;
        next[t as usize] += 1;
    }
    let mut anchor = vec![0.0; n];
    for i in 0..w {
        for j in 0..w {
            let boundary = usize::from(i == 0 || i + 1 == w) + usize::from(j == 0 || j + 1 == w);
            anchor[i * w + j] = boundary as f64;
        }
    }
    LevelGraph {
        n,
        adjacency: CsrAdjacency {
            offsets,
            edges: adj_edges,
            sign,
            other,
        },
        weight: vec![1.0; edges.len()],
        anchor,
        aggregate_of: Vec::new(),
    }
}

struct Timing {
    median_ms: f64,
    min_ms: f64,
}

/// Time `reps` runs of `f` (each ending in a device sync), after one warm-up.
fn time(reps: usize, mut f: impl FnMut()) -> Timing {
    f();
    let mut samples = Vec::with_capacity(reps);
    for _ in 0..reps {
        let t = Instant::now();
        f();
        samples.push(t.elapsed().as_secs_f64() * 1e3);
    }
    samples.sort_by(|a, b| a.partial_cmp(b).unwrap());
    Timing {
        median_ms: samples[samples.len() / 2],
        min_ms: samples[0],
    }
}

fn gbps(bytes: f64, ms: f64) -> f64 {
    bytes / (ms * 1e-3) / 1e9
}

fn main() {
    let probe = probe_gpu();
    println!("probe: {probe}");
    let ctx = match GpuContext::from_env() {
        Ok(ctx) => ctx,
        Err(e) => {
            println!("no usable adapter, nothing to benchmark: {e}");
            return;
        }
    };
    let adapter = ctx.report.clone();
    let nodes = env_list::<usize>("THESEUS_GPU_BENCH_NODES", &[1_000_000, 10_000_000]);
    let reps = env_list::<usize>("THESEUS_GPU_BENCH_REPS", &[10])[0].max(1);
    let wgs = env_list::<u32>(
        "THESEUS_GPU_BENCH_WG",
        &[GpuBackend::DEFAULT_WORKGROUP_SIZE],
    );
    let precisions: Vec<Precision> =
        env_list::<String>("THESEUS_GPU_BENCH_PRECISION", &["f32".to_string()])
            .iter()
            .filter_map(|p| match p.as_str() {
                "f32" => Some(Precision::F32),
                "f64" => Some(Precision::F64),
                _ => None,
            })
            .collect();
    if adapter.software {
        println!("software adapter: timings are not GPU performance data");
    }
    println!(
        "{:>10} {:>4} {:>3} {:>16} {:>10} {:>8} {:>8} {:>10}",
        "nodes", "prec", "wg", "kernel", "median ms", "min ms", "GB/s", "est. MB"
    );

    let mut results = Vec::new();
    let mut ctx = Some(ctx);
    for (wi, &wg) in wgs.iter().enumerate() {
        let ctx_now = match ctx.take() {
            Some(c) => c,
            None => GpuContext::from_env().expect("adapter disappeared"),
        };
        let gpu = match GpuBackend::with_workgroup_size(ctx_now, wg) {
            Ok(gpu) => gpu,
            Err(e) => {
                println!("workgroup size {wg}: {e}");
                continue;
            }
        };
        let _ = wi;
        for &target in &nodes {
            let level = grid_level(target);
            let n = level.n;
            let ne = level.num_edges();
            for &precision in &precisions {
                if !gpu.supports(precision) {
                    println!("{n:>10} {precision:?}: SHADER_F64 not available, skipped");
                    continue;
                }
                let bytes = precision.size_of() as f64;
                let vec_bytes = 3.0 * n as f64 * bytes;
                // apply_graph traffic: x read (3n, once if cached; count the
                // gather as 2·ne extra reads of 3 scalars), y write, weights
                // gathered 2·ne, anchors n, CSR offsets (n+1) + edges 2·ne +
                // other 2·ne as u32.
                let apply_bytes = vec_bytes * 2.0
                    + 2.0 * ne as f64 * (3.0 * bytes + bytes + 8.0)
                    + n as f64 * bytes
                    + (n as f64 + 1.0) * 4.0;
                let dev = match gpu.try_upload_level(&level, precision) {
                    Ok(dev) => dev,
                    Err(e) => {
                        println!("{n:>10} {precision:?}: level upload failed: {e}");
                        continue;
                    }
                };
                let x_host: Vec<f64> = (0..n * 3).map(|i| ((i % 97) as f64) / 97.0 - 0.5).collect();
                let (mut x, mut y) = match (
                    gpu.try_alloc(n * 3, precision),
                    gpu.try_alloc(n * 3, precision),
                ) {
                    (Ok(x), Ok(y)) => (x, y),
                    (Err(e), _) | (_, Err(e)) => {
                        println!("{n:>10} {precision:?}: vector allocation failed: {e}");
                        continue;
                    }
                };
                gpu.pool().reserve_staging(vec_bytes as u64);
                gpu.upload(&x_host, &mut x);
                gpu.sync();

                let t = time(reps, || {
                    gpu.apply_graph(&dev, &x, &mut y);
                    gpu.sync();
                });
                print_row(
                    n,
                    precision,
                    wg,
                    "apply_graph",
                    &t,
                    gbps(apply_bytes, t.median_ms),
                    apply_bytes,
                );
                results.push(json_row(
                    &adapter,
                    n,
                    ne,
                    precision,
                    wg,
                    "apply_graph",
                    &t,
                    apply_bytes,
                ));

                let t = time(reps, || {
                    let _ = gpu.dot3(&x, &y);
                });
                print_row(
                    n,
                    precision,
                    wg,
                    "dot3",
                    &t,
                    gbps(2.0 * vec_bytes, t.median_ms),
                    2.0 * vec_bytes,
                );
                results.push(json_row(
                    &adapter,
                    n,
                    ne,
                    precision,
                    wg,
                    "dot3",
                    &t,
                    2.0 * vec_bytes,
                ));

                // Readback latency: a 12-byte norm-sized readback of a tiny
                // vector (the per-iteration cost of leaving the device), and
                // the full vector download.
                let tiny = gpu.alloc(3, precision);
                let t = time(reps, || {
                    let _ = gpu.norm3(&tiny);
                });
                print_row(n, precision, wg, "readback 3 scalars", &t, 0.0, 3.0 * bytes);
                results.push(json_row(
                    &adapter,
                    n,
                    ne,
                    precision,
                    wg,
                    "readback_scalar",
                    &t,
                    3.0 * bytes,
                ));

                let mut host = vec![0.0; n * 3];
                let t = time(reps, || {
                    gpu.download(&y, &mut host);
                });
                print_row(
                    n,
                    precision,
                    wg,
                    "download vector",
                    &t,
                    gbps(vec_bytes, t.median_ms),
                    vec_bytes,
                );
                results.push(json_row(
                    &adapter,
                    n,
                    ne,
                    precision,
                    wg,
                    "download_vector",
                    &t,
                    vec_bytes,
                ));

                let t = time(reps, || {
                    gpu.upload(&x_host, &mut x);
                    gpu.sync();
                });
                print_row(
                    n,
                    precision,
                    wg,
                    "upload vector",
                    &t,
                    gbps(vec_bytes, t.median_ms),
                    vec_bytes,
                );
                results.push(json_row(
                    &adapter,
                    n,
                    ne,
                    precision,
                    wg,
                    "upload_vector",
                    &t,
                    vec_bytes,
                ));
            }
            println!(
                "{:>10} device bytes after this size: {:.1} MB",
                n,
                gpu.device_bytes() as f64 / 1e6
            );
        }
        let stats = gpu.stats();
        println!(
            "wg {wg}: {} dispatches, {} submits, {} readbacks",
            stats.dispatches, stats.submits, stats.readbacks
        );
        drop(gpu);
    }
    println!();
    println!(
        "bandwidth estimate = bytes each kernel must move at least once; the 'GB/s' column is that\n\
         estimate over the median wall time including submit + sync. Real-GPU numbers are pending\n\
         WS-H's validation machines; software-adapter numbers are not performance data."
    );
    for line in results {
        println!("{line}");
    }
}

fn print_row(
    n: usize,
    precision: Precision,
    wg: u32,
    kernel: &str,
    t: &Timing,
    gbps: f64,
    bytes: f64,
) {
    let gb = if gbps > 0.0 {
        format!("{gbps:8.2}")
    } else {
        format!("{:>8}", "-")
    };
    println!(
        "{n:>10} {:>4} {wg:>3} {kernel:>16} {:>10.3} {:>8.3} {gb} {:>10.1}",
        match precision {
            Precision::F32 => "f32",
            Precision::F64 => "f64",
        },
        t.median_ms,
        t.min_ms,
        bytes / 1e6
    );
}

#[allow(clippy::too_many_arguments)]
fn json_row(
    adapter: &theseus::backend::GpuAdapterReport,
    n: usize,
    ne: usize,
    precision: Precision,
    wg: u32,
    kernel: &str,
    t: &Timing,
    bytes: f64,
) -> String {
    format!(
        "{{\"bench\":\"gpu_kernels\",\"adapter\":{:?},\"backend\":{:?},\"device_type\":{:?},\"software\":{},\
         \"nodes\":{n},\"edges\":{ne},\"precision\":{:?},\"workgroup_size\":{wg},\"kernel\":{kernel:?},\
         \"median_ms\":{:.4},\"min_ms\":{:.4},\"bytes_estimate\":{:.0},\"gbps\":{:.3}}}",
        adapter.name,
        adapter.backend,
        adapter.device_type,
        adapter.software,
        match precision {
            Precision::F32 => "f32",
            Precision::F64 => "f64",
        },
        t.median_ms,
        t.min_ms,
        bytes,
        if t.median_ms > 0.0 { gbps(bytes, t.median_ms) } else { 0.0 }
    )
}

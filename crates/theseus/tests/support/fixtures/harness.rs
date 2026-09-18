//! Helpers shared by the benchmark harnesses (`tests/bench_scale.rs`,
//! `examples/profile_phases.rs`): fixture selection by name at matched edge
//! counts, environment switches, timing statistics, machine/git identification
//! and JSON-lines output (§5.1 of `ITERATIVE_SOLVER_PROGRAM.md`).

use super::grid;
use super::{anisotropic, dome, few_supports, irregular};
use serde_json::{json, Map, Value};
use std::io::Write;
use theseus::types::Problem;

/// Fixture families selectable through `THESEUS_FIXTURE`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FixtureKind {
    Grid,
    Irregular,
    Dome,
    FewSupports,
    Anisotropic,
}

impl FixtureKind {
    pub const ALL: [FixtureKind; 5] = [
        FixtureKind::Grid,
        FixtureKind::Irregular,
        FixtureKind::Dome,
        FixtureKind::FewSupports,
        FixtureKind::Anisotropic,
    ];

    pub fn parse(name: &str) -> Option<Self> {
        match name.trim().to_ascii_lowercase().as_str() {
            "grid" => Some(Self::Grid),
            "irregular" | "irregular-mesh" | "mesh" => Some(Self::Irregular),
            "dome" | "cable-dome" | "cable_dome" => Some(Self::Dome),
            "few-supports" | "few_supports" | "fewsupports" => Some(Self::FewSupports),
            "anisotropic" | "anisotropic-q" | "anisotropic_q" => Some(Self::Anisotropic),
            _ => None,
        }
    }

    /// Name used in tables and JSON.
    pub fn name(self) -> &'static str {
        match self {
            Self::Grid => "grid",
            Self::Irregular => "irregular",
            Self::Dome => "dome",
            Self::FewSupports => "few-supports",
            Self::Anisotropic => "anisotropic",
        }
    }

    /// Reads `THESEUS_FIXTURE` (default `grid`); panics on an unknown name so a
    /// typo does not silently benchmark the wrong fixture.
    pub fn from_env() -> Self {
        match std::env::var("THESEUS_FIXTURE") {
            Ok(v) => Self::parse(&v).unwrap_or_else(|| {
                panic!("THESEUS_FIXTURE={v:?}: expected one of grid, irregular, dome, few-supports, anisotropic")
            }),
            Err(_) => Self::Grid,
        }
    }
}

/// Default seed of the irregular mesh in the harnesses.
pub const IRREGULAR_SEED: u64 = 0x5EED;
/// Default `q*` ratio of the anisotropic fixture in the harnesses
/// (override with `THESEUS_ANISOTROPY`).
pub const ANISOTROPY_RATIO: f64 = 100.0;

/// A built fixture with the parameters that describe it.
pub struct Fixture {
    pub kind: FixtureKind,
    pub problem: Problem,
    /// Generator parameters (grid side, seed, rings, spokes, ratio...).
    pub parameters: Map<String, Value>,
}

/// Builds the fixture of `kind` at the edge count of the `grid_side × grid_side`
/// grid (the grid itself for `Grid`, `FewSupports` and `Anisotropic`; matched
/// sizes for `Irregular` and `Dome`).
pub fn build_fixture(kind: FixtureKind, grid_side: usize) -> Fixture {
    let target_edges = grid::grid_edges(grid_side);
    let mut parameters = Map::new();
    parameters.insert("grid_side".into(), json!(grid_side));
    parameters.insert("target_edges".into(), json!(target_edges));
    let problem = match kind {
        FixtureKind::Grid => grid::make_recoverable_grid_problem(grid_side),
        FixtureKind::Irregular => {
            let side = irregular::side_for_edges(target_edges);
            parameters.insert("side".into(), json!(side));
            parameters.insert("seed".into(), json!(IRREGULAR_SEED));
            parameters.insert("jitter".into(), json!(irregular::JITTER));
            irregular::make_irregular_mesh_problem(side, IRREGULAR_SEED)
        }
        FixtureKind::Dome => {
            let (rings, spokes) = dome::size_for_edges(target_edges);
            parameters.insert("rings".into(), json!(rings));
            parameters.insert("spokes".into(), json!(spokes));
            dome::make_cable_dome_problem(rings, spokes)
        }
        FixtureKind::FewSupports => {
            parameters.insert("supports".into(), json!(2));
            few_supports::make_few_supports_problem(grid_side)
        }
        FixtureKind::Anisotropic => {
            let ratio = std::env::var("THESEUS_ANISOTROPY")
                .ok()
                .and_then(|v| v.parse().ok())
                .unwrap_or(ANISOTROPY_RATIO);
            parameters.insert("ratio".into(), json!(ratio));
            anisotropic::make_anisotropic_q_problem(grid_side, ratio)
        }
    };
    Fixture {
        kind,
        problem,
        parameters,
    }
}

/// Comma-separated list of integers from the environment, or the default.
pub fn env_list(name: &str, default: &[usize]) -> Vec<usize> {
    std::env::var(name)
        .ok()
        .map(|v| {
            v.split(',')
                .filter(|s| !s.trim().is_empty())
                .map(|s| {
                    s.trim()
                        .parse()
                        .unwrap_or_else(|_| panic!("{name}: {s:?} is not an integer"))
                })
                .collect()
        })
        .unwrap_or_else(|| default.to_vec())
}

/// Single integer from the environment, or the default.
pub fn env_usize(name: &str, default: usize) -> usize {
    env_list(name, &[default])[0]
}

/// Requested linear-solver backend (`THESEUS_LINEAR_SOLVER`, default
/// `direct`), lower-cased. Only `direct` is implemented; callers print
/// "not yet available" and skip for anything else.
pub fn linear_solver_from_env() -> String {
    std::env::var("THESEUS_LINEAR_SOLVER")
        .ok()
        .map(|v| v.trim().to_ascii_lowercase())
        .filter(|v| !v.is_empty())
        .unwrap_or_else(|| "direct".into())
}

pub fn is_direct(backend: &str) -> bool {
    matches!(backend, "direct" | "0")
}

/// Thread count: `RAYON_NUM_THREADS` when set, otherwise rayon's detected pool size.
pub fn thread_count() -> usize {
    std::env::var("RAYON_NUM_THREADS")
        .ok()
        .and_then(|v| v.parse().ok())
        .filter(|&n: &usize| n > 0)
        .unwrap_or_else(rayon::current_num_threads)
}

/// Peak resident set size in bytes (`VmHWM` from `/proc/self/status`).
pub fn peak_rss_bytes() -> Option<u64> {
    let status = std::fs::read_to_string("/proc/self/status").ok()?;
    let line = status.lines().find(|l| l.starts_with("VmHWM:"))?;
    let kb: u64 = line.split_whitespace().nth(1)?.parse().ok()?;
    Some(kb * 1024)
}

/// Order statistics of a set of timings.
#[derive(Debug, Clone)]
pub struct Timing {
    pub samples: Vec<f64>,
    pub median: f64,
    pub q25: f64,
    pub q75: f64,
    pub min: f64,
    pub max: f64,
}

impl Timing {
    pub fn from_samples(mut samples: Vec<f64>) -> Self {
        assert!(!samples.is_empty());
        samples.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let q = |p: f64| -> f64 {
            let pos = p * (samples.len() - 1) as f64;
            let lo = pos.floor() as usize;
            let hi = pos.ceil() as usize;
            samples[lo] + (samples[hi] - samples[lo]) * (pos - lo as f64)
        };
        Self {
            median: q(0.5),
            q25: q(0.25),
            q75: q(0.75),
            min: samples[0],
            max: samples[samples.len() - 1],
            samples,
        }
    }

    pub fn iqr(&self) -> f64 {
        self.q75 - self.q25
    }

    pub fn to_json(&self) -> Value {
        json!({
            "median": self.median,
            "iqr": self.iqr(),
            "q25": self.q25,
            "q75": self.q75,
            "min": self.min,
            "max": self.max,
            "n": self.samples.len(),
            "samples": self.samples,
        })
    }
}

/// `git rev-parse HEAD` at run time (with `-dirty` when the tree has
/// uncommitted changes), or `"unknown"`.
pub fn git_sha() -> String {
    let run = |args: &[&str]| {
        std::process::Command::new("git")
            .args(args)
            .output()
            .ok()
            .filter(|o| o.status.success())
            .map(|o| String::from_utf8_lossy(&o.stdout).trim().to_string())
    };
    match run(&["rev-parse", "HEAD"]) {
        Some(sha) if !sha.is_empty() => {
            let dirty = run(&["status", "--porcelain", "--untracked-files=no"])
                .map(|s| !s.is_empty())
                .unwrap_or(false);
            if dirty {
                format!("{sha}-dirty")
            } else {
                sha
            }
        }
        _ => "unknown".into(),
    }
}

/// Lower-case, dash-separated identifier; trademark noise (`(R)`, `(TM)`,
/// `CPU`, `Processor`, the `@ x.xxGHz` suffix) is dropped so the id is short.
/// Keep in sync with `machine_slug` in `scripts/bench_sweep.py`.
fn slug(s: &str) -> String {
    let s = s.split('@').next().unwrap_or(s);
    let mut cleaned = s.to_ascii_lowercase();
    for noise in ["(r)", "(tm)", "(c)", " cpu", " processor", "core "] {
        cleaned = cleaned.replace(noise, " ");
    }
    let mut out = String::new();
    let mut dash = false;
    for c in cleaned.chars() {
        if c.is_ascii_alphanumeric() {
            out.push(c.to_ascii_lowercase());
            dash = false;
        } else if !dash && !out.is_empty() {
            out.push('-');
            dash = true;
        }
    }
    out.trim_end_matches('-').to_string()
}

/// CPU model name (Linux `/proc/cpuinfo`, macOS `sysctl`), or the architecture.
fn cpu_model() -> String {
    if let Ok(info) = std::fs::read_to_string("/proc/cpuinfo") {
        if let Some(line) = info.lines().find(|l| l.starts_with("model name")) {
            if let Some((_, v)) = line.split_once(':') {
                return v.trim().to_string();
            }
        }
    }
    if cfg!(target_os = "macos") {
        if let Ok(o) = std::process::Command::new("sysctl")
            .args(["-n", "machdep.cpu.brand_string"])
            .output()
        {
            let s = String::from_utf8_lossy(&o.stdout).trim().to_string();
            if !s.is_empty() {
                return s;
            }
        }
    }
    std::env::consts::ARCH.to_string()
}

/// Total RAM in whole GiB where discoverable.
fn ram_gib() -> Option<u64> {
    if let Ok(info) = std::fs::read_to_string("/proc/meminfo") {
        let line = info.lines().find(|l| l.starts_with("MemTotal:"))?;
        let kb: u64 = line.split_whitespace().nth(1)?.parse().ok()?;
        return Some(((kb as f64) / (1024.0 * 1024.0)).round() as u64);
    }
    if cfg!(target_os = "macos") {
        let o = std::process::Command::new("sysctl")
            .args(["-n", "hw.memsize"])
            .output()
            .ok()?;
        let bytes: u64 = String::from_utf8_lossy(&o.stdout).trim().parse().ok()?;
        return Some(((bytes as f64) / (1024.0 * 1024.0 * 1024.0)).round() as u64);
    }
    None
}

/// Machine identifier: `THESEUS_MACHINE_ID` when set (the sweep script passes
/// the id it derived, including GPU adapters), else
/// `<os>-<cpu model>-<logical cores>c-<ram>gb` built from what the process can
/// see (no GPU probe in the Rust harness).
pub fn machine_id() -> String {
    if let Ok(id) = std::env::var("THESEUS_MACHINE_ID") {
        if !id.trim().is_empty() {
            return id.trim().to_string();
        }
    }
    let cores = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(0);
    let mut parts = vec![std::env::consts::OS.to_string(), slug(&cpu_model())];
    parts.push(format!("{cores}c"));
    if let Some(gb) = ram_gib() {
        parts.push(format!("{gb}gb"));
    }
    parts.join("-")
}

/// Appends `record` as one line to the file named by `THESEUS_BENCH_JSON`
/// (created when missing). Returns `true` when a line was written.
pub fn append_json_line(record: &Value) -> bool {
    let Ok(path) = std::env::var("THESEUS_BENCH_JSON") else {
        return false;
    };
    if path.trim().is_empty() {
        return false;
    }
    if let Some(parent) = std::path::Path::new(&path).parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent).expect("create THESEUS_BENCH_JSON directory");
        }
    }
    let mut file = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(&path)
        .unwrap_or_else(|e| panic!("open {path}: {e}"));
    let line = serde_json::to_string(record).expect("serialise benchmark record");
    writeln!(file, "{line}").expect("write benchmark record");
    true
}

/// ISO-8601 UTC timestamp (seconds) without pulling in a date crate.
pub fn timestamp_utc() -> String {
    let secs = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_secs() as i64)
        .unwrap_or(0);
    let days = secs.div_euclid(86_400);
    let rem = secs.rem_euclid(86_400);
    let (h, m, s) = (rem / 3600, (rem % 3600) / 60, rem % 60);
    // Civil-from-days (Howard Hinnant's algorithm).
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z.rem_euclid(146_097);
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let y = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let mth = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = if mth <= 2 { y + 1 } else { y };
    format!("{y:04}-{mth:02}-{d:02}T{h:02}:{m:02}:{s:02}Z")
}

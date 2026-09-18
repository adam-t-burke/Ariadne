"""Shared helpers for bench_sweep.py and crossover_fit.py (run via uv, see
those scripts for usage). Loading of runs.jsonl produced by
`crates/theseus/tests/bench_scale.rs`, machine identification and small
formatting utilities. No third-party imports so it can be used from plain
Python as well."""

from __future__ import annotations

import json
import os
import platform
import re
import shutil
import subprocess
from pathlib import Path

FIXTURE_ORDER = ["grid", "irregular", "dome", "few-supports", "anisotropic"]


# ── runs.jsonl ─────────────────────────────────────────────────────────────


def load_runs(path: Path, harness: str = "bench_scale") -> list[dict]:
    """Reads a JSON-lines file; keeps successful `bench_scale` records (those
    without a `status` other than `ok`) and drops malformed lines."""
    runs = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            try:
                rec = json.loads(line)
            except json.JSONDecodeError:
                continue
            if rec.get("harness", harness) != harness:
                continue
            if rec.get("status", "ok") != "ok":
                continue
            if "eval_ms" not in rec or "edges" not in rec:
                continue
            runs.append(rec)
    return runs


def load_failures(path: Path) -> list[dict]:
    """Records appended by the sweep for cells that did not complete."""
    out = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            try:
                rec = json.loads(line)
            except json.JSONDecodeError:
                continue
            if rec.get("status", "ok") != "ok":
                out.append(rec)
    return out


def eval_median(rec: dict) -> float:
    return float(rec["eval_ms"]["median"])


def eval_iqr(rec: dict) -> float:
    return float(rec["eval_ms"].get("iqr", 0.0))


def eval_samples(rec: dict) -> list[float]:
    s = rec["eval_ms"].get("samples")
    return [float(v) for v in s] if s else [eval_median(rec)]


CELL_KEYS = ("backend", "fixture", "grid_side", "threads")


def is_flagged(rec: dict, ratio: float = 0.10) -> bool:
    """IQR above 10% of the median (§5.2 flag)."""
    m = eval_median(rec)
    return m > 0 and eval_iqr(rec) > ratio * m


def select_runs(runs: list[dict]) -> tuple[list[dict], int]:
    """One record per cell. When a cell was measured more than once (a re-run
    after external load disturbed it), keep the run that is not flagged and
    has the smallest evaluation median — under contention the minimum is the
    least-disturbed estimate — and fall back to the least-flagged run.
    Returns the selection and the number of superseded records."""
    cells: dict[tuple, list[dict]] = {}
    for r in runs:
        cells.setdefault(tuple(r.get(k) for k in CELL_KEYS), []).append(r)
    selected = []
    for rs in cells.values():
        rs.sort(key=lambda r: (is_flagged(r), eval_median(r)))
        chosen = dict(rs[0])
        chosen["runs_in_cell"] = len(rs)
        selected.append(chosen)
    return selected, len(runs) - len(selected)


def fixture_sort_key(name: str) -> tuple[int, str]:
    return (FIXTURE_ORDER.index(name) if name in FIXTURE_ORDER else len(FIXTURE_ORDER), name)


def group_by(runs: list[dict], *keys: str) -> dict[tuple, list[dict]]:
    groups: dict[tuple, list[dict]] = {}
    for r in runs:
        groups.setdefault(tuple(r.get(k) for k in keys), []).append(r)
    return groups


# ── machine identification ─────────────────────────────────────────────────

_NOISE = ["(r)", "(tm)", "(c)", " cpu", " processor", "core "]


def slug(text: str) -> str:
    """Mirror of `slug` in tests/support/fixtures/harness.rs."""
    text = text.split("@")[0].lower()
    for noise in _NOISE:
        text = text.replace(noise, " ")
    text = re.sub(r"[^a-z0-9]+", "-", text)
    return text.strip("-")


def _run(cmd: list[str]) -> str | None:
    try:
        out = subprocess.run(cmd, capture_output=True, text=True, timeout=20)
    except (OSError, subprocess.TimeoutExpired):
        return None
    if out.returncode != 0:
        return None
    return out.stdout


def cpu_model() -> str:
    system = platform.system()
    if system == "Linux":
        try:
            with open("/proc/cpuinfo", encoding="utf-8") as fh:
                for line in fh:
                    if line.startswith("model name"):
                        return line.split(":", 1)[1].strip()
        except OSError:
            pass
    elif system == "Darwin":
        out = _run(["sysctl", "-n", "machdep.cpu.brand_string"])
        if out and out.strip():
            return out.strip()
    elif system == "Windows":
        out = _run(["powershell", "-NoProfile", "-Command",
                    "(Get-CimInstance Win32_Processor).Name"])
        if out and out.strip():
            return out.strip().splitlines()[0]
    return platform.processor() or platform.machine()


def physical_cores() -> int | None:
    system = platform.system()
    if system == "Linux":
        out = _run(["lscpu", "-p=CORE,SOCKET"])
        if out:
            cores = {tuple(l.split(",")) for l in out.splitlines() if l and not l.startswith("#")}
            if cores:
                return len(cores)
    elif system == "Darwin":
        out = _run(["sysctl", "-n", "hw.physicalcpu"])
        if out and out.strip().isdigit():
            return int(out.strip())
    elif system == "Windows":
        out = _run(["powershell", "-NoProfile", "-Command",
                    "(Get-CimInstance Win32_Processor | Measure-Object NumberOfCores -Sum).Sum"])
        if out and out.strip().isdigit():
            return int(out.strip())
    return None


def ram_bytes() -> int | None:
    system = platform.system()
    if system == "Linux":
        try:
            with open("/proc/meminfo", encoding="utf-8") as fh:
                for line in fh:
                    if line.startswith("MemTotal:"):
                        return int(line.split()[1]) * 1024
        except OSError:
            pass
    elif system == "Darwin":
        out = _run(["sysctl", "-n", "hw.memsize"])
        if out and out.strip().isdigit():
            return int(out.strip())
    elif system == "Windows":
        out = _run(["powershell", "-NoProfile", "-Command",
                    "(Get-CimInstance Win32_ComputerSystem).TotalPhysicalMemory"])
        if out and out.strip().isdigit():
            return int(out.strip())
    try:
        return os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")
    except (ValueError, OSError, AttributeError):
        return None


def gpu_adapters() -> list[str]:
    """Best-effort GPU adapter names (lspci / system_profiler / CIM)."""
    system = platform.system()
    names: list[str] = []
    if system == "Linux" and shutil.which("lspci"):
        out = _run(["lspci"])
        for line in (out or "").splitlines():
            if re.search(r"VGA|3D controller|Display controller", line):
                names.append(line.split(":", 2)[-1].strip())
        if not names:
            for path in sorted(Path("/sys/class/drm").glob("card[0-9]*/device/uevent")):
                try:
                    txt = path.read_text()
                except OSError:
                    continue
                m = re.search(r"PCI_ID=(\S+)", txt)
                if m:
                    names.append(f"pci {m.group(1)}")
    elif system == "Darwin":
        out = _run(["system_profiler", "SPDisplaysDataType", "-json"])
        if out:
            try:
                data = json.loads(out)
                for item in data.get("SPDisplaysDataType", []):
                    if item.get("sppci_model"):
                        names.append(item["sppci_model"])
            except json.JSONDecodeError:
                pass
    elif system == "Windows":
        out = _run(["powershell", "-NoProfile", "-Command",
                    "(Get-CimInstance Win32_VideoController).Name"])
        for line in (out or "").splitlines():
            if line.strip():
                names.append(line.strip())
    return names


def rust_version() -> str | None:
    """`rustc --version` under the caller's toolchain (RUSTUP_TOOLCHAIN is inherited)."""
    out = _run(["rustc", "--version"])
    return out.strip() if out else None


def git_sha(root: Path) -> str:
    try:
        out = subprocess.run(["git", "rev-parse", "HEAD"], cwd=root, capture_output=True,
                             text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"
    dirty = subprocess.run(["git", "status", "--porcelain", "--untracked-files=no"], cwd=root,
                           capture_output=True, text=True).stdout.strip()
    return f"{out}-dirty" if dirty else out


def os_name() -> str:
    return {"Linux": "linux", "Darwin": "macos", "Windows": "windows"}.get(
        platform.system(), platform.system().lower())


def machine_info(root: Path) -> dict:
    ram = ram_bytes()
    gpus = gpu_adapters()
    info = {
        "os": os_name(),
        "os_release": platform.platform(),
        "cpu_model": cpu_model(),
        "logical_cores": os.cpu_count(),
        "physical_cores": physical_cores(),
        "ram_bytes": ram,
        "ram_gib": round(ram / 2**30) if ram else None,
        "gpu_adapters": gpus,
        "rust_version": rust_version(),
        "python_version": platform.python_version(),
        "hostname": platform.node(),
        "git_sha": git_sha(root),
    }
    info["machine_id"] = machine_id_from_info(info)
    return info


def machine_id_from_info(info: dict) -> str:
    """`<os>-<cpu>-<logical cores>c-<ram>gb[-<gpu>]`; the Rust harness builds
    the same prefix without the GPU part when THESEUS_MACHINE_ID is unset."""
    parts = [info["os"], slug(info["cpu_model"] or "cpu"), f"{info.get('logical_cores') or 0}c"]
    if info.get("ram_gib"):
        parts.append(f"{info['ram_gib']}gb")
    real_gpus = [g for g in info.get("gpu_adapters", [])
                 if not re.search(r"llvmpipe|swiftshader|basic render|virtio|cirrus|vmware|bochs", g, re.I)]
    if real_gpus:
        parts.append(slug(real_gpus[0])[:32].strip("-"))
    return "-".join(p for p in parts if p)


# ── formatting ─────────────────────────────────────────────────────────────


def fmt_ms(v: float) -> str:
    if v >= 10000:
        return f"{v:,.0f}"
    if v >= 100:
        return f"{v:.0f}"
    if v >= 10:
        return f"{v:.1f}"
    return f"{v:.2f}"


def fmt_int(v) -> str:
    return f"{int(v):,}" if v is not None else "n/a"


def fmt_pm(median: float, iqr: float, flag_ratio: float = 0.10) -> str:
    """`median ± iqr`, suffixed with ` !` when the IQR exceeds 10% of the median."""
    flag = " !" if median > 0 and iqr > flag_ratio * median else ""
    return f"{fmt_ms(median)} ± {fmt_ms(iqr)}{flag}"

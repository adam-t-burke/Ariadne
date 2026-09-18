# Benchmark results and reports

Raw measurements and rendered reports of the Theseus scaling benchmarks
(§5 of `crates/theseus/ITERATIVE_SOLVER_PROGRAM.md`). Everything here is
produced by the tooling in `scripts/`; do not edit the JSON by hand.

```
benchmarks/
  results/<machine-id>/<YYYYMMDD>-<sha7>/
    runs.jsonl            one JSON object per benchmark cell (schema below)
    machine.json          OS, CPU model, cores, RAM, GPU adapters, Rust version, git sha
    config.json           sweep configuration (sizes, fixtures, threads, reps, iterations, solvers)
    harness-output.txt    the harness tables, per cell
  reports/
    <machine-id>-<YYYYMMDD>.md          tables (median ± IQR), plots, direct cost-model fit
    <machine-id>-<YYYYMMDD>-eval.svg    evaluation time vs edges, log-log
    <machine-id>-<YYYYMMDD>-total.svg   full-solve time vs edges, log-log
    crossover-<machine-id>.md / .json / .svg   cost models, residuals, crossover (§5.3)
```

## Running a sweep

```sh
export RUSTUP_TOOLCHAIN=1.89.0
uv run --project scripts scripts/bench_sweep.py \
    --sizes 72,160,224,320,448 --fixtures grid,irregular,dome \
    --threads 1,4 --reps 5 --iters 10 --solvers direct
uv run --project scripts scripts/crossover_fit.py \
    benchmarks/results/<machine-id>/<YYYYMMDD>-<sha7>/runs.jsonl
```

`bench_sweep.py` builds `tests/bench_scale.rs` once (`--locked`, release) and
runs one process per cell with `THESEUS_FIXTURE`, `THESEUS_SCALE_GRIDS`,
`THESEUS_SCALE_ITERS`, `THESEUS_BENCH_REPS`, `THESEUS_LINEAR_SOLVER`,
`THESEUS_BENCH_JSON`, `THESEUS_MACHINE_ID` and `RAYON_NUM_THREADS` set. A JSON
config file (`--config`) or CLI flags select the grid; `--report-only <dir>`
re-renders a report from existing results. `crossover_fit.py` fits the §5.3
cost models (direct: `a·n^1.5 + b·n·log n + c`; iterative:
`(d + e·iters)·n + f`) on log-scaled data, reports relative residuals and,
when both backends are present, the crossover edge count with a bootstrap
90% interval.

The same environment switches work on the harness directly:

```sh
THESEUS_FIXTURE=irregular THESEUS_SCALE_GRIDS=160,224 THESEUS_BENCH_REPS=7 \
THESEUS_BENCH_JSON=/tmp/runs.jsonl RAYON_NUM_THREADS=1 \
cargo test --release -p theseus --test bench_scale -- --ignored --nocapture
THESEUS_FIXTURE=dome cargo run --release -p theseus --example profile_phases -- 72 160
```

## `runs.jsonl` record (schema 1, `harness: "bench_scale"`)

| field | meaning |
|---|---|
| `fixture` | `grid`, `irregular`, `dome`, `few-supports`, `anisotropic` |
| `grid_side` | grid side the cell was requested at; non-grid fixtures are generated at that grid's edge count |
| `edges`, `free_nodes`, `nodes`, `fixed_nodes` | realised sizes |
| `backend` | `THESEUS_LINEAR_SOLVER` as given (`direct` is the only implemented value) |
| `adapter` | GPU adapter name (`null` until a GPU backend exists) |
| `threads` | `RAYON_NUM_THREADS`, or rayon's detected pool size |
| `parameters` | generator parameters (`side`, `seed`, `rings`, `spokes`, `ratio`, …), `iterations_budget`, `eval_reps`, `tolerances`, `bounds`, `q_start`, `q_parameterization` |
| `build_ms` | fixture generation time (not part of any benchmark number) |
| `setup_ms` | `FdmCache::new` + symbolic analysis + first factorisation |
| `eval_ms` | `{median, iqr, q25, q75, min, max, n, samples}` of the timed fused evaluations (one warm-up discarded) |
| `evaluations`, `iterations` | fused evaluations and accepted L-BFGS-B iterations of the full solve |
| `linear_solver_iterations` | `null` for `direct`; iterative backends will record forward/adjoint mean/max |
| `total_ms`, `ms_per_iteration`, `non_eval_ms` | full solve wall time; per accepted iteration; total minus `evaluations × eval median` |
| `peak_rss_bytes` | `VmHWM` of the process (Linux), else `null` |
| `device_bytes` | `null` until a GPU backend exists |
| `final_loss`, `termination` | last loss of the trace and the optimizer's termination string |
| `git_sha` | `git rev-parse HEAD` at run time, `-dirty` suffix when the tree has changes, or `unknown` |
| `machine_id` | `<os>-<cpu>-<cores>c-<ram>gb[-<gpu>]` (from `bench_sweep.py`, or derived by the harness) |
| `timestamp_utc`, `schema`, `harness` | provenance |

Cells that time out or fail are appended by the sweep with a `status` field
(`timeout`, `failed`, `unavailable`) instead of the measurement fields.
`examples/profile_phases.rs` appends records with `harness: "profile_phases"`
and a `phases_ms` object (per-phase medians) when `THESEUS_BENCH_JSON` is set.

# External case studies

Structures from other open-source inverse form-finding implementations, run
through the Theseus warm-start pipeline (`docs/warm_start.md`, Section 6.1).

| file | what it does |
|---|---|
| `export_jax_fdm.py` | runs the JAX FDM bundled examples, exports each optimised result as a case, and runs JAX FDM's own L-BFGS-B on the same position fit as the external baseline (`uv run python export_jax_fdm.py`, env from `pyproject.toml` / `uv.lock`) |
| `export_compas_cem.py` | runs the compas_cem examples and the CAD 2023 paper structures, exports the CEM equilibria as cases (separate venv `.venv-cem`; setup and command in the file header) |
| `cases/*.json` | one structure per file: nodes at the reference equilibrium, edges, supports, nodal loads, `q_ref`, box, target(s), and the external tool's results |
| `results/*_summary.md` | the exporters' tables (tool versions, upstream commits, original goals and timings, caveats) |
| `results/rust_external_run.txt` | full output of `warm_start_bench external cases` |
| `results/rust_<case>.json` | per-case rows of that run |
| `upstream_data/` | JAX FDM example meshes that are not in the pip wheel |

Rerun the Rust side from the repository root:

```sh
cargo build --release -p theseus --example warm_start_bench
BENCH_RESULTS=bench/external/results \
  ./target/release/examples/warm_start_bench external bench/external/cases
```

`BENCH_ITERS` caps L-BFGS-B (default 1000), `BENCH_METHODS` selects warm-start
methods. Every case is checked at load time: Theseus' forward solve at `q_ref`
must reproduce the exported geometry.

Case schema: `nodes`, `target`, `loads` are per node (all nodes, 0-based);
`edges`, `q_ref`, `signs`, `bounds.lo/hi` are per edge (`null` bound = that
side unbounded); `fixed` lists supports; `target_original` (optional, per
node, `null` where absent) is the example's own design target; `external`
holds the source tool's runs. Sign convention: `q > 0` tension.

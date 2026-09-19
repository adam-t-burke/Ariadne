# Figures for the warm-start write-up

`render.py` turns the output of `warm_start_bench figures|dense|scale` into
the PNGs under `docs/figures/`, which `docs/warm_start.md` and
`docs/presentation_outline.md` embed.

```bash
cargo build --release -p theseus --example warm_start_bench
B=target/release/examples/warm_start_bench
$B figures bench/figures/data
$B dense  > bench/figures/data/dense.txt
BENCH_METHODS=uniform,s1,gram_sparse,gram_dense,length_ratio,frozen,pipeline,legacy \
$B scale  > bench/figures/data/scale.txt
cd bench/figures && uv run render.py
```

`data/*.json` (geometry, force densities and L-BFGS-B traces of the showcase
cases, ~6 MB) is regenerated in about 20 s and not tracked; `dense.txt` and
`scale.txt` are tracked because the scale run takes ~15 minutes.

The GN-versus-L-BFGS-B residual experiment (`warm_start_bench tradeoff`)
writes `data/tradeoff.json`. Render with `uv run render_tradeoff.py`; the
write-up is `docs/gn_vs_lbfgs.md`.

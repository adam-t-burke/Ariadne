* Shared 4-vCPU cloud VM (no CPU governor exposed). Another agent's single-threaded
  process ran at 100% of one core for most of the sweep and compiled with all cores
  at times; passes 2–4 re-measured cells, and the selection rule above picks the
  least-disturbed run per cell. The grid 708 (1M edges) single-thread evaluation is
  950 ms here against 877 ms in `BENCHMARKS.md` measured on an idle VM.
* 1 and 4 threads give the same times to within noise: the numeric factorisation
  (~75% of an evaluation) runs single-threaded (`Parallelism::None`, see
  `BENCHMARKS.md`), so `RAYON_NUM_THREADS` only touches the O(ne) loops.
* Irregular and dome fixtures at matched edge counts have 0.5–0.6× the free nodes of
  the grid (higher node degree), hence smaller factors and 1.6–1.9× faster evaluations.
* Fits: RMS relative residuals are below 15% in every cell; the maximum residual
  exceeds 15% only at 10k–50k edges on the grid and irregular fixtures, where the
  10k-edge evaluation is faster than the `n^1.5 + n log n` trend (the factor fits in
  cache) — the model has no term for that regime.

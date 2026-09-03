#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "scipy>=1.17,<1.18",
# ]
# ///
"""Benchmark SciPy's public L-BFGS-B wrapper on the reference Rosenbrock case."""

from __future__ import annotations

import argparse
import platform
import time

import numpy as np
import scipy
from scipy.optimize import fmin_l_bfgs_b


def objective_gradient(x: np.ndarray) -> tuple[float, np.ndarray]:
    residual = x[1:] - x[:-1] * x[:-1]
    value = 4.0 * (0.25 * (x[0] - 1.0) ** 2 + np.dot(residual, residual))
    gradient = np.empty_like(x)
    gradient[0] = 2.0 * (x[0] - 1.0) - 16.0 * x[0] * residual[0]
    gradient[1:-1] = (
        8.0 * residual[:-1] - 16.0 * x[1:-1] * residual[1:]
    )
    gradient[-1] = 8.0 * residual[-1]
    return float(value), gradient


def projected_gradient_norm(
    x: np.ndarray, gradient: np.ndarray, lower: np.ndarray, upper: np.ndarray
) -> float:
    projected = np.where(
        gradient < 0.0,
        np.maximum(gradient, x - upper),
        np.minimum(gradient, x - lower),
    )
    return float(np.max(np.abs(projected)))


def solve(
    n: int, m: int, lower: np.ndarray, upper: np.ndarray
) -> tuple[np.ndarray, float, dict[str, object]]:
    x0 = np.full(n, 3.0)
    return fmin_l_bfgs_b(
        objective_gradient,
        x0,
        bounds=np.column_stack((lower, upper)),
        m=m,
        factr=1.0e7,
        pgtol=1.0e-5,
        maxls=20,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=1000)
    parser.add_argument("--m", type=int, default=10)
    parser.add_argument("--repeats", type=int, default=20)
    args = parser.parse_args()
    if args.n < 2 or args.m < 1 or args.repeats < 1:
        parser.error("n must be >= 2; m and repeats must be positive")

    lower = np.full(args.n, -100.0)
    lower[::2] = 1.0
    upper = np.full(args.n, 100.0)

    solve(args.n, args.m, lower, upper)
    started = time.perf_counter_ns()
    result = None
    for _ in range(args.repeats):
        result = solve(args.n, args.m, lower, upper)
    elapsed_ns = time.perf_counter_ns() - started
    assert result is not None
    x, value, info = result
    gradient = np.asarray(info["grad"])
    checksum = float(np.dot(np.arange(1, args.n + 1, dtype=np.float64), x))
    pg = projected_gradient_norm(x, gradient, lower, upper)

    print(
        f"python={platform.python_version()},numpy={np.__version__},"
        f"scipy={scipy.__version__}"
    )
    print(
        f"scipy-lbfgsb,n={args.n},m={args.m},runs={args.repeats},"
        f"microseconds={elapsed_ns / args.repeats / 1.0e3:.3f},"
        f"iterations={info['nit']},evaluations={info['funcalls']},"
        f"warnflag={info['warnflag']},final_f={value:.17e},"
        f"final_pg={pg:.17e},x_checksum={checksum:.17e}"
    )


if __name__ == "__main__":
    main()

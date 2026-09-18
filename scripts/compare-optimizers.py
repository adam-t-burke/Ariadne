#!/usr/bin/env python3
"""Compare the current Theseus optimizers with an isolated pre-Basin revision."""

import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile


def run(*args, cwd, env=None, capture_output=False):
    result = subprocess.run(args, cwd=cwd, env=env, text=True,
                            capture_output=capture_output)
    if result.returncode and capture_output:
        print(result.stdout, end="")
        print(result.stderr, end="")
    result.check_returncode()
    return result


def build(root, env):
    result = run("cargo", "test", "--locked", "-p", "theseus", "--release",
                 "--test", "bench_release", "--test", "optimizer_reference",
                 "--no-run", "--message-format=json", cwd=root, env=env,
                 capture_output=True)
    executables = {}
    for line in result.stdout.splitlines():
        message = json.loads(line)
        if message.get("executable"):
            executables[message["target"]["name"]] = message["executable"]
    return executables


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", default="c8e957c",
                        help="pre-Basin commit (default: c8e957c)")
    args = parser.parse_args()
    root = Path(__file__).resolve().parent.parent
    env = dict(os.environ, RAYON_NUM_THREADS="1")
    with tempfile.TemporaryDirectory(prefix="ariadne-optimizers-") as directory:
        baseline = Path(directory) / "baseline"
        run("git", "worktree", "add", "--detach", str(baseline), args.baseline, cwd=root)
        try:
            tests = Path("crates/theseus/tests")
            (baseline / tests / "support").mkdir(exist_ok=True)
            for name in ("bench_release.rs", "optimizer_reference.rs"):
                shutil.copyfile(root / tests / name, baseline / tests / name)
            shutil.copyfile(root / tests / "support/lbfgsb_reference.rs",
                            baseline / tests / "support/basin_reference.rs")
            # Build both revisions before timing to keep compilation off the CPU.
            print("Building baseline...", flush=True)
            baseline_bins = build(baseline, env)
            print("Building Basin integration...", flush=True)
            basin_bins = build(root, env)
            for label, workspace, executables in (
                ("baseline", baseline, baseline_bins), ("basin", root, basin_bins)
            ):
                print(f"\n{label}: RAYON_NUM_THREADS=1; one warmup, ten samples", flush=True)
                run(executables["optimizer_reference"], "bounded_reference_solutions", cwd=workspace, env=env)
                for target, test in (("bench_release", "bench_optimizer_comparison"),
                                     ("optimizer_reference", "bench_reference_comparison")):
                    run(executables[target], test, "--ignored", "--nocapture", "--test-threads=1",
                        cwd=workspace, env=env)
        finally:
            run("git", "worktree", "remove", "--force", str(baseline), cwd=root)


if __name__ == "__main__":
    main()

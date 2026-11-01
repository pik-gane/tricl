#!/usr/bin/env python3
"""Regenerate all numbers and data behind the figures and tables of the paper (written to paper/data/).

Usage (from the repository root, after building tricl):
    python3 paper/make_data.py [--only NAME ...] [--jobs N]

Parts (names for --only): recovery, macro, atop, benchmark, frames. The alliance part downloads the ATOP and COW
data from the peacesciencer repository into paper/data/atop_raw/ if missing (needs pyreadr and pandas).
Everything is deterministic given the seeds used here (and the libstdc++ random number generator).
"""

import argparse
import concurrent.futures
import csv
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "paper", "data")
BIN = os.path.join(ROOT, "build", "src", "tricl")
PY = sys.executable


def run(cmd, cwd=None, log=None):
    p = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    if log:
        with open(log, "w") as f:
            f.write(p.stdout + p.stderr)
    if p.returncode != 0:
        raise RuntimeError("%s failed:\n%s" % (" ".join(cmd), (p.stdout + p.stderr)[-2000:]))
    return p.stdout


def recovery(jobs):
    """Parameter recovery on the three-entity triangle-closure model, observed for 400 time units."""
    work = tempfile.mkdtemp(prefix="tricl_paper_rec_")
    src = open(os.path.join(ROOT, "tests", "configs", "three_entities.yaml")).read()
    cfg = os.path.join(work, "model.yaml")
    with open(cfg, "w") as f:
        f.write(re.sub(r"\n    t: 40\n", "\n    t: 400\n", src))
    out = os.path.join(DATA, "recovery.json")
    run([PY, os.path.join(ROOT, "python", "tricl_recovery.py"), cfg, "--fit", "A", "K", "D", "BK",
         "--seeds", *[str(s) for s in range(1, 41)], "--workdir", work, "--json", out],
        log=os.path.join(DATA, "recovery.log"))


def macro(jobs):
    """Mean-field approximation versus simulations for three models."""
    cases = [("linear", "tests/configs/mean_field.yaml", 10, 0.5, 20),
             ("si", "config_files/si.yaml", 20, 0.5, 10),
             ("granovetter", "config_files/granovetter_helfmann.yaml", 30, 1.0, 10)]
    def one(case):
        name, cfg, t, dt, n = case
        run([PY, os.path.join(ROOT, "python", "tricl_macro.py"), os.path.join(ROOT, cfg), "--t-max", str(t), "--dt", str(dt),
             "--seeds", *[str(s) for s in range(1, n + 1)], "--out", os.path.join(DATA, "macro_%s.csv" % name)],
            log=os.path.join(DATA, "macro_%s.log" % name))
    with concurrent.futures.ThreadPoolExecutor(jobs) as ex:
        list(ex.map(one, cases))


def atop(jobs):
    """Alliance data: model, fit with standard errors, restricted fits, profile likelihood of b1."""
    raw = os.path.join(DATA, "atop_raw")
    out = os.path.join(DATA, "atop")
    os.makedirs(out, exist_ok=True)
    run([PY, os.path.join(ROOT, "python", "tricl_atop.py"), "--data", raw, "--out", out], log=os.path.join(out, "build.log"))
    shuffled = os.path.join(DATA, "atop_shuffled")
    run([PY, os.path.join(ROOT, "python", "tricl_atop.py"), "--data", raw, "--out", shuffled, "--shuffle", "1"])
    fit = os.path.join(ROOT, "python", "tricl_fit.py")
    allp = ["entry", "exit", "a0_defense", "a1_defense", "d0_defense", "b1_defense"]
    def f(args):
        tag, model_dir, names, starts = args
        cmd = [PY, fit, os.path.join(model_dir, "model.yaml"), os.path.join(model_dir, "events.csv"), "--fit", *names,
               "--json", os.path.join(out, tag + ".json")]
        if starts:
            cmd += ["--start", *starts]
        if tag in ("full", "shuffled"):
            cmd += ["--se"]
        run(cmd, log=os.path.join(out, tag + ".log"))
    tasks = [("full", out, allp, []), ("shuffled", shuffled, allp, []),
             ("no_a1", out, [p for p in allp if p != "a1_defense"], ["a1_defense=1e-9"]),
             ("no_b1", out, [p for p in allp if p != "b1_defense"], ["b1_defense=1e-9"])]
    grid = [-3, -1, -0.3, -0.1, -0.03, -0.01, 0.0, 0.01, 0.02, 0.03, 0.035, 0.04, 0.045, 0.05, 0.06, 0.08, 0.1]
    for b in grid:
        tasks.append(("profile_%g" % b, out, [p for p in allp if p != "b1_defense"], ["b1_defense=%g" % (b if b != 0 else 1e-9)]))
    with concurrent.futures.ThreadPoolExecutor(jobs) as ex:
        list(ex.map(f, tasks))
    with open(os.path.join(out, "profile_b1.csv"), "w", newline="") as fo:
        w = csv.writer(fo)
        w.writerow(["b1", "logl"])
        for b in grid:
            w.writerow([b, json.load(open(os.path.join(out, "profile_%g.json" % b)))["logl"]])


def benchmark(jobs):
    run([PY, os.path.join(ROOT, "tests", "benchmark.py"), "--events", "20000", "--repeat", "3",
         "--out", os.path.join(DATA, "benchmark.csv")], log=os.path.join(DATA, "benchmark.log"))


def frames(jobs):
    """Three frames of the Granovetter model for the illustration figure."""
    work = tempfile.mkdtemp(prefix="tricl_paper_frames_")
    run([BIN, os.path.join(ROOT, "config_files", "granovetter_helfmann.yaml"), "--seed", "3", "--quiet", "--max-t", "40",
         "--links-out", "links.csv", "--entities-out", "entities.csv"], cwd=work)
    fdir = os.path.join(DATA, "frames")
    shutil.rmtree(fdir, ignore_errors=True)
    run([PY, os.path.join(ROOT, "python", "tricl_movie.py"), "entities.csv", "links.csv", "--state", "is", "is not",
         "--t0", "0", "--t1", "40", "--frames", "3", "--size", "4", "4", "--dpi", "150", "--out", fdir], cwd=work)
    shutil.copy(os.path.join(work, "links.csv"), os.path.join(DATA, "frames_links.csv"))
    shutil.copy(os.path.join(work, "entities.csv"), os.path.join(DATA, "frames_entities.csv"))


PARTS = {"recovery": recovery, "macro": macro, "atop": atop, "benchmark": benchmark, "frames": frames}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", nargs="+", default=list(PARTS), choices=list(PARTS))
    ap.add_argument("--jobs", type=int, default=4)
    args = ap.parse_args()
    os.makedirs(DATA, exist_ok=True)
    for name in args.only:
        print("==", name, flush=True)
        PARTS[name](args.jobs)
    return 0


if __name__ == "__main__":
    sys.exit(main())

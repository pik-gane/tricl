#!/usr/bin/env python3
"""Parameter-recovery study: simulate with known metaparameters, re-estimate them, and summarise.

Usage:
    python3 python/tricl_recovery.py CONFIG --fit NAME [NAME ...] [--seeds 1 2 3 ...] [--bin PATH]
                                     [--start-factor F] [--maxiter N] [--workdir DIR] [--json OUT.json]

For each seed, the model in CONFIG is simulated with its metaparameter values (the "true" values), the
performed events are written to a csv file, and the metaparameters NAME ... are re-estimated from these
events with tricl_fit.py, starting from the true values multiplied (alternatingly) by F and 1/F.
The summary reports, per parameter, the true value, the mean and standard deviation of the estimates
across seeds, the mean estimated standard error, and how often the true value lies within
estimate +- 1.96 standard errors.
"""

import argparse
import json
import math
import os
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tricl_fit import TriclModel, fit, standard_errors, read_metaparameters  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("config")
    ap.add_argument("--fit", nargs="+", required=True)
    ap.add_argument("--seeds", nargs="+", type=int, default=[1, 2, 3, 4, 5])
    ap.add_argument("--bin", default=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "build", "src", "tricl"))
    ap.add_argument("--start-factor", type=float, default=1.5)
    ap.add_argument("--maxiter", type=int, default=200)
    ap.add_argument("--workdir", default="")
    ap.add_argument("--json", default="")
    args = ap.parse_args()

    binary = os.path.abspath(args.bin)
    config = os.path.abspath(args.config)
    workdir = args.workdir or tempfile.mkdtemp(prefix="tricl_recovery_")
    os.makedirs(workdir, exist_ok=True)
    config_meta = read_metaparameters(config)
    true = {name: float(config_meta[name]) for name in args.fit}

    results = []
    for k, seed in enumerate(args.seeds):
        events = os.path.join(workdir, "events_seed%d.csv" % seed)
        p = subprocess.run([binary, config, "--summary", "--seed", str(seed), "--events-out", events],
                           cwd=workdir, capture_output=True, text=True)
        if p.returncode != 0:
            print("simulation with seed %d failed:\n%s" % (seed, p.stderr[-1000:]), file=sys.stderr)
            continue
        n_events = json.loads(p.stdout.strip().split("\n")[-1])["events"]
        model = TriclModel(binary, config, events, quiet=True)
        start = {name: true[name] * (args.start_factor if (k + i) % 2 == 0 else 1 / args.start_factor) for i, name in enumerate(args.fit)}
        estimates, logl, grad, n_iter = fit(model, start, args.fit, maxiter=args.maxiter, quiet=True)
        se = standard_errors(model, start, args.fit, estimates)
        results.append({"seed": seed, "events": n_events, "estimates": estimates, "se": se, "logl": logl, "iterations": n_iter})
        print("seed %d: %d events, %d iterations, " % (seed, n_events, n_iter)
              + ", ".join("%s=%.4g (se %.2g)" % (name, estimates[name], se[name]) for name in args.fit), flush=True)

    if not results:
        return 1
    print("\n%-12s %12s %12s %12s %12s %10s" % ("parameter", "true", "mean est.", "sd est.", "mean se", "coverage"))
    summary = {}
    for name in args.fit:
        ests = [r["estimates"][name] for r in results]
        ses = [r["se"][name] for r in results if math.isfinite(r["se"][name])]
        mean = sum(ests) / len(ests)
        sd = math.sqrt(sum((e - mean) ** 2 for e in ests) / max(1, len(ests) - 1))
        mean_se = sum(ses) / len(ses) if ses else float("nan")
        covered = [abs(r["estimates"][name] - true[name]) <= 1.96 * r["se"][name] for r in results if math.isfinite(r["se"][name])]
        coverage = sum(covered) / len(covered) if covered else float("nan")
        summary[name] = {"true": true[name], "mean": mean, "sd": sd, "mean_se": mean_se, "coverage": coverage}
        print("%-12s %12.5g %12.5g %12.3g %12.3g %10.2f" % (name, true[name], mean, sd, mean_se, coverage))
    if args.json:
        with open(args.json, "w") as f:
            json.dump({"true": true, "results": results, "summary": summary}, f, indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())

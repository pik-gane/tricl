#!/usr/bin/env python3
"""Throughput benchmark of tricl: events per second for the example configs at several sizes.

Usage:
    python3 tests/benchmark.py [--bin PATH] [--out FILE.csv] [--events N] [--repeat R] [--quick]

Runs the Helfmann threshold model at N = 100, 1000 and 10000 agents (sparse networks; the largest one has a hub
entity, see issue #4), the 3-blocks model (dense) and the SIR model with social distancing, each for --events
events (default 20000) with a fixed seed, and reports wall time, events per second (net of the initialisation
time, which is measured by a run limited to one event) and the peak memory of the child processes so far (i.e. of
the largest case run up to that point). Results are printed and appended to --out as csv rows with a timestamp and
the git revision, so that the file accumulates a history.
"""

import argparse
import csv
import datetime
import json
import os
import resource
import subprocess
import sys
import tempfile
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

CASES = [  # (name, config, extra arguments)
    ("helfmann N=100", "config_files/granovetter_helfmann.yaml", ["--Na", "20", "--Nc", "60", "--Nn", "20"]),
    ("helfmann N=1000", "config_files/granovetter_helfmann.yaml", ["--Na", "200", "--Nc", "600", "--Nn", "200"]),
    ("helfmann N=10000", "config_files/granovetter_helfmann.yaml", ["--Na", "2000", "--Nc", "6000", "--Nn", "2000"]),
    ("3blocks N=300 (dense)", "config_files/parameters_3blocks.yaml", []),
    ("sir_sd N=1000", "config_files/sir_sd.yaml", []),
]


def run(binary, config, args, cwd):
    """Returns (wall seconds, peak RSS in MB, stdout)."""
    t0 = time.perf_counter()
    p = subprocess.run([binary, os.path.join(ROOT, config), "--summary", "--seed", "1"] + args, cwd=cwd, capture_output=True, text=True)
    secs = time.perf_counter() - t0
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    if p.returncode != 0:
        raise RuntimeError("%s failed:\n%s" % (config, p.stderr[-1000:]))
    return secs, after.ru_maxrss / 1024.0, p.stdout


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", default=os.path.join(ROOT, "build", "src", "tricl"))
    ap.add_argument("--out", default=None, help="csv file to append the results to")
    ap.add_argument("--events", type=int, default=20000)
    ap.add_argument("--repeat", type=int, default=1, help="repetitions (the fastest is reported)")
    ap.add_argument("--quick", action="store_true", help="skip the N=10000 case")
    args = ap.parse_args()
    binary = os.path.abspath(args.bin)
    try:
        rev = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    except Exception:
        rev = ""
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    cwd = tempfile.mkdtemp(prefix="tricl_bench_")
    rows = []
    print("%-24s %10s %10s %12s %10s %10s" % ("case", "events", "init [s]", "run [s]", "events/s", "RSS [MB]"))
    for name, config, extra in CASES:
        if args.quick and "10000" in name:
            continue
        common = extra + ["--quiet"]
        best = None
        for _ in range(args.repeat):
            secs, rss, out = run(binary, config, common + ["--max-events", str(args.events)], cwd)
            if best is None or secs < best[0]:
                best = (secs, rss, out)
        secs, rss, out = best
        summary = json.loads(out.strip().splitlines()[-1])
        n = summary["events"]
        # initialisation time: a run limited to 1 event
        t_init = min(run(binary, config, common + ["--max-events", "1"], cwd)[0] for _ in range(args.repeat))
        rate = n / max(secs - t_init, 1e-9)
        print("%-24s %10d %10.2f %12.2f %10.0f %10.0f" % (name, n, t_init, secs, rate, rss))
        rows.append([stamp, rev, name, n, "%.3f" % t_init, "%.3f" % secs, "%.0f" % rate, "%.0f" % rss])
    if args.out:
        new = not os.path.exists(args.out)
        with open(args.out, "a", newline="") as f:
            w = csv.writer(f)
            if new:
                w.writerow(["timestamp", "revision", "case", "events", "init_seconds", "seconds", "events_per_second", "peak_rss_mb"])
            w.writerows(rows)
        print("appended to", args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())

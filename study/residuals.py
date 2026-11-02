#!/usr/bin/env python3
"""Time-rescaling residuals of contact terminations on the held-out days for each fitted model (on the observed
paths): under a correct model, the integrated termination hazards of the contacts are unit exponential.
Writes study/results/residuals.json (KS statistics, means) and residual samples (residual_samples.npz)."""
import csv
import json
import os
import subprocess
import sys

import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fit  # noqa: E402
import c1  # noqa: E402

DAYS = [5, 6]


def main(models_=("M0", "M1", "M2", "M4", "M3", "M6", "M5", "M7", "M8", "C1")):
    out, samples = {}, {}
    for name in models_:
        p = os.path.join(HERE, "results", "fit_%s.json" % name)
        if not os.path.exists(p):
            continue
        r = json.load(open(p))
        res = []
        if name == "C1":
            m = c1.C1(DAYS, r["age0"])
            th = np.array([r["estimates"][k] for k in m.names])
            for d in DAYS:
                res.append(m.residuals(th, d))
        else:
            args = []
            for k, v in r["estimates"].items():
                args += ["--" + k, repr(float(v))]
            for d in DAYS:
                rf = "/tmp/resid_%s_%d.csv" % (name, d)
                subprocess.run([fit.BIN, os.path.join(HERE, "configs", "%s.yaml" % name), "--quiet", "--events-in", os.path.join(HERE, "data", "events_%d.csv" % d),
                                "--forced-in", os.path.join(HERE, "data", "forced_%d.csv" % d), "--max-t", "540", "--residuals-out", rf] + args,
                               check=True, capture_output=True)
                res.append(np.array([float(x["compensator"]) for x in csv.DictReader(open(rf))]))
        x = np.concatenate(res)
        ks = stats.kstest(x, "expon")
        out[name] = {"n": int(len(x)), "mean": float(x.mean()), "KS": float(ks.statistic),
                     "q50": float(np.median(x)), "P_lt_0.1": float((x < 0.1).mean()), "P_gt_3": float((x > 3).mean())}
        samples[name] = x
        print(name, out[name])
    json.dump(out, open(os.path.join(HERE, "results", "residuals.json"), "w"), indent=1)
    np.savez_compressed(os.path.join(HERE, "results", "residual_samples.npz"), **samples)


if __name__ == "__main__":
    main()

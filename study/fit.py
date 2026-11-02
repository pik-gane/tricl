#!/usr/bin/env python3
"""Maximum-likelihood fits of the contact models M0-M3 with tricl (replay with analytic gradients), summed over days.

Usage: python3 study/fit.py MODEL [--train 2 3 4] [--test 5 6]
Writes study/results/fit_M<MODEL>.json with estimates, standard errors (observed information), the training
log-likelihood, the held-out log-likelihood per test day, and the number of events.
"""
import argparse
import json
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

import numpy as np
from scipy.optimize import minimize

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import models  # noqa: E402

BIN = os.path.join(HERE, "..", "build", "src", "tricl")
DATA = os.path.join(HERE, os.environ.get("STUDY_DATA", "data"))
T_DAY = {2: 300.0}


def run(cfg, day, vals, grad=False, extra=()):
    args = [BIN, cfg, "--summary", "--events-in", os.path.join(DATA, "events_%d.csv" % day),
            "--forced-in", os.path.join(DATA, "forced_%d.csv" % day), "--max-t", str(T_DAY.get(day, 540.0))]
    if grad:
        args.append("--grad")
    for k, v in vals.items():
        args += ["--" + k, repr(float(v))]
    p = subprocess.run(args + list(extra), capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(p.stdout[-500:] + p.stderr[-1000:])
    return json.loads(p.stdout.strip().split("\n")[-1])


def jacobian(cfg, names, vals):
    """d coefficient / d metaparameter, from --dump-parameters (the coefficients are linear in the metaparameters)."""
    def dump(v):
        args = [BIN, cfg, "--dump-parameters"]
        for k, x in v.items():
            args += ["--" + k, repr(float(x))]
        return json.loads(subprocess.run(args, capture_output=True, text=True).stdout.strip().split("\n")[-1])["parameters"]
    base = dump(vals)
    J = {}
    for m in names:
        v = dict(vals); v[m] += 1.0
        d = dump(v)
        J[m] = {lab: d[lab] - base[lab] for lab in base if np.isfinite(base[lab]) and d[lab] != base[lab]}
    return J


class Objective:
    def __init__(self, model, days):
        self.cfg = os.path.join(HERE, "configs", "M%d.yaml" % model)
        self.names = list(models.params(model).keys())
        self.days = days
        p = models.params(model)
        self.J = jacobian(self.cfg, self.names, {m: p[m][2] for m in self.names})
        self.cache = {}

    def evaluate(self, x, days=None):
        key = (tuple(np.round(x, 12)), tuple(days or self.days))
        if key in self.cache:
            return self.cache[key]
        vals = dict(zip(self.names, x))
        with ThreadPoolExecutor(max_workers=4) as ex:
            res = list(ex.map(lambda d: run(self.cfg, d, vals, grad=True), days or self.days))
        L = sum(r["logl"] for r in res)
        g = np.zeros(len(self.names))
        for r in res:
            for i, m in enumerate(self.names):
                g[i] += sum(r["gradient"].get(lab, 0.0) * c for lab, c in self.J[m].items())
        self.cache[key] = (L, g, sum(r["events"] for r in res))
        return self.cache[key]

    def f(self, x):
        L, g, n = self.evaluate(x)
        return -L, -g


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("model", type=int)
    ap.add_argument("--train", type=int, nargs="+", default=[2, 3, 4])
    ap.add_argument("--test", type=int, nargs="+", default=[5, 6])
    ap.add_argument("--tag", default="")
    ap.add_argument("--bound", nargs="*", default=[], help="upper bounds NAME=VALUE (e.g. F_tri=0.5) for a constrained fit")
    ap.add_argument("--start", default="", help="results file to start from (default: the next smaller model)")
    args = ap.parse_args()
    upper = {kv.split("=")[0]: float(kv.split("=")[1]) for kv in args.bound}
    obj = Objective(args.model, args.train)
    p = models.params(args.model)
    x0 = np.array([p[m][2] for m in obj.names])
    # start from the estimates of the next smaller model if available:
    prev = os.path.join(HERE, "results", "fit_M%d%s.json" % ({4: 1, 5: 3, 6: 4, 7: 1, 8: 6}.get(args.model, args.model - 1), args.tag))
    if args.start:
        prev = args.start
    if os.path.exists(prev):
        est = json.load(open(prev))["estimates"]
        x0 = np.array([est.get(m, x0[i]) for i, m in enumerate(obj.names)])
    bounds = [(0.0 if m in models.ATTEMPT_PARAMS else None, upper.get(m)) for m in obj.names]
    x0 = np.array([min(x, upper.get(m, np.inf)) for x, m in zip(x0, obj.names)])
    r = minimize(obj.f, x0, jac=True, method="L-BFGS-B", bounds=bounds, options={"maxiter": 500, "gtol": 1e-5})
    L, g, n = obj.evaluate(r.x)
    # observed information by central finite differences of the analytic gradient:
    k = len(r.x)
    H = np.zeros((k, k))
    for i in range(k):
        h = 1e-4 * max(1.0, abs(r.x[i]))
        xp, xm = r.x.copy(), r.x.copy()
        xp[i] += h; xm[i] -= h
        H[i] = (obj.evaluate(xp)[1] - obj.evaluate(xm)[1]) / (2 * h)
    H = (H + H.T) / 2
    try:
        cov = np.linalg.inv(-H)
        se = np.sqrt(np.maximum(np.diag(cov), 0))
    except np.linalg.LinAlgError:
        se = np.full(k, np.nan)
    test = {}
    for d in args.test:
        Lt, gt, nt = obj.evaluate(r.x, days=[d])
        test[d] = {"logl": Lt, "events": nt}
    out = {"model": args.model, "train_days": args.train, "test_days": args.test, "converged": bool(r.success),
           "message": str(r.message), "iterations": int(r.nit), "logl_train": L, "events_train": n,
           "gradient_norm": float(np.max(np.abs(g))),
           "estimates": dict(zip(obj.names, map(float, r.x))), "standard_errors": dict(zip(obj.names, map(float, se))),
           "test": test, "logl_test": sum(t["logl"] for t in test.values()), "events_test": sum(t["events"] for t in test.values()),
           "n_params": k}
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    json.dump(out, open(os.path.join(HERE, "results", "fit_M%d%s.json" % (args.model, args.tag)), "w"), indent=1)
    print("M%d: logl train %.1f (%d events), test %.1f (%d events), %d params, converged %s, |g| %.2g"
          % (args.model, L, n, out["logl_test"], out["events_test"], k, r.success, out["gradient_norm"]))


if __name__ == "__main__":
    main()

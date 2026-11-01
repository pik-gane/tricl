#!/usr/bin/env python3
"""Maximum-likelihood estimation of tricl metaparameters from an observed event sequence.

Usage:
    python3 python/tricl_fit.py CONFIG EVENTS.csv --fit NAME [NAME ...] [--bin PATH]
                                [--start NAME=VALUE ...] [--maxiter N] [--se] [--json OUT.json] [--quiet]

CONFIG is a tricl config file whose 'dynamics' section uses the metaparameters NAME ... (its initial links
describe the state at t = 0 and its limits:t the end of the observation window), and EVENTS.csv is a complete
event sequence in the format written by tricl --events-out (columns t, event, source, relationship, target,
including the "immediate" events the model performs in response to other events).

The log-likelihood of the events and its gradient w.r.t. the model parameters are computed by
    tricl CONFIG --events-in EVENTS.csv --summary --grad --NAME VALUE ...
and mapped to the metaparameters by the chain rule, using the Jacobian of the model parameters w.r.t. the
metaparameters obtained from finite differences of `tricl CONFIG --dump-parameters`. Metaparameters that only
enter attempt rates are optimised on a log scale so that they stay positive. scipy's L-BFGS-B is used if scipy
is installed, otherwise a simple gradient ascent with backtracking line search.

With --se, standard errors are estimated from the observed information (the negative Hessian of the
log-likelihood, obtained by finite differences of the analytic gradient).

Note that model parameters can only be estimated if their current value is nonzero (for attempt rates: positive),
since angles without effect are not tracked by the simulator.
"""

import argparse
import json
import math
import os
import subprocess
import sys


class TriclModel:
    """Wraps calls of the tricl binary for a config file and an events file."""

    def __init__(self, binary, config, events, extra_args=(), quiet=False):
        self.binary = binary
        self.config = config
        self.events = events
        self.extra_args = list(extra_args)
        self.quiet = quiet
        self.n_evaluations = 0

    def _run(self, args):
        p = subprocess.run([self.binary, self.config] + args + self.extra_args, capture_output=True, text=True)
        if p.returncode != 0:
            raise RuntimeError("tricl failed (exit code %d):\n%s" % (p.returncode, p.stderr[-2000:]))
        for line in reversed(p.stdout.strip().split("\n")):
            if line.startswith("{"):
                return json.loads(line)
        raise RuntimeError("no JSON output from tricl:\n" + p.stdout[-2000:])

    @staticmethod
    def _meta_args(meta):
        args = []
        for name, value in meta.items():
            args += ["--" + name, repr(float(value))]
        return args

    def parameters(self, meta):
        """Model parameters (label -> value) for given metaparameter values."""
        return self._run(["--dump-parameters"] + self._meta_args(meta))["parameters"]

    def loglikelihood(self, meta, gradient=True):
        """Log-likelihood of the events (and its gradient w.r.t. the model parameters, label -> value)."""
        self.n_evaluations += 1
        args = ["--summary", "--events-in", self.events] + (["--grad"] if gradient else []) + self._meta_args(meta)
        result = self._run(args)
        return result["logl"], result.get("gradient", {})

    def jacobian(self, meta, names):
        """d(model parameter) / d(metaparameter) by central finite differences: {label: {name: derivative}}."""
        base = self.parameters(meta)
        jac = {label: {} for label in base}
        for name in names:
            v = float(meta[name])
            h = 1e-6 * max(1.0, abs(v))
            plus = self.parameters(dict(meta, **{name: v + h}))
            minus = self.parameters(dict(meta, **{name: v - h}))
            for label in base:
                if math.isfinite(plus[label]) and math.isfinite(minus[label]):
                    d = (plus[label] - minus[label]) / (2 * h)
                    if d != 0.0:
                        jac[label][name] = d
        return jac


def fit(model, meta, names, maxiter=200, quiet=False, tol=1e-6):
    """Maximise the log-likelihood over the metaparameters `names`, starting from the values in `meta`.

    Returns (estimates dict, final log-likelihood, final gradient dict, number of iterations).
    """
    jac = model.jacobian(meta, names)
    # metaparameters that only enter attempt rates are optimised on a log scale:
    log_scale = {}
    for name in names:
        labels = [label for label in jac if name in jac[label]]
        if not labels:
            raise ValueError("metaparameter %s has no effect on any model parameter" % name)
        log_scale[name] = all(("| base attempt" in label) or ("| attempt via" in label) for label in labels)
        if log_scale[name] and not float(meta[name]) > 0:
            raise ValueError("metaparameter %s enters attempt rates and must therefore start from a positive value" % name)
    current = dict(meta)

    def to_meta(u):
        m = dict(current)
        for i, name in enumerate(names):
            m[name] = math.exp(u[i]) if log_scale[name] else u[i]
        return m

    def from_meta(m):
        return [math.log(m[name]) if log_scale[name] else float(m[name]) for name in names]

    def objective(u):
        m = to_meta(u)
        try:
            logl, grad_raw = model.loglikelihood(m)
        except RuntimeError as e:
            if not quiet:
                print("  evaluation failed at %s: %s" % (m, e), file=sys.stderr)
            return float("inf"), [0.0] * len(names)
        jac_m = model.jacobian(m, names)  # the Jacobian may depend on the metaparameters (nonlinear expressions)
        grad_meta = []
        for name in names:
            g = sum(jac_m[label][name] * grad_raw[label] for label in grad_raw if name in jac_m[label])
            if log_scale[name]:
                g *= m[name]
            grad_meta.append(g)
        if not quiet:
            print("  logl %.6f at %s" % (logl, ", ".join("%s=%.6g" % (n, m[n]) for n in names)), flush=True)
        return -logl, [-g for g in grad_meta]

    u0 = from_meta(current)
    n_iter = 0
    try:
        import numpy as np
        from scipy.optimize import minimize
        res = minimize(lambda u: objective(list(u))[0], np.array(u0), jac=lambda u: np.array(objective(list(u))[1]),
                       method="L-BFGS-B", options={"maxiter": maxiter, "gtol": tol})
        u = list(res.x)
        n_iter = int(res.nit)
    except ImportError:
        # simple gradient ascent with backtracking line search:
        u = list(u0)
        f, g = objective(u)
        step = 1.0
        for n_iter in range(1, maxiter + 1):
            gnorm = math.sqrt(sum(x * x for x in g))
            if gnorm < tol:
                break
            direction = [-x / gnorm for x in g]  # steepest descent of the negative log-likelihood
            while True:
                u_new = [ui + step * di for ui, di in zip(u, direction)]
                f_new, g_new = objective(u_new)
                if f_new < f - 1e-4 * step * gnorm:
                    u, f, g = u_new, f_new, g_new
                    step *= 1.5
                    break
                step *= 0.5
                if step < 1e-12:
                    break
            if step < 1e-12:
                break
    estimates = to_meta(u)
    logl, grad_raw = model.loglikelihood(estimates)
    jac_m = model.jacobian(estimates, names)
    grad_meta = {name: sum(jac_m[label][name] * grad_raw[label] for label in grad_raw if name in jac_m[label]) for name in names}
    return {name: estimates[name] for name in names}, logl, grad_meta, n_iter


def standard_errors(model, meta, names, estimates):
    """Standard errors from the observed information, using finite differences of the analytic gradient."""
    m0 = dict(meta, **estimates)

    def grad_meta(m):
        logl, grad_raw = model.loglikelihood(m)
        jac_m = model.jacobian(m, names)
        return [sum(jac_m[label][name] * grad_raw[label] for label in grad_raw if name in jac_m[label]) for name in names]

    k = len(names)
    hessian = [[0.0] * k for _ in range(k)]
    for j, name in enumerate(names):
        v = float(m0[name])
        h = 1e-4 * max(1.0, abs(v))
        gp = grad_meta(dict(m0, **{name: v + h}))
        gm = grad_meta(dict(m0, **{name: v - h}))
        for i in range(k):
            hessian[i][j] = (gp[i] - gm[i]) / (2 * h)
    # symmetrise and invert the observed information -H (Gauss-Jordan, small matrices only):
    info = [[-(hessian[i][j] + hessian[j][i]) / 2 for j in range(k)] for i in range(k)]
    aug = [row[:] + [1.0 if i == j else 0.0 for j in range(k)] for i, row in enumerate(info)]
    for c in range(k):
        pivot = max(range(c, k), key=lambda r: abs(aug[r][c]))
        if abs(aug[pivot][c]) < 1e-300:
            return {name: float("nan") for name in names}
        aug[c], aug[pivot] = aug[pivot], aug[c]
        p = aug[c][c]
        aug[c] = [x / p for x in aug[c]]
        for r in range(k):
            if r != c and aug[r][c] != 0.0:
                f = aug[r][c]
                aug[r] = [x - f * y for x, y in zip(aug[r], aug[c])]
    cov = [row[k:] for row in aug]
    return {name: (math.sqrt(cov[i][i]) if cov[i][i] > 0 else float("nan")) for i, name in enumerate(names)}


def read_metaparameters(config):
    """Read the metaparameter names and default expressions from a config file (simple line-based parser)."""
    meta = {}
    with open(config) as f:
        lines = f.read().split("\n")
    inside = False
    for line in lines:
        stripped = line.split("#")[0].rstrip()
        if not stripped.strip():
            continue
        if not line[0].isspace():
            inside = stripped.startswith("metaparameters")
            continue
        if inside and ":" in stripped:
            name, expr = stripped.strip().split(":", 1)
            meta[name.strip()] = expr.strip()
    return meta


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("config")
    ap.add_argument("events")
    ap.add_argument("--fit", nargs="+", required=True, help="metaparameters to estimate")
    ap.add_argument("--bin", default=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "build", "src", "tricl"))
    ap.add_argument("--start", nargs="*", default=[], help="starting values NAME=VALUE (default: the config file's values)")
    ap.add_argument("--maxiter", type=int, default=200)
    ap.add_argument("--se", action="store_true", help="also estimate standard errors from the observed information")
    ap.add_argument("--json", default="", help="write the results to this JSON file")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    model = TriclModel(os.path.abspath(args.bin), args.config, args.events, quiet=args.quiet)
    config_meta = read_metaparameters(args.config)
    for name in args.fit:
        if name not in config_meta:
            print("metaparameter %s is not defined in %s" % (name, args.config), file=sys.stderr)
            return 2
    # starting values: numeric metaparameters from the config, overridden by --start:
    start = {}
    for name in args.fit:
        try:
            start[name] = float(config_meta[name])
        except ValueError:
            pass
    for item in args.start:
        name, value = item.split("=", 1)
        start[name] = float(value)
    missing = [name for name in args.fit if name not in start]
    if missing:
        print("please give starting values with --start for: " + " ".join(missing), file=sys.stderr)
        return 2

    logl0, _ = model.loglikelihood(start, gradient=False)
    if not args.quiet:
        print("log-likelihood at start: %.6f" % logl0)
    estimates, logl, grad, n_iter = fit(model, start, args.fit, maxiter=args.maxiter, quiet=args.quiet)
    se = standard_errors(model, start, args.fit, estimates) if args.se else {}

    print("%-12s %14s %14s %14s %14s" % ("parameter", "start", "estimate", "std.error" if args.se else "", "gradient"))
    for name in args.fit:
        print("%-12s %14.6g %14.6g %14s %14.3g" % (name, start[name], estimates[name], ("%.6g" % se[name]) if args.se else "", grad[name]))
    print("log-likelihood: %.6f (start: %.6f), %d iterations, %d likelihood evaluations" % (logl, logl0, n_iter, model.n_evaluations))
    if args.json:
        with open(args.json, "w") as f:
            json.dump({"estimates": estimates, "standard_errors": se, "logl": logl, "logl_start": logl0,
                       "gradient": grad, "iterations": n_iter, "evaluations": model.n_evaluations}, f, indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())

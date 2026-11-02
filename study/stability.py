#!/usr/bin/env python3
"""Stability of fitted network-context models by the pair closure (see closure.py): the time until the contacts
within classes tip into a dense state, starting from an empty network with the timetable of a school day repeated,
and a penalty for constrained fits (fit.py --stability).

The tipping time is continuous in the parameters (it diverges at the boundary of the stable region), unlike the
maximal density, which jumps at the tipping point of the bistable dynamics.
"""
import json
import os
import sys

import numpy as np
from scipy.special import expit
from scipy.stats import poisson

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import models  # noqa: E402
import closure  # noqa: E402

DAY = 3            # timetable of a full day (Tuesday)
REPEAT = 3         # days in a row (a safety margin, since the closure ignores fluctuations)
FACTOR = 3.0       # tipping = contacts within a class exceed FACTOR x the largest number observed per class in a day
_cache = {}


def _expect(a, b, base, tri, deg, att, ceil):
    """E[lam], E[lam n], E[lam m] for n ~ Pois(a), m ~ Pois(2b), lam = 2 (ceil + att n) expit(4 (base + tri n + deg (2n + m)))."""
    nmax = int(a + 8 * np.sqrt(a) + 8)
    mmax = int(2 * b + 8 * np.sqrt(2 * b) + 8)
    n = np.arange(nmax)[:, None]
    m = np.arange(mmax)[None, :]
    pn = poisson.pmf(np.arange(nmax), max(a, 1e-12))[:, None]
    pm = poisson.pmf(np.arange(mmax), max(2 * b, 1e-12))[None, :]
    lam = 2 * (ceil + att * n) * expit(4 * (base + tri * n + deg * (2 * n + m)))
    w = pn * pm * lam
    return w.sum(), (w * n).sum(), (w * m).sum()


def tipping_time(est, data, observed_max, h=0.25):
    """Minimum over classes of the time (minutes) until the pair closure exceeds FACTOR x observed_max[class];
    REPEAT * 540 if no class tips. Vectorized over classes and mark combinations; the numbers of contacts are
    advanced with the exact solution for rates frozen during a step, the two-path sums Z, O explicitly."""
    bl, cls = closure.blocks(est, DAY, data)
    e = est
    names = sorted(bl)
    C = len(names)
    K = max(len(bl[c].combos) for c in names)
    cnt = np.zeros((C, K)); Fc = np.zeros((C, K)); Tc = np.zeros((C, K))
    for i, c in enumerate(names):
        for k, (n_, f_, t_) in enumerate(bl[c].combos):
            cnt[i, k], Fc[i, k], Tc[i, k] = n_, f_, t_
    N = np.array([bl[c].N for c in names], float); M = N * (N - 1) / 2
    cap = np.array([FACTOR * observed_max.get(c, 10) for c in names], float)
    L = np.zeros((C, K)); Z = np.zeros(C); O = np.zeros(C)
    tri_f, deg_f, att_f = e.get("F_tri", 0.0), e["F_deg"], e.get("F_tria", 0.0)
    tri_t, deg_t = e.get("T_tri", 0.0), e["T_deg"]
    T_total = REPEAT * 540.0
    t = 0.0
    while t < T_total:
        tod = t % 540.0
        Fb, Tb = closure.slot_effects(e, DAY, tod)
        Lt = L.sum(axis=1); U = M - Lt
        a1 = np.where(Lt > 1e-9, Z / np.maximum(Lt, 1e-9), 0.0)
        a0 = np.where(U > 1e-9, O / np.maximum(U, 1e-9), 0.0)
        sum_d2 = 2 * (Z + O) + 2 * Lt
        b1 = np.where(Lt > 1e-9, np.maximum((Z + O) / np.maximum(Lt, 1e-9) - a1, 0.0), 0.0)
        b0 = np.where(U > 1e-9, np.maximum((2 * Lt * (N - 1) - sum_d2) / (2 * np.maximum(U, 1e-9)) - a0, 0.0), 0.0)
        res = []
        for a, bb, base, tri, deg, att, ceil in ((a0, b0, Fb + Fc, tri_f, deg_f, att_f, models.A_FORM), (a1, b1, Tb + Tc, tri_t, deg_t, 0.0, models.D_TERM)):
            nmax = int(a.max() + 8 * np.sqrt(a.max()) + 8)
            mmax = int(2 * bb.max() + 8 * np.sqrt(2 * bb.max()) + 8)
            nn = np.arange(nmax); mm = np.arange(mmax)
            pn = poisson.pmf(nn[None, :], np.maximum(a, 1e-12)[:, None])          # C x n
            pm = poisson.pmf(mm[None, :], np.maximum(2 * bb, 1e-12)[:, None])     # C x m
            arg = 4 * (base[:, :, None, None] + tri * nn[None, None, :, None] + deg * (2 * nn[None, None, :, None] + mm[None, None, None, :]))
            lam = 2 * (ceil + att * nn[None, None, :, None]) * expit(arg)         # C x K x n x m
            w = pn[:, None, :, None] * pm[:, None, None, :] * lam
            res.append((w.sum(axis=(2, 3)), (w * nn[None, None, :, None]).sum(axis=(2, 3)), (w * mm[None, None, None, :]).sum(axis=(2, 3))))
        (f0, fn, fm), (g0, gn, gm) = res
        Uc = cnt - L
        dZ = (Uc * 3 * fn - L * 3 * gn).sum(axis=1)
        dO = (Uc * (fm - fn) - L * (gm - gn)).sum(axis=1)
        r = f0 + g0
        Lstar = np.where(r > 0, cnt * f0 / np.maximum(r, 1e-300), L)
        L = Lstar + (L - Lstar) * np.exp(-r * h)
        Z = np.maximum(Z + h * dZ, 0.0); O = np.maximum(O + h * dO, 0.0)
        t += h
        if np.any(L.sum(axis=1) > cap):
            return t
    return T_total


def observed_maxima(data):
    import csv
    cls = dict(data["static"]["class"])
    mx = {}
    for d in (2, 3, 4):
        active = set(); cnt = {}
        ev = closure.read_events(os.path.join(HERE, "data", "events_%d.csv" % d))
        for t, kind, p in ev:
            if cls[p[0]] != cls[p[1]]:
                continue
            (active.add if kind == "establish" else active.discard)(p)
            c = cls[p[0]]
            n = sum(1 for q in active if cls[q[0]] == c)
            mx[c] = max(mx.get(c, 0), n)
    return mx


def penalty(est, data, obs, kappa):
    """kappa (1 - T_tip / T_req)^2 if the closure tips before T_req, else 0; also returns T_tip."""
    T_req = REPEAT * 540.0
    t = tipping_time(est, data, obs)
    return kappa * max(0.0, 1 - t / T_req) ** 2, t


if __name__ == "__main__":
    import time
    data = json.load(open(os.path.join(HERE, "data", "data.json")))
    obs = observed_maxima(data)
    print("observed maxima per class", obs)
    for name in sys.argv[1:] or ["M2"]:
        est = json.load(open(os.path.join(HERE, "results", "fit_%s.json" % name)))["estimates"]
        t0 = time.time()
        print(name, "tipping time (min)", tipping_time(est, data, obs), "%.1fs" % (time.time() - t0))
        e2 = dict(est); e2["F_tri"] = 0.0; e2["F_deg"] = -0.2
        t0 = time.time()
        print(name, "with F_tri=0, F_deg=-0.2:", tipping_time(e2, data, obs), "%.1fs" % (time.time() - t0))


def simulated_tipping(config, est, days=(5,), seeds=(1, 2), wall=300):
    """Stability check by simulation: for each day and seed, simulate the model from an empty network with the
    day's timetable; returns the largest number of simultaneous contacts relative to the largest observed number,
    and whether any run tipped (exceeded FACTOR x the observed maximum, or did not finish within the wall time)."""
    import csv, subprocess
    import fit
    worst, tipped = 0.0, False
    for d in days:
        n = 0; mx_obs = 0
        for t, kind, p in closure.read_events(os.path.join(HERE, "data", "events_%d.csv" % d)):
            n += 1 if kind == "establish" else -1; mx_obs = max(mx_obs, n)
        args = []
        for k, v in est.items():
            args += ["--" + k, repr(float(v))]
        for s in seeds:
            sf = "/tmp/stab_%d_%d.csv" % (d, s)
            out = subprocess.run([fit.BIN, config, "--summary", "--forced-in", os.path.join(HERE, "data", "forced_%d.csv" % d), "--max-t", "540",
                                  "--seed", str(s), "--max-wall", str(wall), "--max-events", "200000", "--stats-out", sf] + args, capture_output=True, text=True).stdout
            summ = json.loads(out.strip().split("\n")[-1])
            mx = max(float(r["student meets student"]) / 2 for r in csv.DictReader(open(sf)))
            worst = max(worst, mx / mx_obs)
            if summ["t"] < 540 - 1e-6 or mx > FACTOR * mx_obs:
                tipped = True
    return worst, tipped

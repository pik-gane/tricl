#!/usr/bin/env python3
"""Mean-field closures for the fitted network-context model M2, compared with tricl simulations of the same model
and with the data, for contacts within classes (93% of all contacts).

Each class is a homogeneous block of N students whose pairs differ by static marks (Facebook friendship, survey
friendship 0/1/2). Rates per pair (both directions) as in M2:
    formation   2 A expit(4 (F(t) + F_c + F_tri n + F_deg (d_i + d_j))),
    termination 2 D expit(4 (T(t) + T_c + T_tri n + T_deg (d_i' + d_j'))),
with n the number of common contacts, d the numbers of contacts (d' not counting the pair's own contact), F(t), T(t)
the timetable effects and F_c, T_c the static effects of the pair's marks.

Simple closure: n ~ Poisson((N - 2) rho^2) and the other contacts of each endpoint ~ Poisson((N - 2) rho),
independently of the pair's state (the closure of tricl_macro.py).

Pair closure: in addition to the numbers of contacts L_c by mark combination c, the sums of common contacts over
linked pairs, Z (= 3 x triangles), and over unlinked pairs, O (open two-paths), are state variables. Given the pair's
state, n ~ Poisson(Z / L) or Poisson(O / U), and d = n + m with m ~ Poisson(b), where the excess degrees b follow
from the degree moments implied by Z + O (= number of two-paths = sum_i d_i (d_i - 1) / 2). The expected changes of
Z and O by an event with n common contacts and other contacts m_i + m_j: formation dZ = 3 n, dO = m_i + m_j - n;
termination the negative.
"""
import csv
import json
import os
import subprocess
import sys

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import expit
from scipy.stats import poisson

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import models  # noqa: E402
import fit  # noqa: E402

NMAX, MMAX = 40, 80


def slot_effects(est, day, t):
    start_h = 12 if day == 2 else 8
    h = start_h + int(t // 60)
    mm = t - (h - start_h) * 60
    F, T = est["F0"], est["T0"]
    if h >= 9:
        F += est.get("F_h%02d" % h, 0); T += est.get("T_h%02d" % h, 0)
    if mm < 20:
        F += est["F_chg"]; T += est["T_chg"]
        if h == 10:
            F += est["F_brk"]; T += est["T_brk"]
    return F, T


def expectations(a, b, base, tri, deg, ceil, linked, att=0.0):
    """E[lam], E[lam n], E[lam (m_i + m_j)] for n ~ Pois(a), m_i + m_j ~ Pois(2 b),
    lam = 2 (ceil + att n) expit(4 (base + tri n + deg S))."""
    n = np.arange(NMAX)[:, None]
    m = np.arange(MMAX)[None, :]
    pn = poisson.pmf(np.arange(NMAX), max(a, 1e-12))[:, None]
    pm = poisson.pmf(np.arange(MMAX), max(2 * b, 1e-12))[None, :]
    S = 2 * n + m  # d_i + d_j (or d_i' + d_j' for linked pairs)
    lam = 2 * (ceil + att * n) * expit(4 * (base + tri * n + deg * S))
    w = pn * pm * lam
    return w.sum(), (w * n).sum(), (w * m).sum()


class Block:
    def __init__(self, N, combos, est, day):
        self.N, self.M = N, N * (N - 1) / 2
        self.combos = combos  # list of (count of pairs, F_c, T_c)
        self.e, self.day = est, day

    def rhs(self, t, y, pair):
        e = self.e
        nc = len(self.combos)
        L = np.maximum(y[:nc], 0)
        Ltot = L.sum()
        Z, O = (max(y[nc], 0), max(y[nc + 1], 0)) if pair else (0, 0)
        U = self.M - Ltot
        Fb, Tb = slot_effects(e, self.day, t)
        N = self.N
        rho = Ltot / self.M
        if pair:
            a1 = Z / Ltot if Ltot > 1e-9 else 0.0
            a0 = O / U if U > 1e-9 else 0.0
            kbar = 2 * Ltot / N
            sum_d2 = 2 * (Z + O) + 2 * Ltot  # sum_i d_i^2
            ex_linked = (Z + O) / Ltot if Ltot > 1e-9 else 0.0  # mean other contacts of an endpoint of a linked pair
            ex_unlinked = (2 * Ltot * (N - 1) - sum_d2) / (2 * U) if U > 1e-9 else 0.0  # mean contacts of an endpoint of an unlinked pair
            b1, b0 = max(ex_linked - a1, 0.0), max(ex_unlinked - a0, 0.0)
        else:
            a0 = a1 = (N - 2) * rho ** 2
            b0 = b1 = max((N - 2) * rho - a0, 0.0)
        dL = np.zeros(nc)
        dZ = dO = 0.0
        for c, (cnt, Fc, Tc) in enumerate(self.combos):
            Uc = cnt - L[c]
            f0, fn, fm = expectations(a0, b0, Fb + Fc, e.get("F_tri", 0.0), e["F_deg"], models.A_FORM, False, e.get("F_tria", 0.0))
            g0, gn, gm = expectations(a1, b1, Tb + Tc, e["T_tri"], e["T_deg"], models.D_TERM, True)
            dL[c] = Uc * f0 - L[c] * g0
            dZ += Uc * 3 * fn - L[c] * 3 * gn
            dO += Uc * (fm - fn) - L[c] * (gm - gn)
        return np.concatenate([dL, [dZ, dO]]) if pair else dL

    def integrate(self, T, pair, grid):
        nc = len(self.combos)
        y = np.zeros(nc + 2 if pair else nc)
        out = []
        edges = sorted(set([0.0, T] + list(np.arange(0, T, 20.0))))
        gi = 0
        for a, b in zip(edges[:-1], edges[1:]):
            ts = [g for g in grid if a <= g <= b]
            sol = solve_ivp(lambda t, yy: self.rhs(t, yy, pair), (a, b), y, t_eval=ts or None, rtol=1e-6, atol=1e-8, method="LSODA")
            for k, tt in enumerate(sol.t):
                if ts and (not out or tt > out[-1][0] + 1e-9):
                    yy = sol.y[:, k]
                    out.append((tt, yy[:nc].sum(), yy[nc] / 3 if pair else None))
            y = sol.y[:, -1]
        return out


def blocks(est, day, data):
    cls = dict(data["static"]["class"])
    fb = {tuple(sorted(p)) for p in data["static"]["facebook"]}
    svy = {}
    for a, b in data["static"]["survey"]:
        k = tuple(sorted((a, b))); svy[k] = svy.get(k, 0) + 1
    members = {}
    for s, c in cls.items():
        members.setdefault(c, []).append(s)
    out = {}
    for c, ss in members.items():
        ss = sorted(ss)
        counts = {}
        for i in range(len(ss)):
            for j in range(i + 1, len(ss)):
                k = (ss[i], ss[j])
                key = ((k in fb), svy.get(k, 0))
                counts[key] = counts.get(key, 0) + 1
        combos = [(n, est["F_cls"] + est["F_fb"] * f + est["F_svy"] * s, est["T_cls"] + est["T_fb"] * f + est["T_svy"] * s)
                  for (f, s), n in sorted(counts.items())]
        out[c] = Block(len(ss), combos, est, day)
    return out, cls


def within_class_counts(events, cls, grid):
    """Numbers of within-class contacts and triangles at the grid times."""
    active = set()
    res, ei = [], 0
    events = [e for e in events if cls[e[2][0]] == cls[e[2][1]]]
    for g in grid:
        while ei < len(events) and events[ei][0] <= g:
            t, ev, p = events[ei]
            (active.add if ev == "establish" else active.discard)(p)
            ei += 1
        nb = {}
        for a, b in active:
            nb.setdefault(a, set()).add(b); nb.setdefault(b, set()).add(a)
        tri = sum(len(nb[a] & nb[b]) for a, b in active) / 3
        res.append((g, len(active), tri))
    return res


def read_events(path):
    out = []
    with open(path) as fh:
        for r in csv.DictReader(fh):
            if r["relationship"] == "meets":
                out.append((float(r["t"]), r["event"], tuple(sorted((r["source"], r["target"])))))
    return out


def main():
    model = sys.argv[1] if len(sys.argv) > 1 else "M2"
    day = 5
    T = 540.0
    est = json.load(open(os.path.join(HERE, "results", "fit_%s.json" % model)))["estimates"]
    data = json.load(open(os.path.join(HERE, "data", "data.json")))
    grid = list(np.arange(2.0, T, 4.0))
    bl, cls = blocks(est, day, data)
    res = {"grid": grid}
    for pair, name in ((False, "simple"), (True, "pair")):
        L = np.zeros(len(grid)); Tr = np.zeros(len(grid))
        for c, b in bl.items():
            out = b.integrate(T, pair, grid)
            L += np.array([o[1] for o in out])
            if pair:
                Tr += np.array([o[2] for o in out])
        res[name] = {"L": L.tolist(), "T": Tr.tolist() if pair else None}
        print(name, "mean within-class contacts", L.mean(), "triangles", Tr.mean() if pair else "-")
    # simulations of the fitted M2 (3 seeds) and the data:
    args = []
    for k, v in est.items():
        args += ["--" + k, repr(float(v))]
    sims = []
    for s in (1, 2, 3):
        ef = "/tmp/closure_sim_%d.csv" % s
        subprocess.run([fit.BIN, os.path.join(HERE, "configs", "%s.yaml" % model), "--quiet", "--forced-in", os.path.join(HERE, "data", "forced_%d.csv" % day),
                        "--max-t", str(T), "--seed", str(s), "--events-out", ef, "--max-wall", "120"] + args, check=True, capture_output=True)
        sims.append(within_class_counts(read_events(ef), cls, grid))
    res["simulation"] = {"L": np.mean([[x[1] for x in s] for s in sims], axis=0).tolist(), "T": np.mean([[x[2] for x in s] for s in sims], axis=0).tolist()}
    dat = within_class_counts(read_events(os.path.join(HERE, "data", "events_%d.csv" % day)), cls, grid)
    res["data"] = {"L": [x[1] for x in dat], "T": [x[2] for x in dat]}
    for k in ("simulation", "data"):
        print(k, "mean within-class contacts", np.mean(res[k]["L"]), "triangles", np.mean(res[k]["T"]))
    json.dump(res, open(os.path.join(HERE, "results", "closure_%s.json" % model), "w"))


if __name__ == "__main__":
    main()

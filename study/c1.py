#!/usr/bin/env python3
"""Competing model C1: static context (as in M1) plus intrinsic, non-Markovian aging, fitted by maximum likelihood.

Rates per pair and direction (the same convention as tricl, so that log-likelihoods are comparable: an event
contributes log(rate per direction), the exposure is the pair rate = 2 x rate per direction):

    formation:   A expit(4 B_F),  B_F = M1 terms + F_mb [met on an earlier day]
                                        + (F_first if the pair has not met today, else F_age log(1 + sigma / sigma0)),
    termination: D expit(4 B_T),  B_T = M1 terms + T_mb [met on an earlier day] + T_age log(1 + tau / tau0),

with tau the age of the current contact and sigma the time since the pair's last contact ended. Far from
saturation, the rates are power laws (1 + tau/tau0)^(4 T_age) and (1 + sigma/sigma0)^(4 F_age), the aging of
contacts and inter-contact intervals in the spirit of the reinforcement models of face-to-face contacts.
The integrals are computed by Gauss-Legendre quadrature in u = log(1 + age / age0) on pieces between the
timetable boundaries, with analytic gradients.

Usage: python3 study/c1.py [--train 2 3 4] [--test 5 6] [--age0 1.0]
"""
import argparse
import json
import os
import sys

import numpy as np
from scipy.optimize import minimize
from scipy.special import expit

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import models  # noqa: E402

GX, GW = np.polynomial.legendre.leggauss(8)
DATA = os.path.join(HERE, os.environ.get("STUDY_DATA", "data"))
STATIC = ["cls", "fb", "svy", "mb"]
SLOTS = models.HOURS[1:] + ["chg", "brk"]


def names():
    n = []
    for e in ("F", "T"):
        n += [e + "0"] + ["%s_%s" % (e, s) for s in SLOTS] + ["%s_%s" % (e, s) for s in STATIC] + [e + "_age"]
    n += ["F_first"]
    return n


def load_day(day, data):
    T = 300.0 if day == 2 else 540.0
    start_h = 12 if day == 2 else 8
    students = data["students"]
    idx = {s: i for i, s in enumerate(students)}
    n = len(students)
    cls = {s: c for s, c in data["static"]["class"]}
    fb = {tuple(sorted(p)) for p in data["static"]["facebook"]}
    svy = {}
    for a, b in data["static"]["survey"]:
        k = tuple(sorted((a, b)))
        svy[k] = svy.get(k, 0) + 1
    mb = set()
    import csv
    with open(os.path.join(DATA, "forced_%d.csv" % day)) as fh:
        for r in csv.DictReader(fh):
            if r["relationship"] == "met before":
                mb.add(tuple(sorted((r["source"], r["target"]))))
    events = []
    with open(os.path.join(DATA, "events_%d.csv" % day)) as fh:
        for r in csv.DictReader(fh):
            events.append((float(r["t"]), r["event"], tuple(sorted((r["source"], r["target"])))))
    pairs = [(students[i], students[j]) for i in range(n) for j in range(i + 1, n)]
    pid = {p: k for k, p in enumerate(pairs)}
    X = np.zeros((len(pairs), len(STATIC)))
    for k, (a, b) in enumerate(pairs):
        X[k] = [cls[a] == cls[b], (a, b) in fb, svy.get((a, b), 0), (a, b) in mb]
    # timetable boundaries and slot covariates (hour index 0..8, hour change, morning break):
    bounds = sorted(set([0.0, T] + [(h - start_h) * 60.0 + o for h in range(start_h, 17) for o in (0.0, 20.0)]))
    bounds = [b for b in bounds if b <= T]

    def slot_of(t):
        h = start_h + int(t // 60)
        m = t - (h - start_h) * 60
        return h, m < 20, (h == 10) and (m < 20)
    # per pair: list of (state, t0, t1, age origin) pieces; events give the event terms
    state = {}  # pair -> (linked?, t_since) ; t_since = start of contact or end of last contact (None if not met today)
    pieces = []  # (pair, kind 0=form 1=term, t0, t1, origin or nan)
    evrows = []  # (pair, kind, t, origin or nan)
    last = {}
    for t, ev, p in events:
        k = pid[p]
        lk, since = state.get(k, (False, None))
        t_prev = last.get(k, 0.0)
        pieces.append((k, 1 if lk else 0, t_prev, t, since))
        evrows.append((k, 1 if lk else 0, t, since))
        state[k] = (not lk, t)
        last[k] = t
    for k in range(len(pairs)):
        lk, since = state.get(k, (False, None))
        pieces.append((k, 1 if lk else 0, last.get(k, 0.0), T, since))
    # split pieces at the timetable boundaries:
    B = np.array(bounds)
    rows = []
    for k, kind, t0, t1, since in pieces:
        if t1 <= t0:
            continue
        cuts = [t0] + [b for b in B if t0 < b < t1] + [t1]
        for a, b in zip(cuts[:-1], cuts[1:]):
            h, chg, brk = slot_of((a + b) / 2)
            rows.append((k, kind, a, b, np.nan if since is None else since, h, chg, brk))
    R = np.array(rows, dtype=float)
    E = np.array([(k, kind, t, np.nan if s is None else s) + slot_of(t) for k, kind, t, s in evrows], dtype=float)
    return {"X": X, "R": R, "E": E, "start_h": start_h, "n_events": len(evrows)}


class C1:
    def __init__(self, days, age0):
        self.data = json.load(open(os.path.join(DATA, "data.json")))
        self.days = {d: load_day(d, self.data) for d in days}
        self.names = names()
        self.ix = {n: i for i, n in enumerate(self.names)}
        self.age0 = age0
        self.A, self.D = models.A_FORM, models.D_TERM

    def _features(self, e, k, kind, h, chg, brk):
        """Feature matrix F (rows x params) without the age/first columns, for kind-specific parameters."""
        p = len(self.names)
        F = np.zeros((len(k), p))
        for kd, pre in ((0, "F"), (1, "T")):
            m = kind == kd
            F[m, self.ix[pre + "0"]] = 1
            for i, hh in enumerate(models.HOURS[1:]):
                F[m & (h == int(hh[1:])), self.ix["%s_%s" % (pre, hh)]] = 1
            F[m & (chg > 0), self.ix[pre + "_chg"]] = 1
            F[m & (brk > 0), self.ix[pre + "_brk"]] = 1
            for i, s in enumerate(STATIC):
                F[m, self.ix["%s_%s" % (pre, s)]] = e["X"][k[m], i]
        return F

    def loglik(self, theta, days=None):
        L, G = 0.0, np.zeros(len(theta))
        for d in (days or self.days):
            e = self.days[d]
            # exposure: integrate over pieces in u = log(1 + age/age0) (age = t - origin), or in t if no origin
            R = e["R"]
            k, kind, t0, t1, orig = R[:, 0].astype(int), R[:, 1].astype(int), R[:, 2], R[:, 3], R[:, 4]
            F = self._features(e, k, kind, R[:, 5], R[:, 6], R[:, 7])
            has = ~np.isnan(orig)
            first = (kind == 0) & ~has
            F[first, self.ix["F_first"]] = 1
            ceil = np.where(kind == 0, self.A, self.D)
            base = F @ theta
            age_par = np.where(kind == 0, theta[self.ix["F_age"]], theta[self.ix["T_age"]])
            age_ix = np.where(kind == 0, self.ix["F_age"], self.ix["T_age"])
            # quadrature nodes:
            o = np.where(has, orig, 0.0)
            u0 = np.where(has, np.log1p((t0 - o) / self.age0), t0)
            u1 = np.where(has, np.log1p((t1 - o) / self.age0), t1)
            half, mid = (u1 - u0) / 2, (u1 + u0) / 2
            for x, w in zip(GX, GW):
                u = mid + half * x
                jac = np.where(has, self.age0 * np.exp(u), 1.0) * half * w  # dt = age0 e^u du
                lage = np.where(has, u, 0.0)
                Bv = base + age_par * lage
                s = expit(4 * Bv)
                rate = 2 * ceil * s
                L -= np.sum(rate * jac)
                gfac = 2 * ceil * 4 * s * (1 - s) * jac
                G -= F.T @ gfac
                np.add.at(G, age_ix, -gfac * lage)
            # event terms:
            Ev = e["E"]
            k, kind, t, orig = Ev[:, 0].astype(int), Ev[:, 1].astype(int), Ev[:, 2], Ev[:, 3]
            F = self._features(e, k, kind, Ev[:, 4], Ev[:, 5], Ev[:, 6])
            has = ~np.isnan(orig)
            F[(kind == 0) & ~has, self.ix["F_first"]] = 1
            lage = np.where(has, np.log1p((t - np.where(has, orig, 0)) / self.age0), 0.0)
            F[kind == 0, self.ix["F_age"]] = lage[kind == 0]
            F[kind == 1, self.ix["T_age"]] = lage[kind == 1]
            Bv = F @ theta
            ceil = np.where(kind == 0, self.A, self.D)
            L += np.sum(np.log(ceil) + np.log(expit(4 * Bv)))
            G += F.T @ (4 * (1 - expit(4 * Bv)))
        return L, G


    def residuals(self, theta, day):
        """Integrated pair termination hazard of each completed contact (time-rescaling residuals)."""
        e = self.days[day]
        R = e["R"]
        k, kind, t0, t1, orig = R[:, 0].astype(int), R[:, 1].astype(int), R[:, 2], R[:, 3], R[:, 4]
        m = kind == 1
        F = self._features(e, k[m], kind[m], R[m, 5], R[m, 6], R[m, 7])
        base = F @ theta
        a = theta[self.ix["T_age"]]
        o = orig[m]
        u0, u1 = np.log1p((t0[m] - o) / self.age0), np.log1p((t1[m] - o) / self.age0)
        half, mid = (u1 - u0) / 2, (u1 + u0) / 2
        comp = np.zeros(m.sum())
        for x, w in zip(GX, GW):
            u = mid + half * x
            comp += 2 * self.D * expit(4 * (base + a * u)) * self.age0 * np.exp(u) * half * w
        # sum the pieces of each contact; only contacts that ended (an event at the end of their last piece):
        ends = {(int(kk), float(tt)) for kk, kd, tt, oo in e["E"][:, :4] if kd == 1}
        key = {}
        for i, (kk, oo, tt) in enumerate(zip(k[m], o, t1[m])):
            key.setdefault((int(kk), float(oo)), [0.0, False])
            key[(int(kk), float(oo))][0] += comp[i]
            if (int(kk), float(tt)) in ends:
                key[(int(kk), float(oo))][1] = True
        return np.array([v for v, ended in key.values() if ended])

    def simulate(self, theta, day, rng):
        """Exact simulation of C1 on a day (pairs are independent): returns a list of (t, event, pair index)."""
        e = self.days[day]
        T = 300.0 if day == 2 else 540.0
        start_h = e["start_h"]
        X = e["X"]
        n = X.shape[0]
        th = theta
        ix = self.ix

        def B(kind, t, origin, linked_before):
            h = start_h + np.floor(t / 60).astype(int)
            mm = t - (h - start_h) * 60
            out = np.where(kind == 0, th[ix["F0"]], th[ix["T0"]]).astype(float)
            for i, hh in enumerate(models.HOURS[1:]):
                hv = int(hh[1:])
                out += np.where((h == hv) & (kind == 0), th[ix["F_" + hh]], 0) + np.where((h == hv) & (kind == 1), th[ix["T_" + hh]], 0)
            chg, brk = mm < 20, (h == 10) & (mm < 20)
            out += np.where(kind == 0, th[ix["F_chg"]] * chg + th[ix["F_brk"]] * brk, th[ix["T_chg"]] * chg + th[ix["T_brk"]] * brk)
            return out, h, mm
        static_F = X @ np.array([th[ix["F_" + s]] for s in STATIC])
        static_T = X @ np.array([th[ix["T_" + s]] for s in STATIC])
        t = np.zeros(n); linked = np.zeros(n, dtype=bool); origin = np.full(n, np.nan)
        active = np.ones(n, dtype=bool)
        events = []

        def rate(t_, kind_, origin_, idx):
            b, h, mm = B(kind_, t_, origin_, None)
            b = b + np.where(kind_ == 0, static_F[idx], static_T[idx])
            has = ~np.isnan(origin_)
            age = np.where(has, np.log1p((t_ - np.where(has, origin_, 0)) / self.age0), 0.0)
            b = b + np.where(kind_ == 0, np.where(has, th[ix["F_age"]] * age, th[ix["F_first"]]), th[ix["T_age"]] * age)
            return 2 * np.where(kind_ == 0, self.A, self.D) * expit(4 * b)
        while active.any():
            idx = np.where(active)[0]
            tt, kind, oo = t[idx], linked[idx].astype(int), origin[idx]
            # end of the current timetable piece (20-minute grid):
            piece_end = np.minimum((np.floor(tt / 20 + 1e-12) + 1) * 20, T)
            eps = 1e-9
            r0 = rate(tt, kind, oo, idx)
            r1 = rate(np.maximum(piece_end - eps, tt), kind, oo, idx)
            bound = np.maximum(np.maximum(r0, r1), 1e-300)  # (rates are monotone in the age within a piece)
            dt = rng.exponential(1.0, len(idx)) / bound
            tn = tt + dt
            cross = tn >= piece_end
            t[idx[cross]] = piece_end[cross]
            fin = idx[cross & (piece_end >= T)]
            active[fin] = False
            prop = ~cross
            if prop.any():
                ip = idx[prop]
                rn = rate(tn[prop], kind[prop], oo[prop], ip)
                acc = rng.random(prop.sum()) * bound[prop] < rn
                t[ip] = tn[prop]
                for j in np.where(acc)[0]:
                    pj = ip[j]
                    events.append((tn[prop][j], "terminate" if linked[pj] else "establish", pj))
                    linked[pj] = not linked[pj]
                    origin[pj] = tn[prop][j]
        events.sort()
        return events


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--train", type=int, nargs="+", default=[2, 3, 4])
    ap.add_argument("--test", type=int, nargs="+", default=[5, 6])
    ap.add_argument("--age0", type=float, nargs="+", default=[0.33, 1.0, 3.0])
    ap.add_argument("--tag", default="")
    args = ap.parse_args()
    best = None
    for a0 in args.age0:
        m = C1(args.train + args.test, a0)
        start = os.path.join(HERE, "results", "fit_M1%s.json" % args.tag)
        x0 = np.zeros(len(m.names))
        if os.path.exists(start):
            est = json.load(open(start))["estimates"]
            for n, v in est.items():
                key = n.replace("svy", "svy")
                if key in m.ix:
                    x0[m.ix[key]] = v
        f = lambda th: tuple(-v for v in m.loglik(th, args.train))
        r = minimize(f, x0, jac=True, method="L-BFGS-B", options={"maxiter": 1000})
        Ltr, g = m.loglik(r.x, args.train)
        Lte = {d: m.loglik(r.x, [d])[0] for d in args.test}
        print("age0 %.2f: logl train %.1f, test %.1f, converged %s, |g| %.2g" % (a0, Ltr, sum(Lte.values()), r.success, np.max(np.abs(g))))
        # standard errors from finite differences of the analytic gradient:
        k = len(r.x); H = np.zeros((k, k))
        for i in range(k):
            h = 1e-4 * max(1.0, abs(r.x[i])); xp, xm = r.x.copy(), r.x.copy(); xp[i] += h; xm[i] -= h
            H[i] = (m.loglik(xp, args.train)[1] - m.loglik(xm, args.train)[1]) / (2 * h)
        H = (H + H.T) / 2
        try:
            se = np.sqrt(np.maximum(np.diag(np.linalg.inv(-H)), 0))
        except np.linalg.LinAlgError:
            se = np.full(k, np.nan)
        res = {"model": "C1", "age0": a0, "logl_train": Ltr, "logl_test": sum(Lte.values()), "test": {d: {"logl": v} for d, v in Lte.items()},
               "estimates": dict(zip(m.names, map(float, r.x))), "standard_errors": dict(zip(m.names, map(float, se))),
               "converged": bool(r.success), "n_params": k, "gradient_norm": float(np.max(np.abs(g))),
               "events_train": sum(m.days[d]["n_events"] for d in args.train), "events_test": sum(m.days[d]["n_events"] for d in args.test)}
        if best is None or Ltr > best["logl_train"]:
            best = res
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    json.dump(best, open(os.path.join(HERE, "results", "fit_C1%s.json" % args.tag), "w"), indent=1)
    print("best age0 %.2f (by training log-likelihood): test %.1f" % (best["age0"], best["logl_test"]))


if __name__ == "__main__":
    main()

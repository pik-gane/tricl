#!/usr/bin/env python3
"""Unobserved heterogeneity: latent affinity classes of pairs for the pair-local models M6 and C1.

In M6 (timetable, static context, met before, decaying impacts of the pair's own contact starts and ends on four time
scales) and C1 (the same with power-law aging instead of impacts), the rates of a pair depend only on its own history,
the static context and the timetable, so the likelihood factorizes over pairs. With K latent classes per pair, which
add offsets a_k to the formation and b_k to the termination probability units (a_0 = b_0 = 0) and have prior
probabilities pi_k, the marginal log-likelihood is exact:

    L = sum_pairs log sum_k pi_k exp(l_pk),

with l_pk the log-likelihood of the pair's events under class k. A pair keeps its class on all days, so the held-out
log-likelihood uses each pair's posterior class probabilities given the training days.

The engine reproduces tricl's log-likelihood of M6 (with the impact threshold of the action types) and the C1 code for
K = 1 (see --validate). Rates are per pair and direction, as in tricl: an event contributes log(rate per direction),
the exposure is the pair rate (twice the rate per direction).

Usage: python3 study/hetero.py MODEL K [--train 2 3 4] [--test 5 6]     (MODEL = M6 or C1)
       python3 study/hetero.py --validate
"""
import argparse
import csv
import json
import os
import sys

import numpy as np
from scipy.optimize import minimize
from scipy.special import expit, logsumexp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import models  # noqa: E402

DATA = os.path.join(HERE, os.environ.get("STUDY_DATA", "data"))
GX, GW = np.polynomial.legendre.leggauss(8)
STATIC = ["cls", "fb", "svy", "mb"]
HOURS = models.HOURS[1:]
DECAYS = models.DECAYS4          # {"d2": 2.0, "d05": 0.5, "d01": 0.1, "d002": 0.02}
DNAMES = list(DECAYS)
DRATES = np.array([DECAYS[d] for d in DNAMES])
THETA_IMPACT = 0.01              # impact threshold of the action types (as in models.py)
AGE0 = 3.0                       # C1 age scale (the best of {1/3, 1, 3} minutes)


def base_names(model):
    n = []
    for e in ("F", "T"):
        n += [e + "0"] + ["%s_%s" % (e, h) for h in HOURS] + ["%s_chg" % e, "%s_brk" % e] + ["%s_%s" % (e, s) for s in STATIC]
    if model == "M6":
        n += ["F_end_" + d for d in DNAMES] + ["T_start_" + d for d in DNAMES]
    else:
        n += ["F_first", "F_age", "T_age"]
    return n


def slot_of(t, start_h):
    h = start_h + int(t // 60)
    m = t - (h - start_h) * 60
    return h, m < 20, (h == 10) and (m < 20)


def load_day(day, data):
    """Rows (pieces of constant regime) of the pairs with events, and the pairs without events."""
    T = 300.0 if day == 2 else 540.0
    start_h = 12 if day == 2 else 8
    students = data["students"]
    n = len(students)
    cls = dict(data["static"]["class"])
    fb = {tuple(sorted(p)) for p in data["static"]["facebook"]}
    svy = {}
    for a, b in data["static"]["survey"]:
        k = tuple(sorted((a, b))); svy[k] = svy.get(k, 0) + 1
    mb = set()
    with open(os.path.join(DATA, "forced_%d.csv" % day)) as fh:
        for r in csv.DictReader(fh):
            if r["relationship"] == "met before":
                mb.add(tuple(sorted((r["source"], r["target"]))))
    pairs = [(students[i], students[j]) for i in range(n) for j in range(i + 1, n)]
    pid = {p: k for k, p in enumerate(pairs)}
    X = np.array([[cls[a] == cls[b], (a, b) in fb, svy.get((a, b), 0), (a, b) in mb] for a, b in pairs], float)
    ev = {}
    with open(os.path.join(DATA, "events_%d.csv" % day)) as fh:
        for r in csv.DictReader(fh):
            ev.setdefault(pid[tuple(sorted((r["source"], r["target"])))], []).append((float(r["t"]), r["event"]))
    bounds = np.array(sorted(set([(h - start_h) * 60.0 + o for h in range(start_h, 17) for o in (0.0, 20.0)] + [T])))
    rows, evrows = [], []
    for p, lst in ev.items():
        # impacts (value at reference time, reference time) of the starts and ends of the pair's contacts, per decay rate:
        Is, ts_ = np.zeros(len(DRATES)), np.zeros(len(DRATES))
        Ie, te_ = np.zeros(len(DRATES)), np.zeros(len(DRATES))
        linked, origin, t_prev = False, np.nan, 0.0

        def imp(I, tr, t):
            v = I * np.exp(-DRATES * (t - tr))
            return np.where(v >= THETA_IMPACT, v, 0.0)  # (an action ceases to exist below the threshold)

        def emit(t0, t1):
            kind = 1 if linked else 0
            I, tr = (Is, ts_) if linked else (Ie, te_)
            # cut points: timetable boundaries and expiries of the relevant impacts
            with np.errstate(divide="ignore"):
                texp = np.where(I > 0, tr + np.log(np.maximum(I, 1e-300) / THETA_IMPACT) / DRATES, -np.inf)
            cuts = sorted(set([t0, t1] + [b for b in bounds if t0 < b < t1] + [x for x in texp if t0 < x < t1]))
            for a, b in zip(cuts[:-1], cuts[1:]):
                if b <= a:
                    continue
                h, chg, brk = slot_of((a + b) / 2, start_h)
                # impact values at a, but zero for impacts that expire before the middle of the piece:
                v = I * np.exp(-DRATES * (a - tr))
                v = np.where((I > 0) & (texp > (a + b) / 2), v, 0.0)
                rows.append((p, kind, a, b, h, chg, brk, origin) + tuple(v))
        for t, e in lst:
            emit(t_prev, t)
            kind = 1 if linked else 0
            h, chg, brk = slot_of(t, start_h)
            vI = imp(Is, ts_, t) if linked else imp(Ie, te_, t)
            evrows.append((p, kind, t, h, chg, brk, origin) + tuple(vI))
            if e == "establish":
                Is = imp(Is, ts_, t) + 1.0; ts_ = np.full(len(DRATES), t)
            else:
                Ie = imp(Ie, te_, t) + 1.0; te_ = np.full(len(DRATES), t)
            linked = not linked
            origin, t_prev = t, t
        emit(t_prev, T)
    R = np.array(rows, float)
    E = np.array(evrows, float)
    quiet = np.array(sorted(set(range(len(pairs))) - set(ev)), int)
    # slot schedule of the day for the pairs without events (constant rates within slots):
    edges = np.concatenate([[0.0], bounds])
    sched = [(b - a,) + slot_of((a + b) / 2, start_h) for a, b in zip(edges[:-1], edges[1:]) if b > a]
    return {"X": X, "R": R, "E": E, "quiet": quiet, "sched": np.array(sched, float), "n_pairs": len(pairs), "n_events": len(E)}


class Hetero:
    def __init__(self, model, K, days):
        self.model, self.K = model, K
        self.data = json.load(open(os.path.join(DATA, "data.json")))
        self.days = {d: load_day(d, self.data) for d in days}
        self.base = base_names(model)
        self.names = self.base + ["a_%d" % k for k in range(1, K)] + ["b_%d" % k for k in range(1, K)] + ["w_%d" % k for k in range(1, K)]
        self.ix = {n: i for i, n in enumerate(self.names)}
        self.nb = len(self.base)

    # --- design ---------------------------------------------------------------------------------------------
    def _design(self, X, kind, h, chg, brk):
        """Static and timetable features (rows x base params), kind-specific."""
        F = np.zeros((len(kind), self.nb))
        for kd, pre in ((0, "F"), (1, "T")):
            m = kind == kd
            if not m.any():
                continue
            F[m, self.ix[pre + "0"]] = 1
            for hh in HOURS:
                F[m & (h == int(hh[1:])), self.ix["%s_%s" % (pre, hh)]] = 1
            F[m & (chg > 0), self.ix[pre + "_chg"]] = 1
            F[m & (brk > 0), self.ix[pre + "_brk"]] = 1
            for i, s in enumerate(STATIC):
                F[m, self.ix["%s_%s" % (pre, s)]] = X[m, i]
        return F

    def _dyn_index(self, kind):
        """Indices of the dynamic parameters (impacts or aging) for each kind."""
        if self.model == "M6":
            return [self.ix["F_end_" + d] for d in DNAMES], [self.ix["T_start_" + d] for d in DNAMES]
        return None

    # --- per-day computation -----------------------------------------------------------------------------------
    def day_terms(self, th, day, W=None):
        """If W is None: matrix l[p, k] of the pair log-likelihoods (pairs of the day x classes).
        Otherwise: gradient sum_pk W[p, k] d l_pk / d theta (without the mixture weights)."""
        D = self.days[day]
        K = self.K
        a = np.concatenate([[0.0], th[[self.ix["a_%d" % k] for k in range(1, K)]]]) if K > 1 else np.zeros(1)
        b = np.concatenate([[0.0], th[[self.ix["b_%d" % k] for k in range(1, K)]]]) if K > 1 else np.zeros(1)
        P = D["n_pairs"]
        L = np.zeros((P, K)) if W is None else None
        G = np.zeros(len(th)) if W is not None else None
        tb = th[:self.nb]
        ceilF, ceilT = models.A_FORM, models.D_TERM
        # pairs without events: only formation exposure, constant within slots
        q = D["quiet"]
        if len(q):
            Xq = D["X"][q]
            for (length, h, chg, brk) in D["sched"]:
                kind = np.zeros(len(q), int)
                F = self._design(Xq, kind, np.full(len(q), h), np.full(len(q), chg), np.full(len(q), brk))
                if self.model == "C1":
                    F[:, self.ix["F_first"]] = 1
                Bq = F @ tb
                for k in range(K):
                    s = expit(4 * (Bq + a[k]))
                    if W is None:
                        L[q, k] -= length * 2 * ceilF * s
                    else:
                        g = -length * 2 * ceilF * 4 * s * (1 - s) * W[q, k]
                        G[:self.nb] += F.T @ g
                        if k > 0:
                            G[self.ix["a_%d" % k]] += g.sum()
        # pieces of pairs with events
        R = D["R"]
        p, kind, t0, t1 = R[:, 0].astype(int), R[:, 1].astype(int), R[:, 2], R[:, 3]
        F = self._design(D["X"][p], kind, R[:, 4], R[:, 5], R[:, 6])
        origin = R[:, 7]
        I0 = R[:, 8:8 + len(DRATES)]
        base = F @ tb
        off = np.where(kind[:, None] == 0, a[None, :], b[None, :])  # rows x K
        ceil = np.where(kind == 0, ceilF, ceilT)
        if self.model == "M6":
            fi, ti = self._dyn_index(None)
            coef = np.where(kind[:, None] == 0, th[fi][None, :], th[ti][None, :])  # rows x 4
            active = (I0 > 0)
            dmax = np.where(active.any(axis=1), (active * DRATES[None, :]).max(axis=1), 0.0)
            # constant rows:
            c = dmax == 0
            # rows with decaying impacts: split into sub-pieces of length <= 1/dmax
            nsub = np.where(c, 0, np.ceil((t1 - t0) * np.maximum(dmax, 1e-12)).astype(int).clip(1, 400))
        for k in range(K):
            if self.model == "M6":
                # constant rows
                s = expit(4 * (base[c] + off[c, k]))
                lam = 2 * ceil[c] * s
                if W is None:
                    np.subtract.at(L[:, k], p[c], lam * (t1[c] - t0[c]))
                else:
                    g = -(t1[c] - t0[c]) * 2 * ceil[c] * 4 * s * (1 - s) * W[p[c], k]
                    G[:self.nb] += F[c].T @ g
                    self._off_grad(G, g, kind[c], k)
                # decaying rows
                idx = np.where(~c)[0]
                if len(idx):
                    rep = np.repeat(idx, nsub[idx])
                    j = np.concatenate([np.arange(m) for m in nsub[idx]])
                    seg = (t1[rep] - t0[rep]) / nsub[rep]
                    sa = j * seg  # offset of the sub-piece from t0
                    for x, w in zip(GX, GW):
                        s_t = sa + seg * (x + 1) / 2
                        imp = I0[rep] * np.exp(-DRATES[None, :] * s_t[:, None])
                        dyn = (coef[rep] * imp).sum(axis=1)
                        sv = expit(4 * (base[rep] + dyn + off[rep, k]))
                        wq = w * seg / 2
                        if W is None:
                            np.subtract.at(L[:, k], p[rep], wq * 2 * ceil[rep] * sv)
                        else:
                            g = -wq * 2 * ceil[rep] * 4 * sv * (1 - sv) * W[p[rep], k]
                            G[:self.nb] += F[rep].T @ g
                            dd = (imp * g[:, None])
                            k0 = kind[rep] == 0
                            G[fi] += dd[k0].sum(axis=0); G[ti] += dd[~k0].sum(axis=0)
                            self._off_grad(G, g, kind[rep], k)
            else:  # C1
                has = ~np.isnan(origin)
                Fx = F.copy()
                Fx[(kind == 0) & ~has, self.ix["F_first"]] = 1
                base1 = Fx @ tb
                agep = np.where(kind == 0, th[self.ix["F_age"]], th[self.ix["T_age"]])
                agei = np.where(kind == 0, self.ix["F_age"], self.ix["T_age"])
                o = np.where(has, origin, 0.0)
                u0 = np.where(has, np.log1p((t0 - o) / AGE0), t0)
                u1 = np.where(has, np.log1p((t1 - o) / AGE0), t1)
                half, mid = (u1 - u0) / 2, (u1 + u0) / 2
                for x, w in zip(GX, GW):
                    u = mid + half * x
                    jac = np.where(has, AGE0 * np.exp(u), 1.0) * half * w
                    la = np.where(has, u, 0.0)
                    sv = expit(4 * (base1 + agep * la + off[:, k]))
                    if W is None:
                        np.subtract.at(L[:, k], p, jac * 2 * ceil * sv)
                    else:
                        g = -jac * 2 * ceil * 4 * sv * (1 - sv) * W[p, k]
                        G[:self.nb] += Fx.T @ g
                        np.add.at(G, agei, g * la)
                        self._off_grad(G, g, kind, k)
        # event terms
        E = D["E"]
        pe, ke, te = E[:, 0].astype(int), E[:, 1].astype(int), E[:, 2]
        Fe = self._design(D["X"][pe], ke, E[:, 3], E[:, 4], E[:, 5])
        ce = np.where(ke == 0, ceilF, ceilT)
        if self.model == "M6":
            fi, ti = self._dyn_index(None)
            Ie = E[:, 7:7 + len(DRATES)]
            coefe = np.where(ke[:, None] == 0, th[fi][None, :], th[ti][None, :])
            Be = Fe @ tb + (coefe * Ie).sum(axis=1)
        else:
            has = ~np.isnan(E[:, 6])
            Fe[(ke == 0) & ~has, self.ix["F_first"]] = 1
            la = np.where(has, np.log1p((te - np.where(has, E[:, 6], 0)) / AGE0), 0.0)
            Be = Fe @ tb + np.where(ke == 0, th[self.ix["F_age"]], th[self.ix["T_age"]]) * la
        for k in range(K):
            Bk = Be + np.where(ke == 0, a[k], b[k])
            if W is None:
                np.add.at(L[:, k], pe, np.log(ce) + np.log(expit(4 * Bk)))
            else:
                g = 4 * (1 - expit(4 * Bk)) * W[pe, k]
                G[:self.nb] += Fe.T @ g
                if self.model == "M6":
                    dd = Ie * g[:, None]
                    G[fi] += dd[ke == 0].sum(axis=0); G[ti] += dd[ke == 1].sum(axis=0)
                else:
                    G[self.ix["F_age"]] += (g * la)[ke == 0].sum(); G[self.ix["T_age"]] += (g * la)[ke == 1].sum()
                self._off_grad(G, g, ke, k)
        return L if W is None else G

    def _off_grad(self, G, g, kind, k):
        if k == 0:
            return
        G[self.ix["a_%d" % k]] += g[kind == 0].sum()
        G[self.ix["b_%d" % k]] += g[kind == 1].sum()

    def logpi(self, th):
        w = np.concatenate([[0.0], th[[self.ix["w_%d" % k] for k in range(1, self.K)]]]) if self.K > 1 else np.zeros(1)
        return w - logsumexp(w)

    def pair_loglik(self, th, days):
        """Summed pair log-likelihoods l[p, k] over days (a pair keeps its class on all days)."""
        return sum(self.day_terms(th, d) for d in days)

    def marginal(self, th, days, grad=False):
        lp = self.logpi(th)
        Lpk = self.pair_loglik(th, days) + lp[None, :]
        tot = logsumexp(Lpk, axis=1)
        L = tot.sum()
        if not grad:
            return L
        R = np.exp(Lpk - tot[:, None])  # responsibilities
        G = sum(self.day_terms(th, d, W=R) for d in days)
        pi = np.exp(lp)
        for k in range(1, self.K):
            G[self.ix["w_%d" % k]] = (R[:, k] - pi[k]).sum()
        return L, G

    def predictive(self, th, train, test):
        """Held-out log-likelihood with the posterior class probabilities of each pair given the training days."""
        lp = self.logpi(th)
        Ltr = self.pair_loglik(th, train) + lp[None, :]
        post = Ltr - logsumexp(Ltr, axis=1)[:, None]
        Lte = self.pair_loglik(th, test)
        return logsumexp(post + Lte, axis=1).sum(), logsumexp(lp[None, :] + Lte, axis=1).sum()


def start_values(h, model, K):
    th = np.zeros(len(h.names))
    src = os.path.join(HERE, "results", "fit_%s.json" % model)
    if os.path.exists(src):
        est = json.load(open(src))["estimates"]
        for n, v in est.items():
            if n in h.ix:
                th[h.ix[n]] = v
    for k in range(1, K):  # spread the classes: higher affinity for higher k
        th[h.ix["a_%d" % k]] = 0.3 * k
        th[h.ix["b_%d" % k]] = -0.1 * k
        th[h.ix["w_%d" % k]] = -1.0 * k
    return th


def fit(model, K, train, test):
    h = Hetero(model, K, sorted(set(train + test)))
    th0 = start_values(h, model, K)
    f = lambda th: tuple(-v for v in h.marginal(th, train, grad=True))
    r = minimize(f, th0, jac=True, method="L-BFGS-B", options={"maxiter": 2000})
    Ltr = h.marginal(r.x, train)
    Lte, Lte_prior = h.predictive(r.x, train, test)
    n_tr = sum(h.days[d]["n_events"] for d in train)
    n_te = sum(h.days[d]["n_events"] for d in test)
    out = {"model": "%s+K%d" % (model, K), "K": K, "logl_train": float(Ltr), "logl_test": float(Lte), "logl_test_prior": float(Lte_prior),
           "events_train": n_tr, "events_test": n_te, "n_params": len(h.names), "converged": bool(r.success), "message": str(r.message),
           "estimates": dict(zip(h.names, map(float, r.x))), "pi": list(map(float, np.exp(h.logpi(r.x))))}
    json.dump(out, open(os.path.join(HERE, "results", "fit_%s_K%d.json" % (model, K)), "w"), indent=1)
    print("%s K=%d: train %.1f, held-out %.1f (posterior classes) / %.1f (prior), %d params, pi %s, converged %s"
          % (model, K, Ltr, Lte, Lte_prior, len(h.names), np.round(out["pi"], 3), r.success), flush=True)
    return out


def validate():
    """K = 1 must reproduce tricl (M6) and the C1 code at their estimates."""
    for model in ("M6", "C1"):
        est = json.load(open(os.path.join(HERE, "results", "fit_%s.json" % model)))["estimates"]
        h = Hetero(model, 1, [3])
        th = np.array([est.get(n, 0.0) for n in h.names])
        L = h.marginal(th, [3])
        if model == "M6":
            import fit as fitmod
            ref = fitmod.run(os.path.join(HERE, "configs", "M6.yaml"), 3, est)["logl"]
        else:
            import c1
            m = c1.C1([3], AGE0)
            ref = m.loglik(np.array([est[k] for k in m.names]), [3])[0]
        print("%s day 3: engine %.6f, reference %.6f, difference %.2e" % (model, L, ref, L - ref))
        # gradient check by finite differences in a few coordinates:
        Lg, G = h.marginal(th, [3], grad=True)
        for n in h.names[:2] + h.names[-3:]:
            i = h.ix[n]; e = 1e-5
            tp, tm = th.copy(), th.copy(); tp[i] += e; tm[i] -= e
            fd = (h.marginal(tp, [3]) - h.marginal(tm, [3])) / (2 * e)
            print("   d/d%-12s analytic %.5f finite difference %.5f" % (n, G[i], fd))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("model", nargs="?")
    ap.add_argument("K", nargs="?", type=int)
    ap.add_argument("--train", type=int, nargs="+", default=[2, 3, 4])
    ap.add_argument("--test", type=int, nargs="+", default=[5, 6])
    ap.add_argument("--validate", action="store_true")
    args = ap.parse_args()
    if args.validate:
        validate()
    else:
        fit(args.model, args.K, args.train, args.test)

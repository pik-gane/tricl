#!/usr/bin/env python3
"""Macroscopic (mean-field) approximation of a tricl model, generated automatically from its config file.

Usage:
    python3 python/tricl_macro.py CONFIG [--bin PATH] [--t-max T] [--dt DT] [--closure poisson|mean]
                                  [--seeds 1 2 ...] [--out FILE.csv] [--inf-rate R] [--meta NAME=VALUE ...]

The model structure is read from `tricl CONFIG --dump-model`. The state variables of the approximation are the
numbers n_L of links of every link type L = (source type, relationship, target type). The rate of an event
(establishment or termination of a link of type L) depends on the numbers of adjacent angles of the influencing
angle types; in the approximation, the number of angles of type (r12, et2, r23) adjacent to a dyad is a Poisson
random variable with mean

    N_et2 * rho(et1, r12, et2) * rho(et2, r23, et3),

independent across angle types, where rho(a, r, b) = n_(a,r,b) / (number of ordered pairs of types a and b) is the
current link density (angles through the identity relationship "=" are links of the same dyad, with mean rho of
the respective link type). Event rates are averaged over these distributions ("poisson" closure, default) or
evaluated at the means ("mean" closure). This yields, per link type,

    dn_L/dt = (number of unlinked ordered pairs) * <establishment rate> - n_L * <termination rate>,

where events of symmetric relationships and of relationships with named inverses change both directed links.

With --seeds, the model is also simulated with tricl for these seeds (writing link counts every DT model time
units via --stats-out), and the mean and standard deviation over seeds are printed next to the approximation.
The approximation neglects all correlations between links (e.g. clustering of the initial network and the
correlations created by triadic closure), so deviations from the simulation measure exactly these effects.
Infinite attempt rates ("immediate" events) are replaced by --inf-rate (default 1e4) with a warning.

Only the Python standard library is needed (scipy is used for the integration if available).
"""

import argparse
import json
import math
import os
import re
import subprocess
import sys
import tempfile


# ------------------------------------------------------------------------------------------------------------
# the sigmoid of tricl (see src/probability.h)

def softplus(x):
    return x + math.log1p(math.exp(-x)) if x > 0 else math.log1p(math.exp(x))


def tail_term(tail, x):
    """T_tail(x) = (1 + tail ln(1 + e^x))^(-1/tail), or 1/(1 + e^x) for tail 0."""
    y = softplus(x)
    return math.exp(-y) if tail == 0 else math.exp(-math.log1p(tail * y) / tail)


def tail_term_complement(tail, x):
    y = softplus(x)
    return -math.expm1(-y) if tail == 0 else -math.expm1(-math.log1p(tail * y) / tail)


def tail2scale(tail):
    """s(tail) = k(tail) / 4 with k(tail) = (1 + tail ln 2)^(-1 - 1/tail), k(0) = 1/2 (normalises the slope at 0 to 1)."""
    return 0.125 if tail == 0 else math.exp(-(1 + 1 / tail) * math.log1p(tail * math.log(2.0))) / 4


def sigmoid(pu, left_tail, right_tail):
    """The sigmoidal function of tricl: T_left(-v)/2 + 1/2 - T_right(v)/2 with v = pu / (s(left) + s(right)); slope 1 at 0."""
    if pu == math.inf:
        return 1.0
    if pu == -math.inf:
        return 0.0
    if left_tail == 0 and right_tail == 0:
        return 1 / (1 + math.exp(-4 * pu)) if pu > -175 else 0.0
    v = pu / (tail2scale(left_tail) + tail2scale(right_tail))
    return (tail_term(left_tail, -v) + tail_term_complement(right_tail, v)) / 2


# ------------------------------------------------------------------------------------------------------------
# reading the model

def with_time_limit(config, t_max):
    """Return the text of a config file with limits:t set to t_max and no limit on the number of events."""
    with open(config) as f:
        lines = f.read().split("\n")
    out, i, n, done = [], 0, len(lines), False
    while i < n:
        line = lines[i]
        out.append(line)
        i += 1
        if re.match(r"^limits\s*:", line):
            indent = "    "
            while i < n and (lines[i].strip() == "" or lines[i].lstrip().startswith("#") or lines[i][0] in " \t"):
                m = re.match(r"^(\s+)\S", lines[i])
                if m:
                    indent = m.group(1)
                if not re.match(r"^\s+(t|events)\s*:", lines[i]):
                    out.append(lines[i])
                i += 1
            out.append(indent + "t: " + repr(t_max))
            out.append(indent + "events: inf")
            done = True
    if not done:
        out += ["limits:", "    t: " + repr(t_max), "    events: inf"]
    return "\n".join(out)


def meta_args(meta):
    args = []
    for name, value in meta.items():
        args += ["--" + name, str(value)]
    return args


def dump_model(binary, config, meta, seed):
    p = subprocess.run([binary, config, "--dump-model", "--seed", str(seed)] + meta_args(meta), capture_output=True, text=True)
    if p.returncode != 0:
        raise SystemExit("tricl --dump-model failed:\n" + p.stderr[-2000:])
    return json.loads(p.stdout.strip().split("\n")[-1])


# ------------------------------------------------------------------------------------------------------------
# the approximation

class MacroModel:
    def __init__(self, model, closure="poisson", inf_rate=1e4, warn=print):
        self.N = model["entity_types"]
        self.inverse = {r: d["inverse"] for r, d in model["relationship_types"].items()}
        self.closure = closure
        self.inf_rate = inf_rate
        # state variables: one per link type (also for link types that only appear in event types):
        self.link_types = [(lt["source"], lt["relationship"], lt["target"]) for lt in model["link_types"]]
        self.n0 = {(lt["source"], lt["relationship"], lt["target"]): float(lt["links"]) for lt in model["link_types"]}
        for evt in model["event_types"]:
            L = (evt["source"], evt["relationship"], evt["target"])
            for key in (L, self.inv(L)):
                if key not in self.n0:
                    self.link_types.append(key)
                    self.n0[key] = 0.0
        self.index = {L: i for i, L in enumerate(self.link_types)}
        # event types:
        self.events = []
        self.stiff = False  # whether infinite rates were replaced by large finite ones
        # mutually exclusive link types on the same dyad: a rule "terminate r immediately when the link r23 of the
        # same dyad exists" (an infinite attempt rate via an angle through the identity relationship) means that
        # establishing r23 removes r. These rules are handled structurally instead of as (stiff) rates:
        self.kills = {}  # (et1, r23, et3) -> set of relationship labels r that are terminated when r23 is established
        for evt in model["event_types"]:
            e = dict(evt)
            e["L"] = (evt["source"], evt["relationship"], evt["target"])
            e["influences"] = [dict(f) for f in evt["influences"]]
            if e["class"] == "terminate":
                for f in e["influences"]:
                    if f["attempt"] == math.inf and f["probunits"] in (0, math.inf):
                        other = None
                        if f["rat12"] == "=" and f["et2"] == e["source"]:
                            other = (e["source"], f["rat23"], e["target"])
                        elif f["rat23"] == "=" and f["et2"] == e["target"]:
                            other = (e["source"], f["rat12"], e["target"])
                        if other is not None:
                            self.kills.setdefault(other, set()).add(e["relationship"])
                            f["attempt"] = 0.0  # handled structurally
            if e["base_attempt"] == math.inf:
                warn("warning: infinite base attempt rate of %s replaced by %g" % (self.evt_label(e), inf_rate))
                e["base_attempt"] = inf_rate
                self.stiff = True
            for f in e["influences"]:
                if f["attempt"] == math.inf:
                    warn("warning: infinite attempt rate of an influence on %s replaced by %g" % (self.evt_label(e), inf_rate))
                    f["attempt"] = inf_rate
                    self.stiff = True
            self.events.append(e)

    def exclusive(self, et1, r, r23, et3):
        """Whether links of types r and r23 never coexist on a dyad (one of them kills the other)."""
        return r in self.kills.get((et1, r23, et3), ()) or r23 in self.kills.get((et1, r, et3), ())

    def inv(self, L):
        et1, r, et3 = L
        return (et3, self.inverse.get(r) or r, et1) if self.inverse.get(r) else L

    @staticmethod
    def evt_label(e):
        return "%s %s %s %s" % (e["class"], e["source"], e["relationship"], e["target"])

    def pairs(self, et1, et3):
        """Number of ordered pairs of entities of the two types (excluding self pairs)."""
        return self.N[et1] * (self.N[et3] - 1) if et1 == et3 else self.N[et1] * self.N[et3]

    def rho(self, n, L):
        """Link density of link type L (0 if the model has no such link type)."""
        if L not in self.index:
            return 0.0
        p = self.pairs(L[0], L[2])
        return min(1.0, max(0.0, n[self.index[L]] / p)) if p > 0 else 0.0

    def angle_mean(self, n, e, infl):
        """Expected number of angles of an influence's type adjacent to a dyad of the event's types."""
        et1, r, et3 = e["L"]
        r12, et2, r23 = infl["rat12"], infl["et2"], infl["rat23"]
        if r12 == "=" or r23 == "=":
            # the angle is a link of another type on the same dyad:
            if r12 == "=":
                other = r23 if et2 == et1 else None
            else:
                other = r12 if et2 == et3 else None
            if other is None:
                return 0.0
            rho_other = self.rho(n, (et1, other, et3))
            if not self.exclusive(et1, r, other, et3):
                return rho_other
            if e["class"] == "terminate":
                return 0.0  # the dyad has r, so it cannot have the exclusive other link
            rho_r = self.rho(n, (et1, r, et3))
            return min(1.0, rho_other / (1 - rho_r)) if rho_r < 1 else 0.0
        n_mid = self.N[et2] - (1 if et2 == et1 else 0) - (1 if et2 == et3 else 0)
        return max(0, n_mid) * self.rho(n, (et1, r12, et2)) * self.rho(n, (et2, r23, et3))

    def rate(self, n, e):
        """Closure-averaged rate of an event type given the state n."""
        mus = [self.angle_mean(n, e, infl) for infl in e["influences"]]
        a0, b0 = e["base_attempt"], e["base_probunits"]
        lt, rt = e["left_tail"], e["right_tail"]
        infl = e["influences"]

        def times(value, count):  # (avoids inf * 0 = nan)
            return 0.0 if count == 0 else value * count

        if self.closure == "mean" or not infl:
            A = a0 + sum(times(f["attempt"], mu) for f, mu in zip(infl, mus))
            B = b0 + sum(times(f["probunits"], mu) for f, mu in zip(infl, mus))
            return A * sigmoid(B, lt, rt)
        # poisson closure: E[(a0 + sum a_j n_j) sigma(b0 + sum b_j n_j)] with independent Poisson n_j.
        # influences without effect on the probunits only enter through their means:
        b_idx = [j for j, f in enumerate(infl) if f["probunits"] != 0 and mus[j] > 0]
        A_rest = a0 + sum(times(infl[j]["attempt"], mus[j]) for j in range(len(infl)) if j not in b_idx)
        if not b_idx:
            return A_rest * sigmoid(b0, lt, rt)
        if len(b_idx) > 3:  # too many dimensions to enumerate
            A = a0 + sum(times(f["attempt"], mu) for f, mu in zip(infl, mus))
            B = b0 + sum(times(f["probunits"], mu) for f, mu in zip(infl, mus))
            return A * sigmoid(B, lt, rt)
        # enumerate the joint distribution of the counts that affect the probunits:
        supports = []
        for j in b_idx:
            mu = mus[j]
            kmax = int(mu + 8 * math.sqrt(mu) + 8)
            probs, p = [], math.exp(-mu)
            for k in range(kmax + 1):
                probs.append(p)
                p *= mu / (k + 1)
            supports.append((j, probs))
        total = 0.0

        def recurse(d, prob, A, B):
            nonlocal total
            if d == len(supports):
                total += prob * A * sigmoid(B, lt, rt)
                return
            j, probs = supports[d]
            for k, pk in enumerate(probs):
                if pk < 1e-15 and k > mus[j]:
                    break
                recurse(d + 1, prob * pk, A + times(infl[j]["attempt"], k), B + times(infl[j]["probunits"], k))

        recurse(0, 1.0, A_rest, b0)
        return total

    def derivative(self, n):
        dn = [0.0] * len(n)
        for e in self.events:
            L = e["L"]
            i = self.index[L]
            r = self.rate(n, e)
            if e["class"] == "establish":
                flow = max(0.0, self.pairs(L[0], L[2]) - n[i]) * r
            elif e["class"] == "terminate":
                flow = -max(0.0, n[i]) * r
            else:
                continue
            dn[i] += flow
            if self.inverse.get(L[1]) is not None:
                # the companion event changes the inverse link as well
                # (for symmetric relationships inv(L) == L: each event changes two directed links):
                dn[self.index[self.inv(L)]] += flow
            if e["class"] == "establish":
                # links of types that are terminated immediately when a link of type L is established:
                rho_L = self.rho(n, L)
                for r in self.kills.get(L, ()):
                    K = (L[0], r, L[2])
                    if K not in self.index:
                        continue
                    # fraction of dyads without L that have K (K and L never coexist):
                    frac = min(1.0, self.rho(n, K) / (1 - rho_L)) if rho_L < 1 else 0.0
                    dn[self.index[K]] -= flow * frac
                    if self.inverse.get(r) is not None and self.inv(K) != K:
                        dn[self.index[self.inv(K)]] -= flow * frac
        return dn

    def integrate(self, t_max, dt, n_steps_per_dt=20):
        """Integrate from the initial link counts; returns (times, states) at multiples of dt.

        Uses scipy's LSODA if available (which handles the stiffness caused by replaced infinite rates),
        otherwise classical Runge-Kutta, or implicit Euler with a numerical Jacobian if the system is stiff.
        """
        n0 = [self.n0[L] for L in self.link_types]
        grid = [k * dt for k in range(int(round(t_max / dt)) + 1)]
        try:
            from scipy.integrate import solve_ivp
            sol = solve_ivp(lambda t, y: self.derivative(list(y)), (0.0, grid[-1]), n0, t_eval=grid,
                            method="LSODA", rtol=1e-7, atol=1e-6)
            if sol.success:
                return list(sol.t), [list(sol.y[:, k]) for k in range(len(sol.t))]
        except ImportError:
            pass
        n = list(n0)
        times, states = [0.0], [list(n)]
        if self.stiff:
            h = dt / max(n_steps_per_dt, 100)
            for t in grid[1:]:
                while times[-1] < t - 1e-12:
                    n = self._implicit_euler_step(n, min(h, t - times[-1]))
                    times.append(min(times[-1] + h, t))
                states.append(list(n))
            return grid, states
        h = dt / n_steps_per_dt
        for t in grid[1:]:
            for _ in range(n_steps_per_dt):
                # classical Runge-Kutta step:
                k1 = self.derivative(n)
                k2 = self.derivative([x + h / 2 * k for x, k in zip(n, k1)])
                k3 = self.derivative([x + h / 2 * k for x, k in zip(n, k2)])
                k4 = self.derivative([x + h * k for x, k in zip(n, k3)])
                n = [x + h / 6 * (a + 2 * b + 2 * c + d) for x, a, b, c, d in zip(n, k1, k2, k3, k4)]
            states.append(list(n))
        return grid, states

    def _implicit_euler_step(self, n, h):
        """One implicit Euler step n_new = n + h f(n_new), solved by Newton's method with a numerical Jacobian."""
        m = len(n)
        x = list(n)
        for _ in range(4):
            f = self.derivative(x)
            g = [xi - ni - h * fi for xi, ni, fi in zip(x, n, f)]
            if max(abs(v) for v in g) < 1e-9 * max(1.0, max(abs(v) for v in x)):
                break
            # Jacobian of g = I - h J_f by forward differences:
            J = [[0.0] * m for _ in range(m)]
            for j in range(m):
                eps = 1e-6 * max(1.0, abs(x[j]))
                xp = list(x)
                xp[j] += eps
                fp = self.derivative(xp)
                for i in range(m):
                    J[i][j] = (1.0 if i == j else 0.0) - h * (fp[i] - f[i]) / eps
            # solve J d = g by Gaussian elimination with partial pivoting:
            A = [row[:] + [g[i]] for i, row in enumerate(J)]
            for c in range(m):
                p = max(range(c, m), key=lambda r: abs(A[r][c]))
                A[c], A[p] = A[p], A[c]
                if abs(A[c][c]) < 1e-300:
                    continue
                for r in range(c + 1, m):
                    fac = A[r][c] / A[c][c]
                    if fac:
                        for k in range(c, m + 1):
                            A[r][k] -= fac * A[c][k]
            d = [0.0] * m
            for c in range(m - 1, -1, -1):
                s = A[c][m] - sum(A[c][k] * d[k] for k in range(c + 1, m))
                d[c] = s / A[c][c] if abs(A[c][c]) >= 1e-300 else 0.0
            x = [xi - di for xi, di in zip(x, d)]
        return x


# ------------------------------------------------------------------------------------------------------------
# simulation for comparison

def simulate_stats(binary, config_text, seeds, t_max, dt, meta, workdir):
    """Run tricl for each seed and return {seed: {t: {link type label: count}}}."""
    config = os.path.join(workdir, "macro_config.yaml")
    with open(config, "w") as f:
        f.write(config_text)
    results = {}
    for seed in seeds:
        stats = os.path.join(workdir, "stats_%d.csv" % seed)
        p = subprocess.run([binary, config, "--summary", "--seed", str(seed), "--stats-out", stats, "--stats-every", repr(dt)] + meta_args(meta),
                           cwd=workdir, capture_output=True, text=True)
        if p.returncode != 0:
            raise SystemExit("simulation with seed %d failed:\n%s" % (seed, p.stderr[-2000:]))
        import csv
        with open(stats) as f:
            rows = list(csv.reader(f))
        header = rows[0]
        results[seed] = {}
        for row in rows[1:]:
            t = float(row[0])
            results[seed][round(t / dt)] = {header[k]: float(row[k]) for k in range(1, len(header))}
    return results


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("config")
    ap.add_argument("--bin", default=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "build", "src", "tricl"))
    ap.add_argument("--t-max", type=float, default=None, help="end of the time interval (default: limits:t of the config)")
    ap.add_argument("--dt", type=float, default=None, help="output time step (default: t_max / 20)")
    ap.add_argument("--closure", choices=["poisson", "mean"], default="poisson")
    ap.add_argument("--seeds", nargs="*", type=int, default=[], help="also simulate with these seeds and compare")
    ap.add_argument("--seed", type=int, default=1, help="seed for the random initial state of the approximation")
    ap.add_argument("--out", default="", help="write the results as csv (columns t, link type, ode, sim_mean, sim_sd)")
    ap.add_argument("--inf-rate", type=float, default=1e4)
    ap.add_argument("--meta", nargs="*", default=[], help="metaparameter overrides NAME=VALUE")
    args = ap.parse_args()

    binary = os.path.abspath(args.bin)
    meta = dict(item.split("=", 1) for item in args.meta)
    model = dump_model(binary, args.config, meta, args.seed)
    t_max = args.t_max if args.t_max is not None else model["t_max"]
    if not (t_max and math.isfinite(t_max)):
        raise SystemExit("please give --t-max (the config has no finite limits:t)")
    dt = args.dt or t_max / 20
    macro = MacroModel(model, closure=args.closure, inf_rate=args.inf_rate)
    times, states = macro.integrate(t_max, dt)
    labels = ["%s %s %s" % L for L in macro.link_types]

    sim = {}
    if args.seeds:
        workdir = tempfile.mkdtemp(prefix="tricl_macro_")
        config_text = with_time_limit(args.config, t_max)
        sim = simulate_stats(binary, config_text, args.seeds, t_max, dt, meta, workdir)

    # output:
    dynamic = [i for i, L in enumerate(macro.link_types) if any(e["L"] in (L, macro.inv(L)) for e in macro.events)]
    rows = []
    print("%8s  %-40s %12s %12s %12s" % ("t", "link type", "approx.", "sim. mean", "sim. sd"))
    for k, (t, n) in enumerate(zip(times, states)):
        for i in dynamic:
            label = labels[i]
            vals = [sim[s][k][label] for s in sim if k in sim[s] and label in sim[s][k]]
            mean = sum(vals) / len(vals) if vals else float("nan")
            sd = math.sqrt(sum((v - mean) ** 2 for v in vals) / (len(vals) - 1)) if len(vals) > 1 else float("nan")
            rows.append((t, label, n[i], mean, sd))
            if k % max(1, len(times) // 10) == 0 or k == len(times) - 1:
                print("%8.3f  %-40s %12.2f %12.2f %12.2f" % (t, label, n[i], mean, sd))
    if args.out:
        with open(args.out, "w") as f:
            f.write("t,link type,ode,sim_mean,sim_sd\n")
            for t, label, ode, mean, sd in rows:
                f.write("%r,%s,%r,%r,%r\n" % (t, '"' + label.replace('"', '""') + '"', ode, mean, sd))
        print("wrote", args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())

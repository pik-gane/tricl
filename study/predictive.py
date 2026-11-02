#!/usr/bin/env python3
"""Predictive checks on the held-out days: simulate each fitted model (tricl models with the timetable as forced
events, C1 with its own exact simulator) and compare statistics not used directly in fitting with the data:
contact durations, inter-contact times of pairs, numbers of contacts, and the sizes of groups (connected components
of the current contact network, as experienced by students). Also time-rescaling residuals of terminations.

Writes study/results/predictive.json and the raw samples used for the figures (study/results/samples_*.npz).
"""
import json
import os
import subprocess
import sys
import tempfile

import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fit  # noqa: E402
import c1  # noqa: E402

DAYS = [5, 6]
SEEDS = [1, 2, 3]
MODELS = ["M0", "M1", "M2", "M4", "M3", "C1"]


def read_events(path):
    import csv
    out = []
    with open(path) as fh:
        for r in csv.DictReader(fh):
            if r["relationship"] == "meets":
                out.append((float(r["t"]), r["event"], tuple(sorted((r["source"], r["target"])))))
    return out


def statistics(events, T):
    """Durations of completed contacts, inter-contact times of pairs, no. of contacts, group sizes sampled per minute."""
    start, last_end = {}, {}
    dur, inter = [], []
    n_contacts = 0
    for t, ev, p in events:
        if ev == "establish":
            n_contacts += 1
            if p in last_end:
                inter.append(t - last_end[p])
            start[p] = t
        else:
            dur.append(t - start.pop(p))
            last_end[p] = t
    # group sizes per minute (size-biased: the size of the group of each student in a group of at least two):
    sizes = []
    active = set()
    ei = 0
    for minute in np.arange(0.5, T, 1.0):
        while ei < len(events) and events[ei][0] <= minute:
            t, ev, p = events[ei]
            if ev == "establish":
                active.add(p)
            else:
                active.discard(p)
            ei += 1
        parent = {}
        def find(x):
            while parent.setdefault(x, x) != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x
        for a, b in active:
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[ra] = rb
        comp = {}
        for x in parent:
            r = find(x)
            comp[r] = comp.get(r, 0) + 1
        for c in comp.values():
            sizes += [c] * c
    return {"dur": np.array(dur), "inter": np.array(inter), "n": n_contacts, "sizes": np.array(sizes)}


def T_of(day):
    return 300.0 if day == 2 else 540.0


def main():
    res = {"data": {}, "models": {}}
    data_stats = {d: statistics(read_events(os.path.join(HERE, "data", "events_%d.csv" % d)), T_of(d)) for d in DAYS}
    pooled = lambda S, key: np.concatenate([S[d][key] for d in DAYS])
    res["data"] = {"contacts": {d: data_stats[d]["n"] for d in DAYS}, "mean_duration": float(pooled(data_stats, "dur").mean()),
                   "mean_group": float(pooled(data_stats, "sizes").mean())}
    samples = {"data_dur": pooled(data_stats, "dur"), "data_inter": pooled(data_stats, "inter"), "data_sizes": pooled(data_stats, "sizes")}
    tmp = tempfile.mkdtemp()
    for name in MODELS:
        fitres = json.load(open(os.path.join(HERE, "results", "fit_%s.json" % name)))
        est = fitres["estimates"]
        sim_stats = []
        resid = []
        if name == "C1":
            m = c1.C1(DAYS, fitres["age0"])
            th = np.array([est[k] for k in m.names])
            students = m.data["students"]
            n = len(students)
            pairs = [(students[i], students[j]) for i in range(n) for j in range(i + 1, n)]
            for d in DAYS:
                resid.append(m.residuals(th, d))
                for s in SEEDS:
                    ev = m.simulate(th, d, np.random.default_rng(100 * d + s))
                    sim_stats.append((d, statistics([(t, e, pairs[k]) for t, e, k in ev], T_of(d))))
        else:
            cfg = os.path.join(HERE, "configs", "%s.yaml" % name)
            args = []
            for k, v in est.items():
                args += ["--" + k, repr(float(v))]
            for d in DAYS:
                forced = os.path.join(HERE, "data", "forced_%d.csv" % d)
                rf = os.path.join(tmp, "res_%s_%d.csv" % (name, d))
                subprocess.run([fit.BIN, cfg, "--quiet", "--events-in", os.path.join(HERE, "data", "events_%d.csv" % d), "--forced-in", forced,
                                "--max-t", str(T_of(d)), "--residuals-out", rf] + args, check=True, capture_output=True)
                import csv
                resid.append(np.array([float(r["compensator"]) for r in csv.DictReader(open(rf))]))
                for s in SEEDS:
                    ef = os.path.join(tmp, "sim_%s_%d_%d.csv" % (name, d, s))
                    subprocess.run([fit.BIN, cfg, "--quiet", "--forced-in", forced, "--max-t", str(T_of(d)), "--seed", str(s),
                                    "--events-out", ef] + args, check=True, capture_output=True)
                    sim_stats.append((d, statistics(read_events(ef), T_of(d))))
        r = np.concatenate(resid)
        dur = np.concatenate([st["dur"] for d, st in sim_stats])
        inter = np.concatenate([st["inter"] for d, st in sim_stats])
        sizes = np.concatenate([st["sizes"] for d, st in sim_stats])
        res["models"][name] = {
            "residual_KS": float(stats.kstest(r, "expon").statistic), "residual_mean": float(r.mean()), "n_residuals": int(len(r)),
            "contacts_per_day": float(np.mean([st["n"] for d, st in sim_stats])),
            "mean_duration": float(dur.mean()), "mean_group": float(sizes.mean()),
            "KS_duration": float(stats.ks_2samp(dur, samples["data_dur"]).statistic),
            "KS_inter": float(stats.ks_2samp(inter, samples["data_inter"]).statistic),
            "KS_group": float(stats.ks_2samp(sizes, samples["data_sizes"]).statistic),
            "P_dur_ge_10min": float((dur >= 10).mean()), "P_inter_le_2min": float((inter <= 2).mean())}
        samples["%s_dur" % name], samples["%s_inter" % name], samples["%s_sizes" % name], samples["%s_resid" % name] = dur, inter, sizes, r
        print(name, json.dumps(res["models"][name]))
    d = samples["data_dur"]; i = samples["data_inter"]
    res["data"].update({"P_dur_ge_10min": float((d >= 10).mean()), "P_inter_le_2min": float((i <= 2).mean())})
    json.dump(res, open(os.path.join(HERE, "results", "predictive.json"), "w"), indent=1)
    np.savez_compressed(os.path.join(HERE, "results", "samples.npz"), **samples)
    print("data", json.dumps(res["data"]))


if __name__ == "__main__":
    main()

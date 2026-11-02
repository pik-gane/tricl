#!/usr/bin/env python3
"""Figures and summary table of the contact study (study/results -> study/figures, study/results/summary.md)."""
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
R = os.path.join(HERE, "results")
F = os.path.join(HERE, "figures")
LABEL = {"M0": "M0 timetable", "M1": "M1 + static context", "M2": "M2 + network context", "M4": "M4 = M1 + memory",
         "M3": "M3 = M2 + memory", "M6": "M6 = M1 + memory (4 scales)", "M5": "M5 = M2 + memory (4 scales)",
         "C1": "C1 = M1 + power-law aging"}
COLORS = {"data": "black", "M1": "#9aa5b1", "M2": "#4c78a8", "M3": "#f58518", "M5": "#e45756", "C1": "#54a24b", "M4": "#b279a2"}


def ccdf(ax, x, **kw):
    x = np.sort(x)
    y = 1 - np.arange(len(x)) / len(x)
    ax.step(x, y, where="post", **kw)


def main(tag=""):
    os.makedirs(F, exist_ok=True)
    rows = []
    for name in ["M0", "M1", "M2", "M4", "M3", "M6", "M5", "C1"]:
        p = os.path.join(R, "fit_%s%s.json" % (name, tag))
        if not os.path.exists(p):
            continue
        r = json.load(open(p))
        rows.append((name, r["n_params"], r["logl_train"], r["logl_test"], r["events_test"]))
    base = dict((n, (lt, lte)) for n, k, lt, lte, ne in rows).get("M1")
    lines = ["| model | parameters | log-likelihood (train) | held-out log-likelihood | held-out per event vs M1 |", "|---|---|---|---|---|"]
    for n, k, lt, lte, ne in rows:
        lines.append("| %s | %d | %.0f | %.0f | %+.3f |" % (LABEL.get(n, n), k, lt, lte, (lte - base[1]) / ne if base else float("nan")))
    open(os.path.join(R, "summary%s.md" % tag), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))
    sp = os.path.join(R, "samples%s.npz" % tag)
    if not os.path.exists(sp):
        return
    S = np.load(sp)
    fig, axes = plt.subplots(1, 4, figsize=(15, 3.6))
    for ax, key, xl in ((axes[0], "dur", "contact duration (min)"), (axes[1], "inter", "inter-contact time of a pair (min)")):
        ccdf(ax, S["data_" + key], color="black", lw=2, label="data (held-out days)")
        for m in ["M1", "M2", "M3", "M5", "C1"]:
            if m + "_" + key in S:
                ccdf(ax, S[m + "_" + key], color=COLORS.get(m), lw=1.2, label=LABEL[m])
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel(xl); ax.set_ylabel("fraction >= x")
    axes[0].legend(fontsize=7)
    ax = axes[2]
    for m in ["data", "M1", "M2", "M3", "M5", "C1"]:
        k = m + "_sizes"
        if k in S:
            v, c = np.unique(S[k], return_counts=True)
            ax.plot(v, c / c.sum(), "o-" if m == "data" else "-", color=COLORS.get(m), ms=3, label=m)
    ax.set_yscale("log"); ax.set_xlabel("size of a student's group"); ax.set_ylabel("fraction"); ax.set_xlim(1.5, 15)
    ax = axes[3]
    for m in ["M1", "M2", "M3", "M5", "C1"]:
        k = m + "_resid"
        if k in S:
            u = np.sort(1 - np.exp(-S[k]))
            ax.plot(np.linspace(0, 1, len(u)), u, color=COLORS.get(m), label=m)
    ax.plot([0, 1], [0, 1], "k:", lw=1)
    ax.set_xlabel("uniform quantile"); ax.set_ylabel("rescaled termination time 1 - exp(-compensator)")
    fig.tight_layout()
    fig.savefig(os.path.join(F, "checks%s.png" % tag), dpi=130)
    print("figure written")


if __name__ == "__main__":
    import sys
    main(sys.argv[1] if len(sys.argv) > 1 else "")

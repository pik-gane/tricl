#!/usr/bin/env python3
"""Generate every figure, table and number of the paper from the study results.

Usage (after the study scripts in study/ have been run, see study/README.md):
    python3 paper/make_figures.py

Reads study/results/*.json, study/results/residual_samples.npz, study/results/closure_M2.json and the event
files in study/data; writes paper/generated/: numbers.tex (one macro per number used in the text), table-*.tex,
fig-*.pdf. Never edit the generated files by hand.
"""
import csv
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
S = os.path.join(ROOT, "study")
RES = os.path.join(S, "results")
GEN = os.path.join(HERE, "generated")
NUMBERS = {}
C = {"data": "black", "M1": "#9aa5b1", "M2": "#4c78a8", "M6": "#b279a2", "M5": "#e45756", "C1": "#54a24b", "simple": "#4c78a8", "pair": "#e45756", "sim": "#f58518"}
plt.rcParams.update({"font.size": 8, "axes.linewidth": 0.6})

LABEL = {"M0": "M0", "M1": "M1", "M2": "M2", "M4": "M4", "M3": "M3", "M6": "M6", "M5": "M5", "C1": "C1"}
DESC = {"M0": "timetable", "M1": "M0 + static context", "M2": "M1 + network",
        "M4": "M1 + memory (2)", "M3": "M2 + memory (2)",
        "M6": "M1 + memory (4)", "M5": "M2 + memory (4)", "C1": "M1 + power-law aging"}


def num(name, value):
    NUMBERS[name] = value


def fmt(x, d=0):
    s = ("{:,.%df}" % d).format(x)
    return s.replace(",", "\\,").replace("-", "$-$")


def load(name, tag=""):
    p = os.path.join(RES, "fit_%s%s.json" % (name, tag))
    return json.load(open(p)) if os.path.exists(p) else None


def contacts(day, data="data"):
    """(start, end or None) of all contacts of a day."""
    start, out = {}, []
    with open(os.path.join(S, data, "events_%d.csv" % day)) as fh:
        for r in csv.DictReader(fh):
            p = (r["source"], r["target"])
            if r["event"] == "establish":
                start[p] = float(r["t"])
            else:
                out.append((start.pop(p), float(r["t"])))
    return out


def data_numbers():
    info = json.load(open(os.path.join(S, "data", "data.json")))["info"]
    info120 = json.load(open(os.path.join(S, "data_gap120", "data.json")))["info"]
    num("Students", "%d" % info["students"])
    num("Classes", "%d" % len(info["classes"]))
    num("Contacts", fmt(info["contacts_total"]))
    num("ContactsGapLong", fmt(info120["contacts_total"]))
    num("FragmentShare", "%d" % round(100 * (1 - info120["contacts_total"] / info["contacts_total"])))
    data = json.load(open(os.path.join(S, "data", "data.json")))
    cls = dict(data["static"]["class"])
    within = total = 0
    durs = []
    for d in (2, 3, 4, 5, 6):
        with open(os.path.join(S, "data", "events_%d.csv" % d)) as fh:
            for r in csv.DictReader(fh):
                if r["event"] == "establish":
                    total += 1
                    within += cls[r["source"]] == cls[r["target"]]
        durs += [e - s for s, e in contacts(d)]
    durs = np.array(durs)
    num("WithinClassShare", "%d" % round(100 * within / total))
    num("MedianDuration", "%.1f" % np.median(durs))
    num("ShareLongContacts", "%.1f" % (100 * np.mean(durs >= 10)))
    num("Facebook", fmt(len(data["static"]["facebook"])))
    num("Survey", fmt(len(data["static"]["survey"])))
    return durs


def table_fits():
    rows = []
    resid = json.load(open(os.path.join(RES, "residuals.json")))
    base = load("M1")
    for m in ["M0", "M1", "M2", "M4", "M3", "M6", "M5", "C1"]:
        r = load(m)
        if r is None:
            continue
        per = (r["logl_test"] - base["logl_test"]) / r["events_test"]
        rows.append("%s & %s & %d & %s & %s & %s & %.3f \\\\" % (m, DESC[m], r["n_params"], fmt(r["logl_train"]), fmt(r["logl_test"]),
                                                               ("$%+.3f$" % per).replace("-", "-"), resid[m]["KS"]))
    with open(os.path.join(GEN, "table-fits.tex"), "w") as fh:
        fh.write("\\begin{ruledtabular}\n\\begin{tabular}{llrrrrr}\n")
        fh.write("model & components & param. & $L$ (training) & $L$ (held out) & $\\Delta L$/event & KS \\\\\n\\colrule\n")
        fh.write("\n".join(rows) + "\n\\end{tabular}\n\\end{ruledtabular}\n")
    r1 = load("M1")
    num("EventsTrain", fmt(r1["events_train"]))
    num("EventsTest", fmt(r1["events_test"]))
    g = {m: load(m)["logl_test"] for m in ["M0", "M1", "M2", "M4", "M3", "M6", "M5", "C1"]}
    n = r1["events_test"]
    num("GainStatic", "%.2f" % ((g["M1"] - g["M0"]) / n))
    num("GainNetwork", "%.2f" % ((g["M2"] - g["M1"]) / n))
    num("GainMemoryAlone", "%.2f" % ((g["M6"] - g["M1"]) / n))
    num("GainMemoryWithNetwork", "%.2f" % ((g["M5"] - g["M2"]) / n))
    num("GainNetworkWithMemory", "%.2f" % ((g["M5"] - g["M6"]) / n))
    num("MemoryAbsorbed", "%d" % round(100 * (1 - (g["M5"] - g["M2"]) / (g["M6"] - g["M1"]))))
    num("MfiveOverCone", fmt(g["M5"] - g["C1"]))
    num("ConeOverMthree", fmt(g["C1"] - g["M3"]))
    num("KSMone", "%.3f" % resid["M1"]["KS"])
    num("KSMfive", "%.3f" % resid["M5"]["KS"])
    num("KSCone", "%.3f" % resid["C1"]["KS"])
    num("ConeAgeScale", "%.0f" % load("C1")["age0"])


def table_estimates():
    m5, c1 = load("M5"), load("C1")
    rows = []
    spec = [("F_cls", "same class"), ("F_fb", "Facebook friends"), ("F_svy", "survey friendship"),
            ("F_tri", "common contact"), ("F_deg", "other current contact"), ("F_mb", "met on earlier day")]
    spec_mem = [("d2", "2"), ("d05", "0.5"), ("d01", "0.1"), ("d002", "0.02")]
    def cell(r, k):
        if r is None or k not in r["estimates"]:
            return "& "
        return "$%.3f$ & $%.3f$" % (r["estimates"][k], r["standard_errors"][k])
    for k, lab in spec:
        rows.append("%s & %s & %s \\\\" % (lab, cell(m5, k), cell(m5, "T" + k[1:])))
    for s, d in spec_mem:
        rows.append("memory, $\\delta = %s$/min & %s & %s \\\\" % (d, cell(m5, "F_end_" + s), cell(m5, "T_start_" + s)))
    with open(os.path.join(GEN, "table-estimates.tex"), "w") as fh:
        fh.write("\\begin{ruledtabular}\n\\begin{tabular}{lrrrr}\n")
        fh.write(" & \\multicolumn{2}{c}{formation} & \\multicolumn{2}{c}{termination} \\\\\neffect & estimate & s.e. & estimate & s.e. \\\\\n\\colrule\n")
        fh.write("\n".join(rows) + "\n\\end{tabular}\n\\end{ruledtabular}\n")
    e = m5["estimates"]
    num("TriFactor", "%.0f" % np.exp(4 * e["F_tri"]))
    num("ClassFactor", "%.0f" % np.exp(4 * e["F_cls"]))
    num("ConeTermExponent", "$%.2f$" % (4 * c1["estimates"]["T_age"]))
    num("ConeFormExponent", "$%.2f$" % (4 * c1["estimates"]["F_age"]))


def table_sensitivity():
    rows = []
    base = load("M1", "_g120")
    for m in ["M0", "M1", "M2", "M4", "M3", "C1"]:
        r = load(m, "_g120")
        if r is None:
            continue
        rows.append("%s & %s & %s & %s \\\\" % (m, fmt(r["logl_test"]), ("$%+.3f$" % ((r["logl_test"] - base["logl_test"]) / r["events_test"])),
                                                  ("$%+.3f$" % ((load(m)["logl_test"] - load("M1")["logl_test"]) / load(m)["events_test"]))))
    with open(os.path.join(GEN, "table-sensitivity.tex"), "w") as fh:
        fh.write("\\begin{ruledtabular}\n\\begin{tabular}{lrrr}\n")
        fh.write("model & $L$ (held out) & $\\Delta L$/event & $\\Delta L$/event, 40-s gaps \\\\\n\\colrule\n")
        fh.write("\n".join(rows) + "\n\\end{tabular}\n\\end{ruledtabular}\n")
    r = {m: load(m, "_g120") for m in ["M1", "M2", "M4"]}
    n = r["M1"]["events_test"]
    num("GainMemoryAloneLong", "%.2f" % ((r["M4"]["logl_test"] - r["M1"]["logl_test"]) / n))
    num("GainNetworkLong", "%.2f" % ((r["M2"]["logl_test"] - r["M1"]["logl_test"]) / n))
    c1l, m3l = load("C1", "_g120"), load("M3", "_g120")
    if c1l is not None and m3l is not None:
        num("ConeLong", "%.2f" % ((c1l["logl_test"] - r["M1"]["logl_test"]) / n))
        num("MthreeLong", "%.2f" % ((m3l["logl_test"] - r["M1"]["logl_test"]) / n))
    n40 = load("M1")["events_test"]
    num("GainMemoryAloneTwo", "%.2f" % ((load("M4")["logl_test"] - load("M1")["logl_test"]) / n40))


def fig_data(durs):
    fig, axes = plt.subplots(1, 2, figsize=(6.8, 2.3))
    ax = axes[0]
    for d, ls in ((3, "-"), (4, "--"), (5, ":"), (6, "-.")):
        st = np.array([s for s, e in contacts(d)])
        h, edges = np.histogram(st, bins=np.arange(0, 541, 20))
        ax.step(8 + edges[:-1] / 60, h / 20, where="post", ls=ls, lw=0.9, color="black", label=["Tue", "Wed", "Thu", "Fri"][d - 3])
    ax.set_xlabel("time of day (h)"); ax.set_ylabel("contacts starting per minute"); ax.legend(fontsize=6, frameon=False)
    ax = axes[1]
    x = np.sort(durs); y = 1 - np.arange(len(x)) / len(x)
    ax.loglog(x, y, color="black", lw=1.2, label="data")
    lam = 1 / durs.mean()
    xx = np.logspace(np.log10(x.min()), np.log10(x.max()), 100)
    ax.loglog(xx, np.exp(-lam * xx), "k:", lw=0.8, label="exponential, same mean")
    ax.set_ylim(1 / len(x), 1.2); ax.set_xlabel("contact duration (min)"); ax.set_ylabel("fraction of contacts $\\geq$ duration")
    ax.legend(fontsize=6, frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(GEN, "fig-data.pdf"))
    plt.close(fig)


def fig_residuals():
    R = np.load(os.path.join(RES, "residual_samples.npz"))
    fig, ax = plt.subplots(figsize=(3.3, 2.6))
    for m in ["M1", "M2", "M6", "M5", "C1"]:
        u = np.sort(1 - np.exp(-R[m]))
        ax.plot(np.linspace(0, 1, len(u)), u, color=C[m], lw=1, label="%s (%s)" % (m, DESC[m]))
    ax.plot([0, 1], [0, 1], "k:", lw=0.8)
    ax.set_xlabel("quantile of the uniform distribution"); ax.set_ylabel("rescaled contact duration $1 - e^{-\\Lambda}$")
    ax.legend(fontsize=5.5, frameon=False, loc="upper left")
    fig.tight_layout()
    fig.savefig(os.path.join(GEN, "fig-residuals.pdf"))
    plt.close(fig)


def fig_closure():
    p = os.path.join(RES, "closure_M2.json")
    if not os.path.exists(p):
        return
    c = json.load(open(p))
    g = np.array(c["grid"])
    fig, ax = plt.subplots(figsize=(3.3, 2.5))
    ax.semilogy(8 + g / 60, np.maximum(c["data"]["L"], 0.5), color="black", lw=1, label="data")
    sim = np.array(c["simulation"]["L"])
    ok = sim > 0
    ax.semilogy(8 + g[ok] / 60, sim[ok], color=C["sim"], lw=1.2, label="simulation of fitted M2")
    ax.semilogy(8 + g / 60, c["simple"]["L"], color=C["simple"], lw=1, ls="--", label="simple closure")
    ax.semilogy(8 + g / 60, c["pair"]["L"], color=C["pair"], lw=1, label="pair closure")
    ax.set_xlabel("time of day (h)"); ax.set_ylabel("contacts within classes")
    ax.legend(fontsize=6, frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(GEN, "fig-closure.pdf"))
    plt.close(fig)
    L = np.array(c["pair"]["L"]); Ls = np.array(c["simple"]["L"])
    hi = 0.5 * max(L.max(), 1)
    num("PairTipMinutes", "%.0f" % g[np.argmax(L > hi)])
    num("SimTipMinutes", "%.0f" % g[np.argmax(sim > hi)])
    num("SimpleTipMinutes", "%.0f" % g[np.argmax(Ls > hi)])
    num("DataWithinContacts", "%.0f" % np.mean(c["data"]["L"]))
    num("DenseWithinContacts", "%.0f" % L[-1])


def main():
    os.makedirs(GEN, exist_ok=True)
    durs = data_numbers()
    table_fits()
    table_estimates()
    table_sensitivity()
    fig_data(durs)
    fig_residuals()
    fig_closure()
    with open(os.path.join(GEN, "numbers.tex"), "w") as fh:
        fh.write("% generated by paper/make_figures.py -- do not edit\n")
        for k in sorted(NUMBERS):
            fh.write("\\newcommand{\\%s}{%s}\n" % (k, NUMBERS[k]))
    for k in sorted(NUMBERS):
        print("%-24s %s" % (k, NUMBERS[k]))


if __name__ == "__main__":
    main()

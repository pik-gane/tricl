#!/usr/bin/env python3
"""Make the figures, tables and the numbers quoted in the text of the paper from paper/data/.

Usage (from the repository root, after paper/make_data.py):
    python3 paper/make_figures.py

Writes paper/generated/: fig-frames.pdf, fig-macro.pdf, fig-atop.pdf, table-benchmark.tex,
table-recovery.tex, table-atop.tex and numbers.tex (LaTeX macros for every number quoted in the text).
Never edit the generated files by hand.
"""

import csv
import json
import math
import os
import re
import statistics
import subprocess
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.image as mpimg  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "paper", "data")
GEN = os.path.join(ROOT, "paper", "generated")
BIN = os.path.join(ROOT, "build", "src", "tricl")

# validated categorical palette (light mode) and ink colours
C1, C2, C3 = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#898781", "#e1e0d9"

plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 9,
    "axes.edgecolor": "#c3c2b7", "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 0.6, "lines.linewidth": 1.5, "legend.frameon": False, "legend.fontsize": 7.5,
    "text.color": INK, "pdf.fonttype": 42,
})

numbers = {}


def num(name, value):
    numbers[name] = value


def fmt_sig(x, digits=2):
    """Format with the given number of significant digits, in LaTeX notation."""
    if x == 0:
        return "0"
    e = int(math.floor(math.log10(abs(x))))
    if -3 <= e <= 3:
        d = max(0, digits - 1 - e)
        return ("%.*f" % (d, x))
    m = x / 10 ** e
    return r"%.*f\times 10^{%d}" % (digits - 1, m, e)


def expit(x):
    return 1 / (1 + math.exp(-x))


# ---------------------------------------------------------------- frames

def fig_frames():
    fdir = os.path.join(DATA, "frames")
    files = sorted(f for f in os.listdir(fdir) if f.endswith(".png"))
    fig, axes = plt.subplots(1, len(files), figsize=(7.2, 2.5))
    for ax, f in zip(axes, files):
        ax.imshow(mpimg.imread(os.path.join(fdir, f)))
        ax.axis("off")
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0, wspace=0.02)
    fig.savefig(os.path.join(GEN, "fig-frames.pdf"))
    plt.close(fig)


# ---------------------------------------------------------------- mean field

def read_macro(name, link_type):
    rows = [r for r in csv.DictReader(open(os.path.join(DATA, "macro_%s.csv" % name))) if r["link type"] == link_type]
    t = [float(r["t"]) for r in rows]
    return t, [float(r["ode"]) for r in rows], [float(r["sim_mean"]) for r in rows], [float(r["sim_sd"]) for r in rows]


def fig_macro():
    cases = [("linear", "agent is active", "agents in state 1", "linear switching"),
             ("si", "agent is infected", "infected agents", "infection on a geometric graph"),
             ("granovetter", "c-agent is active", "active contingent agents", "threshold model")]
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.3))
    for ax, (name, lt, ylabel, title) in zip(axes, cases):
        t, ode, m, sd = read_macro(name, lt)
        ax.fill_between(t, [a - b for a, b in zip(m, sd)], [a + b for a, b in zip(m, sd)], color=C2, alpha=0.12, linewidth=0)
        ax.plot(t, m, "o", color=C2, markersize=3.5, markeredgecolor="white", markeredgewidth=0.6, label="simulation")
        ax.plot(t, ode, "-", color=C1, label="mean field")
        ax.set_title(title, color=INK, loc="left")
        ax.set_xlabel("time")
        ax.set_ylabel(ylabel)
        ax.set_ylim(bottom=0)
        if name == "si":
            num("MacroSiOde", "%.0f" % ode[-1]); num("MacroSiSim", "%.0f" % m[-1])
        if name == "granovetter":
            num("MacroGranOde", "%.0f" % ode[-1]); num("MacroGranSim", "%.0f" % m[-1])
    axes[0].legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(os.path.join(GEN, "fig-macro.pdf"))
    plt.close(fig)


# ---------------------------------------------------------------- recovery

def table_recovery():
    d = json.load(open(os.path.join(DATA, "recovery.json")))
    true = d["true"]
    res = d["results"]
    names = list(true)
    labels = {"A": "$A$", "K": "$K$", "D": "$D$", "BK": "$B_K$"}
    lines = [r"\begin{ruledtabular}", r"\begin{tabular}{lrrrrr}",
             r"parameter & true & mean est. & s.d.\ est. & mean s.e. & coverage \\", r"\colrule"]
    covs = []
    for n in names:
        est = [r["estimates"][n] for r in res]
        se = [r["se"][n] for r in res]
        cov = sum(1 for e, s in zip(est, se) if abs(e - true[n]) <= 1.96 * s) / len(res)
        covs.append(cov)
        lines.append(r"%s & %.3g & %.3f & %.3f & %.3f & %.2f \\" % (labels.get(n, n), true[n], statistics.mean(est),
                                                                 statistics.stdev(est), statistics.mean(se), cov))
    lines += [r"\end{tabular}", r"\end{ruledtabular}"]
    open(os.path.join(GEN, "table-recovery.tex"), "w").write("\n".join(lines) + "\n")
    num("RecoverySeeds", str(len(res)))
    num("RecoveryT", "400")
    num("RecoveryCovMin", "%.2f" % min(covs))
    num("RecoveryCovMax", "%.2f" % max(covs))
    num("RecoveryEvents", "%d" % round(statistics.mean(r["events"] for r in res), -1))


# ---------------------------------------------------------------- benchmark

def table_benchmark():
    rows = list(csv.DictReader(open(os.path.join(DATA, "benchmark.csv"))))
    last = rows[-1]["timestamp"]
    rows = [r for r in rows if r["timestamp"] == last]
    names = {"helfmann N=100": ("Helfmann", "100"), "helfmann N=1000": ("Helfmann", "1\\,000"),
             "helfmann N=10000": ("Helfmann", "10\\,000"), "3blocks N=300 (dense)": ("3-blocks (dense)", "300"),
             "sir_sd N=1000": ("SIR with distancing", "1\\,000")}
    lines = [r"\begin{ruledtabular}", r"\begin{tabular}{llrrr}",
             r"model & agents & events & init. [s] & events/s \\", r"\colrule"]
    for r in rows:
        m, n = names.get(r["case"], (r["case"], ""))
        rate = float(r["events_per_second"])
        lines.append(r"%s & %s & %s & %.2f & %s \\" % (m, n, "{:,}".format(int(r["events"])).replace(",", "\\,"),
                                                      float(r["init_seconds"]), "{:,}".format(int(round(rate, -2))).replace(",", "\\,")))
    lines += [r"\end{tabular}", r"\end{ruledtabular}"]
    open(os.path.join(GEN, "table-benchmark.tex"), "w").write("\n".join(lines) + "\n")


# ---------------------------------------------------------------- alliances

def table_and_fig_atop():
    adir = os.path.join(DATA, "atop")
    full = json.load(open(os.path.join(adir, "full.json")))
    shuf = json.load(open(os.path.join(DATA, "atop", "shuffled.json")))
    est, se = full["estimates"], full["standard_errors"]
    info = open(os.path.join(adir, "build.log")).read()
    m = re.search(r"(\d+) states, (\d+) initial links, (\d+) events, link intervals \{'defense': (\d+)\} \((\d+) censored", info)
    num("AtopStates", m.group(1)); num("AtopEvents", "{:,}".format(int(m.group(3))).replace(",", "\\,"))
    num("AtopIntervals", "{:,}".format(int(m.group(4))).replace(",", "\\,")); num("AtopCensored", "{:,}".format(int(m.group(5))).replace(",", "\\,"))

    a0, a1, d0, b1 = est["a0_defense"], est["a1_defense"], est["d0_defense"], est["b1_defense"]
    num("AtopAzero", fmt_sig(a0, 3)); num("AtopAone", fmt_sig(a1, 3))
    num("AtopFormationRatio", "%d" % round((a0 + a1) / a0))
    # defense is symmetric: both directed links of a pair are possible events with these rates, so a pair forms a
    # pact at twice the rate and a pact dissolves at twice the hazard (paper/sec-appendix.tex, likelihood section):
    num("AtopPairAzero", fmt_sig(2 * a0, 3)); num("AtopPairAone", fmt_sig(2 * a1, 3))
    haz = lambda k: 2 * d0 * expit(-4 * b1 * k)
    num("AtopHazZero", "%.3f" % haz(0)); num("AtopHazFive", "%.3f" % haz(5)); num("AtopHazTwentySix", "%.4f" % haz(26))
    num("AtopLifeZero", "%.0f" % (1 / haz(0)))
    logl = full["logl"]
    no_a1 = json.load(open(os.path.join(adir, "no_a1.json")))["logl"]
    no_b1 = json.load(open(os.path.join(adir, "no_b1.json")))["logl"]
    num("AtopDropAone", "{:,.0f}".format(logl - no_a1).replace(",", "\\,"))
    num("AtopDropBone", "{:,.0f}".format(logl - no_b1).replace(",", "\\,"))
    for key, name in (("a0_defense", "AtopShuffleAzero"), ("a1_defense", "AtopShuffleAone"), ("b1_defense", "AtopShuffleBone")):
        num(name, "%.0f" % (100 * abs(shuf["estimates"][key] - est[key]) / est[key]))

    labels = [("a0_defense", "$a_0$", "formation rate per pair of member states and direction"),
              ("a1_defense", "$a_1$", "additional rate per common ally and direction"),
              ("d0_defense", "$d_0$", "dissolution attempt rate per direction"),
              ("b1_defense", "$b_1$", "probability units per common ally (dissolution)"),
              ("entry", "entry", "entry rate of states"), ("exit", "exit", "exit rate of states")]
    lines = [r"\begin{ruledtabular}", r"\begin{tabular}{llrr}", r"parameter & meaning & estimate & s.e. \\", r"\colrule"]
    for k, sym, meaning in labels:
        lines.append(r"%s & %s & $%s$ & $%s$ \\" % (sym, meaning, fmt_sig(est[k], 3), fmt_sig(se[k], 2)))
    lines += [r"\colrule", r"\multicolumn{2}{l}{log-likelihood} & \multicolumn{2}{r}{$%s$} \\" % ("{:,.1f}".format(logl).replace(",", "\\,")),
              r"\end{tabular}", r"\end{ruledtabular}"]
    open(os.path.join(GEN, "table-atop.tex"), "w").write("\n".join(lines) + "\n")

    # profile
    prof = [(float(r["b1"]), float(r["logl"])) for r in csv.DictReader(open(os.path.join(adir, "profile_b1.csv")))]
    plateau = min(prof, key=lambda p: p[0])[1]
    num("AtopProfileGap", "{:,.0f}".format(logl - plateau).replace(",", "\\,"))

    # pacts in effect per year
    model = open(os.path.join(adir, "model.yaml")).read()
    n = len(re.findall(r"- \[[^,]+, defense, [^\]]+\]", model))
    by_year = {}
    for r in csv.DictReader(open(os.path.join(adir, "events.csv"))):
        if r["relationship"] != "defense":
            continue
        y = 1816 + int(float(r["t"]))
        by_year[y] = by_year.get(y, 0) + (1 if r["event"] == "establish" else -1)
    years, counts = list(range(1816, 2019)), []
    for y in years:
        n += by_year.get(y, 0)
        counts.append(n)

    # timing of one likelihood evaluation
    times = []
    for _ in range(3):
        t0 = time.perf_counter()
        subprocess.run([BIN, os.path.join(adir, "model.yaml"), "--events-in", os.path.join(adir, "events.csv"), "--summary"],
                       capture_output=True, check=True, cwd=adir)
        times.append(time.perf_counter() - t0)
    num("AtopEvalSeconds", "%.2f" % min(times))

    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.3))
    ax = axes[0]
    ax.plot(years, counts, color=C1)
    ax.set_xlabel("year"); ax.set_ylabel("pairs with a defense pledge")
    ax.set_ylim(bottom=0)
    ax = axes[1]
    bx = [p[0] for p in prof]; by = [p[1] - logl for p in prof]
    ax.plot(bx, by, "-o", color=C1, markersize=3.5, markeredgecolor="white", markeredgewidth=0.6)
    ax.plot([b1], [0], "o", color=C2, markersize=5, markeredgecolor="white", markeredgewidth=0.6, zorder=3)
    ax.annotate(r"$\hat b_1$", (b1, 0), xytext=(6, -10), textcoords="offset points", color=INK2, fontsize=7.5)
    ax.set_xscale("symlog", linthresh=0.01)
    ax.set_xticks([-1, -0.1, -0.01, 0, 0.01, 0.1])
    ax.set_xticklabels(["-1", "-0.1", "-0.01", "0", "0.01", "0.1"], fontsize=6.5)
    ax.minorticks_off()
    ax.set_xlabel(r"$b_1$"); ax.set_ylabel("profile log-likelihood")
    ax = axes[2]
    ks = list(range(0, 31))
    ax.plot(ks, [haz(k) for k in ks], color=C1)
    ax.set_xlabel("number of common allies"); ax.set_ylabel("dissolution hazard per year")
    ax.set_ylim(bottom=0)
    fig.tight_layout()
    fig.savefig(os.path.join(GEN, "fig-atop.pdf"))
    plt.close(fig)


def sigmoid_checks():
    out = subprocess.run([os.path.join(ROOT, "build", "tests", "test_sigmoid")], capture_output=True, text=True).stdout
    m = re.search(r"(\d+) of (\d+) checks passed", out)
    num("TestSigmoidChecks", "{:,}".format(int(m.group(2))).replace(",", "\\,"))


def main():
    os.makedirs(GEN, exist_ok=True)
    fig_frames()
    fig_macro()
    table_recovery()
    table_benchmark()
    table_and_fig_atop()
    sigmoid_checks()
    with open(os.path.join(GEN, "numbers.tex"), "w") as f:
        f.write("% generated by paper/make_figures.py -- do not edit\n")
        for k, v in sorted(numbers.items()):
            f.write("\\newcommand{\\%s}{%s}\n" % (k, ("\\ensuremath{%s}" % v) if "\\times" in v else v))
    for k, v in sorted(numbers.items()):
        print("%-22s %s" % (k, v))


if __name__ == "__main__":
    main()

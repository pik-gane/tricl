#!/usr/bin/env python3
"""A first fit of tricl to real data: the formation and dissolution of military alliances (ATOP data).

Usage:
    python3 python/tricl_atop.py [--data DIR] [--out DIR] [--types defense ...] [--from YEAR] [--to YEAR]
                                 [--fit [--bin PATH] [--maxiter N]]

Data: the directed dyad-year alliance data of the Alliance Treaty Obligations and Provisions project (ATOP v5.1,
Leeds et al. 2002) and the Correlates of War state system membership list (v2016), both as shipped in the R package
`peacesciencer` by Steven V. Miller (https://github.com/svmiller/peacesciencer, files data/atop_alliance.rda and
data/cow_states.rda, downloaded into DIR if missing; reading them needs the `pyreadr` package). In the dyad-year
data, a pledge (defense, offense, neutral, nonagg, consul) between two states is either in effect during a year or
not (the pledges are directed, but fewer than 1% of the defense dyad-years lack the mirrored pledge, so a pair is
treated as linked if either state pledges), and every maximal run of consecutive years becomes one link interval: a
link of the symmetric relationship type named after the pledge is established at the start of its first year and
terminated at the start of the year after its last one (links still in effect in the last data year are censored). States enter and leave the international
system; this is encoded, in tricl's usual way, as a link ``[state, is in, system]`` to a single hub entity, which
makes "both states are members of the system" the 2-path ``[~, is in, system, contains, ~]``. Years in which a state
has an alliance but is not a system member (governments in exile, German states before 1816, years after the
end of the membership list) count as membership years, so that no observed event is impossible under the model.

Model (per pledge type T, all rates per year):

    [state, T, state]:
        establish:
            attempt:
                [~, is in, system, contains, ~]: a0_T   # spontaneous formation, member pairs only
                [~, T, state, T, ~]:             a1_T   # additional rate per common T-ally ("ally of an ally")
            success: inf
        terminate:
            attempt:
                base: d0_T                              # dissolution attempt rate (success 1/2 without common allies)
            success:
                tails: 0
                [~, T, state, T, ~]: -b1_T              # each common ally lowers the success probability units
    [state, is in, system]: entry rate `entry` per non-member state, exit rate `exit` per member state

Events of the same year are applied in the order terminations, exits, entries, formations, and formations of the
same year in the order of the state codes (a multilateral treaty appears as many simultaneous dyadic events, so
the later ones already see the earlier ones as common allies; this is a convention).

The script writes OUT/model.yaml (with all states as named entities and the links in effect at the start of the
window as named initial links) and OUT/events.csv, and with --fit runs python/tricl_fit.py on them. The start
values in model.yaml are close to the maximum-likelihood estimates for the defense pledges (see the README);
note that the likelihood has a second, much worse local maximum at negative b1 (a saturated sigmoid), so a fit
started far from the optimum can end there. Check with --shuffle how much the estimates depend on the ordering of
simultaneous events (in the defense data: a0 by about a quarter, a1 and b1 by a few percent).
"""

import argparse
import collections
import csv
import os
import subprocess
import sys
import urllib.request

RAW_URL = "https://raw.githubusercontent.com/svmiller/peacesciencer/master/data/"
FILES = ["atop_alliance.rda", "cow_states.rda"]
PLEDGES = ["defense", "offense", "neutral", "nonagg", "consul"]


def download(data_dir):
    os.makedirs(data_dir, exist_ok=True)
    for f in FILES:
        path = os.path.join(data_dir, f)
        if not os.path.exists(path):
            print("downloading", RAW_URL + f)
            urllib.request.urlretrieve(RAW_URL + f, path)


def load(data_dir):
    """Returns (alliance records, state records) as lists of dicts."""
    import pyreadr
    al = pyreadr.read_r(os.path.join(data_dir, "atop_alliance.rda"))["atop_alliance"]
    st = pyreadr.read_r(os.path.join(data_dir, "cow_states.rda"))["cow_states"]
    return al.to_dict("records"), st.to_dict("records")


def membership_spells(states, alliance_years, last_year):
    """Membership years per state code: the COW spells (extended to last_year if censored at the end of the list)
    plus the years with any alliance; returns {ccode: [(first year, last year), ...]} and {ccode: abbreviation}."""
    years = collections.defaultdict(set)
    abb = {}
    max_end = max(int(r["endyear"]) for r in states)
    for r in states:
        c = int(r["ccode"])
        abb[c] = r["stateabb"]
        end = int(r["endyear"])
        if end == max_end and int(r["endmonth"]) == 12 and int(r["endday"]) == 31:
            end = max(end, last_year)
        years[c].update(range(int(r["styear"]), end + 1))
    for c, ys in alliance_years.items():
        years[c].update(ys)
    spells = {}
    for c, ys in years.items():
        runs, start, prev = [], None, None
        for y in sorted(ys):
            if start is None:
                start = prev = y
            elif y == prev + 1:
                prev = y
            else:
                runs.append((start, prev)); start = prev = y
        runs.append((start, prev))
        spells[c] = runs
    return spells, abb


def intervals(alliances, pledge):
    """Maximal runs of years with the pledge, per unordered pair: {(a, b): [(first, last), ...]}."""
    years = collections.defaultdict(set)
    for r in alliances:
        if r["atop_" + pledge] == 1:
            a, b = int(r["ccode1"]), int(r["ccode2"])
            years[(min(a, b), max(a, b))].add(int(r["year"]))
    out = {}
    for pair, ys in years.items():
        runs, start, prev = [], None, None
        for y in sorted(ys):
            if start is None:
                start = prev = y
            elif y == prev + 1:
                prev = y
            else:
                runs.append((start, prev)); start = prev = y
        runs.append((start, prev))
        out[pair] = runs
    return out


def build(alliances, states, out_dir, pledges=("defense",), first_year=1816, last_year=None, parameters=None, shuffle_seed=None):
    """Writes OUT/model.yaml and OUT/events.csv; returns a dict with counts.

    Events of the same year are ordered by the state codes, or randomly if shuffle_seed is given (a check of how
    much the estimates depend on this convention).
    """
    os.makedirs(out_dir, exist_ok=True)
    if last_year is None:
        last_year = max(int(r["year"]) for r in alliances)
    parameters = dict(parameters or {})
    alliance_years = collections.defaultdict(set)
    for r in alliances:
        if any(r["atop_" + p] == 1 for p in pledges):
            alliance_years[int(r["ccode1"])].add(int(r["year"]))
            alliance_years[int(r["ccode2"])].add(int(r["year"]))
    spells, abb = membership_spells(states, alliance_years, last_year)
    t = lambda year: year - first_year  # model time in years, 0 = start of first_year

    initial, events = [], []  # events: (t, order, event, source, relationship, target)
    n_censored = 0
    for c, runs in spells.items():
        for (s, e) in runs:
            if e < first_year or s > last_year:
                continue
            if s <= first_year:
                initial.append((abb[c], "is in", "system"))
            else:
                events.append((t(s), 2, "establish", abb[c], "is in", "system"))
            if e < last_year:
                events.append((t(e + 1), 1, "terminate", abb[c], "is in", "system"))
    n_links = collections.Counter()
    for p in pledges:
        for (a, b), runs in intervals(alliances, p).items():
            for (s, e) in runs:
                if e < first_year or s > last_year:
                    continue
                n_links[p] += 1
                if s <= first_year:
                    initial.append((abb[a], p, abb[b]))
                else:
                    events.append((t(s), 3, "establish", abb[a], p, abb[b]))
                if e < last_year:
                    events.append((t(e + 1), 0, "terminate", abb[a], p, abb[b]))
                else:
                    n_censored += 1
    if shuffle_seed is None:
        events.sort(key=lambda ev: (ev[0], ev[1], ev[3], ev[5]))
    else:
        import random
        rng = random.Random(shuffle_seed)
        events.sort(key=lambda ev: (ev[0], ev[1], rng.random()))
    with open(os.path.join(out_dir, "events.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["t", "event", "source", "relationship", "target"])
        for ev in events:
            w.writerow([ev[0], ev[2], ev[3], ev[4], ev[5]])

    state_labels = sorted(abb[c] for c in spells)
    lines = ["metadata:",
             "    name: ATOP alliances",
             "    description: formation and dissolution of alliances between states, %d-%d (ATOP v5.1, COW state system)" % (first_year, last_year),
             "    author: generated by python/tricl_atop.py",
             "",
             "metaparameters:"]
    defaults = {"entry": 0.01, "exit": 0.002}
    for p in pledges:
        defaults.update({"a0_" + p: 0.003, "a1_" + p: 0.01, "d0_" + p: 0.1, "b1_" + p: 0.04})
    for name, value in defaults.items():
        lines.append("    %s: %s" % (name, parameters.get(name, value)))
    lines += ["",
              "limits:",
              "    t: %d  # end of %d" % (t(last_year + 1), last_year),
              "",
              "options:",
              "    quiet: true",
              "",
              "entities:",
              "    system:",
              "        - system",
              "    state:"]
    lines += ["        - %s" % s for s in state_labels]
    lines += ["",
              "relationship types:",
              "    is in: contains"]
    lines += ["    %s: symmetric" % p for p in pledges]
    lines += ["",
              "initial links:  # the state at the start of %d" % first_year,
              "    named:"]
    lines += ["        - [%s, %s, %s]" % l for l in initial]
    lines += ["",
              "dynamics:",
              "",
              "    [state, is in, system]:",
              "        establish:",
              "            attempt:",
              "                base: entry  # entry rate per non-member state",
              "            success: inf",
              "        terminate:",
              "            attempt:",
              "                base: exit  # exit rate per member state",
              "            success: inf"]
    for p in pledges:
        lines += ["",
                  "    [state, %s, state]:" % p,
                  "        establish:",
                  "            attempt:",
                  "                [~, is in, system, contains, ~]: a0_%s  # spontaneous formation rate per pair of member states" % p,
                  "                [~, %s, state, %s, ~]: a1_%s  # additional rate per common ally" % (p, p, p),
                  "            success: inf",
                  "        terminate:",
                  "            attempt:",
                  "                base: d0_%s  # dissolution attempt rate; success probability 1/2 without common allies" % p,
                  "            success:",
                  "                tails: 0",
                  "                [~, %s, state, %s, ~]: -b1_%s  # each common ally lowers the success probability units" % (p, p, p)]
    with open(os.path.join(out_dir, "model.yaml"), "w") as f:
        f.write("\n".join(lines) + "\n")
    return {"states": len(state_labels), "initial_links": len(initial), "events": len(events),
            "intervals": dict(n_links), "censored": n_censored, "first_year": first_year, "last_year": last_year,
            "parameters": list(defaults)}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0], formatter_class=argparse.RawDescriptionHelpFormatter,
                                 epilog=__doc__.split("\n", 1)[1])
    ap.add_argument("--data", default="atop_data", help="directory with atop_alliance.rda and cow_states.rda (downloaded if missing)")
    ap.add_argument("--out", default="atop", help="output directory for model.yaml and events.csv (default: atop)")
    ap.add_argument("--types", nargs="+", default=["defense"], choices=PLEDGES, help="pledge types to model (default: defense)")
    ap.add_argument("--from", dest="first_year", type=int, default=1816, help="first year (default 1816)")
    ap.add_argument("--to", dest="last_year", type=int, default=None, help="last year (default: last year in the data)")
    ap.add_argument("--shuffle", type=int, default=None, metavar="SEED", help="order the events of the same year randomly instead of by state codes")
    ap.add_argument("--fit", action="store_true", help="estimate all metaparameters by maximum likelihood (python/tricl_fit.py)")
    ap.add_argument("--bin", default=None, help="tricl binary (default: build/src/tricl next to this repository)")
    ap.add_argument("--maxiter", type=int, default=200)
    args = ap.parse_args(argv)

    download(args.data)
    alliances, states = load(args.data)
    info = build(alliances, states, args.out, args.types, args.first_year, args.last_year, shuffle_seed=args.shuffle)
    print("wrote %s/model.yaml and %s/events.csv: %d states, %d initial links, %d events, link intervals %s (%d censored at %d)"
          % (args.out, args.out, info["states"], info["initial_links"], info["events"], info["intervals"], info["censored"], info["last_year"]))
    if args.fit:
        here = os.path.dirname(os.path.abspath(__file__))
        cmd = [sys.executable, os.path.join(here, "tricl_fit.py"), os.path.join(args.out, "model.yaml"), os.path.join(args.out, "events.csv"),
               "--fit", *info["parameters"], "--se", "--maxiter", str(args.maxiter), "--json", os.path.join(args.out, "fit.json")]
        if args.bin:
            cmd += ["--bin", args.bin]
        print(" ".join(cmd), flush=True)
        return subprocess.call(cmd)
    return 0


if __name__ == "__main__":
    sys.exit(main())

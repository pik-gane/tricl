#!/usr/bin/env python3
"""Test suite for tricl.

Runs the example configs with fixed seeds and reduced event limits (smoke tests),
compares their final state with stored reference values (regression tests),
checks the log-likelihood against an exact formula for a tiny model,
and checks the error handling of the config parser.

Usage:
    python3 tests/run_tests.py [--bin PATH] [--update] [--only SUBSTRING] [--keep]

--bin     path to the tricl binary (default: build/src/tricl relative to the repository root)
--update  regenerate tests/reference.json from the current binary instead of comparing
--only    run only tests whose name contains SUBSTRING
--keep    keep the temporary working directories

Only the Python standard library is needed. Reference values were generated with gcc/libstdc++
on Linux; other standard libraries implement the random distributions differently, so the
regression tests may fail there although the simulator is correct.
"""

import argparse
import csv
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REFERENCE_FILE = os.path.join(ROOT, "tests", "reference.json")

# (config file relative to the repository root, max. no. of events for the test, random seed)
REGRESSION_CASES = [
    ("config_files/parameters_3blocks.yaml", 2000, 1),
    ("config_files/si.yaml", 5000, 1),
    ("config_files/sir.yaml", 20000, 1),
    ("config_files/sir_sd.yaml", 10000, 1),
    ("config_files/granovetter_simple.yaml", 15001, 1),
    ("config_files/granovetter_helfmann.yaml", 20000, 1),
    ("tests/configs/two_entities.yaml", 100000, 1),
    ("tests/configs/three_entities.yaml", 100000, 1),
]

# tolerances for comparing floating point results with the references:
REL_TOL = 1e-9


class Failure(Exception):
    pass


def limited_copy(src, dst, max_events):
    """Copy a config file, setting limits:events to max_events."""
    with open(src) as f:
        lines = f.read().split("\n")
    out = []
    i = 0
    n = len(lines)
    done = False
    while i < n:
        line = lines[i]
        out.append(line)
        i += 1
        if re.match(r"^limits\s*:", line):
            indent = None
            replaced = False
            # the block consists of the following indented (or empty/comment) lines:
            while i < n and (lines[i].strip() == "" or lines[i].lstrip().startswith("#") or lines[i][0] in " \t"):
                block_line = lines[i]
                i += 1
                m = re.match(r"^(\s+)events\s*:", block_line)
                if m:
                    out.append(m.group(1) + "events: " + str(max_events))
                    replaced = True
                    continue
                m = re.match(r"^(\s+)\S", block_line)
                if m and indent is None:
                    indent = m.group(1)
                out.append(block_line)
            if not replaced:
                out.append((indent or "    ") + "events: " + str(max_events))
            done = True
    if not done:
        raise Failure("no 'limits:' section found in " + src)
    with open(dst, "w") as f:
        f.write("\n".join(out))


def run_tricl(binary, config, args, cwd, timeout=900):
    """Run tricl and return (returncode, stdout, stderr, seconds)."""
    t0 = time.time()
    p = subprocess.run([binary, config] + args, cwd=cwd, capture_output=True, text=True, timeout=timeout)
    return p.returncode, p.stdout, p.stderr, time.time() - t0


def parse_summary(stdout):
    """Return the JSON summary that tricl --summary prints as the last line."""
    for line in reversed(stdout.strip().split("\n")):
        if line.startswith("{"):
            return json.loads(line)
    raise Failure("no JSON summary found in output:\n" + stdout[-2000:])


def close(a, b):
    return abs(a - b) <= REL_TOL * max(1.0, abs(a), abs(b))


def check_summary(summary, name):
    if not math.isfinite(summary["logl"]):
        raise Failure("%s: log-likelihood is not finite: %r" % (name, summary["logl"]))
    if summary["total_rate"] < -1e-6:  # tiny negative values are floating point drift
        raise Failure("%s: negative total rate %r" % (name, summary["total_rate"]))


def regression_test(binary, case, references, update, workdir):
    config, max_events, seed = case
    name = "%s (events<=%d, seed %d)" % (config, max_events, seed)
    src = os.path.join(ROOT, config)
    cwd = tempfile.mkdtemp(prefix="tricl_", dir=workdir)
    dst = os.path.join(cwd, os.path.basename(config))
    limited_copy(src, dst, max_events)
    code, out, err, secs = run_tricl(binary, dst, ["--summary", "--seed", str(seed)], cwd)
    if code != 0:
        raise Failure("%s: exit code %d\nstderr:\n%s" % (name, code, err[-2000:]))
    if "WARNING" in err:
        raise Failure("%s: warning on stderr:\n%s" % (name, err[-2000:]))
    summary = parse_summary(out)
    check_summary(summary, name)
    key = "%s|%d|%d" % (config, max_events, seed)
    if update:
        references[key] = summary
        return "updated (%d events, logl %.6f, %.1fs)" % (summary["events"], summary["logl"], secs)
    if key not in references:
        raise Failure("%s: no reference value stored; run with --update" % name)
    ref = references[key]
    for k in ("events", "links", "angles", "seed"):
        if summary[k] != ref[k]:
            raise Failure("%s: %s = %r differs from reference %r" % (name, k, summary[k], ref[k]))
    for k in ("t", "logl", "total_rate"):
        if not close(summary[k], ref[k]):
            raise Failure("%s: %s = %r differs from reference %r" % (name, k, summary[k], ref[k]))
    return "ok (%d events, logl %.6f, %.1fs)" % (summary["events"], summary["logl"], secs)


def exact_logl_test(binary, workdir):
    """Compare the log-likelihood of tests/configs/two_entities.yaml with the exact formula."""
    config = os.path.join(ROOT, "tests", "configs", "two_entities.yaml")
    cwd = tempfile.mkdtemp(prefix="tricl_", dir=workdir)
    events_file = os.path.join(cwd, "events.csv")
    max_t = 50.0
    a_rate = 0.7          # establishment attempt rate per directed pair, success probability 1
    d_rate = 1.3 * 0.5    # termination attempt rate per directed link times expit(0)
    results = []
    for seed in (1, 2, 3):
        code, out, err, secs = run_tricl(binary, config, ["--summary", "--seed", str(seed), "--events-out", events_file], cwd)
        if code != 0:
            raise Failure("exact logl test: exit code %d\n%s" % (code, err[-2000:]))
        summary = parse_summary(out)
        with open(events_file) as f:
            rows = [line.strip().split(",") for line in f.read().strip().split("\n")]
        if rows[0] != ["t", "event", "source", "relationship", "target"]:
            raise Failure("exact logl test: unexpected events csv header %r" % rows[0])
        linked = False
        t_prev = 0.0
        logl = 0.0
        for row in rows[1:]:
            t, event = float(row[0]), row[1]
            total = 2 * d_rate if linked else 2 * a_rate
            rate = d_rate if event == "terminate" else a_rate
            if (event == "establish") == linked:
                raise Failure("exact logl test: event %r in state linked=%r" % (row, linked))
            logl += math.log(rate) - total * (t - t_prev)
            t_prev = t
            linked = (event == "establish")
        total = 2 * d_rate if linked else 2 * a_rate
        logl -= total * (max_t - t_prev)
        if summary["events"] != len(rows) - 1:
            raise Failure("exact logl test: %d events in summary but %d in csv" % (summary["events"], len(rows) - 1))
        if not abs(summary["logl"] - logl) <= 1e-8 * max(1.0, abs(logl)):
            raise Failure("exact logl test (seed %d): tricl reports logl %.12f, exact value is %.12f" % (seed, summary["logl"], logl))
        results.append("seed %d: %d events, logl %.6f" % (seed, summary["events"], logl))
    return "ok (" + "; ".join(results) + ")"


def gradient_test(binary, workdir):
    """Check --events-in replay against the simulation and --grad against finite differences."""
    config = os.path.join(ROOT, "tests", "configs", "three_entities.yaml")
    cwd = tempfile.mkdtemp(prefix="tricl_", dir=workdir)
    events_file = os.path.join(cwd, "events.csv")
    values = {"A": 0.5, "K": 0.8, "BE": 0.2, "BEK": 0.5, "D": 0.7, "B0": 0.3, "BK": -0.4}
    est, term = "establish that thing knows thing", "terminate that thing knows thing"
    labels = {"A": est + " | base attempt", "K": est + " | attempt via knows thing knows",
              "BE": est + " | base probunits", "BEK": est + " | probunits via knows thing knows",
              "D": term + " | base attempt", "B0": term + " | base probunits", "BK": term + " | probunits via knows thing knows"}
    # parameter dump:
    code, out, err, secs = run_tricl(binary, config, ["--dump-parameters"], cwd)
    if code != 0:
        raise Failure("--dump-parameters: exit code %d\n%s" % (code, err[-1000:]))
    dumped = parse_summary(out)["parameters"]
    for m, label in labels.items():
        if label not in dumped or not close(dumped[label], values[m]):
            raise Failure("--dump-parameters: expected %r = %r, got %r" % (label, values[m], dumped.get(label)))
    # simulate with gradient, writing the events:
    code, out, err, secs = run_tricl(binary, config, ["--summary", "--grad", "--seed", "1", "--events-out", events_file], cwd)
    if code != 0:
        raise Failure("simulation with --grad: exit code %d\n%s" % (code, err[-1000:]))
    sim = parse_summary(out)
    if sim["events"] < 20:
        raise Failure("simulation produced only %d events" % sim["events"])
    # replay the same events with gradient:
    code, out, err, secs = run_tricl(binary, config, ["--summary", "--grad", "--events-in", events_file], cwd)
    if code != 0:
        raise Failure("replay with --grad: exit code %d\n%s" % (code, err[-1000:]))
    rep = parse_summary(out)
    if rep["events"] != sim["events"]:
        raise Failure("replay performed %d events, simulation %d" % (rep["events"], sim["events"]))
    if not close(rep["logl"], sim["logl"]):
        raise Failure("replay logl %.12f differs from simulation logl %.12f" % (rep["logl"], sim["logl"]))
    for label in labels.values():
        if label not in rep["gradient"]:
            raise Failure("gradient component %r missing" % label)
        if not close(rep["gradient"][label], sim["gradient"][label]):
            raise Failure("replay gradient %r = %r differs from simulation %r" % (label, rep["gradient"][label], sim["gradient"][label]))
    # replaying must not depend on the seed (no random numbers are drawn in replay mode):
    code, out, err, secs = run_tricl(binary, config, ["--summary", "--grad", "--events-in", events_file, "--seed", "987654"], cwd)
    if code != 0:
        raise Failure("replay with another seed: exit code %d\n%s" % (code, err[-1000:]))
    rep2 = parse_summary(out)
    if rep2["logl"] != rep["logl"] or rep2["gradient"] != rep["gradient"]:
        raise Failure("replay results depend on the seed: logl %r vs %r" % (rep2["logl"], rep["logl"]))
    # finite differences of the replayed log-likelihood w.r.t. each metaparameter:
    results = []
    for m, label in labels.items():
        v = values[m]
        h = 1e-4 * max(1.0, abs(v))
        logls = []
        for x in (v + h, v - h):
            code, out, err, secs = run_tricl(binary, config, ["--summary", "--events-in", events_file, "--" + m, repr(x)], cwd)
            if code != 0:
                raise Failure("replay with --%s %r: exit code %d\n%s" % (m, x, code, err[-1000:]))
            logls.append(parse_summary(out)["logl"])
        fd = (logls[0] - logls[1]) / (2 * h)
        an = rep["gradient"][label]
        if not abs(fd - an) <= 1e-5 * max(1.0, abs(an)):
            raise Failure("d logl / d %s: analytic %.10f, finite difference %.10f" % (m, an, fd))
        results.append("%s: %.4f" % (m, an))
    return "ok (%d events, logl %.4f, gradient %s)" % (sim["events"], sim["logl"], ", ".join(results))


def replay_consistency_test(binary, workdir):
    """Simulate example configs with random initial links, then replay their events: the log-likelihoods must agree."""
    results = []
    for config, max_events, seed in [("config_files/si.yaml", 5000, 3), ("config_files/sir_sd.yaml", 100000, 5)]:
        cwd = tempfile.mkdtemp(prefix="tricl_", dir=workdir)
        dst = os.path.join(cwd, os.path.basename(config))
        limited_copy(os.path.join(ROOT, config), dst, max_events)
        events_file = os.path.join(cwd, "events.csv")
        code, out, err, secs = run_tricl(binary, dst, ["--summary", "--seed", str(seed), "--events-out", events_file], cwd)
        if code != 0:
            raise Failure("%s: simulation exit code %d\n%s" % (config, code, err[-1000:]))
        sim = parse_summary(out)
        code, out, err, secs = run_tricl(binary, dst, ["--summary", "--seed", str(seed), "--events-in", events_file], cwd)
        if code != 0:
            raise Failure("%s: replay exit code %d\n%s" % (config, code, err[-1000:]))
        rep = parse_summary(out)
        if rep["events"] != sim["events"] or rep["links"] != sim["links"] or rep["angles"] != sim["angles"]:
            raise Failure("%s: replay ended with %d events, %d links, %d angles; simulation with %d, %d, %d"
                          % (config, rep["events"], rep["links"], rep["angles"], sim["events"], sim["links"], sim["angles"]))
        if not close(rep["logl"], sim["logl"]):
            raise Failure("%s: replay logl %.12f differs from simulation logl %.12f" % (config, rep["logl"], sim["logl"]))
        results.append("%s: %d events, logl %.3f, replay %.2fs" % (os.path.basename(config), sim["events"], sim["logl"], secs))
    return "ok (" + "; ".join(results) + ")"


def fit_test(binary, workdir):
    """Run the Python maximum-likelihood driver on simulated data of the three-entity model."""
    sys.path.insert(0, os.path.join(ROOT, "python"))
    from tricl_fit import TriclModel, fit, standard_errors
    config = os.path.join(ROOT, "tests", "configs", "three_entities.yaml")
    cwd = tempfile.mkdtemp(prefix="tricl_", dir=workdir)
    events_file = os.path.join(cwd, "events.csv")
    # longer observation window for more events:
    code, out, err, secs = run_tricl(binary, config, ["--summary", "--seed", "7", "--events-out", events_file], cwd)
    if code != 0:
        raise Failure("simulation: exit code %d\n%s" % (code, err[-1000:]))
    model = TriclModel(binary, config, events_file, quiet=True)
    names = ["A", "D"]
    start = {"A": 0.5 * 1.8, "D": 0.7 / 1.8}
    logl0, _ = model.loglikelihood(start, gradient=False)
    estimates, logl, grad, n_iter = fit(model, start, names, maxiter=50, quiet=True)
    if not logl > logl0:
        raise Failure("fit did not improve the log-likelihood (%.6f -> %.6f)" % (logl0, logl))
    gnorm = math.sqrt(sum(g * g for g in grad.values()))
    if not gnorm < 1e-3:
        raise Failure("gradient norm after fitting is %.3g" % gnorm)
    se = standard_errors(model, start, names, estimates)
    for name in names:
        if not (math.isfinite(se[name]) and se[name] > 0):
            raise Failure("standard error of %s is %r" % (name, se[name]))
    return "ok (logl %.3f -> %.3f, %s)" % (logl0, logl, ", ".join("%s=%.3f+-%.3f" % (n, estimates[n], se[n]) for n in names))


def rdf_test(binary, workdir):
    """Round trip: N-Triples -> tricl config skeleton -> simulation -> Turtle-star."""
    cwd = tempfile.mkdtemp(prefix="tricl_", dir=workdir)
    script = os.path.join(ROOT, "python", "tricl_rdf.py")
    data = os.path.join(ROOT, "tests", "data", "people.nt")
    p = subprocess.run([sys.executable, script, "rdf2tricl", data, "--out", cwd, "--name", "people"], capture_output=True, text=True)
    if p.returncode != 0:
        raise Failure("rdf2tricl failed:\n%s" % p.stderr[-1000:])
    with open(os.path.join(cwd, "model.yaml")) as f:
        model = f.read()
    for expected in ['"Person":', '"Organisation":', '- "Alice"', '- "Carol, the third"', '- "dave"', '- "globex"',
                     '"knows": symmetric', '"worksFor": "employs"', 'type: ["Person", "knows", "Person"]',
                     'type: ["Person", "worksFor", "Organisation"]']:
        if expected not in model:
            raise Failure("rdf2tricl: %r missing from model.yaml:\n%s" % (expected, model))
    with open(os.path.join(cwd, "links_Person_knows_Person.csv")) as f:
        rows = [r for r in csv.reader(f)][1:]
    if sorted(rows) != [["Alice", "Bob"], ["Bob", "Carol, the third"], ["Carol, the third", "dave"]]:
        raise Failure("rdf2tricl: unexpected knows links %r (symmetric duplicates and self-links must be dropped)" % rows)
    with open(os.path.join(cwd, "links_Person_worksFor_Organisation.csv")) as f:
        rows = [r for r in csv.reader(f)][1:]
    if sorted(rows) != [["Alice", "ACME"], ["Bob", "ACME"], ["Carol, the third", "globex"]]:
        raise Failure("rdf2tricl: unexpected worksFor links %r (inverse 'employs' must be folded in)" % rows)
    # add dynamics and run tricl on the generated config:
    model = model.replace("    t: 1  # TODO: simulation time", "    t: 5").replace("dynamics: {}", "").rstrip()
    model += """
files:
    gexf: people.gexf.gz
dynamics:
    ["Person", "knows", "Person"]:
        establish:
            attempt:
                base: 0.2
                [~, "knows", "Person", "knows", ~]: 1.0
            success: inf
        terminate:
            attempt: 0.3
            success: inf
    ["Person", "worksFor", "Organisation"]:
        terminate:
            attempt: 0.1
            success: inf
"""
    config = os.path.join(cwd, "model_with_dynamics.yaml")
    with open(config, "w") as f:
        f.write(model)
    code, out, err, secs = run_tricl(binary, config, ["--summary", "--seed", "1"], cwd)
    if code != 0:
        raise Failure("tricl on generated config: exit code %d\n%s" % (code, err[-1000:]))
    summary = parse_summary(out)
    ttl = os.path.join(cwd, "people.ttl")
    p = subprocess.run([sys.executable, script, "tricl2rdf", os.path.join(cwd, "people.gexf.gz"), "--out", ttl], capture_output=True, text=True)
    if p.returncode != 0:
        raise Failure("tricl2rdf failed:\n%s" % p.stderr[-1000:])
    with open(ttl) as f:
        turtle = f.read()
    n_intervals = turtle.count("tricl:from")
    if ":Alice a :Person" not in turtle or "<< :Alice :worksFor :ACME >>" not in turtle or n_intervals < 6:
        raise Failure("tricl2rdf: unexpected output:\n%s" % turtle[:1500])
    return "ok (%d events simulated, %d link intervals exported)" % (summary["events"], n_intervals)


def macro_test(binary, workdir):
    """The macroscopic approximation must be exact for a linear model, and must run on an example with angles."""
    sys.path.insert(0, os.path.join(ROOT, "python"))
    from tricl_macro import MacroModel, dump_model, simulate_stats, with_time_limit
    cwd = tempfile.mkdtemp(prefix="tricl_", dir=workdir)
    config = os.path.join(ROOT, "tests", "configs", "mean_field.yaml")
    model = dump_model(binary, config, {}, 1)
    macro = MacroModel(model, warn=lambda s: None)
    times, states = macro.integrate(10.0, 1.0)
    i = macro.index[("agent", "is", "active")]
    # exact solution of dn/dt = 0.3 (400 - n) - 0.2 n with n(0) = 0:
    exact = [400 * 0.3 / 0.5 * (1 - math.exp(-0.5 * t)) for t in times]
    for t, n, x in zip(times, states, exact):
        if abs(n[i] - x) > 1e-3 * max(1.0, x):
            raise Failure("ODE integration: n(%g) = %.6f, exact %.6f" % (t, n[i], x))
    sim = simulate_stats(binary, with_time_limit(config, 10.0), [1, 2, 3], 10.0, 1.0, {}, cwd)
    label = "agent is active"
    final = [sim[s][10][label] for s in sim]
    mean = sum(final) / len(final)
    if abs(mean - exact[-1]) > 0.1 * exact[-1]:
        raise Failure("simulation mean %.1f at t=10 differs from mean field %.1f by more than 10%%" % (mean, exact[-1]))
    # smoke test on an example with angle influences (SIS on a random geometric graph):
    model = dump_model(binary, os.path.join(ROOT, "config_files", "si.yaml"), {}, 1)
    macro = MacroModel(model, warn=lambda s: None)
    times, states = macro.integrate(20.0, 2.0)
    j = macro.index[("agent", "is", "infected")]
    if not all(math.isfinite(n[j]) and 0 <= n[j] <= 1000 for n in states) or states[-1][j] <= states[0][j]:
        raise Failure("approximation of si.yaml gives implausible infected counts %r" % [n[j] for n in states])
    si_final = states[-1][j]
    # threshold model with mutually exclusive states (immediate events handled structurally):
    config = os.path.join(ROOT, "config_files", "granovetter_helfmann.yaml")
    model = dump_model(binary, config, {}, 1)
    macro = MacroModel(model, warn=lambda s: None)
    if macro.stiff:
        raise Failure("the immediate events of granovetter_helfmann.yaml were not recognised as state switches")
    times, states = macro.integrate(30.0, 5.0)
    k = macro.index[("c-agent", "is", "active")]
    sim = simulate_stats(binary, with_time_limit(config, 30.0), [1, 2, 3], 30.0, 5.0, {}, cwd)
    label = "c-agent is active"
    final = [sim[s][6][label] for s in sim]
    mean_g = sum(final) / len(final)
    if not (math.isfinite(states[-1][k]) and abs(states[-1][k] - mean_g) <= 0.25 * 60):
        raise Failure("granovetter: approximation %.1f vs simulation mean %.1f active c-agents at t=30" % (states[-1][k], mean_g))
    return "ok (linear model: sim. mean %.1f vs mean field %.1f at t=10; si.yaml: %.0f infected at t=20; granovetter: %.1f vs sim. %.1f active at t=30)" % (
        mean, exact[-1], si_final, states[-1][k], mean_g)


def sigmoid_unit_test(binary, workdir):
    """Run the C++ unit tests of the sigmoidal function (built as build/tests/test_sigmoid)."""
    build_dir = os.path.dirname(os.path.dirname(binary))
    test_binary = os.path.join(build_dir, "tests", "test_sigmoid")
    if not os.path.exists(test_binary):
        raise Failure("test binary %s not found (configure and build the whole project)" % test_binary)
    p = subprocess.run([test_binary], capture_output=True, text=True)
    if p.returncode != 0:
        raise Failure("test_sigmoid failed:\n%s" % (p.stdout + p.stderr)[-3000:])
    return "ok (" + p.stdout.strip().split("\n")[-1] + ")"


def error_handling_tests(binary, workdir):
    cwd = tempfile.mkdtemp(prefix="tricl_", dir=workdir)
    results = []
    # undefined symbol in an expression:
    config = os.path.join(ROOT, "tests", "configs", "bad_expression.yaml")
    code, out, err, secs = run_tricl(binary, config, ["--summary"], cwd)
    if code == 0 or "cannot parse expression" not in err or "undefined_symbol" not in err:
        raise Failure("bad expression: expected exit code != 0 and a parse error naming the expression, got %d:\n%s" % (code, err[-1000:]))
    results.append("bad expression rejected")
    # missing csv file:
    config = os.path.join(ROOT, "tests", "configs", "missing_csv.yaml")
    code, out, err, secs = run_tricl(binary, config, ["--summary"], cwd)
    if code == 0 or "this_file_does_not_exist.csv" not in err:
        raise Failure("missing csv: expected exit code != 0 and an error naming the file, got %d:\n%s" % (code, err[-1000:]))
    results.append("missing csv reported")
    # unknown entity type in the draft config:
    config = os.path.join(ROOT, "config_files", "drafts", "identity.yaml")
    code, out, err, secs = run_tricl(binary, config, ["--summary"], cwd)
    if code == 0 or "unknown" not in err:
        raise Failure("draft config: expected exit code != 0 and an 'unknown ...' error, got %d:\n%s" % (code, err[-1000:]))
    results.append("undeclared name reported")
    # one-letter metaparameter given as --X on the command line:
    config = os.path.join(ROOT, "config_files", "granovetter_helfmann.yaml")
    code, out, err, secs = run_tricl(binary, config, ["--summary", "--seed", "1", "--E", "3"], cwd)
    if code != 0:
        raise Failure("--E 3: exit code %d\n%s" % (code, err[-1000:]))
    summary = parse_summary(out)
    if summary["events"] != 6:
        raise Failure("--E 3 should limit the run to 2*3 = 6 events, got %d" % summary["events"])
    results.append("--E works")
    # YAML's .inf literal and an infinite time limit:
    config = os.path.join(ROOT, "config_files", "si.yaml")
    code, out, err, secs = run_tricl(binary, config, ["--summary", "--seed", "1"], cwd)
    if code != 0:
        raise Failure("si.yaml (uses .inf): exit code %d\n%s" % (code, err[-1000:]))
    results.append(".inf accepted")
    return "ok (" + ", ".join(results) + ")"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", default=os.path.join(ROOT, "build", "src", "tricl"))
    ap.add_argument("--update", action="store_true")
    ap.add_argument("--only", default="")
    ap.add_argument("--keep", action="store_true")
    args = ap.parse_args()
    binary = os.path.abspath(args.bin)
    if not os.path.exists(binary):
        print("tricl binary not found at %s (build it first or pass --bin)" % binary)
        return 2

    references = {}
    if os.path.exists(REFERENCE_FILE):
        with open(REFERENCE_FILE) as f:
            references = json.load(f)

    workdir = tempfile.mkdtemp(prefix="tricl_tests_")
    tests = [("sigmoid unit tests", lambda: sigmoid_unit_test(binary, workdir))]
    for case in REGRESSION_CASES:
        tests.append(("regression: " + case[0], lambda case=case: regression_test(binary, case, references, args.update, workdir)))
    tests.append(("exact log-likelihood", lambda: exact_logl_test(binary, workdir)))
    tests.append(("replay and gradient", lambda: gradient_test(binary, workdir)))
    tests.append(("replay of example configs", lambda: replay_consistency_test(binary, workdir)))
    tests.append(("maximum-likelihood fit", lambda: fit_test(binary, workdir)))
    tests.append(("rdf import and export", lambda: rdf_test(binary, workdir)))
    tests.append(("macroscopic approximation", lambda: macro_test(binary, workdir)))
    tests.append(("error handling", lambda: error_handling_tests(binary, workdir)))

    n_failed = 0
    for name, fn in tests:
        if args.only and args.only not in name:
            continue
        try:
            print("%-55s %s" % (name, fn()), flush=True)
        except Failure as e:
            n_failed += 1
            print("%-55s FAILED\n    %s" % (name, str(e).replace("\n", "\n    ")), flush=True)
        except Exception as e:
            n_failed += 1
            print("%-55s ERROR: %r" % (name, e), flush=True)

    if args.update:
        with open(REFERENCE_FILE, "w") as f:
            json.dump(references, f, indent=1, sort_keys=True)
        print("reference values written to", REFERENCE_FILE)
    if not args.keep:
        shutil.rmtree(workdir, ignore_errors=True)
    else:
        print("working directories kept in", workdir)
    print("%d test(s) failed" % n_failed if n_failed else "all tests passed")
    return 1 if n_failed else 0


if __name__ == "__main__":
    sys.exit(main())

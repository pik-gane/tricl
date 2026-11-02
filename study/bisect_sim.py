#!/usr/bin/env python3
"""Largest upper bound on a parameter (with other fixed bounds) for which the constrained maximum-likelihood fit is
stable in simulation (days 3 and 4 of the timetable, no run exceeding 3 x the observed maximum of simultaneous
contacts). Usage: python3 study/bisect_sim.py MODEL PARAM LO HI [FIXED_BOUND ...]  (e.g. F_tri=0.473)"""
import json, os, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import stability as st
model, par, lo, hi = int(sys.argv[1]), sys.argv[2], float(sys.argv[3]), float(sys.argv[4])
fixed = sys.argv[5:]
log = []
def fit_and_check(c):
    tag = "_" + "_".join([b.replace("=", "") for b in fixed] + ["%s%.3f" % (par, c)])
    out = os.path.join(HERE, "results", "fit_M%d%s.json" % (model, tag))
    if not os.path.exists(out):
        subprocess.run([sys.executable, os.path.join(HERE, "fit.py"), str(model), "--tag", tag, "--bound", "%s=%r" % (par, c)] + fixed
                       + ["--start", os.path.join(HERE, "results", "fit_M%d.json" % model)], check=True)
    r = json.load(open(out))
    worst, tipped = st.simulated_tipping(os.path.join(HERE, "configs", "M%d.yaml" % model), r["estimates"], days=(3, 4), seeds=(1,), wall=150)
    log.append({"bound": c, "fixed": fixed, "worst_ratio": worst, "stable": not tipped, "logl_train": r["logl_train"], "logl_test": r["logl_test"], "estimates": r["estimates"]})
    print("bound %s <= %.3f (%s): max/observed %.2f, stable %s, logl train %.1f test %.1f" % (par, c, " ".join(fixed), worst, not tipped, r["logl_train"], r["logl_test"]), flush=True)
    json.dump({"model": model, "param": par, "log": log}, open(os.path.join(HERE, "results", "stable_sim_M%d_%s.json" % (model, par)), "w"), indent=1)
    return not tipped
if not fit_and_check(lo):
    print("lower bound not stable"); sys.exit(1)
for _ in range(5):
    mid = (lo + hi) / 2
    if fit_and_check(mid): lo = mid
    else: hi = mid
print("largest stable bound %s <= %.3f" % (par, lo))

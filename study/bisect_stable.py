#!/usr/bin/env python3
"""Largest upper bound on a destabilizing parameter for which the constrained maximum-likelihood fit is stable
(pair closure: no class tips within three days of the timetable). Usage: python3 study/bisect_stable.py MODEL PARAM LO HI"""
import json, os, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import stability as st
model, par, lo, hi = int(sys.argv[1]), sys.argv[2], float(sys.argv[3]), float(sys.argv[4])
data = json.load(open(os.path.join(HERE, "data", "data.json"))); obs = st.observed_maxima(data)
log = []
def fit_and_check(c):
    tag = "_%s%.3f" % (par, c)
    subprocess.run([sys.executable, os.path.join(HERE, "fit.py"), str(model), "--tag", tag, "--bound", "%s=%r" % (par, c),
                    "--start", os.path.join(HERE, "results", "fit_M%d.json" % model)], check=True)
    r = json.load(open(os.path.join(HERE, "results", "fit_M%d%s.json" % (model, tag))))
    tip = st.tipping_time(r["estimates"], data, obs)
    stable = tip >= st.REPEAT * 540.0
    log.append({"bound": c, "tip": tip, "stable": stable, "logl_train": r["logl_train"], "logl_test": r["logl_test"], "estimates": r["estimates"]})
    print("bound %s <= %.3f: tip %.0f min, stable %s, logl train %.1f test %.1f" % (par, c, tip, stable, r["logl_train"], r["logl_test"]), flush=True)
    return stable
# lo must be stable, hi unstable:
if not fit_and_check(lo):
    print("lower bound not stable"); sys.exit(1)
for _ in range(5):
    mid = (lo + hi) / 2
    if fit_and_check(mid): lo = mid
    else: hi = mid
json.dump({"model": model, "param": par, "largest_stable_bound": lo, "log": log}, open(os.path.join(HERE, "results", "stable_M%d_%s.json" % (model, par)), "w"), indent=1)
print("largest stable bound %s <= %.3f" % (par, lo))

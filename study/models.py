#!/usr/bin/env python3
"""Model configs M0-M3 for the contact study (see PAPER_PLAN.md), written to study/configs/.

All rates are per minute. Contacts are links of the symmetric relationship type "meets". Formation: attempts at the
fixed ceiling rate A_FORM per direction, succeeding with probability expit(4 B), where B is a sum of probability
units (log-linear rates as long as the probability is small). Termination: attempts at the fixed ceiling D_TERM,
success expit(4 B'). Every influence enters both directions of a pair with the same coefficient, so rates are
invariant under reversal (the log-likelihood with event directions differs from the one without by a constant).

M0: timetable only (hour of the day, the first 20 minutes of each hour, the morning break)
M1: + static context (same class, Facebook friends, survey friendship in either direction)
M2: + dynamic network context (common current contacts, current numbers of contacts of both students)
M4 ("M1m"): M1 + the memory terms of M3, without network context
M3: + Markovian memory (met on an earlier day; decaying impacts of the starts and ends of contacts of the pair,
      on two time scales)
"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
A_FORM = 1.0   # formation attempt rate per direction (per minute), the ceiling of the formation rate
D_TERM = 6.0   # termination attempt rate per direction (per minute)
HOURS = ["h%02d" % h for h in range(8, 17)]
DECAYS = {"fast": 0.5, "slow": 0.05}  # decay rates (per minute) of the memory impacts (models 3, 4)
DECAYS4 = {"d2": 2.0, "d05": 0.5, "d01": 0.1, "d002": 0.02}  # four time scales (models 5, 6)


def decays(model):
    return DECAYS4 if model in (5, 6, 8) else DECAYS


ATTEMPT_PARAMS = ("F_tria",)  # parameters that are attempt rates (non-negative)

# initial values of the metaparameters (probability units), refined by the fits:
START = {"F0": -2.0, "T0": -0.6}


def params(model):
    """Metaparameters of a model: name -> (event, cause or 'base', start value)."""
    p = {"F0": ("form", "base", START["F0"]), "T0": ("term", "base", START["T0"])}
    for h in HOURS[1:]:
        p["F_" + h] = ("form", ("leg", h), 0.01)
        p["T_" + h] = ("term", ("leg", h), 0.01)
    for k, hub in (("chg", "hour change"), ("brk", "morning break")):
        p["F_" + k] = ("form", ("leg", hub), 0.01)
        p["T_" + k] = ("term", ("leg", hub), 0.01)
    if model >= 1:
        for k, angles in (("cls", [("in class", "class", "has student")]),
                          ("fb", [("=", "student", "facebook")]),
                          ("svy", [("=", "student", "names"), ("=", "student", "named by")])):
            p["F_" + k] = ("form", ("angles", angles), 0.01)
            p["T_" + k] = ("term", ("angles", angles), 0.01)
    net = model in (2, 3, 5)    # network context
    mem = model in (3, 4, 5, 6, 8)    # memory (model 4: M1 plus memory without network context; 5, 6: as 3, 4 with four time scales)
    if model in (7, 8):
        # saturating network context: common contacts raise the attempt rate additively (per common contact and
        # direction), the numbers of current contacts enter the success probability units (a capacity effect if negative):
        p["F_tria"] = ("form", ("attempt_angles", [("meets", "student", "meets")]), 0.1)
        p["T_tri"] = ("term", ("angles", [("meets", "student", "meets")]), 0.01)
        p["F_deg"] = ("form", ("degree", None), 0.01)
        p["T_deg"] = ("term", ("degree", None), 0.01)
    if net:
        p["F_tri"] = ("form", ("angles", [("meets", "student", "meets")]), 0.01)
        p["T_tri"] = ("term", ("angles", [("meets", "student", "meets")]), 0.01)
        p["F_deg"] = ("form", ("degree", None), 0.01)
        p["T_deg"] = ("term", ("degree", None), 0.01)
    if mem:
        p["F_mb"] = ("form", ("angles", [("=", "student", "met before")]), 0.01)
        p["T_mb"] = ("term", ("angles", [("=", "student", "met before")]), 0.01)
        for s in decays(model):
            p["F_end_" + s] = ("form", ("angles", [("ended " + s, "student", "=")]), 0.01)
            p["T_start_" + s] = ("term", ("angles", [("started " + s, "student", "=")]), 0.01)
    return p


def config(model, data, values=None):
    p = params(model)
    mem = model in (3, 4, 5, 6, 8)
    vals = {k: (values or {}).get(k, v[2]) for k, v in p.items()}
    L = ["metadata:", "    name: contact study model M%d" % model, "", "limits:", "    t: 540", "", "options:", "    seed: 1", "", "metaparameters:"]
    L += ["    %s: %r" % (k, float(v)) for k, v in vals.items()]
    L += ["", "entities:", "    student:"] + ["        - %s" % s for s in data["students"]]
    L += ["    class:"] + ["        - %s" % c for c in data["info"]["classes"]]
    for hub in HOURS + ["hour change", "morning break"]:
        L += ["    %s:" % hub, "        - %s" % hub]
    L += ["", "relationship types:", "    meets: symmetric", "    in class: has student", "    facebook: symmetric",
          "    names: named by", "    met before: symmetric", "    is in: ~"]
    if mem:
        L += ["", "action types:"]
        for s, d in decays(model).items():
            L += ["    started %s:" % s, "        inverse: symmetric", "        decay: %r" % d, "        threshold: 0.01",
                  "        triggered by: [meets, establish]"]
            L += ["    ended %s:" % s, "        inverse: symmetric", "        decay: %r" % d, "        threshold: 0.01",
                  "        triggered by: [meets, terminate]"]
    L += ["", "initial links:", "    named:"]
    L += ["        - [%s, in class, %s]" % (s, c) for s, c in data["static"]["class"]]
    L += ["        - [%s, facebook, %s]" % (a, b) for a, b in data["static"]["facebook"]]
    L += ["        - [%s, names, %s]" % (a, b) for a, b in data["static"]["survey"]]
    L += ["", "dynamics:", "    [student, meets, student]:"]
    for ev, ceiling in (("form", A_FORM), ("term", D_TERM)):
        L += ["        %s:" % ("establish" if ev == "form" else "terminate"), "            attempt:", "                base: %r" % ceiling]
        for k, (e, cause, _) in p.items():
            if (e == ev) and (cause != "base") and (cause[0] == "attempt_angles"):
                for r12, et2, r23 in cause[1]:
                    L.append("                [~, %s, %s, %s, ~]: %s" % (r12, et2, r23, k))
        L += ["            success:", "                tails: 0"]
        for k, (e, cause, _) in p.items():
            if (e != ev) or ((cause != "base") and (cause[0] == "attempt_angles")):
                continue
            if cause == "base":
                L.append("                base: %s" % k)
            elif cause[0] == "leg":
                L.append("                [~, is in, %s]: %s" % (cause[1], k))
            elif cause[0] == "degree":  # current contacts of both students (not counting this pair's own link)
                L.append("                [~, meets, student]: %s" % k)
                L.append("                [student, meets, ~]: %s" % k)
            else:
                for r12, et2, r23 in cause[1]:
                    L.append("                [~, %s, %s, %s, ~]: %s" % (r12, et2, r23, k))
    return "\n".join(L) + "\n"


def main():
    data = json.load(open(os.path.join(HERE, "data", "data.json")))
    os.makedirs(os.path.join(HERE, "configs"), exist_ok=True)
    for m in range(9):
        open(os.path.join(HERE, "configs", "M%d.yaml" % m), "w").write(config(m, data))
        print("M%d: %d metaparameters" % (m, len(params(m))))


if __name__ == "__main__":
    main()

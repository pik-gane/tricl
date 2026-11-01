# Plan for a research paper with tricl (draft, 2026-11-01)

Committed at the author's request.

## 1. Why a new paper

The current draft (`paper/`) describes software: model class, algorithms, tests, one application.
A research paper needs a question whose answer is not known in advance, a comparison with an established competing
model, and results that would change what readers believe.
The draft stays as it is; parts of it (model definition, likelihood, appendices) can become the methods section or
supplementary material.

## 2. Research question

**How much of the apparent memory of human contact dynamics is the effect of the surrounding network?**

Face-to-face contact data show broad (roughly power-law) distributions of contact and inter-contact durations.
The established explanation is intrinsic memory: the longer a contact or a group lasts, the less likely it is to end
("aging" or reinforcement models, which are not Markovian in the network state).
The alternative explanation: the dynamics is Markovian, but in a richer state.
A contact embedded in a group (angles), between students of the same class or between friends (static relationships),
or of a student with many current contacts (legs) has different rates, and a mixture of such rates, changing as
groups form and dissolve, also produces broad durations.
Decaying impacts of past interactions (actions) add memory in a Markovian way, in an extended state.

The answer is decidable with tricl's exact likelihood: fit a hierarchy of Markov models and the competing
non-Markovian model to the same event sequences, and compare them by held-out log-likelihood and by how well they
reproduce statistics that were not used in fitting.

Working title: *Network context explains memory in human contact dynamics* (or, if the result goes the other way,
*Intrinsic memory, not network context, shapes ...*; the plan works for either outcome).

## 3. Data

Primary: **SocioPatterns high school 2013** (Lycée Thiers, Marseille; Netzschleuder `sp_high_school`, reachable from
this environment, 0.5 MB).
329 students in 9 classes, face-to-face proximity at 20-second resolution over 5 days (188,508 contact records),
plus three static networks: Facebook friendships, survey friendships, contact diaries; class membership as metadata.
Contacts become links of a symmetric relationship type `meets` (established when a run of 20-s records starts,
terminated after the last one; gaps of one interval merged, a sensitivity analysis on this choice).
Class membership and friendships become static relationship types, so that angles such as
`[~, member of, class, has member, ~]` (same class) and `[~, =, student, is friend of, ~]` (friends) enter the rates.

Second data set for robustness: **SocioPatterns primary school** (`sp_primary_school`, 2 days, 242 children and
teachers, 20-s resolution) or the **baboon sensor data** (`sp_baboons`, 13 animals, about a month), the latter
being a case without school timetables.

Optional extension with observed acts: **Copenhagen Networks Study** (`copenhagen`): Bluetooth co-location
(5-minute scans), calls and SMS (exact timestamps, as acts with decaying impacts), Facebook friendships.
This would test whether calls and messages trigger meetings and vice versa (cross-excitation between a relationship
and actions).
Its 5-minute resolution of co-location makes event times interval-censored, so it is better as a separate,
second paper.

## 4. The Markov question, addressed explicitly

Contacts in a school are not a homogeneous Markov process: the timetable (lessons, breaks, lunch) is an
exogenous driver.
Two measures:
1. **Exogenous events.** The timetable becomes forced events (all students' links `is in` to a hub entity `break`
   switch at the known times); rates depend on them through legs (`[~, is in, break]`).
   Forced events are applied in simulation and replay at their times, without contribution to the likelihood
   (the likelihood is conditional on them).
   This needs a small code extension (Sec. 7).
2. **A test of the Markov property.** By the time-rescaling theorem, if the model is correct, the integrated rates
   (compensators) between successive events of each pair are independent unit exponentials.
   tricl computes the integrated rates exactly (also for time-varying rates), so Kolmogorov–Smirnov plots of the
   rescaled times are a direct goodness-of-fit test for each model, including the competing one.
   This turns "is it Markov?" from an assumption into a result.

## 5. Models to compare

All models are fitted by maximum likelihood on days 1–3 and evaluated on days 4–5 (held-out log-likelihood per event),
with standard errors from the observed information.

* **M0 (memoryless baseline):** constant rates of forming and ending a contact, per class pair (same class or not).
* **M1 (static context):** M0 plus same class, Facebook friendship and survey friendship (angles via `=`), and the
  timetable (leg to `break`).
* **M2 (dynamic network context):** M1 plus triadic closure and group stabilization (angle
  `[~, meets, student, meets, ~]` in formation and ending) and capacity (legs `[~, meets, student]` and
  `[student, meets, ~]`, i.e. the current numbers of contacts of the two students).
* **M3 (Markovian memory):** M2 plus an action type `has met` with decaying impact that encodes recent contact of a
  pair and enters formation and ending rates (reinforcement).
  Open design point: acts must be observed for the likelihood; the natural observed act is the start of a contact,
  which requires an act triggered by an establishment.
  This needs a small extension (an act that happens immediately with a relationship event, Sec. 7);
  otherwise M3 uses contacts of the previous day as a static relationship (weaker, but needs no extension).
* **C1 (competing, intrinsic memory):** the aging/reinforcement model of face-to-face contacts, in which the rate of
  ending a contact or leaving a group decreases with its current age (Stehlé, Barrat and coworkers; Zhao et al. 2011;
  exact references to be checked against Crossref before use).
  Its likelihood is explicit (rates depend on ages, which are known from the data), so it can be fitted and evaluated
  on the same held-out data; it needs its own short Python implementation.
* **C2 (competing, activity-driven with memory):** activity-driven network with memory (Karsai et al. 2014;
  reference to be checked), the standard temporal-network baseline.

Expected outcome if the hypothesis holds: M2 or M3 reach the held-out likelihood of C1 or exceed it, pass the
time-rescaling test where M0 and M1 fail, and reproduce the contact and inter-contact duration distributions and the
group-size distribution without age-dependent rates.
If C1 wins clearly, the paper reports that network context does not explain the memory, which is also publishable;
quantifying which part of the memory context explains (the gap in held-out likelihood between M0, M2 and C1) is the
main result in either case.

## 6. Higher-order closure (mean-field section)

The simple closure of the current draft (independent Poisson angle counts) is trivial for this purpose.
Plan: generate a **pair approximation** automatically from the config file.
State variables are the densities of links and of pairs of adjacent links (angles) by type; three-link motifs are
closed by the standard Kirkwood or Keeling-type approximation, so that angle counts are correlated with link
presence instead of independent.
For legs, the closure carries the degree distribution (a degree-based pair approximation, as in the
adaptive-network literature).
Use in the paper: predict the stationary contact density, the mean group size and the triangle density of the
fitted M2 from the generated equations, and compare with simulation and data.
The comparison of simple closure, pair closure and simulation shows how much of the group structure the closure
captures; the deviation measures higher-order correlations.
This is a methods contribution in its own right (automatic generation of pair closures for typed multi-relational
rules) and gives an analytic understanding of where the broad durations come from.

## 7. Code extensions needed (in order)

1. Data pipeline for SocioPatterns (`python/tricl_sociopatterns.py`, like `tricl_atop.py`): download, runs of 20-s
   records to establish/terminate events, class and friendship links, timetable from the data (breaks are visible as
   bursts of cross-class contacts) or from the published schedule.
2. Forced (exogenous) events: a file of events applied at given times in simulation and replay, excluded from the
   likelihood. Small change in `pop_next_event` (as for expiries) and in replay.
3. Time-rescaling residuals: output of the integrated rate of each pair between its events (`--residuals-out`), for
   goodness-of-fit plots.
4. Optional: acts triggered by relationship events (for M3), e.g. an `on establish:` key listing acts that happen
   with probability 1 at the same time (no likelihood contribution).
5. Pair-closure generator in `python/tricl_macro.py` (`--closure pair`).
6. Competing models C1 and C2: Python likelihoods on the same event files, and simulation for the statistics.
7. Model comparison script: held-out log-likelihood, AIC/BIC, time-rescaling plots, duration and group-size
   distributions, all into `paper2/generated/` with the same conventions as `paper/` (numbers from a generated file,
   one sentence per line, Crossref-checked references).

Items 1–3 and 5–7 are a few days of work each at most; item 4 needs a short design discussion.

## 8. Figures and tables (target: Physical Review E regular article, about 10 pages)

1. Data and the puzzle: contact and inter-contact duration distributions (broad), compared with M0 (exponential).
2. Model hierarchy and fitted effects (table of estimates with standard errors; held-out log-likelihood per event for
   M0–M3, C1, C2).
3. Time-rescaling (KS) plots for M0, M2, M3, C1.
4. Reproduction of statistics not used in fitting: duration distributions, group sizes, triangle density, per model.
5. Mean field: simple vs pair closure vs simulation vs data for the stationary quantities.
6. (Robustness) second data set; sensitivity to the gap-merging rule and to the timetable encoding.

## 9. Risks

* The timetable may dominate the dynamics; then the forced-event encoding and the time-rescaling test decide whether
  the residual dynamics is Markovian, and the paper reports that.
* The aging model may fit better; see Sec. 5 (publishable either way, the quantification is the result).
* 20-second resolution and detection failures create spurious short gaps and short contacts; the gap-merging rule
  and a minimum duration need a sensitivity analysis.
* Identifiability of triadic closure versus same-class effects (classmates share many contacts); held-out
  likelihood and the class-pair stratification address this.

## 10. Status of the code (as of this plan)

Legs (degree dependence) and actions with decaying impacts (Hawkes-type acts) are implemented on branch
`claude/project-progress-plan-sjcbz2`, exactly in simulation, likelihood and gradient, and tested:
exact stationary distributions of small leg models against long simulations, the Hawkes log-likelihood against an
independent computation, the stationary act rate, replay and gradient checks, and debug-mode recomputation of all
rates after every event.

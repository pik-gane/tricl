# Handover

## State on 2026-11-01

* All work is on branch `claude/project-progress-plan-sjcbz2`, open as PR #10 against `master`
  (https://github.com/pik-gane/tricl/pull/10, "Fixes #1–#9"). Not merged; merging is the author's decision.
* CI: Release and Debug build-and-test passed on the paper commit; the sanitizer job failed on undefined
  behaviour in the hash functions (`src/data_model.h`), fixed in commit `9c291cf` (checked locally with the
  same sanitizer build and commands as the CI job). **First thing to do: check that CI is green on the
  latest head of PR #10**, and fix whatever is red.
* Test suite: `python3 tests/run_tests.py` (18 tests, all passing). Regression references in
  `tests/reference.json` are specific to libstdc++'s random number generator.
* Paper: `paper/` holds a draft for Physical Review E (REVTeX 4.2, 11 pages including appendices, `paper/main.pdf`).
  Open items for the author are listed below.
* Read `CLAUDE.md` for the paper's LaTeX source conventions (one sentence per source line, enforced in CI).

## What the branch contains

See the PR description and the change log at the end of `README.md`. In short:

1. Correctness: config parsing and validation, log-likelihood, rate bookkeeping, a bug that never generated
   random initial links of a symmetric relationship type between entities of different types (also on master).
2. Performance: dense type tables (memory 1 GB → about 10 MB), galloping angle intersection, Gillespie's direct
   method over a sum tree of rates (`src/schedule.h`).
3. The sigmoid: q-exponential form, slope 1 at zero for all tail indices, identical to the old function for
   the default tails [1, 1] (`src/probability.h`, unit tests in `tests/test_sigmoid.cpp`). The author specified it.
4. Estimation: replay mode `--events-in` with exact log-likelihood and analytic gradients (`--grad`),
   `python/tricl_fit.py` (maximum likelihood, standard errors), `python/tricl_recovery.py`.
5. Tools: `python/tricl_macro.py` (mean-field ODEs from a config), `python/tricl_rdf.py` (RDF ↔ tricl),
   `--links-out`/`--entities-out` and `python/tricl_movie.py` (movies without Gephi),
   `python/tricl_atop.py` (alliance data, first real-data fit).
6. Engineering: CMake build types, tests, GitHub Actions (Release, Debug, sanitizers, paper source
   conventions), `CONTRIBUTING.md`, `tests/benchmark.py`, `CITATION.cff`, `--max-events`, `--max-t`, `--max-wall`.

## Paper (`paper/`)

* Build: `python3 paper/make_data.py` (simulations, fits, benchmark → `paper/data/`, about 15 minutes),
  `python3 paper/make_figures.py` (→ `paper/generated/`: figures, tables, `numbers.tex`), `sh paper/build.sh`.
  Every number in the text is a macro from `generated/numbers.tex`; never edit generated files by hand.
* References: `paper/refs.bib`, all 26 entries checked against Crossref (via the WebFetch tool on
  `https://api.crossref.org/works/<DOI>`; direct curl to Crossref is blocked). One DOI (Demirel et al. 2014)
  was wrong in the first draft and is corrected.
* Target venue: Physical Review E. The APS policy on AI tools (updated June 2026) allows substantive AI use
  with disclosure of tool name and version, how it assisted, how it was directed and how the output was
  verified; the statement is in `paper/endmatter.tex`. The APS pages block automated access (Cloudflare);
  the policy was read from search results quoting the APS page and from a Research Information news item of
  18 June 2026. The author should read the policy page once in a browser before submitting.
* Appendices A–E (`paper/sec-appendix.tex`, added 2026-11-01) prove the paper's mathematical claims: sigmoid
  properties, exactness of the summary-event thinning, the log-likelihood with immediate events (and the
  direction-label issue of symmetric/inverse types), sufficient statistics/gradient/concavity/identifiability, and
  exactness of the mean-field equation without angle influences. Lean formalization was discussed and deferred.
  The author should check the proofs.
* **Open items for the author:** insert model names and versions at `[MODELS]` in `paper/endmatter.tex`
  (repository files must not contain model identifiers written by Claude); confirm the AI-use statement (in
  particular "reviewed the changes and the text"); acknowledgments and funding; author list and affiliation
  in `paper/main.tex`; cover letter at submission.
* The model section explains why self-links are excluded (as in the original 2020 code): a self-link would
  make every link of that entity count as a two-path through a third entity of the same type. If the author's
  original motive was different, adjust that sentence in `paper/sec-model.tex`.
* Venues discussed with the author: Physical Review E (recommended), Chaos, Journal of Complex Networks,
  Network Science, Applied Network Science, JASSS. JOSS is probably not eligible yet (its 2026 rules require
  development spread over more than six months and evidence of research use). Conferences: **NetSci 2027
  (Dresden, 17–21 May 2027), abstract deadline 1 December 2026** (one page with a figure; template on
  Overleaf, "netsci2027-abstract-template"); offered to draft it, not started. US conferences in March–August
  2027 were not announced yet (NERCCS and SIAM DS by past pattern).

## Key results (for quick reference)

* Alliance fit (ATOP v5.1 defense pledges, 1816–2018, 217 states, 4514 events): a common ally raises the
  formation rate about 58-fold; pact dissolution hazard 0.108/yr without common allies (2 d0, since both directed
  links of a symmetric pair are possible events; corrected from 0.054 on 2026-11-01); dropping the ally-of-ally
  term costs 9314 log-likelihood units. The profile likelihood in b1 rises again towards a plateau as
  b1 → −∞ (891 units below the maximum), so fits must start at positive b1 (the generated config does).
* Parameter recovery (three-entity model, 40 sequences): coverage of 95% intervals 0.90–1.00.
* Throughput: sparse models 6×10⁴–2.5×10⁵ events/s, dense 3-blocks about 2400 events/s.

## Setting up a fresh container

Container restarts wipe `build/` and installed packages.

```
apt-get install -y --no-install-recommends cmake g++ libyaml-cpp-dev libboost-dev zlib1g-dev graphviz \
    texlive-latex-base texlive-latex-recommended texlive-latex-extra texlive-fonts-recommended \
    texlive-science texlive-publishers lmodern latexmk poppler-utils
python3 -m pip install numpy matplotlib scipy pandas pyreadr networkx pillow   # use python3 -m pip
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j 4
```

`paper/make_data.py` downloads the ATOP and COW data from `raw.githubusercontent.com` (reachable); most
other data hosts (Wikidata, OpenAlex, SNAP, Zenodo, Overleaf) are blocked by the environment's network policy.

## Other context

* A private progress plan (not to be made public, never committed) was given to the author as a file in an
  earlier session; its order of work is done except the optional follow-ups below.
* The repositories `mensch72/responsibility_theory` and `mensch72/responsibility_prima_paper` were attached
  read-only to look up the paper conventions; nothing was changed there.
* Optional technical follow-ups: pair and triple closures for the mean-field generator, estimation from
  network snapshots (simulated moments or data augmentation), rules on single links ("legs", currently
  rejected by the config parser), a pybind11 module, a structural-balance example config, a v0.1 tag with a
  Zenodo DOI after merging.
* Commits end with the attribution lines given by the session; no model identifiers in repository files.

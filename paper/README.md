# Paper: tricl methods paper (draft)

Draft of a methods paper on tricl: the model class, the simulation algorithm, likelihood-based
estimation with gradients, the generated mean-field approximation, and a first application to
military alliances. Generic `article` class, so that it can be posted to arXiv and adapted to the
chosen journal's template later.

## Build

From the repository root, with tricl built in `build/`:

    python3 paper/make_data.py      # simulations, fits, benchmark -> paper/data/   (about 15 minutes)
    python3 paper/make_figures.py   # figures, tables, numbers     -> paper/generated/
    sh paper/build.sh               # -> paper/main.pdf (needs TeX Live with listings, tikz, natbib)

`make_data.py` needs numpy, matplotlib, scipy, pandas and pyreadr; it downloads the ATOP and COW
data from the peacesciencer repository into `paper/data/atop_raw/`. `--only NAME` reruns one part.

## Rules

* Every number in the text comes from `generated/numbers.tex`, every table from `generated/`;
  never edit generated files by hand, rerun the scripts.
* References: only entries checked against Crossref (`refs.bib`; checked 2026-11-01).
* Claims about the code must match the code (e.g. leg rules are rejected by the current version).
* The benchmark depends on the machine; rerun it on the machine whose numbers you want to report.

## Open items for the author

* Author list, affiliation, funding, competing interests (`endmatter.tex`), model names in the
  generative-AI statement ([MODELS]).
* Choice of venue (see the pull request discussion) and conversion to its template.

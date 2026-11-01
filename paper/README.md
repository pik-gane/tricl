# Paper: tricl methods paper (draft)

Draft of a methods paper on tricl for Physical Review E: the model class, the simulation algorithm,
likelihood-based estimation with gradients, the generated mean-field approximation, and a first
application to military alliances. Formatted with REVTeX 4.2 (`aps,pre,reprint`), starting from the
official APS template `apstemplate.tex` of the REVTeX distribution; APS style for cross-references
(Sec., Fig., Eq.) and American spelling. Replace `reprint` by `preprint` for a single-column manuscript.
arXiv accepts the same source.

## Build

From the repository root, with tricl built in `build/`:

    python3 paper/make_data.py      # simulations, fits, benchmark -> paper/data/   (about 15 minutes)
    python3 paper/make_figures.py   # figures, tables, numbers     -> paper/generated/
    sh paper/build.sh               # -> paper/main.pdf (needs TeX Live with REVTeX 4.2, e.g. Debian's texlive-publishers, and listings)

`make_data.py` needs numpy, matplotlib, scipy, pandas and pyreadr; it downloads the ATOP and COW
data from the peacesciencer repository into `paper/data/atop_raw/`. `--only NAME` reruns one part.

## Rules

* Every number in the text comes from `generated/numbers.tex`, every table from `generated/`;
  never edit generated files by hand, rerun the scripts.
* References: only entries checked against Crossref (`refs.bib`; checked 2026-11-01).
* Claims about the code must match the code (e.g. leg rules are rejected by the current version).
* The benchmark depends on the machine; rerun it on the machine whose numbers you want to report.

## Open items for the author

* Author list and affiliation (`main.tex`), acknowledgments and funding (`endmatter.tex`).
* The AI-use statement (`endmatter.tex`) follows the APS policy on AI tools as updated in June 2026,
  which permits substantive use with disclosure of the tool name and version, how it assisted, how
  the authors directed it and how they verified the output. Insert the model names and versions at
  [MODELS] and confirm the statement.
* APS asks for a cover letter, and PACS-free subject areas are chosen at submission.

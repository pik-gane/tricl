# Paper: memory and network context in face-to-face contact dynamics (draft)

Research paper for Physical Review E on the contacts of a high school (SocioPatterns 2013): Markovian models with
network context and decaying memory, fitted with tricl's exact likelihood, against a non-Markovian aging model;
held-out likelihood, time-rescaling residuals, sensitivity, instability of fitted network feedback and a pair closure.
Formatted with REVTeX 4.2 (`aps,pre,reprint`); the REVTeX class files and `main.bbl` are included, so that the paper
builds without TeX Live's `texlive-publishers`. The earlier tricl methods paper is in the git history (commit ebf5431).

## Build

From the repository root, with tricl built in `build/`:

    (run the study scripts, see study/README.md)          # -> study/results/
    python3 paper/make_figures.py                          # figures, tables, numbers -> paper/generated/
    cd paper && pdflatex main && bibtex main && pdflatex main && pdflatex main

The generated files are committed, so only the last line is needed to build the PDF.

## Rules

* LaTeX source conventions: one sentence per source line (also in captions), `\citet` instead of
  written author names, APS abbreviations, few paragraph breaks (details in `../CLAUDE.md`). Enforce with
  `python3 paper/tools/rewrap.py FILE...`; `--check` is run in CI.
* Every number in the text comes from `generated/numbers.tex`, every table from `generated/`;
  never edit generated files by hand, rerun the scripts.
* References: only entries checked against Crossref (`refs.bib`; checked 2026-11-01).
* Claims about the code must match the code.

## Open items for the author

* Author list and affiliation (`main.tex`), acknowledgments and funding (`endmatter.tex`).
* The AI-use statement (`endmatter.tex`) follows the APS policy on AI tools as updated in June 2026,
  which permits substantive use with disclosure of the tool name and version, how it assisted, how
  the authors directed it and how they verified the output. Insert the model names and versions at
  [MODELS] and confirm the statement.
* APS asks for a cover letter, and PACS-free subject areas are chosen at submission.

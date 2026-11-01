# Project notes for Claude

## Paper (`paper/`)

LaTeX source conventions (taken over from the responsibility_theory project; enforced in CI):
* **One sentence per source line**; split a sentence over several lines only at logical units when it is
  overly long. This applies to prose, the abstract and figure/table captions. Run
  `python3 paper/tools/rewrap.py paper/main.tex paper/sec-*.tex paper/endmatter.tex` after editing and check
  that the PDF text is unchanged (`pdftotext`); `--check` lists violating files (the CI runs it).
  Do not start a sentence with a numeral (the tool cannot split there, and it is bad style).
* `\citet` instead of explicitly written author names.
* Abbreviated cross-references in APS style: Sec., Fig., Eq., Ref., Table (full word at the start of a sentence).
* Not too many paragraph breaks.
* American spelling (APS house style), REVTeX 4.2 from the official APS template (`aps,pre,reprint`).
* Every number in the text comes from `paper/generated/numbers.tex`, every table and figure from
  `paper/generated/`, produced by `paper/make_figures.py` from `paper/data/` (`paper/make_data.py`);
  never edit generated files by hand.
* References only from entries checked against Crossref (`paper/refs.bib`).
* Claims about the code must match the code.
* The AI-use statement (`paper/endmatter.tex`) follows the APS policy on AI tools; model names and versions
  are inserted by the author.

## Code

Build with CMake (`build/`), test with `python3 tests/run_tests.py`; regression references are
libstdc++-specific and regenerated with `--update` only after a change that legitimately alters trajectories.

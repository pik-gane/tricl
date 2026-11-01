#!/bin/sh
# Build paper/main.pdf from the generated figures and tables (run paper/make_data.py and
# paper/make_figures.py first if paper/data/ or paper/generated/ is missing or outdated).
set -e
cd "$(dirname "$0")"
pdflatex -interaction=nonstopmode -halt-on-error main >/dev/null
bibtex main >/dev/null
pdflatex -interaction=nonstopmode -halt-on-error main >/dev/null
pdflatex -interaction=nonstopmode -halt-on-error main | grep -E "Output written|undefined" || true

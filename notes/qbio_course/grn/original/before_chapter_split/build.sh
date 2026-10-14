#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
python3 code/make_figures.py
python3 code/topology_analysis.py
python3 code/extended_examples.py
python3 code/positive_feedback_oscillation.py
xelatex -interaction=nonstopmode -halt-on-error main.tex >/dev/null
xelatex -interaction=nonstopmode -halt-on-error main.tex >/dev/null
xelatex -interaction=nonstopmode -halt-on-error main.tex >/dev/null
printf 'Built %s/main.pdf\n' "$PWD"

#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
python3 code/make_figures.py
python3 code/topology_analysis.py
python3 code/extended_examples.py
python3 code/positive_feedback_oscillation.py
python3 code/bistability_details.py
printf 'Updated PNG figures and topology inventory in %s\n' "$PWD"

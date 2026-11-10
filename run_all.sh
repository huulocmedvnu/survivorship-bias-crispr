#!/bin/sh
# Fast reproduction: unit tests, then every figure and figures/numbers.json from the tables in outputs/.
set -e
cd "$(dirname "$0")"
python3 -m pytest -q tests
python3 figures/make_figures.py
python3 figures/make_flowchart.py
echo "figures and numbers.json written to figures/"

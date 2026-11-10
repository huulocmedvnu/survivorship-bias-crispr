#!/bin/bash
# Trimming bounds in the K562 genome-wide screen (day 8) along the day 6 directions.
set -e
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
"$ROOT/rust_scorer/target/release/survb" "$ROOT/data/dosecurve_cache/gwps" "$ROOT/data/intermediate/phaseP_dirs.tsv" \
  "$ROOT/source_data/phaseP_pert.tsv" "$ROOT/source_data/phaseP_bounds.tsv"

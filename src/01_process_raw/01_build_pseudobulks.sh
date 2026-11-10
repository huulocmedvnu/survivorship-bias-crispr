#!/bin/bash
# Per-cell caches (dosecurve 1.0.0) and pseudobulk means (rust_scorer) for the five Replogle and Nadig screens.
# DOSECURVE must point to the dosecurve binary (https://doi.org/10.5281/zenodo.23936124).
set -e
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"; D="$ROOT/data"; S="$ROOT/rust_scorer/target/release/survb"
: "${DOSECURVE:?set DOSECURVE=path/to/dosecurve}"
mkdir -p "$D/pseudobulk"
for s in k562ess gwps rpe1 jurkat hepg2; do
  "$DOSECURVE" cache "$D/h5ad/$s.h5ad" "$D/dosecurve_cache/$s" --col perturbation --control control
  "$S" pb "$D/dosecurve_cache/$s" "$D/pseudobulk/$s"
done

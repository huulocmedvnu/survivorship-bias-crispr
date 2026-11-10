#!/bin/bash
# Download the public raw data into data/ (large: the K562 genome-wide h5ad alone is tens of GB). Only needed for the full reproduction
# (make reproduce-full). The fast reproduction uses source_data/ and needs no download.
# X-Atlas/Orion count files are streamed later by src/01_process_raw/11_orion_stream.py.
set -e
ROOT="$(cd "$(dirname "$0")/.." && pwd)"; D="$ROOT/data"
get() { mkdir -p "$(dirname "$2")"; [ -s "$2" ] || curl -L --fail --retry 5 -C - -o "$2" "$1"; }

# Replogle et al. 2022, via scPerturb (Zenodo 13350497)
Z=https://zenodo.org/records/13350497/files
get "$Z/ReplogleWeissman2022_K562_essential.h5ad?download=1" "$D/h5ad/k562ess.h5ad"
get "$Z/ReplogleWeissman2022_K562_gwps.h5ad?download=1"      "$D/h5ad/gwps.h5ad"
get "$Z/ReplogleWeissman2022_rpe1.h5ad?download=1"           "$D/h5ad/rpe1.h5ad"

# Nadig et al. 2025, GEO GSE264667 (add the perturbation column with dosecurve tools/h5ad_addcol.py before caching)
G=https://ftp.ncbi.nlm.nih.gov/geo/series/GSE264nnn/GSE264667/suppl
get "$G/GSE264667_jurkat_raw_singlecell_01.h5ad" "$D/h5ad/jurkat.h5ad"
get "$G/GSE264667_hepg2_raw_singlecell_01.h5ad"  "$D/h5ad/hepg2.h5ad"

# Drepanos et al. 2026 timecourse (Zenodo 20722064)
for d in 4 7 10 14; do get "https://zenodo.org/records/20722064/files/K562_Day$d.Rds?download=1" "$D/drepanos/K562_Day$d.Rds"; done

# Tahoe-100M metadata and X-Atlas/Orion gene metadata (Hugging Face)
T=https://huggingface.co/datasets/tahoebio/Tahoe-100M/resolve/main/metadata
for f in obs_metadata sample_metadata cell_line_metadata; do get "$T/$f.parquet" "$D/tahoe/$f.parquet"; done
get https://huggingface.co/datasets/Xaira-Therapeutics/X-Atlas-Orion/resolve/main/metadata/gene_metadata.parquet "$D/orion/gene_metadata.parquet"

cat <<'MSG'
Download these files by hand (the portals require a browser session):
  DepMap Public 24Q4: CRISPRGeneEffect.csv, Model.csv                         -> data/depmap/
  PRISM Repurposing 19Q4: secondary-screen-dose-response-curve-parameters.csv -> data/depmap/
      https://depmap.org/portal/data_page
  Replogle et al. 2022 Tables S1 and S2 (mmc1.xlsx, mmc2.xlsx)                -> data/replogle_supp/
      https://doi.org/10.1016/j.cell.2022.05.013
  Dose-response labels regime_v4.tsv (dosecurve 1.0.0, tools/realmap2.py)     -> data/dosecurve/
      https://doi.org/10.5281/zenodo.23936124
MSG

#!/bin/bash
# Rebuild every table in source_data/ from the raw public data under data/ (see README, section 4).
set -e
H="$(cd "$(dirname "$0")" && pwd)"
bash "$H/01_build_pseudobulks.sh"
python3 "$H/02_k562_day6_day8_ratio.py"
python3 "$H/03_k562_gwps_directions.py"
bash "$H/04_k562_gwps_bounds.sh"
python3 "$H/05_k562_gwps_bounds_summary.py"
python3 "$H/06_single_timepoint_bounds.py"
python3 "$H/07_buffered_labels.py"
Rscript "$H/08_drepanos_export.R"
python3 "$H/09_drepanos_timecourse.py"
python3 "$H/10_orion_cell_counts.py"
bash "$H/11_run_orion_stream.sh"
python3 "$H/12_orion_bounds.py" HCT116
python3 "$H/12_orion_bounds.py" HEK293T
python3 "$H/13_tahoe_prism.py"
python3 "$H/14_reference_tables.py"

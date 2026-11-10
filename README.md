# Survivorship bias in single-cell CRISPR perturbation studies

Code and source data for the article "Survivorship bias in single-cell CRISPR perturbation
studies" by Thai Huu Loc Tran and Thi Van Nguyen.

Single-cell CRISPR screens profile only the cells that survive until capture, and for essential
genes the cells with the most complete knockdown divide less or die first. This repository
measures that bias across timepoints and bounds its contribution to the reported responses with
a trimming (Lee) bound. Every number and figure in the article can be rebuilt with one command.

## Quick start for reviewers (under 1 minute)

```
pip install -r environment/requirements.txt
make
```

`make` (the same as `make reproduce-fast`) runs the unit tests, recomputes all 69 numbers quoted in
the Abstract and Results from `source_data/` and compares each with the article, then redraws
every figure. It finishes in about 20 seconds and needs no raw data.

- `output/tables/manuscript_claims_verification.txt`: one line per number, with the value in
  the article, the recomputed value and MATCH or MISMATCH
- `output/tables/numbers.json`: every estimate, confidence interval and P value
- `output/figures/`: Figures 1 to 5 and Figs. S1 and S2 (PNG and PDF)

Without `make`, run `python3 src/02_statistical_tests/verify_manuscript_numbers.py` and the six
scripts in `src/03_figures/`. Scripts resolve paths from the repository root, so they run from
any directory.

## Repository layout

```
README.md  LICENSE  CITATION.cff  Makefile
environment/      requirements.txt (Python 3.12, pinned), r_packages.R (R 4.4)
rust_scorer/      Rust tool survb: pseudobulk means and control-cell scores (sparse x dense product)
src/
  paths.py        repository paths
  analysis.py     every reported estimate, computed from source_data/
  00_download_data.sh
  01_process_raw/ raw public data -> source_data/ (run_full_pipeline.sh runs all steps in order)
  02_statistical_tests/verify_manuscript_numbers.py
  03_figures/     one script per figure
source_data/      per-perturbation tables behind every figure and number (10 MB)
output/           written by make: figures/ and tables/
tests/            unit tests of the bound and the manuscript check
```

## Manuscript to code

| Figure or result | Source data (`source_data/`) | Script | Output (`output/`) |
|---|---|---|---|
| Fig. 1, study design | `phaseQ_drepano.tsv` (readout marks) | `src/03_figures/plot_fig1_schematic.py` | `figures/Figure1.*` |
| Fig. 2a, cells per perturbation | `phaseP_killtest.tsv` | `src/03_figures/plot_fig2_k562_timepoints.py` | `figures/Figure2.*` |
| Fig. 2b, knockdown reported by the authors | `ref_replogle_knockdown.tsv`, `ref_k562_chronos.tsv` | same | same |
| Fig. 2c, day 8 / day 6 response and placebo | `phaseP_killtest.tsv` | same | same |
| Identical guides in both libraries | `phaseP_same_guides.tsv` | `src/02_statistical_tests/verify_manuscript_numbers.py` | `tables/numbers.json` |
| Fig. 3a to c, Drepanos timecourse | `phaseQ_drepano.tsv` | `src/03_figures/plot_fig3_timecourse.py` | `figures/Figure3.*` |
| Fig. 4, explainable fraction in seven screens | `phaseP_bounds_summary.tsv`, `phaseR_{k562ess,rpe1,jurkat,hepg2}_{bounds,pert}.tsv`, `orion_o2_{HCT116,HEK293T}.tsv` | `src/03_figures/plot_fig4_trimming_bounds.py` | `figures/Figure4.*` |
| Fig. 5a and b, buffered genes | `phaseP_dose_labels.tsv` | `src/03_figures/plot_fig5_buffered_genes.py` | `figures/Figure5.*` |
| Fig. S1a and b, Tahoe-100M against PRISM | `tahoe_t1_prism.tsv.gz`, `tahoe_t1_survival.tsv.gz` | `src/03_figures/plot_supp_figS1_S2.py` | `figures/FigureS1.*` |
| Fig. S2, cells per perturbation by essentiality | `ref_gwps_cells.tsv`, `ref_k562_chronos.tsv`, `orion_o1_{HCT116,HEK293T}.tsv` | same | `figures/FigureS2.*` |
| All numbers in the Abstract and Results | all of the above | `src/02_statistical_tests/verify_manuscript_numbers.py` | `tables/manuscript_claims_verification.txt` |

Each table in `source_data/` is written by one script in `src/01_process_raw/`, named in the
table of section "Full reproduction" below.

## Determinism

Bootstrap confidence intervals use 1,000 resamples drawn from one generator,
`numpy.random.default_rng(1)`, in a fixed order (`src/analysis.py`), so every run returns the same
values to the last digit. The two halves of each perturbation are the cells at even and odd
positions in the order of the original count matrix. X-Atlas/Orion cells are assigned to halves by
the lowest bit of the CRC-32 checksum of the cell barcode. The cell pictograms of Fig. 1 use
`random.seed(3)`.

## System requirements and installation

- Platform independent (macOS and Linux)
- Python 3.12 with the packages in `environment/requirements.txt`
- Rust stable toolchain (https://rustup.rs), for the full reproduction only
- R 4.4 with Seurat 5, for one step of the full reproduction only (`environment/r_packages.R`)
- `rsvg-convert` (librsvg) for Fig. 1 (`brew install librsvg` or `apt install librsvg2-bin`)

`make setup` installs the Python packages and builds the Rust scorer, which takes about 2 minutes.

## Full reproduction from the public raw data

`make reproduce-full` rebuilds every table in `source_data/` and then runs `make reproduce-fast`.
It reads per-cell data and is best run on a Linux server. Download the data first with
`make download` (`src/00_download_data.sh`), which also lists the files that must be fetched
from a browser.

### Data

| Data | Identifier | Local path |
|---|---|---|
| Replogle K562 essential (day 6), K562 genome-wide (day 8), RPE1 (day 7) | scPerturb, https://doi.org/10.5281/zenodo.13350497 | `data/h5ad/` |
| Replogle Tables S1 and S2 | https://doi.org/10.1016/j.cell.2022.05.013 | `data/replogle_supp/` |
| Nadig Jurkat and HepG2 | GEO GSE264667 | `data/h5ad/` |
| X-Atlas/Orion HCT116 and HEK293T | https://huggingface.co/datasets/Xaira-Therapeutics/X-Atlas-Orion | `data/orion/` (streamed) |
| Drepanos K562 timecourse | https://doi.org/10.5281/zenodo.20722064 | `data/drepanos/` |
| Tahoe-100M metadata | https://huggingface.co/datasets/tahoebio/Tahoe-100M | `data/tahoe/` |
| DepMap Public 24Q4 Chronos, PRISM Repurposing 19Q4 secondary screen | https://depmap.org/portal/data_page | `data/depmap/` |
| Dose-response labels and per-cell caches (dosecurve 1.0.0) | https://doi.org/10.5281/zenodo.23936124 | `data/dosecurve/`, `data/dosecurve_cache/` |

### Steps

Set `DOSECURVE` to the dosecurve binary, then run `src/01_process_raw/run_full_pipeline.sh`,
which runs:

| Step | Script | Writes |
|---|---|---|
| 01 | `01_build_pseudobulks.sh` (dosecurve cache, `survb pb`) | `data/dosecurve_cache/`, `data/pseudobulk/` |
| 02 | `02_k562_day6_day8_ratio.py` | `phaseP_killtest.tsv` |
| 03 | `03_k562_gwps_directions.py` | `phaseP_pert.tsv`, `data/intermediate/phaseP_dirs.tsv` |
| 04 | `04_k562_gwps_bounds.sh` (`survb`) | `phaseP_bounds.tsv` |
| 05 | `05_k562_gwps_bounds_summary.py` | `phaseP_bounds_summary.tsv` |
| 06 | `06_single_timepoint_bounds.py` (`survb`) | `phaseR_*.tsv` |
| 07 | `07_buffered_labels.py` | `phaseP_dose_labels.tsv` |
| 08 | `08_drepanos_export.R` | `data/drepanos/pb/` |
| 09 | `09_drepanos_timecourse.py` | `phaseQ_drepano.tsv` |
| 10 | `10_orion_cell_counts.py` | `orion_o1_*.tsv` |
| 11 | `11_run_orion_stream.sh` (resumable) | `data/orion/state/` |
| 12 | `12_orion_bounds.py HCT116` and `HEK293T` (`survb`) | `orion_o2_*.tsv` |
| 13 | `13_tahoe_prism.py` | `tahoe_t1_*.tsv.gz` |
| 14 | `14_reference_tables.py` | `ref_*.tsv`, `phaseP_same_guides.tsv` |

`tests/check_survb_reference.py` compares the Rust scorer with NumPy on the K562 genome-wide
cache after step 04.

## Citation and license

See `CITATION.cff`. MIT license, see `LICENSE`.

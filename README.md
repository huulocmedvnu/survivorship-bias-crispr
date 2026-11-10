# Survivorship bias in single-cell CRISPR perturbation studies

Code and source data for the manuscript "Survivorship bias in single-cell CRISPR
perturbation studies" by Thai Huu Loc Tran and Thi Van Nguyen.

Single-cell CRISPR interference screens profile the cells that survive until capture.
For essential genes, cells with the most complete knockdown divide less or die first.
This repository measures that bias across timepoints and bounds its contribution to the
reported transcriptional responses with a trimming (Lee) bound computed on control cells.

- Platform: platform independent (macOS and Linux)
- Languages: Python 3.12, R 4.4, Rust (stable, 2021 edition)
- License: MIT (see `LICENSE`)

Run every command from the repository root.

## Layout

| Folder | Content |
|---|---|
| `scripts/` | Python and R analyses of the Replogle, Nadig and Drepanos screens |
| `tahoe/` | Streaming analyses of X-Atlas/Orion and Tahoe-100M (HTTP range reads, resumable) |
| `tools/survb/` | Rust tool: pseudobulk means (`survb pb`) and the trimming bound (`survb`) |
| `figures/` | `make_figures.py` and `make_flowchart.py`, which read `outputs/` and write the figures and `numbers.json` |
| `outputs/` | Per-perturbation source data from which all figures and reported estimates are computed |
| `tests/` | Unit tests of the bound |

## Requirements

```
pip install -r requirements.txt          # Python 3.12, exact versions used for the manuscript
cargo build --release --locked --manifest-path tools/survb/Cargo.toml
```

R 4.4 with Seurat 5 and Matrix is needed only for `scripts/phaseQ_drepano_export.R`
(`install.packages(c("Seurat", "Matrix"))`). `make_flowchart.py` calls `rsvg-convert`
(librsvg) to write Figure 1 as PNG and PDF.

## 1. Fast reproduction from the source data (about 10 seconds)

```
sh run_all.sh
```

This runs the unit tests and rebuilds Figures 1 to 5, Figs. S1 and S2 and
`figures/numbers.json`, which holds every estimate, confidence interval, regression
coefficient and P value quoted in the manuscript, from the tables in `outputs/`.
Bootstrap confidence intervals use 1,000 resamples drawn from `numpy.random.default_rng(1)`.

## 2. Full reproduction from the public raw data

All input data are public. Place them under `data/` (ignored by git) as below, then run
the steps in order. Steps 2.2 to 2.6 read per-cell data and are best run on a Linux server.

### 2.1 Download

| Data | Source | Local path |
|---|---|---|
| Replogle K562 essential (day 6), K562 genome-wide (day 8), RPE1 (day 7) | scPerturb, https://doi.org/10.5281/zenodo.13350497 | `data/h5ad/` |
| Replogle Tables S1 and S2 (`mmc1.xlsx`, `mmc2.xlsx`) | supplementary data of https://doi.org/10.1016/j.cell.2022.05.013 | `data/replogle_supp/` |
| Nadig Jurkat and HepG2 essential screens | GEO GSE264667 | `data/h5ad/` |
| Drepanos K562 timecourse (`K562_Day{4,7,10,14}.Rds`) | https://doi.org/10.5281/zenodo.20722064 | `ext/drepano/` |
| X-Atlas/Orion HCT116 and HEK293T | https://huggingface.co/datasets/Xaira-Therapeutics/X-Atlas-Orion | streamed by `tahoe/` |
| Tahoe-100M metadata (`obs_metadata.parquet`, `sample_metadata.parquet`, `cell_line_metadata.parquet`) | https://huggingface.co/datasets/tahoebio/Tahoe-100M | `tahoe/meta/` |
| DepMap Public 24Q4 `CRISPRGeneEffect.csv`, `Model.csv` | https://depmap.org/portal/data_page | `data/depmap/` |
| PRISM Repurposing 19Q4 `secondary-screen-dose-response-curve-parameters.csv` | https://depmap.org/portal/data_page | `data/depmap/` |

### 2.2 Per-cell caches and pseudobulks

Caches are built with dosecurve 1.0.0 (https://doi.org/10.5281/zenodo.23936124), which
also provides `h5ad_addcol.py` and `genes_from_var.py` for the GEO files (see the dosecurve
README). Name each h5ad file after its screen:

```
B=path/to/dosecurve/tools/dosecurve/target/release/dosecurve
for s in k562ess gwps rpe1 jurkat hepg2; do
  $B cache data/h5ad/$s.h5ad data/dosecurve_cache/$s --col perturbation --control control
  tools/survb/target/release/survb pb data/dosecurve_cache/$s data/pseudobulk/$s
done
```

Cells of each perturbation keep their order in the original count matrix. The two halves
are the cells at even and odd positions of that order.

### 2.3 Response ratio and trimming bounds in the Replogle and Nadig screens

```
python3 scripts/phaseP_killtest.py                 # day 8 / day 6 response ratio -> outputs/phaseP_killtest.tsv
python3 scripts/phaseP_directions.py               # directions from day 6 -> outputs/phaseP_dirs.tsv, phaseP_pert.tsv
tools/survb/target/release/survb data/dosecurve_cache/gwps outputs/phaseP_dirs.tsv outputs/phaseP_pert.tsv outputs/phaseP_bounds.tsv
python3 scripts/phaseP_bounds_summary.py           # -> outputs/phaseP_bounds_summary.tsv
python3 scripts/phaseR_lines.py                    # K562 day 6, RPE1, Jurkat, HepG2 (calls survb) -> outputs/phaseR_*.tsv
```

`outputs/phaseP_dirs.tsv` (13 MB) is not included and is rebuilt by
`scripts/phaseP_directions.py`.

### 2.4 Dose-response labels

`data/dosecurve/regime_v4.tsv` holds the labels of the three dose-response methods for the
K562 essential screen, written by `tools/realmap2.py` of dosecurve 1.0.0.

```
python3 scripts/phaseP_dose_labels.py              # -> outputs/phaseP_dose_labels.tsv
```

### 2.5 Timecourse

```
Rscript scripts/phaseQ_drepano_export.R            # -> ext/drepano/pb/
python3 scripts/phaseQ_drepano_analysis.py         # -> outputs/phaseQ_drepano.tsv
```

### 2.6 X-Atlas/Orion and Tahoe-100M

```
python3 tahoe/o1_orion_counts.py                   # cells per target -> outputs/orion_o1_*.tsv
sh tahoe/run_o2.sh                                 # streams both lines, resumable
python3 tahoe/o2_finalize.py HCT116                # -> outputs/orion_o2_HCT116.tsv
python3 tahoe/o2_finalize.py HEK293T               # -> outputs/orion_o2_HEK293T.tsv
python3 tahoe/t1_survival.py                       # -> outputs/tahoe_t1_*.tsv.gz
```

X-Atlas/Orion cells are split into halves by the lowest bit of the CRC-32 checksum of the
cell barcode (`zlib.crc32(barcode) & 1`).

### 2.7 Reference tables and figures

```
python3 scripts/export_reference_tables.py        # outputs/ref_*.tsv, outputs/phaseP_same_guides.tsv
sh run_all.sh
```

## License

MIT. See `LICENSE`.

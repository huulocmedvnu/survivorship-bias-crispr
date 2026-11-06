# Survivorship bias in single-cell CRISPR perturbation studies

Code and source data for the manuscript "Survivorship bias in single-cell CRISPR
perturbation studies" by Thai Huu Loc Tran and Thi Van Nguyen.

Single-cell CRISPR interference screens profile the cells that survive until capture.
For essential genes, cells with the most complete knockdown divide less or die first.
This repository measures that bias across timepoints and bounds its contribution to the
reported transcriptional responses with a trimming (Lee) bound computed on control cells.

## Layout

| Folder | Content |
|---|---|
| `scripts/` | Python and R analyses of the Replogle, Nadig, Drepanos and Dixit screens |
| `tahoe/` | Streaming analyses of X-Atlas/Orion and Tahoe-100M (HTTP range reads, resumable) |
| `tools/survb/` | Rust implementation of the trimming bound (`cargo build --release`) |
| `figures/` | `make_figures.py` and `make_flowchart.py`, which read `outputs/` and write the figures and `numbers.json` |
| `outputs/` | Per-perturbation tables from which all figures and reported estimates were computed |
| `tests/` | Unit test of the bound and a check of `survb` against NumPy |

`outputs/phaseP_dirs.tsv` (200-gene direction weights, 13 MB) is not included and is rebuilt
by `scripts/phaseP_directions.py`.

## Reproducing the figures

The figures and every number quoted in the text are rebuilt from the committed tables:

```
pip install -r requirements.txt
python3 figures/make_figures.py
python3 figures/make_flowchart.py
```

## Input data

All input data are public. Scripts expect local copies under `data/`:

- `data/depmap/`: DepMap `CRISPRGeneEffect.csv` and `Model.csv` (https://depmap.org), and PRISM
  secondary screen dose-response curves.
- `data/pseudobulk/`: per-perturbation mean counts per 10,000 for the Replogle K562 essential
  (day 6), K562 genome-wide (day 8) and RPE1 screens (scPerturb), and the Nadig Jurkat and
  HepG2 screens (GEO GSE264667). Each screen is stored as `<screen>.f32` (float32 matrix of
  all cells and two halves of each perturbation), `<screen>.genes.txt` and `<screen>.groups.tsv`.
- `data/dosecurve_cache/`: sparse per-cell count cache written by `dosecurve cache`
  (https://doi.org/10.5281/zenodo.23936124).
- `data/dosecurve/`: dose-response labels from dosecurve.
- X-Atlas/Orion and Tahoe-100M are streamed directly from Hugging Face by the scripts in `tahoe/`.
- The Drepanos timecourse: https://doi.org/10.5281/zenodo.20722064.

## License

MIT. See `LICENSE`.

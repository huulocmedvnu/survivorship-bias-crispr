# Survivorship bias in single-cell CRISPR perturbation studies
# make                 same as make reproduce-fast
# make setup           install Python packages and build the Rust scorer
# make reproduce-fast  verify every manuscript number and redraw every figure from source_data/
# make reproduce-full  rebuild source_data/ from the raw public data, then reproduce-fast
# make package-zip     write Additional_file_2.zip from the committed files

PYTHON ?= python3
FIG := src/03_figures

.PHONY: reproduce-fast setup build-rust test verify figures reproduce-full download package-zip clean

reproduce-fast: test verify figures
	@echo "Done. Figures in output/figures/, verification in output/tables/manuscript_claims_verification.txt"

setup: build-rust
	$(PYTHON) -m pip install -r environment/requirements.txt

build-rust:
	cargo build --release --locked --manifest-path rust_scorer/Cargo.toml

test:
	$(PYTHON) -m pytest -q tests

verify:
	$(PYTHON) src/02_statistical_tests/verify_manuscript_numbers.py

figures:
	$(PYTHON) $(FIG)/plot_fig1_schematic.py
	$(PYTHON) $(FIG)/plot_fig2_k562_timepoints.py
	$(PYTHON) $(FIG)/plot_fig3_timecourse.py
	$(PYTHON) $(FIG)/plot_fig4_trimming_bounds.py
	$(PYTHON) $(FIG)/plot_fig5_buffered_genes.py
	$(PYTHON) $(FIG)/plot_supp_figS1_S2.py

download:
	bash src/00_download_data.sh

reproduce-full: build-rust
	Rscript environment/r_packages.R
	bash src/01_process_raw/run_full_pipeline.sh
	$(MAKE) reproduce-fast

package-zip:
	@git diff --quiet HEAD || { echo "commit your changes first: the zip is built from HEAD"; exit 1; }
	git archive --format=zip --prefix=survivorship-bias-crispr/ -o Additional_file_2.zip HEAD
	@echo "wrote Additional_file_2.zip"

clean:
	rm -rf output Additional_file_2.zip

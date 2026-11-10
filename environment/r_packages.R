# R 4.4 packages for src/01_process_raw/08_drepanos_export.R (reads the Drepanos Seurat objects).
# Run: Rscript environment/r_packages.R
repos <- "https://cloud.r-project.org"
for (p in c("Matrix", "SeuratObject", "Seurat")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = repos)
}
stopifnot(packageVersion("SeuratObject") >= "5.0.0")   # Layers() and LayerData() need SeuratObject 5
cat("R packages ready:", as.character(packageVersion("Seurat")), "\n")

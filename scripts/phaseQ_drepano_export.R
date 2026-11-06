# Phase Q input: per-guide pseudobulk (mean CP10k, all cells and two interleaved halves) for each K562 day of
# Drepanos et al. bioRxiv 2026.07.04.736492 (Zenodo 10.5281/zenodo.20722064), plus cells per guide.
# Run: micromamba run -n psR Rscript scripts/phaseQ_drepano_export.R [days]
# Memory: each sequencing run (layer) is processed on its own, and CP10k scaling is folded into the group
# indicator (entries 1e4 / library size / group size), so the count matrix is never copied.
# Day 14 has one run per modality; cells of the other modality found in a run (stray assignments) are dropped.
suppressMessages({ library(Seurat); library(Matrix) })
ext <- "ext/drepano"; out <- "ext/drepano/pb"; dir.create(out, showWarnings = FALSE)
days <- as.integer(commandArgs(TRUE)); if (!length(days)) days <- c(4, 7, 10, 14)

pb_run <- function(m, md) {
  run <- unique(md$orig.ident); stopifnot(length(run) == 1)
  ok <- !is.na(md$Guide) & md$Guide != "" & md$percent.mt <= 15
  if (grepl("_CRISPR(i|ko)$", run)) ok <- ok & md$modality %in% c("Both", sub(".*_", "", run))
  m <- m[, ok, drop = FALSE]; md <- md[ok, ]
  w <- 1e4 / Matrix::colSums(m)
  grp <- ifelse(grepl("^INTERGENIC", md$Guide), paste0("CONTROL_", run), paste0(md$modality, "|", md$Guide))
  half <- ave(seq_along(grp), grp, FUN = function(i) seq_along(i) %% 2)
  res <- list(n = table(grp))
  for (part in c("all", "h0", "h1")) {
    sel <- if (part == "all") rep(TRUE, length(grp)) else half == as.integer(part == "h1")
    g <- factor(grp[sel]); gi <- as.integer(g); n <- tabulate(gi)
    A <- sparseMatrix(i = which(sel), j = gi, x = w[sel] / n[gi], dims = c(length(grp), nlevels(g)),
                      dimnames = list(NULL, levels(g)))
    res[[part]] <- as.matrix(m %*% A)                               # genes x groups, mean CP10k
  }
  res
}

for (day in days) {
  x <- readRDS(file.path(ext, sprintf("K562_Day%d.Rds", day)))
  md_all <- x@meta.data
  genes <- rownames(x)
  layers <- Layers(x)
  mats <- lapply(layers, function(ly) LayerData(x, assay = "RNA", layer = ly))
  rm(x); gc(verbose = FALSE)
  parts <- list()
  for (k in seq_along(mats)) {
    m <- mats[[k]]; mats[[k]] <- NA
    md <- md_all[colnames(m), ]
    for (run in unique(md$orig.ident)) {
      s <- md$orig.ident == run
      parts[[run]] <- pb_run(m[, s, drop = FALSE], md[s, ])
    }
    rm(m); gc(verbose = FALSE)
  }
  ug <- Reduce(union, lapply(parts, function(p) rownames(p$all)))   # layers may hold different gene sets
  pad <- function(r) { z <- matrix(0, length(ug), ncol(r), dimnames = list(ug, colnames(r))); z[rownames(r), ] <- r; z }
  comb <- function(k) do.call(cbind, lapply(parts, function(p) pad(p[[k]])))
  all <- comb("all"); h0 <- comb("h0"); h1 <- comb("h1")
  stopifnot(!anyDuplicated(colnames(all)))
  h0 <- h0[, intersect(colnames(all), colnames(h0))]; h1 <- h1[, intersect(colnames(all), colnames(h1))]
  n <- do.call(rbind, lapply(parts, function(p) data.frame(group = names(p$n), n = as.integer(p$n))))
  keep <- rowMeans(all) > 0
  write.csv(round(all[keep, ], 5), gzfile(file.path(out, sprintf("day%d_all.csv.gz", day))))
  write.csv(round(h0[keep, ], 5), gzfile(file.path(out, sprintf("day%d_h0.csv.gz", day))))
  write.csv(round(h1[keep, ], 5), gzfile(file.path(out, sprintf("day%d_h1.csv.gz", day))))
  write.csv(n, file.path(out, sprintf("day%d_n.csv", day)), row.names = FALSE)
  cat(sprintf("day %d: %d groups, %d genes, runs: %s\n", day, ncol(all), sum(keep), paste(names(parts), collapse = ", ")))
  rm(parts, all, h0, h1); gc(verbose = FALSE)
}

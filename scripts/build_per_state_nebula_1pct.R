#!/usr/bin/env Rscript
# Companion to scripts/build_per_state_nebula.R, but with min_cell_frac=0.01
# so whole-microglia top hits that are too sparse for the canonical 5%
# threshold (Plac9, Cd276, AI504432, ...) get per-state stats.
#
# The 5% cache remains the published per-state reference because it matches
# what NEBULA actually has decent power on per state; the 1% cache exists
# solely to confirm that the missing whole-microglia top hits are genuinely
# absent per state rather than hidden behind the prevalence filter.
#
# Usage:
#   Rscript scripts/build_per_state_nebula_1pct.R              # all four states
#   Rscript scripts/build_per_state_nebula_1pct.R DAM IFN      # subset
#
# Writes: storage/cache/de_snrnaseq_nebula_per_state_1pct.rds

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})
source("R/helpers.R")

args <- commandArgs(trailingOnly = TRUE)
selected <- if (length(args) > 0) args else NULL

t_start <- Sys.time()
cat("[per-state nebula 1pct] loading microglia_seurat_processed.rds ...\n")
sc <- readRDS("storage/cache/microglia_seurat_processed.rds")
symbol_map <- readRDS("storage/cache/snrnaseq_symbol_map.rds")
cat(sprintf("[per-state nebula 1pct] loaded %d cells x %d genes in %.1f s.\n",
            ncol(sc), nrow(sc),
            as.numeric(difftime(Sys.time(), t_start, units = "secs"))))

cat("[per-state nebula 1pct] adding state labels via label_microglia_states() ...\n")
sc <- label_microglia_states(sc, symbol_map = symbol_map)

states_to_fit <- if (is.null(selected)) levels(sc$state) else selected
cat(sprintf("[per-state nebula 1pct] fitting states: %s\n",
            paste(states_to_fit, collapse = ", ")))

ncore <- max(1L, parallel::detectCores() - 2L)
cat(sprintf("[per-state nebula 1pct] using ncore=%d.\n", ncore))

cache_path <- "storage/cache/de_snrnaseq_nebula_per_state_1pct.rds"
existing <- if (file.exists(cache_path)) readRDS(cache_path) else list()

per_state <- fit_nebula_per_state(
  sc,
  state_col     = "state",
  states        = states_to_fit,
  id_col        = "genotype_batch",
  genotype_col  = "genotype",
  assay         = "RNA",
  layer         = "counts",
  min_cell_frac = 0.01,
  ncore         = ncore,
  symbol_map    = symbol_map
)

for (nm in names(per_state)) existing[[nm]] <- per_state[[nm]]
saveRDS(existing, cache_path)
Sys.chmod(cache_path, mode = "0664")
try(system2("chown", c("rstudio:rstudio", cache_path)), silent = TRUE)

cat("\n==== per-state NEBULA 1% summary (interaction contrast, FDR<0.05) ====\n")
for (st in names(existing)) {
  tbl <- existing[[st]]$top$interaction
  n_sig <- sum(tbl$adj.P.Val < 0.05, na.rm = TRUE)
  min_padj <- suppressWarnings(min(tbl$adj.P.Val, na.rm = TRUE))
  cat(sprintf("  %-15s n_genes=%5d n_sig=%4d min_padj=%.3g\n",
              st, nrow(tbl), n_sig, min_padj))
}

cat(sprintf("\n[per-state nebula 1pct] total wall: %.1f min.\n",
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
cat(sprintf("Cache written: %s (%.1f MB)\n",
            cache_path,
            file.info(cache_path)$size / 1024^2))

#!/usr/bin/env Rscript
# Pre-build the per-substate NEBULA cache outside of the Rmd so the knit
# is fast. Uses the canonical state-labelling helper from R/helpers.R so
# the cache stays in lock-step with the Rmd's substate annotation.
#
# Usage:
#   Rscript scripts/build_per_state_nebula.R              # all four states
#   Rscript scripts/build_per_state_nebula.R DAM IFN      # subset
#
# Writes: storage/cache/de_snrnaseq_nebula_per_state.rds

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})
source("R/helpers.R")

args <- commandArgs(trailingOnly = TRUE)
selected <- if (length(args) > 0) args else NULL

t_start <- Sys.time()
cat("[per-state nebula] loading microglia_seurat_processed.rds ...\n")
sc <- readRDS("storage/cache/microglia_seurat_processed.rds")
symbol_map <- readRDS("storage/cache/snrnaseq_symbol_map.rds")
cat(sprintf("[per-state nebula] loaded %d cells x %d genes in %.1f s.\n",
            ncol(sc), nrow(sc),
            as.numeric(difftime(Sys.time(), t_start, units = "secs"))))

cat("[per-state nebula] adding state labels via label_microglia_states() ...\n")
sc <- label_microglia_states(sc, symbol_map = symbol_map)
cat("\n==== state x genotype ====\n")
print(table(sc$state, sc$genotype))
cat("\n==== state x genotype_batch (min cells per id) ====\n")
tab <- table(sc$state, sc$genotype_batch)
for (s in rownames(tab)) {
  cat(sprintf("  %-15s n=%5d  median=%4d  min=%3d  zero_ids=%d\n",
              s, sum(tab[s, ]), as.integer(median(tab[s, ])),
              min(tab[s, ]), sum(tab[s, ] == 0)))
}

states_to_fit <- if (is.null(selected)) levels(sc$state) else selected
cat(sprintf("\n[per-state nebula] fitting states: %s\n",
            paste(states_to_fit, collapse = ", ")))

ncore <- max(1L, parallel::detectCores() - 2L)
cat(sprintf("[per-state nebula] using ncore=%d.\n", ncore))

cache_path <- "storage/cache/de_snrnaseq_nebula_per_state.rds"
existing <- if (file.exists(cache_path)) readRDS(cache_path) else list()

per_state <- fit_nebula_per_state(
  sc,
  state_col     = "state",
  states        = states_to_fit,
  id_col        = "genotype_batch",
  genotype_col  = "genotype",
  assay         = "RNA",
  layer         = "counts",
  min_cell_frac = 0.05,
  ncore         = ncore,
  symbol_map    = symbol_map
)

# Merge with anything already cached so a partial re-run (e.g. just IFN)
# doesn't wipe the other states.
for (nm in names(per_state)) existing[[nm]] <- per_state[[nm]]
saveRDS(existing, cache_path)
Sys.chmod(cache_path, mode = "0664")
try(system2("chown", c("rstudio:rstudio", cache_path)), silent = TRUE)

cat("\n==== per-state NEBULA summary (interaction contrast, FDR<0.05) ====\n")
for (st in names(existing)) {
  tbl <- existing[[st]]$top$interaction
  n_sig <- sum(tbl$adj.P.Val < 0.05, na.rm = TRUE)
  min_padj <- suppressWarnings(min(tbl$adj.P.Val, na.rm = TRUE))
  cat(sprintf("  %-15s n_genes=%5d n_sig=%4d min_padj=%.3g\n",
              st, nrow(tbl), n_sig, min_padj))
}

cat(sprintf("\n[per-state nebula] total wall: %.1f min.\n",
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
cat(sprintf("Cache written: %s (%.1f MB)\n",
            cache_path,
            file.info(cache_path)$size / 1024^2))

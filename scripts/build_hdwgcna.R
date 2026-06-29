#!/usr/bin/env Rscript
# Pre-build the hdWGCNA cache outside of analysis.Rmd so the knit stays
# fast. The pipeline is monolithic (metacells -> TOM -> modules -> MEs
# -> connectivity) and lives entirely inside the Seurat object's `@misc$
# wgcna` slot.
#
# Outputs:
#   storage/cache/hdwgcna_microglia.rds  - lightweight summary (modules,
#       hub_genes, MEs_cells, module_meta, power_table, soft_power,
#       wgcna_obj). The full Seurat is *not* saved here; metadata needed
#       for downstream DE is bundled separately.
#   storage/cache/hdwgcna_module_de.rds  - per-module DE under the project
#       2x2 factorial design with batch covariate (5 contrasts).
#
# Usage:
#   Rscript scripts/build_hdwgcna.R                 # full pipeline
#   Rscript scripts/build_hdwgcna.R --soft-power 8  # skip soft-power scan
#
# Wall time: ~15-30 minutes on 8 cores.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(hdWGCNA)
  library(WGCNA)
})
source("R/helpers.R")

args <- commandArgs(trailingOnly = TRUE)
soft_power_override <- NULL
if (length(args) > 0L) {
  idx <- which(args == "--soft-power")
  if (length(idx) == 1L && idx < length(args)) {
    soft_power_override <- as.integer(args[idx + 1L])
  }
}

# When --soft-power is provided we treat the run as a sensitivity scan
# and write to suffixed cache filenames so the canonical (auto-power)
# build is preserved. The Rmd's soft-power-sensitivity section reads all
# three caches (canonical + p{N}) and compares module structure.
suffix <- if (is.null(soft_power_override)) "" else
  sprintf("_p%d", soft_power_override)
cache_micro <- sprintf("storage/cache/hdwgcna_microglia%s.rds", suffix)
cache_de    <- sprintf("storage/cache/hdwgcna_module_de%s.rds", suffix)

t_start <- Sys.time()
cat("[hdwgcna] loading microglia_seurat_processed.rds ...\n")
sc <- readRDS("storage/cache/microglia_seurat_processed.rds")
symbol_map <- readRDS("storage/cache/snrnaseq_symbol_map.rds")
cat(sprintf("[hdwgcna] loaded %d cells x %d genes in %.1f s.\n",
            ncol(sc), nrow(sc),
            as.numeric(difftime(Sys.time(), t_start, units = "secs"))))

cat("[hdwgcna] applying label_microglia_states() ...\n")
sc <- label_microglia_states(sc, symbol_map = symbol_map)
cat("\n==== state x genotype (metacell strata) ====\n")
print(table(sc$state, sc$genotype))

cat("\n[hdwgcna] running run_hdwgcna_pipeline() ...\n")
res <- run_hdwgcna_pipeline(
  sc,
  wgcna_name          = "microglia",
  group_by            = c("state", "genotype"),
  metacell_k          = 25,
  metacell_min_cells  = 50,
  metacell_max_shared = 10,
  metacell_reduction  = "harmony",
  gene_fraction       = 0.05,
  batch_var           = "batch",
  network_type        = "signed",
  tom_type            = "signed",
  min_module_size     = 50,
  merge_cut           = 0.2,
  deep_split          = 4,
  soft_power          = soft_power_override,
  verbose             = TRUE
)

cat(sprintf("\n[hdwgcna] pipeline complete. Soft power = %d. ",
            res$soft_power))
cat(sprintf("Modules detected: %d (incl. grey).\n",
            nrow(res$module_meta)))
print(res$module_meta)

# ---- Module DE ----
cat("\n[hdwgcna] fitting per-module DE under 2x2 factorial + batch ...\n")
de_res <- fit_module_de(
  MEs_cells    = res$MEs_cells,
  meta         = res$seurat_obj@meta.data,
  id_col       = "genotype_batch",
  genotype_col = "genotype",
  batch_col    = "batch"
)

cat("\n==== module DE summary (n_sig at FDR<0.05) ====\n")
for (ct in names(de_res$top)) {
  n_sig <- sum(de_res$top[[ct]]$adj.P.Val < 0.05, na.rm = TRUE)
  min_padj <- suppressWarnings(min(de_res$top[[ct]]$adj.P.Val, na.rm = TRUE))
  cat(sprintf("  %-15s n_sig=%2d min_padj=%.3g\n", ct, n_sig, min_padj))
}

# ---- Save lightweight cache ----
# Strip the full Seurat: keep only the wgcna slot + metadata for ME ->
# DE pipeline. This keeps the cache small (~50-100 MB instead of 1+ GB).
sc_meta <- res$seurat_obj@meta.data
sc_wgcna <- res$seurat_obj@misc$wgcna  # contains module assignments, TOM file paths, etc.

light_cache <- list(
  modules     = res$modules,
  hub_genes   = res$hub_genes,
  MEs_cells   = res$MEs_cells,
  module_meta = res$module_meta,
  power_table = res$power_table,
  soft_power  = res$soft_power,
  wgcna_slot  = sc_wgcna,
  cell_meta   = sc_meta[, c("genotype", "genotype_batch", "batch", "state"), drop = FALSE]
)
saveRDS(light_cache, cache_micro)
Sys.chmod(cache_micro, mode = "0664")
try(system2("chown", c("rstudio:rstudio", cache_micro)), silent = TRUE)

saveRDS(de_res, cache_de)
Sys.chmod(cache_de, mode = "0664")
try(system2("chown", c("rstudio:rstudio", cache_de)), silent = TRUE)

cat(sprintf("\n[hdwgcna] caches written:\n  %s (%.1f MB)\n  %s (%.1f MB)\n",
            cache_micro, file.info(cache_micro)$size / 1024^2,
            cache_de,    file.info(cache_de)$size / 1024^2))

cat(sprintf("[hdwgcna] total wall: %.1f min.\n",
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))

#!/usr/bin/env Rscript
# scripts/build_cellchat.R
#
# Per-condition CellChat v2 fits over the 9-category cell_type column
# (Astrocyte / Microglia_{DAM,IFN,homeostatic,proliferative} / Neuronal /
# OPC / Oligodendrocyte / Vascular), then a merged object for the 4
# genotypes (MAPTKI / P301S / NLGF_MAPTKI / NLGF_P301S).
#
# Comparative CCC analysis (rankNet, netVisual_diffInteraction,
# compareInteractions, etc.) all operate on the merged object; the
# per-condition list is the source of truth and is preserved so the
# merge can be re-done with different name orders or subsets without
# re-running the heavy computeCommunProb step.
#
# Inputs:
#   storage/cache/seurat_full_processed.rds   (built by build_seurat_full.R)
#
# Outputs:
#   storage/cache/cellchat_per_condition.rds  (named list of CellChat objs)
#   storage/cache/cellchat_merged.rds         (mergeCellChat output)
#
# Compute notes:
#   - future::plan("multisession", workers = 4) inside computeCommunProb's
#     permutation tests. Set in main(); CellChat picks it up automatically.
#   - With ~250k cells across 4 conditions and 9 cell types, expect
#     ~10-25 min per condition on 4 workers (~40-100 min total).
#   - Peak memory ~20-30 GB while the Seurat is loaded; drops to ~5 GB
#     once we extract data per condition.

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(CellChat)
  library(NMF)
  library(future)
})
source("R/helpers.R")

args        <- commandArgs(trailingOnly = TRUE)
force       <- "--force" %in% args
n_workers   <- 4L  # future workers for CellChat's bootstrap
nboot_perm  <- 100L

per_cond_path <- "storage/cache/cellchat_per_condition.rds"
merged_path   <- "storage/cache/cellchat_merged.rds"

if (file.exists(per_cond_path) && file.exists(merged_path) && !force) {
  message("Both output caches already exist. Pass --force to rebuild.")
  quit(status = 0)
}

t0 <- Sys.time()

# 0. Configure future for CellChat's parallel sections.
#    CellChat's computeCommunProb() calls future::future_sapply without
#    `future.seed=TRUE`, which produces "UNRELIABLE VALUE" warnings on every
#    bootstrap iteration. The permutations themselves are still statistically
#    independent (each worker draws from the host RNG); the warning just notes
#    that the parallel RNG isn't formally reproducible. Silencing it keeps the
#    log readable.
future::plan("multisession", workers = n_workers)
options(future.globals.maxSize  = 16 * 1024^3,   # 16 GB per-worker cap
        future.rng.onMisuse     = "ignore")
set.seed(42L)
message(sprintf("future workers: %d", n_workers))

# 1. Load full processed Seurat
message("[1/3] Loading storage/cache/seurat_full_processed.rds...")
sc <- readRDS("storage/cache/seurat_full_processed.rds")
message(sprintf("       %d genes x %d cells; %d cell types; %d genotypes",
                nrow(sc), ncol(sc),
                nlevels(sc$cell_type), nlevels(sc$genotype)))

# Drop empty cell_type factor levels (just in case)
sc$cell_type <- droplevels(sc$cell_type)

# 2. Per-condition CellChat fits
message("[2/3] Per-condition CellChat fits...")
conditions <- levels(sc$genotype)
cellchat_list <- list()

for (i in seq_along(conditions)) {
  cond <- conditions[i]
  cells_cond <- colnames(sc)[sc$genotype == cond]
  message(sprintf("  [%d/%d] %s: %d cells",
                  i, length(conditions), cond, length(cells_cond)))
  sub <- sc[, cells_cond]
  sub$cell_type <- droplevels(sub$cell_type)

  # Per-cell-type counts for this condition (used to verify min.cells passes)
  ct_counts <- table(sub$cell_type)
  message(sprintf("       cell_type counts: %s",
                  paste(sprintf("%s=%d", names(ct_counts), ct_counts),
                        collapse = ", ")))

  data_input <- GetAssayData(sub, assay = "RNA", layer = "data")
  meta       <- data.frame(group   = as.character(sub$cell_type),
                           samples = as.character(sub$genotype_batch),
                           row.names = colnames(sub))

  cc <- createCellChat(object  = data_input,
                       meta    = meta,
                       group.by = "group")
  cc@DB <- CellChatDB.mouse

  cc <- subsetData(cc)  # restrict to LR genes
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  # Bootstrap-based permutation test for LR probabilities
  cc <- computeCommunProb(cc,
                          type            = "triMean",
                          raw.use         = TRUE,
                          population.size = TRUE,
                          nboot           = nboot_perm)
  cc <- filterCommunication(cc, min.cells = 10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc <- netAnalysis_computeCentrality(cc, slot.name = "netP")

  cellchat_list[[cond]] <- cc

  # Free the per-condition working copy
  rm(sub, data_input, meta, cc); invisible(gc(verbose = FALSE))
}

# 3. Merge for comparative analysis
message("[3/3] Merging the 4 per-condition CellChat objects...")
cellchat_merged <- mergeCellChat(cellchat_list,
                                 add.names = names(cellchat_list))

# Save both
message("Saving caches...")
saveRDS(cellchat_list,   per_cond_path, compress = "gzip")
saveRDS(cellchat_merged, merged_path,   compress = "gzip")
Sys.chmod(per_cond_path, mode = "0644")
Sys.chmod(merged_path,   mode = "0644")

# Provenance
prov_path <- "storage/cache/cellchat_provenance.txt"
sink(prov_path)
cat("CellChat build provenance\n")
cat("=========================\n")
cat("Built at:           ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n", sep = "")
cat("Build script:       scripts/build_cellchat.R\n")
cat("R version:          ", as.character(getRversion()), "\n", sep = "")
cat("CellChat version:   ", as.character(packageVersion("CellChat")), "\n", sep = "")
cat("Future workers:     ", n_workers, "\n", sep = "")
cat("Permutation nboot:  ", nboot_perm, "\n", sep = "")
cat("\n--- Per-condition cell counts and cell-type breakdown ---\n")
for (cond in conditions) {
  cc <- cellchat_list[[cond]]
  cat(sprintf("\n[%s]  %d cells, %d cell types after filterCommunication\n",
              cond, ncol(cc@data.signaling), length(levels(cc@idents))))
  print(table(cc@idents))
}
cat("\n--- Merged object n_LR by condition ---\n")
print(sapply(cellchat_list, function(x) nrow(subsetCommunication(x))))
cat("\n--- Merged object structure ---\n")
print(cellchat_merged)
cat("\n--- Input md5sum ---\n")
cat("seurat_full_processed.rds: ",
    tools::md5sum("storage/cache/seurat_full_processed.rds"), "\n", sep = "")
sink()
Sys.chmod(prov_path, mode = "0644")

dt <- difftime(Sys.time(), t0, units = "mins")
message(sprintf("Wall time: %.1f min", as.numeric(dt)))
message("Done.")

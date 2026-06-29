#!/usr/bin/env Rscript
# scripts/build_seurat_full.R
#
# Build storage/cache/seurat_full_processed.rds: the full snRNAseq Seurat
# object with all broad cell types retained, microglia substate labels
# lifted from the processed microglia cache, and a unified 9-level
# `cell_type` column for downstream cell-cell communication (CellChat /
# MultiNicheNet).
#
# The microglia substate labels live on the post-QC processed cache
# (storage/cache/microglia_seurat_processed.rds). Microglia cells that
# failed substate QC (no entry in the processed cache) are dropped here;
# this keeps the CCC cell-type structure consistent with the rest of
# the project's microglia analyses. The 2 "Unknown" broad_annotation
# cells are also dropped.
#
# Output is normalised RNA (LogNormalize, scale.factor 10000) so both
# CellChat and MultiNicheNet can read the `data` slot directly.
#
# Inputs:
#   storage/data/snrnaseq.rds                   (9 GB; ~286k cells)
#   storage/cache/microglia_seurat_processed.rds (765 MB; ~23k microglia w/ state)
#
# Output:
#   storage/cache/seurat_full_processed.rds      (rstudio:rstudio, gzip)
#
# Wall time: ~10-20 min, peak memory ~22 GB. Run from project root.

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
})
source("R/helpers.R")  # constants$genotype_levels etc.

args     <- commandArgs(trailingOnly = TRUE)
force    <- "--force" %in% args
out_path <- "storage/cache/seurat_full_processed.rds"

if (file.exists(out_path) && !force) {
  message("Output already exists at ", out_path,
          " (pass --force to overwrite). Exiting.")
  quit(status = 0)
}

t0 <- Sys.time()

# 1. Compute microglia state labels (load processed cache, run
#    label_microglia_states(), extract, discard). The state column is
#    not pre-baked into the cache; the Rmd's snrnaseq-substate chunk
#    computes it via the same helper at knit time.
message("[1/9] Loading processed microglia cache for state computation...")
proc <- readRDS("storage/cache/microglia_seurat_processed.rds")
symbol_map <- readRDS("storage/cache/snrnaseq_symbol_map.rds")
message(sprintf("       Processed micro: %d genes x %d cells",
                nrow(proc), ncol(proc)))

message("       Computing state labels via label_microglia_states()...")
proc <- label_microglia_states(proc, symbol_map = symbol_map)
stopifnot("state" %in% colnames(proc@meta.data))
state_lookup <- setNames(as.character(proc$state), colnames(proc))
n_states  <- length(state_lookup)
state_tab <- table(state_lookup)
rm(proc); invisible(gc(verbose = FALSE))
message(sprintf("       Lifted %d microglia state labels.", n_states))
message("       State counts:")
print(state_tab)

# 2. Load full Seurat object
message("[2/9] Loading storage/data/snrnaseq.rds (~9 GB on disk)...")
sc <- readRDS("storage/data/snrnaseq.rds")
message(sprintf("       Full object: %d genes x %d cells",
                nrow(sc), ncol(sc)))

# Sanity checks
required_meta <- c("broad_annotations", "genotype", "batch",
                   "genotype_batch", "sex")
missing_meta <- setdiff(required_meta, colnames(sc@meta.data))
if (length(missing_meta)) {
  stop("Missing metadata columns: ",
       paste(missing_meta, collapse = ", "))
}

# 3. Drop SCT (was computed on all cells; we use RNA + LogNormalize)
message("[3/9] Dropping SCT assay and reductions...")
DefaultAssay(sc) <- "RNA"
sc[["SCT"]]    <- NULL
sc@reductions  <- list()
invisible(gc(verbose = FALSE))

# 4. Drop Unknown cells
n_unknown <- sum(sc$broad_annotations == "Unknown")
sc <- subset(sc, subset = broad_annotations != "Unknown")
message(sprintf("[4/9] Dropped %d 'Unknown' cells.", n_unknown))

# 5. Lift state labels onto matching microglia cells
cells_all <- colnames(sc)
sc$state  <- NA_character_
mg_in     <- intersect(cells_all, names(state_lookup))
sc$state[match(mg_in, cells_all)] <- state_lookup[mg_in]
message(sprintf("[5/9] Lifted state on %d microglia cells (of %d total microglia in broad_annotations).",
                length(mg_in), sum(sc$broad_annotations == "Microglia")))

# 6. Drop microglia that failed substate QC
mg_mask        <- sc$broad_annotations == "Microglia"
mg_no_state    <- mg_mask & is.na(sc$state)
n_drop_qc      <- sum(mg_no_state)
keep_mask      <- !mg_no_state
sc             <- sc[, keep_mask]
message(sprintf("[6/9] Dropped %d microglia cells without substate QC label.",
                n_drop_qc))

# 7. Build the 9-level cell_type column
sc$cell_type <- as.character(sc$broad_annotations)
mg_mask      <- sc$broad_annotations == "Microglia"
sc$cell_type[mg_mask] <- paste0("Microglia_", sc$state[mg_mask])
# Lock factor levels for deterministic ordering in downstream code
ct_levels <- c(sort(setdiff(unique(sc$cell_type), grep("^Microglia_", unique(sc$cell_type), value = TRUE))),
               sort(grep("^Microglia_", unique(sc$cell_type), value = TRUE)))
sc$cell_type <- factor(sc$cell_type, levels = ct_levels)

# Refresh genotype factor to canonical order from constants.R
sc$genotype <- factor(sc$genotype, levels = genotype_levels)

# 8. Ensembl ID -> gene symbol (CellChatDB and NicheNet priors are symbol-keyed).
#    The snrnaseq_symbol_map is a 1:1 bijection (verified: 33,683 rows, no NAs,
#    no duplicate symbols, no duplicate Ensembl IDs), so we can just relabel
#    rownames and rebuild the assay. Doing this once in the cache means
#    build_cellchat.R and build_multinichenet.R can read symbols directly.
message("[7/9] Mapping Ensembl IDs -> gene symbols (CellChat / NicheNet priors are symbol-based)...")
ens2sym <- setNames(symbol_map$symbol, symbol_map$ensembl)
stopifnot(setequal(rownames(sc), names(ens2sym)))      # require 1:1 coverage
counts <- LayerData(sc, assay = "RNA", layer = "counts")
rownames(counts) <- ens2sym[rownames(counts)]
md <- sc@meta.data
sc <- CreateSeuratObject(counts = counts, meta.data = md, assay = "RNA")
# CreateSeuratObject preserves data.frame factor levels, but re-assert
# explicitly so the cache always has the canonical orderings.
sc$cell_type <- factor(as.character(sc$cell_type), levels = ct_levels)
sc$genotype  <- factor(as.character(sc$genotype),  levels = genotype_levels)
rm(counts); invisible(gc(verbose = FALSE))
message(sprintf("       After symbol rename: %d genes x %d cells", nrow(sc), ncol(sc)))

# 9. LogNormalize on RNA
message("[8/9] LogNormalize RNA (scale.factor = 10000)...")
sc <- NormalizeData(sc,
                    normalization.method = "LogNormalize",
                    scale.factor         = 10000,
                    verbose              = FALSE)

# Final summary
message("[9/9] Final object: ", nrow(sc), " genes x ", ncol(sc), " cells")
message("       cell_type counts:")
print(table(sc$cell_type))
message("       genotype x batch:")
print(table(sc$genotype, sc$batch))

# Save
message("Saving to ", out_path, " ...")
saveRDS(sc, out_path, compress = "gzip")
Sys.chmod(out_path, mode = "0644")
sz <- file.size(out_path)
message(sprintf("       File size: %.1f MB", sz / 1024^2))

# Provenance: write a small companion txt with key counts + md5
prov_path <- sub("\\.rds$", "_provenance.txt", out_path)
sink(prov_path, append = FALSE)
cat("seurat_full_processed.rds provenance\n")
cat("====================================\n")
cat("Built at:           ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n", sep = "")
cat("Build script:       scripts/build_seurat_full.R\n")
cat("R version:          ", as.character(getRversion()), "\n", sep = "")
cat("Seurat version:     ", as.character(packageVersion("Seurat")), "\n", sep = "")
cat("\n--- Input md5sums ---\n")
cat("snrnaseq.rds (source):                ", tools::md5sum("storage/data/snrnaseq.rds"), "\n", sep = "")
cat("microglia_seurat_processed.rds (state):", tools::md5sum("storage/cache/microglia_seurat_processed.rds"), "\n", sep = "")
cat("\n--- Cell type counts ---\n")
print(table(sc$cell_type))
cat("\n--- genotype x batch ---\n")
print(table(sc$genotype, sc$batch))
cat("\n--- Final object dimensions ---\n")
cat("Genes: ", nrow(sc), "\n", sep = "")
cat("Cells: ", ncol(sc), "\n", sep = "")
sink()
Sys.chmod(prov_path, mode = "0644")

dt <- difftime(Sys.time(), t0, units = "mins")
message(sprintf("Wall time: %.1f min", as.numeric(dt)))
message("Done.")

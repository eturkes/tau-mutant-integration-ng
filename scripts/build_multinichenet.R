#!/usr/bin/env Rscript
# scripts/build_multinichenet.R
#
# Multi-sample, multi-condition MultiNicheNet analysis over the same 9
# cell-type categories used by build_cellchat.R and the same 5-contrast
# factorial structure used by every other DE step in this project.
#
# Design:
#   sample_id    = genotype_batch (16 levels, one per pseudo-sample)
#   group_id     = genotype       (4 levels)
#   celltype_id  = cell_type      (9 levels, microglia split into 4 substates)
#   batches      = batch          (covariate; nested within genotype_batch but
#                                  spans groups; standard MultiNicheNet usage)
#   contrasts    = the 5 used elsewhere in this project:
#                  tau_alone        = P301S - MAPTKI
#                  nlgf_in_maptki   = NLGF_MAPTKI - MAPTKI
#                  nlgf_in_p301s    = NLGF_P301S - P301S
#                  tau_in_nlgf      = NLGF_P301S - NLGF_MAPTKI
#                  interaction      = (NLGF_P301S - P301S) - (NLGF_MAPTKI - MAPTKI)
#
# Pseudobulk DE under the hood is muscat; ligand-activity scoring uses
# the cached NicheNet mouse prior at storage/data/nichenet/.
#
# Input:
#   storage/cache/seurat_full_processed.rds  (built by build_seurat_full.R)
#   storage/data/nichenet/lr_network_mouse_21122021.rds
#   storage/data/nichenet/ligand_target_matrix_nsga2r_final_mouse.rds
#
# Output:
#   storage/cache/multinichenet_output.rds
#
# Compute: ~30-90 min on 8 cores. Peak memory ~15-25 GB.

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(SingleCellExperiment)
  library(multinichenetr)
  library(nichenetr)
  library(dplyr)
  library(tibble)
})
source("R/helpers.R")

args  <- commandArgs(trailingOnly = TRUE)
force <- "--force" %in% args
out_path <- "storage/cache/multinichenet_output.rds"

if (file.exists(out_path) && !force) {
  message("Output already exists at ", out_path,
          " (pass --force to overwrite). Exiting.")
  quit(status = 0)
}

t0 <- Sys.time()

# 1. Load full processed Seurat -> SingleCellExperiment
message("[1/4] Loading storage/cache/seurat_full_processed.rds...")
sc <- readRDS("storage/cache/seurat_full_processed.rds")

# MultiNicheNet expects the colData fields as named character vectors
sc$sample_id   <- as.character(sc$genotype_batch)
sc$group_id    <- as.character(sc$genotype)
sc$celltype_id <- as.character(sc$cell_type)
sc$batch_id    <- as.character(sc$batch)

# Sanity check: every sample has only one genotype and one batch
sample_check <- sc@meta.data |>
  dplyr::distinct(sample_id, group_id, batch_id)
stopifnot(nrow(sample_check) == length(unique(sc$sample_id)))
message(sprintf("       %d cells, %d samples, %d groups, %d batches, %d cell types",
                ncol(sc), length(unique(sc$sample_id)),
                length(unique(sc$group_id)), length(unique(sc$batch_id)),
                length(unique(sc$celltype_id))))

# Coerce to SingleCellExperiment
message("       Converting to SingleCellExperiment...")
sce <- Seurat::as.SingleCellExperiment(sc, assay = "RNA")
rm(sc); invisible(gc(verbose = FALSE))

# Cache rownames are gene symbols (set in build_seurat_full.R via the
# 1:1 snrnaseq_symbol_map). ~248 entries in the map have ensembl == symbol
# (genes without an MGI symbol), so we use a 5% proportion threshold here
# rather than a strict zero-tolerance check.
n_ens    <- sum(grepl("^ENSMUSG[0-9]+$", rownames(sce)))
prop_ens <- n_ens / nrow(sce)
if (prop_ens > 0.05) {
  stop(sprintf(
    "Expected symbol-keyed cache but %.1f%% of rownames (%d / %d) look like Ensembl IDs. ",
    100 * prop_ens, n_ens, nrow(sce)),
    "Re-run scripts/build_seurat_full.R --force to rebuild with symbols.")
}
message(sprintf("       SCE: %d genes (symbols, %d ENSMUSG-style fall-throughs) x %d cells",
                nrow(sce), n_ens, ncol(sce)))

# 2. Load NicheNet priors
message("[2/4] Loading NicheNet priors (mouse)...")
lr_network <- readRDS("storage/data/nichenet/lr_network_mouse_21122021.rds") |>
  dplyr::rename(from = from, to = to) |>
  dplyr::distinct(from, to)
ligand_target_matrix <- readRDS("storage/data/nichenet/ligand_target_matrix_nsga2r_final_mouse.rds")
message(sprintf("       lr_network: %d unique LR pairs", nrow(lr_network)))
message(sprintf("       ligand_target_matrix: %d targets x %d ligands",
                nrow(ligand_target_matrix), ncol(ligand_target_matrix)))

# 3. Set up contrasts (limma-style linear combinations of group levels)
message("[3/4] Configuring MultiNicheNet contrasts...")
# group levels (must match levels of group_id)
group_levels <- c("MAPTKI", "P301S", "NLGF_MAPTKI", "NLGF_P301S")
# Order is set to match genotype_levels in constants.R

# multinichenet requires contrasts_oi as a single comma-joined string with
# each contrast wrapped in single quotes (parsed at L80 / L82-87 of the
# function source). Building it via paste() keeps the per-contrast format
# readable.
contrasts_oi <- paste(
  c("'P301S-MAPTKI'",
    "'NLGF_MAPTKI-MAPTKI'",
    "'NLGF_P301S-P301S'",
    "'NLGF_P301S-NLGF_MAPTKI'",
    "'(NLGF_P301S-P301S)-(NLGF_MAPTKI-MAPTKI)'"),
  collapse = ","
)
# contrast_tbl$group must contain actual group_id (genotype) levels, not
# semantic labels. For each "A-B" contrast we set group = A (the up-condition);
# for the interaction we use NLGF_P301S as the "divergent" condition of
# interest. multinichenet warns on duplicate groups but accepts them.
contrast_tbl <- tibble::tibble(
  contrast = c(
    "P301S-MAPTKI",
    "NLGF_MAPTKI-MAPTKI",
    "NLGF_P301S-P301S",
    "NLGF_P301S-NLGF_MAPTKI",
    "(NLGF_P301S-P301S)-(NLGF_MAPTKI-MAPTKI)"
  ),
  group = c(
    "P301S",       # tau_alone:        P301S      vs MAPTKI
    "NLGF_MAPTKI", # nlgf_in_maptki:   NLGF_MAPTKI vs MAPTKI
    "NLGF_P301S",  # nlgf_in_p301s:    NLGF_P301S vs P301S
    "NLGF_P301S",  # tau_in_nlgf:      NLGF_P301S vs NLGF_MAPTKI
    "NLGF_P301S"   # interaction:      divergent NLGF effect (NLGF x P301S)
  )
)
message("       Contrasts to fit:")
print(contrast_tbl)

# 4. Run multi_nichenet_analysis (wrapper around the full MultiNicheNet pipeline)
message("[4/4] Running multi_nichenet_analysis() ...")

# Senders / receivers: all cell types as both (full bidirectional analysis)
senders_oi   <- unique(as.character(colData(sce)$celltype_id))
receivers_oi <- unique(as.character(colData(sce)$celltype_id))

# multinichenetr 2.1.0 dropped the prioritizing_weights argument; weights are
# selected through the `scenario` argument ("regular", "lower_DE",
# "no_frac_LR_expr"). "regular" gives the standard tutorial weights.
# Defaults: min_cells = 10, fraction_cutoff = 0.05, logFC_threshold = 0.5.
mnn_output <- multi_nichenet_analysis(
  sce                           = sce,
  celltype_id                   = "celltype_id",
  sample_id                     = "sample_id",
  group_id                      = "group_id",
  batches                       = "batch_id",
  covariates                    = NA,
  lr_network                    = lr_network,
  ligand_target_matrix          = ligand_target_matrix,
  contrasts_oi                  = contrasts_oi,
  contrast_tbl                  = contrast_tbl,
  senders_oi                    = senders_oi,
  receivers_oi                  = receivers_oi,
  scenario                      = "regular",
  min_cells                     = 10,
  fraction_cutoff               = 0.05,
  logFC_threshold               = 0.5,
  p_val_threshold               = 0.05,
  p_val_adj                     = FALSE,  # raw p; multinichenet does its own correction
  empirical_pval                = FALSE,
  top_n_target                  = 250,
  verbose                       = TRUE,
  n.cores                       = as.numeric(parallel::detectCores() - 1L)
)

# Save output
message("Saving MultiNicheNet output to ", out_path, " ...")
saveRDS(mnn_output, out_path, compress = "gzip")
Sys.chmod(out_path, mode = "0644")
sz <- file.size(out_path)
message(sprintf("       File size: %.1f MB", sz / 1024^2))

# Provenance
prov_path <- sub("\\.rds$", "_provenance.txt", out_path)
sink(prov_path)
cat("MultiNicheNet build provenance\n")
cat("==============================\n")
cat("Built at:           ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n", sep = "")
cat("Build script:       scripts/build_multinichenet.R\n")
cat("multinichenetr ver: ", as.character(packageVersion("multinichenetr")), "\n", sep = "")
cat("nichenetr ver:      ", as.character(packageVersion("nichenetr")), "\n", sep = "")
cat("\n--- Output object slot summary ---\n")
cat("Slot names: ", paste(names(mnn_output), collapse = ", "), "\n")
if (!is.null(mnn_output$prioritization_tables)) {
  cat("\nprioritization_table_with_all_scores n_rows: ",
      tryCatch(nrow(mnn_output$prioritization_tables$group_prioritization_tbl),
               error = function(e) NA), "\n")
}
cat("\n--- Input md5sums ---\n")
cat("seurat_full_processed.rds: ",
    tools::md5sum("storage/cache/seurat_full_processed.rds"), "\n", sep = "")
cat("lr_network_mouse_21122021.rds: ",
    tools::md5sum("storage/data/nichenet/lr_network_mouse_21122021.rds"), "\n", sep = "")
cat("ligand_target_matrix_nsga2r_final_mouse.rds: ",
    tools::md5sum("storage/data/nichenet/ligand_target_matrix_nsga2r_final_mouse.rds"), "\n", sep = "")
cat("\n--- Contrast tbl ---\n")
print(contrast_tbl)
sink()
Sys.chmod(prov_path, mode = "0644")

dt <- difftime(Sys.time(), t0, units = "mins")
message(sprintf("Wall time: %.1f min", as.numeric(dt)))
message("Done.")

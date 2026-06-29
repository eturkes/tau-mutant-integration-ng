#!/usr/bin/env Rscript
# Build the LIANA+ cell-cell communication cache for Phase G of the
# mechanism layer (see storage/notes/mechanism_layer_plan.md, sessions
# G2..G5). Heavy compute lives here, OUTSIDE the knit, per the project
# convention in scripts/build_*.R; the Rmd consumer (G3+) `readRDS`s
# the resulting cache directly.
#
# Idempotent: skips the cache write if the output already exists
# unless `--overwrite` is passed. Re-running with no args after a
# successful build is a no-op. Mirrors the
# `scripts/build_kinase_activity_decoupler.R` convention.
#
# Inputs (rds, under storage/cache/):
#   seurat_full_processed.rds   the project-wide Seurat object with
#                               `$genotype`, `$cell_type`, and the RNA
#                               assay's `data` layer (log-normalised).
#
# Python dependencies (in `.venv/`):
#   liana       1.7.1  (multi-method CCC consensus framework;
#                       Dimitrov *et al.* 2024, Nat Cell Biol)
#   anndata     >=0.10 (in-memory single-cell container)
#   scipy       >=1.13 (sparse matrix construction)
#   pandas      >=2.2  (anndata.obs / .var)
#   pip install path: `./.venv/bin/pip install liana==1.7.1`
# The venv was created via `python3 -m venv .venv` and is local to the
# project root; reticulate is configured to use `.venv/bin/python`
# (see `R/ccc_inference.R::.LIANA_VENV` and `.ensure_liana_python()`).
# Symlink note: a venv `python` binary is a symlink chain to the
# system `python3`; we deliberately do NOT `normalizePath()` the venv
# path because following the symlink destroys the venv context (the
# `pyvenv.cfg` marker lives next to the venv symlink, not next to the
# system binary). See `R/ccc_inference.R::.ensure_liana_python()` for
# the bootstrap details.
#
# RETICULATE_PYTHON env var: this script sets RETICULATE_PYTHON BEFORE
# `library(reticulate)` is loaded. Python can only be initialised once
# per R session, so the env var has to be in place at library-load
# time; the bootstrap helper in `.ensure_liana_python()` is the
# session-level safety net.
#
# Output (under storage/cache/):
#   liana_output.rds
#   Nested list with three top-level slots:
#     per_genotype  named list (length 4: MAPTKI, P301S, NLGF_MAPTKI,
#                   NLGF_P301S) of tibbles carrying the full liana_res
#                   schema (13 cols) PLUS `cellphone_padj` (BH-corrected
#                   `cellphone_pvals` within each per-genotype run).
#     per_contrast  named list (length 5: tau_alone, nlgf_in_maptki,
#                   nlgf_in_p301s, tau_in_nlgf, interaction) of tibbles
#                   with columns sender / receiver / ligand / receptor /
#                   lr_interaction / prioritization_score /
#                   n_methods_pass_significance / sign_dir, plus
#                   auxiliary diagnostic cols (per-genotype magnitude
#                   ranks and cellphone_padjs). Zero NAs in the
#                   significance/prio/sign columns by construction
#                   (missing LR pairs imputed as worst-rank / worst-p
#                   after the four-way outer-join; see module header
#                   in `R/ccc_inference.R`).
#     provenance    list of liana / anndata / scipy / python versions,
#                   run timestamp, parameter values, n_pairs_per_genotype,
#                   and the n_cells_per_genotype snapshot.
#
# Runtime: ~6 s for the 2-genotype/3-celltype smoke test on 4000 cells
# (in iterative development). Full run extrapolation: 4 genotypes ×
# 9 cell types × ~17,500 cells each ≈ 5-10 minutes wall on the
# project's Docker container. Memory peak: ~10-20 GB while the full
# Seurat is loaded; drops to ~3-5 GB once subsetting begins.
#
# Phase G is single-modality (snRNAseq only) because CCC inference
# requires cell-resolved expression that only snRNAseq provides in
# this project (GeoMx is spot-level, bulk modalities are not cell-
# resolved). So this script has no analogue of the TF script's
# `modality_top_lists` or the kinase script's `cache_top_lists`; it
# operates directly on the per-genotype Seurat subsets. The cross-tool
# axis (CellChat / MultiNicheNet / liana) is the analogue of the Phase
# E cross-modality axis at G3+.

# RETICULATE_PYTHON must be set BEFORE reticulate loads.
project_root <- "/home/rstudio/tau-mutant-integration-ng"
Sys.setenv(RETICULATE_PYTHON = file.path(project_root, ".venv/bin/python"))

suppressPackageStartupMessages({
  library(reticulate)
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(tibble)
})

args      <- commandArgs(trailingOnly = TRUE)
overwrite <- "--overwrite" %in% args

setwd(project_root)
source("R/helpers.R")  # pulls in R/ccc_inference.R via the project loader

cache_dir <- "storage/cache"
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
out_path  <- file.path(cache_dir, "liana_output.rds")

if (file.exists(out_path) && !overwrite) {
  cat(sprintf("[build_liana_cache] cache exists, skipping: %s\n", out_path))
  cat("Pass --overwrite to rebuild.\n")
  quit(save = "no", status = 0)
}

cat("[build_liana_cache] bootstrapping liana python venv...\n")
.ensure_liana_python()
cat(sprintf("  python = %s\n", reticulate::py_config()$python))
li <- reticulate::import("liana")
cat(sprintf("  liana version = %s\n", li[["__version__"]]))

cat("[build_liana_cache] loading project-wide Seurat...\n")
t0 <- proc.time()
sc <- readRDS(file.path(cache_dir, "seurat_full_processed.rds"))
t1 <- proc.time()
cat(sprintf("  %d genes x %d cells; %d cell types; %d genotypes (loaded in %.1fs)\n",
            nrow(sc), ncol(sc),
            nlevels(sc$cell_type), nlevels(sc$genotype),
            (t1 - t0)[3]))
sc$cell_type <- droplevels(sc$cell_type)
sc$genotype  <- droplevels(sc$genotype)

# Snapshot cells per genotype for provenance / log. Order by the
# canonical project genotype levels (not alphabetical) so the
# provenance log reads MAPTKI -> P301S -> NLGF_MAPTKI -> NLGF_P301S.
n_cells_per_genotype <- table(
  factor(as.character(sc$genotype), levels = .LIANA_GENOTYPES)
)
cat("[build_liana_cache] cells per genotype:\n")
print(n_cells_per_genotype)

# Parameter set. Captured here for the provenance slot and the log;
# kept in sync with the defaults of run_liana_rank_aggregate_for_genotype.
params <- list(
  resource_name = "mouseconsensus",
  use_raw       = FALSE,
  expr_prop     = 0.10,
  min_cells     = 5L,
  n_perms       = 100L,
  seed          = 1L,
  sig_alpha     = 0.10,
  cell_types    = .LIANA_CELL_TYPES,
  genotypes     = .LIANA_GENOTYPES
)
cat("[build_liana_cache] params:\n")
str(params)

cat("\n[build_liana_cache] per-genotype rank_aggregate runs (4 genotypes)...\n")
t_start <- proc.time()
per_genotype <- build_liana_per_genotype_list(
  sc            = sc,
  genotypes     = params$genotypes,
  cell_types    = params$cell_types,
  resource_name = params$resource_name,
  use_raw       = params$use_raw,
  expr_prop     = params$expr_prop,
  min_cells     = params$min_cells,
  n_perms       = params$n_perms,
  seed          = params$seed,
  verbose       = TRUE
)
t_end <- proc.time()
cat(sprintf("[build_liana_cache] all per-genotype runs done in %.1fs\n",
            (t_end - t_start)[3]))

cat("\n[build_liana_cache] deriving per-contrast tibbles (5 contrasts)...\n")
per_contrast <- derive_liana_per_contrast_list(
  per_genotype   = per_genotype,
  contrast_specs = .LIANA_CONTRAST_SPECS,
  sig_alpha      = params$sig_alpha
)
for (cn in names(per_contrast)) {
  pc <- per_contrast[[cn]]
  sig <- sum(pc$n_methods_pass_significance == 1L, na.rm = TRUE)
  cat(sprintf("  %-15s  %5d LR-pair rows  (%d with cellphone_padj<%.2f in primary)\n",
              cn, nrow(pc), sig, params$sig_alpha))
}

cat("\n[build_liana_cache] building provenance...\n")
provenance <- build_liana_provenance(per_genotype = per_genotype, params = params)
provenance$n_cells_per_genotype <- as.integer(n_cells_per_genotype)
names(provenance$n_cells_per_genotype) <- names(n_cells_per_genotype)
str(provenance, max.level = 1)

liana_output <- list(
  per_genotype = per_genotype,
  per_contrast = per_contrast,
  provenance   = provenance
)

cat(sprintf("\n[build_liana_cache] writing cache: %s\n", out_path))
saveRDS(liana_output, out_path, compress = "xz")
Sys.chmod(out_path, mode = "0644")

cat("[build_liana_cache] done.\n")

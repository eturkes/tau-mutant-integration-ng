#!/usr/bin/env Rscript
# Build GO MF and GO CC fgsea result caches + per-contrast TSV exports.
#
# Phase B Step B1 of the pathway-overhaul plan. The GO BP collection is
# already built inside `rmd/07_integration.Rmd`; this script extends the
# project's pathway-collection coverage to GO MF and GO CC so the agnostic
# survey introduced in Phase C (`rmd/12_pathway_survey.Rmd`) can rank
# pathways across three GO sub-ontologies plus the custom curated sets.
#
# Idempotent: skips any cache that already exists unless `--overwrite` is
# passed. Re-running with no args after a successful build is a no-op.
#
# Inputs (rds, all under storage/cache/):
#   de_snrnaseq_nebula.rds  de_geomx.rds  de_proteomics.rds
#   de_phospho.rds          de_phospho_corrected.rds
#
# Outputs (under storage/cache/):
#   msigdb_gomf_mouse.rds   msigdb_gocc_mouse.rds
#   fgsea_gomf_results.rds  fgsea_gocc_results.rds
#
# Outputs (under storage/results/):
#   fgsea_gomf_per_contrast.tsv
#   fgsea_gocc_per_contrast.tsv
#
# Cache schema mirrors `fgsea_gobp_results.rds`:
#   list(modality -> list(contrast -> fgseaResult data.table))
# with modalities = (snrnaseq, geomx, proteomics, phospho, phospho_corr)
# and contrasts   = (nlgf_in_maptki, nlgf_in_p301s, interaction,
#                    tau_alone, tau_in_nlgf).

suppressPackageStartupMessages({
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
overwrite <- "--overwrite" %in% args

setwd("/home/rstudio/tau-mutant-integration-ng")
source("R/helpers.R")

cache_dir   <- "storage/cache"
results_dir <- "storage/results"
dir.create(cache_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# Wrap snRNAseq NEBULA top tables in a fit-shaped list so the shared
# `run_fgsea_per_dataset` helper (which expects `$fit$top[[contrast]]`)
# works without a NEBULA-specific branch. Mirrors `snrnaseq_nebula_shim`
# in rmd/07_integration.Rmd.
shim_nebula <- function(nebula_obj) list(fit = list(top = nebula_obj$top))

# Build per-modality fgsea results for one collection, write the cache and
# the joined-per-contrast TSV. `pathways_loader` is a 0-arg function (e.g.
# `get_gomf`) that returns the pathway-name -> gene-symbols list.
build_collection <- function(label, pathways_loader,
                             fgsea_cache_path, tsv_path,
                             de_snrnaseq_nebula, de_geomx, de_proteomics,
                             de_phospho, de_phospho_corr) {
  if (!overwrite && file.exists(fgsea_cache_path) && file.exists(tsv_path)) {
    message(sprintf("[skip] %s: cache and TSV already present", label))
    return(invisible(NULL))
  }
  message(sprintf("[build] %s", label))
  pathways <- pathways_loader()
  message(sprintf("  collection size: %d gene sets", length(pathways)))

  results <- cache_or_run(fgsea_cache_path, {
    list(
      snrnaseq     = run_fgsea_per_dataset(shim_nebula(de_snrnaseq_nebula),
                                           pathways),
      geomx        = run_fgsea_per_dataset(de_geomx,        pathways),
      proteomics   = run_fgsea_per_dataset(de_proteomics,   pathways),
      phospho      = run_fgsea_per_dataset(de_phospho,      pathways),
      phospho_corr = run_fgsea_per_dataset(de_phospho_corr, pathways)
    )
  }, overwrite = overwrite)

  joint <- join_fgsea_results(results)
  write_tsv_safe(joint, tsv_path)
  message(sprintf("  wrote: %s  (%d pathway-contrast rows)",
                  tsv_path, nrow(joint)))
  invisible(NULL)
}

# Load DE caches once; reuse across collections.
message("loading DE caches")
de_snrnaseq_nebula <- readRDS(file.path(cache_dir, "de_snrnaseq_nebula.rds"))
de_geomx           <- readRDS(file.path(cache_dir, "de_geomx.rds"))
de_proteomics      <- readRDS(file.path(cache_dir, "de_proteomics.rds"))
de_phospho         <- readRDS(file.path(cache_dir, "de_phospho.rds"))
de_phospho_corr    <- readRDS(file.path(cache_dir, "de_phospho_corrected.rds"))

build_collection(
  label              = "GO MF",
  pathways_loader    = get_gomf,
  fgsea_cache_path   = file.path(cache_dir,   "fgsea_gomf_results.rds"),
  tsv_path           = file.path(results_dir, "fgsea_gomf_per_contrast.tsv"),
  de_snrnaseq_nebula = de_snrnaseq_nebula,
  de_geomx           = de_geomx,
  de_proteomics      = de_proteomics,
  de_phospho         = de_phospho,
  de_phospho_corr    = de_phospho_corr
)

build_collection(
  label              = "GO CC",
  pathways_loader    = get_gocc,
  fgsea_cache_path   = file.path(cache_dir,   "fgsea_gocc_results.rds"),
  tsv_path           = file.path(results_dir, "fgsea_gocc_per_contrast.tsv"),
  de_snrnaseq_nebula = de_snrnaseq_nebula,
  de_geomx           = de_geomx,
  de_proteomics      = de_proteomics,
  de_phospho         = de_phospho,
  de_phospho_corr    = de_phospho_corr
)

message("done")

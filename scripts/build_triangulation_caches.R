#!/usr/bin/env Rscript
# Pre-build the snRNAseq triangulation caches outside the Rmd so the knit
# stage can run end-to-end with quick cache reads.

setwd("/home/rstudio/tau-mutant-integration-ng")
suppressPackageStartupMessages({
  source("R/helpers.R")
})

cache_dir <- "storage/cache"

stage <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(stage)) stage <- "all"
cat(sprintf("[build] stage = %s\n", stage))

if (stage %in% c("pb", "all")) {
  message("[build] loading microglia + symbol map")
  sc          <- readRDS(file.path(cache_dir, "microglia_seurat_processed.rds"))
  symbol_map  <- readRDS(file.path(cache_dir, "snrnaseq_symbol_map.rds"))

  sn_tri_pb_path <- file.path(cache_dir, "de_snrnaseq_triangulation_pb.rds")
  cache_or_run(sn_tri_pb_path, overwrite = TRUE, {
    pb_obj <- build_pseudobulk(sc, sample_col = "genotype_batch",
                               covariate_cols = c("genotype", "batch", "sex"))
    meta <- pb_obj$meta
    meta$genotype <- factor(as.character(meta$genotype), levels = genotype_levels)
    meta$batch    <- factor(meta$batch)
    meta$tau      <- as.integer(meta$genotype %in% c("P301S", "NLGF_P301S"))
    meta$nlgf     <- as.integer(meta$genotype %in% c("NLGF_MAPTKI", "NLGF_P301S"))
    meta$tau_nlgf <- meta$tau * meta$nlgf
    fd <- factorial_design(meta, "genotype", "batch")
    message("[build] edger_qlf")
    edger  <- fit_edger_qlf(pb_obj$counts, fd$design, fd$contrasts,
                            min_count = 10, symbol_map = symbol_map)
    message("[build] deseq2")
    deseq2 <- fit_deseq2_pb(pb_obj$counts, meta, fd$design, fd$contrasts,
                            min_count = 10, symbol_map = symbol_map)
    message("[build] dream")
    dream  <- fit_dream_pb(pb_obj$counts, meta, fd$contrasts,
                           min_count = 10, symbol_map = symbol_map)
    list(edger = edger, deseq2 = deseq2, dream = dream)
  })
  message("[build] PB triangulation cache written.")
}

if (stage %in% c("glmmtmb", "all")) {
  if (!exists("sc")) sc <- readRDS(file.path(cache_dir,
                                              "microglia_seurat_processed.rds"))
  if (!exists("symbol_map"))
    symbol_map <- readRDS(file.path(cache_dir, "snrnaseq_symbol_map.rds"))

  sn_tri_tmb_path <- file.path(cache_dir, "de_snrnaseq_glmmtmb.rds")
  cache_or_run(sn_tri_tmb_path, overwrite = TRUE, {
    fit_glmmtmb_microglia(
      sc,
      id_col        = "genotype_batch",
      genotype_col  = "genotype",
      assay         = "RNA",
      layer         = "counts",
      min_cell_frac = 0.01,
      ncore         = max(1L, parallel::detectCores() - 1L),
      symbol_map    = symbol_map
    )
  })
  message("[build] glmmTMB cache written.")
}

message("[build] all stages complete.")

#!/usr/bin/env Rscript
# L3 driver (plan arc L): deconvolve the 91-ROI GeoMx-WTA tissue into
# cell-type / microglial-substate ABUNDANCES against the snRNAseq reference
# and push those abundances through the project's locked 5-contrast 2x2
# factorial. Thin I/O wrapper over R/spatial_decon.R (pure fns); mirrors the
# build_scenic_contrasts.R pattern -- all heavy compute OUT of the knit,
# rmd/21 (L4) is display-only. Idempotent unless --overwrite.
#
# This is the project's first CELL-COMPOSITION readout: every prior arc reads
# expression or activity, none reads HOW MUCH of each cell type sits in tissue.
# Two-stage granularity (user gate, L1): stage-1 6-level (5 broad + pooled
# Microglia = robust total) -> stage-2 sub-deconvolve the 4 substates and
# anchor within-microglia fractions to the stage-1 total.
#
# Heavy step: loads seurat_full_processed.rds (286k cells, ~1.9 G on disk);
# create_profile_matrix densifies the capped per-type subset internally
# (~3.7 G at cap 3000 over the 18,512 common genes). Caveats are carried in
# R/spatial_decon.R's header and re-stated by rmd/21 at every interpretation.
#
# Outputs (storage/results/*.tsv + storage/cache/spatial_decon.rds):
#   spatial_decon_abundance_by_genotype.tsv  per-genotype mean abundance+prop
#   spatial_decon_contrasts.tsv              5-contrast factorial, both layers
#   spatial_decon_spatial_autocorr.tsv       per-slide Moran's I (raw + resid)

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tibble)
})

setwd("/home/rstudio/tau-mutant-integration-ng")
source("R/helpers.R")

args      <- commandArgs(trailingOnly = TRUE)
overwrite <- "--overwrite" %in% args
cc  <- "storage/cache"
res <- "storage/results"
dir.create(res, recursive = TRUE, showWarnings = FALSE)
out_cache <- file.path(cc, "spatial_decon.rds")

## ---- embedded params (state thresholds before applying; L1/L2 locked) ----
CAP        <- 3000L     # cells/type cap before create_profile_matrix densifies
MINCELLNUM <- 15L       # create_profile_matrix per-type min cells
MINGENES   <- 100L      # create_profile_matrix per-gene min expressing cells
OFFSET     <- 1         # log(beta + offset); one-cell pseudocount (L2)
PADJ_CUT   <- 0.10      # project-standard FDR
N_PERM     <- 999L      # Moran's I permutation count
SEED       <- 1L

if (file.exists(out_cache) && !overwrite) {
  cat(sprintf("[build_spatial_deconvolution] cache exists, skipping: %s\n", out_cache))
  cat("Pass --overwrite to rebuild.\n"); quit(save = "no", status = 0)
}

say <- function(...) cat(sprintf(...), "\n")
wtsv <- function(df, name) {
  p <- file.path(res, name)
  write_tsv_safe(as.data.frame(df), p); Sys.chmod(p, "0644")
  say("  wrote %-44s %d rows", name, nrow(df))
}

## ---- load inputs --------------------------------------------------------
say("[build_spatial_deconvolution] loading inputs (seurat_full is heavy)...")
sc    <- readRDS(file.path(cc, "seurat_full_processed.rds"))
geomx <- readRDS("storage/data/geomx.rds")
ct9   <- as.character(sc@meta.data[["cell_type"]])
sc@meta.data$cell_type6 <- ifelse(grepl("^Microglia_", ct9), "Microglia", ct9)
geomx_symbols <- rownames(geomx)
say("  reference %d cells x %d genes | GeoMx %d genes x %d ROIs",
    ncol(sc), nrow(sc), nrow(geomx), ncol(geomx))
say("  9-level cell_type: %s", paste(sort(unique(ct9)), collapse = ", "))
say("  6-level cell_type6: %s",
    paste(sort(unique(sc@meta.data$cell_type6)), collapse = ", "))

# Clean ROI meta (clean names so rmd/21 + the design avoid the "slide name"
# space); rownames = ROI ids = colnames of every abundance matrix below.
md <- geomx@meta.data
roi_meta <- data.frame(
  roi      = colnames(geomx),
  genotype = factor(as.character(md$genotype), levels = genotype_levels),
  slide    = as.character(md[["slide name"]]),
  x        = as.numeric(md[["ROI Coordinate X"]]),
  y        = as.numeric(md[["ROI Coordinate Y"]]),
  nuclei   = as.numeric(md$nuclei),
  row.names = colnames(geomx), stringsAsFactors = FALSE)

## ---- build both reference profiles --------------------------------------
say("[1/6] building reference profiles (cap=%d cells/type)...", CAP)
profile6 <- build_reference_profile(sc, "cell_type6", geomx_symbols,
                                    cap = CAP, seed = SEED,
                                    minCellNum = MINCELLNUM, minGenes = MINGENES)
profile9 <- build_reference_profile(sc, "cell_type", geomx_symbols,
                                    cap = CAP, seed = SEED,
                                    minCellNum = MINCELLNUM, minGenes = MINGENES)
say("  common genes (ref n GeoMx) = %d  [L1 expects 18,512]",
    profile6$qc$n_common_genes)
say("  profile6: %d genes x %d types | kappa=%.1f | max|cor|=%.3f",
    profile6$qc$n_profile_genes, length(profile6$qc$cell_types),
    profile6$qc$condition_number, profile6$qc$max_abs_celltype_cor)
say("  profile9: %d genes x %d types | kappa=%.1f | max|cor|=%.3f  (substate collinearity)",
    profile9$qc$n_profile_genes, length(profile9$qc$cell_types),
    profile9$qc$condition_number, profile9$qc$max_abs_celltype_cor)
rm(sc); invisible(gc())   # free the 1.9 G reference before deconvolution

## ---- background + two deconvolutions ------------------------------------
say("[2/6] deriving Q3-scaled negative-probe background...")
bg <- derive_geomx_background(geomx)
say("  background %d genes x %d ROIs | range [%.3f, %.3f]",
    nrow(bg), ncol(bg), min(bg), max(bg))

say("[3/6] spatialdecon x2 (6-level then 9-level)...")
sd6 <- run_spatialdecon(geomx, profile6$profile, bg)
sd9 <- run_spatialdecon(geomx, profile9$profile, bg)
say("  beta6 %d x %d  resid_rmse=%.3f | beta9 %d x %d  resid_rmse=%.3f | aligned genes=%d",
    nrow(sd6$beta), ncol(sd6$beta), sd6$fit_qc$resid_rmse,
    nrow(sd9$beta), ncol(sd9$beta), sd9$fit_qc$resid_rmse,
    sd6$fit_qc$n_genes_aligned)

## ---- two-stage assembly -------------------------------------------------
say("[4/6] two-stage assembly (anchor substate fractions to stage-1 total)...")
two <- combine_two_stage(sd6$beta, sd9$beta,
                         micro_label = "Microglia", substate_prefix = "Microglia_")
say("  stage1 %d types x %d ROI | stage2 %d substates x %d ROI",
    nrow(two$stage1), ncol(two$stage1), nrow(two$stage2), ncol(two$stage2))
say("  unresolved ROIs (0 microglia in 9-level fit) = %d | consistency Spearman = %.4f",
    two$n_unresolved, two$consistency_spearman)

## ---- factorial contrasts on each layer ----------------------------------
say("[5/6] 5-contrast 2x2 factorial on log-abundance (both layers)...")
ct_stage1 <- fit_abundance_contrasts(two$stage1, roi_meta, offset = OFFSET,
                                     genotype_col = "genotype", slide_col = "slide",
                                     padj_cut = PADJ_CUT) |>
  mutate(layer = "stage1_broad", .before = 1)
ct_stage2 <- fit_abundance_contrasts(two$stage2, roi_meta, offset = OFFSET,
                                     genotype_col = "genotype", slide_col = "slide",
                                     padj_cut = PADJ_CUT) |>
  mutate(layer = "stage2_substate", .before = 1)
contrasts_tbl <- bind_rows(ct_stage1, ct_stage2)
say("  sig (adj.P<%.2f) rows: stage1=%d/%d | stage2=%d/%d | at interaction=%d",
    PADJ_CUT, sum(ct_stage1$sig), nrow(ct_stage1),
    sum(ct_stage2$sig), nrow(ct_stage2),
    sum(contrasts_tbl$sig & contrasts_tbl$contrast == "interaction"))

## ---- spatial autocorrelation (raw + genotype-residualised) --------------
say("[6/6] per-slide Moran's I (%d perms; raw + genotype-residualised)...", N_PERM)
# stage1 is always clean; include stage2 substates only when fully resolved
# (an NA substate column would NaN the per-slide Moran's I).
autocorr_input <- if (two$n_unresolved == 0) rbind(two$stage1, two$stage2) else two$stage1
coords <- as.matrix(roi_meta[, c("x", "y")])
sa_raw <- spatial_autocorrelation(autocorr_input, coords, roi_meta$slide,
                                  covar = NULL, n_perm = N_PERM, seed = SEED)
sa_res <- spatial_autocorrelation(autocorr_input, coords, roi_meta$slide,
                                  covar = as.character(roi_meta$genotype),
                                  n_perm = N_PERM, seed = SEED)
spatial_autocorr <- bind_rows(sa_raw, sa_res)
say("  autocorr rows: %d cell types x %d slides x {raw,resid} = %d",
    nrow(autocorr_input), length(unique(roi_meta$slide)), nrow(spatial_autocorr))
say("  sig positive autocorr (p<0.05): raw=%d | residualised=%d",
    sum(sa_raw$p_perm < 0.05, na.rm = TRUE), sum(sa_res$p_perm < 0.05, na.rm = TRUE))

## ---- abundance-by-genotype summary (all cell types, even-handed) --------
summarise_by_genotype <- function(abund, prop, layer) {
  rows <- list()
  for (k in rownames(abund)) for (g in genotype_levels) {
    ix <- which(roi_meta$genotype == g)
    a <- as.numeric(abund[k, ix]); p <- as.numeric(prop[k, ix])
    rows[[length(rows) + 1L]] <- data.frame(
      layer = layer, cell_type = k, genotype = g,
      n_roi = sum(!is.na(a)),
      mean_abundance = mean(a, na.rm = TRUE), sd_abundance = sd(a, na.rm = TRUE),
      mean_prop = mean(p, na.rm = TRUE), sd_prop = sd(p, na.rm = TRUE),
      stringsAsFactors = FALSE)
  }
  bind_rows(rows)
}
by_genotype <- bind_rows(
  summarise_by_genotype(two$stage1, sd6$prop, "stage1_broad"),
  summarise_by_genotype(two$stage2, two$substate_fractions, "stage2_substate"))

## ---- write results ------------------------------------------------------
say("[build_spatial_deconvolution] writing results...")
wtsv(by_genotype,      "spatial_decon_abundance_by_genotype.tsv")
wtsv(contrasts_tbl,    "spatial_decon_contrasts.tsv")
wtsv(spatial_autocorr, "spatial_decon_spatial_autocorr.tsv")

## ---- assemble cache for rmd/21 ------------------------------------------
spatial_decon <- list(
  profile6 = profile6,                 # list(profile, qc[kappa,max|cor|,cor mat,...])
  profile9 = profile9,
  beta6    = sd6,                       # list(beta, prop, fit_qc, sigmas)
  beta9    = sd9,
  stage1   = two$stage1,               # 6-level abundance matrix
  stage2   = two$stage2,               # substate anchored abundance matrix
  two_stage = two,                     # microglia_total, fractions, n_unresolved, consistency
  fit_qc = list(
    stage6 = sd6$fit_qc, stage9 = sd9$fit_qc,
    n_unresolved = two$n_unresolved,
    consistency_spearman = two$consistency_spearman),
  abundance_contrasts   = contrasts_tbl,
  abundance_by_genotype = by_genotype,
  spatial_autocorr      = spatial_autocorr,
  meta = roi_meta,
  params = list(
    granularity = "two_stage", cap = CAP, minCellNum = MINCELLNUM,
    minGenes = MINGENES, offset = OFFSET, padj_cut = PADJ_CUT,
    n_perm = N_PERM, seed = SEED, design = "genotype + slide",
    background = "Q3-scaled NegGeoMean (NegGeoMean / q_norm_qFactors)",
    cell_counts_used = FALSE,            # nuclei -1 sentinel in 42/91 ROIs
    n_common_genes = profile6$qc$n_common_genes,
    normalize_profile = TRUE, weighted = NA))
saveRDS(spatial_decon, out_cache, compress = "xz")
Sys.chmod(out_cache, "0644")
say("[build_spatial_deconvolution] wrote cache: %s", out_cache)

## ---- validation summary -------------------------------------------------
stopifnot(
  ncol(two$stage1) == 91L, ncol(two$stage2) == 91L,
  nrow(two$stage1) == 6L, nrow(two$stage2) == 4L,
  all(is.finite(two$stage1)),
  profile6$qc$n_common_genes == 18512L)
say("[validation] stage1 6x91 + stage2 4x91 finite; common genes=18512 OK")
say("[build_spatial_deconvolution] done.")

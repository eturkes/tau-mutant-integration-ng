#!/usr/bin/env Rscript
# Build the decoupleR + OmniPath kinase activity cache for Phase F of
# the mechanism layer (see storage/notes/mechanism_layer_plan.md,
# sessions F2..F5). Heavy compute lives here, OUTSIDE the knit, per
# the project convention in scripts/build_*.R; the Rmd consumer in
# session F3+ `readRDS`s the resulting cache directly.
#
# Idempotent: skips the cache write if the output already exists
# unless `--overwrite` is passed. Re-running with no args after a
# successful build is a no-op.
#
# Inputs (rds, under storage/cache/):
#   de_phospho.rds            raw phospho limma DE (14,715 sites,
#                             5 contrasts)
#   de_phospho_corrected.rds  batch-corrected phospho DE (12,821
#                             sites, 5 contrasts)
# Both caches are list(meta, design, contrasts, fit, log_matrix)
# with `fit$fit` = MArrayLM and `fit$top` = named list of
# per-contrast limma top-tables across the canonical five-contrast
# set (nlgf_in_maptki, nlgf_in_p301s, tau_alone, tau_in_nlgf,
# interaction).
#
# Prior: OmniPath kinase-substrate network, fetched via
# `R/kinase_inference.R::build_omnipath_ksn_mouse()` which wraps
# `OmnipathR::enzyme_substrate(organism = 9606)` (the
# `decoupleR::get_ksn_omnipath()` shortcut is broken in the current
# decoupleR/OmnipathR stack -- see kinase_inference.R header for
# details) and mouse-maps via `nichenetr::convert_human_to_mouse_symbols`.
# ~38k mouse edges, ~1.6k unique kinases, ~16k substrate ids;
# columns source / target (`SYMBOL_resPos`) / mor (+1 phosphorylation,
# -1 dephosphorylation).
#
# Output (under storage/cache/):
#   kinase_activity_decoupler.rds
#   Nested list keyed by `cache -> contrast -> tibble(statistic,
#   source, score, p_value)`. Caches (2):
#     phospho_raw
#     phospho_corrected
#   Contrasts (5 per cache):
#     nlgf_in_maptki, nlgf_in_p301s, interaction, tau_alone, tau_in_nlgf
#   Statistics per row: one of (`ulm`, `wsum`, `norm_wsum`,
#   `corr_wsum`, `consensus`). Downstream rankings (F3..F5) should
#   apply the E3 convention: rank on `consensus` for direction /
#   magnitude, gate significance on `ulm`'s p-value (the linear
#   model F-test), BH-adjust within each (cache, contrast) cell.
#
# Runtime: ~6s for the KSN build + mouse map; ~20-25s for each
# (cache x 5 contrasts) decoupleR call; total ~50-60s wall on the
# project's Docker container. The OmnipathR cache for the underlying
# enzyme-substrate table is local-disk-resident after the first
# invocation, so re-runs of this script after the initial pull
# do not re-hit the OmniPath HTTP API.
#
# Phase F is single-modality (phospho only) because phospho is the
# only assay in the project that exposes per-site biology -- bulk
# proteomics and snRNAseq cannot disambiguate phosphorylation events.
# So this script has no analogue of the TF script's
# `modality_top_lists` structure; it iterates over a tighter
# `cache_top_lists` containing two entries (the raw and corrected
# phospho caches). The "cross-cache" axis in F3+ plays the role of
# E3's "cross-modality" axis; the asymmetry is documented in the
# F3 interpretive notes per the plan.

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

args      <- commandArgs(trailingOnly = TRUE)
overwrite <- "--overwrite" %in% args

setwd("/home/rstudio/tau-mutant-integration-ng")
source("R/helpers.R")  # pulls in R/kinase_inference.R via the project loader

cache_dir <- "storage/cache"
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
out_path  <- file.path(cache_dir, "kinase_activity_decoupler.rds")

if (file.exists(out_path) && !overwrite) {
  cat(sprintf("[build_kinase_activity_decoupler] cache exists, skipping: %s\n",
              out_path))
  cat("Pass --overwrite to rebuild.\n")
  quit(save = "no", status = 0)
}

cat("[build_kinase_activity_decoupler] loading phospho DE caches...\n")
de_phospho_raw <- readRDS(file.path(cache_dir, "de_phospho.rds"))
de_phospho_cor <- readRDS(file.path(cache_dir, "de_phospho_corrected.rds"))

cat("[build_kinase_activity_decoupler] building mouse OmniPath KSN...\n")
t0 <- proc.time()
ksn_mouse <- build_omnipath_ksn_mouse()
t1 <- proc.time()
cat(sprintf("  KSN: %d edges, %d kinases, %d substrate ids (built in %.2fs)\n",
            nrow(ksn_mouse),
            length(unique(ksn_mouse$source)),
            length(unique(ksn_mouse$target)),
            (t1 - t0)[3]))

# Per-cache top-tables.
cache_top_lists <- list(
  phospho_raw       = de_phospho_raw$fit$top,
  phospho_corrected = de_phospho_cor$fit$top
)

cat("[build_kinase_activity_decoupler] running decoupleR per cache...\n")
t_start <- proc.time()
kinase_activity_decoupler <- list()
for (cn in names(cache_top_lists)) {
  top_list <- cache_top_lists[[cn]]
  stat_mat <- extract_phospho_stat_matrix(top_list)
  cat(sprintf("  cache %-18s  %5d sites x %d contrasts ... ",
              cn, nrow(stat_mat), ncol(stat_mat)))
  t0 <- proc.time()
  dec_tbl <- run_decoupler_per_cache(
    stat_mat   = stat_mat,
    network    = ksn_mouse,
    statistics = c("ulm", "wsum"),
    minsize    = 5L,
    consensus  = TRUE
  )
  t1 <- proc.time()
  kinase_activity_decoupler[[cn]] <- split_kinase_decoupler_by_contrast(dec_tbl)
  cat(sprintf("done in %.2fs (kinases per contrast: %d, statistics per kinase: %d)\n",
              (t1 - t0)[3],
              length(unique(dec_tbl$source)),
              length(unique(dec_tbl$statistic))))
}
t_end <- proc.time()
cat(sprintf("[build_kinase_activity_decoupler] all caches done in %.2fs\n",
            (t_end - t_start)[3]))

cat(sprintf("[build_kinase_activity_decoupler] writing cache: %s\n", out_path))
saveRDS(kinase_activity_decoupler, out_path, compress = "xz")
Sys.chmod(out_path, mode = "0644")

cat("[build_kinase_activity_decoupler] done.\n")

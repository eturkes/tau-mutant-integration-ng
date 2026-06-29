#!/usr/bin/env Rscript
# Build a per-microglia-substate TF activity cache for Phase I of the
# NF-kB attenuation test (see storage/notes/nfkb_attenuation_plan.md,
# session I2). Heavy compute lives here, OUTSIDE the knit, per the
# project convention in scripts/build_*.R; the Rmd consumer added by
# session I3 `readRDS`s the resulting cache directly.
#
# Distinction from `scripts/build_tf_activity_decoupler.R`: that script
# emits a `modality -> contrast -> tibble` cache covering whole-microglia
# nebula + per-substate nebula (built from the UNFILTERED
# `de_snrnaseq_nebula_per_state.rds`) + GeoMx limma + proteomics limma.
# This script emits a SEPARATE `state -> contrast -> tibble` cache built
# from the 1%-prevalence-FILTERED per-state NEBULA cache
# (`de_snrnaseq_nebula_per_state_1pct.rds`). The 1pct cache has ~2.3x
# more genes per state (homeostatic 11322 vs 5018; DAM 11627 vs 5144;
# IFN 11972 vs 5560; proliferative 11544 vs 4913) which materially
# increases per-state TF target-set coverage and TF inference power.
# The plan locks this 1pct cache as the input for the Phase I per-state
# NF-kB attenuation test because the test asks a per-state-specific
# question (which substate carries the NF-kB suppression?) and benefits
# from the maximum per-state gene coverage.
#
# Idempotent: skips the cache write if the output already exists unless
# `--overwrite` is passed. Re-running with no args after a successful
# build is a no-op.
#
# Inputs (rds, under storage/cache/):
#   de_snrnaseq_nebula_per_state_1pct.rds  per-substate nebula fits
#                                          (homeostatic / DAM / IFN /
#                                          proliferative; same five
#                                          contrasts each; ~11k-12k
#                                          genes per state after the
#                                          1% prevalence filter).
#
# Prior: CollecTRI mouse (decoupleR::get_collectri(organism = "mouse")).
# Cached locally by OmnipathR; ~40k edges, ~1.1k unique TFs, ~6.1k
# targets, columns source / target / mor. Same prior used by the
# whole-microglia / GeoMx / proteomics TF inference layer in
# build_tf_activity_decoupler.R so cross-cache comparisons remain
# meaningful (the prior network is identical; only the per-state input
# matrix differs).
#
# Output (under storage/cache/):
#   per_state_tf_activity.rds
#   Nested list keyed by `state -> contrast -> tibble(statistic, source,
#   score, p_value)`. States (4):
#     homeostatic, DAM, IFN, proliferative
#   Contrasts (5 per state):
#     nlgf_in_maptki, nlgf_in_p301s, interaction, tau_alone, tau_in_nlgf
#   Statistics per row: one of (`ulm`, `wsum`, `norm_wsum`, `corr_wsum`,
#   `consensus`). Downstream rankings should prefer the `consensus`
#   score per decoupleR's recommendation; significance is read from the
#   `ulm` p_value column, BH-adjusted within (state x contrast), with
#   the FDR<0.10 cutoff matching the project's locked TF-inference
#   convention (R/tf_inference.R module header).
#
# Output (under storage/results/):
#   per_state_tf_activity.tsv
#   Long-form table for downstream consumption by the I3 chapter and
#   any external sharing. Columns: state, contrast, tf, score, p_value,
#   padj, sig_fdr10_ulm. Only the `ulm` statistic is written to the TSV
#   because (a) the I3 chapter's NF-kB family table is sourced from
#   `ulm` (it is the canonical frequentist univariate test per the
#   project's locked TF-inference convention), and (b) keeping the TSV
#   single-statistic prevents downstream confusion about which row's
#   score / p / padj triple is the headline. The full multi-statistic
#   cache remains the source of truth for the consensus + wsum scores
#   that the I3 chapter may also need.
#
# Runtime: ~30s per state x 5 contrasts call * 4 states = ~2 min wall,
# plus CollecTRI prior load (~5s on first call, instant from OmnipathR
# cache thereafter). Total well under 3 min on the project's Docker
# container.

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

args      <- commandArgs(trailingOnly = TRUE)
overwrite <- "--overwrite" %in% args

setwd("/home/rstudio/tau-mutant-integration-ng")
source("R/helpers.R")  # pulls in R/tf_inference.R via the project loader

cache_dir   <- "storage/cache"
results_dir <- "storage/results"
dir.create(cache_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

cache_path <- file.path(cache_dir,   "per_state_tf_activity.rds")
tsv_path   <- file.path(results_dir, "per_state_tf_activity.tsv")

if (file.exists(cache_path) && file.exists(tsv_path) && !overwrite) {
  cat(sprintf("[build_per_state_tf_activity] outputs exist, skipping:\n"))
  cat(sprintf("  cache: %s\n", cache_path))
  cat(sprintf("  tsv  : %s\n", tsv_path))
  cat("Pass --overwrite to rebuild.\n")
  quit(save = "no", status = 0)
}

cat("[build_per_state_tf_activity] loading 1pct-filtered per-state NEBULA cache...\n")
neb_per_state <- readRDS(
  file.path(cache_dir, "de_snrnaseq_nebula_per_state_1pct.rds")
)
states <- names(neb_per_state)
cat(sprintf("  states (%d): %s\n", length(states),
            paste(states, collapse = ", ")))

cat("[build_per_state_tf_activity] loading CollecTRI mouse prior...\n")
collectri_mouse <- decoupleR::get_collectri(organism        = "mouse",
                                            split_complexes = FALSE)
cat(sprintf("  CollecTRI mouse: %d edges, %d unique TFs, %d targets\n",
            nrow(collectri_mouse),
            length(unique(collectri_mouse$source)),
            length(unique(collectri_mouse$target))))

cat("[build_per_state_tf_activity] running decoupleR per state...\n")
t_start <- proc.time()
per_state_tf_activity <- list()
for (st in states) {
  top_list <- neb_per_state[[st]]$top
  stat_mat <- extract_de_stat_matrix(top_list, stat_col = "t",
                                     id_col = "symbol")
  cat(sprintf("  state %-15s %5d symbols x %d contrasts ... ",
              st, nrow(stat_mat), ncol(stat_mat)))
  t0 <- proc.time()
  dec_tbl <- run_decoupler_per_modality(
    stat_mat   = stat_mat,
    network    = collectri_mouse,
    statistics = c("ulm", "wsum"),
    minsize    = 5L,
    consensus  = TRUE
  )
  t1 <- proc.time()
  per_state_tf_activity[[st]] <- split_decoupler_by_contrast(dec_tbl)
  cat(sprintf("done in %.2fs (TFs per contrast: %d, statistics per TF: %d)\n",
              (t1 - t0)[3],
              length(unique(dec_tbl$source)),
              length(unique(dec_tbl$statistic))))
}
t_end <- proc.time()
cat(sprintf("[build_per_state_tf_activity] all states done in %.2fs\n",
            (t_end - t_start)[3]))

cat(sprintf("[build_per_state_tf_activity] writing cache: %s\n",
            cache_path))
saveRDS(per_state_tf_activity, cache_path, compress = "xz")
Sys.chmod(cache_path, mode = "0644")

# Emit the long-form TSV. The TSV restricts to the `ulm` statistic
# because it is the canonical frequentist sig test (per R/tf_inference.R)
# and keeps the file readable for downstream consumers. The full
# multi-statistic cache remains the source of truth for the consensus +
# wsum scores; downstream code that needs them reads the cache, not the
# TSV.
cat("[build_per_state_tf_activity] building long-form TSV (ulm only)...\n")
long_rows <- list()
for (st in names(per_state_tf_activity)) {
  for (cn in names(per_state_tf_activity[[st]])) {
    tbl <- per_state_tf_activity[[st]][[cn]]
    ulm_rows <- tbl[tbl$statistic == "ulm", , drop = FALSE]
    if (nrow(ulm_rows) == 0L) next
    padj <- stats::p.adjust(ulm_rows$p_value, method = "BH")
    long_rows[[paste(st, cn, sep = "::")]] <- tibble::tibble(
      state          = st,
      contrast       = cn,
      tf             = ulm_rows$source,
      score          = ulm_rows$score,
      p_value        = ulm_rows$p_value,
      padj           = padj,
      sig_fdr10_ulm  = !is.na(padj) & padj < 0.10
    )
  }
}
long_tbl <- dplyr::bind_rows(long_rows)
long_tbl <- long_tbl[order(long_tbl$state, long_tbl$contrast,
                           long_tbl$padj), , drop = FALSE]
cat(sprintf("  TSV: %d rows (4 states x 5 contrasts x ~%d TFs per slice)\n",
            nrow(long_tbl),
            round(nrow(long_tbl) / (length(states) *
              length(names(per_state_tf_activity[[1]]))))))
n_sig <- sum(long_tbl$sig_fdr10_ulm, na.rm = TRUE)
cat(sprintf("  %d rows pass sig_fdr10_ulm (FDR<0.10 on ulm p_value)\n",
            n_sig))

cat(sprintf("[build_per_state_tf_activity] writing TSV: %s\n", tsv_path))
write.table(long_tbl, tsv_path, sep = "\t", quote = FALSE,
            row.names = FALSE, na = "")
Sys.chmod(tsv_path, mode = "0644")

cat("[build_per_state_tf_activity] done.\n")

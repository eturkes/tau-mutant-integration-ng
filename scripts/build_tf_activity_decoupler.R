#!/usr/bin/env Rscript
# Build the decoupleR + CollecTRI TF activity cache for Phase E of the
# mechanism layer (see storage/notes/mechanism_layer_plan.md, sessions
# E2..E5). Heavy compute lives here, OUTSIDE the knit, per the project
# convention in scripts/build_*.R; the Rmd consumer in session E3+
# `readRDS`s the resulting cache directly.
#
# Idempotent: skips the cache write if the output already exists unless
# `--overwrite` is passed. Re-running with no args after a successful
# build is a no-op.
#
# Inputs (rds, under storage/cache/):
#   de_snrnaseq_nebula.rds                whole-microglia nebula fits
#   de_snrnaseq_nebula_per_state.rds      per-substate nebula fits
#                                         (homeostatic / DAM / IFN /
#                                         proliferative; same five
#                                         contrasts each)
#   de_geomx.rds                          limma GeoMx fit (full pipeline)
#   de_proteomics.rds                     limma proteomics fit
#
# Phospho is deliberately deferred to Phase F because phospho-site stats
# need protein-level collapse + signed re-mapping before they can serve
# as TF-target evidence -- a kinase-side concern, not a TF-side one.
#
# Prior: CollecTRI mouse (decoupleR::get_collectri(organism = "mouse")).
# Cached locally by OmnipathR; ~40k edges, ~1.1k unique TFs, ~6.1k
# targets, columns source / target / mor.
#
# Output (under storage/cache/):
#   tf_activity_decoupler.rds
#   Nested list keyed by `modality -> contrast -> tibble(statistic,
#   source, score, p_value)`. Modalities (7):
#     snrnaseq                          whole-microglia nebula
#     snrnaseq_homeostatic              per-substate nebula
#     snrnaseq_DAM                      per-substate nebula
#     snrnaseq_IFN                      per-substate nebula
#     snrnaseq_proliferative            per-substate nebula
#     geomx                             limma
#     proteomics                        limma
#   Contrasts (5 per modality):
#     nlgf_in_maptki, nlgf_in_p301s, interaction, tau_alone, tau_in_nlgf
#   Statistics per row: one of (`ulm`, `wsum`, `norm_wsum`, `corr_wsum`,
#   `consensus`). Downstream rankings should prefer the `consensus`
#   score per decoupleR's recommendation, with the per-method scores
#   retained for sensitivity checks.
#
# Runtime: ~5-6s per (modality x 5 contrasts) call * 7 modalities =
# ~40-50s wall, plus CollecTRI prior load (~5s on first call, instant
# from OmnipathR cache thereafter). Total well under 2 min on the
# project's Docker container.

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)
overwrite <- "--overwrite" %in% args
# --split-complexes: re-derive the identical cache but with CollecTRI
# complexes split into their single-protein members (NF-kB ->
# Nfkb1/Rela, etc.), writing tf_activity_decoupler_split.rds. Required
# by the causal-network layer (J1 dec.2): CARNIVAL nodes must be single
# proteins in mouse-symbol space, but the default cache keeps complexes
# as UniProt accessions. Everything else (priors-aside) is byte-identical
# to the default build, so the two variants stay directly comparable.
split_complexes <- "--split-complexes" %in% args

setwd("/home/rstudio/tau-mutant-integration-ng")
source("R/helpers.R")  # pulls in R/tf_inference.R via the project loader

cache_dir <- "storage/cache"
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
out_path  <- file.path(cache_dir,
                       if (split_complexes) "tf_activity_decoupler_split.rds"
                       else                 "tf_activity_decoupler.rds")

if (file.exists(out_path) && !overwrite) {
  cat(sprintf("[build_tf_activity_decoupler] cache exists, skipping: %s\n",
              out_path))
  cat("Pass --overwrite to rebuild.\n")
  quit(save = "no", status = 0)
}

cat("[build_tf_activity_decoupler] loading DE caches...\n")
de_snrnaseq_nebula        <- readRDS(file.path(cache_dir,
                                               "de_snrnaseq_nebula.rds"))
de_snrnaseq_nebula_state  <- readRDS(file.path(cache_dir,
                                               "de_snrnaseq_nebula_per_state.rds"))
de_geomx                  <- readRDS(file.path(cache_dir, "de_geomx.rds"))
de_proteomics             <- readRDS(file.path(cache_dir, "de_proteomics.rds"))

cat(sprintf("[build_tf_activity_decoupler] loading CollecTRI mouse prior (split_complexes=%s)...\n",
            split_complexes))
collectri_mouse <- decoupleR::get_collectri(organism        = "mouse",
                                            split_complexes = split_complexes)
cat(sprintf("  CollecTRI mouse: %d edges, %d unique TFs, %d targets\n",
            nrow(collectri_mouse),
            length(unique(collectri_mouse$source)),
            length(unique(collectri_mouse$target))))

# Top-table extractor per modality. The three supported shapes:
#   nebula whole      -> cache$top
#   nebula substate   -> cache[[substate]]$top
#   limma             -> cache$fit$top
# We define a per-modality list-of-top-tables here, then run the same
# decoupleR pipeline across all of them.
modality_top_lists <- list(
  snrnaseq               = de_snrnaseq_nebula$top,
  snrnaseq_homeostatic   = de_snrnaseq_nebula_state$homeostatic$top,
  snrnaseq_DAM           = de_snrnaseq_nebula_state$DAM$top,
  snrnaseq_IFN           = de_snrnaseq_nebula_state$IFN$top,
  snrnaseq_proliferative = de_snrnaseq_nebula_state$proliferative$top,
  geomx                  = de_geomx$fit$top,
  proteomics             = de_proteomics$fit$top
)

cat("[build_tf_activity_decoupler] running decoupleR per modality...\n")
t_start <- proc.time()
tf_activity_decoupler <- list()
for (mn in names(modality_top_lists)) {
  top_list <- modality_top_lists[[mn]]
  stat_mat <- extract_de_stat_matrix(top_list)
  cat(sprintf("  modality %-22s  %5d symbols x %d contrasts ... ",
              mn, nrow(stat_mat), ncol(stat_mat)))
  t0 <- proc.time()
  dec_tbl <- run_decoupler_per_modality(
    stat_mat   = stat_mat,
    network    = collectri_mouse,
    statistics = c("ulm", "wsum"),
    minsize    = 5L,
    consensus  = TRUE
  )
  t1 <- proc.time()
  tf_activity_decoupler[[mn]] <- split_decoupler_by_contrast(dec_tbl)
  cat(sprintf("done in %.2fs (TFs per contrast: %d, statistics per TF: %d)\n",
              (t1 - t0)[3],
              length(unique(dec_tbl$source)),
              length(unique(dec_tbl$statistic))))
}
t_end <- proc.time()
cat(sprintf("[build_tf_activity_decoupler] all modalities done in %.2fs\n",
            (t_end - t_start)[3]))

cat(sprintf("[build_tf_activity_decoupler] writing cache: %s\n", out_path))
saveRDS(tf_activity_decoupler, out_path, compress = "xz")
Sys.chmod(out_path, mode = "0644")

cat("[build_tf_activity_decoupler] done.\n")

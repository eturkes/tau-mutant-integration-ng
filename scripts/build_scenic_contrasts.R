#!/usr/bin/env Rscript
# K4 driver: connect the data-driven SCENIC consensus regulons to the
# project's 5-contrast framework. Thin I/O wrapper over R/scenic.R (pure
# functions); mirrors the build_tf_activity_decoupler.R pattern (heavy
# work outside the knit, results consumed by rmd/20). Idempotent unless
# --overwrite.
#
# Three readouts + diagnostics, all reusing the §14 decoupleR machinery
# (controlled CollecTRI -> SCENIC network swap) and the locked 2x2
# factorial:
#   1. head-to-head  SCENIC vs CollecTRI TF activity on the SAME NEBULA z.
#   2. native AUCell pseudobulk -> 2x2 factorial (logit-transformed means).
#   3. per-substate SCENIC activity at `interaction` (mirrors §18).
#   + target-set Jaccard overlap; + regulon recovery ladder.
#
# Outputs (storage/results/*.tsv + storage/cache/scenic_summary.rds):
#   scenic_regulon_recovery.tsv         recovery ladder (51 regulons + §14/§18 TFs)
#   scenic_vs_collectri_headtohead.tsv  controlled network-swap comparison
#   scenic_aucell_contrasts.tsv         AUCell factorial (all 5 contrasts)
#   scenic_per_substate_interaction.tsv per-substate activity at interaction
#   scenic_target_overlap.tsv           data-driven vs prior target Jaccard

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tibble)
})

setwd("/home/rstudio/tau-mutant-integration-ng")
source("R/helpers.R")

args      <- commandArgs(trailingOnly = TRUE)
overwrite <- "--overwrite" %in% args
cc  <- "storage/cache"
sc  <- file.path(cc, "scenic")
res <- "storage/results"
dir.create(res, recursive = TRUE, showWarnings = FALSE)
out_cache <- file.path(cc, "scenic_summary.rds")
PADJ_CUT  <- 0.10

if (file.exists(out_cache) && !overwrite) {
  cat(sprintf("[build_scenic_contrasts] cache exists, skipping: %s\n", out_cache))
  cat("Pass --overwrite to rebuild.\n"); quit(save = "no", status = 0)
}

say <- function(...) cat(sprintf(...), "\n")
wtsv <- function(df, name) {
  p <- file.path(res, name)
  write_tsv_safe(as.data.frame(df), p); Sys.chmod(p, "0644")
  say("  wrote %-42s %d rows", name, nrow(df))
}

## ---- load inputs --------------------------------------------------------
say("[build_scenic_contrasts] loading inputs...")
reg     <- as.data.frame(fread(file.path(sc, "scenic_regulons.tsv")))
census  <- as.data.frame(fread(file.path(sc, "scenic_recovery_census.tsv")))
neb     <- readRDS(file.path(cc, "de_snrnaseq_nebula.rds"))
tfa     <- readRDS(file.path(cc, "tf_activity_decoupler.rds"))
ps      <- readRDS(file.path(cc, "de_snrnaseq_nebula_per_state_1pct.rds"))
aucell  <- as.data.frame(fread(file.path(sc, "aucell.csv.gz")))
colattr <- as.data.frame(fread(file.path(sc, "microglia_colattrs.tsv")))
mg_genes    <- readLines(file.path(sc, "microglia_genes.txt"))
all_tfs     <- readLines("storage/data/cistarget/allTFs_mm.txt")
candidate_tfs <- intersect(mg_genes, all_tfs)
collectri_net <- decoupleR::get_collectri(organism = "mouse",
                                          split_complexes = FALSE)
say("  regulons=%d edges/%d TFs | candidate TFs=%d | CollecTRI=%d edges",
    nrow(reg), length(unique(reg$TF)), length(candidate_tfs), nrow(collectri_net))

scenic_net <- build_scenic_network(reg, weighted = FALSE)   # sign-only mor (locked)
comparison <- unique(census[, c("TF", "axis")])

## ---- readout 1: controlled head-to-head ---------------------------------
say("[1/5] head-to-head: SCENIC vs CollecTRI on NEBULA z...")
scenic_by_c   <- run_scenic_decoupler(neb$top, scenic_net,
                                      statistics = c("ulm", "wsum"), minsize = 5L)
scenic_long   <- decoupler_activity_long(scenic_by_c, "scenic", PADJ_CUT)
collectri_long <- decoupler_activity_long(tfa$snrnaseq, "collectri", PADJ_CUT)
headtohead    <- build_scenic_headtohead(scenic_long, collectri_long,
                                         scenic_net, comparison, PADJ_CUT)
wtsv(headtohead, "scenic_vs_collectri_headtohead.tsv")
say("  SCENIC regulons sig (padj<%.2f) in >=1 contrast: %d / %d",
    PADJ_CUT, length(unique(scenic_long$source[scenic_long$sig])),
    length(unique(scenic_long$source)))

## ---- readout 2: native AUCell pseudobulk factorial ----------------------
say("[2/5] AUCell pseudobulk -> 2x2 factorial...")
pb <- aucell_to_pseudobulk(aucell, colattr, id_col = "genotype_batch")
say("  pseudobulk %d regulons x %d ids; mean range [%.4f, %.4f]",
    nrow(pb$mat), ncol(pb$mat), min(pb$mat), max(pb$mat))
aucell_contrasts     <- fit_aucell_contrasts(pb$mat, pb$meta, transform = "logit",
                                             padj_cut = PADJ_CUT)
aucell_contrasts_raw <- fit_aucell_contrasts(pb$mat, pb$meta, transform = "none",
                                             padj_cut = PADJ_CUT)
# logit-vs-raw sensitivity: do the two transforms agree on contrast direction?
mrg <- merge(aucell_contrasts[, c("regulon","contrast","t","sig")],
             aucell_contrasts_raw[, c("regulon","contrast","t","sig")],
             by = c("regulon","contrast"), suffixes = c("_logit","_raw"))
sens_rho  <- suppressWarnings(cor(mrg$t_logit, mrg$t_raw, method = "spearman"))
sens_sign <- mean(sign(mrg$t_logit) == sign(mrg$t_raw))
say("  logit vs raw: Spearman(t)=%.3f, sign-agreement=%.1f%%, sig-jaccard=%.2f",
    sens_rho, 100*sens_sign,
    { a <- mrg$sig_logit; b <- mrg$sig_raw
      if (sum(a|b)==0) 1 else sum(a&b)/sum(a|b) })
aucell_contrasts <- aucell_contrasts |>
  mutate(TF = sub("\\(.+$", "", regulon)) |>
  relocate(TF, .after = regulon)
wtsv(aucell_contrasts, "scenic_aucell_contrasts.tsv")
say("  AUCell regulons sig (adj.P<%.2f) in >=1 contrast: %d ; at interaction: %d",
    PADJ_CUT, length(unique(aucell_contrasts$regulon[aucell_contrasts$sig])),
    sum(aucell_contrasts$sig & aucell_contrasts$contrast == "interaction"))

## ---- readout 3: per-substate at interaction -----------------------------
say("[3/5] per-substate SCENIC activity at interaction...")
substate_int <- scenic_substate_activity(ps, scenic_net, contrast = "interaction",
                                         statistics = c("ulm", "wsum"),
                                         minsize = 5L, padj_cut = PADJ_CUT)
substate_int <- substate_int |>
  mutate(is_focus = source %in% c("Rel", "Spi1"))
wtsv(substate_int, "scenic_per_substate_interaction.tsv")
say("  substates covered: %s", paste(sort(unique(substate_int$substate)), collapse=", "))

## ---- readout 4: target overlap ------------------------------------------
say("[4/5] target-set Jaccard (SCENIC vs CollecTRI)...")
overlap <- scenic_collectri_target_overlap(scenic_net, collectri_net)
wtsv(overlap, "scenic_target_overlap.tsv")
say("  shared TFs=%d ; median Jaccard=%.3f ; max=%.3f (%s)",
    nrow(overlap), median(overlap$jaccard, na.rm=TRUE),
    max(overlap$jaccard, na.rm=TRUE), overlap$TF[which.max(overlap$jaccard)])

## ---- recovery ladder ----------------------------------------------------
say("[5/5] recovery ladder...")
recovery <- build_recovery_table(census, scenic_net, candidate_tfs)
wtsv(recovery, "scenic_regulon_recovery.tsv")
say("  recovery_class tally:")
print(table(recovery$recovery_class))

## ---- assemble cache for rmd/20 ------------------------------------------
scenic_summary <- list(
  regulons          = reg,
  scenic_net        = scenic_net,
  headtohead        = headtohead,
  scenic_activity   = scenic_long,        # full 51-TF SCENIC long (all 5 contrasts)
  aucell_contrasts  = aucell_contrasts,
  aucell_pb         = pb,                  # mat + meta (for activity heatmaps)
  substate_int      = substate_int,
  target_overlap    = overlap,
  recovery          = recovery,
  census            = census,
  params = list(
    padj_cut = PADJ_CUT, weighted_mor = FALSE, transform = "logit",
    statistics = c("ulm","wsum"), minsize = 5L,
    recurrence_threshold = 8L, n_runs = 10L,
    n_consensus_regulons = length(unique(reg$TF)),
    aucell_logit_vs_raw = list(spearman_t = sens_rho,
                               sign_agreement = sens_sign),
    candidate_tfs_n = length(candidate_tfs),
    comparison_tfs = comparison)
)
saveRDS(scenic_summary, out_cache, compress = "xz")
Sys.chmod(out_cache, "0644")
say("[build_scenic_contrasts] wrote cache: %s", out_cache)
say("[build_scenic_contrasts] done.")

#!/usr/bin/env Rscript
# N3 driver (plan arc N): CELL-TYPE SPECIFICITY of the tau x amyloid
# interaction. Every prior arc (D..M) reads ONLY microglia; this is the
# project's first CROSS-CELL-TYPE readout. Re-runs the EXACT microglia DE
# treatment (rmd/02b: pseudobulk limma-voom + single-cell NEBULA) SYMMETRICALLY
# over 6 units (5 broad non-microglial types + microglia-as-whole, the 4
# substates pooled) to ask whether the locked interaction -- and the amyloid /
# tau main effects -- are microglia-SPECIFIC or shared. Thin I/O wrapper over
# the pure fns in R/celltype_specificity.R; mirrors build_trajectory.R /
# build_spatial_deconvolution.R -- all heavy compute OUT of the knit, rmd/23
# (N4) is display-only. Idempotent unless --overwrite.
#
# THE CONFOUND (the design is built around it): cell count varies ~25x across
# units (Neuronal ~163k vs the IFN-thin microglia), so pseudobulk precision ->
# power is asymmetric and a specificity claim from NATIVE counts alone is
# confounded. The load-bearing control is the MATCHED regime: downsample every
# unit x replicate to K = min-over-units of the per-replicate minimum (= 289,
# bound by pooled-Microglia; N2-locked) BEFORE refit, so the headline rests on
# power-equalised fits. NATIVE is reported too but flagged power-confounded.
#
# Method depth (N1 gate, user-chosen): the FULL stack -- TWO estimators
# (pseudobulk + NEBULA) x TWO power regimes (NATIVE + MATCHED), all 6 units,
# the locked 5 canonical contrasts. NATIVE NEBULA on the 163k Neuronal unit is
# the compute bottleneck -> each NATIVE NEBULA fit is cached to its own file
# (celltype_specificity_fits/) so a crashed run resumes; units run lightest-
# first so a Neuronal failure still leaves the other 5 cached.
#
# Cross-check (sanity gate, N3-locked): pooled-Microglia from seurat_full must
# reproduce the canonical microglia DE -- assert Spearman(logFC) > 0.90 on
# shared SYMBOLS per contrast for BOTH estimators (report all 5 rho; the
# interaction's ~0.93 is a genuine noisy-difference-of-differences property
# with 0 ambiguous symbols, so a uniform 0.90 gross-mismatch floor is used,
# not a 0.95 hard assert that would fail spuriously).
#
# Readouts (defined BEFORE inspection; all 6 units reported SYMMETRICALLY,
# alphabetical, no privileged ordering -- anti-anchoring guardrail):
#   R1 tally_native   significant-gene counts (FDR<0.05/<0.10, up/down) per
#                     unit x contrast x estimator, NATIVE
#   R2 tally_matched  R1 recomputed on the MATCHED fits (LOAD-BEARING)
#   R3 interaction_concordance  (a) genome-wide Spearman of interaction logFC
#                     microglia-vs-each-unit per estimator; (c) cross-estimator
#                     pseudobulk-vs-NEBULA interaction concordance per unit
#   R4 specificity_class  per-gene microglia_unique / shared /
#                     non_microglial_unique for interaction-sig genes (FDR<0.10)
#   R5 pathway_tally  GO BP fGSEA per unit at interaction -- is the microglial
#                     interaction ENRICHMENT cell-type-restricted? (+ the
#                     microglia-significant pathways' cross-unit NES in the cache)
#
# Outputs (storage/results/*.tsv + storage/cache/celltype_specificity.rds):
#   celltype_specificity_tally_native.tsv             R1
#   celltype_specificity_tally_matched.tsv            R2
#   celltype_specificity_interaction_concordance.tsv  R3
#   celltype_specificity_class.tsv                    R4
#   celltype_specificity_pathway_tally.tsv            R5
#   celltype_specificity_crosscheck.tsv               sanity-gate provenance
# rmd/23 (N4) writes celltype_specificity_verdict.tsv at knit time.
#
# Run: Rscript scripts/build_celltype_specificity.R [--overwrite]

suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(Matrix)
  library(dplyr); library(tibble); library(purrr)
})

setwd("/home/rstudio/tau-mutant-integration-ng")
source("R/helpers.R")

args      <- commandArgs(trailingOnly = TRUE)
overwrite <- "--overwrite" %in% args
cc  <- "storage/cache"
res <- "storage/results"
fit_dir <- file.path(cc, "celltype_specificity_fits")   # per-unit NEBULA resume cache (gitignored)
dir.create(res, recursive = TRUE, showWarnings = FALSE)
dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE)
out_cache <- file.path(cc, "celltype_specificity.rds")

## ---- embedded params (state before applying; N0/N1/N2/N3 locked) --------
SEED          <- 1L
MIN_COUNT     <- 10L      # pseudobulk min_count (rmd/02b)
MIN_CELL_FRAC <- 0.01     # NEBULA gene-prevalence filter (rmd/02b microglia headline)
PADJ_CUT      <- 0.10     # project-standard FDR for the class/pathway tallies
XCHECK_FLOOR  <- 0.90     # microglia cross-check hard floor (N3-locked; report all 5 rho)
HEADLINE_N    <- 20L      # top-N microglia NEBULA interaction symbols for the R3b set
K_EXPECTED    <- 289L     # N2-locked matched depth (bound by pooled-Microglia)
N_CELLS_EXP   <- 286285L
NCORE         <- max(1L, parallel::detectCores() - 2L)
# Canonical UNIT order: ALPHABETICAL, no privileged ordering (guardrail).
UNITS <- c("Astrocyte", "Microglia", "Neuronal", "Oligodendrocyte", "OPC", "Vascular")

if (file.exists(out_cache) && !overwrite) {
  cat(sprintf("[build_celltype_specificity] cache exists, skipping: %s\n", out_cache))
  cat("Pass --overwrite to rebuild.\n"); quit(save = "no", status = 0)
}

say <- function(...) cat(sprintf(...), "\n")
wtsv <- function(df, name) {
  p <- file.path(res, name)
  write_tsv_safe(as.data.frame(df), p); Sys.chmod(p, "0644")
  say("  wrote %-50s %d rows", name, nrow(df))
}
# Slim a wrapper fit to exactly what rmd/23 needs (the per-contrast top tables
# + cell/gene counts), dropping the heavy raw $fit / $design so the cache stays
# lean. The plot fns access $top[[contrast]] / $n_cells, so this is sufficient.
slim_fit <- function(fo) list(top = fo$top, n_cells = fo$n_cells, n_genes = fo$n_genes)
slim_list <- function(l) lapply(l, slim_fit)

## ---- 1. load + label the 6 units ----------------------------------------
say("[build_celltype_specificity] loading seurat_full_processed.rds (~1.9 G)...")
t_load <- Sys.time()
sc <- readRDS(file.path(cc, "seurat_full_processed.rds"))
stopifnot(ncol(sc) == N_CELLS_EXP)
sc$unit <- celltype_unit_labels(sc$cell_type)
stopifnot(setequal(unique(as.character(sc$unit)), UNITS))
unit_counts <- table(factor(as.character(sc$unit), levels = UNITS))
per_unit_min <- min_cells_per_replicate(sc)[UNITS]
K <- min(per_unit_min)
say("  loaded %d cells x %d genes in %.0f s | symbol-keyed (e.g. %s)",
    ncol(sc), nrow(sc), as.numeric(difftime(Sys.time(), t_load, units = "secs")),
    paste(head(rownames(sc), 2), collapse = ", "))
say("  unit cells: %s",
    paste(sprintf("%s=%d", names(unit_counts), as.integer(unit_counts)), collapse = ", "))
say("  per-unit min cells/replicate: %s",
    paste(sprintf("%s=%d", names(per_unit_min), per_unit_min), collapse = ", "))
say("  => K (matched depth) = %d (binding unit: %s)", K,
    names(per_unit_min)[which.min(per_unit_min)])
stopifnot(K == K_EXPECTED)   # N2-locked; a mismatch signals substrate drift

## ---- 2. NATIVE fits: pseudobulk (fast) + NEBULA (heavy, per-unit cached) -
say("[1/6] NATIVE pseudobulk x %d units ...", length(UNITS))
native_pb <- setNames(lapply(UNITS, function(u) {
  fo <- subset_pseudobulk_de(sc, u, min_count = MIN_COUNT, symbol_map = NULL)
  say("  pb NATIVE %-15s %6d cells | %5d genes", u, fo$n_cells, fo$n_genes)
  fo
}), UNITS)

# NEBULA lightest-first so a crash on the 163k Neuronal unit still leaves the
# other 5 cached. Each fit -> its own resume file (skipped if present).
order_light <- names(sort(unit_counts))
say("[2/6] NATIVE NEBULA x %d units (lightest-first: %s) ...",
    length(UNITS), paste(order_light, collapse = " < "))
native_nb <- list()
for (u in order_light) {
  f <- file.path(fit_dir, sprintf("native_nebula_%s.rds", u))
  if (file.exists(f) && !overwrite) {
    native_nb[[u]] <- readRDS(f)
    say("  nebula NATIVE %-15s RESUMED from cache (%d cells)", u, native_nb[[u]]$n_cells)
  } else {
    t0 <- Sys.time()
    say("  nebula NATIVE %-15s %6d cells | fitting on %d cores ...",
        u, as.integer(unit_counts[u]), NCORE)
    fo <- subset_nebula_de(sc, u, min_cell_frac = MIN_CELL_FRAC,
                           ncore = NCORE, symbol_map = NULL)
    saveRDS(fo, f, compress = "xz"); Sys.chmod(f, "0644")
    native_nb[[u]] <- fo
    say("    -> %d genes | %d/%d converged | %.1f min", fo$n_genes,
        sum(fo$fit$convergence == 1, na.rm = TRUE), length(fo$fit$convergence),
        as.numeric(difftime(Sys.time(), t0, units = "mins")))
  }
}
native_nb <- native_nb[UNITS]   # canonical (alphabetical) order

## ---- 3. MATCHED fits: downsample to K, then pseudobulk + NEBULA ----------
say("[3/6] MATCHED regime: downsample to K=%d/replicate (seed=%d) ...", K, SEED)
sc_m <- downsample_balanced(sc, K = K, seed = SEED, units = UNITS)
mcounts <- table(factor(as.character(sc_m$unit), levels = UNITS))
say("  matched object: %d cells | per-unit: %s", ncol(sc_m),
    paste(sprintf("%s=%d", names(mcounts), as.integer(mcounts)), collapse = ", "))
stopifnot(max(table(sc_m$unit, sc_m$genotype_batch)) <= K)

matched_pb <- setNames(lapply(UNITS, function(u)
  subset_pseudobulk_de(sc_m, u, min_count = MIN_COUNT, symbol_map = NULL)), UNITS)
say("  matched pseudobulk: %d units fit", length(matched_pb))
matched_nb <- setNames(lapply(UNITS, function(u) {
  t0 <- Sys.time()
  fo <- subset_nebula_de(sc_m, u, min_cell_frac = MIN_CELL_FRAC,
                         ncore = NCORE, symbol_map = NULL)
  say("  nebula MATCHED %-15s %5d cells | %5d genes | %.1f min", u, fo$n_cells,
      fo$n_genes, as.numeric(difftime(Sys.time(), t0, units = "mins")))
  fo
}), UNITS)

## ---- 4. assemble R1-R4 ---------------------------------------------------
say("[4/6] assembling R1-R4 (tally / concordance / specificity class) ...")
fits <- list(
  native  = list(pseudobulk = native_pb,  nebula = native_nb),
  matched = list(pseudobulk = matched_pb, nebula = matched_nb))
tabs <- assemble_specificity_tables(fits, class_contrast = "interaction",
                                    class_fdr = PADJ_CUT)
tally_native  <- dplyr::filter(tabs$tally, regime == "native")
tally_matched <- dplyr::filter(tabs$tally, regime == "matched")
# Headline read (state BEFORE any verdict prose, N4): interaction breadth per
# unit under the LOAD-BEARING matched regime.
hl <- tally_matched |>
  dplyr::filter(contrast == "interaction") |>
  dplyr::group_by(unit) |>
  dplyr::summarise(nebula_sig10 = sum(n_sig_10[estimator == "nebula"]),
                   pb_sig10 = sum(n_sig_10[estimator == "pseudobulk"]),
                   .groups = "drop") |>
  dplyr::arrange(dplyr::desc(nebula_sig10))
say("  MATCHED interaction sig genes (FDR<0.10) by unit [nebula | pb]:")
for (i in seq_len(nrow(hl)))
  say("    %-15s %4d | %4d", hl$unit[i], hl$nebula_sig10[i], hl$pb_sig10[i])

## ---- 5. microglia cross-check (report all 5 rho; hard floor 0.90) --------
say("[5/6] microglia sanity cross-check vs canonical de_snrnaseq[_nebula] ...")
de_sn <- readRDS(file.path(cc, "de_snrnaseq.rds"))
de_nb <- readRDS(file.path(cc, "de_snrnaseq_nebula.rds"))
xc_pb <- microglia_crosscheck(native_pb$Microglia, de_sn$fit$top, "pseudobulk",
                              min_rho = XCHECK_FLOOR)
xc_nb <- microglia_crosscheck(native_nb$Microglia, de_nb$top, "nebula",
                              min_rho = XCHECK_FLOOR)
crosscheck <- dplyr::bind_rows(xc_pb, xc_nb)
for (i in seq_len(nrow(crosscheck)))
  say("  %-11s %-15s rho=%.3f  n=%d  %s", crosscheck$estimator[i],
      crosscheck$contrast[i], crosscheck$spearman[i], crosscheck$n_shared[i],
      ifelse(crosscheck$passes[i], "PASS", "*** BELOW 0.90 ***"))
stopifnot(all(crosscheck$passes))   # gross-mismatch gate: any contrast <=0.90 fails the build

## ---- 6. R5: GO BP fGSEA per unit at interaction --------------------------
say("[6/6] R5 pathway-level GO BP fGSEA per unit at interaction ...")
gobp <- get_gobp()
set.seed(SEED)   # fgseaMultilevel samples; fix for reproducibility
grid <- expand.grid(regime = c("native", "matched"),
                    estimator = c("pseudobulk", "nebula"),
                    stringsAsFactors = FALSE)
pathway_long <- purrr::pmap_dfr(grid, function(regime, estimator)
  specificity_pathway_fgsea(fits[[regime]][[estimator]], gobp,
                            regime = regime, estimator = estimator))
pathway_tally <- specificity_pathway_tally(pathway_long)
mg_pathway_cross_unit <- microglia_pathway_cross_unit(pathway_long, fdr = PADJ_CUT)
say("  GO BP sets=%d | fgsea rows=%d | microglia-sig pathways tracked=%d",
    length(gobp), nrow(pathway_long),
    length(unique(mg_pathway_cross_unit$pathway)))
mph <- pathway_tally |> dplyr::filter(regime == "matched", estimator == "nebula") |>
  dplyr::arrange(dplyr::desc(n_path_10))
say("  MATCHED/nebula interaction-enriched GO BP pathways (FDR<0.10) by unit:")
for (i in seq_len(nrow(mph)))
  say("    %-15s %4d", mph$unit[i], mph$n_path_10[i])

## ---- headline gene set (R3b; assembled here, plotted in rmd/23) ----------
# DAM activation programme + microglia-state markers (the M-002 progression
# drivers ARE the activation programme) UNION the top-N microglia NEBULA
# interaction symbols (the H-phase headline interaction genes). Reproducible
# + provenance-stamped so rmd/23 is a pure consumer.
dam_programme <- unique(unlist(
  canonical_microglia_markers[c("Microglia", "DAM", "IFN", "Proliferative")],
  use.names = FALSE))
mg_int_top <- de_nb$top$interaction |>
  dplyr::filter(!is.na(symbol), symbol != "") |>
  dplyr::arrange(adj.P.Val) |> head(HEADLINE_N) |> dplyr::pull(symbol)
headline_genes <- unique(c(dam_programme, mg_int_top))
say("  headline gene set: %d symbols (%d DAM/state markers + top-%d microglia interaction)",
    length(headline_genes), length(dam_programme), HEADLINE_N)

## ---- write results -------------------------------------------------------
say("[write] results TSVs + cache ...")
wtsv(tally_native,                "celltype_specificity_tally_native.tsv")
wtsv(tally_matched,               "celltype_specificity_tally_matched.tsv")
wtsv(tabs$interaction_concordance,"celltype_specificity_interaction_concordance.tsv")
wtsv(tabs$specificity_class,      "celltype_specificity_class.tsv")
wtsv(pathway_tally,               "celltype_specificity_pathway_tally.tsv")
wtsv(crosscheck,                  "celltype_specificity_crosscheck.tsv")

## ---- assemble cache for rmd/23 ------------------------------------------
celltype_specificity <- list(
  # slimmed fit grid (top tables only) for the rmd/23 heatmaps + scatter
  fits = list(
    native  = list(pseudobulk = slim_list(native_pb),  nebula = slim_list(native_nb)),
    matched = list(pseudobulk = slim_list(matched_pb), nebula = slim_list(matched_nb))),
  tally = tabs$tally,                                  # R1+R2 stacked
  interaction_concordance = tabs$interaction_concordance,  # R3a + R3c
  specificity_class = tabs$specificity_class,          # R4
  pathway_tally = pathway_tally,                       # R5 headline
  microglia_pathway_cross_unit = mg_pathway_cross_unit,# R5 interpretive
  crosscheck = crosscheck,                             # sanity gate
  headline_genes = headline_genes,                     # R3b curated set
  unit_cell_counts = as.data.frame(unit_counts) |>
    setNames(c("unit", "n_cells")),
  per_unit_min = data.frame(unit = names(per_unit_min),
                            min_cells_per_replicate = as.integer(per_unit_min)),
  params = list(
    K = K, seed = SEED, min_count = MIN_COUNT, min_cell_frac = MIN_CELL_FRAC,
    padj_cut = PADJ_CUT, xcheck_floor = XCHECK_FLOOR, headline_n = HEADLINE_N,
    units = UNITS, n_cells = N_CELLS_EXP, ncore = NCORE,
    design = "~ 0 + genotype + batch", contrasts = specificity_contrasts,
    estimators = c("pseudobulk", "nebula"), regimes = c("native", "matched"),
    headline_provenance = "canonical DAM/state markers + top-N de_snrnaseq_nebula interaction symbols"))
saveRDS(celltype_specificity, out_cache, compress = "xz")
Sys.chmod(out_cache, "0644")
say("[build_celltype_specificity] wrote cache: %s", out_cache)

## ---- validation summary -------------------------------------------------
# 6 units x 5 contrasts x 2 estimators x 2 regimes = 120 tally rows.
stopifnot(
  nrow(tabs$tally) == length(UNITS) * length(specificity_contrasts) * 2L * 2L,
  setequal(unique(tabs$tally$unit), UNITS),
  all(c("regime", "estimator", "unit", "contrast", "n_sig_10") %in% colnames(tabs$tally)),
  all(crosscheck$passes),                              # cross-check floor held
  nrow(pathway_tally) == length(UNITS) * 2L * 2L,      # 6 units x 2 est x 2 regimes
  length(headline_genes) > 0L)
say("[validation] 120 tally rows | 6 units | cross-check >0.90 | %d pathway-tally rows | %d headline genes",
    nrow(pathway_tally), length(headline_genes))
say("[build_celltype_specificity] done.")

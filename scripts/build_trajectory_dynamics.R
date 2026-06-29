#!/usr/bin/env Rscript
# arc-O O3: genome-wide gene-level differential DYNAMICS at the tau x amyloid
# interaction along the clean homeostatic->DAM activation pseudotime.
#
# Question (locked O1): is there a gene-level tau x amyloid interaction in the
# DYNAMICS of expression along pseudotime - genes whose amyloid response is
# localised to a pseudotime window and is tau-modulated - that the STATIC
# matched-power pseudobulk contrast (arc N, ~0 genes) cannot see, naming the
# programme behind M-002 (the aggregate mean-pt synergy)? A null (also ~0) is a
# clean finding: the synergy is aggregate-rate-only.
#
# Readout: tradeSeq NB-GAM (one lineage, conditions = 4 genotypes, nknots = 6)
# with the custom 2x2 difference-of-differences Wald contrast in coefficient
# space (R/trajectory.R O2 helpers; == tradeSeq internal, smoke-test proven).
# Reported side-by-side with the omnibus conditionTest and the associationTest
# (anti-cherry-pick guardrail). Heavy resumable fitGAM; display chapter rmd/22
# (O4) readRDS-loads trajectory_dynamics.rds and recomputes nothing.
#
# Run: Rscript scripts/build_trajectory_dynamics.R [--overwrite] [--refit] [--ncores N]
#   --overwrite  rebuild trajectory_dynamics.rds even if it exists
#   --refit      re-run fitGAM even if the resumable SCE cache exists
#   --ncores N   fitGAM workers (default detectCores()-2)

suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(Matrix)
  library(SummarizedExperiment); library(tradeSeq)
  library(BiocParallel); library(dplyr); library(tibble)
})
setwd("/home/rstudio/tau-mutant-integration-ng")
source("R/helpers.R")

args      <- commandArgs(trailingOnly = TRUE)
overwrite <- "--overwrite" %in% args
refit     <- "--refit" %in% args
ncix      <- which(args == "--ncores")
ncores    <- if (length(ncix)) as.integer(args[ncix + 1L]) else
               max(1L, parallel::detectCores() - 2L)

cc  <- "storage/cache"; res <- "storage/results"
dir.create(res, recursive = TRUE, showWarnings = FALSE)
out_cache <- file.path(cc, "trajectory_dynamics.rds")
sce_cache <- file.path(cc, "trajectory_dynamics_sce.rds")   # resumable fitGAM (gitignored)

## ---- embedded params (locked O1; state before applying) -----------------
NKNOTS       <- 6L      # spline knots per genotype smoother
PADJ_CUT     <- 0.10    # project-standard FDR
SEED         <- 1L
MIN_FRAC     <- 0.01    # gene set = detected in >= 1% of ON-LINEAGE cells (genome-wide, O1 choice 2)
N_TOP_SMOOTH <- 20L     # top interaction genes to store per-genotype smoothers for (O4 plots)
K_SMOOTH     <- 100L    # eval points for the plot smoothers (plot-only; not the contrast)

if (file.exists(out_cache) && !overwrite) {
  cat(sprintf("[build_trajectory_dynamics] cache exists, skipping: %s\n", out_cache))
  cat("Pass --overwrite to rebuild.\n"); quit(save = "no", status = 0)
}
set.seed(SEED)
say <- function(...) cat(sprintf(...), "\n")
wtsv <- function(df, name) {
  p <- file.path(res, name)
  write_tsv_safe(as.data.frame(df), p); Sys.chmod(p, "0644")
  say("  wrote %-44s %d rows", name, nrow(df))
}

## ---- 1. load substrate --------------------------------------------------
say("[build_trajectory_dynamics] loading trajectory.rds + microglia counts + symbol_map ...")
tj    <- readRDS(file.path(cc, "trajectory.rds"))
micro <- readRDS(file.path(cc, "microglia_seurat_processed.rds"))
sm    <- readRDS(file.path(cc, "snrnaseq_symbol_map.rds"))
pc    <- tj$per_cell

counts <- as(GetAssayData(micro, assay = "RNA", layer = "counts"), "CsparseMatrix")
stopifnot(all(pc$cell_id %in% colnames(counts)))
counts <- counts[, pc$cell_id, drop = FALSE]      # align cells to per_cell order
pt        <- pc$pt_clean
genotype  <- pc$genotype
on        <- is.finite(pt)                        # clean homeostatic->DAM lineage (IFN/prolif omitted)
say("cells total=%d  on-lineage(finite pt_clean)=%d", length(pt), sum(on))

cnt_by_geno <- table(factor(genotype[on], levels = genotype_levels))
say("on-lineage cells per genotype (power check):")
print(cnt_by_geno)
on_lineage_counts <- tibble(genotype = names(cnt_by_geno),
                            n_cells  = as.integer(cnt_by_geno))

## ---- 2. locked gene set: genome-wide >= 1%-expressed on-lineage ---------
det      <- Matrix::rowSums(counts[, on, drop = FALSE] > 0)
min_cell <- ceiling(MIN_FRAC * sum(on))
gene_set <- rownames(counts)[det >= min_cell]
say("gene set: detected in >= %d cells (%.0f%% of %d on-lineage) -> %d / %d genes",
    min_cell, 100 * MIN_FRAC, sum(on), length(gene_set), nrow(counts))

## ---- 3. fitGAM (resumable) ----------------------------------------------
if (file.exists(sce_cache) && !refit) {
  say("[fit] loading cached fitGAM SCE: %s (pass --refit to refit)", sce_cache)
  sce <- readRDS(sce_cache)
  stopifnot(setequal(rownames(sce), gene_set))
} else {
  say("[fit] fitGAM on %d genes x %d on-lineage cells, nknots=%d, %d worker(s) ...",
      length(gene_set), sum(on), NKNOTS, ncores)
  BPPARAM <- if (ncores > 1L) MulticoreParam(workers = ncores, RNGseed = SEED) else SerialParam()
  t0  <- proc.time()[3]
  sce <- fit_lineage_gam(counts[gene_set, , drop = FALSE], pt, genotype,
                         nknots = NKNOTS, parallel = ncores > 1L, BPPARAM = BPPARAM)
  say("[fit] done: %d genes x %d cells in %.1f min",
      nrow(sce), ncol(sce), (proc.time()[3] - t0) / 60)
  saveRDS(sce, sce_cache, compress = "gzip"); Sys.chmod(sce_cache, "0644")
  say("[fit] saved resumable SCE -> %s", sce_cache)
}

## ---- 4. contrasts: interaction (targeted) + omnibus + association -------
say("[contrast] 2x2 interaction-dynamics Wald + omnibus conditionTest + associationTest ...")
interaction_tbl <- interaction_dynamics_contrast(sce, symbol_map = sm)
omnibus_tbl     <- omnibus_dynamics(sce, symbol_map = sm)
association_tbl <- association_dynamics(sce, symbol_map = sm)

n_fit <- sum(is.finite(interaction_tbl$waldStat))
n_sig <- sum(interaction_tbl$adj.P.Val < PADJ_CUT, na.rm = TRUE)
say("interaction: %d/%d genes fitted, %d with adj.P < %.2f",
    n_fit, nrow(interaction_tbl), n_sig, PADJ_CUT)
say("omnibus    : %d with adj.P < %.2f",
    sum(omnibus_tbl$adj.P.Val < PADJ_CUT, na.rm = TRUE), PADJ_CUT)
say("association: %d with adj.P < %.2f",
    sum(association_tbl$adj.P.Val < PADJ_CUT, na.rm = TRUE), PADJ_CUT)
say("top interaction genes:")
print(as.data.frame(head(
  interaction_tbl[, c("symbol", "waldStat", "df", "pvalue", "adj.P.Val", "effect_peak")], 12)),
  row.names = FALSE)

## ---- 5. smoothers for top interaction genes (O4 plots) ------------------
top_genes <- head(interaction_tbl$gene[is.finite(interaction_tbl$waldStat)], N_TOP_SMOOTH)
say("[smooth] extracting per-genotype smoothers for top %d interaction genes ...", length(top_genes))
smoothers <- lapply(top_genes, function(g) extract_condition_smoothers(sce, g, K = K_SMOOTH))
names(smoothers) <- top_genes

## ---- 6. static-vs-dynamic honesty join (arc N static + arc M synergy) ---
# arc N: matched-power static interaction in the Microglia unit (NEBULA single-cell
# + pseudobulk), collapses to ~0 genes. arc-N tables are SYMBOL-keyed (built from
# seurat_full_processed with symbol_map = NULL); the dynamics table is Ensembl-keyed
# -> bridge on MGI symbol. Collapse each static table to one row per symbol
# (strongest static signal = min P.Value) so the join cannot duplicate dynamics
# rows; the chapter can then state plainly whether DYNAMICS recovers genes the
# STATIC contrast could not.
cs <- readRDS(file.path(cc, "celltype_specificity.rds"))
collapse_static <- function(df) {
  as_tibble(df) %>% filter(!is.na(symbol)) %>%
    group_by(symbol) %>% slice_min(P.Value, n = 1, with_ties = FALSE) %>% ungroup()
}
stat_nb <- collapse_static(cs$fits$matched$nebula$Microglia$top[["interaction"]])
stat_pb <- collapse_static(cs$fits$matched$pseudobulk$Microglia$top[["interaction"]])
vs_static <- interaction_tbl %>%
  transmute(gene, symbol,
            dyn_waldStat = waldStat, dyn_df = df, dyn_pvalue = pvalue,
            dyn_adj.P.Val = adj.P.Val, dyn_effect_l2 = effect_l2) %>%
  left_join(transmute(stat_nb, symbol,
                      stat_nb_t = t, stat_nb_adj.P.Val = adj.P.Val), by = "symbol") %>%
  left_join(transmute(stat_pb, symbol,
                      stat_pb_t = t, stat_pb_adj.P.Val = adj.P.Val), by = "symbol") %>%
  arrange(dyn_pvalue)

# arc M: aggregate mean-pseudotime synergy (single scalar set; M-002 = the
# progression interaction). Carried for the chapter's vs-M-002 comparison.
arc_m_synergy <- tj$contrasts %>%
  filter(contrast == "interaction",
         measure %in% c("mean_pt", "progression_cf", "composition_cf")) %>%
  transmute(measure, logFC, se, t, P.Value, adj.P.Val, sig)

## ---- 7. assemble + write ------------------------------------------------
params <- list(
  nknots = NKNOTS, padj_cut = PADJ_CUT, seed = SEED, min_frac = MIN_FRAC,
  min_cell = min_cell, gene_set_size = length(gene_set),
  n_top_smoothers = length(top_genes), k_smooth = K_SMOOTH,
  conditions = genotype_levels, lineage = "clean_homeostatic_DAM",
  contrast = "diff_of_diffs_wald_coefficient_space",
  static_comparator = "celltype_specificity matched Microglia interaction (nebula + pseudobulk)",
  ncores = ncores, tradeseq_version = as.character(packageVersion("tradeSeq")))

trajectory_dynamics <- list(
  interaction = interaction_tbl, omnibus = omnibus_tbl, association = association_tbl,
  smoothers = smoothers, vs_static = vs_static, arc_m_synergy = arc_m_synergy,
  gene_set = gene_set, on_lineage_counts = on_lineage_counts,
  n_fitted = n_fit, n_sig_interaction = n_sig, params = params)
saveRDS(trajectory_dynamics, out_cache, compress = "xz"); Sys.chmod(out_cache, "0644")
say("[write] saved %s", out_cache)

wtsv(interaction_tbl, "trajectory_dynamics_interaction.tsv")
wtsv(omnibus_tbl,     "trajectory_dynamics_omnibus.tsv")
wtsv(association_tbl, "trajectory_dynamics_association.tsv")
wtsv(vs_static,       "trajectory_dynamics_vs_static.tsv")

say("[build_trajectory_dynamics] DONE: %d interaction-dynamics genes at FDR < %.2f (of %d fitted).",
    n_sig, PADJ_CUT, n_fit)

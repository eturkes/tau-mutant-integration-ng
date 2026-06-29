#!/usr/bin/env Rscript
# M3 driver (plan arc M): infer a microglial homeostatic -> DAM ACTIVATION
# pseudotime from the snRNAseq and push a per-REPLICATE (genotype_batch, the
# locked 16-id unit) pseudotime summary through the project's locked 5-contrast
# 2x2 factorial -- the project's first DYNAMICS / progression readout (every
# prior arc reads a static quantity: expression DE, pathway/TF/kinase/CCC/
# causal/SCENIC activity, or discrete cell composition). Thin I/O wrapper over
# the pure fns in R/trajectory.R; mirrors build_spatial_deconvolution.R -- all
# heavy compute OUT of the knit, rmd/22 (M4) is display-only. Idempotent unless
# --overwrite.
#
# Maximal-triangulation method set (M1 gate): Slingshot (R, PRIMARY per-lineage
# pseudotime) + scanpy PAGA/DPT + CellRank2 DAM-fate (both velocity-free, via
# the build_trajectory_python.py sub-step) + CytoTRACE2 / entropy / n_genes for
# root validation; tradeSeq differential-dynamics is implemented but GATED behind
# --tradeseq (fitGAM is heavy and not part of the core cache; default OFF).
#
# THE load-bearing piece (guardrail #1): the per-replicate mean-pseudotime shift
# is DECOMPOSED into a between-state COMPOSITION channel (cells reshuffling into
# the already-known larger DAM cluster = the arc-L result) and a within-state
# PROGRESSION channel (genuine activation advancement) so the two are never
# conflated. A null / additive interaction is a valid finding (guardrail #4).
#
# Policy locked BEFORE inspecting which contrast moves (guardrail #8), per M1/M2
# findings (substates are transcriptionally close; cell-cycle pulls a spurious
# intermediate at full harmony dims):
#   * embedding   = harmony truncated to N_DIMS = 15 (M2's choice)
#   * HEADLINE    = the CLEAN homeostatic->DAM lineage (IFN + proliferative
#                   omitted as cycling/secondary side-states): the locked
#                   "lineage of interest", least confounded progression axis.
#   * topology    = an ALL-STATES run (all 4 states, DAM forced terminal) is
#                   ALSO computed + reported (guardrail #7: report all lineages)
#                   and feeds the cross-method concordance.
#   * root        = homeostatic, VALIDATED (not asserted) by 3 potency proxies.
#   * unit/FDR    = genotype_batch (16); adj.P < 0.10; the 5 canonical contrasts.
#
# Outputs (storage/results/*.tsv + storage/cache/trajectory.rds):
#   trajectory_pseudotime_by_genotype.tsv     per-genotype mean of each measure
#   trajectory_contrasts.tsv                  5-contrast factorial, all measures
#   trajectory_progression_decomposition.tsv  per-replicate composition vs progression
#   trajectory_method_concordance.tsv         cross-method pseudotime Spearman
#
# Run: Rscript scripts/build_trajectory.R [--overwrite] [--tradeseq]

suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(Matrix)
  library(data.table); library(dplyr); library(tibble)
})

setwd("/home/rstudio/tau-mutant-integration-ng")
source("R/helpers.R")

args       <- commandArgs(trailingOnly = TRUE)
overwrite  <- "--overwrite" %in% args
run_tradeseq <- "--tradeseq" %in% args
cc  <- "storage/cache"
res <- "storage/results"
export_dir <- file.path(cc, "trajectory")        # R<->Python bundle (gitignored)
dir.create(res, recursive = TRUE, showWarnings = FALSE)
dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)
out_cache <- file.path(cc, "trajectory.rds")

## ---- embedded params (state before applying; M1/M2 locked) --------------
N_DIMS      <- 15L       # harmony dims (M2: full 30 lets cell-cycle intermediate)
N_NEIGHBORS <- 30L       # scanpy neighbours (python lane)
N_DCS       <- 10L       # DPT diffusion components
PADJ_CUT    <- 0.10      # project-standard FDR
SEED        <- 1L
PYTHON      <- ".venv/bin/python"
N_CELLS_EXP <- 26104L

if (file.exists(out_cache) && !overwrite) {
  cat(sprintf("[build_trajectory] cache exists, skipping: %s\n", out_cache))
  cat("Pass --overwrite to rebuild.\n"); quit(save = "no", status = 0)
}

say <- function(...) cat(sprintf(...), "\n")
wtsv <- function(df, name) {
  p <- file.path(res, name)
  write_tsv_safe(as.data.frame(df), p); Sys.chmod(p, "0644")
  say("  wrote %-46s %d rows", name, nrow(df))
}

## ---- 1. load + label states ---------------------------------------------
say("[build_trajectory] loading microglia_seurat_processed + symbol_map ...")
micro      <- readRDS(file.path(cc, "microglia_seurat_processed.rds"))
symbol_map <- readRDS(file.path(cc, "snrnaseq_symbol_map.rds"))
stopifnot(ncol(micro) == N_CELLS_EXP)

micro <- label_microglia_states(micro, symbol_map)       # deterministic; not cached
states <- as.character(micro$state)
meta   <- micro@meta.data
gb     <- as.character(meta$genotype_batch)
cells  <- colnames(micro)
say("  %d cells | states: %s", length(cells),
    paste(sprintf("%s=%d", names(table(states)), as.integer(table(states))),
          collapse = ", "))

# Replicate (genotype_batch) meta: rownames = the 16 ids, cols genotype/batch.
replicate_meta <- meta |>
  dplyr::distinct(genotype_batch, genotype, batch) |>
  dplyr::arrange(genotype_batch) |>
  tibble::remove_rownames() |>
  tibble::column_to_rownames("genotype_batch")
replicate_meta$genotype <- factor(as.character(replicate_meta$genotype),
                                  levels = genotype_levels)
rep_ids <- rownames(replicate_meta)
stopifnot(length(rep_ids) == 16L)
say("  %d replicates (genotype_batch); %d per genotype",
    length(rep_ids), length(rep_ids) / 4L)

## ---- embeddings ---------------------------------------------------------
harmony <- Embeddings(micro, "harmony")
stopifnot(nrow(harmony) == length(cells), ncol(harmony) >= N_DIMS)
harmony15 <- harmony[, seq_len(N_DIMS), drop = FALSE]
umap <- Embeddings(micro, "umap")
say("  harmony %d x %d -> truncated to %d dims | umap %d x %d",
    nrow(harmony), ncol(harmony), N_DIMS, nrow(umap), ncol(umap))

## ---- 2. Slingshot: clean headline + all-states topology -----------------
say("[1/6] Slingshot x2 (clean homeostatic->DAM headline + all-states topology)...")
traj_clean <- build_microglia_trajectory(
  harmony15, states, start_clus = "homeostatic", end_clus = "DAM",
  omit_clusters = c("IFN", "proliferative"), dam_state = "DAM", seed = SEED)
traj_all <- build_microglia_trajectory(
  harmony15, states, start_clus = "homeostatic", end_clus = "DAM",
  omit_clusters = NULL, dam_state = "DAM", seed = SEED)
say("  clean lineages: %s", paste(sapply(traj_clean$lineages, paste, collapse = "->"),
                                  collapse = " | "))
say("  all   lineages: %s", paste(sapply(traj_all$lineages, paste, collapse = "->"),
                                  collapse = " | "))
stopifnot(!is.na(traj_clean$dam_lineage), !is.na(traj_all$dam_lineage))

pt_clean <- traj_clean$pseudotime[, traj_clean$dam_lineage]
w_clean  <- traj_clean$weights[,   traj_clean$dam_lineage]
pt_all   <- traj_all$pseudotime[,   traj_all$dam_lineage]
w_all    <- traj_all$weights[,      traj_all$dam_lineage]
# Monotonicity sanity: state mean pseudotime should rise homeostatic -> DAM.
sm_clean <- tapply(pt_clean, states, mean, na.rm = TRUE)
say("  clean state-mean pt: %s",
    paste(sprintf("%s=%.1f", names(sm_clean), sm_clean), collapse = ", "))
stopifnot(sm_clean["DAM"] > sm_clean["homeostatic"])

## ---- 3. root validation: CytoTRACE2 + entropy + n_genes -----------------
say("[2/6] root-potency validation (CytoTRACE2 + entropy + n_genes)...")
raw_counts <- GetAssayData(micro, assay = "RNA", layer = "counts")
raw_counts <- as(raw_counts, "CsparseMatrix")

# entropy / n_genes are gene-identity-agnostic -> use raw Ensembl counts.
pot_entropy <- cell_potency(raw_counts, "entropy")
pot_ngenes  <- cell_potency(raw_counts, "n_genes")

# CytoTRACE2 needs MGI symbols -> collapse Ensembl rows by summing (SCENIC pattern).
collapse_to_symbols <- function(counts, smap) {
  sym <- smap$symbol[match(rownames(counts), smap$ensembl)]
  keep <- !is.na(sym) & nzchar(sym)
  counts <- counts[keep, , drop = FALSE]; sym <- sym[keep]
  f <- factor(sym, levels = unique(sym))
  inc <- Matrix::fac2sparse(f, to = "d")
  out <- as(inc %*% counts, "CsparseMatrix"); rownames(out) <- levels(f); out
}
pot_cyto <- tryCatch({
  sym_counts <- collapse_to_symbols(raw_counts, symbol_map)
  say("  CytoTRACE2 input: %d symbols x %d cells", nrow(sym_counts), ncol(sym_counts))
  ct <- CytoTRACE2::cytotrace2(
    as.matrix(sym_counts), species = "mouse", is_seurat = FALSE,
    slot_type = "counts", ncores = 4L, seed = SEED)
  score_col <- grep("^CytoTRACE2_Score$", colnames(ct), value = TRUE)
  if (!length(score_col)) score_col <- grep("Score", colnames(ct), value = TRUE)[1]
  stats::setNames(ct[[score_col]], rownames(ct))[cells]
}, error = function(e) { say("  CytoTRACE2 FAILED (%s) -> entropy/n_genes only",
                             conditionMessage(e)); NULL })

potencies <- list(entropy = pot_entropy, n_genes = pot_ngenes)
if (!is.null(pot_cyto)) potencies$cytotrace2 <- pot_cyto
root_validation <- lapply(potencies, function(p)
  validate_root_potency(p[cells], states, root = "homeostatic"))
for (nm in names(root_validation)) {
  rv <- root_validation[[nm]]
  say("  %-11s root rank=%s/%d  most_potent=%s  delta=%.4f", nm,
      rv$root_rank, rv$n_states, rv$root_is_most_potent, rv$delta_to_next)
}

## ---- 4. Python lane: PAGA/DPT + CellRank2 (velocity-free) ----------------
say("[3/6] exporting bundle + invoking scanpy PAGA/DPT + CellRank2 ...")
emb_df <- data.frame(cell_id = cells, harmony15, check.names = FALSE)
colnames(emb_df)[-1] <- paste0("harmony_", seq_len(N_DIMS))
fwrite(emb_df, file.path(export_dir, "embedding.csv"))
fwrite(data.frame(cell_id = cells, state = states,
                  genotype = as.character(meta$genotype),
                  batch = as.character(meta$batch), genotype_batch = gb),
       file.path(export_dir, "obs.csv"))

py_log <- system2(PYTHON,
  c("scripts/build_trajectory_python.py", "--export-dir", export_dir,
    "--n-neighbors", N_NEIGHBORS, "--n-dcs", N_DCS, "--seed", SEED - 1L),
  stdout = TRUE, stderr = TRUE)
cat(paste0("    | ", py_log), sep = "\n")
py_status <- attr(py_log, "status") %||% 0L
py_csv <- file.path(export_dir, "python_pseudotime.csv")
if (!identical(py_status, 0L) && !is.null(attr(py_log, "status")))
  say("  WARNING: python lane exit status %s", py_status)
stopifnot(file.exists(py_csv))
py_pt <- fread(py_csv)
py_pt <- py_pt[match(cells, py_pt$cell_id), ]
dpt <- py_pt$dpt_pseudotime
cellrank_dam_fate <- py_pt$cellrank_dam_fate
say("  DPT finite=%d/%d | CellRank DAM-fate finite=%d/%d",
    sum(is.finite(dpt)), length(dpt),
    sum(is.finite(cellrank_dam_fate)), length(cellrank_dam_fate))

## ---- 5. cross-method concordance ----------------------------------------
say("[4/6] cross-method pseudotime concordance (Spearman)...")
pt_df <- data.frame(slingshot_clean = pt_clean, slingshot_all = pt_all,
                    dpt = dpt, cellrank_dam_fate = cellrank_dam_fate)
concordance <- pseudotime_concordance(pt_df)
print(concordance$pairwise)

## ---- 6. per-replicate summary + decomposition + factorial ---------------
say("[5/6] per-replicate summary + progression/composition decomposition + factorial...")
# Headline = CLEAN homeostatic->DAM lineage (locked lineage of interest); the
# omitted IFN/proliferative carry no on-lineage cells, so restrict within_* to
# the two traversed states (avoids all-NA within_IFN/prolif rows).
pr_clean <- pseudotime_per_replicate(pt_clean, states, gb, weights = w_clean,
                                     dam_state = "DAM",
                                     state_levels = c("homeostatic", "DAM"),
                                     replicate_meta = replicate_meta)
pr_all   <- pseudotime_per_replicate(pt_all, states, gb, weights = w_all,
                                     dam_state = "DAM", replicate_meta = replicate_meta)
decomp   <- decompose_progression_vs_composition(pt_clean, states, gb)

# per-replicate CellRank DAM-fate mean (defined for all cells).
cr_by_rep <- tapply(cellrank_dam_fate, gb, mean, na.rm = TRUE)

# Assemble the measure x replicate matrix (cols = the 16 ids).
S <- pr_clean$by_replicate |>
  dplyr::select(replicate, mean_pt, frac_past, dplyr::starts_with("within_")) |>
  dplyr::left_join(dplyr::select(decomp$by_replicate, replicate,
                                 observed, composition_cf, progression_cf),
                   by = "replicate") |>
  dplyr::left_join(dplyr::transmute(pr_all$by_replicate, replicate,
                                    mean_pt_all = mean_pt), by = "replicate") |>
  dplyr::left_join(data.frame(replicate = names(cr_by_rep),
                              cellrank_dam_fate = as.numeric(cr_by_rep)),
                   by = "replicate") |>
  dplyr::arrange(match(replicate, rep_ids))
stopifnot(identical(S$replicate, rep_ids))

# Fit pseudotime-scale measures and bounded [0,1] measures in SEPARATE eBayes
# calls (mixed scales would distort the shared limma-trend variance prior).
pt_scale <- c("mean_pt", "mean_pt_all", "observed", "composition_cf",
              "progression_cf", grep("^within_", colnames(S), value = TRUE))
bounded  <- c("frac_past", "cellrank_dam_fate")
mat_of <- function(cols) {
  m <- t(as.matrix(S[, cols, drop = FALSE])); colnames(m) <- S$replicate; m
}
ct_pt <- fit_trajectory_contrasts(mat_of(pt_scale), replicate_meta,
                                  transform = "none", padj_cut = PADJ_CUT) |>
  dplyr::mutate(measure_scale = "pseudotime", .before = 1)
ct_bd <- fit_trajectory_contrasts(mat_of(bounded), replicate_meta,
                                  transform = "none", padj_cut = PADJ_CUT) |>
  dplyr::mutate(measure_scale = "bounded_0_1", .before = 1)
contrasts_tbl <- dplyr::bind_rows(ct_pt, ct_bd)
say("  contrast rows: %d | sig (adj.P<%.2f): %d | at interaction: %d",
    nrow(contrasts_tbl), PADJ_CUT, sum(contrasts_tbl$sig),
    sum(contrasts_tbl$sig & contrasts_tbl$contrast == "interaction"))
# Headline read (state BEFORE the verdict prose, M4): which channel carries it.
hl <- contrasts_tbl |>
  dplyr::filter(contrast == "interaction",
                measure %in% c("mean_pt", "composition_cf", "progression_cf"))
say("  interaction by channel: %s",
    paste(sprintf("%s logFC=%.2f adj.P=%.3f", hl$measure, hl$logFC, hl$adj.P.Val),
          collapse = " | "))

## ---- per-genotype summary (all measures, even-handed) -------------------
by_genotype <- S |>
  dplyr::left_join(tibble::rownames_to_column(replicate_meta, "replicate"),
                   by = "replicate") |>
  tidyr::pivot_longer(cols = dplyr::all_of(c(pt_scale, bounded)),
                      names_to = "measure", values_to = "value") |>
  dplyr::group_by(measure, genotype) |>
  dplyr::summarise(n_replicate = sum(is.finite(value)),
                   mean = mean(value, na.rm = TRUE),
                   sd = stats::sd(value, na.rm = TRUE), .groups = "drop")

## ---- per-replicate decomposition table ----------------------------------
decomp_tbl <- decomp$by_replicate |>
  dplyr::left_join(tibble::rownames_to_column(replicate_meta, "replicate"),
                   by = "replicate") |>
  dplyr::arrange(match(replicate, rep_ids)) |>
  dplyr::select(replicate, genotype, batch, observed, composition_cf, progression_cf)

## ---- tradeSeq differential-dynamics (GATED; default OFF) ----------------
differential_dynamics <- NULL
if (run_tradeseq) {
  say("[5b/6] tradeSeq differential-dynamics (clean lineage; conditions=genotype)...")
  differential_dynamics <- tryCatch({
    on <- is.finite(pt_clean)
    panel_syms <- unique(unlist(canonical_microglia_markers))
    panel_ens  <- symbols_to_ensembl(panel_syms, symbol_map)
    hvg <- tryCatch(head(VariableFeatures(micro), 150L), error = function(e) character(0))
    panel <- unique(c(as.character(panel_ens), hvg))
    panel <- panel[panel %in% rownames(raw_counts)]
    say("  fitGAM on %d genes x %d on-lineage cells (nknots=6)...",
        length(panel), sum(on))
    sce <- tradeSeq::fitGAM(
      counts = as.matrix(raw_counts[panel, on, drop = FALSE]),
      pseudotime = matrix(pt_clean[on], ncol = 1),
      cellWeights = matrix(1, nrow = sum(on), ncol = 1),
      conditions = factor(meta$genotype[on], levels = genotype_levels),
      nknots = 6L, verbose = FALSE)
    ct <- tradeSeq::conditionTest(sce, global = TRUE, pairwise = FALSE)
    ct$gene <- rownames(ct); ct$symbol <- symbol_map$symbol[match(ct$gene, symbol_map$ensembl)]
    ct$adj.P.Val <- stats::p.adjust(ct$pvalue, "BH")
    tibble::as_tibble(ct)
  }, error = function(e) { say("  tradeSeq FAILED: %s", conditionMessage(e)); NULL })
} else {
  say("[5b/6] tradeSeq SKIPPED (default; pass --tradeseq to enable).")
}

## ---- write results ------------------------------------------------------
say("[6/6] writing results + cache...")
wtsv(by_genotype,    "trajectory_pseudotime_by_genotype.tsv")
wtsv(contrasts_tbl,  "trajectory_contrasts.tsv")
wtsv(decomp_tbl,     "trajectory_progression_decomposition.tsv")
wtsv(concordance$pairwise, "trajectory_method_concordance.tsv")

## ---- assemble cache for rmd/22 ------------------------------------------
py_prov <- tryCatch(jsonlite::fromJSON(file.path(export_dir, "python_provenance.json")),
                    error = function(e) NULL)
trajectory <- list(
  traj_clean = traj_clean,            # headline homeostatic->DAM Slingshot run
  traj_all   = traj_all,              # all-states topology run
  per_cell = tibble::tibble(
    cell_id = cells, state = states, genotype = as.character(meta$genotype),
    batch = as.character(meta$batch), genotype_batch = gb,
    umap_1 = umap[, 1], umap_2 = umap[, 2],
    pt_clean = pt_clean, pt_all = pt_all, dpt = dpt,
    cellrank_dam_fate = cellrank_dam_fate,
    potency_entropy = pot_entropy[cells], potency_n_genes = pot_ngenes[cells],
    potency_cytotrace2 = if (!is.null(pot_cyto)) pot_cyto[cells] else NA_real_),
  potency = potencies,
  root_validation = root_validation,
  per_replicate = list(clean = pr_clean, all = pr_all),
  decomposition = decomp,
  summary_matrix = S,
  contrasts = contrasts_tbl,
  concordance = concordance,
  by_genotype = by_genotype,
  differential_dynamics = differential_dynamics,
  python_provenance = py_prov,
  replicate_meta = replicate_meta,
  params = list(
    n_dims = N_DIMS, n_neighbors = N_NEIGHBORS, n_dcs = N_DCS,
    padj_cut = PADJ_CUT, seed = SEED, headline = "clean_homeostatic_DAM",
    omit_clusters_clean = c("IFN", "proliferative"),
    dam_threshold = pr_clean$dam_threshold,
    methods = c("slingshot", "paga_dpt", "cellrank2", "cytotrace2", "tradeSeq"),
    tradeseq_run = run_tradeseq, velocity = FALSE,
    design = "~ 0 + genotype + batch", contrasts = names(contrast_definitions)))
saveRDS(trajectory, out_cache, compress = "xz")
Sys.chmod(out_cache, "0644")
say("[build_trajectory] wrote cache: %s", out_cache)

## ---- validation summary -------------------------------------------------
n_int <- contrasts_tbl |>
  dplyr::filter(measure == "mean_pt") |> nrow()
stopifnot(
  all(is.finite(pt_clean[is.finite(pt_clean)])),
  sm_clean["DAM"] > sm_clean["homeostatic"],          # monotone vs state means
  nrow(pr_clean$by_replicate) == 16L,                 # all 16 replicates
  n_int == 5L,                                        # 5 contrasts for mean_pt
  all(c("logFC", "t", "P.Value", "adj.P.Val") %in% colnames(contrasts_tbl)))
say("[validation] pt monotone homeostatic<DAM | 16 replicates | mean_pt has 5 contrasts | cols OK")
say("[build_trajectory] done.")

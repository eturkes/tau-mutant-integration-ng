# Cell-cell communication (CCC) inference helpers used by the Phase G
# mechanism layer (see storage/notes/mechanism_layer_plan.md, sessions
# G2..G5). Wraps the Python LIANA+ framework (Dimitrov *et al.* 2024,
# Nat Cell Biol) via reticulate against a project-local virtual env at
# `.venv/bin/python`. The design mirrors `R/kinase_inference.R` /
# `R/tf_inference.R` so the Phase G outputs read cleanly against the
# Phase E + F outputs at the cross-tool verdict layer in section 16.3
# (vs sections 14.3 and 15.3):
#   - per-(sender, receiver, ligand, receptor) tidy tibbles;
#   - canonical five-contrast set (tau_alone, nlgf_in_maptki,
#     nlgf_in_p301s, tau_in_nlgf, interaction);
#   - never pre-commits to a single CCC method (liana itself aggregates
#     9 methods internally via RobustRankAggregate; G3 then aggregates
#     liana with CellChat and MultiNicheNet);
#   - returns long tidy tibbles with the canonical cross-tool key set
#     (`sender`, `receiver`, `ligand`, `receptor`, `lr_interaction`,
#     `prioritization_score`, `n_methods_pass_significance`, `sign_dir`)
#     so the G3 cross-tool rbind is a one-liner.
#
# Why Python LIANA+ over the R `liana` package. The G2 sub-decision-gate
# resolved against the R port: R `liana` was deprecated in 2024 in
# favour of the Python `liana+` rewrite, which is the SOTA in this space
# (multi-method consensus, mouse-native resources, downstream NMF / MOFA
# integrations not present in the R version). The cost is a reticulate
# bridge against a project-local `.venv`; the design here isolates all
# Python-side state to a small set of helpers and treats the venv as a
# private implementation detail of this module.
#
# Why per-genotype, not per-cell. liana's `rank_aggregate()` does not
# do contrasts internally; it produces a per-condition prioritisation
# of LR pairs. To match the project's 2x2 factorial design and the
# five canonical contrasts used everywhere else in the project (see
# `scripts/build_multinichenet.R`), we fit liana four times -- once per
# genotype -- then derive the five contrasts as deltas of the
# per-genotype `magnitude_rank` (lower-is-better in liana, so we flip
# to `comm_strength = 1 - magnitude_rank` first). This is the same
# pattern CellChat uses (per-condition fits then comparative analysis
# via mergeCellChat / rankNet) and keeps liana's aggregation logic
# untouched within each condition.
#
# scipy / S4 bridge gotcha. reticulate's S4 dispatch does NOT expose
# scipy.sparse method calls cleanly: expressions like `X_csc$T$tocsr()`
# error with "$ operator not defined for this S4 class". The workaround
# adopted in `seurat_subset_to_anndata()` is to pass only raw numeric
# arrays (CSC `data`, `indices`, `indptr`) into Python and run all
# scipy / AnnData construction Python-side via `py_run_string()`. This
# sidesteps the whole class of S4 failures.
#
# Cache shape consumed (per G2): the project-wide Seurat object
#   storage/cache/seurat_full_processed.rds
# with `$genotype` (factor with levels MAPTKI, P301S, NLGF_MAPTKI,
# NLGF_P301S), `$cell_type` (factor with the 9 project cell types
# Astrocyte / Microglia_{DAM,IFN,homeostatic,proliferative} / Neuronal /
# OPC / Oligodendrocyte / Vascular), and an `RNA` assay with a `data`
# layer carrying log-normalised expression (gene symbols on rows, cells
# on columns; dgCMatrix). No Ensembl translation step is needed because
# the project's snRNAseq Seurat is symbol-keyed end-to-end (mirror of
# F2's symbol-keyed kinase-substrate network) and the `mouseconsensus`
# liana resource is symbol-keyed mouse-native.
#
# Cache shape produced: two named lists.
#
#   per_genotype: named list (length 4, keys MAPTKI / P301S /
#   NLGF_MAPTKI / NLGF_P301S) of tibbles carrying the full liana_res
#   schema (13 cols) PLUS `cellphone_padj` (BH-corrected
#   `cellphone_pvals` within each per-genotype run):
#     source, target,
#     ligand_complex, receptor_complex,
#     lr_means, lr_logfc, expr_prod,
#     scaled_weight, spec_weight,
#     cellphone_pvals, cellphone_padj,
#     lrscore, magnitude_rank, specificity_rank
#   Stored for diagnostic / reproducibility / downstream-flexibility
#   reasons; G3 reads only the columns it needs.
#
#   per_contrast: named list (length 5, keys = the five canonical
#   contrast names) of tibbles with exactly the cross-tool aggregation
#   columns required by the G2 contract:
#     sender, receiver,
#     ligand, receptor, lr_interaction,
#     prioritization_score (signed delta of comm_strength,
#                           in approx [-1, 1]),
#     n_methods_pass_significance (0/1: liana cellphone_padj < alpha
#                                  in the contrast's primary genotype),
#     sign_dir (sign of prioritization_score: -1, 0, or +1)
#   Plus auxiliary diagnostic cols (mag_rank_primary, mag_rank_reference,
#   cellphone_padj_primary, cellphone_padj_reference) for downstream
#   sanity-checking and to support G3 visualisations.
#
# Missing-pair imputation. liana excludes LR pairs that fail the
# `expr_prop` (>= 10% of cells in either sender or receiver group) and
# `min_cells` (>= 5 cells per group) gates. A pair present in one
# genotype's run but absent in another is semantically "below detection
# in the other genotype", so missing values after the outer-join are
# imputed:
#   magnitude_rank   -> 1.0 (worst rank, lowest priority)
#   cellphone_padj   -> 1.0 (worst adj p, not significant)
# This produces well-defined per-contrast tibbles with zero NAs in the
# prioritization / significance columns, as required by the G2
# verification gate.
#
# n_methods_pass_significance semantics. Liana internally aggregates 9
# methods via RobustRankAggregate within each genotype's run; at the
# cross-tool level (G3) we treat liana as ONE tool with one verdict per
# pair per contrast, not nine. The verdict is "the per-pair permutation
# test (cellphone_pvals, BH-corrected to cellphone_padj within the
# primary genotype) is < alpha". This is the most honest framing for
# cross-tool aggregation: it does not double-count liana's internal
# multi-method evidence as if it were independent of the other tools.
# alpha default = 0.10 to match the cross-tool sig threshold in G3 / F3.
#
# Why magnitude_rank, not specificity_rank or lrscore. magnitude_rank
# is liana's primary aggregated prioritisation -- the RobustRankAggregate
# combination of magnitude-based scores (lr_means, expr_prod, lrscore,
# etc.). It is the recommended summary score per the liana+ paper for
# "which LR pairs are most strongly communicating in this condition".
# specificity_rank is conceptually orthogonal (how cell-type-specific
# is this interaction); we retain it in per_genotype for downstream G4
# use but do NOT use it for the per-contrast prioritization_score.

# === Constants ===

# Project-local Python virtual env. Relative path; build_liana_cache.R
# runs from project root via setwd().
.LIANA_VENV <- ".venv/bin/python"

# Project-wide canonical genotype levels. Order matches the
# group_levels vector in build_multinichenet.R.
.LIANA_GENOTYPES <- c("MAPTKI", "P301S", "NLGF_MAPTKI", "NLGF_P301S")

# Project-wide canonical cell types (9). NULL means "use whatever
# levels are present in the Seurat subset", which is the production
# default; explicit subsets are only used in unit tests / smoke tests.
.LIANA_CELL_TYPES <- c("Astrocyte",
                       "Microglia_DAM",
                       "Microglia_IFN",
                       "Microglia_homeostatic",
                       "Microglia_proliferative",
                       "Neuronal",
                       "OPC",
                       "Oligodendrocyte",
                       "Vascular")

# Five canonical contrasts. Schema:
#   $kind          one of "simple" (primary - reference) or
#                  "interaction" (double difference).
#   $primary       (simple only) genotype on the "+" side.
#   $reference     (simple only) genotype on the "-" side.
#   $primary_pair  (interaction only) c(plus, minus) for the
#                  primary subcontrast.
#   $reference_pair(interaction only) c(plus, minus) for the
#                  reference subcontrast.
#   $sig_primary   (interaction only) which of the four genotypes
#                  carries the "this LR pair is significant" verdict
#                  for the per-contrast significance column. Set to
#                  NLGF_P301S to mirror the build_multinichenet.R
#                  contrast_tbl$group convention (the divergent
#                  condition under the 2x2 design).
#
# For simple contrasts, the primary genotype also carries the
# significance verdict (sig_primary = primary).
.LIANA_CONTRAST_SPECS <- list(
  tau_alone = list(
    kind = "simple",
    primary = "P301S", reference = "MAPTKI"
  ),
  nlgf_in_maptki = list(
    kind = "simple",
    primary = "NLGF_MAPTKI", reference = "MAPTKI"
  ),
  nlgf_in_p301s = list(
    kind = "simple",
    primary = "NLGF_P301S", reference = "P301S"
  ),
  tau_in_nlgf = list(
    kind = "simple",
    primary = "NLGF_P301S", reference = "NLGF_MAPTKI"
  ),
  interaction = list(
    kind = "interaction",
    primary_pair   = c("NLGF_P301S",  "P301S"),
    reference_pair = c("NLGF_MAPTKI", "MAPTKI"),
    sig_primary    = "NLGF_P301S"
  )
)

# === Reticulate bootstrap ===

# Idempotent setup of the project-local Python venv for liana+ work.
# Safe to call repeatedly; no-op after the first successful call in a
# given R session. Stops with an informative error if liana / anndata /
# scipy.sparse are not importable from the venv.
#
# Symlink note: a Python venv's `bin/python` is a symlink chain to the
# system Python (e.g. `python -> python3 -> /usr/bin/python3`); the
# venv context (site-packages, isolation) is determined by the
# `pyvenv.cfg` file that lives next to the VENV `python` symlink, not
# next to the system `python` it points to. `normalizePath()` follows
# symlinks and so destroys the venv context -- a critical pitfall.
# We construct the absolute path manually (cwd + relative path) so
# the venv symlink is preserved.
#
# Initialization note: Python can only be initialized once per R
# session. If reticulate is already bound to a different Python (e.g.
# because some other package triggered initialization first), we stop
# with an actionable error -- restart R and call this helper before
# any other reticulate operation.
#
# Arguments:
#   venv  path to a python binary, relative to the project root or
#         absolute. Default `.venv/bin/python`.
#
# Returns invisibly NULL.
.ensure_liana_python <- function(venv = .LIANA_VENV) {
  py_path <- if (substr(venv, 1L, 1L) == "/") venv else file.path(getwd(), venv)
  if (!file.exists(py_path)) {
    stop(sprintf("liana python venv binary not found at %s (cwd: %s)",
                 py_path, getwd()))
  }
  Sys.setenv(RETICULATE_PYTHON = py_path)
  reticulate::use_python(py_path, required = TRUE)
  if (reticulate::py_available(initialize = TRUE)) {
    cur_py <- reticulate::py_config()$python
    if (!identical(cur_py, py_path)) {
      stop(sprintf(
        "reticulate is bound to %s but liana needs %s. Restart the R session and call .ensure_liana_python() before any other reticulate operation.",
        cur_py, py_path))
    }
  }
  tryCatch({
    reticulate::import("liana")
    reticulate::import("anndata")
    reticulate::import("scipy.sparse")
  }, error = function(e) {
    stop(sprintf(
      "liana python venv at %s is broken (cannot import liana/anndata/scipy.sparse): %s",
      py_path, conditionMessage(e)))
  })
  invisible(NULL)
}

# === Seurat -> AnnData bridge ===

# Convert a Seurat object (or subset) into an AnnData reference held
# in Python via reticulate. The AnnData is bound to `adata` in the
# Python global namespace; this function also returns the reticulate
# reference for convenience.
#
# All scipy.sparse / AnnData construction happens Python-side via
# `py_run_string()` because reticulate's S4 dispatch doesn't expose
# scipy methods (see header note on the scipy / S4 bridge gotcha).
# We pass only raw arrays (CSC data / indices / indptr from the
# dgCMatrix slots) into Python.
#
# The Seurat dgCMatrix is genes x cells (column-compressed); AnnData
# expects cells x genes. We construct a CSC matrix on the Python side
# matching the dgCMatrix layout, then transpose-and-tocsr.
#
# Arguments:
#   sub       Seurat object or subset.
#   assay     assay name, default "RNA".
#   layer     layer name, default "data" (log-normalised counts).
#   obs_cols  character vector of meta.data columns to carry into
#             `adata.obs`. Factors are coerced to character to avoid
#             pandas categorical-dtype surprises on the liana side.
#
# Returns the `adata` reticulate Python reference.
seurat_subset_to_anndata <- function(sub,
                                     assay = "RNA",
                                     layer = "data",
                                     obs_cols = c("cell_type", "genotype")) {
  .ensure_liana_python()
  mat <- SeuratObject::GetAssayData(sub, assay = assay, layer = layer)
  stopifnot(inherits(mat, "dgCMatrix"))

  py <- reticulate::py
  py$dat     <- mat@x
  py$idx     <- as.integer(mat@i)
  py$ptr     <- as.integer(mat@p)
  py$n_genes <- as.integer(nrow(mat))
  py$n_cells <- as.integer(ncol(mat))
  reticulate::py_run_string("
from scipy.sparse import csc_matrix
X_csc = csc_matrix((dat, idx, ptr), shape=(n_genes, n_cells))
X_csr = X_csc.T.tocsr()
")

  meta <- sub@meta.data
  obs_cols <- intersect(obs_cols, colnames(meta))
  obs_df <- as.data.frame(meta[, obs_cols, drop = FALSE],
                          row.names = colnames(sub),
                          stringsAsFactors = FALSE)
  for (cn in colnames(obs_df)) {
    if (is.factor(obs_df[[cn]])) obs_df[[cn]] <- as.character(obs_df[[cn]])
  }
  var_df <- data.frame(row.names = rownames(sub), stringsAsFactors = FALSE)

  py$obs_py        <- reticulate::r_to_py(obs_df)
  py$var_py        <- reticulate::r_to_py(var_df)
  py$gene_names_py <- reticulate::r_to_py(rownames(sub))
  py$cell_names_py <- reticulate::r_to_py(colnames(sub))
  reticulate::py_run_string("
import anndata as ad
adata = ad.AnnData(X = X_csr, obs = obs_py, var = var_py)
adata.var.index = list(gene_names_py)
adata.obs.index = list(cell_names_py)
")
  reticulate::py$adata
}

# === Per-genotype liana run ===

# Run liana.mt.rank_aggregate on a per-genotype subset of a Seurat
# object. Returns a tibble of the full liana_res schema PLUS a
# BH-corrected `cellphone_padj` column.
#
# Arguments:
#   sc            full Seurat object (carries $genotype and $cell_type).
#   genotype      character, one of the levels of sc$genotype.
#   cell_types    character vector, subset of cell types to include.
#                 Default NULL means "all cell types present in the
#                 genotype subset after droplevels".
#   resource_name liana resource name. Default "mouseconsensus"
#                 (3989 mouse-native LR pairs).
#   use_raw       FALSE because the Seurat layer "data" is already
#                 log-normalised.
#   expr_prop     min fraction of cells in either sender or receiver
#                 group expressing the L or R for the pair to be kept.
#                 Default 0.10 (liana paper default).
#   min_cells     min cells per cell-type group. Default 5L (matches
#                 build_multinichenet.R::min_cells_per_celltype = 5).
#   n_perms       permutations for cellphone_pvals. Default 100L to
#                 match CellChat's nboot_perm = 100L in
#                 build_cellchat.R. Floor on p-value = 1/n_perms.
#   seed          numeric seed for the permutation RNG.
#   verbose       logical; if TRUE, prints liana's per-step messages.
#
# Returns tibble with 14 cols (13 liana_res + cellphone_padj).
run_liana_rank_aggregate_for_genotype <- function(
    sc,
    genotype,
    cell_types    = NULL,
    resource_name = "mouseconsensus",
    use_raw       = FALSE,
    expr_prop     = 0.10,
    min_cells     = 5L,
    n_perms       = 100L,
    seed          = 1L,
    verbose       = FALSE) {
  .ensure_liana_python()
  stopifnot(genotype %in% as.character(sc$genotype))

  keep <- (as.character(sc$genotype) == genotype)
  if (!is.null(cell_types)) {
    keep <- keep & (as.character(sc$cell_type) %in% cell_types)
  }
  sub <- sc[, which(keep)]
  if (is.factor(sub$cell_type))  sub$cell_type  <- droplevels(sub$cell_type)
  if (is.factor(sub$genotype))   sub$genotype   <- droplevels(sub$genotype)
  if (verbose) {
    message(sprintf("  [%s] %d cells, %d cell types: %s",
                    genotype, ncol(sub),
                    length(unique(sub$cell_type)),
                    paste(sort(unique(as.character(sub$cell_type))), collapse = ",")))
  }

  seurat_subset_to_anndata(sub, obs_cols = c("cell_type", "genotype"))

  # `reticulate::py$x <- v` is invalid R syntax (R cannot assign to a
  # namespaced-name slot expression); bind locally first.
  py <- reticulate::py
  py$resource_name <- resource_name
  py$use_raw       <- use_raw
  py$expr_prop     <- expr_prop
  py$min_cells     <- as.integer(min_cells)
  py$n_perms       <- as.integer(n_perms)
  py$seed          <- as.integer(seed)
  py$verbose_py    <- verbose
  reticulate::py_run_string("
import liana.mt as li_mt
li_mt.rank_aggregate(
    adata,
    groupby='cell_type',
    resource_name=resource_name,
    use_raw=use_raw,
    expr_prop=expr_prop,
    min_cells=min_cells,
    n_perms=n_perms,
    seed=seed,
    verbose=verbose_py,
)
res = adata.uns['liana_res']
")
  res <- tibble::as_tibble(reticulate::py$res)

  # BH within this per-genotype run -- cross-tool sig is per-genotype,
  # not per-(genotype, contrast), so this is the right scope.
  res$cellphone_padj <- p.adjust(res$cellphone_pvals, method = "BH")
  res
}

# === Driver: per-genotype list ===

# Iterate over genotypes, returning a named list of liana_res tibbles
# (one entry per genotype). All Python state (the `adata` global) is
# overwritten on each iteration; only the R-side tibbles are retained.
#
# Arguments:
#   sc           full Seurat object.
#   genotypes    character vector. Default `.LIANA_GENOTYPES` (all 4).
#   cell_types   character vector. Default `.LIANA_CELL_TYPES` (all 9).
#                Pass NULL to use "whatever is in each genotype subset".
#   verbose      logical; emits per-genotype timing + size messages.
#   ...          forwarded to `run_liana_rank_aggregate_for_genotype()`.
#
# Returns named list of tibbles, length == length(genotypes).
build_liana_per_genotype_list <- function(
    sc,
    genotypes  = .LIANA_GENOTYPES,
    cell_types = .LIANA_CELL_TYPES,
    verbose    = TRUE,
    ...) {
  out <- vector("list", length(genotypes))
  names(out) <- genotypes
  for (g in genotypes) {
    if (verbose) {
      message(sprintf("[liana] rank_aggregate for genotype: %s", g))
    }
    t0 <- proc.time()
    out[[g]] <- run_liana_rank_aggregate_for_genotype(
      sc = sc, genotype = g, cell_types = cell_types,
      verbose = verbose, ...
    )
    t1 <- proc.time()
    if (verbose) {
      message(sprintf("    -> %d LR-pair rows in %.1fs (cellphone_padj<0.10: %d)",
                      nrow(out[[g]]), (t1 - t0)[3],
                      sum(out[[g]]$cellphone_padj < 0.10, na.rm = TRUE)))
    }
  }
  out
}

# === Per-contrast derivation ===

# Compute the five canonical per-contrast tibbles from a per-genotype
# list. See the module header for the schema and missing-pair
# imputation rules.
#
# Arguments:
#   per_genotype     named list of liana_res tibbles, one per genotype,
#                    as produced by `build_liana_per_genotype_list()`.
#   contrast_specs   list of contrast specifications. Default
#                    `.LIANA_CONTRAST_SPECS` (5 canonical contrasts).
#   sig_alpha        significance threshold for cellphone_padj.
#                    Default 0.10 (project-wide cross-tool / cross-
#                    modality alpha, matches F3 / G3 conventions).
#
# Returns named list of tibbles, one per contrast name.
derive_liana_per_contrast_list <- function(
    per_genotype,
    contrast_specs = .LIANA_CONTRAST_SPECS,
    sig_alpha      = 0.10) {
  stopifnot(is.list(per_genotype))
  stopifnot(all(c("source", "target", "ligand_complex", "receptor_complex",
                  "magnitude_rank", "cellphone_padj") %in%
                colnames(per_genotype[[1]])))

  # Reduce each per-genotype tibble to the per-pair key + the two
  # columns we need downstream. This trims memory before the joins.
  to_pair_tbl <- function(df) {
    tibble::tibble(
      sender    = df$source,
      receiver  = df$target,
      ligand    = df$ligand_complex,
      receptor  = df$receptor_complex,
      mag_rank  = df$magnitude_rank,
      cp_padj   = df$cellphone_padj
    )
  }
  keyed <- lapply(per_genotype, to_pair_tbl)
  key_cols <- c("sender", "receiver", "ligand", "receptor")

  # Helper: outer-join two per-genotype pair tibbles, impute NAs.
  join2 <- function(A, B, suf_a, suf_b) {
    A2 <- dplyr::rename(A,
                        !!paste0("mag_rank", suf_a) := mag_rank,
                        !!paste0("cp_padj",  suf_a) := cp_padj)
    B2 <- dplyr::rename(B,
                        !!paste0("mag_rank", suf_b) := mag_rank,
                        !!paste0("cp_padj",  suf_b) := cp_padj)
    dplyr::full_join(A2, B2, by = key_cols)
  }
  impute_worst <- function(x) {
    x[is.na(x)] <- 1.0
    x
  }

  compute_simple <- function(spec) {
    j <- join2(keyed[[spec$primary]], keyed[[spec$reference]],
               suf_a = "_p", suf_b = "_r")
    j$mag_rank_p <- impute_worst(j$mag_rank_p)
    j$mag_rank_r <- impute_worst(j$mag_rank_r)
    j$cp_padj_p  <- impute_worst(j$cp_padj_p)
    j$cp_padj_r  <- impute_worst(j$cp_padj_r)

    prio <- (1 - j$mag_rank_p) - (1 - j$mag_rank_r)

    tibble::tibble(
      sender                      = j$sender,
      receiver                    = j$receiver,
      ligand                      = j$ligand,
      receptor                    = j$receptor,
      lr_interaction              = paste(j$ligand, j$receptor, sep = "_"),
      prioritization_score        = prio,
      n_methods_pass_significance = as.integer(j$cp_padj_p < sig_alpha),
      sign_dir                    = sign(prio),
      mag_rank_primary            = j$mag_rank_p,
      mag_rank_reference          = j$mag_rank_r,
      cellphone_padj_primary      = j$cp_padj_p,
      cellphone_padj_reference    = j$cp_padj_r
    )
  }

  compute_interaction <- function(spec) {
    # Four-way outer join in two steps.
    P1 <- keyed[[spec$primary_pair[1]]]    # NLGF_P301S
    P0 <- keyed[[spec$primary_pair[2]]]    # P301S
    R1 <- keyed[[spec$reference_pair[1]]]  # NLGF_MAPTKI
    R0 <- keyed[[spec$reference_pair[2]]]  # MAPTKI

    j  <- join2(P1, P0, suf_a = "_p1", suf_b = "_p0")
    j2 <- join2(R1, R0, suf_a = "_r1", suf_b = "_r0")
    j  <- dplyr::full_join(j, j2, by = key_cols)

    for (col in c("mag_rank_p1", "mag_rank_p0",
                  "mag_rank_r1", "mag_rank_r0",
                  "cp_padj_p1",  "cp_padj_p0",
                  "cp_padj_r1",  "cp_padj_r0")) {
      j[[col]] <- impute_worst(j[[col]])
    }

    delta_p <- (1 - j$mag_rank_p1) - (1 - j$mag_rank_p0)
    delta_r <- (1 - j$mag_rank_r1) - (1 - j$mag_rank_r0)
    prio    <- delta_p - delta_r

    sig_col <- switch(spec$sig_primary,
      "NLGF_P301S"  = j$cp_padj_p1,
      "P301S"       = j$cp_padj_p0,
      "NLGF_MAPTKI" = j$cp_padj_r1,
      "MAPTKI"      = j$cp_padj_r0,
      stop("Unknown sig_primary for interaction contrast: ", spec$sig_primary)
    )

    tibble::tibble(
      sender                      = j$sender,
      receiver                    = j$receiver,
      ligand                      = j$ligand,
      receptor                    = j$receptor,
      lr_interaction              = paste(j$ligand, j$receptor, sep = "_"),
      prioritization_score        = prio,
      n_methods_pass_significance = as.integer(sig_col < sig_alpha),
      sign_dir                    = sign(prio),
      delta_in_primary_pair       = delta_p,
      delta_in_reference_pair     = delta_r,
      cellphone_padj_sig_primary  = sig_col
    )
  }

  out <- vector("list", length(contrast_specs))
  names(out) <- names(contrast_specs)
  for (cn in names(contrast_specs)) {
    spec <- contrast_specs[[cn]]
    out[[cn]] <- switch(spec$kind,
      "simple"      = compute_simple(spec),
      "interaction" = compute_interaction(spec),
      stop("Unknown contrast kind: ", spec$kind)
    )
  }
  out
}

# === Provenance metadata ===

# Build a small list documenting the liana run for reproducibility.
# Stored alongside per_genotype / per_contrast in the cache rds.
#
# Arguments:
#   per_genotype  list returned by `build_liana_per_genotype_list()`.
#   params        list of run parameters (resource_name, use_raw,
#                 expr_prop, min_cells, n_perms, seed, sig_alpha,
#                 cell_types).
#
# Returns a list with named slots: liana_version, anndata_version,
# scipy_version, timestamp, params, n_pairs_per_genotype, n_cells_per_genotype.
build_liana_provenance <- function(per_genotype, params) {
  .ensure_liana_python()
  li      <- reticulate::import("liana")
  ad      <- reticulate::import("anndata")
  sp      <- reticulate::import("scipy")
  py_ver  <- reticulate::py_config()$version

  n_pairs <- sapply(per_genotype, nrow)
  list(
    liana_version    = li[["__version__"]],
    anndata_version  = ad[["__version__"]],
    scipy_version    = sp[["__version__"]],
    python_version   = py_ver,
    venv             = .LIANA_VENV,
    timestamp        = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    params           = params,
    n_pairs_per_genotype = n_pairs
  )
}

# ====================================================================
# Phase G3 helpers: cross-tool LR-pair leader board.
# ====================================================================
#
# G3 takes the three CCC tool caches (CellChat per_condition,
# MultiNicheNet group_prioritization_tbl, liana per_contrast from G2)
# and produces a unified cross-tool leader board mirroring the E3 TF
# cross-modality leader board: per-(contrast, LR-pair) rows with
# n_tools_sig / n_tools_consistent_sign / sign_consensus / mean
# normalised score columns, then a per-LR-pair aggregate across the
# five canonical contrasts (leader_score formula identical to E3).
#
# Per-tool significance conventions (a 2026-05-24 G3 pre-build diagnostic
# established each tool's natural scale):
#   - CellChat: per-condition `subsetCommunication()` returns ONLY
#     pairs at pval < 0.05 (the function's own `thresh` default), so
#     CellChat's "sig" vote at any FDR<=0.10 cross-tool gate is the
#     binary "this LR pair exists in CellChat's output for the
#     contrast's primary genotype" -- a meaningful inclusion vote
#     because CellChat already applies its own multiple-testing
#     control through the permutation null at build time.
#   - MultiNicheNet: scaled_p_val_ligand_adapted < alpha AND
#     scaled_p_val_receptor_adapted < alpha. The "_adapted" columns
#     are MNN's BH-adjusted, rank-percentile-scaled DESeq2 p-values
#     for the ligand and receptor genes in the sender / receiver
#     cell types under the contrast; 0 = strongest, 1 = weakest. The
#     AND combination is MNN's recommended "both ends of the L-R
#     pair are DE in the contrast" filter.
#   - liana: cellphone_padj_primary < alpha (BH-adjusted within the
#     per-genotype liana run for simple contrasts; BH within the
#     interaction's primary genotype run for the interaction
#     contrast). Computed at G2 from the per-genotype cellphone_pvals
#     output of liana's rank_aggregate.
#
# Cross-tool sign attribution per-tool:
#   - CellChat: sign(delta) where delta = primary_prob - reference_prob
#     for simple contrasts, double-difference for interaction.
#     subsetCommunication only emits non-zero prob entries; missing
#     pairs are imputed as zero.
#   - MultiNicheNet: sign((scaled_lfc_ligand + scaled_lfc_receptor)/2),
#     i.e. the average direction of regulation across the two ends
#     of the L-R pair.
#   - liana: liana_per_contrast$sign_dir (already signed per the G2
#     contract: sign of the per-contrast prioritization_score, which
#     is a delta of comm strengths).
#
# Receptor heteromer handling. CellChat encodes some receptor complexes
# as underscore-joined symbols (e.g. "NRP1_PLXNA2"). To make the
# cross-tool keys comparable to MultiNicheNet + liana (both of which
# emit single-symbol receptors per row), we EXPLODE such receptors into
# multiple rows during extraction, mirroring the convention used in the
# G3 pre-build diagnostic.
# Ligand complexes in CellChat are exceedingly rare and not currently
# exploded; flag this as a known limitation if such a ligand surfaces
# at the top of a contrast's leader board.
#
# Leader rule (locked at G3 against the diagnostic):
#   per (LR-pair, contrast) cell: leader iff
#     n_tools_sig_consistent_sign >= 2L OR n_tools_sig >= 3L
#   at FDR<0.10 (project-wide cross-tool / cross-modality alpha).
#   This is the OR form mirroring E3 exactly; the AND form
#   `n_tools_sig >= 2 AND n_tools_sig_consistent_sign >= 2` from the
#   G3 plan stub collapses to the same cells in this data because a
#   pair can only be sign-consistent in >= 2 tools if it is also
#   sig in >= 2 tools. The OR form is preferred for E3 parity.
#   Diagnostic counts at FDR<0.10 per contrast:
#     tau_alone      n_2cons_OR_3sig =  70  (n_3t_sig = 0 -- honest
#                                            non-finding for the all-3
#                                            arm; the OR fallback
#                                            covers it)
#     nlgf_in_maptki n_2cons_OR_3sig = 123
#     nlgf_in_p301s  n_2cons_OR_3sig = 128
#     tau_in_nlgf    n_2cons_OR_3sig = 130
#     interaction    n_2cons_OR_3sig = 106
#   Far from the F3 emergency (which yielded 0 cross-cache cells
#   across all five contrasts and forced the user-directed
#   contrast-independent deviation). G3 proceeds with the planned
#   rule.
#
# Per-tool normalised scoring (cross-tool magnitude tiebreak).
# CellChat's delta prob is in approx [-1, +1] (mostly small); MNN's
# prioritization_score is in [0, 1] (signed via lfc separately);
# liana's prioritization_score from G2 is a signed delta in approx
# [-1, +1]. To compare magnitudes across tools we normalise each
# tool's per-contrast distribution to [-1, +1] by dividing by the
# tool's per-contrast max(|score|). The `mean_norm_score` column
# averages these normalised scores over the *significant* per-tool
# cells only (mirror of E3's mean_abs_score over the sig subset).
# Per-tool raw scores are retained alongside the normalised mean for
# transparency.
#
# Schema produced by extract_lr_per_tool() (long form, one row per
# (tool, contrast, sender, receiver, ligand, receptor) cell):
#   tool             "cellchat" / "mnn" / "liana".
#   contrast         one of the 5 canonical contrasts.
#   sender, receiver, ligand, receptor      character keys.
#   lr_interaction   paste0(ligand, "_", receptor) for cross-tool key.
#   score            tool-native magnitude (delta prob for cellchat;
#                    prioritization_score for mnn / liana).
#   padj             tool-native sig statistic at this cell.
#   is_sig           padj < alpha, evaluated at extract time.
#   sign_dir         -1 / 0 / +1.
#
# Schema produced by rank_lr_cross_tool() (one row per
# (contrast, sender, receiver, ligand, receptor) cell):
#   sender, receiver, ligand, receptor      character keys.
#   lr_interaction                          paste0(ligand, "_", receptor).
#   n_tools_present                         integer (1..3).
#   n_tools_sig                             integer (0..3).
#   n_tools_sig_consistent_sign             integer (0..3).
#   sign_consensus                          "+" / "-" / "mixed" / NA.
#   mean_norm_score                         numeric.
#   score_<tool>, padj_<tool>, sig_<tool>   per-tool cells.
#   composite_rank                          integer; sort key.
#
# Schema produced by build_lr_cross_tool_leaderboard() (one row per
# leader LR-pair):
#   lr_interaction, sender, receiver, ligand, receptor      keys.
#   n_contrasts_leader                                      integer.
#   n_contrasts_consistent_sign_ge2                         integer.
#   n_contrasts_sig_ge3                                     integer.
#   max_consistent_sign                                     integer.
#   max_n_tools_sig                                         integer.
#   max_mean_norm_score                                     numeric.
#   dominant_sign                                           "+"/"-"/"mixed"/NA.
#   contrasts_summary                                       pipe-delimited
#                                                           "contrast:n_sig/sign".
#   leader_score                                            composite.
# Sort: leader_score desc, n_contrasts_consistent_sign_ge2 desc,
#       max_consistent_sign desc, max_mean_norm_score desc, lr_interaction asc.

# --------------------------------------------------------------------
# Internal: tool-specific extractors.
# --------------------------------------------------------------------

# CellChat extractor. Takes the per-condition list and the contrast
# specs; for each contrast, derives a per-(sender, receiver, ligand,
# receptor) tibble with delta-prob score (simple) or double-difference
# (interaction), sig boolean (CellChat pval < alpha in the contrast's
# primary genotype), and sign_dir.
#
# Arguments:
#   per_cond         named list of CellChat objects, one per genotype.
#   contrast_specs   named list of contrast specs (see .LIANA_CONTRAST_SPECS).
#   alpha            FDR / pval cut for the per-genotype sig call.
#                    Default 0.10 (project-wide).
#
# Returns long tibble with columns:
#   tool="cellchat", contrast, sender, receiver, ligand, receptor,
#   lr_interaction, score (delta prob), padj (raw pval in primary),
#   is_sig, sign_dir.
.extract_lr_cellchat <- function(per_cond,
                                 contrast_specs = .LIANA_CONTRAST_SPECS,
                                 alpha = 0.10) {
  stopifnot(is.list(per_cond),
            !is.null(names(per_cond)),
            all(c("MAPTKI", "P301S", "NLGF_MAPTKI", "NLGF_P301S")
                %in% names(per_cond)))

  # Per-genotype tidy frame with receptor-complex explosion (same
  # convention as the extraction path above). After
  # explosion, dedupe by (sender, receiver, ligand, receptor): a
  # single LR pair can appear in multiple `subsetCommunication` rows
  # when the same pair is encoded with different `interaction_name_2`
  # downstream-complex annotations (e.g. `Glu-SLC1A2-GLS_GRIK4_GRIA1`
  # and `Glu-SLC1A2-GLS_GRIK4_GRIA2` both reduce to the LR cell
  # `Glu-SLC1A2_GLS -> GRIK4` after the receptor-complex explosion).
  # Without the dedupe, a CellChat-only pair would produce duplicate
  # rows in the long tibble and inflate the cross-tool `n_tools_present`
  # count for that cell. Dedupe rule: keep the row with the strongest
  # CellChat evidence (max prob, lowest pval as tiebreak); this is
  # the CellChat-vignette convention for LR-pair-level aggregation
  # over `interaction_name_2`.
  per_geno <- lapply(per_cond, function(cc) {
    com <- CellChat::subsetCommunication(cc)
    com |>
      dplyr::mutate(
        receptor_parts = strsplit(receptor, "_", fixed = TRUE)
      ) |>
      tidyr::unnest(receptor_parts) |>
      dplyr::transmute(
        sender   = as.character(source),
        receiver = as.character(target),
        ligand   = as.character(ligand),
        receptor = receptor_parts,
        prob     = prob,
        pval     = pval
      ) |>
      dplyr::arrange(dplyr::desc(prob), pval) |>
      dplyr::distinct(sender, receiver, ligand, receptor, .keep_all = TRUE)
  })

  key_cols <- c("sender", "receiver", "ligand", "receptor")
  make_key <- function(df) {
    paste(df$sender, df$receiver, df$ligand, df$receptor, sep = "||")
  }

  out <- vector("list", length(contrast_specs))
  names(out) <- names(contrast_specs)
  for (cn in names(contrast_specs)) {
    spec <- contrast_specs[[cn]]
    prim <- if (spec$kind == "simple") spec$primary else spec$sig_primary
    tbl  <- per_geno[[prim]]
    key_full <- make_key(tbl)
    # Sign / score: delta vs reference (simple) or double-difference
    # (interaction). Imputation rule for missing pairs in the
    # comparator genotypes: prob = 0 (CellChat's subsetCommunication
    # omits pairs with no detected signalling; treat them as zero
    # comm strength, the natural CellChat semantic).
    if (spec$kind == "simple") {
      rtbl <- per_geno[[spec$reference]]
      ref_prob <- rtbl$prob[match(key_full, make_key(rtbl))]
      ref_prob[is.na(ref_prob)] <- 0
      delta <- tbl$prob - ref_prob
    } else {
      pp <- spec$primary_pair    # c(NLGF_P301S, P301S)
      rp <- spec$reference_pair  # c(NLGF_MAPTKI, MAPTKI)
      P0t <- per_geno[[pp[2]]]
      R1t <- per_geno[[rp[1]]]
      R0t <- per_geno[[rp[2]]]
      P0p <- P0t$prob[match(key_full, make_key(P0t))]
      R1p <- R1t$prob[match(key_full, make_key(R1t))]
      R0p <- R0t$prob[match(key_full, make_key(R0t))]
      P0p[is.na(P0p)] <- 0; R1p[is.na(R1p)] <- 0; R0p[is.na(R0p)] <- 0
      delta <- (tbl$prob - P0p) - (R1p - R0p)
    }

    is_sig <- !is.na(tbl$pval) & tbl$pval < alpha
    tibble::tibble(
      tool           = "cellchat",
      contrast       = cn,
      sender         = tbl$sender,
      receiver       = tbl$receiver,
      ligand         = tbl$ligand,
      receptor       = tbl$receptor,
      lr_interaction = paste(tbl$ligand, tbl$receptor, sep = "_"),
      score          = delta,
      padj           = tbl$pval,  # CellChat's pval is already perm-MTC
      is_sig         = is_sig,
      sign_dir       = sign(delta)
    ) -> out[[cn]]
  }
  dplyr::bind_rows(out)
}

# MultiNicheNet extractor. Takes the group_prioritization_tbl and a
# semantic-to-raw contrast map, splits by contrast, and emits the
# canonical tidy schema.
#
# Arguments:
#   mnn_output       the loaded multinichenet_output.rds list.
#   contrast_map     named character vector mapping semantic contrast
#                    names to MNN's raw contrast strings. Default
#                    matches scripts/build_multinichenet.R / rmd/11_ccc.Rmd.
#   alpha            sig cut for both scaled_p_val_*_adapted columns.
#                    Default 0.10.
.extract_lr_multinichenet <- function(mnn_output,
                                      contrast_map = c(
                                        tau_alone       = "P301S-MAPTKI",
                                        nlgf_in_maptki  = "NLGF_MAPTKI-MAPTKI",
                                        nlgf_in_p301s   = "NLGF_P301S-P301S",
                                        tau_in_nlgf     = "NLGF_P301S-NLGF_MAPTKI",
                                        interaction     = "(NLGF_P301S-P301S)-(NLGF_MAPTKI-MAPTKI)"
                                      ),
                                      alpha = 0.10) {
  stopifnot(!is.null(mnn_output$prioritization_tables$group_prioritization_tbl))
  g_tbl <- mnn_output$prioritization_tables$group_prioritization_tbl

  out <- vector("list", length(contrast_map))
  names(out) <- names(contrast_map)
  for (cn in names(contrast_map)) {
    raw <- contrast_map[[cn]]
    sub <- g_tbl[g_tbl$contrast == raw, , drop = FALSE]
    # MNN's adj-p columns are "scaled_p_val_*_adapted" (BH on per-celltype
    # DESeq2 ligand / receptor fits, scaled to a rank percentile in [0, 1];
    # 0 = strongest, 1 = weakest). The AND filter ("both ends of the L-R
    # pair are DE in the contrast") is MNN's recommended sig call.
    pl <- sub[["scaled_p_val_ligand_adapted"]]
    pr <- sub[["scaled_p_val_receptor_adapted"]]
    is_sig <- !is.na(pl) & !is.na(pr) & pl < alpha & pr < alpha
    sign_dir <- sign((sub$scaled_lfc_ligand + sub$scaled_lfc_receptor) / 2)
    sign_dir[is.na(sign_dir)] <- 0L
    out[[cn]] <- tibble::tibble(
      tool           = "mnn",
      contrast       = cn,
      sender         = as.character(sub$sender),
      receiver       = as.character(sub$receiver),
      ligand         = as.character(sub$ligand),
      receptor       = as.character(sub$receptor),
      lr_interaction = paste(as.character(sub$ligand),
                             as.character(sub$receptor), sep = "_"),
      score          = sub$prioritization_score,
      padj           = pmax(pl, pr, na.rm = FALSE),  # worse-of-two summary
      is_sig         = is_sig,
      sign_dir       = sign_dir
    )
  }
  dplyr::bind_rows(out)
}

# Liana extractor. Reads the G2 per_contrast tibbles, which are
# already in canonical schema. Uses cellphone_padj_primary (or
# cellphone_padj_sig_primary for interaction) as the sig statistic.
#
# Arguments:
#   liana_cache      the loaded liana_output.rds list (must contain
#                    `per_contrast`).
#   alpha            sig cut on cellphone_padj. Default 0.10.
.extract_lr_liana <- function(liana_cache, alpha = 0.10) {
  stopifnot(!is.null(liana_cache$per_contrast))
  pc <- liana_cache$per_contrast

  out <- vector("list", length(pc))
  names(out) <- names(pc)
  for (cn in names(pc)) {
    lt <- pc[[cn]]
    padj <- if ("cellphone_padj_primary" %in% colnames(lt)) {
      lt$cellphone_padj_primary
    } else if ("cellphone_padj_sig_primary" %in% colnames(lt)) {
      lt$cellphone_padj_sig_primary
    } else {
      stop("liana per-contrast tibble missing cellphone_padj column for ",
           cn, call. = FALSE)
    }
    is_sig <- !is.na(padj) & padj < alpha
    out[[cn]] <- tibble::tibble(
      tool           = "liana",
      contrast       = cn,
      sender         = lt$sender,
      receiver       = lt$receiver,
      ligand         = lt$ligand,
      receptor       = lt$receptor,
      lr_interaction = lt$lr_interaction,
      score          = lt$prioritization_score,
      padj           = padj,
      is_sig         = is_sig,
      sign_dir       = lt$sign_dir
    )
  }
  dplyr::bind_rows(out)
}

# --------------------------------------------------------------------
# Public: tool extractor dispatcher.
# --------------------------------------------------------------------

# Build the long per-(tool, contrast, lr-pair) tibble by rbinding the
# three tool-specific extractors. The returned tibble is the input to
# `rank_lr_cross_tool()`.
#
# Arguments:
#   cellchat_per_cond   named list of CellChat objects (per genotype).
#   mnn_output          the loaded multinichenet_output.rds list.
#   liana_cache         the loaded liana_output.rds list.
#   alpha               per-tool sig cut. Default 0.10 (project-wide).
#   contrast_specs      contrast specs for CellChat sign / delta
#                       derivation. Default `.LIANA_CONTRAST_SPECS`.
#   mnn_contrast_map    semantic-to-raw map for MNN extraction.
#                       Default matches build_multinichenet.R / rmd/11_ccc.Rmd.
#
# Returns long tibble.
extract_lr_per_tool <- function(cellchat_per_cond,
                                mnn_output,
                                liana_cache,
                                alpha          = 0.10,
                                contrast_specs = .LIANA_CONTRAST_SPECS,
                                mnn_contrast_map = c(
                                  tau_alone       = "P301S-MAPTKI",
                                  nlgf_in_maptki  = "NLGF_MAPTKI-MAPTKI",
                                  nlgf_in_p301s   = "NLGF_P301S-P301S",
                                  tau_in_nlgf     = "NLGF_P301S-NLGF_MAPTKI",
                                  interaction     = "(NLGF_P301S-P301S)-(NLGF_MAPTKI-MAPTKI)"
                                )) {
  cc <- .extract_lr_cellchat(cellchat_per_cond,
                             contrast_specs = contrast_specs,
                             alpha = alpha)
  mn <- .extract_lr_multinichenet(mnn_output,
                                  contrast_map = mnn_contrast_map,
                                  alpha = alpha)
  li <- .extract_lr_liana(liana_cache, alpha = alpha)
  dplyr::bind_rows(cc, mn, li)
}

# --------------------------------------------------------------------
# Public: cross-tool per-contrast ranking.
# --------------------------------------------------------------------

# For a given contrast, build the per-(LR-pair) cross-tool tibble:
# one row per unique (sender, receiver, ligand, receptor) cell across
# the three tools, with cross-tool sig / sign / mean-norm-score
# columns plus per-tool detail cells.
#
# Score normalisation: within each (tool, contrast), divide raw scores
# by max(|score|) in that tool's per-contrast distribution so the
# magnitudes are comparable on [-1, +1]. The `mean_norm_score` column
# is the mean of normalised scores over the per-pair SIGNIFICANT cells
# only (mirror of E3's mean_abs_score over the sig subset); if no cells
# are sig at the pair, mean_norm_score falls back to the mean over
# all tools that have the pair present.
#
# Arguments:
#   lr_long              long tibble from `extract_lr_per_tool()`.
#   contrast             single contrast name.
#   tools                character vector of tool keys to include.
#                        Default `c("cellchat", "mnn", "liana")`.
#
# Returns tibble sorted by composite rank.
rank_lr_cross_tool <- function(lr_long, contrast,
                               tools = c("cellchat", "mnn", "liana")) {
  stopifnot(is.data.frame(lr_long),
            length(contrast) == 1L,
            length(tools) > 0L,
            all(tools %in% unique(lr_long$tool)))
  cn <- contrast
  long <- lr_long |>
    dplyr::filter(.data$contrast == cn, .data$tool %in% tools)

  if (nrow(long) == 0L) {
    return(tibble::tibble(
      sender = character(0), receiver = character(0),
      ligand = character(0), receptor = character(0),
      lr_interaction = character(0),
      n_tools_present = integer(0),
      n_tools_sig = integer(0),
      n_tools_sig_consistent_sign = integer(0),
      sign_consensus = character(0),
      mean_norm_score = numeric(0),
      composite_rank = integer(0)
    ))
  }

  # Per-tool magnitude normalisation within this contrast.
  long <- long |>
    dplyr::group_by(tool) |>
    dplyr::mutate(
      .max_abs = suppressWarnings(max(abs(score), na.rm = TRUE)),
      score_norm = ifelse(is.finite(.max_abs) & .max_abs > 0,
                          score / .max_abs, NA_real_)
    ) |>
    dplyr::ungroup() |>
    dplyr::select(-.max_abs)

  key_cols <- c("sender", "receiver", "ligand", "receptor")

  # Per-(LR-pair) aggregation across tools.
  sig_signed <- ifelse(long$is_sig, long$sign_dir, 0)
  wide <- long |>
    dplyr::mutate(.sig_signed = sig_signed) |>
    dplyr::group_by(sender, receiver, ligand, receptor) |>
    dplyr::summarise(
      lr_interaction              = dplyr::first(lr_interaction),
      n_tools_present             = dplyr::n(),
      n_tools_sig                 = sum(is_sig, na.rm = TRUE),
      n_pos_sig                   = sum(.sig_signed > 0, na.rm = TRUE),
      n_neg_sig                   = sum(.sig_signed < 0, na.rm = TRUE),
      mean_norm_score_sig         = {
        v <- score_norm[is_sig]
        v <- v[is.finite(v)]
        if (length(v) == 0L) NA_real_ else mean(v)
      },
      mean_norm_score_present     = {
        v <- score_norm
        v <- v[is.finite(v)]
        if (length(v) == 0L) NA_real_ else mean(v)
      },
      score_cellchat              = score[tool == "cellchat"][1],
      score_mnn                   = score[tool == "mnn"][1],
      score_liana                 = score[tool == "liana"][1],
      score_norm_cellchat         = score_norm[tool == "cellchat"][1],
      score_norm_mnn              = score_norm[tool == "mnn"][1],
      score_norm_liana            = score_norm[tool == "liana"][1],
      padj_cellchat               = padj[tool == "cellchat"][1],
      padj_mnn                    = padj[tool == "mnn"][1],
      padj_liana                  = padj[tool == "liana"][1],
      sig_cellchat                = any(is_sig[tool == "cellchat"]),
      sig_mnn                     = any(is_sig[tool == "mnn"]),
      sig_liana                   = any(is_sig[tool == "liana"]),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      n_tools_sig_consistent_sign = pmax(n_pos_sig, n_neg_sig),
      sign_consensus = dplyr::case_when(
        n_tools_sig == 0L                 ~ NA_character_,
        n_pos_sig > 0L & n_neg_sig == 0L  ~ "+",
        n_neg_sig > 0L & n_pos_sig == 0L  ~ "-",
        TRUE                              ~ "mixed"
      ),
      mean_norm_score = dplyr::coalesce(mean_norm_score_sig,
                                        mean_norm_score_present)
    )

  wide |>
    dplyr::arrange(
      dplyr::desc(n_tools_sig_consistent_sign),
      dplyr::desc(n_tools_sig),
      dplyr::desc(abs(mean_norm_score))
    ) |>
    dplyr::mutate(composite_rank = dplyr::row_number()) |>
    dplyr::select(
      composite_rank,
      lr_interaction, sender, receiver, ligand, receptor,
      n_tools_present, n_tools_sig, n_tools_sig_consistent_sign,
      sign_consensus, mean_norm_score,
      score_cellchat, score_mnn, score_liana,
      score_norm_cellchat, score_norm_mnn, score_norm_liana,
      padj_cellchat, padj_mnn, padj_liana,
      sig_cellchat, sig_mnn, sig_liana
    )
}

# --------------------------------------------------------------------
# Public: cross-contrast LR-pair leader board.
# --------------------------------------------------------------------

# Aggregate the per-(LR-pair, contrast) cross-tool ranking long tibble
# into a per-LR-pair leader board across the five contrasts. Mirror of
# `build_tf_leader_board()`; leader rule defaults to
# `n_tools_sig_consistent_sign >= 2 OR n_tools_sig >= 3` per the G3
# locked decision.
#
# Arguments:
#   ranking_long           long tibble obtained by rbinding the output
#                          of `rank_lr_cross_tool()` across contrasts,
#                          with an added `contrast` column.
#                          Required columns: lr_interaction, sender,
#                          receiver, ligand, receptor, contrast,
#                          n_tools_sig, n_tools_sig_consistent_sign,
#                          sign_consensus, mean_norm_score.
#   leader_rule            optional vectorised closure
#                          `function(tbl) -> logical(nrow(tbl))`.
#                          Default identical to the TF leader rule.
#
# Returns one-row-per-leader tibble sorted by leader_score desc with
# the same tie-break chain as `build_tf_leader_board`.
build_lr_cross_tool_leaderboard <- function(ranking_long,
                                            leader_rule = NULL) {
  stopifnot(is.data.frame(ranking_long))
  required <- c("lr_interaction", "sender", "receiver", "ligand",
                "receptor", "contrast", "n_tools_sig",
                "n_tools_sig_consistent_sign", "sign_consensus",
                "mean_norm_score")
  missing  <- setdiff(required, names(ranking_long))
  if (length(missing) > 0L) {
    stop(sprintf(
      "build_lr_cross_tool_leaderboard: ranking_long missing columns: %s",
      paste(missing, collapse = ", ")), call. = FALSE)
  }
  if (is.null(leader_rule)) {
    leader_rule <- function(row) {
      row$n_tools_sig_consistent_sign >= 2L |
        row$n_tools_sig                >= 3L
    }
  }

  if (nrow(ranking_long) == 0L) {
    return(tibble::tibble(
      lr_interaction                  = character(0),
      sender                          = character(0),
      receiver                        = character(0),
      ligand                          = character(0),
      receptor                        = character(0),
      n_contrasts_leader              = integer(0),
      n_contrasts_consistent_sign_ge2 = integer(0),
      n_contrasts_sig_ge3             = integer(0),
      max_consistent_sign             = integer(0),
      max_n_tools_sig                 = integer(0),
      max_mean_norm_score             = numeric(0),
      dominant_sign                   = character(0),
      contrasts_summary               = character(0),
      leader_score                    = numeric(0)
    ))
  }

  is_leader <- as.logical(leader_rule(ranking_long))
  if (length(is_leader) != nrow(ranking_long)) {
    stop("build_lr_cross_tool_leaderboard: leader_rule must return a ",
         "logical vector of length nrow(ranking_long).", call. = FALSE)
  }
  is_leader[is.na(is_leader)] <- FALSE
  ranking_long$is_leader <- is_leader

  reduce_dominant_sign <- function(sign_consensus, is_leader) {
    sc <- sign_consensus[is_leader]
    sc <- sc[!is.na(sc)]
    if (length(sc) == 0L) return(NA_character_)
    if (all(sc == "+")) return("+")
    if (all(sc == "-")) return("-")
    "mixed"
  }
  reduce_contrasts_summary <- function(contrast, n_tools_sig,
                                       sign_consensus, is_leader) {
    ix <- which(is_leader)
    if (length(ix) == 0L) return(NA_character_)
    paste(sprintf("%s:%d/%s",
                  contrast[ix],
                  n_tools_sig[ix],
                  ifelse(is.na(sign_consensus[ix]), "NA",
                         sign_consensus[ix])),
          collapse = " | ")
  }

  # The per-pair aggregation key is the LR-pair tuple including the
  # cell-type sender / receiver, NOT the lr_interaction alone --
  # different (sender, receiver) cells of the same lr_interaction are
  # distinct biological signals.
  per_pair <- ranking_long |>
    dplyr::group_by(sender, receiver, ligand, receptor) |>
    dplyr::summarise(
      lr_interaction                  = dplyr::first(lr_interaction),
      n_contrasts_leader              = as.integer(sum(is_leader,
                                                       na.rm = TRUE)),
      n_contrasts_consistent_sign_ge2 = as.integer(sum(
        n_tools_sig_consistent_sign >= 2L, na.rm = TRUE)),
      n_contrasts_sig_ge3             = as.integer(sum(
        n_tools_sig >= 3L, na.rm = TRUE)),
      max_consistent_sign             = as.integer(suppressWarnings(
        max(n_tools_sig_consistent_sign, na.rm = TRUE))),
      max_n_tools_sig                 = as.integer(suppressWarnings(
        max(n_tools_sig, na.rm = TRUE))),
      max_mean_norm_score             = suppressWarnings(
        max(abs(mean_norm_score), na.rm = TRUE)),
      dominant_sign                   = reduce_dominant_sign(
        sign_consensus, is_leader),
      contrasts_summary               = reduce_contrasts_summary(
        contrast, n_tools_sig, sign_consensus, is_leader),
      .groups = "drop"
    ) |>
    dplyr::filter(n_contrasts_leader >= 1L)

  is_neg_inf <- function(x) is.numeric(x) & is.finite(x) == FALSE & x < 0
  per_pair$max_consistent_sign[is_neg_inf(per_pair$max_consistent_sign)] <- 0L
  per_pair$max_n_tools_sig[is_neg_inf(per_pair$max_n_tools_sig)]         <- 0L
  per_pair$max_mean_norm_score[is_neg_inf(per_pair$max_mean_norm_score)] <- NA_real_

  per_pair$leader_score <-
    5 * per_pair$n_contrasts_consistent_sign_ge2 +
    per_pair$n_contrasts_sig_ge3 +
    per_pair$max_consistent_sign / 5

  per_pair |>
    dplyr::arrange(dplyr::desc(leader_score),
                   dplyr::desc(n_contrasts_consistent_sign_ge2),
                   dplyr::desc(max_consistent_sign),
                   dplyr::desc(max_mean_norm_score),
                   lr_interaction, sender, receiver) |>
    dplyr::select(lr_interaction, sender, receiver, ligand, receptor,
                  n_contrasts_leader,
                  n_contrasts_consistent_sign_ge2,
                  n_contrasts_sig_ge3,
                  max_consistent_sign,
                  max_n_tools_sig,
                  max_mean_norm_score,
                  dominant_sign,
                  contrasts_summary,
                  leader_score)
}

# --------------------------------------------------------------------
# Public: per-contrast table formatter.
# --------------------------------------------------------------------

# Format the top-n rows of a per-contrast cross-tool LR ranking as a
# knitr::kable for inline display. Compact display columns drop the
# raw / normalised score and padj triples to keep the table narrow;
# pass `include_per_tool = TRUE` for the full per-tool detail.
#
# Arguments:
#   ranking_tbl       output of `rank_lr_cross_tool()`.
#   contrast          contrast string, used in the default caption.
#   n                 top-n rows. Default 10.
#   include_per_tool  include per-tool score / padj / sig columns.
#                     Default FALSE.
#   caption           optional caption override.
format_lr_ranking_table <- function(ranking_tbl, contrast, n = 10,
                                    include_per_tool = FALSE,
                                    caption = NULL) {
  top <- head(ranking_tbl, n)
  display_cols <- c("composite_rank", "lr_interaction",
                    "sender", "receiver",
                    "n_tools_sig", "n_tools_sig_consistent_sign",
                    "sign_consensus", "mean_norm_score")
  if (isTRUE(include_per_tool)) {
    display_cols <- c(display_cols,
                      "score_norm_cellchat",
                      "score_norm_mnn",
                      "score_norm_liana",
                      "sig_cellchat", "sig_mnn", "sig_liana")
  }
  if (is.null(caption)) {
    caption <- sprintf(
      paste0("Top %d cross-tool LR pairs for contrast '%s' (sorted by ",
             "n_tools_sig_consistent_sign desc, then n_tools_sig desc, ",
             "then |mean_norm_score| desc; per-tool sig at FDR<0.10)."),
      min(n, nrow(top)), contrast)
  }
  knitr::kable(top[, display_cols], digits = 3, caption = caption)
}

# --------------------------------------------------------------------
# Public: LR x tool consistency heatmap.
# --------------------------------------------------------------------

# Plot a per-(LR-pair, tool) heatmap of normalised scores for a given
# contrast across the leader LR-pairs. Mirror of E3's TF x modality
# heatmap. Rows alphabetical so no LR pair is visually privileged;
# cell text marks per-tool significance (`*` = FDR<0.05; `.` = FDR<0.10);
# NA cells (LR pair absent from a tool's per-contrast universe) are grey.
#
# Arguments:
#   ranking_tbl   output of `rank_lr_cross_tool()` for one contrast.
#   contrast      contrast name (for plot title).
#   n             max leaders to display. Default 25.
#   padj_strict   FDR threshold for the strict marker. Default 0.05.
#   padj_cut      FDR threshold for the loose marker. Default 0.10.
plot_lr_cross_tool_heatmap <- function(ranking_tbl, contrast,
                                       n = 25,
                                       padj_strict = 0.05,
                                       padj_cut    = 0.10) {
  tools_order <- c("cellchat", "mnn", "liana")
  top <- head(ranking_tbl, n)
  if (nrow(top) == 0L) {
    grid::grid.newpage()
    grid::grid.text("No LR pairs to display.")
    return(invisible(NULL))
  }
  # Alphabetical row order keyed by lr_interaction + sender + receiver.
  row_key <- sprintf("%s | %s->%s",
                     top$lr_interaction, top$sender, top$receiver)
  row_ord <- order(row_key)
  top    <- top[row_ord, ]
  row_key <- row_key[row_ord]

  score_mat <- as.matrix(top[, paste0("score_norm_", tools_order)])
  padj_mat  <- as.matrix(top[, paste0("padj_",       tools_order)])
  rownames(score_mat) <- row_key
  rownames(padj_mat)  <- row_key
  colnames(score_mat) <- tools_order
  colnames(padj_mat)  <- tools_order

  range_abs <- max(abs(score_mat[is.finite(score_mat)]), na.rm = TRUE)
  if (!is.finite(range_abs) || range_abs == 0) range_abs <- 1
  col_fun <- circlize::colorRamp2(c(-range_abs, 0, range_abs),
                                  c("#3a4cc0", "white", "#b40426"))

  cell_fn <- function(j, i, x, y, width, height, fill) {
    p <- padj_mat[i, j]
    if (is.na(p)) return()
    if (p < padj_strict) {
      grid::grid.text("*", x, y,
                      gp = grid::gpar(fontsize = 11, col = "black"))
    } else if (p < padj_cut) {
      grid::grid.text(".", x, y,
                      gp = grid::gpar(fontsize = 13, col = "black"))
    }
  }

  ComplexHeatmap::Heatmap(
    score_mat,
    name              = "norm score",
    col               = col_fun,
    na_col            = "grey90",
    cluster_rows      = FALSE,
    cluster_columns   = FALSE,
    row_names_side    = "left",
    row_names_gp      = grid::gpar(fontsize = 8),
    column_names_gp   = grid::gpar(fontsize = 10),
    column_names_rot  = 0,
    cell_fun          = cell_fn,
    column_title      = sprintf(
      "Cross-tool LR consistency at '%s' (alphabetical; * FDR<%.2f, . FDR<%.2f)",
      contrast, padj_strict, padj_cut),
    column_title_gp   = grid::gpar(fontsize = 10),
    heatmap_legend_param = list(
      title_gp = grid::gpar(fontsize = 9),
      labels_gp = grid::gpar(fontsize = 8)
    )
  )
}

# ====================================================================
# Phase G4: axis-restricted L-R analysis
# ====================================================================
#
# Mirror of `R/tf_inference.R::score_tf_per_axis()` and
# `R/kinase_inference.R::score_kinase_per_axis()` for the L-R layer.
# The unit of analysis is the per-(sender, receiver, ligand, receptor)
# CCC cell rather than a single regulator symbol; everything else
# parallels E4/F4 so the axis-restricted TSVs at sections 14.2, 15.2,
# and 16.2 are read side-by-side with identical column conventions.
#
# Input contract. The primary input is `ranking_long`, the rbind across
# the five contrasts of `rank_lr_cross_tool()` outputs (built in
# `rmd/15_ccc_mechanism.Rmd` by the G3 chunks as `ccc_ranking_long`).
# Each row is one (contrast, sender, receiver, ligand, receptor) cell
# with the cross-tool aggregates `mean_norm_score`, `n_tools_present`,
# `n_tools_sig`, `n_tools_sig_consistent_sign`, and `sign_consensus`
# already computed. Scoring at the axis layer is therefore a further
# aggregation across axis_contrasts of the already-cross-tool-aggregated
# `mean_norm_score`; the E4/F4 cache shape (modality -> contrast ->
# per-source score) is replaced by the G3 long tibble (contrast ->
# per-cell cross-tool mean_norm_score).
#
# Score source: `mean_norm_score` from G3. This is the mean of per-tool
# normalised scores across the SIGNIFICANT subset of tools at each
# (lr-cell, contrast), falling back to the mean over all PRESENT tools
# when no tool is significant (G3 coalesce). Signed direction (+
# upregulated in primary genotype, - downregulated). Magnitudes are
# normalised to [-1, +1] within each (tool, contrast) by dividing by
# `max(|score|)`, so cross-tool means are scale-comparable even though
# raw CellChat delta-prob (~10^-4), MultiNicheNet prioritization_score
# (~0-1), and LIANA+ prioritization_score (~[-2, +2]) live on disjoint
# scales.
#
# Cross-tool noise filter (locked default min_n_tools_present = 2L).
# The G3 ranking_long contains 232k rows across 47k distinct (s, r, l,
# rec) cells; 87% of rows are single-tool (MultiNicheNet-only, because
# MNN's group_prioritization_tbl contains 41,879 cells per contrast vs
# ~10k for LIANA+ and ~1k for CellChat). MNN's prioritization_score is
# always positive (it's a rank-percentile in [0, 1]); LIANA+'s
# prioritization_score and CellChat's delta_prob are signed. Without a
# baseline cross-tool support filter, the axis ranking is dominated by
# MNN-only cells whose `mean_norm_score` reduces to a single positive
# MNN value -- the top of every axis would be MNN-driven positive
# pseudo-leaders with no cross-tool corroboration. `min_n_tools_present
# = 2L` excludes those MNN-only cells; the result is the natural
# parallel to E4/F4's `min_targets = 5L` baseline data-density filter
# (and to G3's leader rule which itself requires cross-tool agreement
# at >=2 tools). With this filter, every axis still carries
# 3000-6500 cells at mode="either" -- well above the plan's >=5 sanity
# bound.
#
# Hybrid output (axis-mean + per-contrast columns; F4-locked shape).
# Each row carries:
#   * mean_activity_in_axis_contrasts: simple mean of mean_norm_score
#     across (sender, receiver, ligand, receptor, axis_contrast) cells
#     with finite values. For axis 3 (interaction_metabolic, single
#     axis_contrast), this collapses to the per-contrast value by
#     construction; for axes 1/2 (NLGF contrasts, two contrasts each),
#     this averages the two per-contrast cross-tool means.
#   * score_at_<contrast>: per-axis-contrast cross-tool mean (the
#     per-row mean_norm_score from G3). Unique axis_contrasts across
#     all axes are emitted as columns; rows whose axis does not include
#     the contrast carry NA in that column (mirror of the F4 hybrid
#     padding shape).
#
# Mode argument (`"either"` default, `"both"` documented alternative).
# `"either"` (the locked default) qualifies a cell if EITHER the
# ligand OR the receptor lies in the axis gene universe; mirror of
# E3's "TF significant in any modality" permissive logic. `"both"`
# qualifies only cells where the ligand AND receptor BOTH lie in the
# universe (stricter; intersects to a smaller set). Both modes
# populate; with `min_n_tools_present = 2L`, "either" gives 3-6k cells
# per axis and "both" gives 0.4-4k cells per axis (axis 3 is the
# narrowest because its universe excludes most NLGF-axis CCC pairs).
# The plan locks "either" as the primary default; "both" is a
# documented robustness check, not a separate analysis.
#
# Per-pair universe-membership column. `n_lr_ends_in_axis_universe`
# is 1 or 2 depending on the mode filter outcome (= number of {ligand,
# receptor} ends in the axis universe). This is the L-R analogue of
# `n_targets_in_axis_universe` in the TF / kinase axis-restricted
# TSVs (E4 / F4); same role (informative diagnostic on axis-relevance
# strength of each surviving cell), different range (the L-R case
# tops out at 2, whereas TFs / kinases can have dozens to hundreds of
# in-universe targets). The column name is deliberately distinct
# (`_ends` vs `_targets`) so the downstream TSV reader sees the
# semantic shift and does not mis-interpret high-target TFs as
# comparable to two-end LR pairs.
#
# Sorting. Per-axis rows are sorted by `-abs(mean_activity_in_axis_contrasts)`,
# then by `lr_interaction` asc, then `sender` asc, then `receiver` asc.
# Mirror of the TF / kinase sort with the per-cell key disambiguated
# across senders / receivers.

# Filter an L-R tibble to pairs whose ligand and/or receptor lies in
# the supplied axis gene universe.
#
# Arguments:
#   lr_tbl          tibble with character columns `ligand` and
#                   `receptor`. Other columns are preserved. Use either
#                   `ccc_lr_long` (per-(tool, contrast, lr-pair)) or
#                   `ccc_ranking_long` (cross-tool aggregated per-
#                   (contrast, lr-pair)); the filter is column-agnostic.
#   axis_universe   character vector of mouse gene symbols defining the
#                   axis. Same shape as the `universe` slot of an entry
#                   in `build_axis_gene_universe()`'s return.
#   mode            "either" (default) or "both"; see module header.
#
# Returns the filtered tibble with three extra columns:
#   ligand_in_axis_universe       logical.
#   receptor_in_axis_universe     logical.
#   n_lr_ends_in_axis_universe    integer (1 or 2 for kept rows).
restrict_lr_to_universe <- function(lr_tbl, axis_universe,
                                    mode = c("either", "both")) {
  mode <- match.arg(mode)
  stopifnot(is.data.frame(lr_tbl) || tibble::is_tibble(lr_tbl),
            all(c("ligand", "receptor") %in% names(lr_tbl)),
            is.character(axis_universe) || is.null(axis_universe))
  if (length(axis_universe) == 0L) {
    keep <- lr_tbl[FALSE, , drop = FALSE]
    keep$ligand_in_axis_universe    <- logical(0)
    keep$receptor_in_axis_universe  <- logical(0)
    keep$n_lr_ends_in_axis_universe <- integer(0)
    return(tibble::as_tibble(keep))
  }
  uni <- unique(as.character(axis_universe))
  lig_in <- as.character(lr_tbl$ligand)   %in% uni
  rec_in <- as.character(lr_tbl$receptor) %in% uni
  keep_mask <- if (mode == "either") lig_in | rec_in else lig_in & rec_in
  out <- lr_tbl[keep_mask, , drop = FALSE]
  out$ligand_in_axis_universe    <- lig_in[keep_mask]
  out$receptor_in_axis_universe  <- rec_in[keep_mask]
  out$n_lr_ends_in_axis_universe <- as.integer(
    out$ligand_in_axis_universe + out$receptor_in_axis_universe)
  tibble::as_tibble(out)
}

# Score each L-R pair per D2 axis using the G3 cross-tool aggregated
# ranking_long. Hybrid output: axis-mean across axis_contrast cells AND
# per-contrast cross-tool means.
#
# Arguments:
#   ranking_long          rbind-of-per-contrast `rank_lr_cross_tool()`
#                         outputs with a `contrast` column. Required
#                         columns: contrast, sender, receiver, ligand,
#                         receptor, lr_interaction, mean_norm_score,
#                         n_tools_present, sign_consensus.
#   axis_universes        output of `build_axis_gene_universe()` from
#                         `R/tf_inference.R` (gene-level universes per
#                         axis; same input as E4 / F4).
#   mode                  "either" (default) or "both"; see module
#                         header. Forwarded to `restrict_lr_to_universe()`.
#   min_n_tools_present   integer >= 1. Default 2L. Locked default
#                         excludes MNN-only single-tool cells from the
#                         cross-tool ranking (see module header for the
#                         rationale).
#   score_col             character; column of ranking_long carrying the
#                         per-(contrast, cell) signed score to aggregate.
#                         Default "mean_norm_score" (the G3 cross-tool
#                         coalesce of sig and present means).
#
# Returns a tibble with columns:
#   axis                              character; axis name.
#   sender, receiver, ligand, receptor character; per-pair key.
#   lr_interaction                    character; convenience field.
#   mean_activity_in_axis_contrasts   numeric; mean of score_col across
#                                     axis_contrast cells with finite
#                                     score.
#   sd_activity_in_axis_contrasts     numeric; NA when n_cells_used==1.
#   n_cells_used                      integer; count of finite axis-
#                                     contrast cells averaged.
#   n_lr_ends_in_axis_universe        integer; 1 or 2 per the mode filter.
#   leader_rank                       integer; per-axis rank by |mean|
#                                     desc, then lr_interaction,
#                                     sender, receiver asc.
#   score_at_<contrast>               numeric; per-axis-contrast cross-
#                                     tool mean. One column per unique
#                                     axis_contrast across all axes;
#                                     NA where the row's axis does not
#                                     include the contrast (mirror of
#                                     F4 hybrid padding shape).
score_lr_per_axis <- function(ranking_long,
                              axis_universes,
                              mode                = c("either", "both"),
                              min_n_tools_present = 2L,
                              score_col           = "mean_norm_score") {
  mode <- match.arg(mode)
  stopifnot(is.data.frame(ranking_long) || tibble::is_tibble(ranking_long),
            is.list(axis_universes),
            !is.null(names(axis_universes)),
            length(min_n_tools_present) == 1L,
            min_n_tools_present >= 1L)
  required <- c("contrast", "sender", "receiver", "ligand", "receptor",
                "lr_interaction", score_col, "n_tools_present")
  missing  <- setdiff(required, names(ranking_long))
  if (length(missing) > 0L) {
    stop(sprintf(
      "score_lr_per_axis: ranking_long missing columns: %s",
      paste(missing, collapse = ", ")), call. = FALSE)
  }
  rl <- ranking_long[!is.na(ranking_long$n_tools_present) &
                       ranking_long$n_tools_present >= min_n_tools_present, ,
                     drop = FALSE]

  all_contrasts <- unique(unlist(
    lapply(axis_universes, function(a) a$axis_contrasts),
    use.names = FALSE
  ))
  per_contrast_cols <- paste0("score_at_", all_contrasts)

  rows <- list()
  for (ax in names(axis_universes)) {
    aux  <- axis_universes[[ax]]
    uni  <- aux$universe
    ctrs <- aux$axis_contrasts
    if (length(uni) == 0L || length(ctrs) == 0L) next

    # Filter ranking_long to the axis_contrasts; the universe filter is
    # applied per-cell via restrict_lr_to_universe().
    sub <- rl[rl$contrast %in% ctrs, , drop = FALSE]
    if (nrow(sub) == 0L) next
    sub <- restrict_lr_to_universe(sub, uni, mode = mode)
    if (nrow(sub) == 0L) next

    sub$.score <- as.numeric(sub[[score_col]])
    long_df <- sub[is.finite(sub$.score), , drop = FALSE]
    if (nrow(long_df) == 0L) next

    # Per-(sender, receiver, ligand, receptor) aggregate across axis_contrasts.
    per_cell <- long_df |>
      dplyr::group_by(sender, receiver, ligand, receptor) |>
      dplyr::summarise(
        lr_interaction                  = dplyr::first(lr_interaction),
        mean_activity_in_axis_contrasts = mean(.score, na.rm = TRUE),
        sd_activity_in_axis_contrasts   = if (dplyr::n() > 1L) {
          stats::sd(.score, na.rm = TRUE)
        } else {
          NA_real_
        },
        n_cells_used                    = as.integer(sum(is.finite(.score))),
        n_lr_ends_in_axis_universe      = as.integer(max(
          n_lr_ends_in_axis_universe, na.rm = TRUE)),
        .groups = "drop"
      )

    # Hybrid per-contrast columns. Each cell has at most one row per
    # axis_contrast in long_df (ranking_long is per-(contrast, cell)),
    # so the per-contrast value is simply the score at that contrast.
    per_cc <- long_df |>
      dplyr::select(sender, receiver, ligand, receptor, contrast,
                    .score) |>
      tidyr::pivot_wider(names_from  = contrast,
                         values_from = .score,
                         values_fn   = mean)
    # Rename the contrast columns to score_at_<contrast> and pad any
    # missing axis contrasts (defensive; usually pivot_wider creates
    # all the contrasts present in long_df).
    ctr_cols <- setdiff(names(per_cc), c("sender", "receiver",
                                         "ligand", "receptor"))
    names(per_cc)[match(ctr_cols, names(per_cc))] <- paste0("score_at_",
                                                            ctr_cols)
    per_cell <- dplyr::left_join(per_cell, per_cc,
                                 by = c("sender", "receiver",
                                        "ligand", "receptor"))

    # Pad NA columns for any contrast not present in this axis subset
    # so the schema is uniform across axes.
    missing_cols <- setdiff(per_contrast_cols, names(per_cell))
    for (mc in missing_cols) per_cell[[mc]] <- NA_real_

    per_cell <- per_cell[order(-abs(per_cell$mean_activity_in_axis_contrasts),
                               per_cell$lr_interaction,
                               per_cell$sender,
                               per_cell$receiver), , drop = FALSE]
    per_cell$leader_rank <- seq_len(nrow(per_cell))
    per_cell$axis <- ax

    rows[[ax]] <- per_cell[, c("axis", "sender", "receiver",
                               "ligand", "receptor", "lr_interaction",
                               "mean_activity_in_axis_contrasts",
                               "sd_activity_in_axis_contrasts",
                               "n_cells_used",
                               "n_lr_ends_in_axis_universe",
                               "leader_rank",
                               per_contrast_cols), drop = FALSE]
  }
  if (length(rows) == 0L) {
    empty <- tibble::tibble(
      axis                              = character(0),
      sender                            = character(0),
      receiver                          = character(0),
      ligand                            = character(0),
      receptor                          = character(0),
      lr_interaction                    = character(0),
      mean_activity_in_axis_contrasts   = numeric(0),
      sd_activity_in_axis_contrasts     = numeric(0),
      n_cells_used                      = integer(0),
      n_lr_ends_in_axis_universe        = integer(0),
      leader_rank                       = integer(0)
    )
    for (col in per_contrast_cols) empty[[col]] <- numeric(0)
    return(empty)
  }
  tibble::as_tibble(do.call(rbind, rows))
}

# Format the top-n axis-restricted L-R cells for one axis as a
# knitr::kable. Mirror of `format_axis_restricted_table()` in
# `R/tf_inference.R` and `format_axis_restricted_kinase_table()` in
# `R/kinase_inference.R`. Auto-includes per-contrast columns whose
# values are not all NA in the axis subset.
#
# Arguments:
#   axis_tbl     output of `score_lr_per_axis()`.
#   axis_name    one of `unique(axis_tbl$axis)`.
#   n            top-n rows (default 15).
#   caption      optional caption override.
#
# Returns a knitr::kable suitable for `print()` inside an
# `results = 'asis'` chunk.
format_axis_restricted_lr_table <- function(axis_tbl, axis_name,
                                            n = 15, caption = NULL) {
  sub <- axis_tbl[axis_tbl$axis == axis_name, , drop = FALSE]
  if (nrow(sub) == 0L) {
    return(knitr::kable(
      data.frame(message = sprintf(
        "No L-R pairs survive the axis '%s' universe-overlap filter.",
        axis_name)),
      caption = caption %||%
        sprintf("Axis-restricted L-R pairs: %s", axis_name),
      row.names = FALSE))
  }
  sub <- sub[order(sub$leader_rank), , drop = FALSE]
  top <- head(sub, n)

  per_contrast_cols <- grep("^score_at_", names(top), value = TRUE)
  if (length(per_contrast_cols) > 0L) {
    axis_subset <- axis_tbl[axis_tbl$axis == axis_name, , drop = FALSE]
    keep_pc <- vapply(per_contrast_cols, function(col) {
      !all(is.na(axis_subset[[col]]))
    }, logical(1))
    per_contrast_cols <- per_contrast_cols[keep_pc]
  }

  display <- data.frame(
    rank       = top$leader_rank,
    lr_pair    = top$lr_interaction,
    sender     = top$sender,
    receiver   = top$receiver,
    mean_score = round(top$mean_activity_in_axis_contrasts, 3),
    sd_score   = round(top$sd_activity_in_axis_contrasts,   3),
    n_cells    = top$n_cells_used,
    n_ends     = top$n_lr_ends_in_axis_universe,
    stringsAsFactors = FALSE
  )
  for (col in per_contrast_cols) {
    short <- sub("^score_at_", "", col)
    display[[short]] <- round(top[[col]], 3)
  }
  if (is.null(caption)) {
    pc_part <- if (length(per_contrast_cols) > 0L) {
      sprintf(
        paste0(" Hybrid view also shows the per-contrast cross-tool ",
               "mean for each axis-relevant contrast (%s)."),
        paste(sub("^score_at_", "", per_contrast_cols), collapse = ", "))
    } else {
      ""
    }
    caption <- sprintf(
      paste0("Top %d axis-restricted L-R pairs for axis '%s' (sorted by ",
             "|mean cross-tool normalised score| across axis-relevant ",
             "contrasts; pairs filtered to those with at least one ",
             "ligand or receptor in the axis gene universe and at least ",
             "two CCC tools observing the pair at each contrast). ",
             "n_cells is the count of finite axis-contrast cells ",
             "averaged; n_ends is 1 or 2 depending on whether the ",
             "ligand alone, the receptor alone, or both lie in the ",
             "axis universe.%s"),
      nrow(top), axis_name, pc_part)
  }
  knitr::kable(display, caption = caption, row.names = FALSE)
}

# Lollipop chart of the top-n axis-restricted L-R cells for one axis.
# Mirror of `plot_axis_lollipop()` in `R/tf_inference.R` and
# `plot_axis_lollipop_kinase()` in `R/kinase_inference.R`. y-axis
# labels concatenate `lr_interaction | sender -> receiver` so distinct
# cell-type cells of the same L-R biochemistry are distinguishable.
#
# Arguments:
#   axis_tbl     output of `score_lr_per_axis()`.
#   axis_name    one of `unique(axis_tbl$axis)`.
#   n            top-n rows (default 12).
#   title        optional plot title; default derives from axis_name.
#
# Returns a ggplot object.
plot_axis_lollipop_lr <- function(axis_tbl, axis_name, n = 12,
                                  title = NULL) {
  sub <- axis_tbl[axis_tbl$axis == axis_name, , drop = FALSE]
  if (nrow(sub) == 0L) {
    return(ggplot2::ggplot() +
             ggplot2::annotate("text", x = 0, y = 0,
                               label = sprintf(
                                 "No L-R pairs for axis '%s'.", axis_name)) +
             ggplot2::theme_void())
  }
  sub <- sub[order(sub$leader_rank), , drop = FALSE]
  top <- head(sub, n)
  # Re-sort by signed mean so positive- and negative-driver L-R pairs
  # are visually separated. Mirror of plot_axis_lollipop() (TF).
  top <- top[order(top$mean_activity_in_axis_contrasts), , drop = FALSE]
  top$label <- sprintf("%s | %s -> %s",
                       top$lr_interaction, top$sender, top$receiver)
  top$label <- factor(top$label, levels = top$label)
  top$sign  <- ifelse(top$mean_activity_in_axis_contrasts >= 0,
                      "positive", "negative")
  if (is.null(title)) {
    title <- sprintf("Axis-restricted L-R pairs: %s (top %d by |mean score|)",
                     axis_name, nrow(top))
  }
  ggplot2::ggplot(top,
                  ggplot2::aes(x = mean_activity_in_axis_contrasts,
                               y = label,
                               colour = sign)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0,
                                       xend = mean_activity_in_axis_contrasts,
                                       y    = label,
                                       yend = label),
                          linewidth = 0.6) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                        colour = "grey60") +
    ggplot2::scale_colour_manual(values = c(positive = "#b40426",
                                            negative = "#3a4cc0"),
                                 name   = "mean score sign") +
    ggplot2::labs(title = title,
                  x     = "mean cross-tool normalised score (axis contrasts)",
                  y     = NULL) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
                   axis.text.y        = ggplot2::element_text(size = 8))
}

# ====================================================================
# Phase G5: cross-tool + cross-contrast L-R verdict aggregator
# ====================================================================
#
# Mirror of `R/tf_inference.R::build_tf_verdict_table()` and
# `R/kinase_inference.R::build_kinase_verdict_table()`. Thin aggregator
# producing one row per axis from the G4 axis-restricted long tibble
# (output of `score_lr_per_axis()`), with the F4/E5/F5-locked hybrid
# columns (axis-mean PLUS per-contrast min/max) so the L-R verdict
# reads side-by-side with the TF verdict (14.3) and the kinase verdict
# (15.3). The cross-tool consistency column counts top-N cells that
# also reach the G3 cross-tool leader rule at any axis-relevant
# contrast, mirroring F5's "sig in any axis-relevant contrast" pattern
# but substituting the G3 cross-tool leader rule
# (n_tools_sig_consistent_sign >= 2L OR n_tools_sig >= 3L) for F3's
# per-contrast FDR<0.10 ulm gate.
#
# Schema rationale.
#   * `axis` preserves the natural amyloid_activation /
#     synaptic_suppression / interaction_metabolic encounter order
#     rather than alphabetising; matches the 1/2/3 ordering used by
#     section-13 D2, the section-14 TF verdict, and the section-15
#     kinase verdict so all three mechanism-layer verdict TSVs read
#     in parallel.
#   * `top_lr_pairs` is `<lr_interaction> (<sender>-><receiver>)`
#     semicolon-separated. The cell-type pair is essential at the L-R
#     layer: the same biochemistry (e.g. `Apoe_Trem2`) at
#     Microglia_DAM -> Microglia_DAM and Astrocyte -> Microglia_DAM
#     are distinct biological signals, and collapsing to lr_interaction
#     alone would lose this layer of the signal. Semicolon outer
#     separator avoids confusion with commas inside the parenthesised
#     cell-type tag.
#   * `top_lr_signs` is comma-separated +/- per top-N cell, in the
#     same order as `top_lr_pairs`.
#   * `mean_score_range` is the [min, max] of `mean_activity_in_axis_contrasts`
#     across the top-N rows (axis-mean view).
#   * `per_contrast_score_range` is the [min, max] of per-contrast
#     cross-tool means across the top-N x axis_contrasts cells
#     (F4 hybrid view). At axis 3 with one axis_contrast this
#     collapses to the axis-mean range by construction; preserved for
#     schema consistency with axes 1 / 2.
#   * `per_contrast_summary` is "contrast:[min, max]" per axis_contrast,
#     semicolon-separated; same convention as F5 / E5.
#   * `n_top_sig_in_axis_contrasts` is the integer count of top-N LR
#     cells that appear in the G3 cross-tool leader board AND are
#     flagged as a leader at >=1 of the axis-relevant contrasts in
#     the leader board's `contrasts_summary` column. The G3 leader
#     board lists every (s, r, l, rec) cell that reaches the cross-
#     tool leader rule at >=1 of the 5 canonical contrasts; the
#     `contrasts_summary` packs the per-leader-contrast detail as
#     "<contrast>:<n_tools_sig>/<sign_consensus>" entries
#     pipe-separated.
#   * `evidence_summary` is a free-text column filled by the caller
#     via `evidence_summaries`, so the editorial prose lives in the
#     Rmd chunk rather than in code (parallel to F5 / E5).
#
# Expected count behaviour and complementary-views caveat.
# `n_top_sig_in_axis_contrasts` is expected to be small (frequently
# 0/0/0 at top-N = 5) because the G3 cross-contrast leader rule and
# the G4 per-axis-contrast magnitude rule weight different biology:
# G3 rewards pairs with modest scores observed consistently across
# multiple contrasts (e.g. `Sema6a_Plxna4` SOLE 5/5-contrast leader),
# while the G4 axis layer rewards pairs with top-of-normalised-
# distribution scores at the axis's 1-2 contrasts (often
# n_cells_used = 1 single-contrast cells, e.g. `Adam11_Itga4` axis-1
# rank 1 mean 1.000). The complementarity was documented at length in
# the G4 deviation register entry #5; the verdict prose should treat
# `n_top_sig_in_axis_contrasts > 0` as additional cross-tool
# corroboration, not as the dominant evidence signal. A 0/0/0 column
# is the honest reading of "the top-of-axis cells are not in the
# cross-contrast leader board", not a bug.
#
# Arguments:
#   axis_tbl              long tibble from `score_lr_per_axis()`. Must
#                         carry `axis`, `sender`, `receiver`,
#                         `ligand`, `receptor`, `lr_interaction`,
#                         `mean_activity_in_axis_contrasts`,
#                         `leader_rank`, plus the `score_at_<contrast>`
#                         columns for hybrid use.
#   lr_leaderboard        tibble from `build_lr_cross_tool_leaderboard()`.
#                         Must carry `sender`, `receiver`, `ligand`,
#                         `receptor`, `contrasts_summary`.
#   n_top_per_axis        top-N LR cells to summarise per axis.
#                         Default 5 (matches E5 / F5).
#   evidence_summaries    optional named list `axis_name -> string`
#                         containing the prose evidence-summary string
#                         to place in the `evidence_summary` column.
#                         Missing axes get NA_character_.
#
# Returns a data.frame with one row per axis and columns:
#   axis                          character;
#   top_lr_pairs                  semicolon-separated
#                                 "<lr_interaction> (<sender>-><receiver>)";
#   top_lr_signs                  comma-separated +/-;
#   mean_score_range              e.g. "[+0.98, +1.00]" (axis-mean);
#   per_contrast_score_range      e.g. "[+0.97, +1.00]" (hybrid);
#   per_contrast_summary          string like
#                                 "nlgf_in_maptki:[+0.98, +1.00]; nlgf_in_p301s:[+0.97, +1.00]";
#   n_top_sig_in_axis_contrasts   integer (0..n_top_per_axis); count
#                                 of top-N cells that are G3 cross-
#                                 tool leaders at any axis-relevant
#                                 contrast;
#   evidence_summary              character (possibly NA).
build_lr_verdict_table <- function(axis_tbl, lr_leaderboard,
                                   n_top_per_axis = 5L,
                                   evidence_summaries = NULL) {
  stopifnot(is.data.frame(axis_tbl) || tibble::is_tibble(axis_tbl),
            all(c("axis", "sender", "receiver", "ligand", "receptor",
                  "lr_interaction", "mean_activity_in_axis_contrasts",
                  "leader_rank") %in% names(axis_tbl)),
            is.data.frame(lr_leaderboard) ||
              tibble::is_tibble(lr_leaderboard),
            all(c("sender", "receiver", "ligand", "receptor",
                  "contrasts_summary") %in% names(lr_leaderboard)),
            is.numeric(n_top_per_axis), length(n_top_per_axis) == 1L,
            n_top_per_axis >= 1L)

  axes <- unique(as.character(axis_tbl$axis))
  if (is.null(evidence_summaries)) {
    evidence_summaries <- setNames(rep(NA_character_, length(axes)), axes)
  }

  per_contrast_cols <- grep("^score_at_", names(axis_tbl), value = TRUE)

  # Precompute the leader-board 4-key for matching against top-N rows.
  lb_key <- paste(lr_leaderboard$sender, lr_leaderboard$receiver,
                  lr_leaderboard$ligand,  lr_leaderboard$receptor,
                  sep = "|")

  rows <- lapply(axes, function(ax) {
    sub <- axis_tbl[axis_tbl$axis == ax, , drop = FALSE]
    sub <- sub[order(sub$leader_rank), , drop = FALSE]
    top <- head(sub, n_top_per_axis)
    if (nrow(top) == 0L) {
      return(data.frame(
        axis                              = ax,
        top_lr_pairs                      = NA_character_,
        top_lr_signs                      = NA_character_,
        mean_score_range                  = NA_character_,
        per_contrast_score_range          = NA_character_,
        per_contrast_summary              = NA_character_,
        n_top_sig_in_axis_contrasts       = 0L,
        evidence_summary                  = evidence_summaries[[ax]] %||%
                                              NA_character_,
        stringsAsFactors                  = FALSE
      ))
    }

    signs <- ifelse(top$mean_activity_in_axis_contrasts >= 0, "+", "-")
    range_str <- sprintf("[%+0.2f, %+0.2f]",
                         min(top$mean_activity_in_axis_contrasts),
                         max(top$mean_activity_in_axis_contrasts))

    # Hybrid per-contrast view. Identify axis-relevant per-contrast
    # columns by dropping those with all-NA values in the axis subset
    # (mirrors the F5 / E5 derivation). The surviving column names
    # also identify the axis_contrasts for the leader-board cross-
    # reference below.
    pc_keep <- per_contrast_cols[vapply(per_contrast_cols, function(col) {
      !all(is.na(sub[[col]]))
    }, logical(1))]
    axis_contrasts <- sub("^score_at_", "", pc_keep)

    if (length(pc_keep) == 0L) {
      pc_range_str <- NA_character_
      pc_summary   <- NA_character_
    } else {
      pc_values <- unlist(lapply(pc_keep, function(col) top[[col]]),
                          use.names = FALSE)
      pc_values <- pc_values[is.finite(pc_values)]
      if (length(pc_values) == 0L) {
        pc_range_str <- NA_character_
      } else {
        pc_range_str <- sprintf("[%+0.2f, %+0.2f]",
                                min(pc_values), max(pc_values))
      }
      pc_summary <- paste(vapply(pc_keep, function(col) {
        v <- top[[col]]
        v <- v[is.finite(v)]
        if (length(v) == 0L) {
          sprintf("%s:NA", sub("^score_at_", "", col))
        } else {
          sprintf("%s:[%+0.2f, %+0.2f]",
                  sub("^score_at_", "", col),
                  min(v), max(v))
        }
      }, character(1)), collapse = "; ")
    }

    # G5-specific cross-reference column: count of top-N LR cells that
    # appear in the G3 cross-tool leader board AND carry leadership at
    # >=1 axis-relevant contrast in the leader board's contrasts_summary.
    if (length(axis_contrasts) == 0L) {
      n_sig <- 0L
    } else {
      top_key <- paste(top$sender, top$receiver,
                       top$ligand, top$receptor, sep = "|")
      n_sig <- 0L
      for (i in seq_along(top_key)) {
        ix <- match(top_key[i], lb_key)
        if (is.na(ix)) next
        cs <- lr_leaderboard$contrasts_summary[ix]
        if (is.na(cs)) next
        # Each leader contrast in contrasts_summary is encoded as
        # "<contrast>:<n_tools_sig>/<sign_consensus>" entries separated
        # by " | ". A match means the cell is a G3 leader at that
        # contrast under the locked leader rule
        # (n_tools_sig_consistent_sign >= 2L OR n_tools_sig >= 3L).
        is_axis_leader <- any(vapply(axis_contrasts, function(c) {
          grepl(paste0("(^| \\| )", c, ":"), cs, fixed = FALSE)
        }, logical(1)))
        if (is_axis_leader) n_sig <- n_sig + 1L
      }
    }

    top_labels <- sprintf("%s (%s->%s)",
                          top$lr_interaction, top$sender, top$receiver)

    data.frame(
      axis                              = ax,
      top_lr_pairs                      = paste(top_labels, collapse = "; "),
      top_lr_signs                      = paste(signs, collapse = ","),
      mean_score_range                  = range_str,
      per_contrast_score_range          = pc_range_str,
      per_contrast_summary              = pc_summary,
      n_top_sig_in_axis_contrasts       = as.integer(n_sig),
      evidence_summary                  = evidence_summaries[[ax]] %||%
                                            NA_character_,
      stringsAsFactors                  = FALSE
    )
  })
  do.call(rbind, rows)
}

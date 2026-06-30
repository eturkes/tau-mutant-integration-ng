# Activation-trajectory pseudotime (P2-S1). Order microglia along the homeostatic->DAM
# axis with slingshot on the BATCH-CORRECTED harmony embedding, seeded by the substate
# biology. PRIMARY lineage = a forced single clean Homeostatic->DAM curve (clusterLabels =
# the 2 substate super-clusters -> MST on 2 nodes = one edge = one lineage; principal curve
# still fits the full high-dim cell cloud). IFN/Proliferative cells are OMITTED from the
# lineage (off-lineage flag + NA pt), not deleted -> per-unit omitted fraction is reported,
# never hidden. Pseudotime is an ACTIVATION ORDERING (position/extent of advance), NOT
# developmental time or potency -- direction rests on slingshot rooting + DAM-marker
# monotonicity, validated post-hoc, never a potency claim. All non-base calls namespace-
# qualified (targets attaches only `quarto`). Pure helpers unit-tested; the heavy
# orchestrator (build_activation_trajectory) is smoke-tested live on microglia_annotated.

# Marker SCORE-AXIS pseudotime proxy = raw UCell DAM minus Homeostatic. Assumption-light
# (no population z-centring): both are UCell scores on the same rank-based [0,1] scale, so
# their difference is a direct per-cell activation contrast. Used ONLY as a CONCORDANCE
# anchor for the slingshot ordering (it shares the marker SYSTEM that defines the substate
# labels -> catches gross trajectory failure, NOT a statistically independent robustness
# check). Pure.
score_axis_pseudotime <- function(dam_score, homeo_score) {
  stopifnot(is.numeric(dam_score), is.numeric(homeo_score),
            length(dam_score) == length(homeo_score),
            all(is.finite(dam_score)), all(is.finite(homeo_score)))
  dam_score - homeo_score
}

# Smithson-Verkuilen squeeze of a pseudotime vector into the OPEN interval (0,1) for a beta
# GLMM (P2-S2). Two steps: (1) min-max scale finite values to [0,1]; (2) y' = (y*(n-1)+0.5)/n
# with n = number of finite obs -> nudges the 0/1 endpoints inward (~0.5/n) so beta's
# unbounded logit stays finite. NA preserved as NA (off-lineage cells). Degenerate constant
# input -> all 0.5. Pure.
squeeze_unit_interval <- function(x) {
  stopifnot(is.numeric(x))
  ok <- is.finite(x)
  n <- sum(ok)
  stopifnot(n >= 2L)
  rng <- range(x[ok])
  scaled <- if (diff(rng) > 0) (x - rng[1L]) / diff(rng) else rep(0.5, length(x))
  out <- (scaled * (n - 1L) + 0.5) / n
  out[!ok] <- NA_real_
  out
}

# Spearman concordance between the slingshot pseudotime and the marker score-axis, over the
# cells where BOTH are finite (off-lineage pt is NA). Returns list(rho, n) -- rho recorded in
# provenance (concordance is reported, not gated; a low rho is an honest flag, not a build
# failure). Pure.
trajectory_concordance <- function(pt, score_axis) {
  stopifnot(is.numeric(pt), is.numeric(score_axis), length(pt) == length(score_axis))
  ok <- is.finite(pt) & is.finite(score_axis)
  stopifnot(sum(ok) >= 2L)
  list(rho = stats::cor(pt[ok], score_axis[ok], method = "spearman"), n = sum(ok))
}

# Fit slingshot on a reduced-dim embedding + cluster labels, return the pseudotime of the
# lineage TERMINATING in `terminal_clus` (Homeostatic->DAM). matrix in, list out -> pure +
# unit-testable on a synthetic embedding. 2-cluster labels guarantee a single lineage; >2
# (all-retained sensitivity, IFN present) may branch -> the DAM-terminal lineage is selected
# (longest if tied) and off-lineage cells keep slingshot's NA. Rooted at `start_clus` ->
# pseudotime increases away from the homeostatic root. seed set for belt-and-braces
# determinism (the principal curve is RNG-free given embedding + labels).
run_slingshot_lineage <- function(embedding, cluster_labels, start_clus, terminal_clus,
                                   seed = 42L) {
  stopifnot(is.matrix(embedding), !is.null(rownames(embedding)), ncol(embedding) >= 2L,
            all(is.finite(embedding)),
            length(cluster_labels) == nrow(embedding),
            start_clus %in% cluster_labels, terminal_clus %in% cluster_labels)
  old_kind <- RNGkind("Mersenne-Twister")
  on.exit(RNGkind(old_kind[1], old_kind[2], old_kind[3]), add = TRUE)
  set.seed(seed)
  res <- slingshot::slingshot(embedding, clusterLabels = as.character(cluster_labels),
                              start.clus = start_clus)
  lins <- slingshot::slingLineages(res)
  term_ok <- vapply(lins, function(l) l[length(l)] == terminal_clus, logical(1))
  if (!any(term_ok)) {
    stop("run_slingshot_lineage: no lineage terminates in '", terminal_clus,
         "' (lineages: ", paste(vapply(lins, paste, character(1), collapse = "->"),
                                collapse = " | "), ")")
  }
  idx <- which(term_ok)
  if (length(idx) > 1L) idx <- idx[which.max(vapply(lins[idx], length, integer(1)))]
  lin_name <- names(lins)[idx]
  pt_mat <- slingshot::slingPseudotime(res)
  pt <- pt_mat[, lin_name]
  names(pt) <- rownames(pt_mat)
  list(pt = pt, lineage = unname(lins[[idx]]), lineage_name = lin_name,
       n_lineages = length(lins), all_lineages = lapply(lins, unname))
}

# Pinned-stack + run-environment provenance for the trajectory build (mirror of
# reprocess_provenance): pkg versions, seed/RNG, thread snapshot, R version. Pure.
trajectory_provenance <- function(seed) {
  pkgs <- c("slingshot", "princurve", "TrajectoryUtils", "SingleCellExperiment", "UCell")
  list(
    versions = stats::setNames(
      vapply(pkgs, function(p) as.character(utils::packageVersion(p)), character(1)), pkgs),
    seed = seed,
    rng_kind = RNGkind(),
    threads = reprocess_thread_env(),
    r_version = as.character(getRversion())
  )
}

# Orchestrate the activation-trajectory target from microglia_annotated. Reads the cached
# harmony embedding + per-cell UCell scores + substate labels (NO recompute), fits the
# PRIMARY Homeostatic->DAM slingshot pseudotime (dims 1:primary_dims), a FIXED sensitivity
# table (alt dims + an all-cells-retained re-fit), and the marker score-axis concordance.
# Returns a COMPACT target (per-cell frame + per-unit table + lineage/sensitivity/provenance
# lists) -- NEVER the 612MB Seurat (cheap-render invariant; only the report-side qmd must stay
# light, but the target itself is small too). FORK CHOICES (dims, lineage substates) are
# PRE-DECLARED here, not retuned after inspecting contrasts.
build_activation_trajectory <- function(seurat_obj,
                                         primary_dims = 15L,
                                         sensitivity_dims = c(10L, 20L),
                                         lineage_substates = c("Homeostatic", "DAM"),
                                         root_substate = "Homeostatic",
                                         terminal_substate = "DAM",
                                         reduction = "harmony",
                                         substate_col = "microglia_substate",
                                         dam_col = "DAM_UCell",
                                         homeo_col = "Homeostatic_UCell",
                                         unit_col = "genotype_batch",
                                         genotype_col = "genotype",
                                         concordance_floor = 0.5,
                                         seed = 42L) {
  stopifnot(
    inherits(seurat_obj, "Seurat"),
    reduction %in% SeuratObject::Reductions(seurat_obj),
    all(c(substate_col, dam_col, homeo_col, unit_col, genotype_col) %in%
          colnames(seurat_obj@meta.data)),
    is.numeric(primary_dims), length(primary_dims) == 1L,
    all(is.finite(sensitivity_dims)))
  md <- seurat_obj@meta.data
  emb_full <- SeuratObject::Embeddings(seurat_obj, reduction)
  stopifnot(identical(rownames(emb_full), rownames(md)),
            ncol(emb_full) >= max(c(primary_dims, sensitivity_dims)))

  sub <- as.character(md[[substate_col]])
  dam <- as.numeric(md[[dam_col]]); homeo <- as.numeric(md[[homeo_col]])
  on_lineage <- sub %in% lineage_substates
  stopifnot(sum(on_lineage) >= 2L,
            root_substate %in% sub[on_lineage], terminal_substate %in% sub[on_lineage])

  # PRIMARY: forced single Homeostatic->DAM lineage on the on-lineage cells, dims 1:primary.
  emb_on <- emb_full[on_lineage, seq_len(primary_dims), drop = FALSE]
  prim <- run_slingshot_lineage(emb_on, sub[on_lineage],
                                root_substate, terminal_substate, seed)
  pt_raw <- stats::setNames(rep(NA_real_, nrow(md)), rownames(md))
  pt_raw[names(prim$pt)] <- prim$pt

  # marker score-axis (all retained cells) + slingshot-vs-score-axis concordance (on-lineage).
  score_axis_pt <- score_axis_pseudotime(dam, homeo)
  conc <- trajectory_concordance(pt_raw, score_axis_pt)

  # pt01 squeeze (on-lineage cells only -> off-lineage stay NA).
  pt01 <- squeeze_unit_interval(pt_raw)

  # SENSITIVITY: alt dims re-fit (Spearman vs primary on shared on-lineage cells) + an
  # all-cells-retained re-fit (IFN/prolif present -> may branch; take the DAM-terminal
  # lineage) to gauge lineage-conditioning selection effects.
  sens_rows <- list()
  for (d in sensitivity_dims) {
    alt <- run_slingshot_lineage(emb_full[on_lineage, seq_len(d), drop = FALSE],
                                 sub[on_lineage], root_substate, terminal_substate, seed)
    ok <- is.finite(alt$pt) & is.finite(prim$pt[names(alt$pt)])
    sens_rows[[length(sens_rows) + 1L]] <- data.frame(
      variant = paste0("dims_", d),
      spearman_vs_primary = stats::cor(prim$pt[names(alt$pt)][ok], alt$pt[ok],
                                       method = "spearman"),
      n_lineages = alt$n_lineages, n_shared = sum(ok), stringsAsFactors = FALSE)
  }
  allret <- run_slingshot_lineage(emb_full[, seq_len(primary_dims), drop = FALSE],
                                  sub, root_substate, terminal_substate, seed)
  shared <- intersect(names(prim$pt), names(allret$pt)[is.finite(allret$pt)])
  sens_rows[[length(sens_rows) + 1L]] <- data.frame(
    variant = "all_retained", n_lineages = allret$n_lineages, n_shared = length(shared),
    spearman_vs_primary = stats::cor(prim$pt[shared], allret$pt[shared], method = "spearman"),
    stringsAsFactors = FALSE)
  sensitivity <- do.call(rbind, lapply(sens_rows,
    function(r) r[c("variant", "spearman_vs_primary", "n_lineages", "n_shared")]))

  # per-cell compact frame (all retained cells; off-lineage pt NA, score-axis always defined).
  cell_frame <- data.frame(
    cell = rownames(md),
    genotype_batch = as.character(md[[unit_col]]),
    genotype = factor(as.character(md[[genotype_col]]), levels = genotype_levels),
    substate = factor(sub, levels = sort(unique(sub))),
    on_lineage = on_lineage,
    pt_raw = pt_raw, pt01 = pt01, score_axis_pt = score_axis_pt,
    DAM_UCell = dam, Homeostatic_UCell = homeo,
    row.names = NULL, stringsAsFactors = FALSE, check.names = FALSE)

  # per-unit on-lineage / omitted fraction (16 genotype_batch units; conditioning-selection
  # audit -- a genotype-skewed omitted fraction would bias the summaries).
  u <- factor(cell_frame$genotype_batch)
  per_unit <- data.frame(
    genotype_batch = levels(u),
    genotype = vapply(levels(u), function(g) as.character(cell_frame$genotype[u == g][1]),
                      character(1)),
    n_cells = as.integer(table(u)),
    n_on_lineage = as.integer(tapply(cell_frame$on_lineage, u, sum)),
    row.names = NULL, stringsAsFactors = FALSE)
  per_unit$omitted_frac <- 1 - per_unit$n_on_lineage / per_unit$n_cells

  prov <- trajectory_provenance(seed)
  prov$primary_dims <- primary_dims
  prov$sensitivity_dims <- sensitivity_dims
  prov$lineage_substates <- lineage_substates
  prov$root_substate <- root_substate
  prov$terminal_substate <- terminal_substate
  prov$reduction <- reduction
  prov$concordance_rho <- conc$rho
  prov$concordance_n <- conc$n
  # gross-failure floor (0.5 = conventional large-correlation boundary), RECORDED not gated:
  # rho below it (or negative) means the transcriptome ordering disagrees with marker biology
  # -> trajectory suspect. A moderate rho is EXPECTED (slingshot = transcriptome geometry vs
  # score-axis = marker contrast; they are related, not identical), so this is an honest flag.
  prov$concordance_floor <- concordance_floor
  prov$concordant <- isTRUE(conc$rho >= concordance_floor)
  # directional rooting validation (recorded): DAM up, Homeostatic down along pseudotime.
  on0 <- is.finite(pt_raw)
  prov$dam_pt_rho <- stats::cor(pt_raw[on0], dam[on0], method = "spearman")
  prov$homeo_pt_rho <- stats::cor(pt_raw[on0], homeo[on0], method = "spearman")
  prov$n_on_lineage <- sum(on_lineage)
  prov$n_omitted <- sum(!on_lineage)
  prov$omitted_frac_overall <- mean(!on_lineage)

  lineage <- list(primary = prim$lineage, primary_n_lineages = prim$n_lineages,
                  all_retained = allret$all_lineages,
                  all_retained_dam_lineage = allret$lineage)

  # build-time postconditions (the gate's warn=2 unit tests skip this heavy body -> the S1
  # acceptance INVARIANTS are pinned here; a silent regression fails tar_make, gate red):
  on <- is.finite(pt_raw)
  stopifnot(
    prim$n_lineages == 1L,                                   # single clean primary lineage
    identical(prim$lineage, c(root_substate, terminal_substate)),
    sum(on) == sum(on_lineage), !anyNA(pt_raw[on_lineage]),  # every on-lineage cell ordered
    stats::cor(pt_raw[on], dam[on], method = "spearman") > 0, # DAM monotone in pt (rooting OK)
    all(pt01[on] > 0 & pt01[on] < 1),                        # squeeze opened the interval
    is.finite(conc$rho),                                     # concordance recorded
    nrow(cell_frame) == nrow(md), nrow(per_unit) == nlevels(u),
    all(is.finite(sensitivity$spearman_vs_primary)))        # every sensitivity recorded

  list(cell_frame = cell_frame, per_unit = per_unit,
       lineage = lineage, sensitivity = sensitivity, provenance = prov)
}

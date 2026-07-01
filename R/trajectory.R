# Activation-trajectory pseudotime (P2-S1). Order microglia along the homeostatic->DAM
# axis with slingshot on the BATCH-CORRECTED harmony embedding, seeded by the substate
# biology. PRIMARY lineage = a forced single clean Homeostatic->DAM curve (clusterLabels =
# the 2 substate super-clusters -> MST on 2 nodes = one edge = one lineage; principal curve
# still fits the full high-dim cell cloud). IFN/Proliferative cells are OMITTED from the
# lineage (off-lineage flag + NA pt), not deleted -> per-unit omitted fraction is reported,
# never hidden. Pseudotime is an ACTIVATION ORDERING (position/extent of advance), NOT
# developmental time or potency -- direction rests on slingshot rooting + DAM-marker
# rank association, validated post-hoc, never a potency claim. All non-base calls namespace-
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
# GLMM (P2-S3). Two steps: (1) min-max scale finite values to [0,1]; (2) y' = (y*(n-1)+0.5)/n
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
  stopifnot(sum(ok) >= 2L,                             # >=2 finite pairs for a correlation
            length(unique(pt[ok])) >= 2L,              # constant pt or score -> Spearman undefined
            length(unique(score_axis[ok])) >= 2L)      # (fail loud, not a silent NA rho)
  list(rho = stats::cor(pt[ok], score_axis[ok], method = "spearman"), n = sum(ok))
}

# Fit slingshot on a reduced-dim embedding + cluster labels, return the pseudotime of the
# lineage TERMINATING in `terminal_clus` (Homeostatic->DAM). matrix in, list out -> pure +
# unit-testable on a synthetic embedding. 2-cluster labels guarantee a single lineage; >2
# (all-retained sensitivity, IFN present) may branch -> the DAM-terminal lineage is selected
# (longest if tied) and off-lineage cells keep slingshot's NA. Rooted at `start_clus` ->
# pseudotime increases away from the homeostatic root. RNG seed + all three kinds pinned for
# the fit, and the caller's RNG state restored on exit (pure); the curve is RNG-free given inputs.
run_slingshot_lineage <- function(embedding, cluster_labels, start_clus, terminal_clus,
                                   seed = 42L) {
  stopifnot(is.matrix(embedding), !is.null(rownames(embedding)), ncol(embedding) >= 2L,
            all(is.finite(embedding)),
            length(cluster_labels) == nrow(embedding),
            start_clus %in% cluster_labels, terminal_clus %in% cluster_labels,
            all(table(as.character(cluster_labels)) >= 2L)) # >=2 cells/cluster (covariance; singular-cov gotcha)
  # Pin seed + all three RNG kinds for the fit, then RESTORE the caller's RNG state on exit
  # (kind + .Random.seed) so the helper leaves no random-stream side effect downstream.
  old_kind <- RNGkind()
  has_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (has_seed) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    RNGkind(old_kind[1], old_kind[2], old_kind[3])
    if (has_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion",
           sample.kind = "Rejection")
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

# Fail-loud audit of per-unit metadata before the per-unit omitted-fraction summary: every
# cell needs a non-missing replicate unit + a genotype in the declared levels, and each unit
# must carry EXACTLY ONE genotype (the summary takes the first genotype per unit, so mixed
# metadata would silently corrupt the audit). Pure + data.frame-testable (no Seurat needed).
validate_trajectory_units <- function(unit_vec, geno_vec, levels = genotype_levels) {
  stopifnot(is.character(unit_vec), is.character(geno_vec),
            length(unit_vec) == length(geno_vec),
            !anyNA(unit_vec), all(nzchar(unit_vec)),               # no missing/blank unit id
            !anyNA(geno_vec), all(geno_vec %in% levels),          # genotype in declared levels
            all(tapply(geno_vec, unit_vec,                        # one genotype per unit
                       function(g) length(unique(g))) == 1L))
  invisible(TRUE)
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
    primary_dims >= 2L, primary_dims == round(primary_dims),       # positive whole-number dims
    is.numeric(sensitivity_dims), length(sensitivity_dims) >= 1L,
    all(is.finite(sensitivity_dims)), all(sensitivity_dims >= 2L),
    all(sensitivity_dims == round(sensitivity_dims)))
  md <- seurat_obj@meta.data
  emb_full <- SeuratObject::Embeddings(seurat_obj, reduction)
  stopifnot(identical(rownames(emb_full), rownames(md)),
            ncol(emb_full) >= max(c(primary_dims, sensitivity_dims)),
            is.numeric(md[[dam_col]]), is.numeric(md[[homeo_col]])) # score cols numeric, not factor

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

  # fail-loud per-unit metadata audit before the frame (mixed/missing unit or off-levels
  # genotype would silently corrupt the per-unit omitted-fraction summary below).
  unit_vec <- as.character(md[[unit_col]]); geno_vec <- as.character(md[[genotype_col]])
  validate_trajectory_units(unit_vec, geno_vec)

  # per-cell compact frame (all retained cells; off-lineage pt NA, score-axis always defined).
  cell_frame <- data.frame(
    cell = rownames(md),
    genotype_batch = unit_vec,
    genotype = factor(geno_vec, levels = genotype_levels),
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
  # directional rooting (recorded here, GATED in the postconditions): DAM up, Homeostatic down.
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
    prov$dam_pt_rho > 0, prov$homeo_pt_rho < 0,              # DAM rises, Homeo falls along pt (rooting)
    all(pt01[on] > 0 & pt01[on] < 1),                        # squeeze opened the interval
    is.finite(conc$rho),                                     # concordance recorded
    nrow(cell_frame) == nrow(md), nrow(per_unit) == nlevels(u),
    sum(per_unit$n_cells) == nrow(cell_frame),              # per-unit partition covers all cells
    all(is.finite(sensitivity$spearman_vs_primary)))        # every sensitivity recorded

  list(cell_frame = cell_frame, per_unit = per_unit,
       lineage = lineage, sensitivity = sensitivity, provenance = prov)
}

# ============================================================================================
# P2-S2a: per-replicate summary + contrast fit + Kitagawa decomposition (estimation core).
# Collapse per-cell pseudotime to the 16 genotype_batch replicate summaries, then run the 2x2
# factorial interaction as the TRAJECTORY ANALOGUE of the pseudobulk DE -- reuse factorial_design
# + its 5 canonical contrasts (ordinary t at 16 - 7 = 9 residual df, NO eBayes: a handful of
# heterogeneous pseudotime endpoints, not many shrinkable features). The 3-channel Kitagawa /
# Oaxaca shift-share (composition / progression / cross) splits the interaction EXACTLY on the
# additive pt_raw scale. Pure + deterministic-fixture-testable; NO target here (S2b wires it).
# All non-base calls namespace-qualified.

# ONE within-state column sanitizer: S2a builds the within_<state> columns, S2b reads them back
# by the SAME rule -> never hand-case the state label (a verbatim-cased lookup silently misses).
within_state_col <- function(state) paste0("within_", tolower(state))

# Extract the batch label from a genotype_batch id, VECTORISED over its length. Genotypes
# themselves contain "_" (e.g. NLGF_P301S) -> a naive strsplit/regex on "_" mis-splits, and base
# sub() does NOT vectorise over a length>1 pattern (warns + uses pattern[1] -> wrong under
# warn=2). Strip the per-element literal "<genotype>_" prefix instead, then round-trip ASSERT.
derive_batch <- function(genotype_batch, genotype) {
  stopifnot(is.character(genotype_batch),
            length(genotype_batch) == length(genotype))
  geno   <- as.character(genotype)
  prefix <- paste0(geno, "_")                              # per-element literal prefix
  stopifnot(all(startsWith(genotype_batch, prefix)))
  batch <- substring(genotype_batch, nchar(prefix) + 1L)   # vectorised over text + start
  stopifnot(all(nzchar(batch)),
            identical(paste(geno, batch, sep = "_"), genotype_batch))  # round-trip, fail loud
  batch
}

# Collapse on-lineage per-cell pseudotime to per-replicate (genotype_batch) summaries + the
# state x unit composition (pi) and within-state mean (mu) matrices the Kitagawa decomposition
# consumes. Filters to on-lineage cells (finite pt_raw); asserts every cell is in a lineage state
# + the per-unit metadata is clean. dam_onset (pooled DAM median) is PRE-DECLARED -> frac_past =
# the only genuinely [0,1] measure. within_<state> per-unit means stay on the raw pt scale
# (UNtransformed) so they share additive units with the Kitagawa progression channel. within_skip
# flags a state whose smallest unit count is below the floor (S2b drops its within-state measure).
# Pure.
pseudotime_per_replicate <- function(cell_frame, lineage_states, dam_state = "DAM",
                                     min_within = 10L) {
  stopifnot(is.data.frame(cell_frame),
            all(c("genotype_batch", "genotype", "substate", "pt_raw") %in% names(cell_frame)),
            is.character(lineage_states), dam_state %in% lineage_states,
            "min_within must be a single integer >= 2 (within-state sd needs >=2 cells/unit)" =
              (length(min_within) == 1L && is.finite(min_within) && min_within >= 2L))
  cf   <- cell_frame[is.finite(cell_frame$pt_raw), , drop = FALSE]   # on-lineage cells
  stopifnot("no on-lineage cells: every pt_raw is non-finite" = nrow(cf) > 0L)
  sub  <- as.character(cf$substate)
  unit <- as.character(cf$genotype_batch)
  geno <- as.character(cf$genotype)
  stopifnot(all(sub %in% lineage_states))                  # on-lineage cells are all lineage states
  validate_trajectory_units(unit, geno)                    # one geno/unit, geno in levels, no blank unit

  states <- lineage_states
  units  <- sort(unique(unit), method = "radix")           # locale-independent reproducible order
  f_state <- factor(sub,  levels = states)
  f_unit  <- factor(unit, levels = units)

  cnt <- table(f_state, f_unit)                            # state x unit counts
  cnt <- matrix(as.integer(cnt), nrow = length(states), ncol = length(units),
                dimnames = list(states, units))
  stopifnot(all(cnt >= 1L))                                # state present in every unit (Kitagawa needs mu_su)

  pi_mat <- sweep(cnt, 2, colSums(cnt), "/")               # column-normalised composition
  mu_mat <- tapply(cf$pt_raw, list(f_state, f_unit), mean)
  sd_mat <- tapply(cf$pt_raw, list(f_state, f_unit), stats::sd)
  dimnames(mu_mat) <- list(states, units)                  # pin identical dimnames (Kitagawa asserts it)
  dimnames(sd_mat) <- list(states, units)
  pi_bar <- rowSums(cnt) / sum(cnt)                        # CELL-weighted pooled composition
  # per-state pooled mean (states order). as.numeric STRIPS tapply's 1-D `dim` attribute: a 1-D
  # array * matrix is "non-conformable" (both carry dims), whereas a plain vector recycles.
  mu_bar <- stats::setNames(as.numeric(tapply(cf$pt_raw, f_state, mean)[states]), states)
  within_skip <- apply(cnt, 1, function(r) any(r < min_within))

  dam_onset <- stats::median(cf$pt_raw[sub == dam_state])  # PRE-DECLARED progression landmark

  geno_by_unit  <- vapply(units, function(u) geno[unit == u][1L], character(1))
  per_unit <- data.frame(
    genotype_batch = units,
    genotype  = geno_by_unit,
    batch     = derive_batch(units, geno_by_unit),
    n_cells   = as.integer(colSums(cnt)),
    sd_pt     = as.numeric(tapply(cf$pt_raw, f_unit, stats::sd)),
    mean_pt   = as.numeric(tapply(cf$pt_raw, f_unit, mean)),
    median_pt = as.numeric(tapply(cf$pt_raw, f_unit, stats::median)),
    q90       = as.numeric(tapply(cf$pt_raw, f_unit,
                                  function(v) stats::quantile(v, 0.9, names = FALSE))),
    frac_past = as.numeric(tapply(cf$pt_raw, f_unit, function(v) mean(v > dam_onset))),
    row.names = NULL, stringsAsFactors = FALSE, check.names = FALSE)
  for (s in states) per_unit[[within_state_col(s)]] <- as.numeric(mu_mat[s, units])

  list(per_unit = per_unit, states = states, units = units, counts = cnt,
       pi = pi_mat, mu = mu_mat, sd = sd_mat, pi_bar = pi_bar, mu_bar = mu_bar,
       dam_onset = dam_onset, within_skip = within_skip, min_within = min_within)
}

# Ordinary (UN-moderated) t table for one contrast of a contrasts.fit'd limma fit. topTable
# REQUIRES eBayes -> compute the ordinary t by hand (se = sigma * stdev.unscaled, df = residual
# df). For a handful of heterogeneous pseudotime endpoints the moderated/shrunk variance is
# incoherent -> plain OLS t at 9 df. Returns one row per measure (fit feature). Pure.
ordinary_t_table <- function(fit, contrast_name, conf_level = 0.95) {
  stopifnot(contrast_name %in% colnames(fit$coefficients),
            length(conf_level) == 1L, conf_level > 0, conf_level < 1)
  coef <- fit$coefficients[, contrast_name]
  se   <- fit$sigma * fit$stdev.unscaled[, contrast_name]
  df   <- fit$df.residual
  t    <- coef / se
  q    <- stats::qt(1 - (1 - conf_level) / 2, df)
  data.frame(
    measure = rownames(fit$coefficients), contrast = contrast_name,
    coef = coef, se = se, t = t, df = df,
    p_value = 2 * stats::pt(-abs(t), df), ci_l = coef - q * se, ci_r = coef + q * se,
    row.names = NULL, stringsAsFactors = FALSE)
}

# Fit a measures x units matrix against the factorial design + 5 contrasts (ordinary t, NO
# eBayes). limma consumes weights BY POSITION -> assert dimnames match (not just dims) to catch
# row/unit drift. Returns list(fit = contrasts.fit object, top = NAMED-by-contrast list of
# ordinary_t_tables) so top$interaction AND fit$coefficients[, "interaction"] both resolve
# downstream. Do NOT route through fit_limma_log (limma-TREND + eBayes, incoherent here). Pure.
fit_trajectory_contrasts <- function(measure_mat, design, contrasts, weights = NULL,
                                     conf_level = 0.95) {
  stopifnot(is.matrix(measure_mat), !is.null(rownames(measure_mat)),
            !is.null(colnames(measure_mat)),
            identical(colnames(measure_mat), rownames(design)),
            identical(rownames(contrasts), colnames(design)),
            qr(design)$rank == ncol(design),
            all(is.finite(measure_mat)))
  if (!is.null(weights)) {
    stopifnot(is.matrix(weights),
              identical(dimnames(weights), dimnames(measure_mat)),  # weights apply BY POSITION
              all(is.finite(weights)), all(weights > 0))
  }
  cfit <- limma::contrasts.fit(
    limma::lmFit(measure_mat, design = design, weights = weights), contrasts)
  # limma drops the feature rowname on a single-row fit -> restore the measure labels (1:1 with
  # the input rows; lmFit preserves feature count + order) so ordinary_t_table keys them.
  rownames(cfit$coefficients)    <- rownames(measure_mat)
  rownames(cfit$stdev.unscaled)  <- rownames(measure_mat)
  if (!is.null(weights)) {
    # limma::contrasts.fit derives a multi-coefficient contrast's SE from the UNWEIGHTED coef
    # correlation (cov.coefficients = (X'X)^-1) -> under per-feature weights it is exact ONLY for a
    # single-coef contrast or a balanced design (verified to diverge ~13% on an unbalanced one).
    # Override stdev.unscaled with the EXACT per-feature weighted normal-equation value so the
    # ordinary t holds for EVERY contrast independent of balance (the weighted-lmFit coefficients +
    # sigma are already exact WLS). se(contrast c, feature i) = sigma_i * sqrt(c' (X'W_iX)^-1 c).
    su <- vapply(seq_len(nrow(measure_mat)), function(i) {
      xtwx_inv <- chol2inv(chol(crossprod(design, design * weights[i, ])))
      sqrt(diag(crossprod(contrasts, xtwx_inv %*% contrasts)))
    }, numeric(ncol(contrasts)))
    cfit$stdev.unscaled <- matrix(su, nrow(measure_mat), ncol(contrasts), byrow = TRUE,
                                  dimnames = list(rownames(measure_mat), colnames(contrasts)))
  }
  top <- stats::setNames(
    lapply(colnames(contrasts), function(cn) ordinary_t_table(cfit, cn, conf_level)),
    colnames(contrasts))
  list(fit = cfit, top = top)
}

# EXACT 3-channel Kitagawa / Oaxaca shift-share of the per-unit mean pseudotime:
#   mean_pt_r = comp_cf_r + prog_cf_r + cross_r - const, where (state s, unit r)
#   comp_cf = sum_s pi_sr * mu_bar_s  (vary composition, hold within-state means pooled)
#   prog_cf = sum_s pi_bar_s * mu_sr  (vary within-state means, hold composition pooled)
#   cross   = sum_s (pi_sr - pi_bar_s)(mu_sr - mu_bar_s),   const = sum_s pi_bar_s * mu_bar_s.
# Every broadcast (pi_bar/mu_bar down columns) is POSITION-based -> assert the state order is
# identical across pi/mu/pi_bar/mu_bar FIRST (a silent order drift would mis-split). Holds on the
# RAW additive pt scale only (logit/asin break additivity). Returns a per-unit data.frame. Pure.
kitagawa_channels <- function(pi_mat, mu_mat, pi_bar, mu_bar, tol = 1e-8) {
  stopifnot(is.matrix(pi_mat), is.matrix(mu_mat),
            identical(dimnames(pi_mat), dimnames(mu_mat)),
            identical(rownames(pi_mat), names(pi_bar)),
            identical(names(pi_bar), names(mu_bar)),
            all(is.finite(pi_mat)), all(is.finite(mu_mat)),
            all(is.finite(pi_bar)), all(is.finite(mu_bar)),
            all(pi_mat >= 0), all(pi_bar >= 0),                  # pi is a genuine composition ...
            max(abs(colSums(pi_mat) - 1)) < tol,                 # ... columns sum to 1 ...
            abs(sum(pi_bar) - 1) < tol)                          # ... as does the pooled anchor
  const   <- sum(pi_bar * mu_bar)
  mean_pt <- colSums(pi_mat * mu_mat)
  comp_cf <- colSums(pi_mat * mu_bar)
  prog_cf <- colSums(pi_bar * mu_mat)
  cross   <- colSums((pi_mat - pi_bar) * (mu_mat - mu_bar))
  stopifnot(max(abs(mean_pt - (comp_cf + prog_cf + cross - const))) < tol)  # reconstruction
  data.frame(
    genotype_batch = colnames(pi_mat),
    mean_pt = mean_pt, comp_cf = comp_cf, prog_cf = prog_cf, cross = cross, const = const,
    row.names = NULL, stringsAsFactors = FALSE)
}

# Decompose the 2x2 interaction on mean pseudotime into its composition / progression / cross
# channels. The interaction contrast L() is LINEAR + intercept-free -> L(mean_pt) = L(comp_cf) +
# L(prog_cf) + L(cross) EXACTLY (const is unit-constant, annihilated). The exactness needs ONE
# shared per-unit weight vector replicated across the 4 channel-rows (the SAME WLS operator hits
# each row); differing per-row weights generally break it. loadings = the 3 channel
# interaction coefs / L(mean_pt) (NA if |L(mean_pt)| < tol). Cell-weighted pooled anchors =
# PRIMARY; replicate-balanced (rowMeans) anchors = a SENSITIVITY. Pure.
decompose_progression_vs_composition <- function(per_rep, design, contrasts,
                                                 weights = NULL, conf_level = 0.95) {
  tol <- 1e-8
  channel_rows <- c("mean_pt", "comp_cf", "prog_cf", "cross")
  fit_channels <- function(pi_bar, mu_bar) {
    ch <- kitagawa_channels(per_rep$pi, per_rep$mu, pi_bar, mu_bar)
    M4 <- t(as.matrix(ch[, channel_rows]))
    dimnames(M4) <- list(channel_rows, ch$genotype_batch)
    stopifnot(identical(colnames(M4), rownames(design)))
    W4 <- NULL
    if (!is.null(weights)) {
      stopifnot(is.numeric(weights), !is.null(names(weights)),
                all(colnames(M4) %in% names(weights)))
      wv <- weights[colnames(M4)]                          # index by unit -> M4 column order
      stopifnot(all(is.finite(wv)), all(wv > 0))
      W4 <- matrix(wv, nrow(M4), ncol(M4), byrow = TRUE, dimnames = dimnames(M4))
    }
    fit   <- fit_trajectory_contrasts(M4, design, contrasts, weights = W4, conf_level = conf_level)
    L_int <- fit$fit$coefficients[, "interaction"]         # named over channel-rows
    recon <- unname(abs(L_int["mean_pt"] -
                        (L_int["comp_cf"] + L_int["prog_cf"] + L_int["cross"])))
    loadings <- if (abs(L_int["mean_pt"]) < tol)
      stats::setNames(rep(NA_real_, 3L), c("comp_cf", "prog_cf", "cross"))
    else L_int[c("comp_cf", "prog_cf", "cross")] / unname(L_int["mean_pt"])
    list(channels = ch, fit = fit, L_int = L_int, loadings = loadings, recon_resid_max = recon)
  }
  primary  <- fit_channels(per_rep$pi_bar, per_rep$mu_bar)              # cell-weighted anchors
  balanced <- fit_channels(rowMeans(per_rep$pi), rowMeans(per_rep$mu)) # replicate-balanced anchors
  stopifnot(primary$recon_resid_max < tol, balanced$recon_resid_max < tol)
  list(channels = primary$channels, fit = primary$fit, L_int = primary$L_int,
       loadings = primary$loadings, interaction = primary$fit$top$interaction,
       recon_resid_max = primary$recon_resid_max, balanced = balanced)
}

# ============================================================================================
# P2-S2b: progression-interaction inference + orchestrator + the trajectory_progression target.
# The weighted-limma per-replicate-summary interaction (S2a) is the PRIMARY inference; S2b adds a
# Freedman-Lane permutation null (a distribution-light SENSITIVITY) and run_trajectory_progression,
# which wires S1's compact target -> per-replicate summary -> factorial design -> weighted / OLS /
# bounded contrast fits + the 3-channel Kitagawa decomposition + the permutation null, under a
# PRE-REGISTERED primary BH family {progression_cf, within_homeostatic}. Pure-R, NO new dependency
# (glmmTMB is S3). All non-base calls namespace-qualified.

# Freedman-Lane permutation null for the 2x2 interaction coefficient (int_col, default tau_nlgf) of
# a weighted least-squares fit. WLS = OLS on weight-scaled data (row-scale y + design by
# sqrt(weights)); the interaction t is recomputed under each permutation of the REDUCED-model
# (interaction-dropped) WEIGHTED residuals added back to the reduced fit. The weights are ESTIMATED
# from the same per-unit summaries -> exchangeability on the weighted scale is approximate /
# conditional, NOT exact -> a SENSITIVITY cross-check of the parametric ordinary-t, never a
# nominal-exact primary (permuting raw unweighted residuals would be worse still). Xw is fixed
# across permutations -> precompute the pivot-FREE (X'X)^-1 ONCE (its index order matches lm.fit's
# coef order for the full-rank design, sidestepping qr() column-pivoting). RNG-pure: pins the seed +
# all three kinds for the permutations, restores the caller's stream on exit. Returns
# list(t_obs, n_perm, perm_p). Pure.
freedman_lane_interaction <- function(y, design, int_col = "tau_nlgf", weights = NULL,
                                      n_perm = 2000L, seed = 42L) {
  stopifnot(is.numeric(y), is.matrix(design), !is.null(colnames(design)),
            length(y) == nrow(design), all(is.finite(y)), all(is.finite(design)),
            int_col %in% colnames(design), qr(design)$rank == ncol(design),
            length(n_perm) == 1L, is.finite(n_perm), n_perm >= 1L, n_perm == round(n_perm),
            length(seed) == 1L, is.finite(seed), seed == round(seed))
  n <- length(y)
  if (!is.null(weights)) {
    stopifnot(is.numeric(weights), length(weights) == n,
              all(is.finite(weights)), all(weights > 0))
  }
  r  <- if (is.null(weights)) rep(1, n) else sqrt(weights)
  yw <- r * y
  Xw <- r * design                                       # row-scale by sqrt(weights) -> WLS as OLS
  XtXinv <- chol2inv(chol(crossprod(Xw)))                # pivot-FREE (X'X)^-1: index order == lm.fit coefs
  stopifnot(all(is.finite(XtXinv)))
  j  <- match(int_col, colnames(design))
  p  <- ncol(Xw); df <- n - p
  stopifnot(df >= 1L)
  int_t <- function(yv) {                                # full-model WLS interaction t (Xw fixed)
    f <- stats::lm.fit(Xw, yv)
    sigma2 <- sum(f$residuals^2) / df
    unname(f$coefficients[j] / sqrt(sigma2 * XtXinv[j, j]))
  }
  t_obs <- int_t(yw)
  Xw0 <- Xw[, colnames(Xw) != int_col, drop = FALSE]     # reduced (interaction-dropped) weighted fit
  f0  <- stats::lm.fit(Xw0, yw)
  fw0 <- f0$fitted.values; ew0 <- f0$residuals           # exchangeable residuals on the weighted scale
  # RNG-pure: pin seed + all three kinds for the permutations; restore the caller's stream on exit
  # (RNGkind FIRST then assign .Random.seed -- RNGkind reinitialises the stream).
  old_kind <- RNGkind()
  has_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (has_seed) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    RNGkind(old_kind[1], old_kind[2], old_kind[3])
    if (has_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion",
           sample.kind = "Rejection")
  t_star <- vapply(seq_len(n_perm), function(b) int_t(fw0 + ew0[sample.int(n)]), numeric(1))
  list(t_obs = t_obs, n_perm = as.integer(n_perm),
       perm_p = (1 + sum(abs(t_star) >= abs(t_obs))) / (n_perm + 1))
}

# Orchestrate the progression-interaction inference from the COMPACT microglia_trajectory target
# (S1). Collapses on-lineage per-cell pseudotime to the 16 genotype_batch summaries, builds the
# factorial design (9 residual df), and fits the interaction on every endpoint three ways:
#   - weighted: direct measures {mean_pt, median_pt, q90, within_<used>}, inverse-summary-variance
#     weights (mean/median/q90 = n/sd_pt^2 overall precision; within_<state> = state n/sd^2);
#   - ols: the same direct measures unweighted (a weighting-sensitivity);
#   - bounded: frac_past on the logit + asin VST bridges (weights = n_cells; EXPLORATORY).
# Plus the exact 3-channel Kitagawa decomposition (cell-weighted anchors) and a Freedman-Lane
# permutation null on {progression_cf, within_homeostatic, frac_past_logit, mean_pt}, each weighted
# to MATCH its limma fit. PRE-REGISTERED: primary family = BH across {progression_cf,
# within_homeostatic} (composition-robust); everything else is a SEPARATE exploratory BH, and mean_pt
# is flagged composition-conflated. Reads the compact target; pure-R (NO glmmTMB -> S3). Returns the
# per-unit summary, the three fits, the decomposition, the permutation null, the two BH families, and
# provenance (incl. the v1 progression loading ~0.94 / fdr ~0.077 for honest reconciliation).
run_trajectory_progression <- function(microglia_trajectory, min_within = 10L,
                                       n_perm = 2000L, seed = 42L) {
  prov_in        <- microglia_trajectory$provenance
  lineage_states <- prov_in$lineage_substates
  dam_state      <- prov_in$terminal_substate
  root_state     <- prov_in$root_substate
  stopifnot(is.list(microglia_trajectory), "cell_frame" %in% names(microglia_trajectory),
            is.character(lineage_states), length(lineage_states) >= 2L,
            dam_state %in% lineage_states, root_state %in% lineage_states)

  per_rep  <- pseudotime_per_replicate(microglia_trajectory$cell_frame, lineage_states,
                                       dam_state = dam_state, min_within = min_within)
  per_unit <- per_rep$per_unit
  gb       <- per_unit$genotype_batch
  within_skip <- per_rep$within_skip
  used_states <- per_rep$states[!within_skip[per_rep$states]]

  # factorial design over the 16 units (rownames = genotype_batch -> match the measure columns).
  meta <- per_unit[, c("genotype_batch", "genotype", "batch")]
  rownames(meta) <- gb
  assert_complete_crossing(meta, "genotype_batch")            # 4x4 balance, fail loud
  fd <- factorial_design(meta)
  design <- fd$design; contrasts <- fd$contrasts
  units_order <- rownames(design)
  # assert_complete_crossing only checks n_units == prod(OBSERVED levels) -> it passes a complete
  # SUB-rectangle (a dropped batch -> 4x3 = 12 units / 6 resid df). Enforce the DOCUMENTED full 4x4
  # design here: exactly one unit per genotype x batch cell + the 16-unit / 9-resid-df shape. An
  # upstream unit/batch loss (e.g. a genotype_batch with zero on-lineage cells) then fails loud, not
  # silently as a lower-power fit that contradicts the 9-df claim.
  stopifnot("expected the full 16-unit genotype x batch design (9 residual df)" =
              (nrow(design) == 16L && nrow(design) - ncol(design) == 9L),
            "genotype x batch must be crossed exactly once per cell" =
              all(table(meta$genotype, meta$batch) == 1L))

  # direct-measure matrix (measures x units) + the per-endpoint weight matrix (IDENTICAL dimnames --
  # limma applies weights BY POSITION). w_overall = n_cells/sd_pt^2 = a SHARED mean-precision weight
  # reused across mean/median/q90 (heuristic for the two quantiles) + chosen so the Kitagawa channels
  # reconstruct EXACTLY under one shared weight -> progression_cf is a shared-weight DECOMPOSITION
  # test, not channel-specific inverse-variance inference.
  within_used <- unname(vapply(used_states, within_state_col, character(1)))
  direct_rows <- c("mean_pt", "median_pt", "q90", within_used)
  M <- t(as.matrix(per_unit[, direct_rows, drop = FALSE]))
  dimnames(M) <- list(direct_rows, gb)
  w_overall <- stats::setNames(per_unit$n_cells / per_unit$sd_pt^2, gb)
  W <- matrix(NA_real_, nrow(M), ncol(M), dimnames = dimnames(M))
  for (rn in c("mean_pt", "median_pt", "q90")) W[rn, ] <- w_overall[colnames(M)]
  for (s in used_states)
    W[within_state_col(s), ] <-
      per_rep$counts[s, colnames(M)] / per_rep$sd[s, colnames(M)]^2
  within_sd_ok <- vapply(used_states, function(s)
    all(is.finite(per_rep$sd[s, ]) & per_rep$sd[s, ] > 0), logical(1))
  stopifnot("within-state pt sd must be finite + positive (>=2 distinct cells/unit); raise min_within" =
              all(within_sd_ok),
            identical(dimnames(W), dimnames(M)), all(is.finite(W)), all(W > 0))
  fit_weighted <- fit_trajectory_contrasts(M, design, contrasts, weights = W)
  fit_ols      <- fit_trajectory_contrasts(M, design, contrasts, weights = NULL)

  # BOUNDED frac_past on the logit + asin VST bridges (weights = n_cells; EXPLORATORY).
  x <- round(per_unit$frac_past * per_unit$n_cells)
  frac_logit <- log((x + 0.5) / (per_unit$n_cells - x + 0.5))
  frac_asin  <- asin(sqrt(per_unit$frac_past))
  B  <- rbind(frac_past_logit = frac_logit, frac_past_asin = frac_asin)
  colnames(B) <- gb
  WB <- matrix(per_unit$n_cells, nrow(B), ncol(B), byrow = TRUE, dimnames = dimnames(B))
  stopifnot(all(is.finite(WB)), all(WB > 0))
  fit_bounded <- fit_trajectory_contrasts(B, design, contrasts, weights = WB)

  # exact 3-channel decomposition (cell-weighted anchors primary; one shared per-unit weight vector).
  decomposition <- decompose_progression_vs_composition(per_rep, design, contrasts,
                                                        weights = w_overall)

  # Freedman-Lane null, EACH call weighted to MATCH its limma fit.
  fl <- function(y_named, w_named)
    freedman_lane_interaction(y_named[units_order], design, int_col = "tau_nlgf",
                              weights = w_named[units_order], n_perm = n_perm, seed = seed)
  perm <- list(
    mean_pt = fl(stats::setNames(per_unit$mean_pt, gb), w_overall),
    progression_cf = fl(stats::setNames(decomposition$channels$prog_cf,
                                        decomposition$channels$genotype_batch), w_overall),
    frac_past_logit = fl(stats::setNames(frac_logit, gb),
                         stats::setNames(per_unit$n_cells, gb)))
  if (root_state %in% used_states) {
    wcol <- within_state_col(root_state)
    w_within <- stats::setNames(per_rep$counts[root_state, ] / per_rep$sd[root_state, ]^2,
                                colnames(per_rep$counts))
    perm[[wcol]] <- fl(stats::setNames(per_unit[[wcol]], gb), w_within)
  }

  # assemble the per-measure interaction table from the MATCHING fit, attach perm_p, split into the
  # PRE-REGISTERED primary BH family {progression_cf, within_homeostatic} + a separate exploratory BH.
  wint <- fit_weighted$top$interaction
  dint <- decomposition$interaction
  bint <- fit_bounded$top$interaction
  grab <- function(top_int, src, out = src, perm_p = NA_real_) {
    rrow <- top_int[top_int$measure == src,
                    c("coef", "se", "t", "df", "p_value", "ci_l", "ci_r"), drop = FALSE]
    stopifnot(nrow(rrow) == 1L)
    data.frame(measure = out, rrow, perm_p = perm_p, row.names = NULL, stringsAsFactors = FALSE)
  }
  rows <- list(
    grab(dint, "prog_cf", "progression_cf", perm$progression_cf$perm_p),
    grab(wint, "mean_pt", perm_p = perm$mean_pt$perm_p),
    grab(wint, "median_pt"), grab(wint, "q90"),
    grab(dint, "comp_cf"), grab(dint, "cross"),
    grab(bint, "frac_past_logit", perm_p = perm$frac_past_logit$perm_p),
    grab(bint, "frac_past_asin"))
  for (s in used_states) {
    mc <- within_state_col(s)
    pp <- if (!is.null(perm[[mc]])) perm[[mc]]$perm_p else NA_real_
    rows[[length(rows) + 1L]] <- grab(wint, mc, perm_p = pp)
  }
  all_int <- do.call(rbind, rows)
  primary_measures <- c("progression_cf", within_state_col(root_state))
  primary_measures <- primary_measures[primary_measures %in% all_int$measure]
  is_primary <- all_int$measure %in% primary_measures
  primary_family     <- all_int[is_primary, , drop = FALSE]
  exploratory_family <- all_int[!is_primary, , drop = FALSE]
  primary_family$fdr     <- stats::p.adjust(primary_family$p_value, "BH")
  exploratory_family$fdr <- stats::p.adjust(exploratory_family$p_value, "BH")
  rownames(primary_family) <- NULL; rownames(exploratory_family) <- NULL

  provenance <- list(
    limma_version = as.character(utils::packageVersion("limma")),
    r_version = as.character(getRversion()), seed = seed, n_perm = as.integer(n_perm),
    min_within = min_within, used_states = used_states, within_skip = within_skip,
    dam_onset = per_rep$dam_onset, primary_measures = primary_measures,
    planned_primary = c("progression_cf", within_state_col(root_state)),    # the PRE-REGISTERED pair
    primary_within_skipped = !(within_state_col(root_state) %in% primary_measures),  # root too thin -> analyzable subset
    progression_loading = unname(decomposition$loadings["prog_cf"]),
    composition_loading = unname(decomposition$loadings["comp_cf"]),
    cross_loading = unname(decomposition$loadings["cross"]),
    recon_resid_max = decomposition$recon_resid_max,
    v1_progression_loading = 0.94, v1_progression_fdr = 0.077)

  # postconditions: interaction on EVERY measure, exact reconstruction, primary BH present + finite.
  stopifnot(
    all(direct_rows %in% wint$measure),
    all(c("frac_past_logit", "frac_past_asin") %in% bint$measure),
    all(c("mean_pt", "comp_cf", "prog_cf", "cross") %in% dint$measure),
    decomposition$recon_resid_max < 1e-8,
    nrow(primary_family) >= 1L, "progression_cf" %in% primary_family$measure,
    all(is.finite(primary_family$fdr)), all(is.finite(primary_family$p_value)),
    all(is.finite(c(perm$mean_pt$perm_p, perm$progression_cf$perm_p,
                    perm$frac_past_logit$perm_p))))

  list(per_unit = per_unit, counts = per_rep$counts, dam_onset = per_rep$dam_onset,
       within_skip = within_skip, design = design,
       contrasts = list(weighted = fit_weighted, ols = fit_ols, bounded = fit_bounded),
       decomposition = decomposition, permutation = perm,
       primary_family = primary_family, exploratory_family = exploratory_family,
       provenance = provenance)
}

# ============================================================================================
# P2-S3: glmmTMB per-cell pseudotime sensitivity (SUPPORTIVE arm).
# A replication-aware per-cell confirmation that models the FULL bounded (possibly bimodal)
# distribution the 16-unit summary collapses. The weighted-limma summary + Kitagawa decomposition
# (S2a/S2b) is the standalone PRIMARY; this arm is supportive at 16 clusters (asymptotics weak) and
# DEGRADES gracefully (singular RE -> rank-normal LMM -> a RECORDED method="failed"), NEVER blocking.
# On-lock: TMB = a C++ template, NOT Stan. All non-base calls namespace-qualified.
# ============================================================================================

# Evaluate `expr`, capturing + MUFFLING both warnings AND messages (the sccomp lesson, R/composition.R:
# glmmTMB/TMB optimisers can report convergence health via message() carrying a literal "Warning:" too,
# NOT only warning() -> a fresh build would RED the gate's warn=2 / tar_meta / anchored ^Warning: log
# scan). Returns list(value, warnings, messages); value = NULL when `expr` errors (a fit/extraction
# failure -> the caller degrades). The structured health flags below are the authoritative record.
.capture_quietly <- function(expr) {
  warns <- character(0); msgs <- character(0)
  value <- withCallingHandlers(
    tryCatch(expr, error = function(e) { msgs <<- c(msgs, paste0("error: ", conditionMessage(e))); NULL }),
    warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") },
    message = function(m) { msgs  <<- c(msgs,  conditionMessage(m)); invokeRestart("muffleMessage") }
  )
  list(value = value, warnings = warns, messages = msgs)
}

# Health battery for ONE fit's interaction row -> a single tested gate (P2-S3). PURE + deterministic,
# unit-tested directly so the non-convergence / degenerate-SE branches don't hinge on coaxing the
# optimiser. ok = pos-def Hessian & optimiser-converged & finite est & finite POSITIVE se (a zero SE
# gives an infinite z + zero-width CI) & finite z & a valid probability p & a non-singular RE.
.fit_health_ok <- function(pdHess, convergence, est, se, z, p, singular) {
  isTRUE(pdHess) && isTRUE(convergence == 0) &&
    is.finite(est) && is.finite(se) && se > 0 &&
    is.finite(z) && is.finite(p) && p >= 0 && p <= 1 &&
    isFALSE(singular)
}

# Fit ONE glmmTMB model + extract the tau:amyloid interaction Wald row + run the health battery,
# under .capture_quietly. Returns the capture list; $value = a standardized record
# list(method, term, estimate, se, z, p_value, ci_l, ci_r, re_sd, singular, ok) when the fit AND
# extraction succeed, else NULL (any error -> NULL -> the caller degrades). ok = .fit_health_ok()
# (pdHess & converged & finite est & finite POSITIVE se & finite z & valid p & non-singular RE) --
# the SAME tested gate for the beta GLMM and the rank-normal LMM fallback. Wald columns read by
# POSITION (glmmTMB's fixed order Estimate / Std. Error / z value / Pr(>|z|)) so a column-NAME drift
# cannot mis-extract; a positional-integrity guard (z==est/se, p==2*pnorm(-|z|)) catches the converse
# (a column-ORDER change) by degrading, which the method=="glmmTMB_beta" test then flags loud.
.fit_pt_interaction <- function(formula, family, data, method_label) {
  .capture_quietly({
    fit  <- glmmTMB::glmmTMB(formula, data = data, family = family)
    cond <- summary(fit)$coefficients$cond
    term <- intersect(c("tau:amyloid", "amyloid:tau"), rownames(cond))
    stopifnot(length(term) == 1L, ncol(cond) >= 4L)        # interaction row present + unambiguous
    row  <- cond[term, ]
    est  <- unname(row[1L]); se <- unname(row[2L]); z <- unname(row[3L]); p <- unname(row[4L])
    # positional-integrity guard (naming-agnostic): confirm cols 3,4 ARE the Wald z & p, so a future
    # summary-layout change degrades here (NULL -> caller falls back) instead of silently mis-reading.
    if (is.finite(est) && is.finite(se) && se > 0) {
      stopifnot(isTRUE(all.equal(z, est / se, tolerance = 1e-5)),
                isTRUE(all.equal(p, 2 * stats::pnorm(-abs(z)), tolerance = 1e-5)))
    }
    re   <- glmmTMB::VarCorr(fit)$cond$unit                # 1x1 variance matrix for (1|unit)
    re_sd    <- if (is.null(re)) NA_real_ else sqrt(as.numeric(re)[1])
    singular <- !is.finite(re_sd) || re_sd < 1e-4          # collapsed/unidentifiable RE -> degrade
    ok <- .fit_health_ok(fit$sdr$pdHess, fit$fit$convergence, est, se, z, p, singular)
    list(method = method_label, term = term, estimate = est, se = se, z = z, p_value = p,
         ci_l = est - stats::qnorm(0.975) * se, ci_r = est + stats::qnorm(0.975) * se,
         re_sd = re_sd, singular = singular, ok = ok)
  })
}

# Per-cell beta-GLMM sensitivity for the tau:amyloid interaction on bounded pseudotime. Reads the
# COMPACT microglia_trajectory$cell_frame (NOT the 612MB Seurat); on-lineage = finite pt01.
# tau / amyloid = integer 0/1 from genotype (matching factorial_design); batch + unit (genotype_batch)
# as factors. PRIMARY = beta_family() on pt01 (Smithson-Verkuilen-squeezed open (0,1)); DEGRADE to a
# rank-normal LMM (same package, on-lock) on any battery failure; if BOTH fits fail, RECORD
# method="failed" (NA effect) -- a fit/extraction failure NEVER throws (the supportive arm degrades +
# records fail_reason). MALFORMED INPUT (missing cols, non-finite/boundary pt01, unknown genotype,
# broken genotype_batch) fails LOUD via stopifnot -> surfaces an upstream pipeline break, not masks it.
# FIXED batch (de_pb-consistent), (1|unit) random intercept. Returns list(method, term, estimate, se,
# z, p_value, ci_l, ci_r, re_sd, singular, n_cells, n_units, fail_reason, warnings, messages) --
# n_units = genotype_batch clusters present (asymptotics basis, RECORDED not asserted); supportive,
# concordance AND discordance both fine.
glmmtmb_pt_sensitivity <- function(cell_frame, pt_col = "pt01") {
  stopifnot(is.data.frame(cell_frame),
            all(c("genotype_batch", "genotype", pt_col) %in% names(cell_frame)))
  d0   <- cell_frame[is.finite(cell_frame[[pt_col]]), , drop = FALSE]   # on-lineage cells only
  geno <- as.character(d0$genotype)
  stopifnot(all(geno %in% genotype_levels))                            # reject unknown/corrupt genotypes (else silently coded tau=0/amyloid=0)
  dat  <- data.frame(
    pt01    = as.numeric(d0[[pt_col]]),
    tau     = as.integer(geno %in% c("P301S", "NLGF_P301S")),
    amyloid = as.integer(geno %in% c("NLGF_MAPTKI", "NLGF_P301S")),
    batch   = factor(derive_batch(as.character(d0$genotype_batch), geno)),
    unit    = factor(as.character(d0$genotype_batch)),
    stringsAsFactors = FALSE)
  n_cells <- nrow(dat)
  n_units <- nlevels(dat$unit)                                         # genotype_batch clusters present (GLMM asymptotics basis)
  stopifnot(n_cells >= 1L, all(dat$pt01 > 0 & dat$pt01 < 1))            # beta needs OPEN (0,1)

  warns <- character(0); msgs <- character(0)
  collect <- function(cap) { warns <<- c(warns, cap$warnings); msgs <<- c(msgs, cap$messages); cap$value }
  done <- function(rec) list(method = rec$method, term = rec$term, estimate = rec$estimate,
                             se = rec$se, z = rec$z, p_value = rec$p_value, ci_l = rec$ci_l,
                             ci_r = rec$ci_r, re_sd = rec$re_sd, singular = rec$singular,
                             n_cells = n_cells, n_units = n_units, fail_reason = NA_character_,
                             warnings = warns, messages = msgs)

  # PRIMARY: beta GLMM on the bounded squeezed pseudotime.
  beta <- collect(.fit_pt_interaction(pt01 ~ tau * amyloid + batch + (1 | unit),
                                      glmmTMB::beta_family(), dat, "glmmTMB_beta"))
  if (!is.null(beta) && isTRUE(beta$ok)) return(done(beta))

  # DEGRADE: rank-normal LMM (on-lock, SAME package), same interaction extraction + battery.
  dat$rn <- stats::qnorm((rank(dat$pt01) - 0.5) / n_cells)
  lmm <- collect(.fit_pt_interaction(rn ~ tau * amyloid + batch + (1 | unit),
                                     stats::gaussian(), dat, "lmm_ranknorm"))
  if (!is.null(lmm) && isTRUE(lmm$ok)) return(done(lmm))

  # BOTH degraded -> RECORD a failed-supportive result (NA effect), never throw on a FIT failure;
  # carry the best-available attempt diagnostics so the failure is explained, not silent.
  reason <- function(tag, rec) if (is.null(rec)) paste0(tag, ":error")           # fit/extraction threw
                               else if (!is.finite(rec$estimate)) paste0(tag, ":nonestimable")  # interaction dropped (rank-deficient)
                               else if (isTRUE(rec$singular)) paste0(tag, ":singular")          # RE collapsed
                               else paste0(tag, ":nonconverge")                  # !pdHess / non-zero code / bad se|z|p
  last   <- if (!is.null(lmm)) lmm else beta              # most-degraded attempt that still returned a record
  list(method = "failed", term = NA_character_, estimate = NA_real_, se = NA_real_,
       z = NA_real_, p_value = NA_real_, ci_l = NA_real_, ci_r = NA_real_,
       re_sd = if (is.null(last)) NA_real_ else last$re_sd,
       singular = if (is.null(last)) NA else last$singular,
       n_cells = n_cells, n_units = n_units,
       fail_reason = paste(reason("beta", beta), reason("lmm", lmm)),
       warnings = warns, messages = msgs)
}

# ============================================================================================
# P2-S4: compact report-data extraction for the trajectory chapter (keeps the gate render cheap).
# Bundle EVERYTHING _trajectory.qmd plots/tabulates from the three trajectory targets into one
# small target, so the force-rendered report (hence EVERY scripts/check.sh run) tar_loads a single
# compact object -- never the 612MB Seurat. All three inputs are ALREADY compact (microglia_trajectory
# ~3.3MB in memory; the two inference targets are small), so no heavy object is read here. Pure: no
# RNG, no I/O. The per-cell plotting frame is asserted render-clean by construction (finite pt on
# on-lineage cells, finite score-axis, no missing genotype/substate) so the qmd never trips a ggplot
# missing-value warning under warn=2 (which would red the gate).
# ============================================================================================
trajectory_report_data <- function(microglia_trajectory, trajectory_progression,
                                   trajectory_glmm_sensitivity) {
  stopifnot(
    is.list(microglia_trajectory),
    all(c("cell_frame", "per_unit", "sensitivity", "provenance") %in% names(microglia_trajectory)),
    is.list(trajectory_progression),
    all(c("per_unit", "contrasts", "primary_family",
          "exploratory_family", "provenance") %in% names(trajectory_progression)),
    is.list(trajectory_glmm_sensitivity),
    # FULL 13-name glmm set: the [c(...)] row-subset below NA-FILLS any missing name (a silent
    # false-green that surfaces only mid-render) -> require every name the subset pulls, up front.
    all(c("method", "term", "estimate", "se", "z", "p_value", "ci_l", "ci_r", "re_sd",
          "singular", "n_cells", "n_units", "fail_reason") %in%
          names(trajectory_glmm_sensitivity)))
  tcf <- microglia_trajectory$cell_frame
  tp  <- trajectory_progression
  mp  <- microglia_trajectory$provenance
  stopifnot(
    is.data.frame(tcf),
    all(c("genotype", "substate", "on_lineage", "pt_raw", "score_axis_pt") %in% names(tcf)),
    # NESTED fields the body reads (mirror microglia_report_data: guard EVERY field the body pulls,
    # not just the top-level containers) -> a malformed input fails HERE, never as a silent NULL that
    # only breaks mid-render.
    is.list(tp$contrasts), "weighted" %in% names(tp$contrasts),
    is.list(tp$contrasts$weighted), "top" %in% names(tp$contrasts$weighted),   # weighted_top (5 contrasts)
    # per_unit + sensitivity columns the conditioning + robustness panels read:
    is.data.frame(microglia_trajectory$per_unit),
    all(c("genotype", "n_cells", "n_on_lineage", "omitted_frac") %in%
          names(microglia_trajectory$per_unit)),
    is.data.frame(microglia_trajectory$sensitivity),
    all(c("variant", "spearman_vs_primary") %in% names(microglia_trajectory$sensitivity)),
    # provenance source fields the inline prose pulls (assembled into out$provenance below):
    all(c("primary_dims", "lineage_substates", "root_substate", "terminal_substate",
          "concordance_rho", "concordance_floor", "concordant", "dam_pt_rho", "homeo_pt_rho",
          "omitted_frac_overall") %in% names(mp)),
    all(c("dam_onset", "used_states", "within_skip", "primary_measures", "planned_primary",
          "primary_within_skipped", "progression_loading", "composition_loading", "cross_loading",
          "recon_resid_max", "v1_progression_loading", "v1_progression_fdr", "n_perm", "seed",
          "limma_version", "r_version") %in% names(tp$provenance)))

  # slim per-cell plotting frame (off-lineage cells keep NA pt; the qmd filters on_lineage for the
  # pseudotime panels and uses the always-defined score-axis for the concordance panel).
  cell_frame <- data.frame(
    genotype      = tcf$genotype,                          # factor over genotype_levels (preserved)
    substate      = tcf$substate,
    on_lineage    = tcf$on_lineage,
    pt_raw        = tcf$pt_raw,
    score_axis_pt = tcf$score_axis_pt,
    row.names = NULL, stringsAsFactors = FALSE)

  # interaction across EVERY measure (primary BH family + the separate exploratory BH), tagged with
  # the family so the forest can colour/group them; each row carries coef / SE / 95% CI / ordinary-t
  # p / Freedman-Lane perm_p / BH FDR. The pt_raw-scale measures {mean_pt, comp_cf, progression_cf,
  # cross, within_*} share additive units (one forest axis); frac_past_* are transformed bridges (flagged).
  interaction <- rbind(
    cbind(family = "primary",     tp$primary_family,     stringsAsFactors = FALSE),
    cbind(family = "exploratory", tp$exploratory_family, stringsAsFactors = FALSE))
  rownames(interaction) <- NULL

  # the 3-channel decomposition is NOT re-bundled: the qmd draws its loadings figure + prose from the
  # provenance loadings (composition/progression/cross_loading, guarded finite below) and the per-channel
  # coefs from the comp_cf/progression_cf/cross rows of `interaction` above -> a `decomposition` field
  # would only duplicate those two live sources (dead figure-shaped output, codex 955).
  list(
    cell_frame   = cell_frame,
    interaction  = interaction,
    weighted_top = tp$contrasts$weighted$top,             # named-by-contrast per-measure tables (5 contrasts)
    per_unit     = tp$per_unit,                           # 16-unit summary (mean_pt, frac_past, within_*, n_cells, sd_pt)
    lineage_per_unit = microglia_trajectory$per_unit,     # per-unit n_on_lineage + omitted_frac (conditioning audit)
    sensitivity  = microglia_trajectory$sensitivity,      # dims {10,20} + all-retained robustness vs primary
    glmm = trajectory_glmm_sensitivity[c("method", "term", "estimate", "se", "z", "p_value",
                                         "ci_l", "ci_r", "re_sd", "singular", "n_cells",
                                         "n_units", "fail_reason")],  # supportive per-cell arm (warnings/messages dropped)
    provenance = list(
      # slingshot trajectory build (microglia_trajectory):
      primary_dims       = mp$primary_dims,
      lineage_substates  = mp$lineage_substates,
      root_substate      = mp$root_substate,
      terminal_substate  = mp$terminal_substate,
      concordance_rho    = mp$concordance_rho,
      concordance_floor  = mp$concordance_floor,
      concordant         = mp$concordant,
      dam_pt_rho         = mp$dam_pt_rho,
      homeo_pt_rho       = mp$homeo_pt_rho,
      omitted_frac_overall = mp$omitted_frac_overall,
      # progression inference (trajectory_progression):
      dam_onset          = tp$provenance$dam_onset,
      used_states        = tp$provenance$used_states,
      within_skip        = tp$provenance$within_skip,
      primary_measures   = tp$provenance$primary_measures,
      planned_primary    = tp$provenance$planned_primary,
      primary_within_skipped = tp$provenance$primary_within_skipped,
      progression_loading = tp$provenance$progression_loading,
      composition_loading = tp$provenance$composition_loading,
      cross_loading      = tp$provenance$cross_loading,
      recon_resid_max    = tp$provenance$recon_resid_max,
      v1_progression_loading = tp$provenance$v1_progression_loading,
      v1_progression_fdr = tp$provenance$v1_progression_fdr,
      n_perm             = tp$provenance$n_perm,
      seed               = tp$provenance$seed,
      limma_version      = tp$provenance$limma_version,
      r_version          = tp$provenance$r_version)) -> out

  # render-cleanliness postconditions (mirror microglia_report_data): the qmd must never hit a ggplot
  # missing-value warning (warn=2 -> render error) nor a sprintf/%d type error. every value the qmd
  # inline-formats or feeds a geom is present + finite; existence is asserted BEFORE is.finite() so a
  # dropped column fails HERE, never vacuously (all(is.finite(NULL)) == TRUE lets a missing col slip).
  on   <- cell_frame$on_lineage
  it   <- out$interaction                                 # qmd irow()/interaction-table source
  gl   <- out$glmm                                        # qmd supportive-arm sentence + provenance line
  pv   <- out$provenance                                  # qmd inline-formats many scalars + feeds 2 geoms
  lpu  <- out$lineage_per_unit                            # qmd conditioning audit: aggregate n_cells/n_on_lineage by genotype
  fin1 <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)   # finite length-1 numeric
  int1 <- function(x) fin1(x) && x == round(x)            # + integer-VALUED (R %d accepts whole doubles -> no coercion)
  str1 <- function(x) is.character(x) && length(x) == 1L && !is.na(x)    # non-NA length-1 string
  # every measure _trajectory.qmd indexes via irow() (match into interaction$measure) + every canonical
  # contrast mp_ctr() pulls from weighted_top -> assert both present + UNIQUE so a dropped/renamed/
  # duplicated measure or contrast fails at the extractor, not as an empty/first-of-dup inline value.
  irow_measures <- c("mean_pt", "comp_cf", "progression_cf", "cross", "within_homeostatic",
                     "within_dam", "median_pt", "q90", "frac_past_logit", "frac_past_asin")
  canonical_contrasts <- c("tau_alone", "nlgf_in_maptki", "nlgf_in_p301s", "tau_in_nlgf", "interaction")
  perm_inlined <- c("mean_pt", "progression_cf")          # rows whose perm_p the qmd prints via round() (headline stat)
  stopifnot(
    nrow(cell_frame) >= 1L, is.logical(on), !anyNA(on), any(on),  # logical mask (cf[cf$on_lineage,]) -> no NA, >= 1 on-lineage
    !anyNA(cell_frame$genotype), !anyNA(cell_frame$substate),
    all(is.finite(cell_frame$pt_raw[on])),                # every on-lineage cell ordered
    all(is.finite(cell_frame$score_axis_pt)),             # score-axis defined for all cells
    all(c("measure", "family", "coef", "ci_l", "ci_r", "p_value", "fdr") %in% names(it)),  # cols exist BEFORE finite
    all(irow_measures %in% it$measure),                   # every qmd-indexed interaction measure present ...
    all(vapply(irow_measures, function(m) sum(it$measure == m) == 1L, logical(1))),  # ... exactly once (match 1st-of-dup)
    # finite coef/CI/p/fdr across BOTH families is an INTENTIONAL build-fatal data-quality gate: a
    # zero-variance endpoint (constant per-unit composition -> comp_cf/cross se=0 -> NaN p/fdr) halts
    # the report rather than tabling NaN. Validated non-degenerate on the current data; if real data
    # ever hits a degenerate EXPLORATORY endpoint, S4b adds graceful "-" formatting + exempts it here.
    all(is.finite(it$coef)), all(is.finite(it$fdr)),
    all(is.finite(it$ci_l)), all(is.finite(it$ci_r)), all(is.finite(it$p_value)),
    "perm_p" %in% names(it),                              # table col it$perm_p (NA-tolerant there); inlined rows:
    all(is.finite(it$perm_p[it$measure %in% perm_inlined])),  # mean_pt+prog_cf prose p
    all(canonical_contrasts %in% names(out$weighted_top)),  # 5 canonical contrasts for mp_ctr()
    all(vapply(canonical_contrasts, function(cn) {         # each mp_ctr(cn) mean_pt row feeds p_ctr geom_pointrange
      w <- out$weighted_top[[cn]]                          # (coef/CI) + the two amyloid rows' p_value the prose (fmt_p)
      is.data.frame(w) && all(c("measure", "coef", "ci_l", "ci_r", "p_value") %in% names(w)) &&
        sum(w$measure == "mean_pt") == 1L &&
        all(is.finite(unlist(w[w$measure == "mean_pt", c("coef", "ci_l", "ci_r", "p_value")])))
    }, logical(1))),
    is.data.frame(out$per_unit), nrow(out$per_unit) >= 1L,     # qmd prints nrow(trd$per_unit) (%d)
    nrow(out$sensitivity) >= 1L, all(is.finite(out$sensitivity$spearman_vs_primary)),  # qmd min(...) -> Inf+warn if empty
    # lineage-conditioning audit (lpu): qmd aggregates n_cells/n_on_lineage by genotype + prints
    # omitted% -> finite positive counts, on-lineage bounded [0, n_cells], known non-missing genotype
    # (an NA/Inf count or unknown genotype -> NaN omitted% or a dropped-level warn under warn=2).
    nrow(lpu) >= 1L, all(is.finite(lpu$n_cells)), all(lpu$n_cells > 0),
    all(is.finite(lpu$n_on_lineage)), all(lpu$n_on_lineage >= 0 & lpu$n_on_lineage <= lpu$n_cells),
    !anyNA(lpu$genotype), all(as.character(lpu$genotype) %in% genotype_levels),
    # glmm supportive arm: method picks the sentence branch; the NOT-failed branch inline-formats
    # effect/CI/p + re_sd (finite -> never prints "NA"/"Inf"); n_cells/n_units print unconditionally.
    length(gl$method) == 1L, gl$method %in% c("glmmTMB_beta", "lmm_ranknorm", "failed"),
    fin1(gl$n_cells), int1(gl$n_units), length(gl$re_sd) == 1L, length(gl$singular) == 1L,
    gl$method == "failed" || (fin1(gl$estimate) && fin1(gl$ci_l) && fin1(gl$ci_r) &&
                              fin1(gl$p_value) && fin1(gl$re_sd)),
    gl$method != "failed" || length(gl$fail_reason) == 1L,   # failed branch prints fail_reason instead
    # provenance: fields the qmd feeds a geom (dam_onset -> geom_hline; *_loading -> geom_col) or inline-
    # formats must be finite (NA = warn=2 red / "NA" text); dims/perm/seed integer-valued for %d; the
    # substate/version strings non-NA; concordant a scalar flag; the rest (assembled, not inline-read) non-NULL.
    all(vapply(pv[c("dam_onset", "composition_loading", "progression_loading", "cross_loading",
                    "dam_pt_rho", "homeo_pt_rho", "concordance_rho", "concordance_floor",
                    "recon_resid_max", "omitted_frac_overall",
                    "v1_progression_loading", "v1_progression_fdr")], fin1, logical(1))),  # last two inline-printed
    # the qmd prints the reconstruction residual + claims the three loadings "sum to one by
    # construction" -> back both at the report layer (run_trajectory_progression gates recon < 1e-8;
    # this catches a corrupted provenance copy where residual + loadings diverge).
    pv$recon_resid_max < 1e-6,
    abs(pv$composition_loading + pv$progression_loading + pv$cross_loading - 1) < 1e-5,
    all(vapply(pv[c("primary_dims", "n_perm", "seed")], int1, logical(1))),
    all(vapply(pv[c("root_substate", "terminal_substate", "limma_version", "r_version")], str1, logical(1))),
    is.logical(pv$concordant), length(pv$concordant) == 1L, !is.na(pv$concordant),
    !any(vapply(pv, is.null, logical(1))))                # every assembled provenance field non-NULL
  out
}

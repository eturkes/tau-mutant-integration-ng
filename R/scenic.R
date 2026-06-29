# --------------------------------------------------------------------
# SCENIC regulon contrast modelling (plan arc K, step K4). Connects the
# data-driven consensus regulons (storage/cache/scenic/, built by the
# pyscenic GRNBoost2 -> cisTarget -> AUCell pipeline in K2/K3) to the
# project's locked 5-contrast framework via three complementary readouts
# plus a target-overlap diagnostic. Pure functions only; the I/O driver
# is scripts/build_scenic_contrasts.R and the display chapter is rmd/20.
#
# Design principle (controlled network swap): readout 1 reuses the §14
# decoupleR machinery from R/tf_inference.R VERBATIM (same NEBULA z, same
# ulm/wsum/consensus statistics, same minsize=5, same FDR<0.10), changing
# ONLY the prior network -- data-driven SCENIC regulons in place of the
# CollecTRI prior. Any difference in TF activity is then attributable to
# the network, not the scoring. Readout 2 scores the native AUCell
# activity through the LOCKED 2x2 factorial (R/design.R + R/de_pb.R).
# Readout 3 mirrors §18's per-substate analysis on the per-state NEBULA z.
#
# Caveat carried throughout (set in K3): the consensus regulons are ALL
# activating (sign "+"); the SCENIC decoupleR net therefore has no
# repressing edges, unlike signed CollecTRI. A §14 repressing call with
# no SCENIC counterpart is a structural consequence of that, not a
# discordance. Stated wherever signs are compared.
# --------------------------------------------------------------------

# Build a decoupleR-compatible prior network from the long regulon table
# (storage/cache/scenic/scenic_regulons.tsv: TF, sign, target, recurrence,
# mean_importance). Mirrors the CollecTRI tibble shape (source/target/mor)
# so it is a drop-in for run_decoupler_per_modality().
#
#   weighted  FALSE (default, locked K4 choice) -> mor = sign as +/-1, the
#             minimal swap that keeps decoupleR's mode-of-regulation
#             semantics identical to CollecTRI (which also carries +/-1
#             mor). TRUE -> mor = sign * mean_importance, a GRNBoost2-
#             weighted variant retained for sensitivity only.
#
# Returns a tibble(source, target, mor).
build_scenic_network <- function(regulons, weighted = FALSE) {
  stopifnot(is.data.frame(regulons),
            all(c("TF", "sign", "target") %in% names(regulons)))
  s <- ifelse(regulons$sign == "+", 1, -1)
  mor <- if (isTRUE(weighted)) {
    stopifnot("mean_importance" %in% names(regulons))
    s * as.numeric(regulons$mean_importance)
  } else {
    s
  }
  tibble::tibble(source = as.character(regulons$TF),
                 target = as.character(regulons$target),
                 mor    = mor)
}

# Thin wrapper: run the §14 decoupleR pipeline on one cache's top-table
# list using the SCENIC network. Returns the per-contrast split list
# (contrast -> tibble(statistic, source, score, p_value)), identical in
# shape to a single modality of the tf_activity_decoupler cache.
run_scenic_decoupler <- function(top_list, scenic_net,
                                 statistics = c("ulm", "wsum"),
                                 minsize = 5L, consensus = TRUE) {
  stat_mat <- extract_de_stat_matrix(top_list)          # stat_col="t", id_col="symbol"
  dec <- run_decoupler_per_modality(stat_mat, network = scenic_net,
                                    statistics = statistics,
                                    minsize = as.integer(minsize),
                                    consensus = consensus)
  split_decoupler_by_contrast(dec)
}

# Tidy a per-contrast decoupler split list into a long activity table with
# BH-adjusted significance, computed WITHIN each contrast across that
# network's full scored-TF set (mirrors §14's rank_tfs_cross_modality FDR
# scope so each network carries its native multiple-testing burden). The
# significance statistic is ulm (the §14 convention); the reported
# magnitude/sign is the cross-method consensus.
#
# Returns tibble(network, source, contrast, ulm_score, ulm_p, padj,
#                cons_score, sig).
decoupler_activity_long <- function(by_contrast, network_label,
                                    padj_cut = 0.10) {
  if (length(by_contrast) == 0L) {
    return(tibble::tibble(network = character(), source = character(),
                          contrast = character(), ulm_score = numeric(),
                          ulm_p = numeric(), padj = numeric(),
                          cons_score = numeric(), sig = logical()))
  }
  rows <- lapply(names(by_contrast), function(cn) {
    tb  <- by_contrast[[cn]]
    ulm <- tb[tb$statistic == "ulm", c("source", "score", "p_value"), drop = FALSE]
    names(ulm) <- c("source", "ulm_score", "ulm_p")
    ulm$padj <- p.adjust(ulm$ulm_p, method = "BH")
    con <- tb[tb$statistic == "consensus", c("source", "score"), drop = FALSE]
    names(con) <- c("source", "cons_score")
    out <- merge(ulm, con, by = "source", all.x = TRUE)
    out$contrast <- cn
    out
  })
  res <- dplyr::bind_rows(rows)
  res$network <- network_label
  res$sig <- !is.na(res$padj) & res$padj < padj_cut
  tibble::as_tibble(res[, c("network", "source", "contrast", "ulm_score",
                            "ulm_p", "padj", "cons_score", "sig")])
}

# Controlled head-to-head: SCENIC vs CollecTRI activity for every TF in
# the union of the SCENIC regulon set and a supplied comparison-TF set
# (the §14 verdict + §18 NF-kB TFs), per contrast. For a TF absent from a
# network the network columns are NA and *_sig is FALSE. sign_concordant
# compares the consensus-score signs (NA if either side absent / zero);
# both_sig is the strict agreement flag.
#
# Returns tibble keyed (TF, contrast) with paired scenic_* / collectri_*
# columns + in_scenic_regulons, scenic_n_targets, is_comparison_tf, axis.
build_scenic_headtohead <- function(scenic_long, collectri_long,
                                    scenic_net, comparison = NULL,
                                    padj_cut = 0.10) {
  scenic_tfs <- sort(unique(scenic_net$source))
  comp_tfs   <- if (is.null(comparison)) character() else unique(comparison$TF)
  universe   <- sort(unique(c(scenic_tfs, comp_tfs)))
  contrasts  <- sort(unique(c(scenic_long$contrast, collectri_long$contrast)))
  grid <- expand.grid(TF = universe, contrast = contrasts,
                      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

  pick <- function(long, pref) {
    sub <- long[, c("source", "contrast", "cons_score", "ulm_score",
                    "ulm_p", "padj", "sig")]
    names(sub) <- c("TF", "contrast",
                    paste0(pref, c("_cons", "_ulm_score", "_ulm_p",
                                   "_padj", "_sig")))
    sub
  }
  out <- merge(grid, pick(scenic_long, "scenic"),
               by = c("TF", "contrast"), all.x = TRUE)
  out <- merge(out, pick(collectri_long, "collectri"),
               by = c("TF", "contrast"), all.x = TRUE)

  out$scenic_sig[is.na(out$scenic_sig)]       <- FALSE
  out$collectri_sig[is.na(out$collectri_sig)] <- FALSE

  ntar <- tapply(scenic_net$target, scenic_net$source, length)
  out$in_scenic_regulons <- out$TF %in% scenic_tfs
  out$scenic_n_targets   <- as.integer(ntar[out$TF])

  ss <- sign(out$scenic_cons); cs <- sign(out$collectri_cons)
  out$sign_concordant <- ifelse(is.na(ss) | is.na(cs) | ss == 0 | cs == 0,
                                NA, ss == cs)
  out$both_sig <- out$scenic_sig & out$collectri_sig

  out$is_comparison_tf <- out$TF %in% comp_tfs
  if (!is.null(comparison)) {
    ax <- tapply(comparison$axis, comparison$TF,
                 function(z) paste(sort(unique(z)), collapse = ";"))
    out$axis <- as.character(ax[out$TF])
  } else {
    out$axis <- NA_character_
  }
  out <- out[order(out$TF, out$contrast), ]
  tibble::as_tibble(out)
}

# Pseudobulk the AUCell per-cell matrix to the 16 genotype_batch ids by
# arithmetic mean per regulon (the locked replicate unit). aucell: a data
# frame with a cell-id first column + one column per "TF(+)" regulon.
# colattrs: CellID + genotype + batch + genotype_batch.
#
# Returns list(mat = regulon x 16 numeric matrix, meta = per-id data frame
# rows = the 16 ids with genotype + batch, rownames = id). Columns of mat
# are ordered to match rownames(meta).
aucell_to_pseudobulk <- function(aucell, colattrs, id_col = "genotype_batch") {
  cell_col <- names(aucell)[1]
  reg_cols <- setdiff(names(aucell), cell_col)
  key <- colattrs[match(aucell[[cell_col]], colattrs$CellID), ]
  stopifnot(!anyNA(key$CellID))
  ids <- key[[id_col]]
  id_levels <- sort(unique(ids))
  mat <- vapply(reg_cols, function(rc) {
    tapply(aucell[[rc]], ids, mean)[id_levels]
  }, numeric(length(id_levels)))
  mat <- t(mat)                                   # regulon x id
  colnames(mat) <- id_levels
  meta <- unique(data.frame(id = ids,
                            genotype = key$genotype,
                            batch = key$batch,
                            stringsAsFactors = FALSE))
  meta <- meta[match(id_levels, meta$id), ]
  rownames(meta) <- meta$id
  list(mat = mat[, rownames(meta), drop = FALSE], meta = meta)
}

# Fit the LOCKED 2x2 factorial + batch design to the pseudobulk AUCell
# matrix and extract all 5 contrasts per regulon via limma. transform:
# "logit" (default) variance-stabilises the bounded [0,1] AUCell means;
# K4 verified every pseudobulk mean lies strictly interior (~0.005-0.097)
# so qlogis needs no epsilon clamp. "none" keeps raw means (sensitivity).
#
# Returns tibble(regulon, contrast, logFC, t, P.Value, adj.P.Val, sig,
#                transform), adj.P.Val being limma's BH across regulons
# within each contrast.
fit_aucell_contrasts <- function(pb_mat, meta, transform = c("logit", "none"),
                                 padj_cut = 0.10) {
  transform <- match.arg(transform)
  mat <- switch(transform,
                logit = {
                  stopifnot(all(pb_mat > 0 & pb_mat < 1))
                  qlogis(pb_mat)
                },
                none  = pb_mat)
  fd <- factorial_design(meta, genotype_col = "genotype",
                         batch_col = "batch", add_batch = TRUE)
  fit <- fit_limma_log(mat[, rownames(fd$design), drop = FALSE],
                       group = meta$genotype,
                       design = fd$design, contrasts = fd$contrasts)
  rows <- lapply(names(fit$top), function(cn) {
    tt <- fit$top[[cn]]
    data.frame(regulon = tt$feature, contrast = cn,
               logFC = tt$logFC, t = tt$t, P.Value = tt$P.Value,
               adj.P.Val = tt$adj.P.Val, stringsAsFactors = FALSE)
  })
  res <- dplyr::bind_rows(rows)
  res$sig <- res$adj.P.Val < padj_cut
  res$transform <- transform
  tibble::as_tibble(res)
}

# Per-substate SCENIC activity at one contrast (default "interaction"),
# mirroring §18. per_state_cache: the de_snrnaseq_nebula_per_state_1pct
# list (substate -> $top -> contrast). Runs the same decoupleR swap per
# substate, BH within substate-contrast, returns the chosen contrast.
#
# Returns tibble(substate, source, contrast, ulm_score, ulm_p, padj,
#                cons_score, sig).
scenic_substate_activity <- function(per_state_cache, scenic_net,
                                     contrast = "interaction",
                                     statistics = c("ulm", "wsum"),
                                     minsize = 5L, padj_cut = 0.10) {
  rows <- lapply(names(per_state_cache), function(st) {
    by_c <- run_scenic_decoupler(per_state_cache[[st]]$top, scenic_net,
                                 statistics = statistics, minsize = minsize)
    long <- decoupler_activity_long(by_c, network_label = "scenic",
                                    padj_cut = padj_cut)
    sub <- long[long$contrast == contrast, , drop = FALSE]
    sub$substate <- st
    sub
  })
  res <- dplyr::bind_rows(rows)
  res <- res[, c("substate", "source", "contrast", "ulm_score", "ulm_p",
                 "padj", "cons_score", "sig")]
  tibble::as_tibble(res)
}

# Jaccard overlap of target sets for every TF shared between the SCENIC
# regulon network and the CollecTRI prior -- how much of each data-driven
# regulon the literature prior already encoded. Both nets are source/
# target tibbles. Returns tibble(TF, n_scenic_targets, n_collectri_targets,
# n_shared, n_union, jaccard) for the shared TFs, descending jaccard.
scenic_collectri_target_overlap <- function(scenic_net, collectri_net) {
  s_by <- split(scenic_net$target, scenic_net$source)
  c_by <- split(collectri_net$target, collectri_net$source)
  shared <- intersect(names(s_by), names(c_by))
  rows <- lapply(shared, function(tf) {
    a <- unique(s_by[[tf]]); b <- unique(c_by[[tf]])
    inter <- length(intersect(a, b)); uni <- length(union(a, b))
    data.frame(TF = tf, n_scenic_targets = length(a),
               n_collectri_targets = length(b), n_shared = inter,
               n_union = uni, jaccard = if (uni > 0) inter / uni else NA_real_,
               stringsAsFactors = FALSE)
  })
  res <- dplyr::bind_rows(rows)
  res <- res[order(-res$jaccard), ]
  tibble::as_tibble(res)
}

# Recovery ladder for the full regulon set + the §14/§18 comparison TFs.
# Classifies each TF by why it did or did not enter the locked >=8/10
# consensus, distinguishing structural non-recovery (below the >=1%
# expression floor -> never a candidate regulator) from motif pruning
# (a candidate with adjacencies but no motif-supported regulon) from
# stochastic near-misses. census: scenic_recovery_census.tsv (the
# comparison-TF audit from K3). scenic_net: the consensus network.
# candidate_tfs: TFs eligible as GRN regulators (filtered genes ∩ allTFs).
#
# Returns tibble(TF, axis, is_comparison_tf, in_consensus,
#                n_consensus_targets, max_runs_present, is_candidate,
#                recovery_class).
build_recovery_table <- function(census, scenic_net, candidate_tfs) {
  ntar <- tapply(scenic_net$target, scenic_net$source, length)
  consensus_tfs <- names(ntar)

  cen <- as.data.frame(census, stringsAsFactors = FALSE)
  cen_axis <- tapply(cen$axis, cen$TF,
                     function(z) paste(sort(unique(z)), collapse = ";"))
  cen_max  <- tapply(cen$max_runs_present, cen$TF, max)
  cen_pres <- tapply(cen$regulon_present_ge_thr, cen$TF, max)
  cen_inco <- tapply(cen$in_consensus, cen$TF, max)

  universe <- sort(unique(c(consensus_tfs, cen$TF)))
  classify <- function(tf) {
    in_con <- tf %in% consensus_tfs ||
      (!is.na(cen_inco[tf]) && cen_inco[tf] == 1)
    is_comp <- tf %in% cen$TF
    if (in_con) return(if (is_comp) "recovered_comparison" else "recovered_novel")
    pres <- if (!is.na(cen_pres[tf])) cen_pres[tf] else 0
    mx   <- if (!is.na(cen_max[tf]))  cen_max[tf]  else 0
    if (pres == 1) return("present_targets_subthreshold")   # regulon formed, <5 stable targets
    if (mx >= 6)   return("near_miss")                       # 6-7/10 runs
    if (mx >= 1)   return("weak")                            # 1-5/10 runs
    if (tf %in% candidate_tfs) return("motif_pruned")        # candidate, 0 regulons
    "expression_floor"                                       # below the >=1% filter
  }
  rows <- lapply(universe, function(tf) {
    data.frame(
      TF = tf,
      axis = as.character(if (!is.na(cen_axis[tf])) cen_axis[tf] else NA_character_),
      is_comparison_tf = tf %in% cen$TF,
      in_consensus = tf %in% consensus_tfs,
      n_consensus_targets = if (tf %in% consensus_tfs) as.integer(ntar[tf]) else 0L,
      max_runs_present = if (!is.na(cen_max[tf])) as.integer(cen_max[tf])
                         else if (tf %in% consensus_tfs) NA_integer_ else 0L,
      is_candidate = tf %in% candidate_tfs,
      recovery_class = classify(tf),
      stringsAsFactors = FALSE)
  })
  res <- dplyr::bind_rows(rows)
  res <- res[order(!res$in_consensus, res$is_comparison_tf == FALSE, res$TF), ]
  tibble::as_tibble(res)
}

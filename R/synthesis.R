# P5 synthesis helpers. This layer is read-only: it compresses the existing compact
# report targets into claim/status rows and deliberately avoids new inference machinery.

synthesis_status_levels <- function() {
  c("core_supported", "corroborated", "focused_support",
    "not_supported", "not_earned", "open_caveat")
}

synthesis_contrasts <- function() {
  pairwise <- names(contrast_definitions)
  c("tau_alone", setdiff(pairwise, "tau_alone"), "interaction")
}

synthesis_focus_contrasts <- function() {
  c("interaction", "nlgf_in_maptki", "nlgf_in_p301s", "tau_in_nlgf")
}

.synthesis_require_cols <- function(x, cols, label) {
  if (!is.data.frame(x)) stop(label, " must be a data frame", call. = FALSE)
  missing <- setdiff(cols, names(x))
  if (length(missing)) {
    stop(label, " missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  TRUE
}

.synthesis_finite_cols <- function(x, cols, label) {
  .synthesis_require_cols(x, cols, label)
  for (col in cols) {
    if (!is.numeric(x[[col]]) && !is.integer(x[[col]])) {
      stop(label, " column is not numeric: ", col, call. = FALSE)
    }
    if (!all(is.finite(x[[col]]))) {
      stop(label, " column has non-finite values: ", col, call. = FALSE)
    }
  }
  TRUE
}

.synthesis_one_row <- function(x, idx, label) {
  if (!is.logical(idx) || length(idx) != nrow(x)) {
    stop("bad row index for ", label, call. = FALSE)
  }
  hit <- which(idx %in% TRUE)
  if (length(hit) != 1L) {
    stop("expected exactly one ", label, "; found ", length(hit), call. = FALSE)
  }
  x[hit, , drop = FALSE]
}

.synthesis_fmt_p <- function(p) {
  if (!is.finite(p)) return("NA")
  formatC(p, format = "e", digits = 2)
}

.synthesis_fmt_num <- function(x, digits = 2L) {
  if (!is.finite(x)) return("NA")
  formatC(x, format = "f", digits = digits)
}

synthesis_validate_evidence_table <- function(evidence_table) {
  required <- c("claim_id", "axis", "status", "direction", "evidence",
                "primary_sources", "supporting_sources", "caveat", "report_anchor")
  .synthesis_require_cols(evidence_table, required, "synthesis evidence_table")
  forbidden <- c("support", "contradict", "net_score", "evidence_score", "model_score",
                 "model_margin", "contest_margin", "ledger_id", "arc_id", "claim_score")
  bad <- intersect(forbidden, names(evidence_table))
  if (length(bad)) {
    stop("synthesis evidence_table has ledger-like columns: ",
         paste(bad, collapse = ", "), call. = FALSE)
  }
  if (!nrow(evidence_table)) stop("synthesis evidence_table is empty", call. = FALSE)
  if (anyNA(evidence_table$claim_id) || any(evidence_table$claim_id == "")) {
    stop("synthesis evidence_table has blank claim_id", call. = FALSE)
  }
  if (anyDuplicated(evidence_table$claim_id)) {
    stop("synthesis evidence_table has duplicated claim_id", call. = FALSE)
  }
  bad_status <- setdiff(unique(evidence_table$status), synthesis_status_levels())
  if (length(bad_status)) {
    stop("synthesis evidence_table has invalid status: ",
         paste(bad_status, collapse = ", "), call. = FALSE)
  }
  char_cols <- setdiff(required, "status")
  for (col in char_cols) {
    if (!is.character(evidence_table[[col]]) ||
        anyNA(evidence_table[[col]]) ||
        any(evidence_table[[col]] == "")) {
      stop("synthesis evidence_table has blank/non-character column: ", col, call. = FALSE)
    }
  }
  TRUE
}

.synthesis_evidence_row <- function(claim_id, axis, status, direction, evidence,
                                    primary_sources, supporting_sources, caveat,
                                    report_anchor) {
  data.frame(
    claim_id = claim_id,
    axis = axis,
    status = status,
    direction = direction,
    evidence = evidence,
    primary_sources = primary_sources,
    supporting_sources = supporting_sources,
    caveat = caveat,
    report_anchor = report_anchor,
    stringsAsFactors = FALSE
  )
}

.synthesis_status_summary <- function(evidence_table) {
  tab <- table(factor(evidence_table$status, levels = synthesis_status_levels()))
  data.frame(status = names(tab), n = as.integer(tab), row.names = NULL,
             stringsAsFactors = FALSE)
}

synthesis_report_data <- function(microglia_report, trajectory_report, mechanism_report,
                                  crossmodality_report, alpha = 0.10) {
  stopifnot(is.list(microglia_report), is.list(trajectory_report),
            is.list(mechanism_report), is.list(crossmodality_report),
            is.numeric(alpha), length(alpha) == 1L, is.finite(alpha),
            alpha > 0, alpha < 1)
  contrasts <- synthesis_contrasts()
  focus <- synthesis_focus_contrasts()

  # P1 microglia compact anchor.
  stopifnot("n_cells" %in% names(microglia_report),
            is.list(microglia_report$provenance),
            "substate_table" %in% names(microglia_report$provenance),
            is.numeric(microglia_report$n_cells),
            length(microglia_report$n_cells) == 1L,
            is.finite(microglia_report$n_cells),
            microglia_report$n_cells > 0)
  substate_table <- as.matrix(microglia_report$provenance$substate_table)
  stopifnot(all(genotype_levels %in% rownames(substate_table)),
            all(c("Homeostatic", "DAM") %in% colnames(substate_table)))
  storage.mode(substate_table) <- "double"
  stopifnot(all(is.finite(substate_table)), all(substate_table >= 0))
  substate_table <- substate_table[genotype_levels, , drop = FALSE]
  microglia_summary <- data.frame(
    genotype = rownames(substate_table),
    n_cells = rowSums(substate_table),
    homeostatic = substate_table[, "Homeostatic"],
    dam = substate_table[, "DAM"],
    dam_fraction = substate_table[, "DAM"] / pmax(rowSums(substate_table), 1),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  stopifnot(all(is.finite(microglia_summary$dam_fraction)))

  # P1/P2/P3 anchors carried by the mechanism and trajectory report bundles.
  comp <- mechanism_report$composition_anchor
  .synthesis_require_cols(comp, c("method", "contrast", "substate", "prop_ratio", "t",
                                  "p_value", "fdr_global"),
                          "mechanism_report$composition_anchor")
  .synthesis_finite_cols(comp, c("prop_ratio", "t", "p_value", "fdr_global"),
                         "mechanism_report$composition_anchor")
  dam_maptki <- .synthesis_one_row(comp, comp$contrast == "nlgf_in_maptki" &
                                     comp$substate == "DAM", "DAM nlgf_in_maptki composition")
  dam_p301s <- .synthesis_one_row(comp, comp$contrast == "nlgf_in_p301s" &
                                    comp$substate == "DAM", "DAM nlgf_in_p301s composition")
  homeo_maptki <- .synthesis_one_row(comp, comp$contrast == "nlgf_in_maptki" &
                                       comp$substate == "Homeostatic",
                                     "Homeostatic nlgf_in_maptki composition")
  homeo_p301s <- .synthesis_one_row(comp, comp$contrast == "nlgf_in_p301s" &
                                      comp$substate == "Homeostatic",
                                    "Homeostatic nlgf_in_p301s composition")
  dam_interaction <- .synthesis_one_row(comp, comp$contrast == "interaction" &
                                          comp$substate == "DAM", "DAM composition interaction")

  amyloid_dam_supported <- all(c(dam_maptki$t, dam_p301s$t) > 0) &&
    all(c(dam_maptki$fdr_global, dam_p301s$fdr_global) < alpha) &&
    all(c(homeo_maptki$t, homeo_p301s$t) < 0) &&
    all(c(homeo_maptki$fdr_global, homeo_p301s$fdr_global) < alpha)
  dam_interaction_supported <- dam_interaction$t > 0 && dam_interaction$fdr_global < alpha
  if (!amyloid_dam_supported) {
    stop("amyloid-to-DAM support anchor failed current synthesis contract", call. = FALSE)
  }
  if (!dam_interaction_supported) {
    stop("DAM composition interaction anchor failed current synthesis contract", call. = FALSE)
  }

  traj <- trajectory_report$interaction
  .synthesis_require_cols(traj, c("family", "measure", "coef", "p_value", "perm_p",
                                  "fdr", "ci_l", "ci_r"),
                          "trajectory_report$interaction")
  .synthesis_finite_cols(traj, c("coef", "p_value", "fdr", "ci_l", "ci_r"),
                         "trajectory_report$interaction")
  mean_pt <- .synthesis_one_row(traj, traj$measure == "mean_pt", "trajectory mean_pt")
  comp_cf <- .synthesis_one_row(traj, traj$measure == "comp_cf", "trajectory comp_cf")
  progression_cf <- .synthesis_one_row(traj, traj$measure == "progression_cf",
                                       "trajectory progression_cf")
  within_homeo <- .synthesis_one_row(traj, traj$measure == "within_homeostatic",
                                     "trajectory within_homeostatic")
  if (!(comp_cf$fdr < alpha)) {
    stop("trajectory comp_cf anchor is not significant under the synthesis contract",
         call. = FALSE)
  }
  if (progression_cf$fdr < alpha || within_homeo$fdr < alpha) {
    stop("trajectory progression-beyond-composition anchor changed support status",
         call. = FALSE)
  }
  stopifnot(is.list(trajectory_report$provenance),
            all(c("composition_loading", "progression_loading", "cross_loading") %in%
                  names(trajectory_report$provenance)))
  trajectory_provenance <- trajectory_report$provenance[c("composition_loading",
                                                          "progression_loading",
                                                          "cross_loading")]
  stopifnot(all(vapply(trajectory_provenance, function(x)
    is.numeric(x) && length(x) == 1L && is.finite(x), logical(1))))

  tf <- mechanism_report$tf_highlights
  .synthesis_require_cols(tf, c("population", "source", "contrast", "score",
                                "p_value", "fdr"),
                          "mechanism_report$tf_highlights")
  .synthesis_finite_cols(tf, c("score", "p_value", "fdr"),
                         "mechanism_report$tf_highlights")
  myc <- .synthesis_one_row(tf, tf$population == "whole_microglia" &
                              tf$source == "Myc" &
                              tf$contrast == "interaction",
                            "whole Myc interaction TF")
  if (!(myc$score < 0 && myc$fdr < alpha)) {
    stop("whole Myc interaction TF anchor is not supported under the synthesis contract",
         call. = FALSE)
  }

  stopifnot(is.list(mechanism_report$nfkb),
            is.data.frame(mechanism_report$nfkb$table),
            is.list(mechanism_report$nfkb$verdict))
  nf <- mechanism_report$nfkb$verdict
  stopifnot(all(c("status", "supported", "primary_discordant", "primary_negative",
                  "primary_positive", "alpha") %in% names(nf)),
            nf$status %in% c("supported", "not_supported", "discordant"),
            is.logical(nf$supported), length(nf$supported) == 1L, !is.na(nf$supported))
  if (isTRUE(nf$supported)) {
    stop("NF-kB attenuation anchor changed to supported", call. = FALSE)
  }
  .synthesis_require_cols(mechanism_report$nfkb$table,
                          c("test", "score", "p_value", "primary_family_fdr",
                            "primary_test", "detail"),
                          "mechanism_report$nfkb$table")
  nf_primary <- mechanism_report$nfkb$table[
    mechanism_report$nfkb$table$primary_test %in% TRUE, , drop = FALSE]
  if (nrow(nf_primary) != 2L) {
    stop("mechanism_report$nfkb$table expected two primary rows", call. = FALSE)
  }

  kinase <- mechanism_report$kinase$table
  .synthesis_require_cols(kinase, c("source", "contrast", "score", "p_value", "fdr",
                                    "significant", "run_index_supports",
                                    "run_order_confounded"),
                          "mechanism_report$kinase$table")
  .synthesis_finite_cols(kinase, c("score", "p_value", "fdr"),
                         "mechanism_report$kinase$table")
  gsk <- kinase[kinase$source == "Gsk3b", , drop = FALSE]
  missing_gsk <- setdiff(contrasts, gsk$contrast)
  if (length(missing_gsk)) {
    stop("Gsk3b kinase rows missing contrasts: ", paste(missing_gsk, collapse = ", "),
         call. = FALSE)
  }
  gsk_interaction <- .synthesis_one_row(kinase, kinase$source == "Gsk3b" &
                                          kinase$contrast == "interaction",
                                        "Gsk3b interaction kinase")
  gsk_tau_nlgf <- .synthesis_one_row(kinase, kinase$source == "Gsk3b" &
                                       kinase$contrast == "tau_in_nlgf",
                                     "Gsk3b tau_in_nlgf kinase")
  gsk_focal_supported <- any(gsk_interaction$significant %in% TRUE,
                             gsk_tau_nlgf$significant %in% TRUE,
                             gsk_interaction$fdr < alpha,
                             gsk_tau_nlgf$fdr < alpha)
  if (gsk_focal_supported) {
    stop("Gsk3b focal kinase anchor changed to supported", call. = FALSE)
  }

  # P4 compact cross-modality anchors.
  stopifnot(is.list(crossmodality_report$geomx),
            is.list(crossmodality_report$bulk),
            is.list(crossmodality_report$clearance),
            is.list(crossmodality_report$pathway))
  geomx_sens <- crossmodality_report$geomx$sensitivity
  .synthesis_require_cols(geomx_sens, c("fit", "status", "duplicate_correlation_used",
                                        "reason"),
                          "crossmodality_report$geomx$sensitivity")
  bulk_run <- crossmodality_report$bulk$run_index
  .synthesis_require_cols(bulk_run, c("layer", "contrast", "status", "reason",
                                      "n_primary_sig", "n_lost_or_flipped"),
                          "crossmodality_report$bulk$run_index")
  .synthesis_finite_cols(bulk_run, c("n_primary_sig", "n_lost_or_flipped"),
                         "crossmodality_report$bulk$run_index")
  if (!all(focus %in% bulk_run$contrast)) {
    stop("crossmodality_report$bulk$run_index missing focus contrasts", call. = FALSE)
  }

  spatial_decon <- crossmodality_report$clearance$spatial_decon
  verdict <- crossmodality_report$clearance$verdict
  stopifnot(is.list(spatial_decon), all(c("status", "action", "reasons") %in% names(spatial_decon)),
            spatial_decon$status %in% c("earned", "defer", "blocked"),
            is.character(spatial_decon$action), length(spatial_decon$action) == 1L,
            !is.na(spatial_decon$action),
            is.character(spatial_decon$reasons), !anyNA(spatial_decon$reasons),
            is.list(verdict), all(c("status", "ccc_called") %in% names(verdict)),
            verdict$status %in% c("earned", "not_earned"),
            is.logical(verdict$ccc_called), length(verdict$ccc_called) == 1L,
            !is.na(verdict$ccc_called))
  if (identical(spatial_decon$status, "earned") || isTRUE(verdict$ccc_called)) {
    stop("SpatialDecon/full CCC anchor is now earned; revise synthesis", call. = FALSE)
  }
  pairs <- crossmodality_report$clearance$pair_support
  .synthesis_require_cols(pairs, c("pair", "contrast", "n_sides_measured",
                                   "modalities_measured",
                                   "coherent_supported_modalities",
                                   "n_coherent_supported_modalities",
                                   "microglia_strong", "status"),
                          "crossmodality_report$clearance$pair_support")
  earned_pairs <- pairs[pairs$status == "earned", , drop = FALSE]

  pathway_axis <- crossmodality_report$pathway$axis_summary
  .synthesis_require_cols(pathway_axis, c("collection", "set", "axis", "contrast",
                                          "n_modalities_present", "n_modalities_sig",
                                          "mixed_sign", "consistent_direction",
                                          "rank_score"),
                          "crossmodality_report$pathway$axis_summary")
  .synthesis_finite_cols(pathway_axis, c("n_modalities_present", "n_modalities_sig",
                                         "rank_score"),
                         "crossmodality_report$pathway$axis_summary")
  amyloid_axis <- pathway_axis[pathway_axis$contrast %in% c("nlgf_in_maptki", "nlgf_in_p301s") &
                                 pathway_axis$axis %in% c("DAM", "synaptic", "clearance"),
                               , drop = FALSE]
  if (!nrow(amyloid_axis)) {
    stop("crossmodality pathway axis summary missing amyloid DAM/synaptic/clearance rows",
         call. = FALSE)
  }
  if (!any(amyloid_axis$n_modalities_sig > 0)) {
    stop("crossmodality amyloid-axis rows have no significant modality support",
         call. = FALSE)
  }

  bulk_primary_sig <- sum(bulk_run$n_primary_sig, na.rm = TRUE)
  bulk_lost <- sum(bulk_run$n_lost_or_flipped, na.rm = TRUE)
  earned_text <- if (nrow(earned_pairs)) {
    paste(paste0(earned_pairs$pair, " in ", earned_pairs$contrast), collapse = "; ")
  } else {
    "no clearance pair passed the earned CCC-lite rule"
  }
  spatial_text <- paste(spatial_decon$reasons, collapse = "; ")

  evidence_table <- do.call(rbind, list(
    .synthesis_evidence_row(
      "amyloid_dam_activation", "microglia_state", "core_supported",
      "amyloid increases DAM and lowers homeostatic microglia in both tau backgrounds",
      paste0("DAM propeller rows are positive in both amyloid contrasts (FDR ",
             .synthesis_fmt_p(max(dam_maptki$fdr_global, dam_p301s$fdr_global)), "); ",
             "homeostatic rows move oppositely."),
      "microglia_report; mechanism_report$composition_anchor",
      "crossmodality_report$pathway",
      "Substate labels are cluster-primary and snRNAseq under-detects part of the DAM programme.",
      "#sec-microglia"
    ),
    .synthesis_evidence_row(
      "tau_amyloid_dam_composition", "interaction", "core_supported",
      "the tau-amyloid interaction is expressed mainly as extra DAM-cell composition",
      paste0("DAM composition interaction FDR ", .synthesis_fmt_p(dam_interaction$fdr_global),
             "; trajectory comp_cf effect ", .synthesis_fmt_num(comp_cf$coef),
             " (FDR ", .synthesis_fmt_p(comp_cf$fdr), ")."),
      "mechanism_report$composition_anchor; trajectory_report$interaction",
      "crossmodality_report",
      "This is a composition result, not a per-cell acceleration claim.",
      "#sec-trajectory"
    ),
    .synthesis_evidence_row(
      "progression_beyond_composition", "interaction", "not_supported",
      "no statistically supported further advance beyond composition is detected",
      paste0("progression_cf effect ", .synthesis_fmt_num(progression_cf$coef),
             " with FDR ", .synthesis_fmt_p(progression_cf$fdr),
             "; within-homeostatic FDR ", .synthesis_fmt_p(within_homeo$fdr), "."),
      "trajectory_report$interaction",
      "trajectory_report$provenance",
      "This is absence of supported evidence, not proof of no effect.",
      "#sec-trajectory"
    ),
    .synthesis_evidence_row(
      "myc_rna_interaction", "mechanism", "focused_support",
      "RNA TF activity supports a Myc-linked interaction signal",
      paste0("whole-microglia Myc ULM score ", .synthesis_fmt_num(myc$score),
             " with FDR ", .synthesis_fmt_p(myc$fdr), "."),
      "mechanism_report$tf_highlights",
      "mechanism_report$pathway_project",
      "TF activity is inferred from pseudobulk RNA ranks, not direct protein activity.",
      "#sec-mechanism"
    ),
    .synthesis_evidence_row(
      "nfkb_attenuation", "mechanism",
      "not_supported",
      "NF-kB attenuation does not pass the primary gate",
      paste0("NF-kB gate status is ", nf$status, " across the two primary rows."),
      "mechanism_report$nfkb",
      "mechanism_report$pathway_project",
      "Discordant or non-significant RNA target/TF rows are not averaged into support.",
      "#sec-mechanism"
    ),
    .synthesis_evidence_row(
      "gsk3b_kinase", "mechanism",
      "not_supported",
      "Gsk3b is covered but not recovered as a focal rebuilt kinase signal",
      paste0("Gsk3b interaction FDR ", .synthesis_fmt_p(gsk_interaction$fdr),
             "; tau-in-NLGF FDR ", .synthesis_fmt_p(gsk_tau_nlgf$fdr), "."),
      "mechanism_report$kinase",
      "crossmodality_report$bulk",
      "The kinase layer is 24M bulk hippocampus, not microglia-sorted.",
      "#sec-mechanism"
    ),
    .synthesis_evidence_row(
      "crossmodality_amyloid_axes", "cross_modality", "corroborated",
      "GeoMx and bulk layers mainly corroborate amyloid-response and synaptic-clearance axes",
      paste0("Amyloid DAM/synaptic/clearance axis rows cover ",
             length(unique(amyloid_axis$axis)), " axes; bulk primary significant rows total ",
             bulk_primary_sig, "."),
      "crossmodality_report$pathway; crossmodality_report$geomx",
      "crossmodality_report$bulk; crossmodality_report$divergence",
      "Interaction evidence outside microglia composition remains smaller and mixed by modality.",
      "#sec-crossmodality"
    ),
    .synthesis_evidence_row(
      "clearance_axis", "cross_modality",
      if (nrow(earned_pairs)) "focused_support" else "not_earned",
      if (nrow(earned_pairs)) "measured clearance-axis support is focused" else
        "measured clearance-axis pair support is not earned",
      earned_text,
      "crossmodality_report$clearance$pair_support",
      "microglia RNA; GeoMx spatial; bulk hippocampus",
      "CCC-lite is a measured-anchor rule and is not a full cell-cell communication model.",
      "#sec-crossmodality"
    ),
    .synthesis_evidence_row(
      "spatial_decon_full_ccc", "unearned", "not_earned",
      "SpatialDecon abundance and full CCC are not earned in this rebuild",
      paste0("SpatialDecon status is ", spatial_decon$status,
             "; ccc_called is ", verdict$ccc_called, "."),
      "crossmodality_report$clearance$spatial_decon",
      "crossmodality_report$clearance$verdict",
      spatial_text,
      "#sec-crossmodality"
    ),
    .synthesis_evidence_row(
      "bulk_run_index_sensitivity", "caveat", "open_caveat",
      "bulk proteome/phospho signals are sensitive to genotype-blocked run order",
      paste0(bulk_lost, " of ", bulk_primary_sig,
             " primary-significant bulk rows are lost or flipped under run-index sensitivity."),
      "crossmodality_report$bulk$run_index",
      "crossmodality_report$bulk$run_index",
      "Bulk hippocampus rows are corroborative unless run-index support is explicit.",
      "#sec-crossmodality"
    )
  ))
  rownames(evidence_table) <- NULL
  synthesis_validate_evidence_table(evidence_table)

  status_summary <- .synthesis_status_summary(evidence_table)
  open <- evidence_table[evidence_table$status %in% c("not_supported", "not_earned", "open_caveat"),
                         c("claim_id", "axis", "status", "direction", "caveat", "report_anchor"),
                         drop = FALSE]
  rownames(open) <- NULL

  out <- list(
    headline = c(
      "Amyloid drives a microglial homeostatic-to-DAM activation programme.",
      "Mutant tau modulates the amyloid response mainly through DAM-cell composition, not supported further activation-axis progression.",
      "Mechanism support is asymmetric: Myc is supported in RNA, while NF-kB attenuation and Gsk3b are not recovered.",
      "Cross-modality evidence corroborates the amyloid-response and synaptic-clearance axes; SpatialDecon/full CCC remain unearned."
    ),
    evidence_table = evidence_table,
    status_summary = status_summary,
    open_questions = open,
    source_highlights = list(
      microglia_substates = microglia_summary,
      trajectory = traj[traj$measure %in% c("mean_pt", "comp_cf", "progression_cf",
                                            "within_homeostatic"),
                        , drop = FALSE],
      mechanism = list(
        myc = myc,
        nfkb_verdict = nf,
        gsk3b_focal = rbind(gsk_interaction, gsk_tau_nlgf)
      ),
      crossmodality = list(
        geomx_sensitivity = geomx_sens,
        bulk_run_index = bulk_run,
        earned_pairs = earned_pairs,
        spatial_decon = spatial_decon
      )
    ),
    provenance = list(
      source_targets = c("microglia_report", "trajectory_report",
                         "mechanism_report", "crossmodality_report"),
      alpha = alpha,
      statuses = synthesis_status_levels(),
      contrasts = contrasts,
      focus_contrasts = focus,
      trajectory_loadings = trajectory_provenance,
      r_version = paste(R.version$major, R.version$minor, sep = "."),
      report_contract = paste(
        "compact read-only synthesis over existing P1-P4 report targets;",
        "no v1 ledger, contest scoring, new modality, or heavy report-time read"
      )
    )
  )

  stopifnot(length(out$headline) >= 3L,
            nrow(out$evidence_table) >= 8L,
            nrow(out$evidence_table) <= 14L,
            nrow(out$status_summary) == length(synthesis_status_levels()),
            sum(out$status_summary$n) == nrow(out$evidence_table),
            nrow(out$open_questions) >= 1L,
            !("cell_frame" %in% names(out)),
            identical(out$provenance$source_targets,
                      c("microglia_report", "trajectory_report",
                        "mechanism_report", "crossmodality_report")))
  out
}

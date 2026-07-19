# P7.5: execute the frozen DAM-occupancy robustness grid. Every perturbation is
# one-at-a-time from the reference; E1 beta-binomial probability occupancy defines the verdict.

occupancy_sweep_estimators <- function() {
  c("beta_binomial", "empirical_logit", "simple_proportion")
}

occupancy_sweep_frozen_interaction <- function() {
  c(
    e1_estimate = 0.17411412574101759,
    e1_ci_l = 0.095412497552753661,
    e1_ci_r = 0.25281575392928152,
    e1_fdr_zero = 1.8129871062254304e-05,
    e1_fdr_minimum = 0.081166440287495756,
    e2_coef = 0.93405922366503336,
    e2_ci_l = 0.16587418318000835,
    e2_ci_r = 1.7022442641500584,
    e2_fdr = 0.028067318488813262,
    e2_perm_p = 0.021000000000000001,
    e3_estimate = 0.17333675165569132,
    e3_ci_l = 0.042849082730726323,
    e3_ci_r = 0.30382442058065628,
    e3_fdr = 0.018544617034855316
  )
}

occupancy_sweep_interaction_row <- function(x, label) {
  stopifnot(is.data.frame(x), "contrast" %in% names(x))
  out <- x[x$contrast == "interaction", , drop = FALSE]
  if (nrow(out) != 1L) {
    stop(label, " must contain exactly one interaction row", call. = FALSE)
  }
  out
}

occupancy_sweep_interaction_vector <- function(fit) {
  stopifnot(
    is.list(fit), is.list(fit$diagnostics),
    identical(fit$diagnostics$estimators$beta_binomial$status, "ok"),
    identical(fit$diagnostics$estimators$empirical_logit$status, "ok"),
    identical(fit$diagnostics$estimators$simple_proportion$status, "ok")
  )
  e1 <- occupancy_sweep_interaction_row(
    fit$probability_contrasts, "E1 probability contrasts"
  )
  e2 <- occupancy_sweep_interaction_row(fit$empirical_logit, "E2 empirical logit")
  e3 <- occupancy_sweep_interaction_row(fit$simple_proportion, "E3 simple proportion")
  stopifnot(
    is.data.frame(fit$permutation), nrow(fit$permutation) == 1L,
    identical(fit$permutation$contrast, "interaction")
  )
  c(
    e1_estimate = unname(e1$estimate),
    e1_ci_l = unname(e1$ci_l),
    e1_ci_r = unname(e1$ci_r),
    e1_fdr_zero = unname(e1$fdr_zero),
    e1_fdr_minimum = unname(e1$fdr_minimum),
    e2_coef = unname(e2$coef),
    e2_ci_l = unname(e2$ci_l),
    e2_ci_r = unname(e2$ci_r),
    e2_fdr = unname(e2$fdr),
    e2_perm_p = unname(fit$permutation$perm_p),
    e3_estimate = unname(e3$coef),
    e3_ci_l = unname(e3$ci_l),
    e3_ci_r = unname(e3$ci_r),
    e3_fdr = unname(e3$fdr)
  )
}

occupancy_sweep_max_abs_diff <- function(x, y) {
  stopifnot(identical(names(x), names(y)), all(is.finite(x)), all(is.finite(y)))
  max(abs(unname(x) - unname(y)))
}

occupancy_sweep_spec_row <- function(
    variant_id, axis, setting, is_reference = FALSE,
    primary_resolution = 0.4, id_floor = 0.15, mglike_floor = 0.30,
    tol = 0.10, amb_floor = 0.10,
    dropped_unit = NA_character_, dropped_batch = NA_character_) {
  data.frame(
    variant_id = variant_id,
    axis = axis,
    setting = setting,
    is_reference = is_reference,
    primary_resolution = primary_resolution,
    id_floor = id_floor,
    mglike_floor = mglike_floor,
    tol = tol,
    amb_floor = amb_floor,
    dropped_unit = dropped_unit,
    dropped_batch = dropped_batch,
    stringsAsFactors = FALSE
  )
}

occupancy_sweep_grid <- function(units, batches) {
  units <- sort(unique(as.character(units)), method = "radix")
  batches <- sort(unique(as.character(batches)), method = "radix")
  stopifnot(length(units) == 16L, length(batches) == 4L)

  rows <- c(
    list(occupancy_sweep_spec_row(
      "reference", "reference", "reference labeling", is_reference = TRUE
    )),
    lapply(c(0.2, 0.3, 0.5, 0.6, 0.8), function(value) {
      occupancy_sweep_spec_row(
        paste0("resolution_", sprintf("%.1f", value)), "resolution",
        paste0("primary_resolution=", sprintf("%.1f", value)),
        primary_resolution = value
      )
    }),
    lapply(c(0.10, 0.20), function(value) {
      occupancy_sweep_spec_row(
        paste0("id_floor_", sprintf("%.2f", value)), "id_floor",
        paste0("id_floor=", sprintf("%.2f", value)), id_floor = value
      )
    }),
    lapply(c(0.20, 0.40), function(value) {
      occupancy_sweep_spec_row(
        paste0("mglike_floor_", sprintf("%.2f", value)), "mglike_floor",
        paste0("mglike_floor=", sprintf("%.2f", value)), mglike_floor = value
      )
    }),
    list(occupancy_sweep_spec_row(
      "no_prune", "no_prune", "no pruning", id_floor = 0, mglike_floor = 0
    )),
    lapply(c(0.05, 0.15), function(value) {
      occupancy_sweep_spec_row(
        paste0("tol_", sprintf("%.2f", value)), "tol",
        paste0("tol=", sprintf("%.2f", value)), tol = value
      )
    }),
    lapply(c(0.05, 0.15), function(value) {
      occupancy_sweep_spec_row(
        paste0("amb_floor_", sprintf("%.2f", value)), "amb_floor",
        paste0("amb_floor=", sprintf("%.2f", value)), amb_floor = value
      )
    }),
    lapply(units, function(value) {
      occupancy_sweep_spec_row(
        paste0("LOU_", value), "leave_one_unit_out", paste0("drop ", value),
        dropped_unit = value
      )
    }),
    lapply(batches, function(value) {
      occupancy_sweep_spec_row(
        paste0("LOBO_", value), "leave_one_batch_out", paste0("drop ", value),
        dropped_batch = value
      )
    })
  )
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  expected <- c(
    reference = 1L, resolution = 5L, id_floor = 2L, mglike_floor = 2L,
    no_prune = 1L, tol = 2L, amb_floor = 2L,
    leave_one_unit_out = 16L, leave_one_batch_out = 4L
  )
  observed <- table(factor(out$axis, levels = names(expected)))
  stopifnot(
    nrow(out) == 35L, sum(!out$is_reference) == 34L,
    !anyDuplicated(out$variant_id), identical(as.integer(observed), unname(expected))
  )
  out
}

occupancy_sweep_resolution_object <- function(microglia_processed, resolution) {
  stopifnot(
    inherits(microglia_processed, "Seurat"), length(resolution) == 1L,
    is.finite(resolution), resolution > 0
  )
  obj <- microglia_processed
  resolution_col <- paste0("SCT_snn_res.", format(resolution, trim = TRUE, scientific = FALSE))
  if (!resolution_col %in% colnames(obj@meta.data)) {
    if (!length(SeuratObject::Graphs(obj))) {
      obj <- Seurat::FindNeighbors(
        obj, reduction = "harmony", dims = seq_len(20L), verbose = FALSE
      )
    }
    obj <- Seurat::FindClusters(
      obj, resolution = resolution, algorithm = 1L,
      random.seed = 42L, verbose = FALSE
    )
  }
  stopifnot(resolution_col %in% colnames(obj@meta.data))
  obj@meta.data[["microglia_clusters"]] <- factor(obj@meta.data[[resolution_col]])
  SeuratObject::Idents(obj) <- "microglia_clusters"
  obj
}

occupancy_sweep_membership_evaluation <- function(obj, expected_units) {
  stopifnot(inherits(obj, "Seurat"))
  required <- c("genotype_batch", "genotype", "batch", "microglia_subpopulation")
  stopifnot(all(required %in% colnames(obj@meta.data)))
  md <- obj@meta.data[, required, drop = FALSE]
  unit_map <- unique(md[c("genotype_batch", "genotype", "batch")])
  stopifnot(
    !anyDuplicated(unit_map$genotype_batch),
    identical(
      as.character(unit_map$genotype_batch),
      paste(as.character(unit_map$genotype), as.character(unit_map$batch), sep = "_")
    )
  )
  present_units <- sort(unique(as.character(md$genotype_batch)), method = "radix")
  expected_units <- sort(unique(as.character(expected_units)), method = "radix")
  missing_units <- setdiff(expected_units, present_units)
  extra_units <- setdiff(present_units, expected_units)
  expected_reason <- character()
  if (length(missing_units)) {
    expected_reason <- c(
      expected_reason,
      paste0("labeling removed every retained cell from expected unit(s): ",
             paste(missing_units, collapse = ", "))
    )
  }
  if (length(extra_units)) {
    expected_reason <- c(
      expected_reason,
      paste0("labeling introduced unexpected unit(s): ", paste(extra_units, collapse = ", "))
    )
  }

  prune <- obj@misc$microglia_prune
  dropped <- if (is.list(prune) && !is.null(prune$dropped)) {
    as.character(prune$dropped)
  } else {
    character()
  }
  base <- list(
    fit = NULL,
    membership_status = if (length(expected_reason)) "non_estimable" else "pending",
    membership_reason = if (length(expected_reason)) {
      paste(expected_reason, collapse = " | ")
    } else {
      NA_character_
    },
    n_cells = nrow(md),
    n_units = length(present_units),
    n_batches = length(unique(as.character(md$batch))),
    n_dropped_clusters = length(dropped),
    dropped_clusters = if (length(dropped)) paste(dropped, collapse = ",") else ""
  )
  if (length(expected_reason)) return(base)

  fit <- occupancy_from_membership(
    membership = md$microglia_subpopulation,
    meta = md[c("genotype_batch", "genotype", "batch")],
    estimators = occupancy_sweep_estimators(),
    margin = 0.10, n_perm = 9999L, seed = 614L
  )
  base$fit <- fit
  base$membership_status <- fit$membership_diagnostics$status
  base$membership_reason <- fit$membership_diagnostics$reason
  base$n_units <- fit$membership_diagnostics$n_units
  base
}

occupancy_sweep_annotation_evaluation <- function(
    microglia_processed, symbol_map, expected_units,
    primary_resolution = 0.4, id_floor = 0.15, mglike_floor = 0.30,
    tol = 0.10, amb_floor = 0.10, require_no_prune = FALSE) {
  obj <- if (identical(primary_resolution, 0.4)) {
    microglia_processed
  } else {
    occupancy_sweep_resolution_object(microglia_processed, primary_resolution)
  }
  annotated <- annotate_microglia(
    obj, symbol_map,
    id_floor = id_floor, mglike_floor = mglike_floor,
    tol = tol, amb_floor = amb_floor
  )
  if (require_no_prune) {
    stopifnot(
      is.list(annotated@misc$microglia_prune),
      length(annotated@misc$microglia_prune$dropped) == 0L,
      annotated@misc$microglia_prune$n_dropped == 0L
    )
  }
  out <- occupancy_sweep_membership_evaluation(annotated, expected_units)
  rm(annotated, obj)
  invisible(gc(FALSE))
  out
}


occupancy_sweep_safe_annotation_evaluation <- function(...) {
  tryCatch(
    occupancy_sweep_annotation_evaluation(...),
    error = function(e) {
      reason <- paste0(
        "variant labeling/construction failed: ", conditionMessage(e)
      )
      list(
        fit = NULL,
        membership_status = "non_estimable",
        membership_reason = reason,
        n_cells = NA_integer_,
        n_units = NA_integer_,
        n_batches = NA_integer_,
        n_dropped_clusters = NA_integer_,
        dropped_clusters = NA_character_
      )
    }
  )
}

occupancy_sweep_coverage_evaluation <- function(unit_coverage, unit_meta) {
  stopifnot(
    is.data.frame(unit_coverage), is.data.frame(unit_meta),
    identical(as.character(unit_coverage$genotype_batch), rownames(unit_meta))
  )
  fit <- occupancy_family(
    unit_coverage, unit_meta,
    estimators = occupancy_sweep_estimators(),
    margin = 0.10, n_perm = 9999L, seed = 614L
  )
  list(
    fit = fit,
    membership_status = "reference_aggregate",
    membership_reason = NA_character_,
    n_cells = sum(unit_coverage$n_retained),
    n_units = nrow(unit_coverage),
    n_batches = length(unique(as.character(unit_meta$batch))),
    n_dropped_clusters = NA_integer_,
    dropped_clusters = NA_character_
  )
}

occupancy_sweep_estimator_diagnostic <- function(evaluation, estimator) {
  if (is.null(evaluation$fit)) {
    stopifnot(
      identical(evaluation$membership_status, "non_estimable"),
      length(evaluation$membership_reason) == 1L,
      !is.na(evaluation$membership_reason), nzchar(evaluation$membership_reason)
    )
    return(list(status = "estimator_failed", reason = evaluation$membership_reason))
  }
  diagnostic <- evaluation$fit$diagnostics$estimators[[estimator]]
  stopifnot(
    is.list(diagnostic), diagnostic$status %in% c("ok", "estimator_failed"),
    length(diagnostic$reason) == 1L
  )
  list(status = diagnostic$status, reason = diagnostic$reason)
}

occupancy_sweep_variant_row <- function(spec, evaluation) {
  stopifnot(is.data.frame(spec), nrow(spec) == 1L, is.list(evaluation))
  e1_diagnostic <- occupancy_sweep_estimator_diagnostic(evaluation, "beta_binomial")
  e2_diagnostic <- occupancy_sweep_estimator_diagnostic(evaluation, "empirical_logit")
  e3_diagnostic <- occupancy_sweep_estimator_diagnostic(evaluation, "simple_proportion")

  e1 <- e2 <- e3 <- NULL
  permutation <- NULL
  if (identical(e1_diagnostic$status, "ok")) {
    e1 <- occupancy_sweep_interaction_row(
      evaluation$fit$probability_contrasts, "E1 probability contrasts"
    )
  }
  if (identical(e2_diagnostic$status, "ok")) {
    e2 <- occupancy_sweep_interaction_row(
      evaluation$fit$empirical_logit, "E2 empirical logit"
    )
    permutation <- evaluation$fit$permutation
    stopifnot(
      is.data.frame(permutation), nrow(permutation) == 1L,
      identical(permutation$contrast, "interaction")
    )
  }
  if (identical(e3_diagnostic$status, "ok")) {
    e3 <- occupancy_sweep_interaction_row(
      evaluation$fit$simple_proportion, "E3 simple proportion"
    )
  }

  data.frame(
    spec,
    n_cells = as.integer(evaluation$n_cells),
    n_units = as.integer(evaluation$n_units),
    n_batches = as.integer(evaluation$n_batches),
    n_dropped_clusters = as.integer(evaluation$n_dropped_clusters),
    dropped_clusters = evaluation$dropped_clusters,
    membership_status = evaluation$membership_status,
    membership_reason = evaluation$membership_reason,
    e1_status = e1_diagnostic$status,
    e1_reason = e1_diagnostic$reason,
    e1_estimate = if (is.null(e1)) NA_real_ else unname(e1$estimate),
    e1_ci_l = if (is.null(e1)) NA_real_ else unname(e1$ci_l),
    e1_ci_r = if (is.null(e1)) NA_real_ else unname(e1$ci_r),
    e1_fdr_zero = if (is.null(e1)) NA_real_ else unname(e1$fdr_zero),
    e1_fdr_minimum = if (is.null(e1)) NA_real_ else unname(e1$fdr_minimum),
    e2_status = e2_diagnostic$status,
    e2_reason = e2_diagnostic$reason,
    e2_coef = if (is.null(e2)) NA_real_ else unname(e2$coef),
    e2_ci_l = if (is.null(e2)) NA_real_ else unname(e2$ci_l),
    e2_ci_r = if (is.null(e2)) NA_real_ else unname(e2$ci_r),
    e2_fdr = if (is.null(e2)) NA_real_ else unname(e2$fdr),
    e2_perm_p = if (is.null(permutation)) NA_real_ else unname(permutation$perm_p),
    e3_status = e3_diagnostic$status,
    e3_reason = e3_diagnostic$reason,
    e3_estimate = if (is.null(e3)) NA_real_ else unname(e3$coef),
    e3_ci_l = if (is.null(e3)) NA_real_ else unname(e3$ci_l),
    e3_ci_r = if (is.null(e3)) NA_real_ else unname(e3$ci_r),
    e3_fdr = if (is.null(e3)) NA_real_ else unname(e3$fdr),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

occupancy_sweep_failure_inventory <- function(variant_table) {
  estimators <- c(e1 = "beta_binomial", e2 = "empirical_logit", e3 = "simple_proportion")
  rows <- list()
  for (prefix in names(estimators)) {
    failed <- variant_table[[paste0(prefix, "_status")]] == "estimator_failed"
    if (!any(failed)) next
    reason <- variant_table[[paste0(prefix, "_reason")]][failed]
    stopifnot(!anyNA(reason), all(nzchar(reason)))
    rows[[length(rows) + 1L]] <- data.frame(
      variant_id = variant_table$variant_id[failed],
      axis = variant_table$axis[failed],
      estimator = unname(estimators[[prefix]]),
      reason = reason,
      membership_reason = variant_table$membership_reason[failed],
      stringsAsFactors = FALSE
    )
  }
  if (!length(rows)) {
    return(data.frame(
      variant_id = character(), axis = character(), estimator = character(),
      reason = character(), membership_reason = character(),
      stringsAsFactors = FALSE
    ))
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

occupancy_sweep_concordance <- function(variant_table, estimator) {
  stopifnot(estimator %in% c("e2", "e3"))
  status_col <- paste0(estimator, "_status")
  estimate_col <- if (identical(estimator, "e2")) "e2_coef" else "e3_estimate"
  comparable <- variant_table$e1_status == "ok" & variant_table[[status_col]] == "ok"
  concordant <- rep(NA, nrow(variant_table))
  concordant[comparable] <-
    sign(variant_table[[estimate_col]][comparable]) ==
    sign(variant_table$e1_estimate[comparable])
  list(
    values = concordant,
    summary = list(
      n_comparable = sum(comparable),
      n_concordant = sum(concordant[comparable]),
      n_discordant = sum(!concordant[comparable]),
      discordant_variants = variant_table$variant_id[comparable & !concordant]
    )
  )
}

run_occupancy_robustness <- function(
    microglia_processed, microglia_annotated, microglia_state_substrate, symbol_map) {
  stopifnot(
    inherits(microglia_processed, "Seurat"),
    inherits(microglia_annotated, "Seurat"),
    is.list(microglia_state_substrate),
    is.data.frame(microglia_state_substrate$unit_coverage),
    is.data.frame(microglia_state_substrate$unit_meta),
    is.data.frame(symbol_map)
  )
  unit_coverage <- microglia_state_substrate$unit_coverage
  unit_meta <- microglia_state_substrate$unit_meta
  expected_units <- as.character(unit_coverage$genotype_batch)
  expected_batches <- as.character(unit_meta$batch)
  stopifnot(
    nrow(unit_coverage) == 16L, nrow(unit_meta) == 16L,
    identical(expected_units, rownames(unit_meta)),
    identical(
      expected_units,
      paste(as.character(unit_meta$genotype), as.character(unit_meta$batch), sep = "_")
    )
  )

  # Reference oracle: aggregate path -> frozen interaction; cached membership -> exact aggregate;
  # default re-annotation -> exact cached labels and the same occupancy family.
  substrate_evaluation <- occupancy_sweep_coverage_evaluation(unit_coverage, unit_meta)
  substrate_interaction <- occupancy_sweep_interaction_vector(substrate_evaluation$fit)
  frozen_interaction <- occupancy_sweep_frozen_interaction()
  frozen_max_abs_diff <- occupancy_sweep_max_abs_diff(
    substrate_interaction, frozen_interaction
  )

  cached_md <- microglia_annotated@meta.data
  cached_aggregation <- membership_to_unit_coverage(
    cached_md$microglia_subpopulation,
    cached_md[c("genotype_batch", "genotype", "batch")]
  )
  cached_evaluation <- occupancy_sweep_membership_evaluation(
    microglia_annotated, expected_units
  )
  cached_interaction <- occupancy_sweep_interaction_vector(cached_evaluation$fit)
  membership_max_abs_diff <- occupancy_sweep_max_abs_diff(
    cached_interaction, substrate_interaction
  )

  regenerated <- annotate_microglia(
    microglia_processed, symbol_map,
    id_floor = 0.15, mglike_floor = 0.30, tol = 0.10, amb_floor = 0.10
  )
  regenerated_cells_identical <- identical(
    colnames(regenerated), colnames(microglia_annotated)
  )
  regenerated_labels_identical <- identical(
    regenerated@meta.data$microglia_subpopulation,
    microglia_annotated@meta.data$microglia_subpopulation
  )
  reference_evaluation <- occupancy_sweep_membership_evaluation(
    regenerated, expected_units
  )
  regenerated_interaction <- occupancy_sweep_interaction_vector(reference_evaluation$fit)
  regenerated_max_abs_diff <- occupancy_sweep_max_abs_diff(
    regenerated_interaction, substrate_interaction
  )

  stopifnot(
    frozen_max_abs_diff <= 1e-8,
    identical(round(unname(substrate_interaction["e1_estimate"]), 7L), 0.1741141),
    membership_max_abs_diff <= 1e-8,
    identical(cached_aggregation$diagnostics$status, "estimable"),
    identical(cached_aggregation$unit_coverage, unit_coverage),
    identical(cached_aggregation$unit_meta, unit_meta),
    regenerated_cells_identical,
    regenerated_labels_identical,
    regenerated_max_abs_diff <= 1e-8
  )

  grid <- occupancy_sweep_grid(expected_units, expected_batches)
  rows <- vector("list", nrow(grid))
  rows[[1L]] <- occupancy_sweep_variant_row(grid[1L, , drop = FALSE], reference_evaluation)
  rm(regenerated, cached_md, cached_evaluation, reference_evaluation)
  invisible(gc(FALSE))

  for (i in seq.int(2L, nrow(grid))) {
    spec <- grid[i, , drop = FALSE]
    evaluation <- switch(
      spec$axis,
      resolution = occupancy_sweep_safe_annotation_evaluation(
        microglia_processed, symbol_map, expected_units,
        primary_resolution = spec$primary_resolution
      ),
      id_floor = occupancy_sweep_safe_annotation_evaluation(
        microglia_processed, symbol_map, expected_units,
        id_floor = spec$id_floor
      ),
      mglike_floor = occupancy_sweep_safe_annotation_evaluation(
        microglia_processed, symbol_map, expected_units,
        mglike_floor = spec$mglike_floor
      ),
      no_prune = occupancy_sweep_safe_annotation_evaluation(
        microglia_processed, symbol_map, expected_units,
        id_floor = 0, mglike_floor = 0, require_no_prune = TRUE
      ),
      tol = occupancy_sweep_safe_annotation_evaluation(
        microglia_processed, symbol_map, expected_units,
        tol = spec$tol
      ),
      amb_floor = occupancy_sweep_safe_annotation_evaluation(
        microglia_processed, symbol_map, expected_units,
        amb_floor = spec$amb_floor
      ),
      leave_one_unit_out = {
        keep <- as.character(unit_coverage$genotype_batch) != spec$dropped_unit
        uc <- unit_coverage[keep, , drop = FALSE]
        um <- unit_meta[as.character(uc$genotype_batch), , drop = FALSE]
        occupancy_sweep_coverage_evaluation(uc, um)
      },
      leave_one_batch_out = {
        keep <- as.character(unit_coverage$batch) != spec$dropped_batch
        uc <- droplevels(unit_coverage[keep, , drop = FALSE])
        um <- droplevels(unit_meta[as.character(uc$genotype_batch), , drop = FALSE])
        occupancy_sweep_coverage_evaluation(uc, um)
      },
      stop("unknown occupancy robustness axis: ", spec$axis, call. = FALSE)
    )
    rows[[i]] <- occupancy_sweep_variant_row(spec, evaluation)
    rm(evaluation)
    invisible(gc(FALSE))
  }

  variant_table <- do.call(rbind, rows)
  rownames(variant_table) <- NULL
  stopifnot(
    nrow(variant_table) == 35L,
    sum(!variant_table$is_reference) == 34L,
    !anyDuplicated(variant_table$variant_id),
    all(variant_table$e1_status %in% c("ok", "estimator_failed")),
    all(variant_table$e2_status %in% c("ok", "estimator_failed")),
    all(variant_table$e3_status %in% c("ok", "estimator_failed"))
  )

  e2_concordance <- occupancy_sweep_concordance(variant_table, "e2")
  e3_concordance <- occupancy_sweep_concordance(variant_table, "e3")
  variant_table$e2_sign_concordant <- e2_concordance$values
  variant_table$e3_sign_concordant <- e3_concordance$values

  e1_ok <- variant_table$e1_status == "ok" & is.finite(variant_table$e1_estimate) &
    is.finite(variant_table$e1_fdr_zero) & is.finite(variant_table$e1_fdr_minimum)
  e1_pass <- rep(FALSE, nrow(variant_table))
  e1_pass[e1_ok] <- variant_table$e1_estimate[e1_ok] > 0 &
    variant_table$e1_fdr_zero[e1_ok] <= 0.05
  verdict <- if (all(e1_pass)) "ROBUST-POSITIVE" else "FRAGILE"

  direction_tipping <- e1_ok & !variant_table$is_reference &
    variant_table$e1_estimate <= 0
  zero_null_tipping <- e1_ok & !variant_table$is_reference &
    variant_table$e1_fdr_zero > 0.05
  tipping_rows <- direction_tipping | zero_null_tipping
  tipping_set <- variant_table[
    tipping_rows,
    c("variant_id", "axis", "setting", "e1_estimate", "e1_fdr_zero"),
    drop = FALSE
  ]
  tipping_set$direction_tipping <- direction_tipping[tipping_rows]
  tipping_set$zero_null_tipping <- zero_null_tipping[tipping_rows]
  rownames(tipping_set) <- NULL

  margin_resolved <- variant_table[
    e1_ok & variant_table$e1_fdr_minimum <= 0.05,
    c("variant_id", "axis", "setting", "e1_estimate", "e1_fdr_minimum"),
    drop = FALSE
  ]
  within_margin <- variant_table[
    e1_ok & abs(variant_table$e1_estimate) <= 0.10,
    c("variant_id", "axis", "setting", "e1_estimate", "e1_fdr_minimum"),
    drop = FALSE
  ]
  rownames(margin_resolved) <- rownames(within_margin) <- NULL

  failure_inventory <- occupancy_sweep_failure_inventory(variant_table)
  e1_values <- variant_table$e1_estimate[e1_ok]
  minimum_fdr_values <- variant_table$e1_fdr_minimum[e1_ok]
  stopifnot(length(e1_values) >= 1L, length(minimum_fdr_values) >= 1L)

  reference_diagnostics <- list(
    frozen_interaction = frozen_interaction,
    substrate_interaction = substrate_interaction,
    substrate_vs_frozen_max_abs_diff = frozen_max_abs_diff,
    cached_membership_vs_substrate_max_abs_diff = membership_max_abs_diff,
    cached_membership_aggregation_identical = TRUE,
    reannotation_cells_identical = regenerated_cells_identical,
    reannotation_labels_identical = regenerated_labels_identical,
    reannotation_vs_substrate_max_abs_diff = regenerated_max_abs_diff
  )
  out <- list(
    schema = "p7_occupancy_robustness_v1",
    reference_diagnostics = reference_diagnostics,
    variant_table = variant_table,
    e1_range = stats::setNames(range(e1_values), c("min", "max")),
    tipping_set = tipping_set,
    margin = list(
      fdr_minimum_range = stats::setNames(
        range(minimum_fdr_values), c("min", "max")
      ),
      resolved_variants = margin_resolved,
      within_margin_variants = within_margin
    ),
    estimator_failed = failure_inventory,
    cross_estimator_concordance = list(
      e2_empirical_logit = e2_concordance$summary,
      e3_simple_proportion = e3_concordance$summary
    ),
    verdict = verdict,
    audit = list(
      protocol = ".agent/p7_dam_occupancy_prereg.md frozen at 7755c9f",
      reference_settings = c(
        primary_resolution = 0.4, id_floor = 0.15, mglike_floor = 0.30,
        tol = 0.10, amb_floor = 0.10
      ),
      margin = 0.10,
      n_perm = 9999L,
      occupancy_seed = 614L,
      clustering_seed = 42L,
      estimator_names = occupancy_sweep_estimators(),
      n_variants = nrow(variant_table),
      n_nonreference = sum(!variant_table$is_reference),
      n_estimator_attempts = 3L * nrow(variant_table),
      n_estimator_failed = nrow(failure_inventory),
      optional_ucell_rescoring = "not run; preregistered optional secondary",
      e1_failure_blocks_robust_positive = TRUE,
      parent_isolated = FALSE,
      in_memory_bytes = NA_real_,
      serialized_bytes = NA_real_,
      max_target_bytes = 256 * 1024
    )
  )

  stopifnot(!state_substrate_contains_parent(out))
  out$audit$parent_isolated <- TRUE
  out$audit$in_memory_bytes <- as.numeric(utils::object.size(out))
  out$audit$serialized_bytes <- as.numeric(length(qs2::qs_serialize(out)))
  stopifnot(
    !state_substrate_contains_parent(out),
    out$audit$in_memory_bytes <= out$audit$max_target_bytes,
    out$audit$serialized_bytes <= out$audit$max_target_bytes,
    out$audit$n_variants == 35L,
    out$audit$n_nonreference == 34L,
    out$audit$n_estimator_attempts == 105L
  )
  out
}

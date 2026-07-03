# Inline figure-data contracts for the post-report visual-density pass. These
# builders return compact, qmd-ready data frames/lists; they do not draw plots.

.fig_contrasts <- function() {
  c("tau_alone", "nlgf_in_maptki", "nlgf_in_p301s", "tau_in_nlgf", "interaction")
}

.fig_focus_contrasts <- function() {
  c("interaction", "nlgf_in_maptki", "nlgf_in_p301s", "tau_in_nlgf")
}

figure_manifest <- function(chapter = NULL) {
  out <- data.frame(
    figure_id = c(
      "fig-story-core",
      "fig-story-mechanism-crossmodality",
      "fig-microglia-umap-substate",
      "fig-microglia-score-triptych",
      "fig-microglia-unit-composition",
      "fig-microglia-score-distribution",
      "fig-microglia-whole-volcano",
      "fig-microglia-substate-audit",
      "fig-microglia-substate-volcano",
      "fig-trajectory-pt-density",
      "fig-trajectory-unit-pt-dam",
      "fig-trajectory-kitagawa-forest",
      "fig-trajectory-audit",
      "fig-mechanism-project-pathway",
      "fig-mechanism-go-dotplot",
      "fig-mechanism-tf-focus",
      "fig-mechanism-kinase-heatmap",
      "fig-crossmodality-amyloid-response",
      "fig-crossmodality-synaptic-clearance",
      "fig-crossmodality-interaction-boundary",
      "fig-crossmodality-geomx-volcano",
      "fig-crossmodality-geomx-sensitivity",
      "fig-crossmodality-phospho-correction"
    ),
    chapter = c(
      rep("story", 2),
      rep("microglia", 7),
      rep("trajectory", 4),
      rep("mechanism", 4),
      rep("crossmodality", 6)
    ),
    target = c(
      rep("story_figures", 2),
      rep("microglia_figures", 7),
      rep("trajectory_figures", 4),
      rep("mechanism_figures", 4),
      rep("crossmodality_figures", 6)
    ),
    slot = c(
      "core_story",
      "mechanism_crossmodality",
      "umap_by_substate",
      "score_triptych",
      "unit_composition",
      "score_distribution",
      "whole_de_volcano",
      "substate_fit_audit",
      "substate_de_volcano",
      "pt_density",
      "unit_pt_vs_dam",
      "kitagawa_forest",
      "trajectory_audit",
      "project_pathway_heatmap",
      "go_top_dotplot",
      "tf_focus",
      "kinase_heatmap",
      "amyloid_response_plate",
      "synaptic_clearance_plate",
      "interaction_boundary_plate",
      "geomx_volcano",
      "geomx_sensitivity",
      "phospho_raw_corrected"
    ),
    stringsAsFactors = FALSE
  )
  stopifnot(!anyDuplicated(out$figure_id), !any(grepl("_", out$figure_id, fixed = TRUE)))
  if (!is.null(chapter)) out <- out[out$chapter %in% chapter, , drop = FALSE]
  rownames(out) <- NULL
  out
}

.fig_story_dam_response <- function(composition_results) {
  stopifnot(is.list(composition_results), !is.null(composition_results$counts))
  comp <- .fig_unit_composition(composition_results$counts)
  dam <- comp[comp$substate == "DAM", , drop = FALSE]
  dam$tau_background <- factor(ifelse(grepl("P301S", as.character(dam$genotype)),
                                      "P301S", "MAPTKI"),
                               levels = c("MAPTKI", "P301S"))
  dam$amyloid <- factor(ifelse(startsWith(as.character(dam$genotype), "NLGF"),
                               "NLGF+", "NLGF-"),
                        levels = c("NLGF-", "NLGF+"))
  means <- stats::aggregate(proportion ~ tau_background + amyloid, data = dam, FUN = mean)
  names(means)[3L] <- "mean"
  se <- stats::aggregate(proportion ~ tau_background + amyloid, data = dam,
                         FUN = function(x) if (length(x) > 1L) stats::sd(x) / sqrt(length(x)) else 0)
  names(se)[3L] <- "se"
  means <- merge(means, se, by = c("tau_background", "amyloid"), all.x = TRUE, sort = FALSE)
  means$lo <- pmax(0, means$mean - means$se)
  means$hi <- pmin(1, means$mean + means$se)
  .fig_assert_finite(dam, c("n_cells", "unit_total", "proportion"), "story DAM unit response")
  .fig_assert_finite(means, c("mean", "se", "lo", "hi"), "story DAM mean response")
  list(unit = dam, mean = means)
}

.fig_story_de_counts <- function(pb_de_microglia) {
  stopifnot(is.list(pb_de_microglia), is.list(pb_de_microglia$top))
  alpha <- pb_de_microglia$thresholds$fdr %||% 0.10
  lfc <- pb_de_microglia$thresholds$lfc %||% 0
  counts <- .fig_bind(lapply(names(pb_de_microglia$top), function(cn) {
    tt <- pb_de_microglia$top[[cn]]
    .fig_require_cols(tt, c("logFC", "adj.P.Val"), paste0("pb_de_microglia$top$", cn))
    sig <- is.finite(tt$adj.P.Val) & tt$adj.P.Val < alpha &
      is.finite(tt$logFC) & abs(tt$logFC) > lfc
    data.frame(
      contrast = cn,
      n = sum(sig),
      up = sum(sig & tt$logFC > 0),
      down = sum(sig & tt$logFC < 0),
      stringsAsFactors = FALSE
    )
  }))
  signed <- rbind(
    data.frame(contrast = counts$contrast, direction = "up", n = counts$up,
               n_signed = counts$up, stringsAsFactors = FALSE),
    data.frame(contrast = counts$contrast, direction = "down", n = counts$down,
               n_signed = -counts$down, stringsAsFactors = FALSE)
  )
  signed$direction <- factor(signed$direction, levels = c("down", "up"))
  stage <- data.frame(n_interaction_stageR = NA_real_,
                      median_abs_lfc = NA_real_,
                      stringsAsFactors = FALSE)
  if (is.list(pb_de_microglia$stageR) && is.matrix(pb_de_microglia$stageR$stage_padj) &&
      "interaction" %in% colnames(pb_de_microglia$stageR$stage_padj) &&
      "interaction" %in% names(pb_de_microglia$top)) {
    hit <- pb_de_microglia$stageR$stage_padj[, "interaction"] <= pb_de_microglia$stageR$alpha
    hit[is.na(hit)] <- FALSE
    tt <- pb_de_microglia$top$interaction
    .fig_require_cols(tt, c("gene", "logFC"), "pb_de_microglia$top$interaction")
    idx <- match(rownames(pb_de_microglia$stageR$stage_padj)[hit], tt$gene)
    l <- abs(tt$logFC[idx])
    stage$n_interaction_stageR <- sum(hit)
    stage$median_abs_lfc <- if (length(l)) stats::median(l, na.rm = TRUE) else NA_real_
  }
  .fig_assert_finite(counts, c("n", "up", "down"), "story DE counts")
  .fig_assert_finite(signed, c("n", "n_signed"), "story signed DE counts")
  list(counts = counts, signed = signed, stageR = stage)
}

.fig_story_trajectory <- function(trajectory_report) {
  stopifnot(is.list(trajectory_report), is.data.frame(trajectory_report$interaction))
  keep <- c("mean_pt", "comp_cf", "progression_cf", "within_homeostatic")
  x <- trajectory_report$interaction[trajectory_report$interaction$measure %in% keep, , drop = FALSE]
  .fig_require_cols(x, c("measure", "coef", "ci_l", "ci_r", "p_value", "fdr"),
                    "trajectory_report$interaction")
  x$measure_label <- c(
    mean_pt = "mean position",
    comp_cf = "composition",
    progression_cf = "progression",
    within_homeostatic = "within-homeostatic"
  )[x$measure]
  x$measure_label <- factor(x$measure_label,
                            levels = rev(c("mean position", "composition",
                                           "progression", "within-homeostatic")))
  x$supported <- is.finite(x$fdr) & x$fdr < 0.10
  .fig_assert_finite(x, c("coef", "ci_l", "ci_r", "p_value", "fdr"),
                     "story trajectory decomposition")
  x
}

.fig_story_mechanism <- function(mechanism_report, alpha = 0.10) {
  stopifnot(is.list(mechanism_report), is.data.frame(mechanism_report$tf_highlights),
            is.list(mechanism_report$nfkb), is.list(mechanism_report$kinase))
  tf <- mechanism_report$tf_highlights
  .fig_require_cols(tf, c("population", "source", "contrast", "score", "fdr"),
                    "mechanism_report$tf_highlights")
  myc <- tf[tf$source == "Myc" & tf$contrast == "interaction" &
              tf$population %in% c("whole_microglia", "DAM"), , drop = FALSE]
  stopifnot(nrow(myc) >= 1L)
  myc <- data.frame(
    track = "Myc TF",
    item = ifelse(myc$population == "whole_microglia", "whole MG", "DAM"),
    contrast = myc$contrast,
    score = myc$score,
    fdr = myc$fdr,
    supported = is.finite(myc$fdr) & myc$fdr < alpha,
    stringsAsFactors = FALSE
  )

  nf <- mechanism_report$nfkb$table
  .fig_require_cols(nf, c("test", "score", "primary_test", "primary_family_fdr"),
                    "mechanism_report$nfkb$table")
  nf <- nf[nf$primary_test %in% TRUE & nf$contrast == "interaction", , drop = FALSE]
  stopifnot(nrow(nf) >= 1L)
  nf <- data.frame(
    track = "NF-kB gate",
    item = c(tf_family = "TF family", target_gsea = "target GSEA")[nf$test],
    contrast = "interaction",
    score = nf$score,
    fdr = nf$primary_family_fdr,
    supported = FALSE,
    stringsAsFactors = FALSE
  )

  kin <- mechanism_report$kinase$table
  .fig_require_cols(kin, c("source", "contrast", "score", "fdr", "significant"),
                    "mechanism_report$kinase$table")
  gsk <- kin[kin$source == "Gsk3b" & kin$contrast %in% c("interaction", "tau_in_nlgf"),
             , drop = FALSE]
  stopifnot(nrow(gsk) >= 1L)
  gsk <- data.frame(
    track = "Gsk3b kinase",
    item = ifelse(gsk$contrast == "interaction", "interaction", "tau in NLGF"),
    contrast = gsk$contrast,
    score = gsk$score,
    fdr = gsk$fdr,
    supported = gsk$significant %in% TRUE,
    stringsAsFactors = FALSE
  )
  out <- rbind(myc, nf, gsk)
  out$track <- factor(out$track, levels = c("Myc TF", "NF-kB gate", "Gsk3b kinase"))
  out$neg_log10_fdr <- -log10(pmax(out$fdr, 1e-300))
  stopifnot(!anyNA(out$item), !anyNA(out$track))
  .fig_assert_finite(out, c("score", "fdr", "neg_log10_fdr"), "story mechanism focus")
  out
}

.fig_story_pathways <- function(crossmodality_figures, axes = c("DAM", "synaptic", "clearance",
                                                                "antigen_presentation", "NFkB")) {
  stopifnot(is.list(crossmodality_figures),
            is.data.frame(crossmodality_figures$four_modality_pathways))
  x <- crossmodality_figures$four_modality_pathways
  .fig_require_cols(x, c("axis", "contrast", "n_modalities_present", "n_modalities_sig",
                         "n_evidence_groups_sig", "direction", "rank_score"),
                    "crossmodality_figures$four_modality_pathways")
  x <- x[x$axis %in% axes, , drop = FALSE]
  stopifnot(nrow(x) > 0L)
  x$axis_label <- c(
    DAM = "DAM",
    synaptic = "synaptic",
    clearance = "clearance",
    antigen_presentation = "antigen presentation",
    NFkB = "NF-kB"
  )[x$axis]
  x$axis_label <- factor(x$axis_label,
                         levels = rev(c("DAM", "synaptic", "clearance",
                                        "antigen presentation", "NF-kB")))
  .fig_assert_finite(x, c("n_modalities_present", "n_modalities_sig",
                          "n_evidence_groups_sig", "rank_score"),
                     "story pathway axes")
  x
}

.fig_story_clearance <- function(crossmodality_report) {
  stopifnot(is.list(crossmodality_report), is.list(crossmodality_report$clearance),
            is.data.frame(crossmodality_report$clearance$pair_support))
  x <- crossmodality_report$clearance$pair_support
  .fig_require_cols(x, c("pair", "contrast", "status", "coherent_supported_modalities"),
                    "crossmodality_report$clearance$pair_support")
  x$n_supported_modalities <- vapply(x$coherent_supported_modalities, function(z) {
    length(.fig_semicolon_tokens(z))
  }, integer(1), USE.NAMES = FALSE)
  x$supported <- x$status == "earned"
  x$pair_label <- gsub("_", "-", x$pair, fixed = TRUE)
  .fig_assert_finite(x, "n_supported_modalities", "story clearance support")
  x
}

.fig_story_clearance_effects <- function(crossmodality_figures) {
  stopifnot(is.list(crossmodality_figures),
            is.data.frame(crossmodality_figures$axis_effect_spine))
  x <- .fig_axis_effect_plot_rows(
    crossmodality_figures$axis_effect_spine,
    axes = "clearance",
    contrasts = c("nlgf_in_maptki", "nlgf_in_p301s"),
    feature_ids = c("clearance:Apoe", "clearance:Trem2"),
    modality_classes = names(.fig_four_modality_classes())
  )
  stopifnot(nrow(x) > 0L)
  x
}

story_figure_data <- function(qc_figures, composition_results, pb_de_microglia,
                              trajectory_report, mechanism_report,
                              crossmodality_report, crossmodality_figures,
                              alpha = 0.10) {
  stopifnot(is.list(qc_figures), is.list(composition_results),
            is.list(pb_de_microglia), is.list(trajectory_report),
            is.list(mechanism_report), is.list(crossmodality_report),
            is.list(crossmodality_figures))
  sample_counts <- qc_figures$study_design$sample_counts
  .fig_require_cols(sample_counts, c("genotype", "modality", "n"),
                    "qc_figures$study_design$sample_counts")
  dam <- .fig_story_dam_response(composition_results)
  de <- .fig_story_de_counts(pb_de_microglia)
  traj <- .fig_story_trajectory(trajectory_report)
  mech <- .fig_story_mechanism(mechanism_report, alpha)
  path <- .fig_story_pathways(crossmodality_figures)
  clr <- .fig_story_clearance(crossmodality_report)
  clr_effects <- .fig_story_clearance_effects(crossmodality_figures)
  out <- list(
    manifest = figure_manifest("story"),
    sample_counts = sample_counts,
    dam_response = dam,
    de_counts = de,
    trajectory = traj,
    mechanism = mech,
    pathway_axes = path,
    clearance = clr,
    clearance_effects = clr_effects,
    core_story = list(dam_response = dam, de_counts = de, trajectory = traj),
    mechanism_crossmodality = list(mechanism = mech, clearance = clr,
                                   clearance_effects = clr_effects),
    provenance = list(
      source_targets = c("qc_figures", "composition_results", "pb_de_microglia",
                         "trajectory_report", "mechanism_report",
                         "crossmodality_report", "crossmodality_figures"),
      alpha = alpha,
      contract = "compact story-plate figure data; no new inference"
    )
  )
  .fig_assert_finite(out$sample_counts, "n", "story sample counts")
  out
}

visual_reduction_slot_map <- function(disposition = NULL) {
  pinned <- figure_manifest()
  out <- data.frame(
    manifest_slot = pinned$figure_id,
    target = pinned$target,
    slot = pinned$slot,
    disposition = "figure;caption",
    stringsAsFactors = FALSE
  )
  retained <- data.frame(
    manifest_slot = c(
      "fig-qc-depth",
      "fig-qc-fractions",
      "fig-microglia-umap",
      "fig-microglia-composition-shift",
      "fig-microglia-composition-forest",
      "fig-microglia-amyloid-volcano",
      "fig-microglia-dropout-audit",
      "fig-trajectory-pseudotime-shift",
      "fig-trajectory-decomposition",
      "fig-trajectory-concordance",
      "fig-mechanism-tf-interaction",
      "collapsed-qc-audit",
      "collapsed-microglia-audit",
      "collapsed-trajectory-audit",
      "section-nav"
    ),
    target = c(
      "qc_figures",
      "qc_figures",
      "microglia_report",
      "microglia_figures",
      "microglia_figures",
      "microglia_figures",
      "microglia_report",
      "trajectory_figures",
      "trajectory_figures",
      "trajectory_figures",
      "mechanism_figures",
      "qc_figures",
      "microglia_report",
      "trajectory_figures",
      "report"
    ),
    slot = c(
      "depth_distribution",
      "fraction_distribution",
      "cell_frame",
      "composition_shift",
      "composition_forest",
      "amyloid_volcano",
      "prune",
      "pseudotime_shift",
      "decomposition",
      "concordance",
      "tf_interaction",
      "audit_notes",
      "prune",
      "trajectory_audit",
      "heading"
    ),
    disposition = c(
      "figure;caption",
      "figure;caption",
      "caption",
      "figure;caption",
      "figure;caption",
      "caption",
      "caption",
      "figure;caption",
      "figure;caption",
      "figure;caption",
      "caption",
      "collapsed_audit",
      "collapsed_audit",
      "collapsed_audit",
      "keep"
    ),
    stringsAsFactors = FALSE
  )
  out <- rbind(out, retained)
  stopifnot(!anyDuplicated(out$manifest_slot))
  if (!is.null(disposition)) {
    keep <- vapply(out$disposition, function(x) {
      any(disposition %in% strsplit(x, ";", fixed = TRUE)[[1L]])
    }, logical(1))
    out <- out[keep, , drop = FALSE]
  }
  rownames(out) <- NULL
  out
}

visual_slot_coverage <- function(manifest, disposition = c("figure", "schematic"),
                                 slot_map = visual_reduction_slot_map()) {
  .fig_require_cols(manifest, c("disposition", "target_slot"), "prose replacement manifest")
  selected <- manifest[manifest$disposition %in% disposition, , drop = FALSE]
  if (nrow(selected) > 0L) {
    split_slots <- strsplit(selected$target_slot, ";", fixed = TRUE)
    slots <- sort(unique(unlist(split_slots, use.names = FALSE)), method = "radix")
    slots <- slots[nzchar(slots)]
  } else {
    slots <- character()
  }
  missing <- setdiff(slots, slot_map$manifest_slot)
  data.frame(
    n_slots = length(slots),
    n_missing = length(missing),
    missing = paste(missing, collapse = ";"),
    stringsAsFactors = FALSE
  )
}

.fig_require_cols <- function(x, cols, label) {
  if (!is.data.frame(x)) stop(label, " must be a data.frame", call. = FALSE)
  miss <- setdiff(cols, names(x))
  if (length(miss)) {
    stop(label, " missing columns: ", paste(miss, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

.fig_bind <- function(xs) {
  xs <- Filter(function(x) is.data.frame(x) && nrow(x) > 0L, xs)
  if (!length(xs)) return(data.frame())
  out <- do.call(rbind, xs)
  rownames(out) <- NULL
  out
}

.fig_object_dim <- function(x, label) {
  d <- tryCatch(dim(x), error = function(e) NULL)
  if ((is.null(d) || length(d) < 2L) && is.list(x) && !is.null(x$dims)) d <- x$dims
  if (is.null(d) || length(d) < 2L || !all(is.finite(d[1:2]))) {
    stop(label, " dimensions unavailable", call. = FALSE)
  }
  as.integer(d[1:2])
}

.fig_meta_data <- function(x, label) {
  if (is.data.frame(x)) return(x)
  if (is.list(x) && is.data.frame(x$meta.data)) return(x$meta.data)
  if (isS4(x) && "meta.data" %in% methods::slotNames(x)) return(x@meta.data)
  stop(label, " metadata unavailable", call. = FALSE)
}

.fig_metric_histogram <- function(df, metrics, labels, transform = identity,
                                  bins = 50L, lower = NULL) {
  pieces <- lapply(names(metrics), function(metric) {
    value <- transform(as.numeric(df[[metric]]))
    data.frame(
      genotype = df$genotype,
      metric = labels[[metric]],
      value = value,
      stringsAsFactors = FALSE
    )
  })
  long <- .fig_bind(pieces)
  .fig_histogram(long, "value", c("genotype", "metric"), bins = bins, lower = lower)
}

qc_figure_data <- function(microglia_seurat_raw, geomx, proteomics, phospho, sample_key) {
  md <- .fig_meta_data(microglia_seurat_raw, "microglia_seurat_raw")
  gd <- .fig_meta_data(geomx, "geomx")
  .fig_require_cols(md, c("genotype", "batch", "genotype_batch", "nCount_RNA",
                          "nFeature_RNA", "percent_mt", "percent_ribo",
                          "percent_malat1", "percent_contam"),
                    "microglia_seurat_raw@meta.data")
  .fig_require_cols(gd, "genotype", "geomx metadata")
  .fig_require_cols(sample_key, "genotype", "sample_key")
  microglia_dim <- .fig_object_dim(microglia_seurat_raw, "microglia_seurat_raw")
  geomx_dim <- .fig_object_dim(geomx, "geomx")
  proteomics_dim <- .fig_object_dim(proteomics, "proteomics")
  phospho_dim <- .fig_object_dim(phospho, "phospho")

  modality_table <- data.frame(
    modality = c("snRNAseq microglia", "GeoMx spatial", "proteomics PTM",
                 "phosphoproteomics PTM", "proteomics sample key"),
    feature_axis = c("genes", "genes", "sites", "sites", "rows"),
    n_features = c(microglia_dim[1L], geomx_dim[1L], proteomics_dim[1L],
                   phospho_dim[1L], nrow(sample_key)),
    sample_axis = c("nuclei", "AOIs", "columns", "columns", "24M samples"),
    n_samples = c(microglia_dim[2L], geomx_dim[2L], proteomics_dim[2L],
                  phospho_dim[2L], nrow(sample_key)),
    stringsAsFactors = FALSE
  )
  modality_table$shape <- paste(modality_table$n_features, modality_table$n_samples, sep = " x ")

  genotype_batch <- as.data.frame(table(
    genotype = factor(md$genotype, levels = genotype_levels),
    batch = md$batch
  ), stringsAsFactors = FALSE)
  genotype_batch$genotype <- factor(as.character(genotype_batch$genotype), levels = genotype_levels)
  genotype_batch$Freq <- as.numeric(genotype_batch$Freq)
  genotype_batch <- genotype_batch[order(genotype_batch$batch, genotype_batch$genotype,
                                         method = "radix"), , drop = FALSE]

  geomx_genotype <- as.data.frame(table(
    genotype = factor(gd$genotype, levels = genotype_levels)
  ), stringsAsFactors = FALSE)
  names(geomx_genotype)[2L] <- "n_aoi"
  geomx_genotype$n_aoi <- as.numeric(geomx_genotype$n_aoi)

  design_grid <- data.frame(
    genotype = genotype_levels,
    tau_background = ifelse(grepl("P301S", genotype_levels), "P301S", "MAPTKI"),
    amyloid = ifelse(startsWith(genotype_levels, "NLGF"), "NLGF+", "NLGF-"),
    label = c("MAPTKI", "P301S", "NLGF\nMAPTKI", "NLGF\nP301S"),
    stringsAsFactors = FALSE
  )
  design_grid$genotype <- factor(design_grid$genotype, levels = genotype_levels)
  design_grid$tau_background <- factor(design_grid$tau_background,
                                       levels = c("MAPTKI", "P301S"))
  design_grid$amyloid <- factor(design_grid$amyloid, levels = c("NLGF-", "NLGF+"))

  sn_counts <- stats::aggregate(Freq ~ genotype, data = genotype_batch, FUN = sum)
  names(sn_counts)[2L] <- "n"
  sn_counts$modality <- "snRNAseq nuclei"
  geo_counts <- geomx_genotype
  names(geo_counts)[2L] <- "n"
  geo_counts$modality <- "GeoMx AOIs"
  sample_genotype <- factor(as.character(sample_key$genotype), levels = genotype_levels)
  bulk_counts <- as.data.frame(table(genotype = sample_genotype), stringsAsFactors = FALSE)
  names(bulk_counts)[2L] <- "n"
  bulk_counts$n <- as.numeric(bulk_counts$n)
  bulk_counts$modality <- "24M bulk runs"
  sample_counts <- rbind(sn_counts[c("genotype", "modality", "n")],
                         geo_counts[c("genotype", "modality", "n")],
                         bulk_counts[c("genotype", "modality", "n")])
  sample_counts$genotype <- factor(as.character(sample_counts$genotype), levels = genotype_levels)
  sample_counts$modality <- factor(sample_counts$modality,
                                   levels = c("snRNAseq nuclei", "GeoMx AOIs", "24M bulk runs"))

  depth_metrics <- c(nCount_RNA = "total counts", nFeature_RNA = "detected genes")
  frac_metrics <- c(percent_mt = "mitochondrial %", percent_ribo = "ribosomal %",
                    percent_malat1 = "MALAT1 %", percent_contam = "contamination %")
  md$genotype <- factor(as.character(md$genotype), levels = genotype_levels)
  depth_distribution <- .fig_metric_histogram(md, depth_metrics, depth_metrics,
                                              transform = function(x) log10(pmax(x, 1)),
                                              bins = 50L, lower = 0)
  fraction_distribution <- .fig_metric_histogram(md, frac_metrics, frac_metrics,
                                                 bins = 50L, lower = 0)

  qc_metrics <- c("nCount_RNA", "nFeature_RNA", names(frac_metrics))
  bounds <- data.frame(
    metric = qc_metrics,
    lower = c(1, 1, 0, 0, 0, 0),
    upper = c(Inf, microglia_dim[1L], 100, 100, 100, 100),
    rail = c("UMIs >= 1", "genes >= 1 and <= assayed genes",
             rep("percentage 0-100", 4L)),
    stringsAsFactors = FALSE
  )
  bounds$obs_min <- vapply(qc_metrics, function(cl) min(md[[cl]], na.rm = TRUE), numeric(1))
  bounds$obs_max <- vapply(qc_metrics, function(cl) max(md[[cl]], na.rm = TRUE), numeric(1))
  bounds$n_na <- vapply(qc_metrics, function(cl) sum(is.na(md[[cl]])), integer(1))
  bounds$within <- bounds$n_na == 0 & bounds$obs_min >= bounds$lower & bounds$obs_max <= bounds$upper
  stopifnot(all(bounds$within),
            all(md$nFeature_RNA <= md$nCount_RNA),
            all(genotype_batch$Freq > 0),
            nrow(sample_key) == 16L,
            !anyNA(sample_genotype),
            !anyNA(sample_counts$genotype),
            all(table(sample_genotype) == nrow(sample_key) %/% length(genotype_levels)),
            setequal(as.character(unique(gd$genotype)), genotype_levels))

  out <- list(
    manifest = visual_reduction_slot_map("figure")[visual_reduction_slot_map("figure")$target == "qc_figures", ],
    study_design = list(genotype_grid = design_grid,
                        sample_counts = sample_counts),
    modality_table = modality_table,
    geomx_genotype = geomx_genotype,
    genotype_batch = genotype_batch,
    depth_distribution = depth_distribution,
    fraction_distribution = fraction_distribution,
    metric_bounds = bounds,
    audit_notes = data.frame(
      note = c("modality dimensions", "complete genotype-batch grid", "definitional metric bounds"),
      status = "pass",
      stringsAsFactors = FALSE
    ),
    provenance = list(
      source_targets = c("microglia_seurat_raw", "geomx", "proteomics", "phospho", "sample_key"),
      contract = "compact QC figure data; report layer need not read raw modality objects for QC visuals"
    )
  )
  .fig_assert_finite(out$modality_table, c("n_features", "n_samples"), "qc modality_table")
  .fig_assert_finite(out$study_design$sample_counts, "n", "qc study_design sample_counts")
  .fig_assert_finite(out$genotype_batch, "Freq", "qc genotype_batch")
  .fig_assert_finite(out$geomx_genotype, "n_aoi", "qc geomx_genotype")
  .fig_assert_finite(out$depth_distribution, c("x_mid", "n"), "qc depth_distribution")
  .fig_assert_finite(out$fraction_distribution, c("x_mid", "n"), "qc fraction_distribution")
  .fig_assert_finite(out$metric_bounds, c("lower", "obs_min", "obs_max", "n_na"), "qc metric_bounds")
  out
}

.fig_assert_nonempty <- function(x, label) {
  if (!is.data.frame(x) || nrow(x) < 1L) {
    stop(label, " is empty", call. = FALSE)
  }
  invisible(TRUE)
}

.fig_assert_finite <- function(x, cols, label) {
  .fig_require_cols(x, cols, label)
  bad <- cols[!vapply(cols, function(cl) all(is.finite(x[[cl]])), logical(1))]
  if (length(bad)) {
    stop(label, " has non-finite numeric columns: ", paste(bad, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

.fig_midpoints <- function(breaks) {
  (utils::head(breaks, -1L) + utils::tail(breaks, -1L)) / 2
}

.fig_breaks <- function(x, n = 50L, symmetric = FALSE, lower = NULL, upper_q = 0.995) {
  x <- x[is.finite(x)]
  if (!length(x)) stop("cannot bin an all-non-finite vector", call. = FALSE)
  if (symmetric) {
    hi <- as.numeric(stats::quantile(abs(x), upper_q, names = FALSE, na.rm = TRUE))
    hi <- max(hi, 1e-6)
    return(seq(-hi, hi, length.out = n + 1L))
  }
  lo <- if (is.null(lower)) min(x, na.rm = TRUE) else lower
  hi <- as.numeric(stats::quantile(x, upper_q, names = FALSE, na.rm = TRUE))
  hi <- max(hi, lo + 1e-6)
  seq(lo, hi, length.out = n + 1L)
}

.fig_clip <- function(x, breaks) {
  pmin(pmax(x, min(breaks)), max(breaks))
}

.fig_grid_bins <- function(x, y, group, x_bins = 60L, y_bins = 50L,
                           symmetric_x = TRUE, lower_y = 0) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]; group <- group[ok]
  xb <- .fig_breaks(x, n = x_bins, symmetric = symmetric_x)
  yb <- .fig_breaks(y, n = y_bins, symmetric = FALSE, lower = lower_y)
  bx <- cut(.fig_clip(x, xb), breaks = xb, include.lowest = TRUE, labels = FALSE)
  by <- cut(.fig_clip(y, yb), breaks = yb, include.lowest = TRUE, labels = FALSE)
  raw <- data.frame(group = group, bin_x = bx, bin_y = by, n = 1L,
                    stringsAsFactors = FALSE)
  out <- stats::aggregate(n ~ group + bin_x + bin_y, data = raw, FUN = sum)
  xm <- .fig_midpoints(xb); ym <- .fig_midpoints(yb)
  out$x_mid <- xm[out$bin_x]
  out$y_mid <- ym[out$bin_y]
  out$x_min <- xb[out$bin_x]
  out$x_max <- xb[out$bin_x + 1L]
  out$y_min <- yb[out$bin_y]
  out$y_max <- yb[out$bin_y + 1L]
  out[order(out$group, out$bin_x, out$bin_y, method = "radix"), , drop = FALSE]
}

.fig_histogram <- function(df, value_col, group_cols, bins = 50L, lower = NULL) {
  .fig_require_cols(df, c(value_col, group_cols), paste0("histogram(", value_col, ")"))
  x <- df[[value_col]]
  ok <- is.finite(x)
  df <- df[ok, , drop = FALSE]
  x <- x[ok]
  br <- .fig_breaks(x, n = bins, lower = lower, symmetric = FALSE)
  b <- cut(.fig_clip(x, br), breaks = br, include.lowest = TRUE, labels = FALSE)
  raw <- cbind(df[group_cols], data.frame(bin = b, n = 1L, stringsAsFactors = FALSE))
  out <- stats::aggregate(n ~ ., data = raw, FUN = sum)
  mids <- .fig_midpoints(br)
  out$x_mid <- mids[out$bin]
  out$x_min <- br[out$bin]
  out$x_max <- br[out$bin + 1L]
  rownames(out) <- NULL
  out
}

.fig_symbol_labels <- function(tt, symbol_map = NULL) {
  if ("symbol" %in% names(tt)) {
    lab <- as.character(tt$symbol)
  } else if ("gene" %in% names(tt) && is.data.frame(symbol_map)) {
    .fig_require_cols(symbol_map, c("ensembl", "symbol"), "symbol_map")
    map <- stats::setNames(as.character(symbol_map$symbol), as.character(symbol_map$ensembl))
    lab <- unname(map[as.character(tt$gene)])
    lab[is.na(lab) | lab == ""] <- as.character(tt$gene)[is.na(lab) | lab == ""]
  } else if ("site_id" %in% names(tt)) {
    lab <- as.character(tt$site_id)
  } else if ("feature" %in% names(tt)) {
    lab <- as.character(tt$feature)
  } else {
    lab <- paste0("row", seq_len(nrow(tt)))
  }
  lab[is.na(lab) | lab == ""] <- paste0("row", which(is.na(lab) | lab == ""))
  lab
}

.fig_volcano_data <- function(top_list, contrasts = .fig_contrasts(), symbol_map = NULL,
                              alpha = 0.10, n_label = 10L, x_bins = 70L, y_bins = 55L) {
  stopifnot(is.list(top_list), all(contrasts %in% names(top_list)))
  pieces <- lapply(contrasts, function(cn) {
    tt <- top_list[[cn]]
    .fig_require_cols(tt, c("logFC", "P.Value", "adj.P.Val"), paste0("top$", cn))
    lab <- .fig_symbol_labels(tt, symbol_map)
    x <- as.numeric(tt$logFC)
    p <- pmax(as.numeric(tt$P.Value), 1e-300)
    f <- as.numeric(tt$adj.P.Val)
    ok <- is.finite(x) & is.finite(p) & is.finite(f)
    y <- -log10(p[ok])
    sig <- f[ok] < alpha
    bins <- .fig_grid_bins(x[ok], y, group = rep(cn, sum(ok)),
                           x_bins = x_bins, y_bins = y_bins,
                           symmetric_x = TRUE, lower_y = 0)
    names(bins)[names(bins) == "group"] <- "contrast"
    raw_lab <- data.frame(
      contrast = cn,
      label = lab[ok],
      effect = x[ok],
      neg_log10_p = y,
      neg_log10_fdr = -log10(pmax(f[ok], 1e-300)),
      p_value = p[ok],
      fdr = f[ok],
      significant = sig,
      stringsAsFactors = FALSE
    )
    raw_lab$direction <- ifelse(raw_lab$fdr < alpha & raw_lab$effect > 0, "up",
                                ifelse(raw_lab$fdr < alpha & raw_lab$effect < 0,
                                       "down", "not significant"))
    raw_lab$direction <- factor(raw_lab$direction,
                                levels = c("down", "not significant", "up"))
    ord <- order(raw_lab$fdr, raw_lab$p_value, -abs(raw_lab$effect),
                 raw_lab$label, method = "radix", na.last = TRUE)
    labels <- raw_lab[utils::head(ord, n_label), , drop = FALSE]
    counts <- data.frame(
      contrast = cn,
      n_features = nrow(raw_lab),
      n_fdr_0_10 = sum(raw_lab$fdr < alpha),
      n_up_fdr_0_10 = sum(raw_lab$fdr < alpha & raw_lab$effect > 0),
      n_down_fdr_0_10 = sum(raw_lab$fdr < alpha & raw_lab$effect < 0),
      stringsAsFactors = FALSE
    )
    list(points = raw_lab, bins = bins, labels = labels, counts = counts)
  })
  out <- list(
    points = .fig_bind(lapply(pieces, `[[`, "points")),
    bins = .fig_bind(lapply(pieces, `[[`, "bins")),
    labels = .fig_bind(lapply(pieces, `[[`, "labels")),
    counts = .fig_bind(lapply(pieces, `[[`, "counts"))
  )
  .fig_assert_nonempty(out$points, "volcano points")
  .fig_assert_nonempty(out$bins, "volcano bins")
  .fig_assert_nonempty(out$labels, "volcano labels")
  .fig_assert_finite(out$points, c("effect", "neg_log10_p", "neg_log10_fdr", "p_value", "fdr"),
                     "volcano points")
  .fig_assert_finite(out$bins, c("x_mid", "y_mid", "n"), "volcano bins")
  .fig_assert_finite(out$labels, c("effect", "neg_log10_p", "neg_log10_fdr", "p_value", "fdr"),
                     "volcano labels")
  out
}

.fig_unit_composition <- function(counts) {
  mat <- as.matrix(counts)
  stopifnot(nrow(mat) > 0L, ncol(mat) > 0L, all(is.finite(mat)), all(mat >= 0))
  rn <- rownames(mat) %||% paste0("unit", seq_len(nrow(mat)))
  cn <- colnames(mat) %||% paste0("substate", seq_len(ncol(mat)))
  dimnames(mat) <- list(rn, cn)
  out <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  names(out) <- c("genotype_batch", "substate", "n_cells")
  out$n_cells <- as.numeric(out$n_cells)
  out$genotype <- vapply(out$genotype_batch, function(u) {
    hit <- genotype_levels[startsWith(u, paste0(genotype_levels, "_"))]
    if (length(hit)) hit[1L] else NA_character_
  }, character(1))
  out$batch <- mapply(function(u, g) {
    if (is.na(g)) NA_character_ else substring(u, nchar(g) + 2L)
  }, out$genotype_batch, out$genotype, USE.NAMES = FALSE)
  totals <- stats::aggregate(n_cells ~ genotype_batch, data = out, FUN = sum)
  names(totals)[2L] <- "unit_total"
  out <- merge(out, totals, by = "genotype_batch", all.x = TRUE, sort = FALSE)
  out$proportion <- ifelse(out$unit_total > 0, out$n_cells / out$unit_total, NA_real_)
  out$genotype <- factor(out$genotype, levels = genotype_levels)
  out <- out[order(match(out$genotype, genotype_levels), out$genotype_batch, out$substate,
                   method = "radix", na.last = TRUE), , drop = FALSE]
  rownames(out) <- NULL
  .fig_assert_finite(out, c("n_cells", "unit_total", "proportion"), "unit composition")
  stopifnot(!anyNA(out$genotype))
  out
}

.fig_score_triptych <- function(cf, score_cols) {
  .fig_require_cols(cf, c("umap_1", "umap_2", score_cols), "microglia cell_frame")
  pieces <- lapply(score_cols, function(sc) {
    df <- data.frame(umap_1 = cf$umap_1, umap_2 = cf$umap_2, score = cf[[sc]],
                     score_name = sub("_UCell_z$", "", sc), stringsAsFactors = FALSE)
    ok <- is.finite(df$umap_1) & is.finite(df$umap_2) & is.finite(df$score)
    df <- df[ok, , drop = FALSE]
    xb <- .fig_breaks(df$umap_1, n = 80L, symmetric = FALSE)
    yb <- .fig_breaks(df$umap_2, n = 80L, symmetric = FALSE)
    bx <- cut(.fig_clip(df$umap_1, xb), xb, include.lowest = TRUE, labels = FALSE)
    by <- cut(.fig_clip(df$umap_2, yb), yb, include.lowest = TRUE, labels = FALSE)
    raw <- data.frame(score_name = df$score_name, bin_x = bx, bin_y = by,
                      score_sum = df$score, n = 1L, stringsAsFactors = FALSE)
    ag <- stats::aggregate(cbind(score_sum, n) ~ score_name + bin_x + bin_y,
                           data = raw, FUN = sum)
    ag$score_mean <- ag$score_sum / ag$n
    xm <- .fig_midpoints(xb); ym <- .fig_midpoints(yb)
    ag$umap_1 <- xm[ag$bin_x]
    ag$umap_2 <- ym[ag$bin_y]
    ag[c("score_name", "bin_x", "bin_y", "umap_1", "umap_2", "score_mean", "n")]
  })
  out <- .fig_bind(pieces)
  .fig_assert_finite(out, c("umap_1", "umap_2", "score_mean", "n"), "score triptych")
  out
}

.fig_composition_forest <- function(composition_results) {
  stopifnot(is.list(composition_results))
  if (!is.null(composition_results$sccomp) && nrow(composition_results$sccomp)) {
    x <- composition_results$sccomp
    .fig_require_cols(x, c("contrast", "substate", "c_effect", "c_lower", "c_upper", "c_fdr"),
                      "composition_results$sccomp")
    out <- data.frame(
      contrast = x$contrast,
      substate = x$substate,
      effect = x$c_effect,
      ci_l = x$c_lower,
      ci_r = x$c_upper,
      fdr = x$c_fdr,
      method = "sccomp",
      scale = "logit proportion",
      stringsAsFactors = FALSE
    )
  } else {
    x <- composition_results$propeller_logit
    .fig_require_cols(x, c("contrast", "substate", "t", "fdr_global"),
                      "composition_results$propeller_logit")
    out <- data.frame(
      contrast = x$contrast,
      substate = x$substate,
      effect = x$t,
      ci_l = x$t,
      ci_r = x$t,
      fdr = x$fdr_global,
      method = "propeller_logit",
      scale = "t-statistic",
      stringsAsFactors = FALSE
    )
  }
  .fig_assert_finite(out, c("effect", "ci_l", "ci_r", "fdr"), "composition_forest")
  out
}

microglia_figure_data <- function(microglia_report, composition_results,
                                  pb_de_microglia, pb_de_substate, symbol_map,
                                  alpha = 0.10) {
  stopifnot(is.list(microglia_report), is.list(composition_results),
            is.list(pb_de_microglia), is.list(pb_de_substate), is.data.frame(symbol_map))
  cf <- microglia_report$cell_frame
  score_cols <- c("Homeostatic_UCell_z", "DAM_UCell_z", "MHC_APC_UCell_z")
  .fig_require_cols(cf, c("umap_1", "umap_2", "genotype", "substate", score_cols),
                    "microglia_report$cell_frame")
  umap <- cf[c("umap_1", "umap_2", "genotype", "substate")]
  score_long <- .fig_bind(lapply(score_cols, function(sc) {
    data.frame(genotype = cf$genotype, substate = cf$substate,
               score_name = sub("_UCell_z$", "", sc), score = cf[[sc]],
               stringsAsFactors = FALSE)
  }))
  score_dist <- .fig_histogram(score_long, "score", c("genotype", "substate", "score_name"),
                               bins = 50L)

  comp <- .fig_unit_composition(composition_results$counts)
  genotype_comp <- stats::aggregate(n_cells ~ genotype + substate, data = comp, FUN = sum)
  genotype_total <- stats::aggregate(n_cells ~ genotype, data = genotype_comp, FUN = sum)
  names(genotype_total)[2L] <- "genotype_total"
  genotype_comp <- merge(genotype_comp, genotype_total, by = "genotype", all.x = TRUE, sort = FALSE)
  genotype_comp$proportion <- genotype_comp$n_cells / pmax(genotype_comp$genotype_total, 1)
  .fig_require_cols(composition_results$concordance,
                    c("contrast", "substate", "dir_logit", "sig_logit",
                      "dir_asin", "sig_asin", "dir_concordant", "sig_concordant", "flag"),
                    "composition_results$concordance")
  concordance <- composition_results$concordance
  composition_forest <- .fig_composition_forest(composition_results)

  whole_volcano <- .fig_volcano_data(pb_de_microglia$top, .fig_contrasts(), symbol_map,
                                     alpha = alpha, n_label = 8L)
  amyloid_volcano <- lapply(whole_volcano, function(x) {
    if (is.data.frame(x) && "contrast" %in% names(x)) {
      x[x$contrast == "nlgf_in_p301s", , drop = FALSE]
    } else {
      x
    }
  })
  substate_volcano_pieces <- lapply(names(pb_de_substate$per_substate), function(st) {
    res <- pb_de_substate$per_substate[[st]]
    if (!is.list(res) || !identical(res$status, "fit") || !"top" %in% names(res)) return(NULL)
    v <- .fig_volcano_data(res$top, .fig_contrasts(), symbol_map,
                           alpha = alpha, n_label = 5L)
    v$bins$substate <- st
    v$labels$substate <- st
    v$counts$substate <- st
    v
  })
  substate_volcano_pieces <- Filter(Negate(is.null), substate_volcano_pieces)
  substate_volcano <- list(
    bins = .fig_bind(lapply(substate_volcano_pieces, `[[`, "bins")),
    labels = .fig_bind(lapply(substate_volcano_pieces, `[[`, "labels")),
    counts = .fig_bind(lapply(substate_volcano_pieces, `[[`, "counts"))
  )

  cc <- as.matrix(pb_de_substate$cell_counts)
  cell_counts <- as.data.frame(as.table(cc), stringsAsFactors = FALSE)
  names(cell_counts) <- c("substate", "genotype_batch", "n_cells")
  cell_counts$n_cells <- as.numeric(cell_counts$n_cells)
  fit_status <- .fig_bind(lapply(names(pb_de_substate$per_substate), function(st) {
    res <- pb_de_substate$per_substate[[st]]
    data.frame(
      substate = st,
      status = as.character(res$status %||% NA_character_),
      n_cells = as.numeric(res$n_cells %||% NA_real_),
      units = as.numeric(res$units %||% NA_real_),
      reason = as.character(res$reason %||% NA_character_),
      stringsAsFactors = FALSE
    )
  }))

  out <- list(
    manifest = figure_manifest("microglia"),
    umap_by_substate = umap,
    score_triptych = .fig_score_triptych(cf, score_cols),
    unit_composition = comp,
    score_distribution = score_dist,
    composition_shift = list(genotype_substate = genotype_comp,
                             score_distribution = score_dist),
    composition_forest = composition_forest,
    amyloid_volcano = amyloid_volcano,
    whole_de_volcano = whole_volcano,
    substate_fit_audit = list(cell_counts = cell_counts, fit_status = fit_status),
    substate_de_volcano = substate_volcano,
    provenance = list(
      source_targets = c("microglia_report", "composition_results",
                         "pb_de_microglia", "pb_de_substate", "symbol_map"),
      alpha = alpha,
      contract = "compact inline microglia figure data; no heavy Seurat object"
    )
  )

  .fig_assert_finite(out$umap_by_substate, c("umap_1", "umap_2"), "umap_by_substate")
  .fig_assert_nonempty(out$score_triptych, "score_triptych")
  .fig_assert_nonempty(out$score_distribution, "score_distribution")
  .fig_assert_nonempty(out$composition_forest, "composition_forest")
  .fig_assert_nonempty(out$substate_de_volcano$bins, "substate_de_volcano bins")
  out
}

trajectory_figure_data <- function(trajectory_report, composition_results, alpha = 0.10) {
  stopifnot(is.list(trajectory_report), is.list(composition_results))
  cf <- trajectory_report$cell_frame
  .fig_require_cols(cf, c("genotype", "substate", "on_lineage", "pt_raw", "score_axis_pt"),
                    "trajectory cell_frame")
  pt <- cf[cf$on_lineage %in% TRUE & is.finite(cf$pt_raw), , drop = FALSE]
  density <- .fig_histogram(pt, "pt_raw", c("genotype", "substate"), bins = 55L, lower = 0)
  concordance <- .fig_grid_bins(pt$pt_raw, pt$score_axis_pt, group = rep("all", nrow(pt)),
                                x_bins = 60L, y_bins = 60L,
                                symmetric_x = FALSE, lower_y = NULL)
  names(concordance)[names(concordance) == "group"] <- "panel"

  comp <- .fig_unit_composition(composition_results$counts)
  dam <- comp[comp$substate == "DAM", c("genotype_batch", "proportion"), drop = FALSE]
  names(dam)[2L] <- "dam_fraction"
  unit <- merge(trajectory_report$per_unit, dam, by = "genotype_batch", all.x = TRUE, sort = FALSE)
  .fig_require_cols(unit, c("genotype_batch", "genotype", "batch", "mean_pt",
                            "within_homeostatic", "within_dam", "dam_fraction"),
                    "trajectory unit_pt_vs_dam")

  measures <- c("mean_pt", "comp_cf", "progression_cf", "cross",
                "within_homeostatic", "within_dam")
  forest <- .fig_bind(lapply(names(trajectory_report$weighted_top), function(cn) {
    w <- trajectory_report$weighted_top[[cn]]
    .fig_require_cols(w, c("measure", "coef", "ci_l", "ci_r", "p_value"), paste0("weighted_top$", cn))
    w <- w[w$measure %in% measures, , drop = FALSE]
    data.frame(contrast = cn, w, stringsAsFactors = FALSE)
  }))
  forest$significant <- is.finite(forest$p_value) & forest$p_value < alpha
  mean_position <- forest[forest$measure == "mean_pt", , drop = FALSE]

  interaction <- trajectory_report$interaction
  .fig_require_cols(interaction, c("measure", "coef", "ci_l", "ci_r", "p_value", "fdr", "perm_p"),
                    "trajectory_report$interaction")
  decomp_measures <- c("mean_pt", "comp_cf", "progression_cf", "cross",
                       "within_homeostatic", "within_dam")
  decomposition <- interaction[interaction$measure %in% decomp_measures, , drop = FALSE]
  decomposition$significant <- is.finite(decomposition$fdr) & decomposition$fdr < alpha

  audit <- list(
    sensitivity = trajectory_report$sensitivity,
    lineage_per_unit = trajectory_report$lineage_per_unit,
    glmm = as.data.frame(trajectory_report$glmm, stringsAsFactors = FALSE),
    provenance = as.data.frame(trajectory_report$provenance[c("concordance_rho",
                                                              "omitted_frac_overall",
                                                              "composition_loading",
                                                              "progression_loading",
                                                              "cross_loading")],
                               stringsAsFactors = FALSE)
  )

  out <- list(
    manifest = figure_manifest("trajectory"),
    pt_density = density,
    pseudotime_shift = list(density = density, mean_position = mean_position),
    unit_pt_vs_dam = unit,
    decomposition = decomposition,
    kitagawa_forest = forest,
    concordance = concordance,
    trajectory_audit = audit,
    provenance = list(
      source_targets = c("trajectory_report", "composition_results"),
      alpha = alpha,
      contract = "compact trajectory figure data; composition joined only from compact counts"
    )
  )
  .fig_assert_finite(out$pt_density, c("x_mid", "n"), "pt_density")
  .fig_assert_finite(out$unit_pt_vs_dam, c("mean_pt", "dam_fraction"), "unit_pt_vs_dam")
  .fig_assert_finite(out$kitagawa_forest, c("coef", "ci_l", "ci_r", "p_value"), "kitagawa_forest")
  .fig_assert_finite(out$decomposition, c("coef", "ci_l", "ci_r", "p_value", "fdr"), "decomposition")
  .fig_assert_finite(out$concordance, c("x_mid", "y_mid", "n"), "trajectory concordance")
  out
}

mechanism_figure_data <- function(mechanism_report, alpha = 0.10) {
  stopifnot(is.list(mechanism_report))
  project <- mechanism_report$pathway_project
  go_top <- mechanism_report$pathway_go_top
  tf <- mechanism_report$tf_highlights
  kinase <- mechanism_report$kinase$table
  .fig_require_cols(project, c("pathway", "population", "contrast", "NES", "fdr"),
                    "mechanism pathway_project")
  .fig_require_cols(go_top, c("pathway", "population", "contrast", "NES", "fdr", "size"),
                    "mechanism pathway_go_top")
  .fig_require_cols(tf, c("population", "source", "contrast", "score", "fdr", "selection"),
                    "mechanism tf_highlights")
  .fig_require_cols(kinase, c("source", "contrast", "score", "fdr", "significant",
                              "run_index_score", "run_index_fdr", "run_index_supports",
                              "include_reason"),
                    "mechanism kinase")
  out <- list(
    manifest = figure_manifest("mechanism"),
    project_pathway_heatmap = project,
    go_top_dotplot = go_top,
    tf_interaction = tf[tf$contrast == "interaction", , drop = FALSE],
    tf_focus = tf,
    kinase_heatmap = kinase,
    provenance = list(
      source_targets = "mechanism_report",
      alpha = alpha,
      contract = "compact mechanism figure data from the guarded report bundle"
    )
  )
  .fig_assert_finite(out$project_pathway_heatmap, c("NES", "fdr"), "project_pathway_heatmap")
  .fig_assert_finite(out$go_top_dotplot, c("NES", "fdr", "size"), "go_top_dotplot")
  .fig_assert_finite(out$tf_focus, c("score", "fdr"), "tf_focus")
  .fig_assert_finite(out$kinase_heatmap, c("score", "fdr", "run_index_score"), "kinase_heatmap")
  out
}

.fig_sensitivity_counts <- function(primary_top, sensitivity, contrasts, alpha = 0.10) {
  stopifnot(is.list(primary_top), is.list(sensitivity), all(contrasts %in% names(primary_top)))
  fits <- names(sensitivity)
  .fig_bind(lapply(fits, function(fit) {
    sens <- sensitivity[[fit]]
    .fig_bind(lapply(contrasts, function(cn) {
      p <- primary_top[[cn]]
      .fig_require_cols(p, c("symbol", "logFC", "adj.P.Val"), paste0("primary$", cn))
      if (!is.list(sens) || !identical(sens$status %||% NA_character_, "fit") ||
          !"top" %in% names(sens) || !cn %in% names(sens$top)) {
        return(data.frame(fit = fit, contrast = cn, status = as.character(sens$status %||% NA_character_),
                          n_primary_sig = sum(p$adj.P.Val < alpha, na.rm = TRUE),
                          n_sensitivity_sig = NA_real_, n_lost = NA_real_,
                          n_gained = NA_real_, n_sign_flip = NA_real_,
                          stringsAsFactors = FALSE))
      }
      s <- sens$top[[cn]]
      .fig_require_cols(s, c("symbol", "logFC", "adj.P.Val"), paste0(fit, "$", cn))
      m <- merge(p[c("symbol", "logFC", "adj.P.Val")],
                 s[c("symbol", "logFC", "adj.P.Val")],
                 by = "symbol", suffixes = c("_primary", "_sensitivity"))
      ps <- is.finite(m$adj.P.Val_primary) & m$adj.P.Val_primary < alpha
      ss <- is.finite(m$adj.P.Val_sensitivity) & m$adj.P.Val_sensitivity < alpha
      data.frame(
        fit = fit, contrast = cn, status = as.character(sens$status %||% "fit"),
        n_primary_sig = sum(ps),
        n_sensitivity_sig = sum(ss),
        n_lost = sum(ps & !ss),
        n_gained = sum(!ps & ss),
        n_sign_flip = sum(ps & ss & sign(m$logFC_primary) != sign(m$logFC_sensitivity)),
        stringsAsFactors = FALSE
      )
    }))
  }))
}

.fig_best_by_id <- function(x, id_col, effect_col = "logFC", fdr_col = "adj.P.Val") {
  .fig_require_cols(x, c(id_col, effect_col, fdr_col), "best_by_id input")
  x <- x[!is.na(x[[id_col]]) & x[[id_col]] != "" &
           is.finite(x[[effect_col]]) & is.finite(x[[fdr_col]]), , drop = FALSE]
  x <- x[order(x[[id_col]], x[[fdr_col]], -abs(x[[effect_col]]),
               method = "radix", na.last = TRUE), , drop = FALSE]
  x[!duplicated(x[[id_col]]), , drop = FALSE]
}

.fig_phospho_correction <- function(raw, corrected, contrasts, alpha = 0.10) {
  stopifnot(is.list(raw), is.list(corrected), all(contrasts %in% names(raw)),
            all(contrasts %in% names(corrected)))
  pieces <- lapply(contrasts, function(cn) {
    r <- .fig_best_by_id(raw[[cn]], "site_id")
    c <- .fig_best_by_id(corrected[[cn]], "site_id")
    m <- merge(r[c("site_id", "logFC", "adj.P.Val")],
               c[c("site_id", "logFC", "adj.P.Val")],
               by = "site_id", suffixes = c("_raw", "_corrected"))
    sig <- (m$adj.P.Val_raw < alpha) | (m$adj.P.Val_corrected < alpha)
    bins <- .fig_grid_bins(m$logFC_raw, m$logFC_corrected, group = rep(cn, nrow(m)),
                           x_bins = 70L, y_bins = 70L, symmetric_x = TRUE, lower_y = NULL)
    names(bins)[names(bins) == "group"] <- "contrast"
    ord <- order(pmin(m$adj.P.Val_raw, m$adj.P.Val_corrected),
                 -abs(m$logFC_raw) - abs(m$logFC_corrected),
                 m$site_id, method = "radix", na.last = TRUE)
    labels <- m[utils::head(ord, 12L), , drop = FALSE]
    points <- data.frame(
      contrast = cn,
      site_id = m$site_id,
      raw_logFC = m$logFC_raw,
      corrected_logFC = m$logFC_corrected,
      raw_fdr = m$adj.P.Val_raw,
      corrected_fdr = m$adj.P.Val_corrected,
      status = ifelse(m$adj.P.Val_raw < alpha & m$adj.P.Val_corrected < alpha, "both",
                      ifelse(m$adj.P.Val_raw < alpha, "raw only",
                             ifelse(m$adj.P.Val_corrected < alpha, "corrected only",
                                    "not significant"))),
      stringsAsFactors = FALSE
    )
    points$status <- factor(points$status,
                            levels = c("not significant", "raw only",
                                       "corrected only", "both"))
    data <- data.frame(
      contrast = cn,
      site_id = labels$site_id,
      raw_logFC = labels$logFC_raw,
      corrected_logFC = labels$logFC_corrected,
      raw_fdr = labels$adj.P.Val_raw,
      corrected_fdr = labels$adj.P.Val_corrected,
      stringsAsFactors = FALSE
    )
    counts <- data.frame(
      contrast = cn, n_sites = nrow(m), n_any_fdr_0_10 = sum(sig),
      n_raw_only = sum(m$adj.P.Val_raw < alpha & !(m$adj.P.Val_corrected < alpha)),
      n_corrected_only = sum(!(m$adj.P.Val_raw < alpha) & m$adj.P.Val_corrected < alpha),
      n_both = sum(m$adj.P.Val_raw < alpha & m$adj.P.Val_corrected < alpha),
      stringsAsFactors = FALSE
    )
    list(points = points, bins = bins, labels = data, counts = counts)
  })
  out <- list(
    points = .fig_bind(lapply(pieces, `[[`, "points")),
    bins = .fig_bind(lapply(pieces, `[[`, "bins")),
    labels = .fig_bind(lapply(pieces, `[[`, "labels")),
    counts = .fig_bind(lapply(pieces, `[[`, "counts"))
  )
  .fig_assert_finite(out$points, c("raw_logFC", "corrected_logFC", "raw_fdr", "corrected_fdr"),
                     "phospho correction points")
  .fig_assert_finite(out$bins, c("x_mid", "y_mid", "n"), "phospho correction bins")
  .fig_assert_finite(out$labels, c("raw_logFC", "corrected_logFC", "raw_fdr", "corrected_fdr"),
                     "phospho correction labels")
  out
}

.fig_four_modality_classes <- function() {
  c(
    snRNAseq_microglia = "snRNAseq microglia",
    GeoMx_spatial = "GeoMx spatial",
    bulk_proteome = "bulk proteome",
    bulk_phosphoproteome = "bulk phosphoproteome"
  )
}

.fig_semicolon_tokens <- function(x) {
  if (length(x) != 1L || is.na(x) || !nzchar(x)) return(character())
  z <- unlist(strsplit(as.character(x), ";", fixed = TRUE), use.names = FALSE)
  z[nzchar(z)]
}

.fig_best_four_modality_evidence <- function(crossmodality_table, contrasts,
                                             alpha = 0.10,
                                             modalities = names(.fig_four_modality_classes())) {
  stopifnot(is.list(crossmodality_table), is.data.frame(crossmodality_table$evidence),
            is.numeric(alpha), length(alpha) == 1L, alpha > 0, alpha < 1)
  ev <- crossmodality_table$evidence
  .fig_require_cols(ev, c("symbol", "contrast", "modality_class", "modality_group",
                          "effect", "statistic", "fdr", "sign", "missingness_reason"),
                    "crossmodality_table$evidence")
  ev <- ev[ev$contrast %in% contrasts &
             ev$modality_class %in% modalities &
             !is.na(ev$symbol) & ev$symbol != "" &
             is.na(ev$missingness_reason), , drop = FALSE]
  stopifnot(nrow(ev) > 0L, all(contrasts %in% ev$contrast),
            all(modalities %in% ev$modality_class))
  fdr_ord <- ev$fdr
  fdr_ord[!is.finite(fdr_ord)] <- Inf
  stat_abs <- abs(ev$statistic)
  stat_abs[!is.finite(stat_abs)] <- -Inf
  ev <- ev[order(ev$contrast, ev$modality_class, ev$symbol, fdr_ord, -stat_abs,
                 ev$modality_group, method = "radix", na.last = TRUE), , drop = FALSE]
  key <- interaction(ev$contrast, ev$modality_class, ev$symbol, drop = TRUE, lex.order = TRUE)
  out <- ev[!duplicated(key), , drop = FALSE]
  out$significant <- is.finite(out$fdr) & out$fdr < alpha
  out$best_sign <- ifelse(is.finite(out$sign) & out$sign != 0L, as.integer(out$sign),
                          ifelse(out$effect > 0, 1L, ifelse(out$effect < 0, -1L, 0L)))
  out[, c("symbol", "contrast", "modality_class", "modality_group", "effect",
          "statistic", "fdr", "significant", "best_sign"), drop = FALSE]
}

.fig_four_modality_counts <- function(crossmodality_table, contrasts, alpha = 0.10) {
  modalities <- .fig_four_modality_classes()
  best <- .fig_best_four_modality_evidence(crossmodality_table, contrasts, alpha,
                                           names(modalities))
  split_key <- interaction(best$contrast, best$modality_class, drop = TRUE, lex.order = TRUE)
  rows <- lapply(split(best, split_key), function(d) {
    data.frame(
      contrast = d$contrast[1L],
      modality_class = d$modality_class[1L],
      n_symbols = nrow(d),
      n_sig = sum(d$significant),
      n_up_sig = sum(d$significant & d$best_sign > 0L),
      n_down_sig = sum(d$significant & d$best_sign < 0L),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  grid <- expand.grid(contrast = contrasts, modality_class = names(modalities),
                      stringsAsFactors = FALSE)
  out <- merge(grid, out, by = c("contrast", "modality_class"), all.x = TRUE,
               sort = FALSE)
  count_cols <- c("n_symbols", "n_sig", "n_up_sig", "n_down_sig")
  out[count_cols] <- lapply(out[count_cols], function(x) {
    x[is.na(x)] <- 0L
    as.integer(x)
  })
  out$modality <- unname(modalities[out$modality_class])
  out$signed_balance <- ifelse(out$n_sig > 0L,
                               (out$n_up_sig - out$n_down_sig) / out$n_sig,
                               0)
  out$log_n_sig <- log10(out$n_sig + 1)
  out[order(match(out$contrast, contrasts), match(out$modality_class, names(modalities)),
            method = "radix"), , drop = FALSE]
}

.fig_four_modality_pathways <- function(crossmodality_report, contrasts) {
  stopifnot(is.list(crossmodality_report), is.list(crossmodality_report$pathway),
            is.data.frame(crossmodality_report$pathway$axis_summary))
  x <- crossmodality_report$pathway$axis_summary
  .fig_require_cols(x, c("axis", "contrast", "n_modalities_present", "n_modalities_sig",
                         "n_evidence_groups_sig", "n_positive_modalities",
                         "n_negative_modalities", "mixed_sign", "consistent_direction",
                         "rank_score"),
                    "crossmodality_report$pathway$axis_summary")
  x <- x[x$contrast %in% contrasts, , drop = FALSE]
  stopifnot(nrow(x) > 0L, all(contrasts %in% x$contrast))
  x$direction <- ifelse(x$mixed_sign, "mixed", as.character(x$consistent_direction))
  x$direction[is.na(x$direction) | !nzchar(x$direction)] <- "none"
  x$direction <- factor(x$direction, levels = c("positive", "negative", "mixed", "none"))
  x$modality_fraction <- x$n_modalities_sig / pmax(1, x$n_modalities_present)
  x[order(match(x$contrast, contrasts), -x$n_modalities_sig, -x$n_modalities_present,
          -x$rank_score, x$axis, method = "radix", na.last = TRUE), , drop = FALSE]
}

.fig_four_modality_symbols <- function(crossmodality_report, contrasts,
                                       top_per_contrast = 6L) {
  stopifnot(is.list(crossmodality_report), is.list(crossmodality_report$divergence),
            is.data.frame(crossmodality_report$divergence$axis_symbols),
            top_per_contrast >= 1L)
  modalities <- .fig_four_modality_classes()
  x <- crossmodality_report$divergence$axis_symbols
  .fig_require_cols(x, c("symbol", "contrast", "axis", "n_modalities_present",
                         "n_modalities_sig", "min_fdr", "max_abs_statistic",
                         "modalities_present", "modalities_sig", "direction_call",
                         "rank_score"),
                    "crossmodality_report$divergence$axis_symbols")
  x <- x[x$contrast %in% contrasts, , drop = FALSE]
  stopifnot(nrow(x) > 0L, all(contrasts %in% x$contrast))
  x$n_four_modalities_present <- vapply(x$modalities_present, function(z) {
    sum(names(modalities) %in% .fig_semicolon_tokens(z))
  }, integer(1), USE.NAMES = FALSE)
  x$n_four_modalities_sig <- vapply(x$modalities_sig, function(z) {
    sum(names(modalities) %in% .fig_semicolon_tokens(z))
  }, integer(1), USE.NAMES = FALSE)
  pieces <- lapply(contrasts, function(cn) {
    d <- x[x$contrast == cn, , drop = FALSE]
    d <- d[order(-d$n_four_modalities_present, -d$n_four_modalities_sig, -d$rank_score,
                 d$min_fdr, -d$max_abs_statistic, d$symbol,
                 method = "radix", na.last = TRUE), , drop = FALSE]
    utils::head(d, top_per_contrast)
  })
  picked <- do.call(rbind, pieces)
  rows <- lapply(seq_len(nrow(picked)), function(i) {
    d <- picked[i, , drop = FALSE]
    present <- .fig_semicolon_tokens(d$modalities_present)
    sig <- .fig_semicolon_tokens(d$modalities_sig)
    data.frame(
      symbol = d$symbol,
      contrast = d$contrast,
      axis = d$axis,
      direction_call = d$direction_call,
      min_fdr = d$min_fdr,
      n_modalities_present = d$n_four_modalities_present,
      n_modalities_sig = d$n_four_modalities_sig,
      modality_class = names(modalities),
      modality = unname(modalities),
      status = ifelse(names(modalities) %in% sig, "FDR < 0.10",
                      ifelse(names(modalities) %in% present, "measured", "not observed")),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$status <- factor(out$status, levels = c("not observed", "measured", "FDR < 0.10"))
  out$symbol_axis <- ifelse(is.na(out$axis) | out$axis == "other" | out$axis == "",
                            out$symbol,
                            paste0(out$symbol, " (", out$axis, ")"))
  out[order(match(out$contrast, contrasts), -out$n_modalities_present,
            -out$n_modalities_sig, out$symbol, match(out$modality_class, names(modalities)),
            method = "radix", na.last = TRUE), , drop = FALSE]
}

.fig_axis_effect_axes <- function() {
  c("DAM", "antigen_presentation", "synaptic", "clearance",
    "interaction_boundary", "mechanism_boundary")
}

.fig_axis_effect_contrasts <- function() {
  c("nlgf_in_maptki", "nlgf_in_p301s", "tau_in_nlgf", "interaction")
}

.fig_axis_effect_modalities <- function() {
  c(
    .fig_four_modality_classes(),
    TF_activity = "TF activity",
    kinase_activity = "kinase activity",
    crossmodality_pathway = "cross-modality pathway",
    clearance_pair = "clearance pair"
  )
}

.fig_axis_effect_standard_axis <- function(axis) {
  out <- rep(NA_character_, length(axis))
  z <- as.character(axis)
  out[z %in% "DAM" | grepl("DAM", z, fixed = TRUE)] <- "DAM"
  out[z %in% c("antigen_presentation", "MHC_APC") |
        grepl("MHC_APC", z, fixed = TRUE)] <- "antigen_presentation"
  out[z %in% "synaptic" | grepl("synaptic", z, fixed = TRUE)] <- "synaptic"
  out[z %in% c("clearance", "complement") |
        grepl("clearance", z, fixed = TRUE) |
        grepl("complement", z, fixed = TRUE)] <- "clearance"
  out[z %in% c("NFkB", "gsk3b")] <- "mechanism_boundary"
  out
}

.fig_axis_effect_anchor_symbols <- function() {
  rows <- data.frame(
    axis = c(rep("DAM", 3), "antigen_presentation",
             rep("synaptic", 3), rep("clearance", 4),
             rep("mechanism_boundary", 4)),
    symbol = c("Apoe", "Cst7", "Spp1", "Cd74",
               "Syn1", "Syp", "Snap25",
               "Apoe", "Trem2", "Cd74", "Spp1",
               "Myc", "Nfkb1", "Rela", "Gsk3b"),
    stringsAsFactors = FALSE
  )
  rows$feature_id <- paste(rows$axis, rows$symbol, sep = ":")
  rows$feature_label <- rows$symbol
  rows$feature_type <- "symbol"
  rows$selection_rule <- "predeclared_anchor_symbol"
  rows$selection_rank <- seq_len(nrow(rows))
  rows
}

.fig_axis_effect_top_symbols <- function(crossmodality_report, top_per_axis = 2L) {
  stopifnot(is.list(crossmodality_report), is.list(crossmodality_report$divergence),
            is.data.frame(crossmodality_report$divergence$axis_symbols),
            top_per_axis >= 0L)
  if (top_per_axis == 0L) return(.fig_axis_effect_anchor_symbols()[0, ])
  x <- crossmodality_report$divergence$axis_symbols
  .fig_require_cols(x, c("symbol", "axis", "n_modalities_present",
                         "n_modalities_sig", "min_fdr", "max_abs_statistic",
                         "rank_score"),
                    "crossmodality_report$divergence$axis_symbols")
  x$axis <- .fig_axis_effect_standard_axis(x$axis)
  x <- x[!is.na(x$axis) & x$axis %in% .fig_axis_effect_axes() &
           !is.na(x$symbol) & nzchar(x$symbol), , drop = FALSE]
  if (!nrow(x)) return(.fig_axis_effect_anchor_symbols()[0, ])
  split_key <- interaction(x$axis, x$symbol, drop = TRUE, lex.order = TRUE)
  agg <- do.call(rbind, lapply(split(x, split_key), function(d) {
    data.frame(
      axis = d$axis[1L],
      symbol = d$symbol[1L],
      n_modalities_present = max(d$n_modalities_present, na.rm = TRUE),
      n_modalities_sig = max(d$n_modalities_sig, na.rm = TRUE),
      min_fdr = suppressWarnings(min(d$min_fdr, na.rm = TRUE)),
      max_abs_statistic = max(d$max_abs_statistic, na.rm = TRUE),
      rank_score = max(d$rank_score, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  agg$min_fdr[!is.finite(agg$min_fdr)] <- Inf
  anchors <- .fig_axis_effect_anchor_symbols()
  anchor_key <- paste(anchors$axis, anchors$symbol, sep = "\r")
  agg <- agg[!(paste(agg$axis, agg$symbol, sep = "\r") %in% anchor_key), , drop = FALSE]
  if (!nrow(agg)) return(anchors[0, ])
  agg <- agg[order(match(agg$axis, .fig_axis_effect_axes()),
                   -agg$n_modalities_present, -agg$n_modalities_sig,
                   agg$min_fdr, -agg$max_abs_statistic, -agg$rank_score,
                   agg$symbol, method = "radix", na.last = TRUE), , drop = FALSE]
  pieces <- lapply(.fig_axis_effect_axes(), function(ax) {
    utils::head(agg[agg$axis == ax, , drop = FALSE], top_per_axis)
  })
  out <- do.call(rbind, pieces)
  if (!nrow(out)) return(anchors[0, ])
  out$feature_id <- paste(out$axis, out$symbol, sep = ":")
  out$feature_label <- out$symbol
  out$feature_type <- "symbol"
  out$selection_rule <- "top_ranked_axis_symbol"
  out$selection_rank <- 100L + seq_len(nrow(out))
  out[, c("axis", "symbol", "feature_id", "feature_label", "feature_type",
          "selection_rule", "selection_rank"), drop = FALSE]
}

.fig_axis_effect_selection <- function(crossmodality_report) {
  modalities <- .fig_axis_effect_modalities()
  assay_modalities <- names(.fig_four_modality_classes())
  symbols <- rbind(.fig_axis_effect_anchor_symbols(),
                   .fig_axis_effect_top_symbols(crossmodality_report))
  symbol_rows <- do.call(rbind, lapply(seq_len(nrow(symbols)), function(i) {
    d <- symbols[i, , drop = FALSE]
    data.frame(
      axis = d$axis,
      feature_id = d$feature_id,
      feature_label = d$feature_label,
      feature_type = d$feature_type,
      modality_class = assay_modalities,
      modality = unname(modalities[assay_modalities]),
      lookup_symbol = d$symbol,
      lookup_source = "evidence_symbol",
      applicable = TRUE,
      selection_rule = d$selection_rule,
      selection_rank = d$selection_rank,
      stringsAsFactors = FALSE
    )
  }))

  pathway_rows <- data.frame(
    axis = .fig_axis_effect_axes(),
    feature_id = paste0("pathway:", .fig_axis_effect_axes()),
    feature_label = paste(gsub("_", " ", .fig_axis_effect_axes(), fixed = TRUE),
                          "axis"),
    feature_type = "pathway",
    modality_class = "crossmodality_pathway",
    modality = unname(modalities["crossmodality_pathway"]),
    lookup_symbol = NA_character_,
    lookup_source = "pathway_axis",
    applicable = TRUE,
    selection_rule = "fixed_axis_pathway_summary",
    selection_rank = 300L + seq_along(.fig_axis_effect_axes()),
    stringsAsFactors = FALSE
  )

  pairs <- data.frame(
    axis = "clearance",
    feature_id = paste0("pair:", c("Apoe_Trem2", "App_Cd74", "Pros1_Mertk")),
    feature_label = gsub("_", "-", c("Apoe_Trem2", "App_Cd74", "Pros1_Mertk"),
                         fixed = TRUE),
    feature_type = "pair",
    modality_class = "clearance_pair",
    modality = unname(modalities["clearance_pair"]),
    lookup_symbol = c("Apoe_Trem2", "App_Cd74", "Pros1_Mertk"),
    lookup_source = "clearance_pair",
    applicable = TRUE,
    selection_rule = "predeclared_clearance_pair",
    selection_rank = 400L + seq_len(3L),
    stringsAsFactors = FALSE
  )

  mech <- data.frame(
    axis = "mechanism_boundary",
    feature_id = c("boundary:Myc_TF", "boundary:Nfkb1_TF",
                   "boundary:Rela_TF", "boundary:Gsk3b_kinase"),
    feature_label = c("Myc TF", "Nfkb1 TF", "Rela TF", "Gsk3b kinase"),
    feature_type = "boundary",
    symbol = c("Myc", "Nfkb1", "Rela", "Gsk3b"),
    applicable_modality = c("TF_activity", "TF_activity", "TF_activity",
                            "kinase_activity"),
    selection_rank = 500L + seq_len(4L),
    stringsAsFactors = FALSE
  )
  mech_rows <- do.call(rbind, lapply(seq_len(nrow(mech)), function(i) {
    d <- mech[i, , drop = FALSE]
    data.frame(
      axis = d$axis,
      feature_id = d$feature_id,
      feature_label = d$feature_label,
      feature_type = d$feature_type,
      modality_class = c("TF_activity", "kinase_activity"),
      modality = unname(modalities[c("TF_activity", "kinase_activity")]),
      lookup_symbol = d$symbol,
      lookup_source = "mechanism_boundary",
      applicable = c("TF_activity", "kinase_activity") == d$applicable_modality,
      selection_rule = "predeclared_mechanism_boundary",
      selection_rank = d$selection_rank,
      stringsAsFactors = FALSE
    )
  }))

  spatial_rows <- data.frame(
    axis = "interaction_boundary",
    feature_id = "boundary:SpatialDecon_abundance",
    feature_label = "SpatialDecon abundance",
    feature_type = "boundary",
    modality_class = assay_modalities,
    modality = unname(modalities[assay_modalities]),
    lookup_symbol = NA_character_,
    lookup_source = "spatial_decon_boundary",
    applicable = assay_modalities == "GeoMx_spatial",
    selection_rule = "blocked_spatial_abundance_boundary",
    selection_rank = 600L,
    stringsAsFactors = FALSE
  )

  out <- rbind(symbol_rows, pathway_rows, pairs, mech_rows, spatial_rows)
  out$selection_key <- paste(out$axis, out$feature_id, out$modality_class,
                             sep = "\r")
  stopifnot(!anyDuplicated(out$selection_key),
            all(out$axis %in% .fig_axis_effect_axes()),
            all(out$modality_class %in% names(modalities)),
            all(out$feature_type %in% c("symbol", "pathway", "pair", "boundary")))
  out[order(match(out$axis, .fig_axis_effect_axes()), out$selection_rank,
            out$feature_id, match(out$modality_class, names(modalities)),
            method = "radix"), , drop = FALSE]
}

.fig_axis_effect_best_evidence <- function(crossmodality_table, contrasts, alpha = 0.10) {
  stopifnot(is.list(crossmodality_table), is.data.frame(crossmodality_table$evidence),
            is.numeric(alpha), length(alpha) == 1L, alpha > 0, alpha < 1)
  ev <- crossmodality_table$evidence
  .fig_require_cols(ev, c("symbol", "contrast", "modality_class", "modality_group",
                          "feature_type", "effect", "statistic", "fdr", "sign",
                          "feature_id", "missingness_reason"),
                    "crossmodality_table$evidence")
  ev <- ev[ev$contrast %in% contrasts &
             ev$modality_class %in% names(.fig_axis_effect_modalities()) &
             !is.na(ev$symbol) & nzchar(ev$symbol) &
             is.na(ev$missingness_reason), , drop = FALSE]
  stopifnot(nrow(ev) > 0L)
  fdr_ord <- ev$fdr
  fdr_ord[!is.finite(fdr_ord)] <- Inf
  stat_abs <- abs(ev$statistic)
  stat_abs[!is.finite(stat_abs)] <- -Inf
  ev <- ev[order(ev$symbol, ev$contrast, ev$modality_class,
                 fdr_ord, -stat_abs, ev$modality_group, ev$feature_id,
                 method = "radix", na.last = TRUE), , drop = FALSE]
  key <- interaction(ev$symbol, ev$contrast, ev$modality_class,
                     drop = TRUE, lex.order = TRUE)
  out <- ev[!duplicated(key), , drop = FALSE]
  out$lookup_key <- paste(out$symbol, out$contrast, out$modality_class,
                          sep = "\r")
  out$supported <- is.finite(out$fdr) & out$fdr < alpha
  out
}

.fig_axis_effect_best_pathway <- function(crossmodality_report, contrasts) {
  x <- crossmodality_report$pathway$axis_summary
  .fig_require_cols(x, c("axis", "contrast", "n_modalities_present",
                         "n_modalities_sig", "mixed_sign",
                         "consistent_direction", "rank_score"),
                    "crossmodality_report$pathway$axis_summary")
  x$axis <- .fig_axis_effect_standard_axis(x$axis)
  x <- x[!is.na(x$axis) & x$contrast %in% contrasts, , drop = FALSE]
  if (!nrow(x)) return(x)
  x <- x[order(x$axis, x$contrast, -x$n_modalities_sig,
               -x$n_modalities_present, -x$rank_score,
               method = "radix", na.last = TRUE), , drop = FALSE]
  key <- interaction(x$axis, x$contrast, drop = TRUE, lex.order = TRUE)
  out <- x[!duplicated(key), , drop = FALSE]
  out$lookup_key <- paste(out$axis, out$contrast, sep = "\r")
  out
}

.fig_axis_effect_scale <- function(source_feature_type) {
  if (source_feature_type %in% c("tf_activity", "kinase_activity")) {
    "activity_score"
  } else if (source_feature_type %in% c("gene", "protein_group",
                                        "phosphosite_gene",
                                        "phosphosite_gene_corrected")) {
    "log2_effect"
  } else {
    source_feature_type
  }
}

.fig_axis_effect_direction <- function(effect) {
  if (!is.finite(effect) || effect == 0) {
    "none"
  } else if (effect > 0) {
    "positive"
  } else {
    "negative"
  }
}

.fig_axis_effect_spine <- function(crossmodality_report, crossmodality_table,
                                   alpha = 0.10) {
  contrasts <- .fig_axis_effect_contrasts()
  selection <- .fig_axis_effect_selection(crossmodality_report)
  evidence <- .fig_axis_effect_best_evidence(crossmodality_table, contrasts, alpha)
  pathways <- .fig_axis_effect_best_pathway(crossmodality_report, contrasts)
  pairs <- crossmodality_report$clearance$pair_support
  .fig_require_cols(pairs, c("pair", "contrast", "status",
                             "n_coherent_supported_modalities"),
                    "crossmodality_report$clearance$pair_support")
  pairs$lookup_key <- paste(pairs$pair, pairs$contrast, sep = "\r")
  spatial_status <- crossmodality_report$clearance$spatial_decon$status %||% NA_character_

  grid <- merge(selection, data.frame(contrast = contrasts, stringsAsFactors = FALSE),
                by = NULL, sort = FALSE)
  grid <- grid[order(match(grid$axis, .fig_axis_effect_axes()), grid$selection_rank,
                     grid$feature_id, match(grid$contrast, contrasts),
                     match(grid$modality_class, names(.fig_axis_effect_modalities())),
                     method = "radix"), , drop = FALSE]
  rows <- lapply(seq_len(nrow(grid)), function(i) {
    d <- grid[i, , drop = FALSE]
    base <- data.frame(
      axis = d$axis,
      feature_id = d$feature_id,
      feature_label = d$feature_label,
      feature_type = d$feature_type,
      modality_class = d$modality_class,
      modality = d$modality,
      contrast = d$contrast,
      effect = NA_real_,
      effect_scale = "not_estimated",
      fdr = NA_real_,
      support_status = "not_observed",
      direction = "not_observed",
      measured_state = "not_observed",
      source_slot = d$lookup_source,
      source_feature_type = NA_character_,
      source_feature_id = NA_character_,
      selection_rank = d$selection_rank,
      selection_rule = d$selection_rule,
      stringsAsFactors = FALSE
    )
    if (!isTRUE(d$applicable)) {
      base$support_status <- "not_applicable"
      base$direction <- "not_applicable"
      base$measured_state <- "not_applicable"
      base$effect_scale <- "not_applicable"
      return(base)
    }

    if (d$lookup_source %in% c("evidence_symbol", "mechanism_boundary")) {
      key <- paste(d$lookup_symbol, d$contrast, d$modality_class, sep = "\r")
      j <- match(key, evidence$lookup_key)
      if (!is.na(j)) {
        e <- evidence[j, , drop = FALSE]
        base$effect <- e$effect
        base$effect_scale <- .fig_axis_effect_scale(e$feature_type)
        base$fdr <- e$fdr
        base$support_status <- ifelse(e$supported, "supported",
                                      "measured_not_supported")
        base$direction <- .fig_axis_effect_direction(e$effect)
        base$measured_state <- "measured"
        base$source_slot <- paste0("crossmodality_table$evidence:", e$modality_group)
        base$source_feature_type <- e$feature_type
        base$source_feature_id <- e$feature_id
      }
      return(base)
    }

    if (d$lookup_source == "pathway_axis") {
      key <- paste(d$axis, d$contrast, sep = "\r")
      j <- match(key, pathways$lookup_key)
      if (!is.na(j)) {
        p <- pathways[j, , drop = FALSE]
        base$effect <- p$rank_score
        base$effect_scale <- "axis_rank_score"
        base$support_status <- ifelse(p$n_modalities_sig > 0,
                                      "supported", "measured_not_supported")
        dir <- ifelse(p$mixed_sign, "mixed", as.character(p$consistent_direction))
        base$direction <- ifelse(is.na(dir) | !nzchar(dir), "none", dir)
        base$measured_state <- "measured"
        base$source_slot <- "crossmodality_report$pathway$axis_summary"
        base$source_feature_type <- "pathway_axis"
        base$source_feature_id <- p$axis
      }
      return(base)
    }

    if (d$lookup_source == "clearance_pair") {
      key <- paste(d$lookup_symbol, d$contrast, sep = "\r")
      j <- match(key, pairs$lookup_key)
      if (!is.na(j)) {
        p <- pairs[j, , drop = FALSE]
        base$effect <- p$n_coherent_supported_modalities
        base$effect_scale <- "supported_modality_count"
        base$support_status <- ifelse(p$status == "earned", "earned", "not_earned")
        base$direction <- "none"
        base$measured_state <- "measured"
        base$source_slot <- "crossmodality_report$clearance$pair_support"
        base$source_feature_type <- "clearance_pair"
        base$source_feature_id <- p$pair
      }
      return(base)
    }

    if (d$lookup_source == "spatial_decon_boundary") {
      if (identical(spatial_status, "blocked")) {
        base$support_status <- "blocked"
        base$direction <- "blocked"
        base$measured_state <- "blocked"
        base$effect_scale <- "blocked"
        base$source_slot <- "crossmodality_report$clearance$spatial_decon"
        base$source_feature_type <- "spatial_decon_abundance"
        base$source_feature_id <- "geomx_abundance_de"
      } else {
        base$support_status <- "not_observed"
        base$direction <- "not_observed"
        base$measured_state <- "not_observed"
        base$effect_scale <- "not_estimated"
        base$source_slot <- "crossmodality_report$clearance$spatial_decon"
      }
      return(base)
    }
    base
  })

  out <- do.call(rbind, rows)
  out$axis <- factor(out$axis, levels = .fig_axis_effect_axes())
  out$contrast <- factor(out$contrast, levels = contrasts)
  out$modality_class <- factor(out$modality_class,
                               levels = names(.fig_axis_effect_modalities()))
  out$measured_state <- factor(out$measured_state,
                               levels = c("measured", "not_observed",
                                          "blocked", "not_applicable"))
  out$support_status <- factor(out$support_status,
                               levels = c("supported", "earned",
                                          "measured_not_supported",
                                          "not_earned", "not_observed",
                                          "blocked", "not_applicable"))
  out$direction <- factor(out$direction,
                          levels = c("positive", "negative", "mixed", "none",
                                     "not_observed", "blocked", "not_applicable"))
  out <- out[order(as.integer(out$axis), out$selection_rank, out$feature_id,
                   as.integer(out$contrast), as.integer(out$modality_class),
                   method = "radix"), , drop = FALSE]
  rownames(out) <- NULL

  expected <- as.vector(outer(selection$selection_key, contrasts, paste, sep = "\r"))
  observed <- paste(as.character(out$axis), out$feature_id,
                    as.character(out$modality_class), as.character(out$contrast),
                    sep = "\r")
  stopifnot(!anyDuplicated(observed),
            setequal(expected, observed),
            any(out$measured_state == "blocked"),
            any(out$measured_state == "not_observed"),
            any(out$measured_state == "not_applicable"))
  .fig_assert_finite(out[out$measured_state == "measured", , drop = FALSE],
                     "effect", "axis_effect_spine measured effects")
  list(spine = out, selection = selection)
}

.fig_axis_effect_plot_rows <- function(spine, axes = NULL, contrasts = NULL,
                                       feature_types = NULL, feature_ids = NULL,
                                       modality_classes = NULL,
                                       keep_not_applicable = FALSE) {
  .fig_require_cols(spine, c("axis", "feature_id", "feature_label", "feature_type",
                             "modality_class", "contrast", "effect", "effect_scale",
                             "fdr", "support_status", "direction", "measured_state",
                             "selection_rank", "selection_rule"),
                    "axis_effect_spine")
  x <- spine
  if (!is.null(axes)) {
    x <- x[as.character(x$axis) %in% axes, , drop = FALSE]
  }
  if (!is.null(contrasts)) {
    x <- x[as.character(x$contrast) %in% contrasts, , drop = FALSE]
  }
  if (!is.null(feature_types)) {
    x <- x[as.character(x$feature_type) %in% feature_types, , drop = FALSE]
  }
  if (!is.null(feature_ids)) {
    x <- x[as.character(x$feature_id) %in% feature_ids, , drop = FALSE]
  }
  if (!is.null(modality_classes)) {
    x <- x[as.character(x$modality_class) %in% modality_classes, , drop = FALSE]
  }
  if (!keep_not_applicable) {
    x <- x[as.character(x$measured_state) != "not_applicable", , drop = FALSE]
  }
  stopifnot(nrow(x) > 0L)

  axis_labs <- c(
    DAM = "DAM",
    antigen_presentation = "antigen presentation",
    synaptic = "synaptic",
    clearance = "clearance",
    interaction_boundary = "interaction boundary",
    mechanism_boundary = "mechanism boundary"
  )
  contrast_labs <- c(
    nlgf_in_maptki = "amyloid on MAPTKI",
    nlgf_in_p301s = "amyloid on P301S",
    tau_in_nlgf = "tau in NLGF",
    interaction = "interaction"
  )
  modality_labs <- .fig_axis_effect_modalities()

  x$axis_chr <- as.character(x$axis)
  x$contrast_chr <- as.character(x$contrast)
  x$modality_chr <- as.character(x$modality_class)
  x$axis_label_chr <- unname(axis_labs[x$axis_chr])
  x$axis_label_chr[is.na(x$axis_label_chr)] <- gsub("_", " ", x$axis_chr[is.na(x$axis_label_chr)],
                                                    fixed = TRUE)
  x$contrast_label_chr <- unname(contrast_labs[x$contrast_chr])
  x$contrast_label_chr[is.na(x$contrast_label_chr)] <- x$contrast_chr[is.na(x$contrast_label_chr)]
  x$modality_label_chr <- unname(modality_labs[x$modality_chr])
  x$modality_label_chr[is.na(x$modality_label_chr)] <- x$modality_chr[is.na(x$modality_label_chr)]
  x$feature_label_chr <- gsub("_", " ", as.character(x$feature_label), fixed = TRUE)
  x$feature_label_chr <- gsub(" axis$", "", x$feature_label_chr)
  x$plot_effect <- ifelse(is.finite(x$effect), x$effect, 0)
  x$effect_abs <- ifelse(is.finite(x$effect), abs(x$effect), NA_real_)
  x$neg_log10_fdr <- ifelse(is.finite(x$fdr), -log10(pmax(x$fdr, 1e-300)), NA_real_)
  x$measured_state_chr <- as.character(x$measured_state)
  support_chr <- as.character(x$support_status)
  x$plot_status_chr <- ifelse(
    support_chr %in% c("supported", "earned"),
    "supported/earned",
    ifelse(x$measured_state_chr == "measured", "measured, not supported",
           ifelse(x$measured_state_chr == "blocked", "blocked",
                  ifelse(x$measured_state_chr == "not_observed", "not observed",
                         "not applicable")))
  )
  x$effect_sign_chr <- ifelse(as.character(x$direction) %in% c("positive", "negative", "mixed"),
                              as.character(x$direction), "none")

  axis_order <- .fig_axis_effect_axes()
  if (!is.null(axes)) axis_order <- axis_order[axis_order %in% axes]
  contrast_order <- .fig_axis_effect_contrasts()
  if (!is.null(contrasts)) contrast_order <- contrast_order[contrast_order %in% contrasts]
  modality_order <- names(.fig_axis_effect_modalities())
  if (!is.null(modality_classes)) {
    modality_order <- modality_order[modality_order %in% modality_classes]
  }
  feature_order <- unique(x$feature_label_chr[order(match(x$axis_chr, axis_order),
                                                    x$selection_rank, x$feature_id,
                                                    method = "radix", na.last = TRUE)])

  x$axis_label <- factor(x$axis_label_chr, levels = unname(axis_labs[axis_order]))
  x$contrast_label <- factor(x$contrast_label_chr, levels = unname(contrast_labs[contrast_order]))
  x$modality_label <- factor(x$modality_label_chr, levels = unname(modality_labs[modality_order]))
  x$feature_label_plot <- factor(x$feature_label_chr, levels = rev(feature_order))
  x$plot_status <- factor(x$plot_status_chr,
                          levels = c("supported/earned", "measured, not supported",
                                     "not observed", "blocked", "not applicable"))
  x$effect_sign <- factor(x$effect_sign_chr, levels = c("negative", "positive", "mixed", "none"))
  x <- x[order(match(x$axis_chr, axis_order), x$selection_rank, x$feature_id,
               match(x$contrast_chr, contrast_order),
               match(x$modality_chr, modality_order),
               method = "radix", na.last = TRUE), , drop = FALSE]
  rownames(x) <- NULL
  measured <- x[x$measured_state_chr == "measured", , drop = FALSE]
  if (nrow(measured)) {
    .fig_assert_finite(measured, "plot_effect", "axis effect plot measured rows")
  }
  x
}

.fig_crossmodality_amyloid_response_plate <- function(spine) {
  effects <- .fig_axis_effect_plot_rows(
    spine,
    axes = c("DAM", "antigen_presentation"),
    contrasts = c("nlgf_in_maptki", "nlgf_in_p301s"),
    feature_types = "symbol",
    modality_classes = names(.fig_four_modality_classes())
  )
  stopifnot(any(effects$measured_state_chr == "measured"),
            any(effects$axis_chr == "DAM"),
            any(effects$axis_chr == "antigen_presentation"))
  list(effects = effects)
}

.fig_crossmodality_synaptic_clearance_plate <- function(spine) {
  effects <- .fig_axis_effect_plot_rows(
    spine,
    axes = c("synaptic", "clearance"),
    contrasts = c("nlgf_in_maptki", "nlgf_in_p301s"),
    feature_types = "symbol",
    modality_classes = names(.fig_four_modality_classes())
  )
  pairs <- .fig_axis_effect_plot_rows(
    spine,
    axes = "clearance",
    contrasts = .fig_axis_effect_contrasts(),
    feature_types = "pair",
    modality_classes = "clearance_pair"
  )
  stopifnot(any(effects$axis_chr == "synaptic"),
            any(effects$axis_chr == "clearance"),
            any(pairs$feature_type == "pair"))
  list(effects = effects, pairs = pairs)
}

.fig_crossmodality_interaction_boundary_plate <- function(spine) {
  feature_ids <- c(
    "DAM:Apoe", "DAM:Cst7",
    "antigen_presentation:Cd74",
    "synaptic:Syn1", "synaptic:Syp",
    "clearance:Apoe", "clearance:Trem2",
    "boundary:Myc_TF", "boundary:Nfkb1_TF", "boundary:Rela_TF",
    "boundary:Gsk3b_kinase", "boundary:SpatialDecon_abundance"
  )
  effects <- .fig_axis_effect_plot_rows(
    spine,
    axes = c("DAM", "antigen_presentation", "synaptic", "clearance",
             "interaction_boundary", "mechanism_boundary"),
    contrasts = "interaction",
    feature_ids = feature_ids,
    modality_classes = names(.fig_axis_effect_modalities())
  )
  stopifnot(any(effects$axis_chr == "mechanism_boundary"),
            any(effects$measured_state_chr == "blocked"),
            any(effects$measured_state_chr == "measured"))
  list(effects = effects)
}

crossmodality_figure_data <- function(crossmodality_report, geomx_de, bulk_omics_summary,
                                      phospho_de_24m, phospho_corrected_24m,
                                      crossmodality_table = NULL,
                                      alpha = 0.10) {
  stopifnot(is.list(crossmodality_report), is.list(geomx_de),
            is.list(bulk_omics_summary), is.list(phospho_de_24m),
            is.list(phospho_corrected_24m))
  if (is.null(crossmodality_table)) {
    stop("crossmodality_figure_data needs crossmodality_table for four-modality figures",
         call. = FALSE)
  }
  contrasts <- .fig_focus_contrasts()
  geomx_volcano <- .fig_volcano_data(geomx_de$primary$top, contrasts,
                                     alpha = alpha, n_label = 8L)
  geomx_sens <- .fig_sensitivity_counts(geomx_de$primary$top, geomx_de$sensitivity,
                                        contrasts, alpha = alpha)
  geomx_counts <- crossmodality_report$geomx$counts
  bulk_counts <- crossmodality_report$bulk$significant_counts
  .fig_require_cols(geomx_counts, c("contrast", "n_up_fdr_0_10", "n_down_fdr_0_10", "n_fdr_0_10"),
                    "crossmodality geomx counts")
  .fig_require_cols(bulk_counts, c("layer", "contrast", "n_fdr_0_10"), "crossmodality bulk counts")
  geomx_counts$n_sig <- geomx_counts$n_fdr_0_10
  geomx_counts$n_up <- geomx_counts$n_up_fdr_0_10
  geomx_counts$n_down <- geomx_counts$n_down_fdr_0_10
  bulk_counts$n_sig <- bulk_counts$n_fdr_0_10
  phospho_correction <- .fig_phospho_correction(phospho_de_24m$top,
                                                phospho_corrected_24m$top,
                                                contrasts, alpha = alpha)
  four_counts <- .fig_four_modality_counts(crossmodality_table, contrasts, alpha)
  four_pathways <- .fig_four_modality_pathways(crossmodality_report, contrasts)
  four_symbols <- .fig_four_modality_symbols(crossmodality_report, contrasts)
  axis_effect <- .fig_axis_effect_spine(crossmodality_report, crossmodality_table,
                                        alpha)
  amyloid_response <- .fig_crossmodality_amyloid_response_plate(axis_effect$spine)
  synaptic_clearance <- .fig_crossmodality_synaptic_clearance_plate(axis_effect$spine)
  interaction_boundary <- .fig_crossmodality_interaction_boundary_plate(axis_effect$spine)

  out <- list(
    manifest = figure_manifest("crossmodality"),
    axis_effect_spine = axis_effect$spine,
    axis_effect_selection = axis_effect$selection,
    amyloid_response_plate = amyloid_response,
    synaptic_clearance_plate = synaptic_clearance,
    interaction_boundary_plate = interaction_boundary,
    four_modality_counts = four_counts,
    four_modality_pathways = four_pathways,
    four_modality_symbols = four_symbols,
    geomx_counts = geomx_counts,
    geomx_volcano = geomx_volcano,
    geomx_sensitivity = geomx_sens,
    bulk_counts = bulk_counts,
    phospho_raw_corrected = phospho_correction,
    provenance = list(
      source_targets = c("crossmodality_report", "geomx_de",
                         "bulk_omics_summary", "phospho_de_24m",
                         "phospho_corrected_24m", "crossmodality_table"),
      alpha = alpha,
      axis_effect_axes = .fig_axis_effect_axes(),
      axis_effect_contrasts = .fig_axis_effect_contrasts(),
      axis_effect_states = levels(axis_effect$spine$measured_state),
      contract = "compact cross-modality figure data; heavy top tables reduced to points and bins"
    )
  )
  .fig_assert_finite(out$axis_effect_spine[
    out$axis_effect_spine$measured_state == "measured", , drop = FALSE],
    "effect", "axis_effect_spine measured effects")
  .fig_assert_finite(out$geomx_sensitivity,
                     c("n_primary_sig", "n_sensitivity_sig", "n_lost", "n_gained", "n_sign_flip"),
                     "geomx_sensitivity")
  .fig_assert_finite(out$four_modality_counts,
                     c("n_symbols", "n_sig", "n_up_sig", "n_down_sig",
                       "signed_balance", "log_n_sig"),
                     "four_modality_counts")
  .fig_assert_finite(out$four_modality_pathways,
                     c("n_modalities_present", "n_modalities_sig",
                       "n_evidence_groups_sig", "modality_fraction"),
                     "four_modality_pathways")
  .fig_assert_finite(out$four_modality_symbols,
                     c("min_fdr", "n_modalities_present", "n_modalities_sig"),
                     "four_modality_symbols")
  .fig_assert_finite(out$geomx_counts, c("n_up", "n_down", "n_sig"), "geomx_counts")
  .fig_assert_finite(out$bulk_counts, "n_sig", "bulk_counts")
  out
}

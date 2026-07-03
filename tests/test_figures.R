options(warn = 2)
source("R/constants.R")
source("R/utils.R")
source("R/figures.R")
source("tests/helpers.R")

contrasts <- c("tau_alone", "nlgf_in_maptki", "nlgf_in_p301s", "tau_in_nlgf", "interaction")
focus <- c("interaction", "nlgf_in_maptki", "nlgf_in_p301s", "tau_in_nlgf")

manifest <- figure_manifest()
stopifnot(nrow(manifest) == 25L,
          !anyDuplicated(manifest$figure_id),
          all(grepl("^fig-", manifest$figure_id)),
          !any(grepl("_", manifest$figure_id, fixed = TRUE)),
          all(c("microglia", "trajectory", "mechanism", "crossmodality") %in%
                manifest$chapter))
cat("ok - figure_manifest pins 25 hyphenated inline figure ids\n")

symbol_map <- data.frame(
  ensembl = paste0("ENSMUSG", sprintf("%08d", 1:80)),
  symbol = paste0("Gene", 1:80),
  stringsAsFactors = FALSE
)

make_top_gene <- function(n = 40L) {
  stats::setNames(lapply(seq_along(contrasts), function(i) {
    p <- pmin(0.95, seq(0.001, 0.70, length.out = n) + i * 0.001)
    data.frame(
      gene = symbol_map$ensembl[seq_len(n)],
      logFC = seq(-2.2, 2.2, length.out = n) + (i - 3L) * 0.08,
      P.Value = p,
      adj.P.Val = pmin(0.99, p * 2),
      stringsAsFactors = FALSE
    )
  }), contrasts)
}

make_top_symbol <- function(n = 50L) {
  stats::setNames(lapply(seq_along(contrasts), function(i) {
    p <- pmin(0.95, seq(0.001, 0.80, length.out = n) + i * 0.001)
    data.frame(
      symbol = paste0("Sym", seq_len(n)),
      logFC = seq(-2.5, 2.5, length.out = n) + i * 0.05,
      P.Value = p,
      adj.P.Val = pmin(0.99, p * 1.8),
      stringsAsFactors = FALSE
    )
  }), contrasts)
}

make_top_site <- function(n = 45L, shift = 0) {
  stats::setNames(lapply(seq_along(contrasts), function(i) {
    p <- pmin(0.95, seq(0.001, 0.75, length.out = n) + i * 0.001)
    data.frame(
      feature = paste0("row", seq_len(n)),
      site_id = paste0("Kin", seq_len(n), "_S", 100 + seq_len(n)),
      logFC = seq(-1.8, 1.8, length.out = n) + shift + i * 0.03,
      P.Value = p,
      adj.P.Val = pmin(0.99, p * 2),
      stringsAsFactors = FALSE
    )
  }), contrasts)
}

units16 <- as.vector(outer(genotype_levels, sprintf("batch%02d", 1:4), paste, sep = "_"))
counts <- matrix(rep(c(25, 35, 4), length.out = length(units16) * 3L),
                 nrow = length(units16), ncol = 3L,
                 dimnames = list(units16, c("Homeostatic", "DAM", "IFN")))

cf <- data.frame(
  umap_1 = seq(-4, 4, length.out = 48),
  umap_2 = sin(seq_len(48) / 3),
  genotype = factor(rep(genotype_levels, length.out = 48), levels = genotype_levels),
  substate = factor(rep(c("Homeostatic", "DAM", "IFN"), length.out = 48),
                    levels = c("Homeostatic", "DAM", "IFN")),
  Homeostatic_UCell_z = cos(seq_len(48) / 5),
  DAM_UCell_z = sin(seq_len(48) / 4),
  MHC_APC_UCell_z = cos(seq_len(48) / 7),
  stringsAsFactors = FALSE
)
composition_results <- list(
  counts = counts,
  concordance = expand.grid(contrast = contrasts,
                            substate = c("Homeostatic", "DAM", "IFN"),
                            stringsAsFactors = FALSE)
)
composition_results$concordance$dir_logit <- rep(c(-1, 1, 1), length.out = nrow(composition_results$concordance))
composition_results$concordance$sig_logit <- rep(c(FALSE, TRUE), length.out = nrow(composition_results$concordance))
composition_results$concordance$dir_asin <- composition_results$concordance$dir_logit
composition_results$concordance$sig_asin <- composition_results$concordance$sig_logit
composition_results$concordance$dir_concordant <- TRUE
composition_results$concordance$sig_concordant <- TRUE
composition_results$concordance$flag <- FALSE
composition_results$propeller_logit <- expand.grid(contrast = contrasts,
                                                   substate = c("Homeostatic", "DAM", "IFN"),
                                                   stringsAsFactors = FALSE)
composition_results$propeller_logit$t <- seq(-3, 3, length.out = nrow(composition_results$propeller_logit))
composition_results$propeller_logit$fdr_global <- seq(0.01, 0.5, length.out = nrow(composition_results$propeller_logit))

pb_de_microglia <- list(top = make_top_gene())
cell_counts <- t(counts)
pb_de_substate <- list(
  cell_counts = cell_counts,
  per_substate = list(
    Homeostatic = list(status = "fit", n_cells = 800, units = 16, top = make_top_gene()),
    DAM = list(status = "fit", n_cells = 900, units = 16, top = make_top_gene()),
    IFN = list(status = "skip", n_cells = 30, units = 12, reason = "thin"),
    Proliferative = list(status = "skip", n_cells = 0, units = 0, reason = "empty")
  )
)

mf <- microglia_figure_data(list(cell_frame = cf), composition_results,
                            pb_de_microglia, pb_de_substate, symbol_map)
stopifnot(all(figure_manifest("microglia")$slot %in% names(mf)),
          all(c("summary_board", "composition_shift", "composition_forest", "amyloid_volcano") %in% names(mf)),
          nrow(mf$umap_by_substate) == nrow(cf),
          nrow(mf$score_triptych) > 0L,
          nrow(mf$unit_composition) == length(units16) * 3L,
          nrow(mf$whole_de_volcano$bins) > 0L,
          nrow(mf$substate_de_volcano$labels) > 0L)
cat("ok - microglia_figure_data builds every microglia slot from compact inputs\n")

weighted_top <- stats::setNames(lapply(seq_along(contrasts), function(i) {
  data.frame(
    measure = c("mean_pt", "comp_cf", "progression_cf", "cross",
                "within_homeostatic", "within_dam"),
    coef = c(1, 0.8, -0.2, 0.4, 0.1, 0.2) + i / 10,
    ci_l = c(1, 0.8, -0.2, 0.4, 0.1, 0.2) + i / 10 - 0.2,
    ci_r = c(1, 0.8, -0.2, 0.4, 0.1, 0.2) + i / 10 + 0.2,
    p_value = seq(0.01, 0.30, length.out = 6),
    stringsAsFactors = FALSE
  )
}), contrasts)
trajectory_report <- list(
  cell_frame = data.frame(genotype = cf$genotype, substate = cf$substate,
                          on_lineage = TRUE, pt_raw = seq(1, 5, length.out = nrow(cf)),
                          score_axis_pt = seq(-1, 1, length.out = nrow(cf)),
                          stringsAsFactors = FALSE),
  per_unit = data.frame(genotype_batch = units16,
                        genotype = rep(genotype_levels, each = 4),
                        batch = rep(sprintf("batch%02d", 1:4), times = 4),
                        n_cells = 50, sd_pt = 1,
                        mean_pt = seq(2, 5, length.out = 16),
                        median_pt = seq(2, 5, length.out = 16),
                        q90 = seq(3, 6, length.out = 16),
                        frac_past = seq(0.1, 0.8, length.out = 16),
                        within_homeostatic = seq(1, 4, length.out = 16),
                        within_dam = seq(3, 6, length.out = 16),
                        stringsAsFactors = FALSE),
  interaction = data.frame(
    family = c("exploratory", "exploratory", "primary", "exploratory",
               "primary", "exploratory"),
    measure = c("mean_pt", "comp_cf", "progression_cf", "cross",
                "within_homeostatic", "within_dam"),
    coef = c(1.1, 1.4, -0.2, -0.1, 0.05, 0.1),
    ci_l = c(0.4, 0.6, -0.5, -0.4, -0.2, -0.1),
    ci_r = c(1.8, 2.2, 0.1, 0.2, 0.3, 0.3),
    p_value = c(0.01, 0.005, 0.2, 0.5, 0.7, 0.4),
    perm_p = c(0.04, NA, 0.18, NA, 0.6, NA),
    fdr = c(0.04, 0.02, 0.4, 0.6, 0.7, 0.5),
    stringsAsFactors = FALSE
  ),
  weighted_top = weighted_top,
  lineage_per_unit = data.frame(genotype_batch = units16,
                                genotype = rep(genotype_levels, each = 4),
                                n_cells = 60, n_on_lineage = 55,
                                omitted_frac = 5 / 60,
                                stringsAsFactors = FALSE),
  sensitivity = data.frame(variant = c("dims_10", "dims_20"),
                           spearman_vs_primary = c(0.95, 0.97),
                           n_lineages = 1, n_shared = 100,
                           stringsAsFactors = FALSE),
  glmm = list(method = "glmmTMB_beta", estimate = 0.1, p_value = 0.02),
  provenance = list(concordance_rho = 0.7, omitted_frac_overall = 0.05,
                    composition_loading = 1.1, progression_loading = -0.3,
                    cross_loading = 0.2)
)
tf <- trajectory_figure_data(trajectory_report, composition_results)
stopifnot(all(figure_manifest("trajectory")$slot %in% names(tf)),
          all(c("pseudotime_shift", "decomposition", "concordance", "logic_board") %in% names(tf)),
          nrow(tf$pt_density) > 0L,
          all(is.finite(tf$unit_pt_vs_dam$dam_fraction)),
          all(c("comp_cf", "progression_cf", "cross") %in% tf$kitagawa_forest$measure))
cat("ok - trajectory_figure_data joins DAM fraction and Kitagawa rows\n")

mechanism_report <- list(
  pathway_project = expand.grid(pathway = c("DAM", "Homeostatic", "NFkB_Activated_Targets"),
                                population = c("whole_microglia", "DAM"),
                                contrast = contrasts,
                                stringsAsFactors = FALSE),
  pathway_go_top = expand.grid(pathway = paste0("GO_", 1:4),
                               population = c("whole_microglia", "DAM"),
                               contrast = contrasts,
                               stringsAsFactors = FALSE),
  tf_highlights = expand.grid(population = c("whole_microglia", "DAM"),
                              source = c("Myc", "Nfkb1", "Rela"),
                              contrast = contrasts,
                              stringsAsFactors = FALSE),
  nfkb = list(table = expand.grid(population = c("whole_microglia", "DAM"),
                                  contrast = c("interaction", "tau_in_nlgf"),
                                  test = c("target_gsea", "tf_family"),
                                  stringsAsFactors = FALSE),
              verdict = list(status = "discordant", supported = FALSE)),
  kinase = list(table = expand.grid(source = c("Gsk3b", "Camk2a"),
                                    contrast = contrasts,
                                    stringsAsFactors = FALSE))
)
mechanism_report$pathway_project$NES <- seq(-2, 2, length.out = nrow(mechanism_report$pathway_project))
mechanism_report$pathway_project$fdr <- seq(0.01, 0.5, length.out = nrow(mechanism_report$pathway_project))
mechanism_report$pathway_go_top$NES <- seq(-2, 2, length.out = nrow(mechanism_report$pathway_go_top))
mechanism_report$pathway_go_top$fdr <- seq(0.01, 0.5, length.out = nrow(mechanism_report$pathway_go_top))
mechanism_report$pathway_go_top$size <- 20
mechanism_report$tf_highlights$score <- seq(-3, 3, length.out = nrow(mechanism_report$tf_highlights))
mechanism_report$tf_highlights$fdr <- seq(0.01, 0.5, length.out = nrow(mechanism_report$tf_highlights))
mechanism_report$tf_highlights$selection <- "key"
mechanism_report$nfkb$table$score <- seq(-1, 1, length.out = nrow(mechanism_report$nfkb$table))
mechanism_report$nfkb$table$p_value <- 0.2
mechanism_report$nfkb$table$primary_test <- mechanism_report$nfkb$table$population == "whole_microglia"
mechanism_report$nfkb$table$supportive_only <- !mechanism_report$nfkb$table$primary_test
mechanism_report$nfkb$table$primary_family_fdr <- 0.4
mechanism_report$kinase$table$score <- seq(-2, 2, length.out = nrow(mechanism_report$kinase$table))
mechanism_report$kinase$table$fdr <- seq(0.02, 0.6, length.out = nrow(mechanism_report$kinase$table))
mechanism_report$kinase$table$significant <- mechanism_report$kinase$table$fdr < 0.10
mechanism_report$kinase$table$run_index_score <- mechanism_report$kinase$table$score / 2
mechanism_report$kinase$table$run_index_fdr <- pmin(0.9, mechanism_report$kinase$table$fdr * 2)
mechanism_report$kinase$table$run_index_supports <- FALSE
mechanism_report$kinase$table$include_reason <- ifelse(mechanism_report$kinase$table$source == "Gsk3b",
                                                       "gsk3b_carry", "significant")
mecf <- mechanism_figure_data(mechanism_report)
stopifnot(all(figure_manifest("mechanism")$slot %in% names(mecf)),
          all(c("status_board", "project_sets", "tf_interaction") %in% names(mecf)),
          nrow(mecf$project_pathway_heatmap) > 0L,
          nrow(mecf$nfkb_discordance$table) > 0L,
          any(mecf$kinase_heatmap$source == "Gsk3b"))
cat("ok - mechanism_figure_data exposes supported and unsupported mechanism rows\n")

geomx_top <- make_top_symbol()
geomx_de <- list(
  primary = list(top = geomx_top),
  sensitivity = list(
    unblocked = list(status = "fit", top = make_top_symbol()),
    collapsed_bio_unit = list(status = "fit", top = make_top_symbol())
  )
)
crossmodality_report <- list(
  geomx = list(counts = data.frame(contrast = contrasts,
                                   n_features = 100,
                                   n_fdr_0_05 = 1:5,
                                   n_fdr_0_10 = 2:6,
                                   n_up_fdr_0_10 = 1:5,
                                   n_down_fdr_0_10 = 1,
                                   stringsAsFactors = FALSE)),
  bulk = list(
    run_index = expand.grid(layer = c("proteome", "phospho_raw", "phospho_corrected"),
                            contrast = contrasts, stringsAsFactors = FALSE),
    significant_counts = expand.grid(layer = c("proteome", "phospho_raw", "phospho_corrected"),
                                     contrast = contrasts, stringsAsFactors = FALSE),
    anchor_rows = data.frame(layer = "proteome", contrast = rep(focus, each = 4),
                             feature = paste0("F", 1:16), site_id = NA_character_,
                             symbols = paste0("A", 1:16),
                             anchor_symbols = paste0("A", 1:16),
                             anchor_class = rep(c("clearance", "synaptic"), 8),
                             logFC = seq(-2, 2, length.out = 16),
                             t = seq(-3, 3, length.out = 16),
                             P.Value = seq(0.01, 0.5, length.out = 16),
                             adj.P.Val = seq(0.02, 0.8, length.out = 16),
                             stringsAsFactors = FALSE)
  ),
  clearance = list(spatial_decon = list(status = "blocked", action = "attempted",
                                        unresolved_aoi_n = 4L),
                   pair_support = expand.grid(pair = c("Apoe_Trem2", "App_Cd74"),
                                              contrast = focus,
                                              stringsAsFactors = FALSE)),
  divergence = list(axis_symbols = expand.grid(symbol = paste0("S", 1:6),
                                               contrast = focus,
                                               stringsAsFactors = FALSE)),
  pathway = list(axis_summary = expand.grid(axis = c("DAM", "synaptic", "NFkB"),
                                            contrast = focus,
                                            stringsAsFactors = FALSE))
)
crossmodality_report$bulk$run_index$status <- "fit"
crossmodality_report$bulk$run_index$n_primary_sig <- rep(0:4, length.out = nrow(crossmodality_report$bulk$run_index))
crossmodality_report$bulk$run_index$n_lost_or_flipped <- rep(0:3, length.out = nrow(crossmodality_report$bulk$run_index))
crossmodality_report$bulk$significant_counts$n_features <- 100
crossmodality_report$bulk$significant_counts$n_fdr_0_05 <- rep(0:2, length.out = nrow(crossmodality_report$bulk$significant_counts))
crossmodality_report$bulk$significant_counts$n_fdr_0_10 <- rep(1:3, length.out = nrow(crossmodality_report$bulk$significant_counts))
crossmodality_report$bulk$significant_counts$n_up_fdr_0_10 <- 1
crossmodality_report$bulk$significant_counts$n_down_fdr_0_10 <- 0
crossmodality_report$clearance$pair_support$n_sides_measured <- 2
crossmodality_report$clearance$pair_support$modalities_measured <- "GeoMx_spatial;snRNAseq_microglia"
crossmodality_report$clearance$pair_support$coherent_supported_modalities <- ""
crossmodality_report$clearance$pair_support$n_coherent_supported_modalities <- 0
crossmodality_report$clearance$pair_support$microglia_strong <- FALSE
crossmodality_report$clearance$pair_support$status <- "not_earned"
crossmodality_report$divergence$axis_symbols$axis <- rep(c("DAM", "synaptic"), length.out = nrow(crossmodality_report$divergence$axis_symbols))
crossmodality_report$divergence$axis_symbols$axis_member <- TRUE
crossmodality_report$divergence$axis_symbols$n_modalities_present <- 3
crossmodality_report$divergence$axis_symbols$n_modalities_sig <- rep(0:2, length.out = nrow(crossmodality_report$divergence$axis_symbols))
crossmodality_report$divergence$axis_symbols$min_fdr <- seq(0.01, 0.7, length.out = nrow(crossmodality_report$divergence$axis_symbols))
crossmodality_report$divergence$axis_symbols$rank_score <- seq_len(nrow(crossmodality_report$divergence$axis_symbols))
crossmodality_report$pathway$axis_summary$n_modalities_present <- 3
crossmodality_report$pathway$axis_summary$n_modalities_sig <- rep(0:2, length.out = nrow(crossmodality_report$pathway$axis_summary))
crossmodality_report$pathway$axis_summary$rank_score <- seq_len(nrow(crossmodality_report$pathway$axis_summary))

cmf <- crossmodality_figure_data(crossmodality_report, geomx_de, list(),
                                 list(top = make_top_site()),
                                 list(top = make_top_site(shift = 0.2)))
stopifnot(all(figure_manifest("crossmodality")$slot %in% names(cmf)),
          all(c("status_board", "geomx_counts", "bulk_counts", "axis_heatmap") %in% names(cmf)),
          nrow(cmf$geomx_volcano$bins) > 0L,
          nrow(cmf$geomx_sensitivity) == length(focus) * 2L,
          nrow(cmf$phospho_raw_corrected$labels) > 0L,
          all(is.finite(cmf$bulk_run_index$loss_fraction)))
cat("ok - crossmodality_figure_data bins heavy contrasts into compact plotting tables\n")

qc_md <- data.frame(
  genotype = rep(genotype_levels, each = 8),
  batch = rep(rep(sprintf("batch%02d", 1:4), each = 2), times = 4),
  genotype_batch = rep(units16, each = 2),
  nCount_RNA = seq(1000, 5000, length.out = 32),
  nFeature_RNA = seq(300, 1200, length.out = 32),
  percent_mt = seq(1, 8, length.out = 32),
  percent_ribo = seq(5, 18, length.out = 32),
  percent_malat1 = seq(0.5, 3, length.out = 32),
  percent_contam = seq(0.2, 4, length.out = 32),
  stringsAsFactors = FALSE
)
qc <- qc_figure_data(
  list(meta.data = qc_md, dims = c(2000L, nrow(qc_md))),
  list(meta.data = data.frame(genotype = rep(genotype_levels, c(3, 3, 4, 4)),
                              stringsAsFactors = FALSE),
       dims = c(120L, 14L)),
  data.frame(a = 1:4, b = 5:8),
  data.frame(a = 1:5, b = 6:10),
  data.frame(sample = seq_len(16))
)
stopifnot(all(c("modality_table", "genotype_batch", "depth_distribution",
                "fraction_distribution", "metric_bounds") %in% names(qc)),
          nrow(qc$genotype_batch) == 16L,
          all(qc$metric_bounds$within))
cat("ok - qc_figure_data builds compact QC visual slots\n")

rv <- report_visual_data(data.frame(component = "R", version = "test"),
                         qc, mf, tf, mecf, cmf)
stopifnot("report_spine_schematic" %in% names(rv),
          !"synthesis" %in% rv$report_spine_schematic$nodes$node,
          all(c("inputs", "qc", "microglia", "trajectory", "mechanism",
                "crossmodality", "environment") %in% rv$report_spine_schematic$nodes$node),
          rv$source_target_contract$n_manifest_slots[
            rv$source_target_contract$target == "report_visuals"] > 0L)
cat("ok - report_visual_data builds overview spine without synthesis chapter dependency\n")

replacement_manifest <- utils::read.delim(".agent/prose_replacement_manifest.tsv",
                                          stringsAsFactors = FALSE, check.names = FALSE)
coverage <- visual_slot_coverage(replacement_manifest)
stopifnot(coverage$n_missing == 0L)
cat("ok - visual slot map covers every figure/schematic prose-replacement disposition\n")

bad <- mf
bad$whole_de_volcano$bins$x_mid[1] <- Inf
expect_error(.fig_assert_finite(bad$whole_de_volcano$bins, c("x_mid"), "bad bins"), "finite")
cat("ok - figure finite guards fail loud on malformed geom inputs\n")

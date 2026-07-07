options(warn = 2)
source("R/constants.R")
source("R/utils.R")
source("R/design.R")
source("R/figures.R")
source("tests/helpers.R")

contrasts <- c("tau_alone", "nlgf_in_maptki", "nlgf_in_p301s", "tau_in_nlgf", "interaction")

manifest <- figure_manifest()
stopifnot(nrow(manifest) == 11L,
          !anyDuplicated(manifest$figure_id),
          all(grepl("^fig-", manifest$figure_id)),
          !any(grepl("_", manifest$figure_id, fixed = TRUE)),
          all(c("microglia", "trajectory") %in% manifest$chapter))
cat("ok - figure_manifest pins 11 hyphenated inline figure ids\n")

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
          all(c("composition_shift", "composition_forest", "amyloid_volcano") %in% names(mf)),
          nrow(mf$umap_by_substate) == nrow(cf),
          nrow(mf$score_triptych) > 0L,
          nrow(mf$unit_composition) == length(units16) * 3L,
          nrow(mf$whole_de_volcano$points) > 0L,
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
          all(c("pseudotime_shift", "decomposition", "concordance") %in% names(tf)),
          nrow(tf$pt_density) > 0L,
          all(is.finite(tf$unit_pt_vs_dam$dam_fraction)),
          all(c("comp_cf", "progression_cf", "cross") %in% tf$kitagawa_forest$measure))
cat("ok - trajectory_figure_data joins DAM fraction and Kitagawa rows\n")

replacement_manifest <- utils::read.delim(".agent/prose_replacement_manifest.tsv",
                                          stringsAsFactors = FALSE, check.names = FALSE)
coverage <- visual_slot_coverage(replacement_manifest)
stopifnot(coverage$n_missing == 0L)
cat("ok - visual slot map covers every figure/schematic prose-replacement disposition\n")

bad <- mf
bad$whole_de_volcano$bins$x_mid[1] <- Inf
expect_error(.fig_assert_finite(bad$whole_de_volcano$bins, c("x_mid"), "bad bins"), "finite")
cat("ok - figure finite guards fail loud on malformed geom inputs\n")

# --- modality_logfc_scatter_data: y=nlgf_in_maptki, x=nlgf_in_p301s, key-aligned, empirical off-diagonal labels ---
ens5 <- symbol_map$ensembl[1:5]
# snRNAseq: nlgf_in_p301s rows deliberately REVERSED -> exercises the match()-by-key alignment.
bulk_meta <- data.frame(
  genotype = factor(rep(genotype_levels, each = 2), levels = genotype_levels),
  run_index = 1:8,
  row.names = paste0("s", 1:8),
  stringsAsFactors = FALSE
)
pb <- list(top = list(
  nlgf_in_maptki = data.frame(gene = ens5,       logFC = c(1, 2, 3, 4, 5),      stringsAsFactors = FALSE),
  nlgf_in_p301s  = data.frame(gene = rev(ens5),  logFC = c(50, 40, 30, 20, 10), stringsAsFactors = FALSE)))
gx <- list(primary = list(top = list(
  nlgf_in_maptki = data.frame(symbol = paste0("G", 1:4), logFC = c(0.1, 0.2, 0.3, 0.4),     stringsAsFactors = FALSE),
  nlgf_in_p301s  = data.frame(symbol = paste0("G", 1:4), logFC = c(-0.1, -0.2, -0.3, -0.4), stringsAsFactors = FALSE))),
  spatial = list(aoi = data.frame(
    slide = factor(rep(c("slide1", "slide2"), each = 4)),
    genotype = factor(rep(genotype_levels, times = 2), levels = genotype_levels),
    x_coord = rep(c(0, 1, 0, 1), times = 2),
    y_coord = rep(c(0, 0, 1, 1), times = 2),
    signed_response_score = seq(-1, 1, length.out = 8),
    score_abs = abs(seq(-1, 1, length.out = 8)),
    stringsAsFactors = FALSE)))
pr <- list(top = list(
  nlgf_in_maptki = data.frame(feature = paste0("PG", 1:3), gene_first = c("Apoe", "Trem2", ""),
                              gene_symbols = c("Apoe;Trem2", "Trem2", ""),
                              logFC = c(1, 2, 3), P.Value = c(0.01, 0.02, 0.20),
                              adj.P.Val = c(0.03, 0.05, 0.30), stringsAsFactors = FALSE),
  nlgf_in_p301s  = data.frame(feature = paste0("PG", 1:3), gene_first = c("Apoe", "Trem2", ""),
                              gene_symbols = c("Apoe;Trem2", "Trem2", ""),
                              logFC = c(4, 8, 6), P.Value = c(0.001, 0.02, 0.50),
                              adj.P.Val = c(0.004, 0.05, 0.60), stringsAsFactors = FALSE)),
  matrix = matrix(seq(1, 24), nrow = 3, dimnames = list(paste0("PG", 1:3), rownames(bulk_meta))),
  meta = bulk_meta)
ph <- list(top = list(
  nlgf_in_maptki = data.frame(feature = paste0("row", 1:5, "|k"),
                              site_id = c("Mapt_S404", NA, "", "Sfr1_T102", "Lcp1_S5"),
                              gene = c("Mapt", "Apoe", "hMapt", NA, "Lcp1_S5"),
                              logFC = c(0.5, -0.5, 1.0, -3.0, -4.0),
                              P.Value = c(0.01, 0.02, 0.03, 0.004, 0.002),
                              adj.P.Val = c(0.03, 0.04, 0.05, 0.01, 0.008),
                              stringsAsFactors = FALSE),
  nlgf_in_p301s  = data.frame(feature = paste0("row", 1:5, "|k"),
                              site_id = c("Mapt_S404", NA, "", "Sfr1_T102", "Lcp1_S5"),
                              gene = c("Mapt", "Apoe", "hMapt", NA, "Lcp1_S5"),
                              logFC = c(1.5, -1.5, 2.0, 4.0, 5.0),
                              P.Value = c(0.01, 0.04, 0.03, 0.002, 0.001),
                              adj.P.Val = c(0.03, 0.08, 0.06, 0.006, 0.004),
                              stringsAsFactors = FALSE)),
  matrix = matrix(seq(1, 40), nrow = 5, dimnames = list(paste0("row", 1:5, "|k"), rownames(bulk_meta))),
  meta = bulk_meta)
group_sets <- list(
  `Immune activation` = c("Gene1", "Gene5", "Apoe", "Trem2", "Lcp1", "Sfr1"),
  `Neuronal axis` = c("Mapt", "Gene2")
)
ms <- modality_logfc_scatter_data(pb, symbol_map, gx, pr, ph,
                                  group_gene_sets = group_sets,
                                  offdiag_tail_quantile = 0.6,
                                  group_min_genes = 1L,
                                  group_max_groups = 3L)
sn  <- ms$panels$snRNAseq$data
sn_i1 <- match(ens5[1], sn$feature); sn_i5 <- match(ens5[5], sn$feature)
pr_d <- ms$panels$Proteome$data; ph_d <- ms$panels$Phospho$data
stopifnot(
  identical(ms$order, c("snRNAseq", "GeoMx", "Proteome", "Phospho")),
  all(ms$order %in% names(ms$panels)),
  # AXIS INVARIANT: y is the tau-KO amyloid effect (nlgf_in_maptki), x the mutant-tau one (nlgf_in_p301s)
  sn$y[sn_i1] == 1, sn$x[sn_i1] == 10, sn$y[sn_i5] == 5, sn$x[sn_i5] == 50,
  sn$interaction[sn_i1] == sn$x[sn_i1] - sn$y[sn_i1],
  sn$abs_interaction[sn_i1] == abs(sn$interaction[sn_i1]),
  setequal(sn$label, symbol_map$symbol[1:5]),                    # Ensembl gene keys -> symbol labels
  setequal(sn$gene_symbols, symbol_map$symbol[1:5]),
  pr_d$label[match("PG1", pr_d$feature)] == "Apoe",             # gene_first label
  pr_d$gene_symbols[match("PG1", pr_d$feature)] == "Apoe;Trem2", # all group genes feed group scoring
  pr_d$label[match("PG3", pr_d$feature)] == "PG3",              # blank gene_first -> feature id
  nrow(ph_d) == 4L,                                             # duplicate Mapt phosphosites -> one protein point
  ph_d$label[match("phospho_protein:Mapt", ph_d$feature)] == "Mapt",
  ph_d$gene_symbols[match("phospho_protein:Mapt", ph_d$feature)] == "Mapt",
  ph_d$n_phosphosite[match("phospho_protein:Mapt", ph_d$feature)] == 2L,
  ph_d$y[match("phospho_protein:Mapt", ph_d$feature)] == mean(c(0.5, 1.0)),
  ph_d$x[match("phospho_protein:Mapt", ph_d$feature)] == mean(c(1.5, 2.0)),
  ph_d$gene_symbols[match("phospho_protein:Sfr1", ph_d$feature)] == "Sfr1", # blank gene -> parent from site_id
  ph_d$gene_symbols[match("phospho_protein:Lcp1", ph_d$feature)] == "Lcp1", # site-like gene -> parent substitute
  ms$provenance$y_contrast == "nlgf_in_maptki",
  ms$provenance$x_contrast == "nlgf_in_p301s",
  ms$provenance$n_features[["snRNAseq"]] == 5L,
  ms$provenance$n_features[["Phospho"]] == 4L,
  ms$provenance$phospho_site_features == 5L,
  ms$provenance$phospho_parent_proteins == 4L,
  is.data.frame(ms$descriptive$GeoMx$aoi),
  is.data.frame(ms$descriptive$Proteome$pca),
  is.data.frame(ms$descriptive$Proteome$volcano),
  is.data.frame(ms$descriptive$Phospho$volcano),
  is.data.frame(ms$descriptive$Phospho$heatmap),
  nrow(ms$descriptive$Proteome$pca) == nrow(bulk_meta),
  nrow(ms$descriptive$Phospho$heatmap) > 0L)
cat("ok - modality_logfc_scatter_data maps y/x to the two amyloid contrasts and labels each modality\n")

gs <- ms$groups$summary
gg <- ms$groups$selected_genes
panel_cutoffs <- modality_scatter_panel_cutoffs(ms$panels, ms$order, tail_quantile = 0.6)
cutoff_by_row <- unname(panel_cutoffs[as.character(gg$modality)])
stopifnot(
  is.data.frame(gs), nrow(gs) > 0L,
  all(c("modality", "group", "group_label_plot", "n_gene", "score_maptki",
        "score_p301s", "delta", "top_genes", "top_features") %in% names(gs)),
  is.factor(gs$modality), is.factor(gs$group_label_plot),
  any(as.character(gs$group) == "Immune activation"),
  any(gg$gene_symbol == "Gene5" & as.character(gg$modality) == "snRNAseq"),
  all(c("score_feature", "score_label") %in% names(gg)),
  any(gg$gene_symbol == "Lcp1" & gg$score_label == "Lcp1" &
        as.character(gg$modality) == "Phospho"),
  all(!grepl("_[A-Za-z][0-9]", gg$score_label[as.character(gg$modality) == "Phospho"])),
  any(as.character(gg$modality) == "GeoMx"),                   # within-method cutoff, not suppressed by bulk spread
  all(gg$offdiag_distance >= cutoff_by_row),
  all(gs$n_feature <= 3L),
  all.equal(gs$delta, gs$score_p301s - gs$score_maptki, tolerance = 1e-12) == TRUE,
  ms$groups$provenance$group_set_source == "custom functional groups",
  ms$groups$provenance$offdiag_tail_quantile == 0.6,
  all.equal(ms$provenance$offdiag_cutoff, panel_cutoffs, tolerance = 1e-12) == TRUE,
  grepl("within-method", ms$provenance$offdiag_cutoff_source, fixed = TRUE),
  all.equal(ms$groups$provenance$offdiag_cutoff, panel_cutoffs, tolerance = 1e-12) == TRUE,
  grepl("within-method", ms$groups$provenance$offdiag_cutoff_source, fixed = TRUE),
  all(vapply(ms$order, function(m) {
    identical(attr(ms$panels[[m]]$data, "offdiag_cutoff"), unname(panel_cutoffs[[m]]))
  }, logical(1))),
  ms$groups$provenance$min_genes == 1L)
cat("ok - modality_logfc_scatter_data scores within-method off-diagonal labels after phospho parent-protein collapse\n")

# non-finite logFC rows are dropped (finite filter); a missing amyloid contrast fails loud
pb_na <- list(top = list(
  nlgf_in_maptki = data.frame(gene = ens5[1:3], logFC = c(1, NA, 3), stringsAsFactors = FALSE),
  nlgf_in_p301s  = data.frame(gene = ens5[1:3], logFC = c(4, 5, 6),  stringsAsFactors = FALSE)))
ms_na <- modality_logfc_scatter_data(pb_na, symbol_map, gx, pr, ph,
                                     group_gene_sets = group_sets,
                                     offdiag_tail_quantile = 0.6,
                                     group_min_genes = 1L,
                                     group_max_groups = 3L)
stopifnot(nrow(ms_na$panels$snRNAseq$data) == 2L,
          ms_na$provenance$n_features[["snRNAseq"]] == 2L)
pb_missing <- list(top = list(nlgf_in_maptki = pb$top$nlgf_in_maptki))
expect_error(modality_logfc_scatter_data(pb_missing, symbol_map, gx, pr, ph,
                                         group_gene_sets = group_sets,
                                         group_min_genes = 1L))
# a duplicated key in the x-contrast table (setequal ignores multiplicity) must fail loud, not first-match
pb_dupx <- list(top = list(
  nlgf_in_maptki = data.frame(gene = ens5[1:2],                    logFC = c(1, 2),    stringsAsFactors = FALSE),
  nlgf_in_p301s  = data.frame(gene = c(ens5[1], ens5[2], ens5[2]), logFC = c(4, 5, 6), stringsAsFactors = FALSE)))
expect_error(modality_logfc_scatter_data(pb_dupx, symbol_map, gx, pr, ph,
                                         group_gene_sets = group_sets,
                                         group_min_genes = 1L))
cat("ok - modality_logfc_scatter_data drops non-finite pairs, fails loud on missing / duplicated-key contrasts\n")

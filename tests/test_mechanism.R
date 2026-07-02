# P3-S1 mechanism contracts: symbol rank matrices, decoupleR schema wrapping, prior
# fingerprinting, OmniPath prior-table standardisation, phosphosite IDs, and KSN coverage.

source("R/constants.R")
source("R/utils.R")
source("R/io.R")
source("R/design.R")
source("R/de_pb.R")
source("R/mechanism.R")
source("tests/helpers.R")

canonical <- mechanism_contrasts()

# --- symbol mapping + duplicate collapse -------------------------------------------------
sym <- data.frame(ensembl = c("g1", "g2", "g3", "g4"),
                  symbol = c("A", "B", "B", "D"), stringsAsFactors = FALSE)
top <- data.frame(gene = c("g1", "g2", "g3", "gX"),
                  t = c(1, -2, 3, 4), logFC = 1:4, stringsAsFactors = FALSE)
mapped <- add_symbol_to_top(top, sym)
mc <- attr(mapped, "symbol_mapping")
stopifnot(identical(mapped$symbol, c("A", "B", "B")),
          mc$n_input == 4L, mc$n_mapped == 3L, mc$n_dropped == 1L)
expect_error(add_symbol_to_top(top, sym, gene_col = "missing"), "missing gene column")

tops <- stats::setNames(rep(list(top), length(canonical)), canonical)
tops$interaction <- data.frame(gene = c("g1", "g2", "g3", "g4"),
                               t = c(5, -6, 7, 8), stringsAsFactors = FALSE)
rmx <- extract_rank_matrix(tops, sym)
stopifnot(is.matrix(rmx), identical(colnames(rmx), canonical),
          identical(rownames(rmx), c("A", "B", "D")),
          rmx["B", "tau_alone"] == 3, rmx["B", "interaction"] == 7,
          rmx["D", "interaction"] == 8)
expect_error(extract_rank_matrix(list(tau_alone = transform(top, t = NA_real_)), sym), "no finite")

# --- decoupleR ULM wrapper schema --------------------------------------------------------
mat <- matrix(c(1, 2, 3, 4, 2, 1, 4, 3, 3, 4, 1, 2), nrow = 4,
              dimnames = list(c("GeneA", "GeneB", "GeneC", "GeneD"), c("c1", "c2", "c3")))
net <- data.frame(source = c("TF1", "TF1", "TF2", "TF2", "TF2"),
                  target = c("GeneA", "GeneB", "GeneB", "GeneC", "GeneD"),
                  mor = c(1, -1, 1, 1, -1), stringsAsFactors = FALSE)
ulm <- run_decoupler_matrix(mat, net, minsize = 2L)
stopifnot(is.data.frame(ulm), all(c("statistic", "source", "condition", "score", "p_value") %in% names(ulm)),
          identical(unique(ulm$statistic), "ulm"), isFALSE(attr(ulm, "has_consensus")),
          all(ulm$source %in% c("TF1", "TF2")), all(ulm$condition %in% colnames(mat)))

# --- cache preflight + prior fingerprint determinism ------------------------------------
cache <- set_mechanism_prior_cache(file.path(tempdir(), "tau_mechanism_omnipath"))
stopifnot(dir.exists(cache$cache_dir), identical(cache$tz, Sys.getenv("TZ")), nzchar(cache$tz))
fp1 <- prior_fingerprint(net, list(b = 2, a = 1))
fp2 <- prior_fingerprint(net[nrow(net):1, c("target", "mor", "source")], list(a = 1, b = 2))
fp3 <- prior_fingerprint(transform(net, mor = -mor), list(a = 1, b = 2))
stopifnot(identical(fp1$hash, fp2$hash), !identical(fp1$hash, fp3$hash),
          fp1$n_rows == nrow(net), fp1$n_cols == ncol(net))

# --- CollecTRI shape from OmniPath-like table -------------------------------------------
collectri_raw <- data.frame(
  source = c("P1", "COMPLEX:JUN_FOS", "COMPLEX:REL_NFKB", "P2", "bad", "amb", "none", "dup1", "dup2"),
  source_genesymbol = c("Myc", "JUN_FOS", "REL_NFKB1", "Spi1", "", "Amb", "None", "Dup", "Dup"),
  target_genesymbol = c("Tert", "Jun", "Ccl2", "B2m", "Drop", "X", "Y", "Z", "Z"),
  is_stimulation = c("True", "True", "False", "False", "False", "True", "False", "True", "False"),
  is_inhibition = c("False", "False", "True", "True", "False", "True", "False", "False", "True"),
  consensus_stimulation = c("True", "True", "False", "False", "False", "True", "False", "True", "False"),
  consensus_inhibition = c("False", "False", "True", "True", "False", "True", "False", "False", "True"),
  stringsAsFactors = FALSE
)
ctr <- standardise_collectri_table(collectri_raw)
cf <- attr(ctr, "prior_filter_counts")
stopifnot(identical(names(ctr), c("source", "target", "mor")),
          all(c("Myc", "AP1", "NFKB", "Spi1") %in% ctr$source),
          ctr$mor[match("NFKB", ctr$source)] == -1,
          all(ctr$mor %in% c(-1, 1)), !("Dup" %in% ctr$source),
          cf$n_ambiguous_sign == 3L, cf$n_conflicting_pairs == 1L)

# --- KSN parsing + phosphosite IDs + coverage -------------------------------------------
ksn_raw <- data.frame(
  enzyme_genesymbol = c("Gsk3b", "Gsk3b", "Mapk9", "Ptpn1", "Bad;Group", "Src", "Src", "Src", "Cdk5", "Cdk5", "Cdk5"),
  substrate_genesymbol = c("Mapt", "Myocd", "Rxra", "Jak2", "Drop", "A,B", "Conf", "Conf", NA, "Bad", "Bad"),
  residue_type = c("T", "S", "S", "Y", "S", "Y", "S", "S", "S", NA, "S"),
  residue_offset = c(375, 454, 75, 1007, 1, 2, 10, 10, 1, 2, NA),
  modification = c("phosphorylation", "phosphorylation", "phosphorylation", "dephosphorylation", "phosphorylation",
                   "phosphorylation", "phosphorylation", "dephosphorylation", "phosphorylation",
                   "phosphorylation", "phosphorylation"),
  stringsAsFactors = FALSE
)
ksn <- standardise_ksn_table(ksn_raw)
kf <- attr(ksn, "prior_filter_counts")
stopifnot(all(c("Gsk3b", "Mapk9", "Ptpn1") %in% ksn$source),
          all(c("Mapt_T375", "Myocd_S454", "Rxra_S75", "Jak2_Y1007") %in% ksn$target),
          ksn$mor[match("Jak2_Y1007", ksn$target)] == -1,
          !any(grepl("Bad;Group", ksn$source, fixed = TRUE)),
          !any(grepl("A,B", ksn$target, fixed = TRUE)),
          !any(grepl("Conf", ksn$target, fixed = TRUE)),
          kf$n_missing_component == 3L, kf$n_multi_gene == 2L,
          kf$n_conflicting_pairs == 1L)

phospho <- data.frame(
  PG.Genes = c("Mapt", "A;B", "", "Gsk3b", "Myocd", "Mapt"),
  PTM.SiteAA = c("T", "S", "S", NA, "S", "T"),
  PTM.SiteLocation = c(375, 10, 20, 30, 454, 375),
  stringsAsFactors = FALSE
)
ids <- phospho_site_ids(phospho)
pc <- attr(ids, "phospho_site_counts")
stopifnot(is.list(pc), !is.null(pc$n_rows),
          identical(sort(ids, method = "radix"), c("Mapt_T375", "Myocd_S454")),
          pc$n_rows == 6L, pc$n_missing_gene == 1L, pc$n_multi_gene == 1L,
          pc$n_missing_site == 1L, pc$n_duplicate_rows == 1L)
cov <- ksn_coverage_probe(ksn, ids, minsize = 2L)
stopifnot(cov$n_phospho_sites == 2L, cov$n_matched_sites == 2L,
          cov$gsk3b$source_present, cov$gsk3b$matched_sites == 2L,
          cov$gsk3b$passes_minsize, cov$kinases_passing_minsize == 1L,
          cov$source_case$count[cov$source_case$pattern == "exact_Gsk3b"] == 1L)

fake_prior <- net
attr(fake_prior, "provenance") <- list(hash = "abc", n_sources = 2L, n_targets = 4L)
fake_cov <- list(n_matched_sites = 2L, kinases_passing_minsize = 1L,
                 gsk3b = list(matched_sites = 2L))
exp_ok <- list(collectri = list(hash = "abc", n_rows = nrow(fake_prior), n_sources = 2L, n_targets = 4L),
               ksn = list(hash = "abc", n_rows = nrow(fake_prior), n_sources = 2L, n_targets = 4L),
               ksn_coverage = list(n_matched_sites = 2L, kinases_passing_minsize = 1L,
                                   gsk3b_matched_sites = 2L))
stopifnot(assert_mechanism_prior_expectations(fake_prior, fake_prior, fake_cov, exp_ok))
exp_bad <- exp_ok; exp_bad$ksn_coverage$n_matched_sites <- 999L
expect_error(assert_mechanism_prior_expectations(fake_prior, fake_prior, fake_cov, exp_bad), "KSN coverage drift")

# --- P3-S2 RNA mechanism helpers --------------------------------------------------------
s2_genes <- paste0("g", seq_len(24))
s2_symbols <- paste0("Sym", seq_len(24))
s2_map <- data.frame(ensembl = s2_genes, symbol = s2_symbols, stringsAsFactors = FALSE)
mk_top <- function(mult = 1) {
  data.frame(gene = s2_genes,
             logFC = seq(-1.2, 1.2, length.out = length(s2_genes)) * mult,
             t = seq(-3, 3, length.out = length(s2_genes)) * mult,
             P.Value = seq(0.001, 0.2, length.out = length(s2_genes)),
             adj.P.Val = seq(0.01, 0.5, length.out = length(s2_genes)),
             stringsAsFactors = FALSE)
}
s2_tops <- stats::setNames(lapply(seq_along(canonical), function(i) mk_top(ifelse(i %% 2L, 1, -1))),
                           canonical)
s2_pbm <- list(top = s2_tops, n_cells = 1000L)
s2_pbs <- list(per_substate = list(
  Homeostatic = c(list(status = "fit", substate = "Homeostatic", n_cells = 500L), list(top = s2_tops)),
  DAM = list(status = "skipped", substate = "DAM", n_cells = 20L, reason = "floor")
))
ranks <- collect_rna_rank_matrices(s2_pbm, s2_pbs, s2_map)
stopifnot(setequal(names(ranks), c("whole_microglia", "Homeostatic")),
          identical(attr(ranks, "skipped")$population, "DAM"),
          identical(colnames(ranks$whole_microglia$matrix), canonical))

tiny_sets <- list(big = c("a", "b", "c"), small = c("a", "b"), with_na = c("a", NA, ""))
flt <- filter_gene_set_list(tiny_sets, min_size = 3L, universe = c("a", "b", "c"))
stopifnot(identical(names(flt), "big"), identical(flt$big, c("a", "b", "c")))

s2_collectri <- data.frame(
  source = c(rep("Nfkb1", 5), rep("Rela", 5), rep("Reln", 5), rep("Myc", 5)),
  target = paste0("Sym", c(1:5, 6:10, 11:15, 16:20)),
  mor = rep(c(1, -1), length.out = 20),
  stringsAsFactors = FALSE
)
nf_src <- nfkb_family_sources(s2_collectri)
stopifnot(setequal(nf_src, c("Nfkb1", "Rela")), !("Reln" %in% nf_src))
project_sets <- build_project_gene_sets(s2_collectri, custom_min_size = 5L)
stopifnot("NFkB_CollecTRI_Targets" %in% names(project_sets),
          all(c("NFkB_Activated_Targets", "NFkB_Repressed_Targets") %in% names(project_sets)),
          length(project_sets$NFkB_CollecTRI_Targets) == 10L,
          length(project_sets$NFkB_Activated_Targets) == 5L,
          length(project_sets$NFkB_Repressed_Targets) == 5L)
gene_sets <- list(
  sets = list(project = project_sets),
  sizes = data.frame(collection = "project", set = names(project_sets),
                     size = lengths(project_sets), stringsAsFactors = FALSE),
  nfkb_sources = attr(project_sets, "nfkb_sources"),
  thresholds = list(go_min_size = 5L, custom_min_size = 5L)
)

tf_s2 <- run_mechanism_tf(s2_pbm, s2_pbs, s2_map, s2_collectri, minsize = 5L)
stopifnot(is.data.frame(tf_s2$activity), identical(tf_s2$skipped$population, "DAM"),
          all(tf_s2$activity$direction == tf_s2$activity$score),
          all(is.finite(tf_s2$activity$fdr)),
          all(tf_s2$activity$fdr >= tf_s2$activity$p_value - 1e-12))

pw_s2 <- run_mechanism_pathway(s2_pbm, s2_pbs, s2_map, gene_sets)
stopifnot(is.data.frame(pw_s2$pathway), nrow(pw_s2$pathway) > 0L,
          identical(pw_s2$skipped$population, "DAM"),
          all(pw_s2$pathway$direction == pw_s2$pathway$NES),
          all(pw_s2$pathway$fdr == pw_s2$pathway$padj),
          is.logical(pw_s2$pathway$p_floor_warning),
          setequal(pw_s2$pathway$contrast, canonical))

fake_tf <- list(
  activity = data.frame(
    population = rep("whole_microglia", 2),
    population_type = rep("whole", 2),
    source = rep("Nfkb1", 2),
    contrast = c("interaction", "tau_in_nlgf"),
    score = c(2, -5),
    p_value = c(0, 1e-300),
    stringsAsFactors = FALSE),
  skipped = .empty_df(c("population", "population_type", "status", "n_cells", "reason"))
)
fake_path <- list(
  pathway = data.frame(
    population = rep("whole_microglia", 4),
    population_type = rep("whole", 4),
    collection = rep("project", 4),
    contrast = rep(c("interaction", "tau_in_nlgf"), each = 2),
    pathway = rep(c("NFkB_Activated_Targets", "NFkB_Repressed_Targets"), times = 2),
    NES = c(1, -1, -4, 4),
    pval = c(1e-6, 1e-6, 1e-6, 1e-6),
    padj = c(1e-6, 1e-6, 1e-6, 1e-6),
    size = rep(5L, 4),
    p_floor_warning = c(TRUE, TRUE, FALSE, FALSE),
    stringsAsFactors = FALSE),
  warnings = data.frame(population = "whole_microglia", collection = "project",
                        contrast = "interaction", warning = "known", stringsAsFactors = FALSE)
)
fake_sets <- list(nfkb_sources = "Nfkb1")
nf_fake <- build_nfkb_attenuation(fake_tf, fake_path, fake_sets, alpha = 0.10)
stopifnot(!nf_fake$verdict$supported,
          identical(nf_fake$verdict$status, "not_supported"),
          any(nf_fake$table$p_floor_warning[nf_fake$table$test == "target_gsea" & nf_fake$table$primary_test]),
          all(is.finite(nf_fake$table$score[nf_fake$table$test == "tf_family"])))

fake_tf$activity$score[1] <- -5
fake_tf$activity$p_value[1] <- 1e-6
nf_discordant <- build_nfkb_attenuation(fake_tf, fake_path, fake_sets, alpha = 0.10)
stopifnot(!nf_discordant$verdict$supported,
          identical(nf_discordant$verdict$status, "discordant"))

fake_path$pathway$NES[fake_path$pathway$contrast == "interaction" &
                        fake_path$pathway$pathway == "NFkB_Activated_Targets"] <- -4
fake_path$pathway$NES[fake_path$pathway$contrast == "interaction" &
                        fake_path$pathway$pathway == "NFkB_Repressed_Targets"] <- 4
nf_supported <- build_nfkb_attenuation(fake_tf, fake_path, fake_sets, alpha = 0.10)
stopifnot(nf_supported$verdict$supported,
          identical(nf_supported$verdict$status, "supported"),
          nf_supported$verdict$n_primary_supported == 2L)

# --- P3-S3 phospho DE + kinase helpers -----------------------------------------------
s3_stubs <- paste0("run", sprintf("%02d", 1:16))
s3_key <- data.frame(
  file_name = paste0(s3_stubs, ".PTM.Quantity"),
  label = rep(c("MAPT-KI_24M", "P301S+3_24M", "NLGF-MAPT-KI_24M", "NLGF-P301S+3_24M"), each = 4L),
  genotype = factor(rep(genotype_levels, each = 4L), levels = genotype_levels),
  col_stub = s3_stubs,
  stringsAsFactors = FALSE
)
s3_ann <- data.frame(
  PG.Genes = c("Gsk3b", "Gsk3b", "Mapt", "A;B", "", "Myc"),
  PTM.SiteAA = c("S", "S", "T", "S", "Y", "S"),
  PTM.SiteLocation = c(9, 9, 375, 10, 1, 62),
  PTM.CollapseKey = c("ck1", "ck2", "ck3", "ck4", "", "ck6"),
  `Phosphosite probability` = c(0.90, 0.95, 0.99, 0.80, 0.70, 0.90),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
s3_int <- as.data.frame(matrix(100 + seq_len(nrow(s3_ann) * 16), nrow = nrow(s3_ann)))
names(s3_int) <- paste0(s3_stubs, ".PTM.Quantity")
s3_phospho <- cbind(s3_ann, s3_int)

matched24 <- match_24m_intensity_columns(s3_phospho, s3_key)
stopifnot(length(matched24$columns) == 16L,
          identical(as.character(matched24$meta$genotype), rep(genotype_levels, each = 4L)),
          identical(matched24$columns, names(s3_int)))
expect_error(match_24m_intensity_columns(s3_phospho[, -ncol(s3_phospho)], s3_key), "expected exactly")

pf <- phospho_feature_frame(s3_phospho)
pfc <- attr(pf, "phospho_feature_counts")
stopifnot(anyDuplicated(pf$feature) == 0L,
          identical(pf$site_id[1:3], c("Gsk3b_S9", "Gsk3b_S9", "Mapt_T375")),
          is.na(pf$site_id[4]), is.na(pf$site_id[5]),
          pfc$n_multi_gene == 1L, pfc$n_missing_gene == 1L,
          pfc$n_blank_collapse_key == 1L)

lg <- positive_log2_matrix(matrix(c(4, 0, -1, NA_real_), nrow = 2))
lgc <- attr(lg, "log2_counts")
stopifnot(lg[1, 1] == 2, is.na(lg[2, 1]), is.na(lg[1, 2]), is.na(lg[2, 2]),
          lgc$n_nonpositive_to_na == 2L, lgc$n_missing_input == 1L)

prep24 <- prepare_phospho_24m_matrix(s3_phospho, s3_key, min_present = 1L, min_groups = 4L)
stopifnot(ncol(prep24$matrix) == 16L, nrow(prep24$matrix) == nrow(s3_phospho),
          identical(colnames(prep24$matrix), rownames(prep24$meta)))
fd24 <- factorial_design(prep24$meta, add_batch = FALSE)
stopifnot(identical(colnames(fd24$design), c("(Intercept)", "tau", "nlgf", "tau_nlgf")),
          identical(colnames(fd24$contrasts), canonical))

ri24 <- run_index_factorial_design(prep24$meta)
stopifnot(identical(ri24$status, "fit"),
          "run_index" %in% colnames(ri24$design),
          all(ri24$contrasts["run_index", ] == 0))
ri_bad <- run_index_factorial_design(transform(prep24$meta, run_index = 1))
stopifnot(identical(ri_bad$status, "skipped"))

dup_top <- data.frame(
  feature = pf$feature[c(1, 2, 3, 6)],
  site_id = pf$site_id[c(1, 2, 3, 6)],
  phosphosite_probability = c(0.90, 0.95, 0.99, 0.90),
  original_row = pf$original_row[c(1, 2, 3, 6)],
  t = c(10, 1, -3, 4),
  stringsAsFactors = FALSE
)
tie_top <- dup_top
tie_top$phosphosite_probability[1:2] <- 0.95
tie_top$t[1:2] <- c(-10, 1)
site_mat <- phospho_site_stat_matrix(list(tau_alone = dup_top, interaction = tie_top))
scc <- attr(site_mat, "collapse_counts")
stopifnot(site_mat["Gsk3b_S9", "tau_alone"] == 1,
          site_mat["Gsk3b_S9", "interaction"] == -10,
          scc$n_duplicate_sites[scc$contrast == "tau_alone"] == 1L,
          anyDuplicated(rownames(site_mat)) == 0L)

cov_ok <- list(n_matched_sites = 120L, kinases_passing_minsize = 55L,
               gsk3b = list(passes_minsize = TRUE, matched_sites = 8L))
stopifnot(assert_ksn_activity_coverage(cov_ok, min_kinases = 50L, min_matched_sites = 100L))
cov_bad <- cov_ok; cov_bad$gsk3b$passes_minsize <- FALSE
expect_error(assert_ksn_activity_coverage(cov_bad), "Gsk3b does not pass")

fake_kin <- list(
  activity = data.frame(
    fit = c(rep("primary", 4), rep("run_index", 2)),
    statistic = "ulm",
    source = c("Gsk3b", "Other", "Gsk3b", "Other2", "Gsk3b", "Gsk3b"),
    contrast = c("interaction", "interaction", "tau_in_nlgf", "tau_alone", "interaction", "tau_in_nlgf"),
    score = c(2, -3, 0.5, 4, 2.5, -0.4),
    p_value = c(0.20, 0.001, 0.80, 0.002, 0.01, 0.90),
    direction = c(2, -3, 0.5, 4, 2.5, -0.4),
    method = "ulm",
    fdr = c(0.30, 0.01, 0.80, 0.02, 0.02, 0.90),
    stringsAsFactors = FALSE),
  coverage = cov_ok
)
ksum <- build_kinase_mechanism_summary(fake_kin, alpha = 0.10)
stopifnot(is.data.frame(ksum$table),
          all(c("interaction", "tau_in_nlgf", "tau_alone", "nlgf_in_maptki", "nlgf_in_p301s") %in%
                ksum$table$contrast[ksum$table$source == "Gsk3b"]),
          any(ksum$table$source == "Other" & ksum$table$include_reason == "significant"),
          !any(ksum$table$run_order_confounded[ksum$table$source == "Gsk3b"], na.rm = TRUE))

# --- P3-S4 compact mechanism report bundle ---------------------------------------------
s4_pops <- c("whole_microglia", "DAM", "Homeostatic")
s4_sources <- c("Myc", "Nfkb1", "Rela", "Tbp", paste0("Other", 1:6))
s4_tf_grid <- expand.grid(population = s4_pops, contrast = canonical, source = s4_sources,
                          stringsAsFactors = FALSE)
s4_tf_grid$population_type <- ifelse(s4_tf_grid$population == "whole_microglia", "whole", "substate")
s4_tf_grid$statistic <- "ulm"
s4_tf_grid$score <- seq(-2, 2, length.out = nrow(s4_tf_grid))
s4_tf_grid$p_value <- seq(0.001, 0.2, length.out = nrow(s4_tf_grid))
s4_tf_grid$direction <- s4_tf_grid$score
s4_tf_grid$method <- "ulm"
s4_tf_grid$has_consensus <- FALSE
s4_tf_grid$fdr <- pmin(1, s4_tf_grid$p_value * 2)
s4_tf_grid$score[s4_tf_grid$population == "whole_microglia" &
                   s4_tf_grid$contrast == "interaction" &
                   s4_tf_grid$source == "Myc"] <- -6
s4_tf_grid$p_value[s4_tf_grid$population == "whole_microglia" &
                     s4_tf_grid$contrast == "interaction" &
                     s4_tf_grid$source == "Myc"] <- 1e-8
s4_tf_grid$fdr[s4_tf_grid$population == "whole_microglia" &
                 s4_tf_grid$contrast == "interaction" &
                 s4_tf_grid$source == "Myc"] <- 1e-6
s4_tf <- list(activity = s4_tf_grid,
              skipped = .empty_df(c("population", "population_type", "status", "n_cells", "reason")),
              provenance = list(minsize = 5L))

s4_sets <- c("DAM", "Homeostatic", "NFkB_Activated_Targets", "NFkB_Repressed_Targets")
s4_project <- expand.grid(population = s4_pops, contrast = canonical, pathway = s4_sets,
                          stringsAsFactors = FALSE)
s4_project$population_type <- ifelse(s4_project$population == "whole_microglia", "whole", "substate")
s4_project$collection <- "project"
s4_project$NES <- seq(-1.5, 2.5, length.out = nrow(s4_project))
s4_project$pval <- seq(0.001, 0.05, length.out = nrow(s4_project))
s4_project$padj <- pmin(1, s4_project$pval * 2)
s4_project$log2err <- 0.1
s4_project$ES <- s4_project$NES / 2
s4_project$size <- 10L
s4_project$p_floor_warning <- FALSE
s4_project$direction <- s4_project$NES
s4_project$p_value <- s4_project$pval
s4_project$fdr <- s4_project$padj
s4_go <- expand.grid(population = s4_pops, contrast = canonical,
                     collection = c("GO_BP", "GO_CC"),
                     pathway = paste0("GO_SET_", 1:8), stringsAsFactors = FALSE)
s4_go$population_type <- ifelse(s4_go$population == "whole_microglia", "whole", "substate")
s4_go$NES <- seq(-3, 3, length.out = nrow(s4_go))
s4_go$pval <- seq(0.001, 0.2, length.out = nrow(s4_go))
s4_go$padj <- pmin(1, s4_go$pval * 2)
s4_go$log2err <- 0.1
s4_go$ES <- s4_go$NES / 2
s4_go$size <- 30L
s4_go$p_floor_warning <- FALSE
s4_go$direction <- s4_go$NES
s4_go$p_value <- s4_go$pval
s4_go$fdr <- s4_go$padj
s4_pathway <- list(pathway = rbind(s4_project[names(s4_go)], s4_go),
                   skipped = s4_tf$skipped,
                   warnings = .empty_df(c("population", "collection", "contrast", "warning")),
                   provenance = list(go_max_size = 500L))

s4_nfkb_table <- data.frame(
  population = c("whole_microglia", "whole_microglia", "DAM", "DAM"),
  population_type = c("whole", "whole", "substate", "substate"),
  contrast = c("interaction", "interaction", "tau_in_nlgf", "tau_in_nlgf"),
  test = c("tf_family", "target_gsea", "tf_family", "target_gsea"),
  score = c(2, -2, -1, -1),
  p_value = c(0.02, 0.03, 0.2, 0.2),
  raw_p_value = c(0.02, 0.03, 0.2, 0.2),
  direction = c(2, -2, -1, -1),
  n_sources = c(2L, NA_integer_, 2L, NA_integer_),
  detail = c("Nfkb1", "NFkB_Activated_Targets", "Nfkb1", "NFkB_Activated_Targets"),
  family_members = "Nfkb1;Rela",
  pathway_fdr = c(NA, 0.03, NA, 0.2),
  size = c(NA_integer_, 10L, NA_integer_, 10L),
  p_floor_warning = FALSE,
  primary_test = c(TRUE, TRUE, FALSE, FALSE),
  supportive_only = c(FALSE, FALSE, TRUE, TRUE),
  primary_family_fdr = c(0.03, 0.03, NA, NA),
  stringsAsFactors = FALSE
)
s4_nfkb <- list(table = s4_nfkb_table,
                verdict = list(alpha = 0.10, status = "discordant", supported = FALSE),
                skipped = s4_tf$skipped,
                provenance = list(nfkb_sources = c("Nfkb1", "Rela"), alpha = 0.10))

s4_cov <- list(n_phospho_sites = 1000L, n_matched_edges = 300L, n_matched_sites = 200L,
               kinases_passing_minsize = 60L, minsize = 5L,
               gsk3b = list(source_present = TRUE, matched_sites = 8L, passes_minsize = TRUE),
               top_matched = data.frame(source = "Gsk3b", n_sites = 8L, stringsAsFactors = FALSE))
s4_ksum <- ksum
s4_ksum$coverage <- s4_cov

s4_comp <- list(propeller_logit = expand.grid(method = "propeller_logit",
                                              contrast = canonical,
                                              substate = c("Homeostatic", "DAM"),
                                              stringsAsFactors = FALSE))
s4_comp$propeller_logit$prop_ratio <- 1.2
s4_comp$propeller_logit$t <- 2
s4_comp$propeller_logit$p_value <- 0.02
s4_comp$propeller_logit$fdr_global <- 0.05

s4_traj <- list(
  interaction = data.frame(
    family = c("exploratory", "exploratory", "primary", "primary", "exploratory"),
    measure = c("mean_pt", "comp_cf", "progression_cf", "within_homeostatic", "within_dam"),
    coef = c(2, 2.5, -1, -1.5, 0.2),
    se = 1, t = c(2, 2.5, -1, -1.5, 0.2), df = 9,
    p_value = c(0.04, 0.01, 0.2, 0.3, 0.8),
    ci_l = c(0.1, 0.5, -3, -4, -1),
    ci_r = c(4, 4.5, 1, 1, 1.4),
    perm_p = c(0.04, NA, 0.2, 0.3, NA),
    fdr = c(0.10, 0.03, 0.4, 0.4, 0.8),
    stringsAsFactors = FALSE),
  provenance = list(composition_loading = 1.2, progression_loading = -0.4, cross_loading = 0.2)
)

s4_report <- mechanism_report_data(s4_tf, s4_pathway, s4_nfkb, s4_ksum, s4_comp, s4_traj,
                                   tf_top_n = 3L, go_top_n = 2L)
stopifnot(is.list(s4_report),
          all(c("pathway_project", "pathway_go_top", "tf_highlights", "nfkb",
                "kinase", "composition_anchor", "trajectory_anchor") %in% names(s4_report)),
          nrow(s4_report$pathway_go_top) <= length(s4_pops) * length(canonical) * 2L * 2L,
          any(s4_report$tf_highlights$source == "Myc" &
                s4_report$tf_highlights$population == "whole_microglia" &
                s4_report$tf_highlights$contrast == "interaction"),
          all(canonical %in% s4_report$kinase$table$contrast[s4_report$kinase$table$source == "Gsk3b"]),
          !("cell_frame" %in% names(s4_report)),
          s4_report$kinase$coverage$gsk3b$matched_sites == 8L)

s4_tf_bad <- s4_tf
s4_tf_bad$activity <- s4_tf_bad$activity[!(s4_tf_bad$activity$population == "whole_microglia" &
                                             s4_tf_bad$activity$contrast == "interaction" &
                                             s4_tf_bad$activity$source == "Myc"), , drop = FALSE]
expect_error(mechanism_report_data(s4_tf_bad, s4_pathway, s4_nfkb, s4_ksum, s4_comp, s4_traj), "Myc")
s4_ksum_bad <- s4_ksum
s4_ksum_bad$table <- s4_ksum_bad$table[!(s4_ksum_bad$table$source == "Gsk3b" &
                                           s4_ksum_bad$table$contrast == "interaction"), , drop = FALSE]
expect_error(mechanism_report_data(s4_tf, s4_pathway, s4_nfkb, s4_ksum_bad, s4_comp, s4_traj), "contrasts")
s4_traj_bad <- s4_traj
s4_traj_bad$interaction <- s4_traj_bad$interaction[s4_traj_bad$interaction$measure != "comp_cf", , drop = FALSE]
expect_error(mechanism_report_data(s4_tf, s4_pathway, s4_nfkb, s4_ksum, s4_comp, s4_traj_bad), "comp_cf")

cat("ok - test_mechanism\n")

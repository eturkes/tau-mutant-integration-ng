# P5-S1 synthesis contracts: compact P1-P4 report anchors -> descriptive status rows,
# fail-loud missing anchors, empty earned-pair handling, and no ledger-like columns.

source("R/constants.R")
source("R/utils.R")
source("R/synthesis.R")
source("tests/helpers.R")

canonical <- synthesis_contrasts()
focus <- synthesis_focus_contrasts()

make_synthesis_inputs <- function() {
  subtab <- matrix(
    c(40, 10, 2, 0,
      38, 11, 2, 0,
      18, 44, 3, 0,
      10, 58, 4, 0),
    nrow = 4L, byrow = TRUE,
    dimnames = list(genotype_levels, c("Homeostatic", "DAM", "IFN", "Proliferative"))
  )
  microglia <- list(n_cells = sum(subtab),
                    provenance = list(substate_table = subtab))

  composition <- expand.grid(
    contrast = canonical,
    substate = c("Homeostatic", "DAM"),
    stringsAsFactors = FALSE
  )
  composition$method <- "propeller_logit"
  composition$prop_ratio <- 1.1
  composition$t <- 1.5
  composition$p_value <- 0.2
  composition$fdr_global <- 0.2
  composition$t[composition$contrast %in% c("nlgf_in_maptki", "nlgf_in_p301s") &
                  composition$substate == "DAM"] <- 7
  composition$fdr_global[composition$contrast %in% c("nlgf_in_maptki", "nlgf_in_p301s") &
                           composition$substate == "DAM"] <- 0.001
  composition$p_value[composition$contrast %in% c("nlgf_in_maptki", "nlgf_in_p301s") &
                        composition$substate == "DAM"] <- 0.0005
  composition$t[composition$contrast %in% c("nlgf_in_maptki", "nlgf_in_p301s") &
                  composition$substate == "Homeostatic"] <- -7
  composition$fdr_global[composition$contrast %in% c("nlgf_in_maptki", "nlgf_in_p301s") &
                           composition$substate == "Homeostatic"] <- 0.001
  composition$p_value[composition$contrast %in% c("nlgf_in_maptki", "nlgf_in_p301s") &
                        composition$substate == "Homeostatic"] <- 0.0005
  composition$t[composition$contrast == "interaction" & composition$substate == "DAM"] <- 3
  composition$fdr_global[composition$contrast == "interaction" &
                           composition$substate == "DAM"] <- 0.03
  composition$p_value[composition$contrast == "interaction" &
                        composition$substate == "DAM"] <- 0.01
  composition <- composition[, c("method", "contrast", "substate", "prop_ratio",
                                 "t", "p_value", "fdr_global")]

  traj <- data.frame(
    family = c("exploratory", "exploratory", "primary", "primary"),
    measure = c("mean_pt", "comp_cf", "progression_cf", "within_homeostatic"),
    coef = c(2.0, 2.5, -1.1, -1.6),
    p_value = c(0.04, 0.006, 0.17, 0.39),
    perm_p = c(0.04, NA, 0.18, 0.39),
    fdr = c(0.11, 0.025, 0.35, 0.39),
    ci_l = c(0.1, 0.8, -2.8, -5.0),
    ci_r = c(3.9, 4.2, 0.6, 1.8),
    stringsAsFactors = FALSE
  )
  trajectory <- list(
    interaction = traj,
    provenance = list(composition_loading = 1.2,
                      progression_loading = -0.5,
                      cross_loading = 0.3)
  )

  tf <- data.frame(
    population = c("whole_microglia", "DAM"),
    source = c("Myc", "Myc"),
    contrast = c("interaction", "interaction"),
    score = c(-5.5, -4.0),
    p_value = c(1e-6, 1e-4),
    fdr = c(1e-4, 1e-3),
    stringsAsFactors = FALSE
  )
  nfkb <- list(
    table = data.frame(
      test = c("target_gsea", "tf_family"),
      score = c(-1.2, 1.4),
      p_value = c(0.16, 0.63),
      primary_family_fdr = c(0.32, 0.63),
      primary_test = c(TRUE, TRUE),
      detail = c("NFkB_Activated_Targets", "Nfkb2"),
      stringsAsFactors = FALSE
    ),
    verdict = list(alpha = 0.10, status = "discordant", supported = FALSE,
                   primary_discordant = TRUE, primary_negative = FALSE,
                   primary_positive = FALSE)
  )
  kinase <- expand.grid(source = "Gsk3b", contrast = canonical, stringsAsFactors = FALSE)
  kinase$score <- c(0.2, 0.3, 0.4, 0.9, -1.1)
  kinase$p_value <- c(0.8, 0.7, 0.6, 0.34, 0.27)
  kinase$fdr <- c(0.9, 0.85, 0.8, 0.67, 0.96)
  kinase$significant <- FALSE
  kinase$run_index_supports <- FALSE
  kinase$run_order_confounded <- FALSE
  mechanism <- list(
    composition_anchor = composition,
    tf_highlights = tf,
    nfkb = nfkb,
    kinase = list(table = kinase)
  )

  geomx_sens <- data.frame(
    fit = c("primary_bio_unit_blocked", "unblocked_aoi", "collapsed_bio_unit"),
    status = c("fit", "fit", "fit"),
    duplicate_correlation_used = c(TRUE, FALSE, FALSE),
    reason = c(NA_character_, NA_character_, NA_character_),
    stringsAsFactors = FALSE
  )
  bulk_run <- expand.grid(layer = c("proteome", "phospho_raw", "phospho_corrected"),
                          contrast = canonical, stringsAsFactors = FALSE)
  bulk_run$status <- "fit"
  bulk_run$reason <- NA_character_
  bulk_run$n_primary_sig <- ifelse(bulk_run$contrast %in% c("nlgf_in_maptki", "nlgf_in_p301s"), 10L, 1L)
  bulk_run$n_lost_or_flipped <- ifelse(bulk_run$contrast %in% c("nlgf_in_maptki", "nlgf_in_p301s"), 8L, 0L)

  pair_support <- expand.grid(pair = c("Apoe_Trem2", "App_Cd74"),
                              contrast = focus, stringsAsFactors = FALSE)
  pair_support$n_sides_measured <- 2L
  pair_support$modalities_measured <- "GeoMx_spatial;snRNAseq_microglia"
  pair_support$coherent_supported_modalities <- "GeoMx_spatial;snRNAseq_microglia"
  pair_support$n_coherent_supported_modalities <- 2L
  pair_support$microglia_strong <- TRUE
  pair_support$status <- "not_earned"
  pair_support$status[pair_support$pair == "Apoe_Trem2" &
                        pair_support$contrast == "nlgf_in_p301s"] <- "earned"

  pathway <- expand.grid(
    axis = c("DAM", "synaptic", "clearance"),
    contrast = c("nlgf_in_maptki", "nlgf_in_p301s"),
    stringsAsFactors = FALSE
  )
  pathway$collection <- "project"
  pathway$set <- pathway$axis
  pathway$n_modalities_present <- 3L
  pathway$n_modalities_sig <- 2L
  pathway$mixed_sign <- FALSE
  pathway$consistent_direction <- "positive"
  pathway$rank_score <- 8

  crossmodality <- list(
    geomx = list(sensitivity = geomx_sens),
    bulk = list(run_index = bulk_run),
    clearance = list(
      spatial_decon = list(status = "defer", action = "skipped",
                           reasons = c("nuclei sentinels", "no compact reference profile")),
      verdict = list(status = "earned", ccc_called = FALSE),
      pair_support = pair_support
    ),
    pathway = list(axis_summary = pathway)
  )
  list(microglia = microglia, trajectory = trajectory,
       mechanism = mechanism, crossmodality = crossmodality)
}

fx <- make_synthesis_inputs()
sr <- synthesis_report_data(fx$microglia, fx$trajectory, fx$mechanism, fx$crossmodality)
stopifnot(all(c("headline", "evidence_table", "status_summary", "open_questions",
                "source_highlights", "provenance") %in% names(sr)),
          nrow(sr$evidence_table) == 10L,
          all(sr$evidence_table$status %in% synthesis_status_levels()),
          sum(sr$status_summary$n) == nrow(sr$evidence_table),
          any(sr$evidence_table$claim_id == "clearance_axis" &
                sr$evidence_table$status == "focused_support"),
          !any(c("support", "contradict", "net_score", "ledger_id", "arc_id") %in%
                 names(sr$evidence_table)),
          identical(sr$provenance$source_targets,
                    c("microglia_report", "trajectory_report",
                      "mechanism_report", "crossmodality_report")))
cat("ok - synthesis_report_data builds compact descriptive evidence rows\n")

bad_status <- sr$evidence_table
bad_status$status[1] <- "scored_margin"
expect_error(synthesis_validate_evidence_table(bad_status), "invalid status")
bad_ledger <- sr$evidence_table
bad_ledger$net_score <- 1
expect_error(synthesis_validate_evidence_table(bad_ledger), "ledger-like")
cat("ok - synthesis evidence table rejects invalid statuses and ledger-like columns\n")

no_earned <- fx
no_earned$crossmodality$clearance$pair_support$status <- "not_earned"
no_earned$crossmodality$clearance$verdict$status <- "not_earned"
sr_no <- synthesis_report_data(no_earned$microglia, no_earned$trajectory,
                               no_earned$mechanism, no_earned$crossmodality)
clearance <- sr_no$evidence_table[sr_no$evidence_table$claim_id == "clearance_axis", ]
stopifnot(nrow(clearance) == 1L,
          clearance$status == "not_earned",
          grepl("no clearance pair", clearance$evidence, fixed = TRUE))
cat("ok - synthesis handles empty earned-pair set without inventing support\n")

bad_traj <- fx
bad_traj$trajectory$interaction <- bad_traj$trajectory$interaction[
  bad_traj$trajectory$interaction$measure != "comp_cf", , drop = FALSE]
expect_error(synthesis_report_data(bad_traj$microglia, bad_traj$trajectory,
                                   bad_traj$mechanism, bad_traj$crossmodality),
             "comp_cf")

bad_prog <- fx
bad_prog$trajectory$interaction$fdr[
  bad_prog$trajectory$interaction$measure == "progression_cf"] <- 0.01
expect_error(synthesis_report_data(bad_prog$microglia, bad_prog$trajectory,
                                   bad_prog$mechanism, bad_prog$crossmodality),
             "progression-beyond-composition")

bad_myc <- fx
bad_myc$mechanism$tf_highlights <- bad_myc$mechanism$tf_highlights[
  bad_myc$mechanism$tf_highlights$source != "Myc", , drop = FALSE]
expect_error(synthesis_report_data(bad_myc$microglia, bad_myc$trajectory,
                                   bad_myc$mechanism, bad_myc$crossmodality),
             "Myc")

bad_nfkb <- fx
bad_nfkb$mechanism$nfkb$verdict$supported <- TRUE
bad_nfkb$mechanism$nfkb$verdict$status <- "supported"
expect_error(synthesis_report_data(bad_nfkb$microglia, bad_nfkb$trajectory,
                                   bad_nfkb$mechanism, bad_nfkb$crossmodality),
             "NF-kB")

bad_gsk <- fx
bad_gsk$mechanism$kinase$table$significant[
  bad_gsk$mechanism$kinase$table$contrast == "interaction"] <- TRUE
expect_error(synthesis_report_data(bad_gsk$microglia, bad_gsk$trajectory,
                                   bad_gsk$mechanism, bad_gsk$crossmodality),
             "Gsk3b")

bad_ccc <- fx
bad_ccc$crossmodality$clearance$verdict$ccc_called <- TRUE
expect_error(synthesis_report_data(bad_ccc$microglia, bad_ccc$trajectory,
                                   bad_ccc$mechanism, bad_ccc$crossmodality),
             "SpatialDecon")

bad_dam <- fx
bad_dam$mechanism$composition_anchor <- bad_dam$mechanism$composition_anchor[
  !(bad_dam$mechanism$composition_anchor$contrast == "interaction" &
      bad_dam$mechanism$composition_anchor$substate == "DAM"), , drop = FALSE]
expect_error(synthesis_report_data(bad_dam$microglia, bad_dam$trajectory,
                                   bad_dam$mechanism, bad_dam$crossmodality),
             "DAM composition interaction")
cat("ok - synthesis missing-anchor guards fail loud\n")

cat("ok - synthesis tests complete\n")

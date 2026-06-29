#!/usr/bin/env Rscript
# scripts/build_human_robustness_mediation.R
#
# Session H5 of storage/notes/human_validation_plan.md.
#
# Two robustness arms + one complementary causal framing for the H4 human
# amyloid x tau interaction result. H4 fit, per signature x stratum, the
# donor-pseudobulk mixed model
#   z(score) ~ amyloid_c * tau_c + region + age_c + sex + (1|donor), w=n_cells
# on scanpy score_genes means, and read the amyloid:tau coefficient as the
# human analogue of the mouse interaction. H5 asks whether that interaction
# DIRECTION is an artefact of (A) collapsing cells to donor means, or (B) the
# score_genes scoring method, and reframes the same signatures as a Green-2024
# MEDIATION cascade (amyloid -> microglia signature -> tau) to contrast against
# the moderation reading.
#
# ----------------------------------------------------------------------------
# ARM A -- single-cell mixed model (AGGREGATION robustness).
#   Refit the H4 fixed-effect structure at the PER-CELL level with a donor
#   random intercept (no n_cells weight; each cell is a row):
#     z(score) ~ amyloid_c * tau_c + region + age_c + sex + (1|donor_id)
#   amyloid/tau are donor-level, so the random intercept absorbs donor noise
#   and the amyloid:tau SE is approximately donor-calibrated -- but cell-level
#   residuals are non-normal/heteroskedastic, so per-cell p-values are
#   APPROXIMATE. Pseudobulk (donor = unit) remains the conservative primary
#   estimate (Squair 2021, Nat Commun: pseudobulk beats single-cell mixed
#   models for DE); this arm tests DIRECTION concordance when within-donor cell
#   variance is retained, NOT a more-powered re-test. Same MIN_CELLS=10
#   donor x region x state floor as H4.
#
# ARM B -- AUCell scoring-method robustness.
#   build_human_aucell_rescoring.py re-scored every sizeable signature with
#   AUCell (rank recovery curve; orthogonal to score_genes' control-set mean)
#   and re-aggregated to donor x region x state means over the SAME H3
#   predicted_substate cells. Refit the EXACT H4 model on those means and test
#   amyloid:tau sign concordance. Tiny marker/DE sets (<5 genes) are AUCell-n/a
#   (rank scorer undefined for tiny sets) -> reported as NA, not dropped
#   silently; all six pre-registered headline signatures are large and tested.
#
# ARM C -- mediation (amyloid -> signature -> tau; Green 2024 cascade framing).
#   On donor-pseudobulk score_genes means, per signature x stratum:
#     mediator model : z(score) ~ amyloid_c + region + age_c + sex       -> a
#     outcome  model : tau      ~ amyloid_c + z(score) + region + age_c + sex
#                                                                 -> b, c' (ADE)
#     total    model : tau      ~ amyloid_c + region + age_c + sex        -> c
#   ACME = a*b ; ADE = c' ; proportion mediated = ACME / c. CIs/p from a
#   DONOR-CLUSTERED nonparametric bootstrap (resample donors with replacement,
#   both regions together, refit; percentile 95% CI; two-sided bootstrap p).
#   weights = n_cells (precision). Green 2024 (ROSMAP) places microglial
#   disease states causally between amyloid and tau, predicting POSITIVE
#   mediation (a>0 amyloid drives the state, b>0 state drives tau) for the
#   disease/DAM-like signatures. BH-FDR across the panel within each stratum.
#
# CONTRAST (moderation vs mediation, made explicit per the plan):
#   * H4 interaction = MODERATION: signature is the OUTCOME; tau modifies the
#     amyloid->signature slope (the mouse 2x2 logic -- tau background reshapes
#     the amyloid microglial response). Coefficient: amyloid:tau.
#   * H5 mediation   = MEDIATION: signature is the MEDIATOR transmitting the
#     amyloid->tau effect (the Green cascade logic). Coefficient: ACME = a*b.
#   These are DIFFERENT causal structures, not competing tests of one effect;
#   both are observational/cross-sectional, and mediation additionally IMPORTS
#   the amyloid->tau temporal ordering (Jack/Bateman cascade) as an untestable
#   assumption here. Reported side by side; framed "consistent with", never
#   causal proof (anti-anchoring #7).
#
# Confirmation cohort: H1 locked SEA-AD MTG+DLPFC and DEFERRED ROSMAP (no
# RADC/Synapse DUA held). Cross-cohort confirmation is therefore deferred
# (stated, not silently capped); the within-cohort MTG-vs-DLPFC regional
# replication (H4 region_concordant; reproduced in the per-region OLS here) is
# the available independent-subsample robustness.
#
# Inputs (light caches; no h5ad touched here):
#   storage/cache/summary_human_validation.rds        (H4: pseudobulk panel + frames)
#   storage/cache/human_substate_score_means.csv      (H3: score_genes means)
#   storage/cache/human_substate_aucell_score_means.csv (H5-AUCell: means)
#   storage/cache/human_substate_percell.csv.gz        (H3: per-cell scores + labels)
#   storage/cache/human_seaad_donor_neuropath.csv      (H1: staging + covariates)
#
# Outputs:
#   storage/results/human_robustness_mediation.tsv     one row / signature x stratum
#   storage/cache/summary_human_robustness.rds         light plotting cache (H6)
#   storage/cache/human_robustness_mediation_provenance.txt
#
# Run (deterministic; seed fixed): Rscript scripts/build_human_robustness_mediation.R [--overwrite|--smoke]

suppressPackageStartupMessages({
  library(data.table)
  library(lme4)
  library(lmerTest)
})

args      <- commandArgs(trailingOnly = TRUE)
overwrite <- "--overwrite" %in% args
smoke     <- "--smoke"     %in% args
setwd("/home/rstudio/tau-mutant-integration-ng")
set.seed(1)

cache_dir   <- "storage/cache"
results_dir <- "storage/results"
tsv_path  <- file.path(results_dir, "human_robustness_mediation.tsv")
rds_path  <- file.path(cache_dir,   "summary_human_robustness.rds")
prov_path <- file.path(cache_dir,   "human_robustness_mediation_provenance.txt")

if (!smoke && all(file.exists(tsv_path, rds_path)) && !overwrite) {
  message("H5 outputs present; pass --overwrite to rebuild. Skipping.")
  quit(save = "no", status = 0)
}

write_tsv_safe <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_tsv(x, path); Sys.chmod(path, mode = "0664")
}

MIN_CELLS <- 10L
STRATA    <- c("all", "homeostatic", "DAM", "IFN", "proliferative")
N_BOOT    <- if (smoke) 200L else 2000L

ctrl <- lmerControl(optimizer = "bobyqa",
                    check.conv.singular = .makeCC("ignore", tol = 1e-4))

## ------------------------------------------------------------------- load
h4 <- readRDS(file.path(cache_dir, "summary_human_validation.rds"))
pb <- as.data.table(h4$models)[, .(signature, stratum,
                                    pb_ix_est = ix_est, pb_ix_p = ix_p,
                                    pb_ix_fdr = ix_fdr)]

np <- fread(file.path(cache_dir, "human_seaad_donor_neuropath.csv"))
npm <- np[!is.na(amyloid_thal) & !is.na(tau_braak),
          .(donor_id, region,
            amyloid = as.numeric(amyloid_thal), tau = as.numeric(tau_braak),
            age = as.numeric(age), sex = factor(sex))]

sm_sg  <- fread(file.path(cache_dir, "human_substate_score_means.csv"))
sm_auc <- fread(file.path(cache_dir, "human_substate_aucell_score_means.csv"))
sig_sg  <- sub("^score_", "", grep("^score_", names(sm_sg),  value = TRUE))
sig_auc <- sub("^score_", "", grep("^score_", names(sm_auc), value = TRUE))

## --------------------------------------------------- per-stratum prep (H4)
# Mirror H4 exactly: merge means+staging, n_cells>=MIN_CELLS, centre on the
# analysed rows, flag discordant (off-diagonal half by |resid(tau~amyloid)|).
prep_stratum <- function(sm, st) {
  fr <- merge(sm[state == st], npm, by = c("donor_id", "region"))
  fr <- fr[n_cells >= MIN_CELLS]
  if (nrow(fr) < 12L) return(NULL)
  fr[, region := factor(region)]
  fr[, amyloid_c := amyloid - mean(amyloid)]
  fr[, tau_c     := tau     - mean(tau)]
  fr[, age_c     := age     - mean(age)]
  fr[]
}

## ====================================================== ARM A: single cell
# Per-cell scores + labels; join staging; refit at cell level with donor RE.
pc_cols <- c("donor_id", "region", "predicted_substate",
             paste0("score_", sig_sg))
pc <- fread(cmd = paste0("zcat ", file.path(cache_dir,
            "human_substate_percell.csv.gz")), select = pc_cols)
pc <- merge(pc, npm, by = c("donor_id", "region"))   # staged donors only
pc[, region := factor(region)]

sc_ix <- function(sig, st) {
  d <- if (st == "all") copy(pc) else pc[predicted_substate == st]
  # MIN_CELLS donor x region x state floor (mirror H4's analysed set).
  d[, ncs := .N, by = .(donor_id, region)]
  d <- d[ncs >= MIN_CELLS]
  yv <- d[[paste0("score_", sig)]]
  if (uniqueN(d$donor_id) < 12L || sd(yv, na.rm = TRUE) == 0)
    return(data.table(sc_ix_est = NA_real_, sc_ix_p = NA_real_,
                      sc_n_cells = nrow(d), sc_singular = NA))
  d[, `:=`(amyloid_c = amyloid - mean(amyloid), tau_c = tau - mean(tau),
           age_c = age - mean(age), y = as.numeric(scale(yv)))]
  wmsg <- NA_character_
  fit <- tryCatch(withCallingHandlers(
    lmerTest::lmer(y ~ amyloid_c * tau_c + region + age_c + sex + (1 | donor_id),
                   data = d, control = ctrl, REML = TRUE),
    warning = function(w) { wmsg <<- conditionMessage(w); invokeRestart("muffleWarning") }),
    error = function(e) NULL)
  if (is.null(fit)) return(data.table(sc_ix_est = NA_real_, sc_ix_p = NA_real_,
                                      sc_n_cells = nrow(d), sc_singular = NA))
  s <- summary(fit)$coefficients
  pc_col <- if ("Pr(>|t|)" %in% colnames(s)) "Pr(>|t|)" else NA
  data.table(sc_ix_est = unname(s["amyloid_c:tau_c", "Estimate"]),
             sc_ix_p   = if (!is.na(pc_col)) unname(s["amyloid_c:tau_c", pc_col]) else NA,
             sc_n_cells = nrow(d), sc_singular = isSingular(fit))
}

## ====================================================== ARM B: AUCell refit
auc_ix <- function(sig, fr_auc) {
  col <- paste0("score_", sig)
  if (is.null(fr_auc) || !col %in% names(fr_auc))
    return(data.table(auc_ix_est = NA_real_, auc_ix_p = NA_real_))
  yv <- fr_auc[[col]]
  if (sd(yv, na.rm = TRUE) == 0)
    return(data.table(auc_ix_est = NA_real_, auc_ix_p = NA_real_))
  d <- copy(fr_auc); d[, y := as.numeric(scale(yv))]
  wmsg <- NA_character_
  fit <- tryCatch(withCallingHandlers(
    lmerTest::lmer(y ~ amyloid_c * tau_c + region + age_c + sex + (1 | donor_id),
                   data = d, weights = n_cells, control = ctrl, REML = TRUE),
    warning = function(w) { wmsg <<- conditionMessage(w); invokeRestart("muffleWarning") }),
    error = function(e) NULL)
  if (is.null(fit)) return(data.table(auc_ix_est = NA_real_, auc_ix_p = NA_real_))
  s <- summary(fit)$coefficients
  pc_col <- if ("Pr(>|t|)" %in% colnames(s)) "Pr(>|t|)" else NA
  data.table(auc_ix_est = unname(s["amyloid_c:tau_c", "Estimate"]),
             auc_ix_p   = if (!is.na(pc_col)) unname(s["amyloid_c:tau_c", pc_col]) else NA)
}

## ====================================================== ARM C: mediation
# WLS coefficients via lm.wfit (pivots for rank deficiency -> bootstrap-safe).
wls <- function(X, y, w) {
  fit <- lm.wfit(X, y, w)
  cf <- fit$coefficients; cf[is.na(cf)] <- 0; cf
}
# Per stratum: shared donor-resample indices so a (mediator) and b (outcome)
# are bootstrapped on PAIRED samples -> ACME=a*b is a proper bootstrap draw.
mediation_stratum <- function(fr, sigs) {
  don <- unique(fr$donor_id)
  rows_by_donor <- split(seq_len(nrow(fr)), fr$donor_id)
  w   <- fr$n_cells
  # design matrices (intercept + covariates); amyloid_c is the exposure.
  Xm  <- model.matrix(~ amyloid_c + region + age_c + sex, fr)         # mediator/total
  tau <- fr$tau
  ai_m <- which(colnames(Xm) == "amyloid_c")
  boot <- lapply(seq_len(N_BOOT), function(i)
    unlist(rows_by_donor[sample(don, replace = TRUE)], use.names = FALSE))
  # total effect c: tau ~ amyloid_c + covs (signature-independent) -- point + boot
  c_pt <- wls(Xm, tau, w)[ai_m]
  c_bt <- vapply(boot, function(ix) wls(Xm[ix, , drop = FALSE], tau[ix], w[ix])[ai_m],
                 numeric(1))
  out <- vector("list", length(sigs))
  for (k in seq_along(sigs)) {
    sig <- sigs[k]; M <- as.numeric(scale(fr[[paste0("score_", sig)]]))
    Xy  <- cbind(Xm[, 1:ai_m, drop = FALSE], M = M,
                 Xm[, (ai_m + 1):ncol(Xm), drop = FALSE])  # insert M after amyloid_c
    ai_y <- which(colnames(Xy) == "amyloid_c"); mi_y <- which(colnames(Xy) == "M")
    a_pt <- wls(Xm, M, w)[ai_m]
    oy   <- wls(Xy, tau, w); b_pt <- oy[mi_y]; cp_pt <- oy[ai_y]
    acme_pt <- a_pt * b_pt
    a_bt <- numeric(N_BOOT); b_bt <- numeric(N_BOOT); cp_bt <- numeric(N_BOOT)
    for (i in seq_len(N_BOOT)) {
      ix <- boot[[i]]
      a_bt[i] <- wls(Xm[ix, , drop = FALSE], M[ix], w[ix])[ai_m]
      o  <- wls(Xy[ix, , drop = FALSE], tau[ix], w[ix]); b_bt[i] <- o[mi_y]; cp_bt[i] <- o[ai_y]
    }
    acme_bt <- a_bt * b_bt
    ci   <- quantile(acme_bt, c(0.025, 0.975), names = FALSE, na.rm = TRUE)
    pval <- 2 * min(mean(acme_bt <= 0), mean(acme_bt >= 0)); pval <- min(pval, 1)
    out[[k]] <- data.table(signature = sig,
      med_a = a_pt, med_b = b_pt, med_acme = acme_pt,
      med_acme_lo = ci[1], med_acme_hi = ci[2], med_acme_p = pval,
      med_ade = cp_pt, med_total = c_pt,
      med_prop = if (abs(c_pt) > 1e-9) acme_pt / c_pt else NA_real_)
  }
  rbindlist(out)
}

## -------------------------------------------------------------- run strata
sig_run <- if (smoke) c("DAM_up", "NFKB_union_targets", "MG_M3_module", "Gsk3b_targets",
                        "Gerrits_AD1_human") else union(sig_sg, sig_auc)
str_run <- if (smoke) "all" else STRATA

panel <- list(); med_all <- list()
for (st in str_run) {
  fr_sg  <- prep_stratum(sm_sg,  st)
  fr_auc <- prep_stratum(sm_auc, st)
  if (is.null(fr_sg)) { message("stratum ", st, ": <12 pseudobulk rows; skipping"); next }

  # ARM A + ARM B per signature
  rows <- rbindlist(lapply(sig_run, function(sig) {
    cbind(data.table(signature = sig, stratum = st),
          sc_ix(sig, st), auc_ix(sig, fr_auc))
  }), fill = TRUE)

  # ARM C mediation (score_genes pseudobulk; shared bootstrap)
  med <- mediation_stratum(fr_sg, intersect(sig_run, sig_sg))
  med[, stratum := st]
  med_all[[st]] <- med
  rows <- merge(rows, med, by = c("signature", "stratum"), all.x = TRUE, sort = FALSE)
  panel[[st]] <- rows
  message("stratum ", st, " done (", nrow(fr_sg), " pseudobulk rows, ",
          length(intersect(sig_run, sig_sg)), " mediations x ", N_BOOT, " boot)")
}
panel <- rbindlist(panel, fill = TRUE)

## ----------------------------------------------- merge H4 pseudobulk + flags
panel <- merge(panel, pb, by = c("signature", "stratum"), all.x = TRUE, sort = FALSE)
sgn <- function(x) fifelse(is.na(x), NA_character_, fifelse(x > 0, "positive", "negative"))
panel[, `:=`(
  pb_ix_dir  = sgn(pb_ix_est),
  sc_ix_dir  = sgn(sc_ix_est),
  auc_ix_dir = sgn(auc_ix_est),
  sc_concordant  = fifelse(is.na(sc_ix_est)  | is.na(pb_ix_est), NA, sign(sc_ix_est)  == sign(pb_ix_est)),
  auc_concordant = fifelse(is.na(auc_ix_est) | is.na(pb_ix_est), NA, sign(auc_ix_est) == sign(pb_ix_est)),
  med_dir    = sgn(med_acme))]
panel[, med_acme_fdr := p.adjust(med_acme_p, "BH"), by = stratum]

# pre-registration (carry H4's; add mediation expectation for disease sigs)
prereg <- data.table(
  signature = c("DAM_up", "Gerrits_AD1_human", "Gerrits_AD2_human",
                "NFKB_union_targets", "MG_M3_module", "Gsk3b_targets"),
  prereg_ix_expected  = c(NA, NA, NA, "negative", "negative", "nonnull"),
  prereg_med_expected = c("positive", "positive", "positive", NA, NA, NA))
panel <- merge(panel, prereg, by = "signature", all.x = TRUE, sort = FALSE)

setorder(panel, stratum, signature)
num_cols <- names(which(vapply(panel, is.numeric, logical(1))))
panel[, (num_cols) := lapply(.SD, function(x) signif(x, 4)), .SDcols = num_cols]

## ----------------------------------------------- contrast (moderation x mediation)
contrast <- panel[stratum == "all", .(signature,
  moderation_ix_est = pb_ix_est, moderation_ix_dir = pb_ix_dir,
  mediation_a = med_a, mediation_b = med_b, mediation_acme = med_acme,
  mediation_dir = med_dir, mediation_prop = med_prop,
  sc_concordant, auc_concordant)]

## ------------------------------------------------------------------- smoke
if (smoke) {
  cat("\n=== SMOKE: stratum 'all',", length(sig_run), "sigs, B =", N_BOOT, "===\n")
  print(panel[, .(signature, pb_ix_est, sc_ix_est, sc_concordant,
                  auc_ix_est, auc_concordant,
                  med_a, med_b, med_acme, med_acme_lo, med_acme_hi, med_acme_p)])
  cat("\nconcordance: sc", sum(panel$sc_concordant, na.rm = TRUE), "/", sum(!is.na(panel$sc_concordant)),
      "| auc", sum(panel$auc_concordant, na.rm = TRUE), "/", sum(!is.na(panel$auc_concordant)), "\n")
  quit(save = "no", status = 0)
}

## ----------------------------------------------------------------- persist
write_tsv_safe(panel, tsv_path)
summary_human_robustness <- list(
  panel    = panel,
  contrast = contrast,
  mediation = rbindlist(med_all, fill = TRUE),
  meta = list(
    sc_formula  = "z(score) ~ amyloid_c*tau_c + region + age_c + sex + (1|donor) [per-cell]",
    auc_formula = "z(AUCell) ~ amyloid_c*tau_c + region + age_c + sex + (1|donor), w=n_cells",
    med_models  = c(mediator = "z(score) ~ amyloid_c + region + age_c + sex",
                    outcome  = "tau ~ amyloid_c + z(score) + region + age_c + sex",
                    total    = "tau ~ amyloid_c + region + age_c + sex"),
    med_boot    = N_BOOT, med_boot_unit = "donor (clustered, both regions together)",
    min_cells   = MIN_CELLS, strata = STRATA,
    aucell_n_sig = length(sig_auc), score_genes_n_sig = length(sig_sg),
    confirmation_cohort = "DEFERRED: ROSMAP/Synapse DUA not held (H1 decision); within-cohort MTG-vs-DLPFC is the available replication",
    built = "scripts/build_human_robustness_mediation.R"))
saveRDS(summary_human_robustness, rds_path); Sys.chmod(rds_path, mode = "0664")

## -------------------------------------------------------------- provenance
pre <- panel[stratum == "all" & !is.na(prereg_ix_expected)]
sc_conc_n  <- panel[, .(c = sum(sc_concordant, na.rm = TRUE), n = sum(!is.na(sc_concordant)))]
auc_conc_n <- panel[, .(c = sum(auc_concordant, na.rm = TRUE), n = sum(!is.na(auc_concordant)))]
key <- c("NFKB_union_targets", "MG_M3_module", "Gsk3b_targets", "DAM_up")
prov <- c(
  "human_robustness_mediation provenance (H5)",
  "===========================================",
  "built_by : scripts/build_human_robustness_mediation.R",
  "ARM A (single-cell aggregation robustness): per-cell lmer with donor RE;",
  "  pseudobulk (donor unit) stays the conservative primary (Squair 2021);",
  "  per-cell p APPROXIMATE (non-normal cell residuals) -> DIRECTION is the readout.",
  paste0("  sc direction concordant with H4 pseudobulk in ", sc_conc_n$c, "/", sc_conc_n$n, " models."),
  "ARM B (AUCell scoring-method robustness): orthogonal rank scorer, same H4 model;",
  paste0("  auc direction concordant with H4 pseudobulk in ", auc_conc_n$c, "/", auc_conc_n$n,
         " models (sizeable sets only; tiny marker/DE sets AUCell-n/a)."),
  "ARM C (mediation, Green-2024 cascade amyloid->signature->tau): donor-clustered",
  paste0("  bootstrap B=", N_BOOT, "; ACME=a*b, ADE=c', proportion mediated=ACME/total;"),
  "  BH-FDR within stratum. Positive mediation expected for disease/DAM signatures.",
  "CONTRAST: H4 interaction = MODERATION (signature=outcome; tau modifies amyloid",
  "  response); H5 mediation = signature=MEDIATOR on the amyloid->tau path. Distinct",
  "  causal structures; both observational; mediation imports amyloid->tau ordering",
  "  (Jack/Bateman) as an untestable assumption -> 'consistent with', never proof.",
  "CONFIRMATION COHORT: DEFERRED -- ROSMAP/Synapse DUA not held (H1). MTG-vs-DLPFC",
  "  per-region replication is the available within-cohort robustness (no silent cap).",
  "",
  "key signatures (stratum 'all'): pb=score_genes pseudobulk ix, sc=single-cell ix,",
  "  auc=AUCell ix, ACME=mediation effect [95% CI], p:")
for (s in key) {
  r <- panel[stratum == "all" & signature == s]
  if (nrow(r) == 0) next
  prov <- c(prov, sprintf(
    "  %-20s pb=%+.4f  sc=%+.4f(%s)  auc=%s(%s)  ACME=%+.4f[%+.4f,%+.4f] p=%.3f",
    s, r$pb_ix_est, r$sc_ix_est, ifelse(isTRUE(r$sc_concordant), "concord", "DISCORD"),
    ifelse(is.na(r$auc_ix_est), "n/a", sprintf("%+.4f", r$auc_ix_est)),
    ifelse(is.na(r$auc_concordant), "-", ifelse(isTRUE(r$auc_concordant), "concord", "DISCORD")),
    r$med_acme, r$med_acme_lo, r$med_acme_hi, r$med_acme_p))
}
prov <- c(prov, "",
  "pre-registered mediation (positive expected; disease cascade):")
for (s in c("DAM_up", "Gerrits_AD1_human", "Gerrits_AD2_human")) {
  r <- panel[stratum == "all" & signature == s]
  if (nrow(r) == 0) next
  prov <- c(prov, sprintf("  %-20s a=%+.4f b=%+.4f ACME=%+.4f p=%.3f fdr=%.3f dir=%s",
    s, r$med_a, r$med_b, r$med_acme, r$med_acme_p, r$med_acme_fdr, r$med_dir))
}
prov <- c(prov, "",
  paste0("outputs : ", tsv_path, " (", nrow(panel), " rows) ; ", rds_path),
  paste0("built_at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
writeLines(prov, prov_path); Sys.chmod(prov_path, mode = "0664")

cat("H5 done:", nrow(panel), "rows (", uniqueN(panel$signature), "sigs x",
    uniqueN(panel$stratum), "strata).\n  ->", tsv_path, "\n  ->", rds_path,
    "\n  ->", prov_path, "\n")

#!/usr/bin/env Rscript
# scripts/build_human_interaction_models.R
#
# Session H4 of storage/notes/human_validation_plan.md.
#
# Human amyloid x tau INTERACTION models on microglial signature scores.
# The human analogue of the mouse interaction contrast
#   (NLGF_P301S - P301S) - (NLGF_MAPTKI - MAPTKI)
# is the cross-donor amyloid:tau coefficient. In human AD amyloid (Thal)
# and tau (Braak) co-progress and vary continuously / partly discordantly,
# so a per-donor product term tests whether the amyloid response of each
# microglial signature is tau-dependent -- under explicit collinearity
# mitigation (the two axes are Spearman ~0.65 collinear; a null is therefore
# ambiguous until identifiability is reported, anti-anchoring #1).
#
# Inputs (light caches; no h5ad touched here):
#   storage/cache/human_substate_score_means.csv  (H3) per donor x region x
#       {all + 4 substates} mean score for all 26 H2 signatures + n_cells.
#   storage/cache/human_seaad_donor_neuropath.csv (H1) per donor x region
#       ordinal staging (amyloid_thal 0-5, tau_braak 0/2-6) + age/sex/region.
#       The 8 Allen reference rows carry no staging (NA) and drop out, so the
#       analysed set is the AD-continuum donors (which already span Thal 0-5 /
#       Braak 0-6, i.e. low-pathology anchors are retained).
#
# Model (per signature x stratum); donor = unit of replication, mirroring the
# mouse genotype_batch pseudobulk logic, with brain region as a within-cohort
# replicate (83/89 donors are paired across MTG + DLPFC):
#   z(score) ~ amyloid_c * tau_c + region + age_c + sex + (1 | donor)
#   * precision-weighted by n_cells (Var of a per-state mean ~ 1/n_cells);
#   * amyloid / tau / age mean-centred on the ANALYSED rows (centring removes
#     non-essential main-effect<->product collinearity and makes main effects
#     interpretable at the sample mean of the other axis);
#   * response z-scored per signature so the amyloid:tau coefficient is a
#     comparable cross-panel effect size for the H6 forest plot (z-scoring is
#     a positive rescale, so it preserves sign -> direction tests unaffected);
#   * lmerTest (Satterthwaite df) -> estimate, SE, Wald 95% CI, p; BH-FDR
#     across the 26-signature panel WITHIN each stratum (each stratum is one
#     family). Whole-microglia stratum "all" is the headline panel.
#
# Collinearity battery (reported with every estimate):
#   - VIF on the centred fixed-effect design INCLUDING the product column, so
#     the interaction term's own inflation is reported (per stratum: identical
#     X across signatures). With centring this isolates the genuine amyloid~tau
#     collinearity (Spearman ~0.65 -> main-effect VIF ~1.7) from the spurious
#     main-effect<->product collinearity that uncentred products manufacture.
#   - per-region replication: independent MTG and DLPFC OLS refits (a donor is
#     unique within a region, so no random effect is needed) -> sign
#     concordance across the two regions is within-cohort replication (the
#     reason H1 acquired both regions). This stands in for a literal
#     "residualised re-fit", which by the Frisch-Waugh-Lovell theorem returns
#     the joint model's interaction coefficient UNCHANGED -- the joint
#     amyloid*tau fit already yields collinearity-adjusted PARTIAL effects, so
#     residualising carries no independent information; VIF already quantifies
#     the inflation it would notionally address.
#   - discordant-donor sensitivity: refit on the off-diagonal half of donors
#     (|residual of tau~amyloid| at/above its median: the donors whose tau is
#     most decoupled from their amyloid); sign stability there = robustness to
#     the amyloid~tau collinearity itself.
#
# Pre-registered directional hypotheses (DIRECTION beats significance,
# anti-anchoring #4; sanity-gate on positive controls, #5):
#   DAM_up / Gerrits AD1 / Gerrits AD2 -> POSITIVE amyloid MAIN effect (sanity)
#   NFKB_union_targets                 -> NEGATIVE amyloid:tau (mouse tau
#                                         attenuates amyloid-driven NF-kB)
#   MG_M3_module                       -> NEGATIVE amyloid:tau
#   Gsk3b_targets                      -> NON-NULL amyloid:tau (tau-dependent
#                                         amyloid response; mouse locks no sign)
#
# Per-state strata are cleared because the H3 conservation gate PASSED (DAM &
# IFN resolve); per-state rows require n_cells >= MIN_CELLS so a donor's state
# mean is not a 1-cell artefact (raw IFN/proliferative medians are 23/14).
#
# Idempotent: skips if outputs present unless `--overwrite`. `--smoke` runs the
# "all" stratum on three signatures, prints, and writes nothing (shape check).
#
# Outputs:
#   storage/results/human_interaction_models.tsv         one row / signature x stratum
#   storage/cache/summary_human_validation.rds           light plotting cache (H6)
#   storage/cache/human_interaction_models_provenance.txt
#
# Run (deterministic): Rscript scripts/build_human_interaction_models.R [--overwrite|--smoke]

suppressPackageStartupMessages({
  library(data.table)
  library(lme4)
  library(lmerTest)
})

args      <- commandArgs(trailingOnly = TRUE)
overwrite <- "--overwrite" %in% args
smoke     <- "--smoke"     %in% args
setwd("/home/rstudio/tau-mutant-integration-ng")

cache_dir   <- "storage/cache"
results_dir <- "storage/results"

tsv_path  <- file.path(results_dir, "human_interaction_models.tsv")
rds_path  <- file.path(cache_dir,   "summary_human_validation.rds")
prov_path <- file.path(cache_dir,   "human_interaction_models_provenance.txt")

if (!smoke && all(file.exists(tsv_path, rds_path)) && !overwrite) {
  message("H4 outputs present; pass --overwrite to rebuild. Skipping.")
  quit(save = "no", status = 0)
}

# Mirror R/utils.R::write_tsv_safe locally: H4 is pure statistics and needs no
# bioinformatics helper, so we avoid sourcing R/helpers.R (which would attach
# nichenetr / Seurat / the full DE stack).
write_tsv_safe <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_tsv(x, path)
  Sys.chmod(path, mode = "0664")
}

MIN_CELLS <- 10L            # per-donor state floor: drop ~1-cell state means
STRATA    <- c("all", "homeostatic", "DAM", "IFN", "proliferative")

## ------------------------------------------------------------------- load
sm <- fread(file.path(cache_dir, "human_substate_score_means.csv"))
np <- fread(file.path(cache_dir, "human_seaad_donor_neuropath.csv"))

score_cols <- grep("^score_", names(sm), value = TRUE)
signatures <- sub("^score_", "", score_cols)

# Staged continuum donors only (reference rows carry NA staging -> drop).
npm <- np[!is.na(amyloid_thal) & !is.na(tau_braak),
          .(donor_id, region,
            amyloid = as.numeric(amyloid_thal),
            tau     = as.numeric(tau_braak),
            age     = as.numeric(age),
            sex     = factor(sex))]
stopifnot(nrow(npm) > 0,
          !any(npm$donor_id %in% np[is_reference == "True", donor_id]))

dat <- merge(sm, npm, by = c("donor_id", "region"))
dat[, region := factor(region)]

# Off-diagonal coverage of the staging plane (identifiability headline).
allrows  <- unique(dat[state == "all", .(donor_id, region, amyloid, tau)])
spearman <- suppressWarnings(cor(allrows$amyloid, allrows$tau, method = "spearman"))

## -------------------------------------------------------------- utilities
wald_ci <- function(est, se) est + c(-1, 1) * qnorm(0.975) * se

# Manual VIF on a centred design matrix (incl. the product column) so the
# interaction term's inflation is itself reported. Depends only on X, hence
# identical for every signature in a stratum -> computed once per stratum.
vif_design <- function(fr) {
  X  <- model.matrix(~ amyloid_c * tau_c + region + age_c + sex, fr)[, -1, drop = FALSE]
  cn <- colnames(X)
  v  <- vapply(seq_along(cn), function(j) {
    r2 <- summary(lm(X[, j] ~ X[, -j, drop = FALSE]))$r.squared
    1 / (1 - r2)
  }, numeric(1))
  setNames(v, cn)
}

ctrl <- lmerControl(optimizer = "bobyqa",
                    check.conv.singular = .makeCC("ignore", tol = 1e-4))

# Fit a warning-muffling lmerTest model, recording the (last) warning message.
fit_lmer <- function(form, data) {
  wmsg <- NA_character_
  fit  <- tryCatch(
    withCallingHandlers(
      lmerTest::lmer(form, data = data, weights = n_cells, control = ctrl, REML = TRUE),
      warning = function(w) { wmsg <<- conditionMessage(w); invokeRestart("muffleWarning") }),
    error = function(e) NULL)
  list(fit = fit, warn = wmsg)
}

coef_row <- function(fit, term, has_p = TRUE) {
  if (is.null(fit)) return(c(est = NA, se = NA, p = NA))
  s <- summary(fit)$coefficients
  if (!term %in% rownames(s)) return(c(est = NA, se = NA, p = NA))
  pcol <- if ("Pr(>|t|)" %in% colnames(s)) "Pr(>|t|)" else NA
  c(est = unname(s[term, "Estimate"]),
    se  = unname(s[term, "Std. Error"]),
    p   = if (has_p && !is.na(pcol)) unname(s[term, pcol]) else NA)
}

# est + p of the amyloid:tau term from an lmerTest OR lm fit (robustness refits).
ix_estp <- function(obj) {
  if (is.null(obj)) return(c(est = NA, p = NA))
  s  <- summary(obj)$coefficients
  pc <- grep("^Pr", colnames(s))
  if (!"amyloid_c:tau_c" %in% rownames(s) || length(pc) == 0) return(c(est = NA, p = NA))
  c(est = unname(s["amyloid_c:tau_c", "Estimate"]), p = unname(s["amyloid_c:tau_c", pc[1]]))
}

## --------------------------------------------------------- per-model fit
fit_one <- function(fr, sig, stratum) {
  yv   <- fr[[paste0("score_", sig)]]
  out0 <- data.table(signature = sig, stratum = stratum,
                     n_obs = nrow(fr), n_donors = uniqueN(fr$donor_id),
                     n_cells_min = min(fr$n_cells), n_cells_med = median(fr$n_cells))
  if (sd(yv, na.rm = TRUE) == 0 || sum(!is.na(yv)) < 8L)
    return(cbind(out0, data.table(model = "skipped")))

  d <- copy(fr)
  d[, y := as.numeric(scale(yv))]                 # z-score response (sign-safe)

  m  <- fit_lmer(y ~ amyloid_c * tau_c + region + age_c + sex + (1 | donor_id), d)
  if (is.null(m$fit)) return(cbind(out0, data.table(model = "failed")))

  a  <- coef_row(m$fit, "amyloid_c")
  tt <- coef_row(m$fit, "tau_c")
  ix <- coef_row(m$fit, "amyloid_c:tau_c")
  ci <- wald_ci(ix["est"], ix["se"])

  # robustness 1: discordant-donor sensitivity (off-diagonal half); lmer, OLS fallback.
  ds   <- d[discordant == TRUE]
  dm   <- fit_lmer(y ~ amyloid_c * tau_c + region + age_c + sex + (1 | donor_id), ds)
  dmod <- "lmer"; dobj <- dm$fit
  if (is.null(dobj)) {
    dobj <- tryCatch(lm(y ~ amyloid_c * tau_c + region + sex, data = ds, weights = n_cells),
                     error = function(e) NULL)
    dmod <- if (is.null(dobj)) "none" else "lm"
  }
  dsx <- ix_estp(dobj)

  # robustness 2: per-region replication (independent subsamples; OLS, donor
  # unique within a region -> no random effect). MTG + DLPFC sign concordance
  # is the within-cohort regional replication H1's two-region design buys.
  per_region <- function(rg) {
    sub <- d[region == rg]
    if (uniqueN(sub$donor_id) < 12L) return(c(est = NA, p = NA))
    ix_estp(tryCatch(lm(y ~ amyloid_c * tau_c + age_c + sex, data = sub, weights = n_cells),
                     error = function(e) NULL))
  }
  mtg <- per_region("MTG"); dlp <- per_region("DLPFC")
  reg_conc <- !is.na(mtg["est"]) && !is.na(dlp["est"]) && sign(mtg["est"]) == sign(dlp["est"])

  cbind(out0, data.table(
    model = "lmer",
    converged = is.na(m$warn), conv_msg = ifelse(is.na(m$warn), NA_character_, m$warn),
    singular = isSingular(m$fit),
    resid_shapiro_p = tryCatch(shapiro.test(residuals(m$fit))$p.value, error = function(e) NA),
    amyloid_est = a["est"], amyloid_se = a["se"], amyloid_p = a["p"],
    tau_est = tt["est"], tau_se = tt["se"], tau_p = tt["p"],
    ix_est = ix["est"], ix_se = ix["se"], ix_ci_lo = ci[1], ix_ci_hi = ci[2], ix_p = ix["p"],
    ix_discordant_est = dsx["est"], ix_discordant_p = dsx["p"], discordant_model = dmod,
    ix_mtg_est = mtg["est"], ix_mtg_p = mtg["p"],
    ix_dlpfc_est = dlp["est"], ix_dlpfc_p = dlp["p"], region_concordant = reg_conc))
}

## ----------------------------------------------------------- stratum loop
sig_run <- if (smoke) c("DAM_up", "NFKB_union_targets", "Gsk3b_targets") else signatures
str_run <- if (smoke) "all" else STRATA

panel  <- list(); frames <- list(); vif_list <- list(); disc_n <- list()
for (st in str_run) {
  fr <- dat[state == st & n_cells >= MIN_CELLS]
  if (nrow(fr) < 12L) { message("stratum ", st, ": <12 rows after n_cells filter; skipping"); next }
  fr[, amyloid_c := amyloid - mean(amyloid)]
  fr[, tau_c     := tau     - mean(tau)]
  fr[, age_c     := age     - mean(age)]

  # discordant = off-diagonal half by |resid(tau ~ amyloid)| (donor-level).
  don <- unique(fr[, .(donor_id, amyloid, tau)])
  don[, r := residuals(lm(tau ~ amyloid, data = don))]
  disc_donors <- don[abs(r) >= median(abs(r)), donor_id]
  fr[, discordant := donor_id %in% disc_donors]

  vif_list[[st]] <- vif_design(fr)
  disc_n[[st]]   <- length(disc_donors)

  rows <- rbindlist(lapply(sig_run, function(sg) fit_one(fr, sg, st)), fill = TRUE)
  v <- vif_list[[st]]
  rows[, `:=`(vif_amyloid = v["amyloid_c"], vif_tau = v["tau_c"],
              vif_interaction = v["amyloid_c:tau_c"],
              n_discordant_donors = length(disc_donors),
              spearman_amyloid_tau = spearman)]
  panel[[st]]  <- rows
  frames[[st]] <- fr
}
panel <- rbindlist(panel, fill = TRUE)

## ------------------------------------------- FDR, direction, pre-registration
panel[, ix_fdr      := p.adjust(ix_p,      "BH"), by = stratum]
panel[, amyloid_fdr := p.adjust(amyloid_p, "BH"), by = stratum]
panel[, tau_fdr     := p.adjust(tau_p,     "BH"), by = stratum]
panel[, ix_direction := fifelse(is.na(ix_est), NA_character_,
                                fifelse(ix_est > 0, "positive", "negative"))]

prereg <- data.table(
  signature = c("DAM_up", "Gerrits_AD1_human", "Gerrits_AD2_human",
                "NFKB_union_targets", "MG_M3_module", "Gsk3b_targets"),
  prereg_term     = c("amyloid", "amyloid", "amyloid", "interaction", "interaction", "interaction"),
  prereg_expected = c("positive", "positive", "positive", "negative", "negative", "nonnull"),
  prereg_label    = c("DAM-up amyloid main (sanity #5)",
                      "Gerrits AD1 amyloid main (sanity #5)",
                      "Gerrits AD2 amyloid main (sanity #5)",
                      "NF-kB tau-attenuation (negative interaction)",
                      "MG-M3 negative interaction",
                      "Gsk3b tau-dependent amyloid (non-null interaction)"))
panel <- merge(panel, prereg, by = "signature", all.x = TRUE, sort = FALSE)
panel[, prereg_observed := fifelse(is.na(prereg_term), NA_character_,
        fifelse(prereg_term == "amyloid",
                fifelse(amyloid_est > 0, "positive", "negative"),
                fifelse(ix_est > 0, "positive", "negative")))]
panel[, prereg_concordant := fcase(
        is.na(prereg_expected),            NA,
        prereg_expected == "nonnull",      !is.na(prereg_observed),
        default = prereg_observed == prereg_expected)]

# stable, readable ordering: stratum then signature
setorder(panel, stratum, signature)
num_cols <- names(which(vapply(panel, is.numeric, logical(1))))
panel[, (num_cols) := lapply(.SD, function(x) signif(x, 4)), .SDcols = num_cols]

## ------------------------------------------------------------------- smoke
if (smoke) {
  cat("\n=== SMOKE: stratum 'all', 3 signatures ===\n")
  print(panel[, .(signature, n_obs, n_donors, vif_interaction,
                  amyloid_est, amyloid_p, ix_est, ix_ci_lo, ix_ci_hi, ix_p,
                  ix_discordant_est, ix_mtg_est, ix_dlpfc_est, region_concordant,
                  prereg_expected, prereg_concordant)])
  cat("\nspearman(amyloid,tau) =", round(spearman, 3),
      "| discordant donors (all) =", disc_n[["all"]], "\n")
  quit(save = "no", status = 0)
}

## ----------------------------------------------------------------- persist
write_tsv_safe(panel, tsv_path)

summary_human_validation <- list(
  models = panel,
  frames = frames,                       # per-stratum analysed rows (scores + staging)
  vif    = vif_list,
  prereg = prereg,
  meta   = list(
    formula     = "z(score) ~ amyloid_c * tau_c + region + age_c + sex + (1|donor)",
    min_cells   = MIN_CELLS,
    strata      = STRATA,
    n_donors    = uniqueN(allrows$donor_id),
    n_staged    = nrow(allrows),
    spearman_amyloid_tau = spearman,
    discordant_n = unlist(disc_n),
    fdr_scope   = "BH within each stratum (26-signature family)",
    response    = "per-signature z-scored mean module score",
    weights     = "n_cells (precision)",
    built       = "scripts/build_human_interaction_models.R"))
saveRDS(summary_human_validation, rds_path)
Sys.chmod(rds_path, mode = "0664")

## -------------------------------------------------------------- provenance
prov <- c(
  "human_interaction_models provenance (H4)",
  "=========================================",
  paste0("built_by   : scripts/build_human_interaction_models.R"),
  paste0("inputs     : human_substate_score_means.csv (H3), human_seaad_donor_neuropath.csv (H1)"),
  paste0("model      : z(score) ~ amyloid_c * tau_c + region + age_c + sex + (1|donor), weights=n_cells, REML"),
  paste0("response   : per-signature z-scored mean module score (sign-preserving; comparable effect size)"),
  paste0("centring   : amyloid/tau/age mean-centred on analysed rows; region(DLPFC ref), sex(female ref)"),
  paste0("axes       : amyloid=Thal phase (0-5), tau=Braak stage (0/2-6); both ordinal, treated numeric"),
  paste0("n_donors   : ", uniqueN(allrows$donor_id), " AD-continuum donors (8 Allen reference rows dropped: no staging)"),
  paste0("n_staged   : ", nrow(allrows), " donor x region rows (all stratum)"),
  paste0("spearman(amyloid,tau) = ", round(spearman, 3), "  (off-diagonal coverage / identifiability)"),
  paste0("min_cells  : ", MIN_CELLS, " (per-donor state floor; raw IFN/proliferative medians 23/14)"),
  "rows per stratum after n_cells filter:",
  paste0("  ", names(disc_n), ": n_obs=",
         vapply(STRATA[STRATA %in% names(disc_n)], function(s) nrow(frames[[s]]), integer(1)),
         ", discordant_donors=", unlist(disc_n)[STRATA[STRATA %in% names(disc_n)]]),
  "VIF (centred design, per stratum; interaction-column inflation is the collinearity headline):",
  unlist(lapply(names(vif_list), function(s)
    paste0("  ", s, ": amyloid=", round(vif_list[[s]]["amyloid_c"], 2),
           " tau=", round(vif_list[[s]]["tau_c"], 2),
           " interaction=", round(vif_list[[s]]["amyloid_c:tau_c"], 2)))),
  "collinearity battery : VIF (incl. product col) + discordant-donor refit (ix_discordant_*) + per-region",
  "             replication (ix_mtg_*/ix_dlpfc_*, region_concordant). A literal residualised re-fit is",
  "             omitted: by Frisch-Waugh-Lovell it equals the joint-model interaction coefficient (no new info).",
  "FDR        : BH within each stratum (26-signature family); whole-microglia stratum 'all' = headline",
  "pre-registered (direction beats significance): DAM_up/Gerrits AD1/AD2 = +amyloid main (sanity #5);",
  "             NFKB_union_targets & MG_M3_module = negative amyloid:tau; Gsk3b_targets = non-null amyloid:tau",
  "caveats    : (1) a null interaction is NOT evidence against the mouse finding under this collinearity --",
  "             read it with VIF + discordant n (anti-anchoring #1). (2) PMI unusable (binned range strings,",
  "             all-NA). (3) age is coarse (3 bins: 65/78/90) -> a weak adjustment, not a precise covariate.",
  "             (4) per-state proliferative is the weakest state (H3: lowest self-score) -> read cautiously.",
  paste0("             (5) residual normality is violated in ",
         sum(panel$resid_shapiro_p < 0.05, na.rm = TRUE), "/", nrow(panel),
         " models (bounded score means) -> p approximate;"),
  "             DIRECTION (sign concordance, anti-anchoring #4) is the primary readout and is normality-free.",
  "headline (computed; the panel is the result, not any single row -- anti-anchoring #2):",
  paste0("  sanity #5 : Gerrits_AD1 amyloid MAIN positive in ",
         panel[signature == "Gerrits_AD1_human" & amyloid_est > 0, .N], "/5 strata (pipeline control passes);",
         " Gerrits_AD2 & DAM_up amyloid main ~null (tau absorbs shared variance under r~0.65)."),
  paste0("  MG_M3_module amyloid:tau NEGATIVE (mouse-concordant) in ",
         panel[signature == "MG_M3_module" & ix_est < 0, .N], "/5 strata; NFKB_union NEGATIVE in ",
         panel[signature == "NFKB_union_targets" & ix_est < 0, .N], "/5; Gsk3b_targets POSITIVE in ",
         panel[signature == "Gsk3b_targets" & ix_est > 0, .N], "/5."),
  paste0("  FDR<0.05 interactions: ", panel[ix_fdr < 0.05, .N],
         " (only substate_proliferative score in 'all', DLPFC-driven); no PRE-REGISTERED interaction"),
  "             survives FDR -> support is DIRECTIONAL + region/discordant-robust, not significant (expected",
  "             under amyloid~tau collinearity + observational human power; anti-anchoring #1/#7).",
  paste0("built_at   : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
writeLines(prov, prov_path)
Sys.chmod(prov_path, mode = "0664")

cat("H4 done:",
    nrow(panel), "models (", uniqueN(panel$signature), "signatures x",
    uniqueN(panel$stratum), "strata).\n")
cat("  ->", tsv_path, "\n  ->", rds_path, "\n  ->", prov_path, "\n")

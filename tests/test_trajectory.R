# P2-S1 trajectory pure-helper tests: score-axis, Smithson-Verkuilen squeeze, Spearman
# concordance, slingshot lineage extraction (forced single H->D + branched DAM-terminal),
# provenance. Synthetic deterministic embedding (helpers.R make_trajectory_embedding; no RNG).
# Run from project root: Rscript tests/test_trajectory.R (gate runs it under options(warn=2)).
source("R/constants.R")
source("R/microglia.R")     # reprocess_thread_env (trajectory_provenance dep)
source("R/design.R")        # factorial_design (S2a contrast-fit + decomposition tests)
source("R/de_pb.R")         # assert_complete_crossing (S2b orchestrator dep)
source("R/trajectory.R")
source("tests/helpers.R")
suppressMessages(library(slingshot))
invisible(loadNamespace("glmmTMB"))   # S3 arm: LOAD the namespace WITHOUT attaching -> glmmtmb_pt_sensitivity uses glmmTMB:: (an accidental unqualified prod call still fails, matching the target env)

# --- score_axis_pseudotime ---------------------------------------------------------------
stopifnot(isTRUE(all.equal(score_axis_pseudotime(c(0.3, 0.5), c(0.1, 0.2)), c(0.2, 0.3))))
expect_error(score_axis_pseudotime(c(0.3, 0.5), 0.1))            # length mismatch
expect_error(score_axis_pseudotime(c(0.3, NA), c(0.1, 0.2)))     # non-finite
cat("ok - score_axis_pseudotime\n")

# --- squeeze_unit_interval (Smithson-Verkuilen) ------------------------------------------
sq <- squeeze_unit_interval(c(0, 5, 10))
stopifnot(all(sq > 0 & sq < 1),                                  # opened the interval
          all(diff(sq) > 0),                                     # monotone preserved
          isTRUE(all.equal(sq, c(0.5, 1.5, 2.5) / 3)))           # (y*(n-1)+0.5)/n, n=3
sqn <- squeeze_unit_interval(c(0, NA, 10))
stopifnot(is.na(sqn[2]), isTRUE(all.equal(sqn[c(1, 3)], c(0.25, 0.75))))  # NA preserved, n=2
stopifnot(isTRUE(all.equal(squeeze_unit_interval(c(5, 5, 5)), rep(0.5, 3))))  # constant -> 0.5
expect_error(squeeze_unit_interval(c(1, NA)))                    # < 2 finite obs
cat("ok - squeeze_unit_interval\n")

# --- trajectory_concordance ---------------------------------------------------------------
stopifnot(isTRUE(all.equal(trajectory_concordance(1:10, 1:10)$rho, 1)))
stopifnot(isTRUE(all.equal(trajectory_concordance(1:10, 10:1)$rho, -1)))
cc <- trajectory_concordance(c(NA, 2:10), 1:10)                  # complete pairs only
stopifnot(cc$n == 9L, isTRUE(all.equal(cc$rho, 1)))
expect_error(trajectory_concordance(c(1, NA), c(NA, 1)))         # < 2 complete pairs
expect_error(trajectory_concordance(rep(2, 5), 1:5))             # constant pt -> Spearman undefined
expect_error(trajectory_concordance(1:5, rep(2, 5)))             # constant score -> Spearman undefined
cat("ok - trajectory_concordance\n")

# --- run_slingshot_lineage: forced single Homeostatic->DAM --------------------------------
fx <- make_trajectory_embedding(n_per = 80L)
sl <- run_slingshot_lineage(fx$embedding, fx$labels, "Homeostatic", "DAM")
stopifnot(sl$n_lineages == 1L,
          identical(sl$lineage, c("Homeostatic", "DAM")),
          !anyNA(sl$pt),
          identical(names(sl$pt), rownames(fx$embedding)),
          mean(sl$pt[fx$labels == "Homeostatic"]) < mean(sl$pt[fx$labels == "DAM"]))  # direction
cat("ok - run_slingshot_lineage single H->D\n")

# --- run_slingshot_lineage: branched, DAM-terminal selected, IFN -> NA --------------------
fb <- make_trajectory_embedding(n_per = 80L, with_ifn = TRUE)
sb <- run_slingshot_lineage(fb$embedding, fb$labels, "Homeostatic", "DAM")
stopifnot(sb$n_lineages >= 2L,                                   # branched
          sb$lineage[length(sb$lineage)] == "DAM",              # DAM-terminal lineage chosen
          all(is.na(sb$pt[fb$labels == "IFN"])),                # IFN off the DAM lineage
          !anyNA(sb$pt[fb$labels == "DAM"]))                    # all DAM cells ordered on it
expect_error(run_slingshot_lineage(fx$embedding, fx$labels, "Homeostatic", "Homeostatic"),
             "no lineage terminates")                                          # root is no leaf
expect_error(run_slingshot_lineage(fx$embedding, fx$labels, "Homeostatic", "Nope"))  # absent cluster
emb_nr <- fx$embedding; rownames(emb_nr) <- NULL
expect_error(run_slingshot_lineage(emb_nr, fx$labels, "Homeostatic", "DAM"))  # needs rownames
cat("ok - run_slingshot_lineage branched\n")

# --- trajectory_provenance ----------------------------------------------------------------
pv <- trajectory_provenance(42L)
stopifnot(is.list(pv), "slingshot" %in% names(pv$versions), pv$seed == 42L,
          length(pv$rng_kind) == 3L, nzchar(pv$r_version))
cat("ok - trajectory_provenance\n")

# --- validate_trajectory_units (per-unit metadata fail-loud audit) ------------------------
g1 <- genotype_levels[1]; g2 <- genotype_levels[2]
stopifnot(validate_trajectory_units(c("u1", "u1", "u2"), c(g1, g1, g2)))  # one geno/unit -> ok
expect_error(validate_trajectory_units(c("u1", NA), c(g1, g2)))           # missing unit id
expect_error(validate_trajectory_units(c("u1", ""), c(g1, g2)))           # blank unit id
expect_error(validate_trajectory_units(c("u1", "u2"), c(g1, "ZZZ")))      # genotype off-levels
expect_error(validate_trajectory_units(c("u1", "u1"), c(g1, g2)))         # mixed geno in one unit
cat("ok - validate_trajectory_units\n")

# =========================================================================================
# P2-S2a: estimation core (per-replicate summary + contrast fit + Kitagawa decomposition).
# =========================================================================================

# --- derive_batch (vectorised literal-prefix extraction; genotypes contain "_") ------------
gb <- c("MAPTKI_batch01", "NLGF_P301S_batch03", "P301S_batch02")
gv <- c("MAPTKI", "NLGF_P301S", "P301S")
stopifnot(identical(derive_batch(gb, gv), c("batch01", "batch03", "batch02")),  # VECTOR call
          identical(derive_batch("NLGF_MAPTKI_batch04", "NLGF_MAPTKI"), "batch04"))
expect_error(derive_batch("MAPTKI_batch01", "P301S"))            # prefix mismatch -> fail loud
expect_error(derive_batch("MAPTKI_", "MAPTKI"))                  # empty batch -> fail loud
cat("ok - derive_batch\n")

# --- pseudotime_per_replicate (per-unit summary + within-state floor) ----------------------
cf8 <- make_trajectory_cell_frame(per_state = 8L)
pr  <- pseudotime_per_replicate(cf8, lineage_states = c("Homeostatic", "DAM"))
stopifnot(nrow(pr$per_unit) == 16L, length(pr$units) == 16L,     # 4 geno x 4 batch
          identical(pr$states, c("Homeostatic", "DAM")),
          all(dim(pr$counts) == c(2L, 16L)), all(pr$counts == 8L),
          all(abs(colSums(pr$pi) - 1) < 1e-12),                  # composition columns sum to 1
          isTRUE(all.equal(sum(pr$pi_bar), 1)),
          all(c("within_homeostatic", "within_dam") %in% names(pr$per_unit)),
          all(pr$per_unit$frac_past >= 0 & pr$per_unit$frac_past <= 1),
          all(pr$per_unit$batch %in% sprintf("batch%02d", 1:4)),  # derived (no batch col in frame)
          setequal(pr$per_unit$genotype, genotype_levels),
          all(abs(pr$per_unit$within_dam - pr$per_unit$within_homeostatic - 3) < 1e-9))  # DAM base 4 vs Homeo 1
prf <- pseudotime_per_replicate(cf8, c("Homeostatic", "DAM"), min_within = 20L)  # 8 < 20
stopifnot(all(prf$within_skip))                                  # every state below floor -> skipped
cat("ok - pseudotime_per_replicate\n")

# --- ordinary_t_table / fit_trajectory_contrasts vs manual OLS -----------------------------
meta16 <- make_meta16()
fd <- factorial_design(meta16)
gi <- match(meta16$genotype, genotype_levels)
bi <- match(meta16$batch, paste0("batch0", 1:4))
y  <- 5 + 0.3 * gi + 0.5 * (gi == 4L) + ((gi * bi) %% 5L) * 0.07  # real interaction + NON-additive resid
M  <- matrix(y, nrow = 1L, dimnames = list("y", rownames(fd$design)))
ft <- fit_trajectory_contrasts(M, fd$design, fd$contrasts)
ols <- stats::lm.fit(fd$design, y)
p <- ncol(fd$design); dfres <- nrow(fd$design) - p
sigma2 <- sum(ols$residuals^2) / dfres
XtXinv <- chol2inv(chol(crossprod(fd$design)))
j <- match("tau_nlgf", colnames(fd$design))
b_int <- ols$coefficients[["tau_nlgf"]]                         # interaction == tau_nlgf coef
se_int <- sqrt(sigma2 * XtXinv[j, j])
stopifnot(sigma2 > 1e-8,                                        # design does NOT saturate -> t well-defined
          isTRUE(all.equal(ft$top$interaction$coef, b_int)),
          isTRUE(all.equal(ft$top$interaction$se,   se_int)),
          isTRUE(all.equal(ft$top$interaction$t,    b_int / se_int)),
          ft$top$interaction$df == dfres,
          identical(names(ft$top), c("tau_alone", "nlgf_in_maptki", "nlgf_in_p301s",
                                     "tau_in_nlgf", "interaction")))
cat("ok - fit_trajectory_contrasts vs OLS\n")

# --- kitagawa_channels (exact 3-channel shift-share identity) ------------------------------
prog_pr <- pseudotime_per_replicate(make_trajectory_cell_frame(per_state = 8L),
                                    c("Homeostatic", "DAM"))
kc <- kitagawa_channels(prog_pr$pi, prog_pr$mu, prog_pr$pi_bar, prog_pr$mu_bar)
stopifnot(max(abs(kc$mean_pt - (kc$comp_cf + kc$prog_cf + kc$cross - kc$const))) < 1e-8,
          max(abs(kc$comp_cf - mean(kc$comp_cf))) < 1e-8,        # pure progression -> comp_cf flat
          max(abs(kc$cross)) < 1e-8)                             # ... + zero cross term
comp_pr <- pseudotime_per_replicate(
  make_trajectory_cell_frame(per_state = 8L, dam_extra = 6L,
    adv = c(MAPTKI = 0, P301S = 0, NLGF_MAPTKI = 0, NLGF_P301S = 0)), c("Homeostatic", "DAM"))
kc2 <- kitagawa_channels(comp_pr$pi, comp_pr$mu, comp_pr$pi_bar, comp_pr$mu_bar)
stopifnot(max(abs(kc2$mean_pt - (kc2$comp_cf + kc2$prog_cf + kc2$cross - kc2$const))) < 1e-8,
          max(abs(kc2$prog_cf - mean(kc2$prog_cf))) < 1e-8,      # pure composition -> prog_cf flat
          max(abs(kc2$cross)) < 1e-8)
cat("ok - kitagawa_channels\n")

# --- decompose_progression_vs_composition (interaction loadings + reconstruction) ----------
dec_design <- function(per_rep) {                               # factorial design over the per_rep units
  m <- per_rep$per_unit[, c("genotype_batch", "genotype", "batch")]
  rownames(m) <- m$genotype_batch
  factorial_design(m)
}
dfp <- dec_design(prog_pr)
dec <- decompose_progression_vs_composition(prog_pr, dfp$design, dfp$contrasts)
stopifnot(all(c("channels", "fit", "L_int", "loadings", "interaction",
                "recon_resid_max", "balanced") %in% names(dec)),
          dec$recon_resid_max < 1e-8,
          isTRUE(all.equal(unname(dec$L_int[["mean_pt"]]), 0.4)),  # pure within-state interaction
          isTRUE(all.equal(unname(dec$loadings[["prog_cf"]]), 1)),
          abs(dec$loadings[["comp_cf"]]) < 1e-8)
dfc <- dec_design(comp_pr)
dec2 <- decompose_progression_vs_composition(comp_pr, dfc$design, dfc$contrasts)
stopifnot(dec2$recon_resid_max < 1e-8,
          isTRUE(all.equal(unname(dec2$loadings[["comp_cf"]]), 1)),
          abs(dec2$loadings[["prog_cf"]]) < 1e-8)
cat("ok - decompose_progression_vs_composition\n")

# --- weighted contrasts: EXACT per-feature WLS SE on an UNBALANCED design ------------------
# limma::contrasts.fit approximates a multi-coef contrast's SE under per-row weights (it reuses the
# unweighted coef correlation); fit_trajectory_contrasts overrides it with the exact value. An
# UNBALANCED design is required to distinguish the two -- on a balanced design both are exact.
Xw <- cbind(b0 = 1, g = c(0, 0, 0, 1, 1, 1), x = c(0, 1, 2, 0, 1, 3)); rownames(Xw) <- paste0("u", 1:6)
Cw <- cbind(g_only = c(0, 1, 0), sum_gx = c(0, 1, 1)); rownames(Cw) <- colnames(Xw)
Mw <- rbind(m1 = c(2, 3, 5, 4, 6, 9), m2 = c(1, 1, 2, 3, 5, 8)); colnames(Mw) <- rownames(Xw)
Ww <- rbind(m1 = c(1, 1, 1, 1, 1, 20), m2 = c(20, 1, 1, 1, 1, 1)); colnames(Ww) <- rownames(Xw)
ftw <- fit_trajectory_contrasts(Mw, Xw, Cw, weights = Ww)
for (i in 1:2) {
  o   <- stats::lm.wfit(Xw, Mw[i, ], Ww[i, ])
  sig <- sqrt(sum(Ww[i, ] * o$residuals^2) / o$df.residual)
  XtWXinv <- chol2inv(chol(crossprod(Xw, Xw * Ww[i, ])))
  for (cn in colnames(Cw))
    stopifnot(abs(ftw$top[[cn]]$se[i] -
                  sig * sqrt(as.numeric(t(Cw[, cn]) %*% XtWXinv %*% Cw[, cn]))) < 1e-9)
  stopifnot(max(abs(ftw$fit$coefficients[i, ] - drop(crossprod(Cw, o$coefficients)))) < 1e-9)
}
cat("ok - fit_trajectory_contrasts weighted SE exact (unbalanced)\n")

# --- fail-loud contract guards -------------------------------------------------------------
ls2 <- c("Homeostatic", "DAM")
cf_na <- make_trajectory_cell_frame(); cf_na$pt_raw <- NA_real_
expect_error(pseudotime_per_replicate(cf_na, ls2), "no on-lineage cells")          # all pt non-finite
expect_error(pseudotime_per_replicate(make_trajectory_cell_frame(), ls2, min_within = 0L), "min_within")
expect_error(pseudotime_per_replicate(make_trajectory_cell_frame(), ls2, min_within = NA_integer_), "min_within")
M_nc <- M; colnames(M_nc) <- NULL
expect_error(fit_trajectory_contrasts(M_nc, fd$design, fd$contrasts), "colnames(measure_mat)")
expect_error(kitagawa_channels(prog_pr$pi * 2, prog_pr$mu, prog_pr$pi_bar, prog_pr$mu_bar))  # cols sum to 2
cat("ok - fail-loud contract guards\n")

# =========================================================================================
# P2-S2b: progression-interaction inference (Freedman-Lane + orchestrator).
# =========================================================================================

# --- freedman_lane_interaction: null / signal / determinism / RNG-purity / weighted ---------
fdm <- factorial_design(make_meta16())
v   <- as.numeric(1:16)                                    # deterministic, NO RNG
e   <- stats::lm.fit(fdm$design, v)$residuals              # design-ORTHOGONAL -> tau_nlgf coef ~ 0
fl_null <- freedman_lane_interaction(e, fdm$design, n_perm = 999L)
stopifnot(abs(fl_null$t_obs) < 1e-6,                       # orthogonal residual -> ~0 interaction t
          fl_null$perm_p > 0.9)                            # ... so almost all |t*| >= |t_obs|
fl_sig <- freedman_lane_interaction(2 * fdm$design[, "tau_nlgf"] + 0.1 * e, fdm$design, n_perm = 999L)
stopifnot(fl_sig$perm_p < 0.05)                            # strong interaction -> small perm_p
y_det <- 2 * fdm$design[, "tau_nlgf"] + 0.1 * e
stopifnot(identical(freedman_lane_interaction(y_det, fdm$design, n_perm = 999L, seed = 7L)$perm_p,
                    freedman_lane_interaction(y_det, fdm$design, n_perm = 999L, seed = 7L)$perm_p))
set.seed(123); a1 <- runif(1)                              # RNG-purity: FL leaves the caller's stream
invisible(freedman_lane_interaction(y_det, fdm$design, n_perm = 50L))   # ... untouched (on.exit restore)
a2 <- runif(1)
set.seed(123); b1 <- runif(1); b2 <- runif(1)
stopifnot(identical(a1, b1), identical(a2, b2))
flw <- freedman_lane_interaction(y_det, fdm$design, weights = 1 + (1:16), n_perm = 99L)  # weighted path
stopifnot(is.finite(flw$t_obs), flw$perm_p >= 0, flw$perm_p <= 1)
# t_obs matches an INDEPENDENT weighted-lm oracle (lm.wfit internally, no sqrt-scaling) -> validates
# the WLS-as-OLS reduction + the pivot-free chol2inv coefficient indexing
wv  <- 1 + (1:16)
dd  <- as.data.frame(fdm$design); names(dd) <- make.names(names(dd)); dd$RESP <- y_det
lmw <- stats::lm(RESP ~ . - 1, data = dd, weights = wv)   # y_det near-collinear -> |t| huge; use RELATIVE tol
stopifnot(isTRUE(all.equal(flw$t_obs, summary(lmw)$coefficients["tau_nlgf", "t value"], tolerance = 1e-6)))
# RNG-purity edge cases: fabricates NO caller seed when none existed, and restores a NON-default kind
if (exists(".Random.seed", envir = .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
invisible(freedman_lane_interaction(y_det, fdm$design, n_perm = 10L))
stopifnot(!exists(".Random.seed", envir = .GlobalEnv))     # no .Random.seed fabricated in the caller
k0 <- RNGkind(); RNGkind("L'Ecuyer-CMRG")
invisible(freedman_lane_interaction(y_det, fdm$design, n_perm = 10L))
stopifnot(identical(RNGkind()[1], "L'Ecuyer-CMRG"))        # kind restored, not left Mersenne-Twister
do.call(RNGkind, as.list(k0))                              # restore the test's own RNG kind
cat("ok - freedman_lane_interaction\n")

# --- run_trajectory_progression: STRUCTURE on a NON-additive fixture (sigma > 0) -------------
# jitter > 0 breaks the saturated design's zero residual -> finite t (no limma zero-variance warning
# under warn=2); assert RETURN STRUCTURE + exact reconstruction only (inferential values come from
# the S2a component tests + the LIVE smoke + the gate's fresh tar_make on real data).
fake_traj <- list(
  cell_frame = make_trajectory_cell_frame(per_state = 12L, jitter = 0.3),
  provenance = list(lineage_substates = c("Homeostatic", "DAM"),
                    root_substate = "Homeostatic", terminal_substate = "DAM"))
rp <- run_trajectory_progression(fake_traj, n_perm = 199L)
stopifnot(
  all(c("per_unit", "counts", "dam_onset", "within_skip", "design", "contrasts",
        "decomposition", "permutation", "primary_family", "exploratory_family",
        "provenance") %in% names(rp)),
  all(c("weighted", "ols", "bounded") %in% names(rp$contrasts)),
  all(c("mean_pt", "median_pt", "q90", "within_homeostatic", "within_dam") %in%
        rp$contrasts$weighted$top$interaction$measure),    # interaction on every direct measure
  all(c("frac_past_logit", "frac_past_asin") %in% rp$contrasts$bounded$top$interaction$measure),
  all(c("mean_pt", "comp_cf", "prog_cf", "cross") %in% rp$decomposition$interaction$measure),
  rp$decomposition$recon_resid_max < 1e-8,                 # exact reconstruction on noisy data
  "progression_cf" %in% rp$primary_family$measure,
  "within_homeostatic" %in% rp$primary_family$measure,     # used (per_state 12 >= floor 10)
  nrow(rp$primary_family) == 2L,
  all(is.finite(rp$primary_family$p_value)), all(is.finite(rp$primary_family$fdr)),
  all(c("within_homeostatic", "within_dam") %in% names(rp$per_unit)),   # within_<lc> wire through
  all(is.finite(c(rp$permutation$mean_pt$perm_p, rp$permutation$progression_cf$perm_p,
                  rp$permutation$frac_past_logit$perm_p,
                  rp$permutation$within_homeostatic$perm_p))),
  # primary + exploratory FDR are adjusted SEPARATELY (each its own BH, not one pooled p.adjust)
  isTRUE(all.equal(rp$primary_family$fdr, p.adjust(rp$primary_family$p_value, "BH"))),
  isTRUE(all.equal(rp$exploratory_family$fdr, p.adjust(rp$exploratory_family$p_value, "BH"))),
  rp$provenance$primary_within_skipped == FALSE,           # root state retained on this fixture
  identical(rp$provenance$planned_primary, c("progression_cf", "within_homeostatic")))
cat("ok - run_trajectory_progression structure (non-additive fixture)\n")

# --- glmmtmb_pt_sensitivity (P2-S3 supportive arm) -------------------------------------------
# HEALTH GATE unit: .fit_health_ok() is the single branch point of the degrade cascade -> exercise
# each FALSE arm directly with synthetic inputs (no fragile optimiser coaxing for !pdHess / non-
# convergence / degenerate SE), so the cascade's logic is tested independent of glmmTMB's numerics.
ok_args <- list(pdHess = TRUE, convergence = 0, est = 0.3, se = 0.1, z = 3, p = 0.1, singular = FALSE)
do_ok   <- function(...) do.call(.fit_health_ok, modifyList(ok_args, list(...)))
stopifnot(
  isTRUE(do_ok()),                                         # all-healthy -> TRUE
  isFALSE(do_ok(pdHess = FALSE)),                          # non-pos-def Hessian
  isFALSE(do_ok(convergence = 1L)),                        # optimiser did not converge
  isFALSE(do_ok(se = 0)), isFALSE(do_ok(se = NA_real_)),   # zero / non-finite SE -> infinite z, zero-width CI
  isFALSE(do_ok(z = Inf)),                                 # non-finite z
  isFALSE(do_ok(p = NA_real_)), isFALSE(do_ok(p = 1.5)),   # invalid probability
  isFALSE(do_ok(singular = TRUE)), isFALSE(do_ok(singular = NA)))  # collapsed / unknown RE
cat("ok - .fit_health_ok degrade-gate branches\n")

# SUCCESS path: a non-additive fixture (jitter > 0 -> non-singular unit RE) yields a RECORDED,
# finite tau:amyloid interaction via the beta GLMM -> extraction + CI + term resolution exercised.
# Deterministic (no RNG in the fixture / nlminb optimiser; beta wins under the pinned snapshot).
gs <- glmmtmb_pt_sensitivity(make_trajectory_cell_frame(per_state = 12L, jitter = 0.3))
stopifnot(
  gs$method == "glmmTMB_beta",                            # the non-singular RE fits the beta GLMM (NOT a degrade)
  identical(gs$term, "tau:amyloid"),
  all(is.finite(c(gs$estimate, gs$se, gs$z, gs$p_value, gs$ci_l, gs$ci_r, gs$re_sd))),
  gs$ci_l < gs$ci_r, gs$p_value >= 0, gs$p_value <= 1,
  isFALSE(gs$singular), gs$n_cells == 384L, gs$n_units == 16L,  # 16 units x (12 Homeo + 12 DAM)
  is.na(gs$fail_reason),                                  # a real fit records no failure
  is.character(gs$warnings), is.character(gs$messages))
cat("ok - glmmtmb_pt_sensitivity (beta success path)\n")

# DEGRADE via SINGULAR RE: the default fixture's identical within-genotype units collapse the
# (1|unit) variance -> both the beta GLMM and the rank-normal LMM are singular -> method="failed".
gd <- glmmtmb_pt_sensitivity(make_trajectory_cell_frame(per_state = 12L))
stopifnot(gd$method == "failed",
          all(is.na(c(gd$estimate, gd$se, gd$p_value, gd$ci_l, gd$ci_r))),
          gd$n_units == 16L, grepl("singular", gd$fail_reason))  # 16 units present, both arms singular -> recorded
cat("ok - glmmtmb_pt_sensitivity (singular -> failed)\n")

# DEGRADE via NON-ESTIMABLE interaction (+ the never-RAISE contract): drop both amyloid genotypes ->
# amyloid + tau:amyloid columns are all-zero -> glmmTMB drops them as rank-deficient (a captured,
# muffled MESSAGE, NOT an exception) -> the interaction row is absent -> est NA -> BOTH arms degrade
# -> method="failed", never raised. fail_reason records "nonestimable"; the drop message is captured.
cf_noamy <- make_trajectory_cell_frame(per_state = 12L)
cf_noamy <- cf_noamy[cf_noamy$genotype %in% c("MAPTKI", "P301S"), ]
gn <- glmmtmb_pt_sensitivity(cf_noamy)
stopifnot(gn$method == "failed", is.na(gn$estimate), is.na(gn$p_value),
          gn$n_cells == 192L, gn$n_units == 8L,            # 8 units x 24 cells (subset is fine, recorded)
          grepl("nonestimable", gn$fail_reason),           # interaction dropped (rank-deficient), recorded
          any(grepl("rank-deficient", gn$messages)))       # the muffled column-drop message captured
cat("ok - glmmtmb_pt_sensitivity (non-estimable -> failed, never raises)\n")

# FAIL-LOUD on an UNKNOWN genotype: a corrupt / typo label must NOT be silently coded tau=0/amyloid=0
# (that would bias the interaction towards control). The input contract throws -> surfaces the
# upstream break (genotype_batch kept prefix-consistent so the genotype-validity guard is the gate).
cf_bad <- make_trajectory_cell_frame(per_state = 12L)
cf_bad$genotype <- as.character(cf_bad$genotype)
hit <- cf_bad$genotype == "MAPTKI"
cf_bad$genotype[hit]       <- "BOGUS"
cf_bad$genotype_batch[hit] <- sub("^MAPTKI_", "BOGUS_", cf_bad$genotype_batch[hit])
stopifnot(inherits(try(glmmtmb_pt_sensitivity(cf_bad), silent = TRUE), "try-error"))
cat("ok - glmmtmb_pt_sensitivity (unknown genotype -> fail loud)\n")

# --- P2-S4a: trajectory_report_data compact extractor ---------------------------------------
# REUSE the S2b/S3 non-additive fixture (jitter > 0 -> a real beta GLMM + sigma > 0), extended into
# a microglia_trajectory stub carrying the per_unit / sensitivity / provenance fields the extractor
# reads. score_axis_pt is the one S1 cell_frame column the fixture omits -> add it (= pt01) so the
# extractor's score-axis guard has an input. Assert the bundle carries every documented field + every
# qmd-indexed measure/contrast, and that the HARDENED guards fail loud on a dropped required name.
s4_cf <- make_trajectory_cell_frame(per_state = 12L, jitter = 0.3)
# Vary composition per unit so comp_cf + cross are non-degenerate: a constant 12/12 split gives them
# ZERO residual variance -> NaN p -> NaN fdr (the extractor's finite-fdr postcondition, valid on real
# data where composition varies). Drop the last k%%3 DAM cells (0/1/2) per unit -> every unit keeps
# >= 10 DAM (the min_within floor) and the DAM fraction varies within + across genotypes.
s4_cf <- do.call(rbind, Map(function(part, k) {
  d <- k %% 3L
  if (d == 0L) part else part[-utils::tail(which(part$substate == "DAM"), d), , drop = FALSE]
}, split(s4_cf, s4_cf$genotype_batch), seq_len(16L)))
rownames(s4_cf) <- NULL
s4_cf$score_axis_pt <- s4_cf$pt01                            # S1 column the fixture omits; extractor reads it
s4_units <- unique(s4_cf[c("genotype_batch", "genotype")])
s4_n     <- as.integer(table(s4_cf$genotype_batch)[s4_units$genotype_batch])
s4_mt <- list(
  cell_frame  = s4_cf,
  per_unit    = data.frame(genotype_batch = s4_units$genotype_batch, genotype = s4_units$genotype,
                           n_cells = s4_n, n_on_lineage = s4_n, omitted_frac = 0,
                           stringsAsFactors = FALSE),
  sensitivity = data.frame(variant = c("dims10", "dims20", "all_retained"),
                           spearman_vs_primary = c(0.99, 0.997, 0.99),
                           n_lineages = 1L, n_shared = 100L, stringsAsFactors = FALSE),
  provenance  = list(primary_dims = 15L, lineage_substates = c("Homeostatic", "DAM"),
                     root_substate = "Homeostatic", terminal_substate = "DAM",
                     concordance_rho = 0.62, concordance_floor = 0.5, concordant = TRUE,
                     dam_pt_rho = 0.56, homeo_pt_rho = -0.40, omitted_frac_overall = 0.02))
s4_tp   <- run_trajectory_progression(s4_mt, n_perm = 199L)
s4_glmm <- glmmtmb_pt_sensitivity(s4_mt$cell_frame)
trd <- trajectory_report_data(s4_mt, s4_tp, s4_glmm)
stopifnot(
  all(c("cell_frame", "interaction", "weighted_top", "decomposition", "per_unit",
        "lineage_per_unit", "sensitivity", "glmm", "provenance") %in% names(trd)),
  all(c("mean_pt", "comp_cf", "progression_cf", "cross", "within_homeostatic", "within_dam",
        "median_pt", "q90", "frac_past_logit", "frac_past_asin") %in% trd$interaction$measure),
  all(is.finite(trd$interaction$coef)), all(is.finite(trd$interaction$ci_l)),
  all(is.finite(trd$interaction$ci_r)), all(is.finite(trd$interaction$p_value)),
  all(is.finite(trd$interaction$fdr)),
  "perm_p" %in% names(trd$interaction),                      # qmd table/prose col (extractor guards existence + inlined rows)
  all(c("comp_cf", "prog_cf", "cross") %in% names(trd$decomposition$loadings)),
  all(c("tau_alone", "nlgf_in_maptki", "nlgf_in_p301s", "tau_in_nlgf", "interaction") %in%
        names(trd$weighted_top)),
  setequal(names(trd$glmm), c("method", "term", "estimate", "se", "z", "p_value", "ci_l", "ci_r",
                              "re_sd", "singular", "n_cells", "n_units", "fail_reason")),
  !any(vapply(trd$provenance, is.null, logical(1))))         # every assembled provenance field non-NULL
cat("ok - trajectory_report_data (fields + qmd-indexed measures + finite inference)\n")

# HARDENED guards FAIL LOUD on a dropped now-required name (else a silent NA-fill breaks mid-render):
expect_error(trajectory_report_data(s4_mt, s4_tp, modifyList(s4_glmm, list(n_cells = NULL))))  # 13-name glmm set
s4_tp_bad <- s4_tp; s4_tp_bad$contrasts$weighted <- NULL                                        # weighted_top gone
expect_error(trajectory_report_data(s4_mt, s4_tp_bad, s4_glmm))
s4_mt_bad <- s4_mt; s4_mt_bad$provenance$dam_pt_rho <- NULL                                      # a prose provenance source
expect_error(trajectory_report_data(s4_mt_bad, s4_tp, s4_glmm))
cat("ok - trajectory_report_data hardened guards fail loud on malformed input\n")

# perm_p column + each canonical weighted_top mean_pt row are DIRECT qmd render inputs (it$perm_p table +
# int_*$perm_p prose; mp_ctr(cn) coef/CI feed the p_ctr geom_pointrange) -> guard them too, else a dropped
# column/row breaks mid-render, not at the extractor.
s4_tp_np <- s4_tp; s4_tp_np$primary_family$perm_p <- NULL; s4_tp_np$exploratory_family$perm_p <- NULL
expect_error(trajectory_report_data(s4_mt, s4_tp_np, s4_glmm))                                   # interaction$perm_p gone
s4_tp_wm <- s4_tp
s4_tp_wm$contrasts$weighted$top$tau_alone <-
  s4_tp_wm$contrasts$weighted$top$tau_alone[s4_tp_wm$contrasts$weighted$top$tau_alone$measure != "mean_pt", ]
expect_error(trajectory_report_data(s4_mt, s4_tp_wm, s4_glmm))                                   # weighted_top mean_pt row gone
cat("ok - trajectory_report_data guards perm_p col + weighted mean_pt render inputs\n")

cat("ALL trajectory tests passed\n")

# P2-S1 trajectory pure-helper tests: score-axis, Smithson-Verkuilen squeeze, Spearman
# concordance, slingshot lineage extraction (forced single H->D + branched DAM-terminal),
# provenance. Synthetic deterministic embedding (helpers.R make_trajectory_embedding; no RNG).
# Run from project root: Rscript tests/test_trajectory.R (gate runs it under options(warn=2)).
source("R/constants.R")
source("R/microglia.R")     # reprocess_thread_env (trajectory_provenance dep)
source("R/design.R")        # factorial_design (S2a contrast-fit + decomposition tests)
source("R/trajectory.R")
source("tests/helpers.R")
suppressMessages(library(slingshot))

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

cat("ALL trajectory tests passed\n")

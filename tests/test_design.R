# S3 acceptance: the 5-contrast machinery. Verifies (1) factorial_design reproduces the exact
# coefficient weights for all 5 canonical contrasts on a synthetic 16-row meta, and (2) the
# treatment-coded factorial and the cell-means parameterisation estimate IDENTICAL contrast
# values for arbitrary responses (they span the same column space) -- a property check across
# several deterministic response vectors, stronger than a single fixed-data comparison.

source("R/constants.R")
source("R/design.R")
source("tests/helpers.R")

meta16 <- make_meta16()
fd <- factorial_design(meta16)            # add_batch = TRUE, batch present -> ~ tau+nlgf+tau_nlgf+batch
cm <- fd$contrasts

# --- (1) exact coefficient weights -------------------------------------------------------
canonical <- c("tau_alone", "nlgf_in_maptki", "nlgf_in_p301s", "tau_in_nlgf", "interaction")
stopifnot(identical(colnames(cm), canonical))               # canonical names AND order
stopifnot(all(c("tau", "nlgf", "tau_nlgf") %in% rownames(cm)))

exp_w <- list(
  tau_alone      = c(tau = 1),
  nlgf_in_maptki = c(nlgf = 1),
  nlgf_in_p301s  = c(nlgf = 1, tau_nlgf = 1),
  tau_in_nlgf    = c(tau = 1, tau_nlgf = 1),
  interaction    = c(tau_nlgf = 1)
)
for (nm in canonical) {
  col  <- cm[, nm]
  want <- stats::setNames(rep(0, length(col)), rownames(cm))
  want[names(exp_w[[nm]])] <- exp_w[[nm]]
  stopifnot(all(col == want),                               # every coefficient exactly as expected
            sum(col != 0) == length(exp_w[[nm]]))           # and no stray non-zero weights
}
# intercept + batch coefficients are never weighted (genotype contrasts are intercept/batch-free)
stopifnot(all(cm[grepl("^batch|Intercept", rownames(cm)), ] == 0))

# --- (2) factorial == cell-means for ANY response ----------------------------------------
# cell-means design from clean factor variables (column names must be syntactically valid for
# makeContrasts), genotype columns renamed to bare levels so makeContrasts can reference them.
cm_frame <- data.frame(geno  = factor(meta16$genotype, levels = genotype_levels),
                       batch = factor(meta16$batch))
dc <- model.matrix(~ 0 + geno + batch, data = cm_frame)
colnames(dc)[1:4] <- genotype_levels
cmc <- make_contrast_matrix(dc)
stopifnot(identical(colnames(cmc), canonical))

# Prove equality for EVERY response (not a handful of samples): each contrast estimate is linear in
# y, est = L %*% y with L = t(C) %*% (X'X)^-1 X' (X full column rank). Identical estimator maps L ->
# identical contrast values for ALL y in R^16. (check.attributes = FALSE: the two designs carry
# different sample rownames -- 1..16 vs s1..s16 -- but identical row order, so values align.)
pinv <- function(X) solve(crossprod(X), t(X))               # (X'X)^-1 X', X full column rank
Lf   <- t(fd$contrasts)[canonical, ] %*% pinv(fd$design)
Lc   <- t(cmc)[canonical, ]          %*% pinv(dc)
stopifnot(isTRUE(all.equal(Lf, Lc, tolerance = 1e-9, check.attributes = FALSE)))

# spot-check the contrast VALUES against hand computation on group means (y = group index)
y <- c(rep(10, 4), rep(13, 4), rep(20, 4), rep(31, 4))      # MAPTKI, P301S, NLGF_MAPTKI, NLGF_P301S
bf <- stats::coef(stats::lm.fit(fd$design, y))[rownames(fd$contrasts)]
est <- drop(t(fd$contrasts) %*% bf)
stopifnot(isTRUE(all.equal(unname(est["tau_alone"]),      13 - 10)),
          isTRUE(all.equal(unname(est["nlgf_in_maptki"]), 20 - 10)),
          isTRUE(all.equal(unname(est["nlgf_in_p301s"]),  31 - 13)),
          isTRUE(all.equal(unname(est["tau_in_nlgf"]),    31 - 20)),
          isTRUE(all.equal(unname(est["interaction"]),    (31 - 13) - (20 - 10))))

# make_contrast_matrix fails loud if a genotype level column is absent
expect_error(make_contrast_matrix(dc[, -1, drop = FALSE]), "levels")
# factorial_design fails loud on an out-of-level genotype
bad <- meta16; bad$genotype[1] <- "WT"
expect_error(factorial_design(bad), "geno")
# ... and when batch adjustment is requested (add_batch=TRUE default) but no batch column exists
expect_error(factorial_design(meta16[, "genotype", drop = FALSE]), "batch_col")

cat("ok - test_design\n")

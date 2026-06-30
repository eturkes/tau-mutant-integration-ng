# P2-S1 trajectory pure-helper tests: score-axis, Smithson-Verkuilen squeeze, Spearman
# concordance, slingshot lineage extraction (forced single H->D + branched DAM-terminal),
# provenance. Synthetic deterministic embedding (helpers.R make_trajectory_embedding; no RNG).
# Run from project root: Rscript tests/test_trajectory.R (gate runs it under options(warn=2)).
source("R/constants.R")
source("R/microglia.R")     # reprocess_thread_env (trajectory_provenance dep)
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

cat("ALL trajectory tests passed\n")

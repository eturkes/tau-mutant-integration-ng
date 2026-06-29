# S3 acceptance: substate-composition machinery. The contrast WIRING itself is proven in
# test_design.R (factorial_design == cell-means, exact weights); here we test what composition.R
# adds: the per-sample x substate count table (incl. dropping globally-empty levels + the covariate-
# constancy guard), the cross-method concordance/discordance logic, and an end-to-end propeller smoke
# on a synthetic count table with a KNOWN amyloid -> DAM shift (direction must come out right). sccomp
# is backend-gated (CmdStan) -> never invoked here; we only check the gate predicate is a logical.

source("R/constants.R")
source("R/utils.R")
source("R/design.R")
source("R/composition.R")
source("tests/helpers.R")

canonical <- c("tau_alone", "nlgf_in_maptki", "nlgf_in_p301s", "tau_in_nlgf", "interaction")

# Deterministic cell-level metadata: 16 genotype_batch samples; amyloid (NLGF*) genotypes are
# DAM-dominant, non-amyloid Homeostatic-dominant; mild per-batch jitter gives nonzero across-replicate
# variance (else propeller's eBayes sees zero residual variance). Substate factor declares the full
# 6-level set so Proliferative/ambiguous/unassigned are empty -> must be dropped. No RNG.
make_cell_meta <- function() {
  mix <- list(MAPTKI = c(40L, 8L, 3L), P301S = c(42L, 7L, 3L),         # c(Homeostatic, DAM, IFN)
              NLGF_MAPTKI = c(12L, 36L, 3L), NLGF_P301S = c(10L, 40L, 3L))
  jit <- c(batch01 = 0L, batch02 = 3L, batch03 = -2L, batch04 = 1L)    # perturb DAM count per batch
  sub_levels <- c("Homeostatic", "DAM", "IFN", "Proliferative", "ambiguous", "unassigned")
  rows <- list()
  for (gt in names(mix)) for (b in names(jit)) {
    n <- mix[[gt]]; nH <- n[1]; nD <- max(1L, n[2] + jit[[b]]); nI <- n[3]
    subs <- rep(c("Homeostatic", "DAM", "IFN"), times = c(nH, nD, nI))
    gb <- paste(gt, b, sep = "_")
    rows[[gb]] <- data.frame(genotype = gt, batch = b, genotype_batch = gb,
                             microglia_substate = factor(subs, levels = sub_levels),
                             stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, rows); rownames(out) <- NULL
  out$genotype       <- factor(out$genotype, levels = genotype_levels)
  out$batch          <- factor(out$batch)
  out$genotype_batch <- factor(out$genotype_batch)
  out
}
meta <- make_cell_meta()

# --- composition_counts: shapes, empty-level drop, canonical order, count correctness -------------
cc <- composition_counts(meta)
stopifnot(nrow(cc$counts) == 16L, ncol(cc$counts) == 3L,
          identical(cc$present_groups, c("Homeostatic", "DAM", "IFN")),          # canonical biological order
          identical(sort(cc$dropped_groups, method = "radix"),
                    c("Proliferative", "ambiguous", "unassigned")),              # globally-empty -> dropped
          identical(rownames(cc$sample_meta), rownames(cc$counts)),
          all(rowSums(cc$counts) > 0))
# count sums match a direct tabulation for every sample x group
tab_ref <- table(meta$genotype_batch, meta$microglia_substate)[rownames(cc$counts), cc$present_groups]
stopifnot(identical(unname(cc$counts), unname(matrix(as.integer(tab_ref), nrow = 16L))))
# proportions sum to 1 per sample; long form is sample x group rows with joined covariates
stopifnot(isTRUE(all.equal(unname(rowSums(cc$proportions)), rep(1, 16L))),
          nrow(cc$long) == 16L * 3L,
          all(c("sample", "cell_group", "count", "genotype", "batch") %in% names(cc$long)))
# sample_meta covariates are the 4-genotype x 4-batch crosswalk
stopifnot(setequal(cc$sample_meta$genotype, genotype_levels),
          all(table(cc$sample_meta$genotype) == 4L))

# --- covariate-constancy guard: genotype varying within a sample -> fail loud --------------------
bad <- meta
bad$genotype[bad$genotype_batch == "MAPTKI_batch01"][1] <- "P301S"
expect_error(composition_counts(bad), "length(v) == 1L")

# --- run_propeller: structure + KNOWN direction (amyloid raises DAM, lowers Homeostatic) ---------
pl <- run_propeller(cc$per_cell, cc$sample_meta, "logit")
stopifnot(is.data.frame(pl), nrow(pl) == 5L * 3L,
          setequal(pl$contrast, canonical), setequal(pl$substate, cc$present_groups),
          all(c("method", "contrast", "substate", "prop_ratio", "t", "p_value", "fdr_contrast") %in% names(pl)),
          all(pl$method == "propeller_logit"))
get1 <- function(d, ct, su, col) d[[col]][d$contrast == ct & d$substate == su]
stopifnot(get1(pl, "nlgf_in_maptki", "DAM", "t") > 0,                  # amyloid -> DAM up
          get1(pl, "nlgf_in_maptki", "Homeostatic", "t") < 0,         # ... Homeostatic down
          get1(pl, "nlgf_in_maptki", "DAM", "prop_ratio") > 1)        # ratio direction agrees
pa <- run_propeller(cc$per_cell, cc$sample_meta, "asin")              # asin transform also runs
stopifnot(nrow(pa) == 15L, all(pa$method == "propeller_asin"),
          get1(pa, "nlgf_in_maptki", "DAM", "t") > 0)                 # same direction under asin

# --- composition_concordance: agreement vs discordance flagging ----------------------------------
mk <- function(method, dir, sig) data.frame(                          # one synthetic row, controllable sign/sig
  method = method, contrast = "nlgf_in_maptki", substate = "DAM",
  t = dir, c_effect = dir, prop_ratio = 2, p_value = 0.001,
  fdr_contrast = if (sig) 0.01 else 0.5, c_fdr = if (sig) 0.01 else 0.5, stringsAsFactors = FALSE)
con_ok  <- composition_concordance(mk("propeller_logit", 1, TRUE), mk("propeller_asin", 1, TRUE),
                                   mk("sccomp", 1, TRUE))
stopifnot(nrow(con_ok) == 1L, isFALSE(con_ok$flag), isTRUE(con_ok$dir_concordant), isTRUE(con_ok$sig_concordant))
con_dir <- composition_concordance(mk("propeller_logit", 1, TRUE), mk("propeller_asin", -1, TRUE),  # sign clash
                                   mk("sccomp", 1, TRUE))
stopifnot(isTRUE(con_dir$flag), isFALSE(con_dir$dir_concordant))
con_sig <- composition_concordance(mk("propeller_logit", 1, TRUE), mk("propeller_asin", 1, FALSE),  # significance clash
                                   mk("sccomp", 1, TRUE))
stopifnot(isTRUE(con_sig$flag), isFALSE(con_sig$sig_concordant))
con_noSC <- composition_concordance(mk("propeller_logit", 1, TRUE), mk("propeller_asin", 1, TRUE))  # sccomp absent
stopifnot(isFALSE(con_noSC$flag), is.null(con_noSC$dir_sccomp))       # only propeller transforms compared

# --- sccomp gate predicate is a logical (never runs the backend in tests) ------------------------
stopifnot(is.logical(sccomp_backend_ready()), length(sccomp_backend_ready()) == 1L)

cat("ok - test_composition\n")

# Microglia substate COMPOSITION across the 5 canonical contrasts (P1-S3). Tests whether the
# amyloid -> DAM shift (the S2 headline) is statistically supported per genotype contrast, INCLUDING
# the tau x amyloid interaction. Sample unit = genotype_batch (16; 4 per genotype, fully crossed
# with 4 batches). Consumes microglia_annotated (S2); produces the composition_results target.
#
# METHOD STACK (reproducibility-tiered -- the S3 decision):
#   PRIMARY  (locked, reproducible from the pinned P3M snapshot): propeller (speckle) on the LOGIT
#     transform + limma robust eBayes. Batch is a FIXED design covariate -- consistent with the de_pb
#     DE modality. CELL-MEANS parameterisation (~ 0 + genotype + batch) via make_contrast_matrix()
#     (R/design.R): speckle's PropRatio needs per-genotype mean coefficients, and that helper's 5
#     contrasts are proven == the factorial form in tests/test_design.R, so the wiring isn't re-derived.
#   SENSITIVITY 1 (locked): propeller on the ASIN transform (Phipson 2022: asin is always real-
#     valued, robust where a proportion approaches 0/1; logit is the slightly more powerful default).
#   SENSITIVITY 2 (OPTIONAL, OFF-lock): sccomp Bayesian beta-binomial, ~ 0 + genotype + (1 | batch)
#     with a RANDOM batch intercept (few-level batch regularised by priors -- the random-vs-fixed
#     asymmetry vs propeller/limma is intentional). Cell-means form -> colon-free contrasts (no
#     backtick hazard). Runs ONLY when the CmdStan backend is provisioned (scripts/install-cmdstan.sh);
#     otherwise recorded as skipped so a fresh clone still builds propeller-only and stays green. Its
#     HMC is seed-reproducible on one platform but NOT bitwise-locked across CmdStan builds.
#   DISCORDANCE RULE (pre-declared): propeller-logit stands as THE call. Where asin or sccomp differ
#     in effect SIGN or significance for a (contrast, substate), FLAG and report it -- never average.
#
# Substate levels carrying 0 cells genome-wide (here Proliferative/ambiguous/unassigned) are DROPPED:
# an all-zero group has no variance (breaks the propeller logit and sccomp's likelihood). The drop is
# data-driven (recorded in provenance), so the function adapts if a future re-run populates a level.

# ---- pure: per-sample x substate cell-count table + aligned sample-level covariates --------------
# Aggregates cells to (sample_col x group_col); drops globally-empty group levels; orders the kept
# groups canonically (microglia_substate_levels first, any extra observed levels after). Returns the
# shapes each downstream method needs + the matched sample covariates (asserted constant within a
# sample -> fail loud, mirroring build_pseudobulk's constancy guard). Pure: no RNG, no I/O.
composition_counts <- function(meta, sample_col = "genotype_batch",
                               group_col = "microglia_substate",
                               covariate_cols = c("genotype", "batch")) {
  meta <- as.data.frame(meta)
  stopifnot(sample_col %in% names(meta), group_col %in% names(meta),
            length(covariate_cols) >= 1L, all(covariate_cols %in% names(meta)))
  declared <- if (is.factor(meta[[group_col]])) levels(meta[[group_col]]) else
    sort(unique(as.character(meta[[group_col]])), method = "radix")
  obs <- unique(as.character(meta[[group_col]]))
  ordered_levels <- c(intersect(microglia_substate_levels, obs),   # canonical biological order first
                      sort(setdiff(obs, microglia_substate_levels), method = "radix"))  # then any extras
  sample <- factor(as.character(meta[[sample_col]]))               # drops unused sample levels
  group  <- factor(as.character(meta[[group_col]]), levels = ordered_levels)
  stopifnot(!anyNA(sample), !anyNA(group))
  present <- levels(group)
  dropped <- setdiff(declared, present)
  stopifnot(length(present) >= 2L)                                 # composition needs >= 2 groups

  per_cell <- data.frame(sample = sample, cell_group = group)      # one row per CELL (propeller input)
  tab    <- table(sample = sample, cell_group = group)             # [sample x group]
  counts <- matrix(as.integer(tab), nrow = nrow(tab), dimnames = dimnames(tab))
  props  <- counts / rowSums(counts)
  samples <- rownames(counts)

  cov_list <- lapply(covariate_cols, function(cc) {                # one value per sample; fail loud if it varies
    vapply(samples, function(s) {
      v <- unique(as.character(meta[[cc]][as.character(meta[[sample_col]]) == s]))
      stopifnot(length(v) == 1L)                                   # covariate varies within a sample -> design ill-defined
      v
    }, character(1))
  })
  names(cov_list) <- covariate_cols
  sample_meta <- data.frame(cov_list, row.names = samples, stringsAsFactors = FALSE, check.names = FALSE)

  long <- as.data.frame(tab, stringsAsFactors = FALSE)             # cols: sample, cell_group, Freq
  names(long)[names(long) == "Freq"] <- "count"
  long <- cbind(long, sample_meta[as.character(long$sample), covariate_cols, drop = FALSE])
  rownames(long) <- NULL

  stopifnot(ncol(counts) == length(present), all(rowSums(counts) > 0),
            identical(rownames(sample_meta), samples),
            nrow(long) == length(samples) * length(present))
  list(counts = counts, proportions = props, per_cell = per_cell, long = long,
       sample_meta = sample_meta, present_groups = present, dropped_groups = dropped)
}

# ---- propeller: PRIMARY (logit) + sensitivity (asin). CELL-MEANS design (~ 0 + genotype + batch):
# propeller.ttest derives PropRatio by raising the per-GENOTYPE mean-proportion coefficients to the
# contrast powers, so it REQUIRES cell-means (one coefficient = one genotype's mean proportion); a
# treatment/factorial design makes PropRatio a meaningless ratio of effect coefficients AND collapses
# speckle's apply() whenever a contrast loads a single coefficient (e.g. tau_alone). Batch enters as
# fixed adjustment columns (zero contrast weight): the t-statistic/p-value (fit on the TRANSFORMED
# props with the FULL design) is batch-adjusted, while PropRatio is speckle's SEPARATE refit of the
# untransformed proportions over ONLY the contrasted genotype columns (batch dropped) -> a ratio of
# MARGINAL genotype mean proportions, exact under the balanced crossed design. NB for `interaction`
# PropRatio is a raw-scale ratio-of-ratios whose SIGN need not match the primary logit t -> read
# direction from t, never from prop_ratio, for that contrast. make_contrast_matrix() supplies the SAME
# 5 contrasts proven == the factorial form in test_design.R. propeller.ttest tests ONE contrast vector
# per call -> loop the 5 columns. Returns a tidy long df (one row per contrast x substate). ---------
run_propeller <- function(per_cell, sample_meta, transform = c("logit", "asin")) {
  transform <- match.arg(transform)
  geno  <- factor(as.character(sample_meta$genotype), levels = genotype_levels)
  batch <- factor(as.character(sample_meta$batch))
  stopifnot(!anyNA(geno), nlevels(batch) >= 2L)                    # >=2 batches -> estimable fixed adjustment
  stopifnot(all(table(geno, batch) == 1L))                        # balanced fully-crossed design (1 sample / genotype x batch): PropRatio's MARGINAL-mean reading assumes it -> fail loud on imbalance (e.g. a genotype_batch with no microglia)
  design <- stats::model.matrix(~ 0 + geno + batch)
  colnames(design) <- sub("^geno", "", colnames(design))          # genotype cols -> bare levels (makeContrasts needs level-named cols)
  rownames(design) <- rownames(sample_meta)
  stopifnot(qr(design)$rank == ncol(design))                      # full rank -> all 5 contrasts estimable
  cm <- make_contrast_matrix(design)                              # (design-cols x 5); batch rows = 0
  props <- speckle::getTransformedProps(clusters = per_cell$cell_group,
                                        sample = per_cell$sample, transform = transform)
  samp_order <- colnames(props$Counts)
  stopifnot(setequal(samp_order, rownames(design)))               # design rows <-> proportion-matrix samples
  design <- design[samp_order, , drop = FALSE]
  parts <- lapply(colnames(cm), function(cn) {
    out <- speckle::propeller.ttest(prop.list = props, design = design, contrasts = cm[, cn],
                                    robust = TRUE, trend = FALSE, sort = FALSE)
    data.frame(method = paste0("propeller_", transform), contrast = cn, substate = rownames(out),
               prop_ratio = out$PropRatio, t = out$Tstatistic, p_value = out$P.Value,
               fdr_contrast = out$FDR, stringsAsFactors = FALSE)   # fdr_contrast = BH within this contrast
  })
  res <- do.call(rbind, parts)
  rownames(res) <- NULL
  res
}

# ---- sccomp backend gate: TRUE only when the PROJECT-LOCAL (off-lock) CmdStan + cmdstanr stack is
# provisioned (the same tools/ dirs _targets.R keys its .libPaths/CMDSTAN preflight on) AND instantiate
# confirms a usable CmdStan. Requiring the project-local dirs blocks a GLOBAL cmdstanr/CmdStan (common
# on a Stan user's machine, or a stray CMDSTAN env) from silently activating the arm on a fresh clone
# -> the "fresh clone builds propeller-only" guarantee holds. Paths are relative to the project root
# (targets + the test harness both run there). ---------------------------------------------------
sccomp_backend_ready <- function() {
  isTRUE(tryCatch(
    dir.exists("tools/rlib-stan") &&                              # project-local cmdstanr lib (scripts/install-cmdstan.sh)
      length(Sys.glob(file.path("tools", "cmdstan", "cmdstan-*"))) >= 1L &&  # project-local compiled CmdStan tree
      requireNamespace("sccomp", quietly = TRUE) &&
      requireNamespace("instantiate", quietly = TRUE) &&
      instantiate::stan_cmdstan_exists(),
    error = function(e) FALSE))
}

# ---- sccomp (optional Bayesian cross-check). Cell-means ~ 0 + genotype -> design columns
# genotype<level> (colon-free; the 5 contrasts built from genotype_levels exactly mirror the
# canonical definitions). Random batch intercept. HMC, fixed seed. Records HMC convergence health two
# ways, both RECORDED not gated (the OPTIONAL non-locked arm must not fail the zero-fault gate on a
# sampler note; the LOCKED propeller path keeps full warning strictness): any R warnings -> the
# `warnings` attr; the final (outlier-removed) cmdstanr fit's diagnostic_summary() -- divergent
# transitions / treedepth / E-BFMI -- -> the `diagnostics` attr (the real signal, since sccomp emits
# divergences via message() not warning(), and per-contrast c_rhat is structurally all-NA: sccomp
# computes rhat/ESS for base design params only). Tidy long df (one row per contrast x substate). ----
run_sccomp <- function(long, seed = 42L, cores = 4L) {
  g <- paste0("genotype", genotype_levels)
  contrasts <- c(
    tau_alone      = paste(g[2], "-", g[1]),
    nlgf_in_maptki = paste(g[3], "-", g[1]),
    nlgf_in_p301s  = paste(g[4], "-", g[2]),
    tau_in_nlgf    = paste(g[4], "-", g[3]),
    interaction    = sprintf("(%s - %s) - (%s - %s)", g[4], g[2], g[3], g[1])
  )
  long <- as.data.frame(long)
  long$genotype <- factor(as.character(long$genotype), levels = genotype_levels)
  long$batch    <- factor(as.character(long$batch))
  stopifnot(!anyNA(long$genotype), nlevels(long$batch) >= 2L, all(long$count >= 0))

  warns <- character(0)
  msgs  <- character(0)
  test <- withCallingHandlers(
    {
      # `cores` is sccomp's single real parallelism knob (a fit_model formal): it caps the chain
      # count (find_optimal_number_of_chains %>% min(cores)) and sets parallel_chains/threads_per_chain
      # in the internal mod$sample() call. Pass cores ALONE: `parallel_chains` is not a formal -> it
      # falls into `...` -> duplicates that hardcoded mod$sample(parallel_chains=) arg (error); `chains`
      # IS a fit_model formal (would bind + override, not collide) but is redundant once cores governs it.
      fit <- sccomp::sccomp_estimate(
        long, formula_composition = ~ 0 + genotype + (1 | batch), formula_variability = ~ 1,
        sample = "sample", cell_group = "cell_group", abundance = "count",
        inference_method = "hmc", cores = cores, mcmc_seed = seed, verbose = FALSE)
      fit <- sccomp::sccomp_remove_outliers(fit, cores = cores, mcmc_seed = seed, verbose = FALSE)
      sccomp::sccomp_test(fit, contrasts = contrasts)
    },
    warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") },
    # cmdstanr reports sampler health (e.g. "Warning: N of M transitions ended with a divergence") via
    # message() with a literal "Warning:" prefix -- NOT R's warning() -> the warning handler never saw
    # these, so they reached stderr -> the tee'd build log -> the gate's anchored ^Warning: scan would
    # RED a FRESH composition_results rebuild (the off-lock arm must be RECORDED, never gating; the gate
    # only stays green today because it leaves composition_results cached). Capture into `msgs` and muffle
    # so they stay out of the log; the structured diagnostic_summary() below is the authoritative record.
    message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") }
  )

  # HMC convergence of the FINAL (outlier-removed) fit. withCallingHandlers evaluates its block in
  # THIS env, so `fit` leaks out here. pass_fit=TRUE stores the cmdstanr fit on attr(fit, "fit").
  # diagnostic_summary(quiet=TRUE) RETURNS the per-chain divergent / max-treedepth / E-BFMI numbers
  # WITHOUT printing cmdstanr's note (we keep only the structured values; the printed note would
  # otherwise reach the build log and trip the gate's warning scan). Total by construction: the outer
  # tryCatch plus the structural `ok` check make any surprise (fit absent, API drift, a malformed or
  # short/non-finite-count summary) degrade to NULL -- never a misleading all-zero record, never a crash.
  diagnostics <- tryCatch(local({
    f <- attr(fit, "fit")
    if (is.null(f)) return(NULL)
    ds <- f$diagnostic_summary(quiet = TRUE)
    div <- ds$num_divergent; td <- ds$num_max_treedepth; eb <- ds$ebfmi       # per-chain vectors
    nchain <- length(div)
    ok <- is.numeric(div) && is.numeric(td) && is.numeric(eb) && nchain >= 1L &&
      length(td) == nchain && length(eb) == nchain &&
      all(is.finite(div)) && all(is.finite(td))         # counts must be finite; E-BFMI may be NaN on a pathological chain -> surfaced via min, not dropped
    if (!ok) return(NULL)
    ndraws <- tryCatch(as.integer(f$metadata()$iter_sampling) * nchain, error = function(e) NA_integer_)
    list(n_chains = nchain, n_post_warmup_draws = ndraws,
         num_divergent = sum(div), divergent_by_chain = div,
         num_max_treedepth = sum(td), ebfmi_min = min(eb))
  }), error = function(e) NULL)

  test <- as.data.frame(test)
  # sccomp labels `parameter` with the contrast NAME when contrasts is named, else the expression.
  rev_map <- stats::setNames(names(contrasts), unname(contrasts))
  contrast <- ifelse(test$parameter %in% names(contrasts), test$parameter, rev_map[test$parameter])
  out <- data.frame(method = "sccomp", contrast = contrast, substate = test$cell_group,
                    c_effect = test$c_effect, c_lower = test$c_lower, c_upper = test$c_upper,
                    c_pH0 = test$c_pH0, c_fdr = test$c_FDR, stringsAsFactors = FALSE)
  rownames(out) <- NULL
  stopifnot(nrow(out) >= 1L, !anyNA(out$contrast))                 # nonempty + every row mapped to a canonical contrast
  attr(out, "warnings") <- warns
  attr(out, "diagnostics") <- diagnostics                         # final-fit HMC health (divergences / treedepth / E-BFMI), recorded
  attr(out, "messages") <- msgs                                   # muffled sccomp/cmdstanr sampler notes -- kept out of the gate log, recorded for traceability
  out
}

# ---- cross-method concordance keyed by (contrast, substate). Direction = sign of the propeller
# t-statistic / sccomp c_effect (both = effect on the logit proportion, same contrast coding).
# Significance = FDR < alpha (propeller within-contrast FDR; sccomp c_FDR). Flags any (contrast,
# substate) where available methods disagree on sign or significance -> the discordance to report. ---
composition_concordance <- function(propeller_logit, propeller_asin, sccomp_df = NULL, alpha = 0.05) {
  key <- function(d) paste(d$contrast, d$substate, sep = "||")
  base <- propeller_logit[, c("contrast", "substate")]
  bkey <- key(base)
  stopifnot(!anyDuplicated(bkey))                                 # base keys unique -> 1:1 matching
  matched <- function(d, what) {                                  # every base (contrast,substate) MUST be present in `d`
    i <- match(bkey, key(d))                                      # else a missing/partial/empty method row would drop to NA
    if (anyNA(i))                                                 # and false-green as "no disagreement" downstream
      stop(sprintf("composition_concordance: method '%s' is missing %d (contrast, substate) row(s)",
                   what, sum(is.na(i))), call. = FALSE)
    d[i, ]
  }
  m_logit <- matched(propeller_logit, "propeller_logit")
  m_asin  <- matched(propeller_asin,  "propeller_asin")
  out <- data.frame(
    contrast = base$contrast, substate = base$substate,
    dir_logit = sign(m_logit$t), sig_logit = m_logit$fdr_contrast < alpha,
    dir_asin  = sign(m_asin$t),  sig_asin  = m_asin$fdr_contrast  < alpha,
    stringsAsFactors = FALSE)
  if (!is.null(sccomp_df)) {
    m_sc <- matched(sccomp_df, "sccomp")
    out$dir_sccomp <- sign(m_sc$c_effect)
    out$sig_sccomp <- m_sc$c_fdr < alpha
  }
  dir_cols <- grep("^dir_", names(out), value = TRUE)
  sig_cols <- grep("^sig_", names(out), value = TRUE)
  out$dir_concordant <- apply(out[dir_cols], 1L, function(x) { x <- x[!is.na(x)]; length(unique(x)) <= 1L })
  out$sig_concordant <- apply(out[sig_cols], 1L, function(x) { x <- x[!is.na(x)]; length(unique(x)) <= 1L })
  out$flag <- !(out$dir_concordant & out$sig_concordant)          # TRUE -> discordance to report
  rownames(out) <- NULL
  out
}

# ---- orchestrator -> the composition_results target. propeller always runs (locked primary +
# sensitivity). sccomp runs iff its PROJECT-LOCAL backend is provisioned: backend ABSENT -> recorded
# skip (a fresh clone stays green); backend PRESENT but sccomp errors STRUCTURALLY -> LOUD (so a broken
# optional arm cannot hide behind a green gate) unless allow_sccomp_failure downgrades it to a skip.
# Sampler-note WARNINGS never fail the build (recorded in provenance) -- only a hard error does. -----
test_composition <- function(microglia_annotated, seed = 42L, cores = NULL, alpha = 0.05,
                             allow_sccomp_failure = FALSE) {
  meta <- microglia_annotated@meta.data
  cc <- composition_counts(meta)

  prop_logit <- run_propeller(cc$per_cell, cc$sample_meta, "logit")
  prop_asin  <- run_propeller(cc$per_cell, cc$sample_meta, "asin")
  prop_logit$fdr_global <- stats::p.adjust(prop_logit$p_value, method = "BH")  # across all 5x|groups| logit tests

  if (is.null(cores)) cores <- max(1L, parallel::detectCores() - 1L)
  sccomp_df <- NULL
  sccomp_status <- "skipped: project-local CmdStan backend absent (run scripts/install-cmdstan.sh for the Bayesian arm)"
  sccomp_warnings <- character(0)
  sccomp_diagnostics <- NULL
  sccomp_messages <- character(0)
  if (sccomp_backend_ready()) {
    sccomp_df <-
      if (allow_sccomp_failure)
        tryCatch(run_sccomp(cc$long, seed = seed, cores = cores),
                 error = function(e) { sccomp_status <<- paste0("error: ", conditionMessage(e)); NULL })
      else run_sccomp(cc$long, seed = seed, cores = cores)        # backend provisioned -> a structural failure is loud
    if (!is.null(sccomp_df)) {
      sccomp_warnings <- attr(sccomp_df, "warnings") %||% character(0)
      sccomp_diagnostics <- attr(sccomp_df, "diagnostics")
      sccomp_messages <- attr(sccomp_df, "messages") %||% character(0)
      attr(sccomp_df, "warnings") <- NULL
      attr(sccomp_df, "diagnostics") <- NULL
      attr(sccomp_df, "messages") <- NULL
      div_note <- if (!is.null(sccomp_diagnostics))                # HMC health caveat surfaced in status (recorded, never gated)
        sprintf("; HMC %d/%s divergent, %d max-treedepth, min E-BFMI %.2f across %d chains",
                sccomp_diagnostics$num_divergent,
                ifelse(is.na(sccomp_diagnostics$n_post_warmup_draws), "?", sccomp_diagnostics$n_post_warmup_draws),
                sccomp_diagnostics$num_max_treedepth, sccomp_diagnostics$ebfmi_min, sccomp_diagnostics$n_chains)
      else ""
      sccomp_status <- sprintf("ran (hmc, outlier-removed; %d warning(s) recorded%s)",
                               length(sccomp_warnings), div_note)
    }
  }

  concordance <- composition_concordance(prop_logit, prop_asin, sccomp_df, alpha = alpha)

  # headline: DAM under the amyloid contrasts (the S2 amyloid -> DAM claim, now tested)
  amyloid_contrasts <- c("nlgf_in_maptki", "nlgf_in_p301s", "interaction")
  headline <- prop_logit[prop_logit$substate == "DAM" & prop_logit$contrast %in% amyloid_contrasts,
                         c("contrast", "substate", "prop_ratio", "t", "p_value", "fdr_global")]
  rownames(headline) <- NULL

  # off-lock backend identity, RECORDED (the reproducibility blocker is non-pinned, so traceability =
  # capturing the resolved versions/path) -- only meaningful when sccomp actually ran.
  sccomp_backend <- if (!is.null(sccomp_df)) tryCatch(c(
    cmdstanr     = as.character(utils::packageVersion("cmdstanr")),
    cmdstan      = as.character(cmdstanr::cmdstan_version()),
    cmdstan_path = normalizePath(cmdstanr::cmdstan_path())),
    error = function(e) c(cmdstanr = NA_character_, cmdstan = NA_character_, cmdstan_path = NA_character_)) else NULL

  provenance <- list(
    seed = seed, alpha = alpha, sample_unit = "genotype_batch", n_samples = nrow(cc$counts),
    present_groups = cc$present_groups, dropped_groups = cc$dropped_groups,
    method_tiers = c(primary = "propeller_logit", sensitivity1 = "propeller_asin",
                     sensitivity2 = "sccomp (optional, off-lock)"),
    batch_handling = "fixed design covariate in propeller/limma; random (1|batch) intercept in sccomp",
    discordance_rule = "propeller-logit stands; asin/sccomp sign/significance disagreements are flagged, never averaged",
    prop_ratio_semantics = "raw MARGINAL genotype mean-proportion ratio (speckle refits on the contrasted genotype columns, batch dropped); for `interaction` a raw-scale ratio-of-ratios whose sign may differ from the primary logit t -> direction comes from t/p_value",
    sccomp_status = sccomp_status, sccomp_warnings = sccomp_warnings,
    sccomp_messages = sccomp_messages, sccomp_diagnostics = sccomp_diagnostics, sccomp_backend = sccomp_backend,
    versions = c(speckle = as.character(utils::packageVersion("speckle")),
                 sccomp = tryCatch(as.character(utils::packageVersion("sccomp")), error = function(e) NA_character_)))

  list(counts = cc$counts, proportions = cc$proportions,
       propeller_logit = prop_logit, propeller_asin = prop_asin,
       sccomp = sccomp_df, sccomp_status = sccomp_status,
       concordance = concordance, headline = headline, provenance = provenance)
}

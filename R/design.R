# 2x2 factorial design + the 5 canonical contrasts, shared by every DE method (de_pb.R
# now; P1 single-cell DE later). Two equivalent parameterisations, BOTH keyed by the same
# 5 names (constants.R contrast_definitions + the difference-of-differences `interaction`):
#   - factorial_design()    : treatment-coded ~ tau + nlgf + tau_nlgf [+ batch]; coefficients
#                             mapped to the 5 contrasts by hand.
#   - make_contrast_matrix(): cell-means form (~ 0 + genotype[+ batch]); contrasts expressed
#                             as differences of per-genotype means via limma::makeContrasts.
# They span the same column space -> identical estimable contrast values for ANY response
# (asserted in tests/test_design.R). Names: tau_alone, nlgf_in_maptki, nlgf_in_p301s,
# tau_in_nlgf, interaction.

# Cell-means contrast matrix. `design` must carry one column per genotype level
# (~ 0 + genotype, columns named by level); extra columns (e.g. batch) are allowed and take
# zero weight. ALL design column names must be syntactically valid R names (limma::makeContrasts
# requirement) -> rename genotype columns to bare levels and build batch from a named factor.
# Returns a (design-cols x 5) contrast matrix, columns in canonical order.
make_contrast_matrix <- function(design, levels = genotype_levels) {
  stopifnot(all(levels %in% colnames(design)))
  limma::makeContrasts(
    tau_alone      = P301S       - MAPTKI,
    nlgf_in_maptki = NLGF_MAPTKI - MAPTKI,
    nlgf_in_p301s  = NLGF_P301S  - P301S,
    tau_in_nlgf    = NLGF_P301S  - NLGF_MAPTKI,
    interaction    = (NLGF_P301S - P301S) - (NLGF_MAPTKI - MAPTKI),
    levels         = design
  )
}

# Treatment-coded 2x2 factorial. tau = P301S background, nlgf = amyloid, tau_nlgf = their
# interaction. genotype is refactored to canonical levels first; any value outside
# genotype_levels -> NA -> fail loud (no silent dropped level). Batch is included iff add_batch
# (default TRUE), which then REQUIRES batch_col present -> fail loud otherwise (no silent omit);
# a batch-free modality (e.g. GeoMx) MUST pass add_batch=FALSE.
# Returns list(design, contrasts): `contrasts` is (design-cols x 5), mapping coefficients to
# the 5 canonical contrasts; intercept + batch rows stay 0 (genotype contrasts are
# intercept/batch-free).
factorial_design <- function(meta, genotype_col = "genotype",
                             batch_col = "batch", add_batch = TRUE) {
  stopifnot(genotype_col %in% names(meta))
  geno <- factor(as.character(meta[[genotype_col]]), levels = genotype_levels)
  stopifnot(!anyNA(geno))
  tau  <- as.integer(geno %in% c("P301S", "NLGF_P301S"))
  nlgf <- as.integer(geno %in% c("NLGF_MAPTKI", "NLGF_P301S"))
  df <- data.frame(tau = tau, nlgf = nlgf, tau_nlgf = tau * nlgf,
                   row.names = rownames(meta))
  if (add_batch) {
    stopifnot(!is.null(batch_col), batch_col %in% names(meta))   # asked for batch -> column MUST exist (no silent omit)
    df$batch <- factor(meta[[batch_col]])
    stopifnot(!anyNA(df$batch),          # missing batch -> model.matrix drops rows -> rowname shape error
              nlevels(df$batch) >= 2L)   # a single-level batch would be silently dropped by model.matrix
    design <- stats::model.matrix(~ tau + nlgf + tau_nlgf + batch, data = df)
  } else {
    design <- stats::model.matrix(~ tau + nlgf + tau_nlgf, data = df)
  }
  rownames(design) <- rownames(meta)
  stopifnot(qr(design)$rank == ncol(design))   # full rank -> all 5 contrasts estimable (else limma only warns)

  cn <- c("tau_alone", "nlgf_in_maptki", "nlgf_in_p301s", "tau_in_nlgf", "interaction")
  cm <- matrix(0, nrow = ncol(design), ncol = length(cn),
               dimnames = list(colnames(design), cn))
  cm["tau",      "tau_alone"]      <- 1
  cm["nlgf",     "nlgf_in_maptki"] <- 1
  cm["nlgf",     "nlgf_in_p301s"]  <- 1
  cm["tau_nlgf", "nlgf_in_p301s"]  <- 1
  cm["tau",      "tau_in_nlgf"]    <- 1
  cm["tau_nlgf", "tau_in_nlgf"]    <- 1
  cm["tau_nlgf", "interaction"]    <- 1
  list(design = design, contrasts = cm)
}

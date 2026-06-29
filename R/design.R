# Design matrix and contrast builders shared by all DE methods. The
# contrast names are canonical across the project; downstream consumers
# rely on exactly these five labels.

# Names of contrasts in the limma sense (numerator - denominator).
make_contrast_matrix <- function(design, levels = genotype_levels) {
  stopifnot(all(levels %in% colnames(design)))
  cm <- limma::makeContrasts(
    nlgf_in_maptki = NLGF_MAPTKI - MAPTKI,
    nlgf_in_p301s  = NLGF_P301S  - P301S,
    tau_alone      = P301S       - MAPTKI,
    tau_in_nlgf    = NLGF_P301S  - NLGF_MAPTKI,
    interaction    = (NLGF_P301S - P301S) - (NLGF_MAPTKI - MAPTKI),
    levels         = design
  )
  cm
}

# Build a 2x2 factorial design matrix + 5-contrast matrix from a sample
# metadata frame. Mirrors NEBULA's parameterisation so cross-method betas
# are directly comparable.
factorial_design <- function(meta, genotype_col = "genotype",
                             batch_col = "batch", add_batch = TRUE) {
  geno <- factor(as.character(meta[[genotype_col]]), levels = genotype_levels)
  tau  <- as.integer(geno %in% c("P301S", "NLGF_P301S"))
  nlgf <- as.integer(geno %in% c("NLGF_MAPTKI", "NLGF_P301S"))
  df <- data.frame(tau = tau, nlgf = nlgf, tau_nlgf = tau * nlgf,
                   row.names = rownames(meta))
  if (add_batch && !is.null(batch_col) && batch_col %in% names(meta)) {
    df$batch <- factor(meta[[batch_col]])
    design <- model.matrix(~ tau + nlgf + tau_nlgf + batch, data = df)
  } else {
    design <- model.matrix(~ tau + nlgf + tau_nlgf, data = df)
  }
  rownames(design) <- rownames(meta)

  cn <- c("tau_alone", "nlgf_in_maptki", "interaction",
          "nlgf_in_p301s", "tau_in_nlgf")
  cm <- matrix(0, nrow = ncol(design), ncol = length(cn),
               dimnames = list(colnames(design), cn))
  cm["tau",      "tau_alone"]      <- 1
  cm["nlgf",     "nlgf_in_maptki"] <- 1
  cm["tau_nlgf", "interaction"]    <- 1
  cm["nlgf",     "nlgf_in_p301s"]  <- 1
  cm["tau_nlgf", "nlgf_in_p301s"]  <- 1
  cm["tau",      "tau_in_nlgf"]    <- 1
  cm["tau_nlgf", "tau_in_nlgf"]    <- 1

  list(design = design, contrasts = cm, meta = df)
}

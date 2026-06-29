# Single-cell DE methods: NEBULA (whole-microglia + per-substate) and
# glmmTMB. Both fit a 2x2 factorial NB GLMM with a (1 | id_col) random
# intercept. Returns are list-shaped to match the pseudobulk helpers in
# de_pb.R so downstream cross-method consolidation is uniform.

# Fit NEBULA (negative-binomial mixed model) on a Seurat microglia object.
#
# Uses a 2x2 factorial parameterisation of the four genotypes so that three
# of the five canonical contrasts (tau_alone, nlgf_in_maptki, interaction)
# are direct model coefficients with NEBULA-reported SE and p-value. The
# remaining two (nlgf_in_p301s, tau_in_nlgf) are computed as linear
# combinations of those betas with variances assembled from the per-gene
# covariance matrix (NEBULA's `covariance = TRUE` output).
#
# The random intercept is on `id_col` (genotype_batch, 16 levels), which
# absorbs batch and sex effects (both fully nested within genotype_batch in
# this study, so they cannot enter as fixed effects).
#
# Returns: list(top = named list of 5 tibbles, fit = raw nebula output,
#               n_genes, n_cells, design).
fit_nebula_microglia <- function(seurat_obj,
                                 id_col       = "genotype_batch",
                                 genotype_col = "genotype",
                                 assay        = "RNA",
                                 layer        = "counts",
                                 min_cell_frac = 0.01,
                                 ncore        = max(1L, parallel::detectCores() - 2L),
                                 symbol_map   = NULL,
                                 nebula_kwargs = list()) {
  stopifnot(id_col %in% colnames(seurat_obj@meta.data))
  stopifnot(genotype_col %in% colnames(seurat_obj@meta.data))

  meta <- seurat_obj@meta.data
  ord  <- order(as.character(meta[[id_col]]))
  meta <- meta[ord, , drop = FALSE]
  cells <- rownames(meta)

  counts <- GetAssayData(seurat_obj, assay = assay, layer = layer)[, cells, drop = FALSE]

  min_cells <- ceiling(min_cell_frac * ncol(counts))
  keep <- Matrix::rowSums(counts > 0) >= min_cells
  counts <- counts[keep, , drop = FALSE]
  message(sprintf("[nebula] retaining %d / %d genes expressed in >= %d cells (%.1f%%).",
                  nrow(counts), length(keep), min_cells, 100 * min_cell_frac))

  geno <- factor(as.character(meta[[genotype_col]]), levels = genotype_levels)
  tau  <- as.integer(geno %in% c("P301S", "NLGF_P301S"))
  nlgf <- as.integer(geno %in% c("NLGF_MAPTKI", "NLGF_P301S"))
  desmat <- cbind(
    `(Intercept)` = 1L,
    tau           = tau,
    nlgf          = nlgf,
    tau_nlgf      = tau * nlgf
  )
  # NEBULA's Rcpp layer requires the design matrix to be double-typed; an
  # integer-storage matrix triggers "Wrong R type for mapped matrix".
  storage.mode(desmat) <- "double"

  sid <- as.integer(factor(as.character(meta[[id_col]]),
                           levels = sort(unique(as.character(meta[[id_col]])))))
  off <- log(Matrix::colSums(counts))

  base_args <- c(list(
    id         = sid,
    pred       = desmat,
    offset     = off,
    model      = "NBGMM",
    method     = "LN",
    covariance = TRUE,
    ncore      = 1L,
    verbose    = FALSE
  ), nebula_kwargs)
  message(sprintf("[nebula] fitting on %d cells x %d genes x %d subjects (chunked over %d workers) ...",
                  ncol(counts), nrow(counts), length(unique(sid)), ncore))
  t0 <- Sys.time()

  # NEBULA's internal future-based parallelism serialises workers via
  # md5sum, which hits R's long-vector limit (2^31) on a count matrix of
  # this size. Workaround: split genes into chunks and run NEBULA with
  # ncore = 1 inside each fork-based BiocParallel worker — the dgCMatrix
  # is shared copy-on-write so we sidestep serialisation entirely.
  fit <- if (ncore <= 1L) {
    do.call(nebula::nebula, c(list(count = counts), base_args))
  } else {
    chunks <- split(seq_len(nrow(counts)),
                    cut(seq_len(nrow(counts)), ncore, labels = FALSE))
    bp <- BiocParallel::MulticoreParam(workers = ncore, progressbar = FALSE)
    chunk_fits <- BiocParallel::bplapply(chunks, function(idx) {
      do.call(nebula::nebula,
              c(list(count = counts[idx, , drop = FALSE]), base_args))
    }, BPPARAM = bp)
    summary_df  <- do.call(rbind, lapply(chunk_fits, `[[`, "summary"))
    cov_df      <- do.call(rbind, lapply(chunk_fits, `[[`, "covariance"))
    re_df       <- do.call(rbind, lapply(chunk_fits, `[[`, "random_effect"))
    od_df       <- do.call(rbind, lapply(chunk_fits, `[[`, "overdispersion"))
    conv_vec    <- unlist(lapply(chunk_fits, `[[`, "convergence"), use.names = FALSE)
    algo_vec    <- unlist(lapply(chunk_fits, `[[`, "algorithm"),   use.names = FALSE)
    # NEBULA numbers gene_id locally per chunk; reset to the global index.
    summary_df$gene_id <- seq_len(nrow(summary_df))
    list(summary = summary_df, covariance = cov_df,
         random_effect = re_df, overdispersion = od_df,
         convergence = conv_vec, algorithm = algo_vec)
  }
  message(sprintf("[nebula] done in %.1f s. summary rows: %d.",
                  as.numeric(difftime(Sys.time(), t0, units = "secs")),
                  nrow(fit$summary)))

  top <- assemble_nebula_top(fit, symbol_map = symbol_map)

  list(fit = fit, top = top, design = desmat,
       n_genes = nrow(counts), n_cells = ncol(counts))
}

# Build the five per-contrast top tibbles from a raw NEBULA fit list
# (output of `nebula::nebula()` or the chunked equivalent assembled in
# `fit_nebula_microglia`). Pulled out so that the cached fit can be
# reprocessed (e.g. after adjusting BH filtering rules) without re-running
# the ~25 min model fit.
#
# Args:
#   fit:        named list with `summary`, `covariance`, `convergence`.
#   symbol_map: optional data frame with `ensembl` and `symbol` columns.
#
# The BH adjustment is computed only over genes with convergence == 1.
# Non-converged genes keep their raw P-value but receive adj.P.Val = NA.
# NEBULA convergence codes: 1 = converged, 0 = test-only, -10 = too few
# cells per subject, -50 = failed.
assemble_nebula_top <- function(fit, symbol_map = NULL) {
  s   <- fit$summary
  cov <- fit$covariance
  stopifnot(!is.null(s), !is.null(cov))
  # NEBULA's `logFC_*` are on natural-log scale (log link); convert to log2
  # for direct comparability with limma-voom log2FC. Variances scale by k^2.
  k    <- 1 / log(2)
  k2   <- k * k
  beta_tau  <- s$logFC_tau        * k
  beta_nlgf <- s$logFC_nlgf       * k
  beta_int  <- s$logFC_tau_nlgf   * k

  # Packed lower-triangular column-major layout for a 4-predictor model:
  #   cov_1 = (1,1)             intercept var
  #   cov_2 = (2,1), cov_3 = (3,1), cov_4 = (4,1)
  #   cov_5 = (2,2)             tau var
  #   cov_6 = (3,2)             cov(tau, nlgf)
  #   cov_7 = (4,2)             cov(tau, tau_nlgf)
  #   cov_8 = (3,3)             nlgf var
  #   cov_9 = (4,3)             cov(nlgf, tau_nlgf)
  #   cov_10 = (4,4)            tau_nlgf var
  var_tau      <- cov$cov_5  * k2
  var_nlgf     <- cov$cov_8  * k2
  var_int      <- cov$cov_10 * k2
  cov_tau_int  <- cov$cov_7  * k2
  cov_nlgf_int <- cov$cov_9  * k2

  gene_ids <- s$gene
  stopifnot(!is.null(gene_ids))
  symb <- if (!is.null(symbol_map)) {
    symbol_map$symbol[match(gene_ids, symbol_map$ensembl)]
  } else NA_character_

  conv <- if (!is.null(fit$convergence)) fit$convergence else rep(NA_integer_, length(gene_ids))

  make_tbl <- function(est, var) {
    se   <- sqrt(pmax(var, 0))
    z    <- est / se
    p    <- 2 * pnorm(-abs(z))
    ok   <- !is.na(conv) & conv == 1
    padj <- rep(NA_real_, length(p))
    if (any(ok)) padj[ok] <- p.adjust(p[ok], method = "BH")
    tibble::tibble(
      gene        = gene_ids,
      symbol      = symb,
      logFC       = est,
      se          = se,
      t           = z,
      P.Value     = p,
      adj.P.Val   = padj,
      convergence = conv
    )
  }

  list(
    nlgf_in_maptki = make_tbl(beta_nlgf, var_nlgf),
    tau_alone      = make_tbl(beta_tau,  var_tau),
    interaction    = make_tbl(beta_int,  var_int),
    nlgf_in_p301s  = make_tbl(beta_nlgf + beta_int,
                              var_nlgf + var_int + 2 * cov_nlgf_int),
    tau_in_nlgf    = make_tbl(beta_tau  + beta_int,
                              var_tau  + var_int + 2 * cov_tau_int)
  )
}

# Per-substate NEBULA wrapper: fit the whole-microglia model independently
# within each microglia state so we can ask "which state carries the
# divergence signal for gene X?". Returns a named list keyed by state, with
# each element matching the shape of `fit_nebula_microglia()`.
#
# Args:
#   seurat_obj   : Seurat object containing at least all four canonical states.
#   state_col    : meta.data column carrying the state factor.
#   states       : character vector of states to fit. Defaults to all levels
#                  of `seurat_obj@meta.data[[state_col]]` that have non-zero
#                  cells in every level of `id_col`.
#   id_col       : grouping factor for the NEBULA random intercept.
#   genotype_col : factor giving genotype.
#   min_cell_frac: per-state expression-prevalence threshold for gene
#                  inclusion. 5% by default (higher than the whole-microglia
#                  1% because substates are sparser).
#   ncore        : chunk-based BiocParallel workers; passed through.
#   symbol_map   : optional ENSEMBL->symbol map.
#
# The fitter skips states whose subject coverage is degenerate (any id with
# zero cells) and reports which states it actually fit.
fit_nebula_per_state <- function(seurat_obj,
                                 state_col     = "state",
                                 states        = NULL,
                                 id_col        = "genotype_batch",
                                 genotype_col  = "genotype",
                                 assay         = "RNA",
                                 layer         = "counts",
                                 min_cell_frac = 0.05,
                                 ncore         = max(1L, parallel::detectCores() - 2L),
                                 symbol_map    = NULL,
                                 nebula_kwargs = list()) {
  stopifnot(state_col %in% colnames(seurat_obj@meta.data))
  state_vec <- seurat_obj@meta.data[[state_col]]
  if (is.null(states)) states <- levels(factor(state_vec))

  id_vec <- as.character(seurat_obj@meta.data[[id_col]])
  id_levels <- sort(unique(id_vec))

  results <- list()
  for (st in states) {
    cells_in_state <- which(as.character(state_vec) == st)
    if (length(cells_in_state) == 0L) {
      message(sprintf("[per-state nebula] skipping '%s': 0 cells.", st))
      next
    }
    id_in_state <- table(factor(id_vec[cells_in_state], levels = id_levels))
    zero_ids <- names(id_in_state)[id_in_state == 0L]
    if (length(zero_ids) > 0L) {
      message(sprintf("[per-state nebula] skipping '%s': empty id(s) %s.",
                      st, paste(zero_ids, collapse = ", ")))
      next
    }
    message(sprintf("[per-state nebula] fitting '%s' on %d cells across %d ids ...",
                    st, length(cells_in_state), length(id_in_state)))
    sub <- subset(seurat_obj, cells = colnames(seurat_obj)[cells_in_state])
    t0 <- Sys.time()
    results[[st]] <- fit_nebula_microglia(
      sub,
      id_col        = id_col,
      genotype_col  = genotype_col,
      assay         = assay,
      layer         = layer,
      min_cell_frac = min_cell_frac,
      ncore         = ncore,
      symbol_map    = symbol_map,
      nebula_kwargs = nebula_kwargs
    )
    results[[st]]$state <- st
    message(sprintf("[per-state nebula] '%s' done in %.1f min.",
                    st, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
  results
}

# glmmTMB single-cell NB GLMM with random intercept (1 | genotype_batch).
#
# Per-gene fit:
#   count ~ tau + nlgf + tau:nlgf + offset(log_lib_size) + (1 | genotype_batch)
# family = nbinom2 (variance = mu + mu^2 / phi). This mirrors NEBULA's NBGMM
# but uses Laplace approximation (lme4/TMB) instead of NEBULA's LN method.
#
# Gene-wise fits are parallelised via BiocParallel::MulticoreParam (fork-based
# on Linux, so the count matrix is shared copy-on-write — no serialisation).
# Each worker fits one gene at a time; failed fits are captured as NAs.
fit_glmmtmb_microglia <- function(seurat_obj,
                                  id_col       = "genotype_batch",
                                  genotype_col = "genotype",
                                  assay        = "RNA",
                                  layer        = "counts",
                                  min_cell_frac = 0.01,
                                  ncore        = max(1L, parallel::detectCores() - 2L),
                                  symbol_map   = NULL,
                                  gene_subset  = NULL) {
  stopifnot(requireNamespace("glmmTMB", quietly = TRUE))
  stopifnot(id_col %in% colnames(seurat_obj@meta.data))
  stopifnot(genotype_col %in% colnames(seurat_obj@meta.data))

  meta <- seurat_obj@meta.data
  ord  <- order(as.character(meta[[id_col]]))
  meta <- meta[ord, , drop = FALSE]
  cells <- rownames(meta)
  counts <- GetAssayData(seurat_obj, assay = assay, layer = layer)[, cells, drop = FALSE]

  min_cells <- ceiling(min_cell_frac * ncol(counts))
  keep <- Matrix::rowSums(counts > 0) >= min_cells
  counts <- counts[keep, , drop = FALSE]
  message(sprintf("[glmmTMB] retaining %d / %d genes expressed in >= %d cells.",
                  nrow(counts), length(keep), min_cells))

  if (!is.null(gene_subset)) {
    in_set <- intersect(gene_subset, rownames(counts))
    counts <- counts[in_set, , drop = FALSE]
    message(sprintf("[glmmTMB] restricted to %d gene_subset entries.",
                    nrow(counts)))
  }

  geno <- factor(as.character(meta[[genotype_col]]), levels = genotype_levels)
  df_cell <- data.frame(
    tau            = as.integer(geno %in% c("P301S", "NLGF_P301S")),
    nlgf           = as.integer(geno %in% c("NLGF_MAPTKI", "NLGF_P301S")),
    genotype_batch = factor(as.character(meta[[id_col]])),
    log_lib        = log(Matrix::colSums(counts) + 1)
  )

  fit_one <- function(i) {
    df_cell$count <- as.numeric(counts[i, ])
    out <- tryCatch({
      m <- glmmTMB::glmmTMB(
        count ~ tau + nlgf + tau:nlgf + offset(log_lib) + (1 | genotype_batch),
        family = glmmTMB::nbinom2,
        data   = df_cell,
        control = glmmTMB::glmmTMBControl(parallel = 1L)
      )
      s  <- summary(m)$coefficients$cond
      vc <- vcov(m)$cond
      list(
        b_tau   = s["tau",       "Estimate"],
        b_nlgf  = s["nlgf",      "Estimate"],
        b_int   = s["tau:nlgf",  "Estimate"],
        se_tau  = s["tau",       "Std. Error"],
        se_nlgf = s["nlgf",      "Std. Error"],
        se_int  = s["tau:nlgf",  "Std. Error"],
        v_tau   = vc["tau",      "tau"],
        v_nlgf  = vc["nlgf",     "nlgf"],
        v_int   = vc["tau:nlgf", "tau:nlgf"],
        c_tau_int  = vc["tau",  "tau:nlgf"],
        c_nlgf_int = vc["nlgf", "tau:nlgf"],
        converged  = !is.null(m$fit$convergence) && m$fit$convergence == 0
      )
    }, error = function(e) {
      list(b_tau = NA_real_, b_nlgf = NA_real_, b_int = NA_real_,
           se_tau = NA_real_, se_nlgf = NA_real_, se_int = NA_real_,
           v_tau = NA_real_, v_nlgf = NA_real_, v_int = NA_real_,
           c_tau_int = NA_real_, c_nlgf_int = NA_real_,
           converged = FALSE)
    })
    out
  }

  bp <- if (ncore > 1L) {
    BiocParallel::MulticoreParam(workers = ncore, progressbar = FALSE)
  } else BiocParallel::SerialParam()

  message(sprintf("[glmmTMB] fitting %d cells x %d genes on %d worker(s) ...",
                  ncol(counts), nrow(counts), ncore))
  t0 <- Sys.time()
  fits <- BiocParallel::bplapply(seq_len(nrow(counts)), fit_one, BPPARAM = bp)
  message(sprintf("[glmmTMB] done in %.1f min.",
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))

  fits_df <- do.call(rbind, lapply(fits, function(x) {
    as.data.frame(x, stringsAsFactors = FALSE)
  }))
  fits_df$gene <- rownames(counts)
  fits_df$symbol <- if (!is.null(symbol_map)) {
    symbol_map$symbol[match(fits_df$gene, symbol_map$ensembl)]
  } else NA_character_

  make_tbl <- function(est, var) {
    se   <- sqrt(pmax(var, 0))
    z    <- est / se
    p    <- 2 * pnorm(-abs(z))
    ok   <- isTRUE_vec(fits_df$converged) & is.finite(p)
    padj <- rep(NA_real_, length(p))
    if (any(ok)) padj[ok] <- p.adjust(p[ok], method = "BH")
    tibble::tibble(
      gene        = fits_df$gene,
      symbol      = fits_df$symbol,
      logFC       = est / log(2),                # ln -> log2
      se          = se / log(2),
      t           = z,
      P.Value     = p,
      adj.P.Val   = padj,
      convergence = as.integer(fits_df$converged)
    )
  }

  list(
    fit       = fits_df,
    top       = list(
      tau_alone      = make_tbl(fits_df$b_tau,                          fits_df$v_tau),
      nlgf_in_maptki = make_tbl(fits_df$b_nlgf,                         fits_df$v_nlgf),
      interaction    = make_tbl(fits_df$b_int,                          fits_df$v_int),
      nlgf_in_p301s  = make_tbl(fits_df$b_nlgf + fits_df$b_int,
                                fits_df$v_nlgf + fits_df$v_int + 2 * fits_df$c_nlgf_int),
      tau_in_nlgf    = make_tbl(fits_df$b_tau  + fits_df$b_int,
                                fits_df$v_tau  + fits_df$v_int + 2 * fits_df$c_tau_int)
    ),
    n_genes  = nrow(fits_df),
    n_cells  = ncol(counts),
    design   = c("tau", "nlgf", "tau:nlgf", "(1|genotype_batch)")
  )
}

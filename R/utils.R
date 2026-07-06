# Small shared utilities. (v1's cache_or_run is dropped -- the targets store owns caching;
# its 0664 chmod was a shared-rocker artefact, obsolete in the single-user setup.)

# NULL-coalescing infix (rlang-style; kept local to avoid the dep). Used by R/plot.R.
`%||%` <- function(a, b) if (is.null(a)) b else a

# Shared Figure 6 off-diagonal rule: finite x/y rows whose absolute distance from y=x is in the
# modality's empirical extreme tail. One row per display label; Figure 7 starts from this same set.
modality_scatter_label_rows <- function(df, n_label = NULL, label_col = "label",
                                        tail_quantile = 0.998,
                                        robust_mad_min = 6) {
  stopifnot(is.data.frame(df), all(c("x", "y", label_col) %in% names(df)),
            is.null(n_label) || (is.numeric(n_label) && length(n_label) == 1L && n_label >= 1L),
            is.numeric(tail_quantile), length(tail_quantile) == 1L,
            tail_quantile > 0, tail_quantile < 1,
            is.numeric(robust_mad_min), length(robust_mad_min) == 1L,
            robust_mad_min >= 0)
  x <- df[is.finite(df$x) & is.finite(df$y), , drop = FALSE]
  if (!nrow(x)) {
    attr(x, "offdiag_cutoff") <- NA_real_
    attr(x, "offdiag_quantile_cutoff") <- NA_real_
    attr(x, "offdiag_robust_cutoff") <- NA_real_
    return(x)
  }
  x$offdiag_distance <- abs(x$y - x$x)
  quantile_cutoff <- as.numeric(stats::quantile(x$offdiag_distance, probs = tail_quantile,
                                                names = FALSE, type = 8))
  robust_cutoff <- stats::median(x$offdiag_distance) +
    robust_mad_min * stats::mad(x$offdiag_distance, constant = 1.4826)
  cutoff <- max(quantile_cutoff, robust_cutoff, 0, na.rm = TRUE)
  if (!is.finite(cutoff) || length(unique(x$offdiag_distance)) < 2L) {
    out <- x[0, , drop = FALSE]
    attr(out, "offdiag_cutoff") <- cutoff
    attr(out, "offdiag_quantile_cutoff") <- quantile_cutoff
    attr(out, "offdiag_robust_cutoff") <- robust_cutoff
    return(out)
  }
  x <- x[x$offdiag_distance >= cutoff & x$offdiag_distance > 0, , drop = FALSE]
  if (!nrow(x)) {
    attr(x, "offdiag_cutoff") <- cutoff
    attr(x, "offdiag_quantile_cutoff") <- quantile_cutoff
    attr(x, "offdiag_robust_cutoff") <- robust_cutoff
    return(x)
  }
  tie <- if ("feature" %in% names(x)) as.character(x$feature) else as.character(seq_len(nrow(x)))
  ord <- order(-x$offdiag_distance, as.character(x[[label_col]]), tie, method = "radix")
  x <- x[ord, , drop = FALSE]
  x <- x[!duplicated(x[[label_col]]), , drop = FALSE]
  if (!is.null(n_label)) x <- utils::head(x, as.integer(n_label))
  x$scatter_label_rank <- seq_len(nrow(x))
  rownames(x) <- NULL
  attr(x, "offdiag_cutoff") <- cutoff
  attr(x, "offdiag_quantile_cutoff") <- quantile_cutoff
  attr(x, "offdiag_robust_cutoff") <- robust_cutoff
  x
}

# Write a tibble/data.frame to TSV, creating the parent directory if absent. For result
# tables emitted outside the targets store (e.g. figures/tables a report reads).
write_tsv_safe <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_tsv(x, path)
  invisible(path)
}

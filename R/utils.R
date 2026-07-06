# Small shared utilities. (v1's cache_or_run is dropped -- the targets store owns caching;
# its 0664 chmod was a shared-rocker artefact, obsolete in the single-user setup.)

# NULL-coalescing infix (rlang-style; kept local to avoid the dep). Used by R/plot.R.
`%||%` <- function(a, b) if (is.null(a)) b else a

# Shared Figure 6 off-diagonal utilities. Production Figure 6/7 passes one pooled absolute
# |x-y| cutoff learned across modalities; the panel-local tail is a fallback for tests/ad hoc plots.
modality_scatter_pooled_cutoff <- function(panels, order = names(panels), tail_quantile = 0.998) {
  stopifnot(is.list(panels), is.character(order), length(order) >= 1L,
            all(order %in% names(panels)),
            is.numeric(tail_quantile), length(tail_quantile) == 1L,
            tail_quantile > 0, tail_quantile < 1)
  dist <- unlist(lapply(order, function(m) {
    d <- panels[[m]]
    if (is.list(d) && is.data.frame(d$data)) d <- d$data
    stopifnot(is.data.frame(d), all(c("x", "y") %in% names(d)))
    abs(d$y[is.finite(d$x) & is.finite(d$y)] -
          d$x[is.finite(d$x) & is.finite(d$y)])
  }), use.names = FALSE)
  dist <- dist[is.finite(dist) & dist > 0]
  if (length(dist) < 2L || length(unique(dist)) < 2L) return(NA_real_)
  as.numeric(stats::quantile(dist, probs = tail_quantile, names = FALSE, type = 8))
}

modality_scatter_label_rows <- function(df, n_label = NULL, label_col = "label",
                                        tail_quantile = 0.998,
                                        cutoff = NULL,
                                        cutoff_source = NULL) {
  stopifnot(is.data.frame(df), all(c("x", "y", label_col) %in% names(df)),
            is.null(n_label) || (is.numeric(n_label) && length(n_label) == 1L && n_label >= 1L),
            is.numeric(tail_quantile), length(tail_quantile) == 1L,
            tail_quantile > 0, tail_quantile < 1,
            is.null(cutoff) || (is.numeric(cutoff) && length(cutoff) == 1L))
  attr_cutoff <- attr(df, "offdiag_cutoff", exact = TRUE)
  attr_source <- attr(df, "offdiag_cutoff_source", exact = TRUE)
  if (is.null(cutoff) && !is.null(attr_cutoff)) cutoff <- attr_cutoff
  if (is.null(cutoff_source) && !is.null(attr_source)) cutoff_source <- attr_source
  cutoff_source <- cutoff_source %||% if (is.null(cutoff)) "panel_tail_quantile" else "provided"
  x <- df[is.finite(df$x) & is.finite(df$y), , drop = FALSE]
  if (!nrow(x)) {
    attr(x, "offdiag_cutoff") <- NA_real_
    attr(x, "offdiag_quantile_cutoff") <- NA_real_
    attr(x, "offdiag_cutoff_source") <- cutoff_source
    return(x)
  }
  x$offdiag_distance <- abs(x$y - x$x)
  quantile_cutoff <- NA_real_
  if (is.null(cutoff)) {
    quantile_cutoff <- as.numeric(stats::quantile(x$offdiag_distance, probs = tail_quantile,
                                                  names = FALSE, type = 8))
    cutoff <- quantile_cutoff
  }
  if (!is.finite(cutoff) || length(unique(x$offdiag_distance)) < 2L) {
    out <- x[0, , drop = FALSE]
    attr(out, "offdiag_cutoff") <- cutoff
    attr(out, "offdiag_quantile_cutoff") <- quantile_cutoff
    attr(out, "offdiag_cutoff_source") <- cutoff_source
    return(out)
  }
  x <- x[x$offdiag_distance >= cutoff & x$offdiag_distance > 0, , drop = FALSE]
  if (!nrow(x)) {
    attr(x, "offdiag_cutoff") <- cutoff
    attr(x, "offdiag_quantile_cutoff") <- quantile_cutoff
    attr(x, "offdiag_cutoff_source") <- cutoff_source
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
  attr(x, "offdiag_cutoff_source") <- cutoff_source
  x
}

# Write a tibble/data.frame to TSV, creating the parent directory if absent. For result
# tables emitted outside the targets store (e.g. figures/tables a report reads).
write_tsv_safe <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_tsv(x, path)
  invisible(path)
}

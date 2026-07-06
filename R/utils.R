# Small shared utilities. (v1's cache_or_run is dropped -- the targets store owns caching;
# its 0664 chmod was a shared-rocker artefact, obsolete in the single-user setup.)

# NULL-coalescing infix (rlang-style; kept local to avoid the dep). Used by R/plot.R.
`%||%` <- function(a, b) if (is.null(a)) b else a

# Shared Figure 6 label rule: finite x/y rows, ranked by distance from y=x, one row per display
# label. Figure 7 starts from this visible-label set, then substitutes parent genes for phosphosites.
modality_scatter_label_rows <- function(df, n_label = 12L, label_col = "label") {
  stopifnot(is.data.frame(df), all(c("x", "y", label_col) %in% names(df)),
            is.numeric(n_label), length(n_label) == 1L, n_label >= 1L)
  x <- df[is.finite(df$x) & is.finite(df$y), , drop = FALSE]
  if (!nrow(x)) return(x)
  tie <- if ("feature" %in% names(x)) as.character(x$feature) else as.character(seq_len(nrow(x)))
  ord <- order(-abs(x$y - x$x), as.character(x[[label_col]]), tie, method = "radix")
  x <- x[ord, , drop = FALSE]
  x <- x[!duplicated(x[[label_col]]), , drop = FALSE]
  x <- utils::head(x, as.integer(n_label))
  x$scatter_label_rank <- seq_len(nrow(x))
  rownames(x) <- NULL
  x
}

# Write a tibble/data.frame to TSV, creating the parent directory if absent. For result
# tables emitted outside the targets store (e.g. figures/tables a report reads).
write_tsv_safe <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_tsv(x, path)
  invisible(path)
}

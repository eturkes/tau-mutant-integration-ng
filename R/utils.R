# Small shared utilities. (v1's cache_or_run is dropped -- the targets store owns caching;
# its 0664 chmod was a shared-rocker artefact, obsolete in the single-user setup.)

# NULL-coalescing infix (rlang-style; kept local to avoid the dep). Used by R/plot.R.
`%||%` <- function(a, b) if (is.null(a)) b else a

# Write a tibble/data.frame to TSV, creating the parent directory if absent. For result
# tables emitted outside the targets store (e.g. figures/tables a report reads).
write_tsv_safe <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_tsv(x, path)
  invisible(path)
}

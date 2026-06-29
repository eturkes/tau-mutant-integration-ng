# Tiny utilities used everywhere: caching, NULL-fallback operator, safe
# TSV writes, NA-safe TRUE check.

# Cache helper: run an expression, save the result, return loaded value on re-run.
cache_or_run <- function(path, expr, overwrite = FALSE) {
  if (!overwrite && file.exists(path)) {
    message(sprintf("[cache hit ] %s", path))
    return(readRDS(path))
  }
  message(sprintf("[cache miss] %s -> computing", path))
  val <- force(expr)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(val, path)
  Sys.chmod(path, mode = "0664")
  val
}

# NULL-fallback operator (rlang-style); kept local to avoid the rlang dep.
`%||%` <- function(a, b) if (is.null(a)) b else a

# Safe TSV writer: creates the parent dir, sets group-readable perms.
write_tsv_safe <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_tsv(x, path)
  Sys.chmod(path, mode = "0664")
}

# Treat NA convergence as FALSE.
isTRUE_vec <- function(x) !is.na(x) & as.logical(x)

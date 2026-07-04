source("rv/scripts/rvr.R")
source("rv/scripts/activate.R")

# Fail loud if rv did not activate its library (rv missing from PATH, `rv info` failed, or
# R version != rproject.toml r_version -> activate.R safe-modes to a temp lib). Without this,
# base R silently keeps the global/site libraries -> wrong package versions, a repro hole.
# Stop non-interactive runs (pipeline, smoke tests); warn when interactive so a session
# still opens for debugging. Re-add if `rv` regenerates .Rprofile.
local({
  if (!any(grepl("rv/library", .libPaths(), fixed = TRUE))) {
    msg <- "rv library not active: run scripts/install-rv.sh, put ~/.local/bin on PATH, and match rproject.toml r_version to this R."
    if (interactive()) warning(msg, call. = FALSE) else stop(msg, call. = FALSE)
  }
})

# Redirect OmnipathR's API logs out of the repo root. Its logger defaults to writing
# `./omnipathr-log/` at the working dir; setting the logdir option before the package
# loads (it reads this at load) keeps probe/download logs in the gitignored cache tree
# instead of littering the working tree. Re-add if `rv` regenerates .Rprofile.
local({
  d <- file.path(getwd(), "storage", "cache", "omnipathr-log")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  options(omnipathr.logdir = d)
})

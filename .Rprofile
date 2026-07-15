source("rv/scripts/rvr.R")
source("rv/scripts/activate.R")

# Fail loud if rv did not activate its library (rv missing from PATH, `rv info` failed, or
# R version != rproject.toml r_version -> activate.R safe-modes to a temp lib). Without this,
# base R silently keeps the global/site libraries -> wrong package versions, a repro hole.
# Stop non-interactive runs (pipeline, smoke tests); warn when interactive so a session
# still opens for debugging. Re-add if `rv` regenerates .Rprofile.
local({
  if (!any(grepl("rv/library", .libPaths(), fixed = TRUE))) {
    msg <- paste(
      "rv library not active: run scripts/bootstrap/rv.sh, put ~/.local/bin on PATH,",
      "and match rproject.toml r_version to this R."
    )
    if (interactive()) warning(msg, call. = FALSE) else stop(msg, call. = FALSE)
  }
})

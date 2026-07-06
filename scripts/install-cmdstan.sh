#!/usr/bin/env bash
# OPTIONAL Stan backend for sccomp -- the Bayesian composition CROSS-CHECK (P1-S3). propeller
# (the reproducible PRIMARY) needs NONE of this; R/composition.R degrades to propeller-only when
# this backend is absent, so a fresh clone still builds and scripts/check.sh stays green WITHOUT
# ever running this script. Run it only to add the sccomp sensitivity arm.
#
# OFF rv.lock BY NECESSITY -- the recorded reproducibility blocker: cmdstanr lives on the Stan
# r-universe (not the pinned P3M snapshot) and CmdStan is a compiled C++ tree rv cannot lock. Both
# land project-local under tools/ (gitignored + read-economy skip), in a library SEPARATE from rv/library
# so `rv sync` never prunes them. Idempotent. The CmdStan version is RECORDED (provenance), not
# bitwise-pinned: the Bayesian arm is explicitly non-locked; the locked guarantee is propeller.
# Toolchain: build-essential (g++/make) -- already a project sysdep (scripts/install-sysdeps.sh).
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
STAN_LIB="${ROOT}/tools/rlib-stan"     # cmdstanr (+ any missing CRAN deps from P3M) -- OFF rv/library
CMDSTAN_DIR="${ROOT}/tools/cmdstan"    # compiled CmdStan tree (cmdstan-<ver>/)
mkdir -p "${STAN_LIB}" "${CMDSTAN_DIR}"

# 1. cmdstanr -> separate project-local lib. cmdstanr is ONLY on the Stan r-universe; its CRAN deps
#    (posterior/checkmate/processx/R6/data.table/...) resolve from P3M (listed first -> pinned),
#    cmdstanr itself from r-universe (listed last -> the only source). rv/library stays on .libPaths
#    so deps already present there are reused, not duplicated into STAN_LIB.
Rscript -e '
  lib <- "'"${STAN_LIB}"'"; .libPaths(c(lib, .libPaths()))
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    install.packages("cmdstanr", lib = lib,
                     repos = c(getOption("repos"), stan = "https://stan-dev.r-universe.dev"))
  }
  stopifnot(requireNamespace("cmdstanr", quietly = TRUE))
  cat("cmdstanr:", as.character(utils::packageVersion("cmdstanr")), "\n")
'

# 2. CmdStan compiled from source, project-local. cmdstanr picks the latest compatible CmdStan
#    (override by exporting CMDSTAN_VERSION); we RECORD the resolved version. Skips if already built.
Rscript -e '
  lib <- "'"${STAN_LIB}"'"; .libPaths(c(lib, .libPaths()))
  dir <- "'"${CMDSTAN_DIR}"'"
  cmdstanr::check_cmdstan_toolchain(fix = FALSE)
  if (length(Sys.glob(file.path(dir, "cmdstan-*"))) == 0L) {
    ver <- Sys.getenv("CMDSTAN_VERSION"); if (!nzchar(ver)) ver <- NULL
    cmdstanr::install_cmdstan(dir = dir, cores = parallel::detectCores(),
                              overwrite = FALSE, version = ver)
  }
  path <- Sys.glob(file.path(dir, "cmdstan-*"))[1]
  stopifnot(!is.na(path), dir.exists(path))
  cat("CMDSTAN:", path, "\n")
'

RESOLVED="$(ls -d "${CMDSTAN_DIR}"/cmdstan-* 2>/dev/null | head -n1)"
echo "Stan backend ready: ${RESOLVED}"
echo "(_targets.R prepends ${STAN_LIB} to .libPaths + sets CMDSTAN to this path when both exist)"

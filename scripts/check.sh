#!/usr/bin/env bash
# Quality gate (P0-S5). Fail-loud, zero-fault. Enforces MORE than tar_make's exit code:
# tar_make() returns 0 even when targets CAPTURE warnings, so beyond the exit code we
#   (a) loop the unit tests (warn=2 promotes any stray R warning to an error),
#   (b) assert tar_meta error+warnings are all NA across every target,
#   (c) grep the captured render log for Quarto / pandoc / knitr warnings (knitr runs the
#       .qmd in a SEPARATE R process, so those warnings never reach tar_meta -> log scan catches them).
# Any fault -> non-zero exit. A clean store makes tar_make a fast no-op; this still re-runs (a)-(c).
#
# Run:  scripts/check.sh                      full: rv sync + uv sync + tests + pipeline + enforce
#       CHECK_SKIP_SYNC=1 scripts/check.sh    skip env sync for fast local iteration
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

LOG="$(mktemp)"
trap 'rm -f "$LOG"' EXIT
trap 'echo "check.sh: GATE FAILED -> ${BASH_COMMAND}" >&2' ERR

# 1. Reproducible environment (idempotent; both fast when already satisfied). Their own
#    stdout is NOT captured into $LOG, so it cannot trip the render-log scan in 4c.
if [ "${CHECK_SKIP_SYNC:-0}" != "1" ]; then
  echo "== rv sync =="; rv sync
  echo "== uv sync =="; uv sync
fi

# 2. Unit tests: each tests/test_*.R is a fail-loud stopifnot script that prints "ok - <name>".
#    warn=2 turns any stray R warning into an error. Per-file Rscript -> isolation + a non-zero
#    exit on the first failure (set -e aborts the loop).
echo "== unit tests =="
shopt -s nullglob
tests=(tests/test_*.R)
shopt -u nullglob
[ "${#tests[@]}" -gt 0 ] || { echo "check.sh: no tests/test_*.R found" >&2; exit 1; }
for t in "${tests[@]}"; do
  Rscript -e 'options(warn = 2); source(commandArgs(TRUE)[1])' "$t"
done

# 3. Pipeline: tar_make() builds/refreshes the DAG (.Rprofile activates rv). Capture combined
#    output (incl. the Quarto render log) for the 4c scan. A HARD target error makes tar_make
#    exit non-zero -> pipefail + set -e abort here.
echo "== tar_make =="
Rscript -e 'targets::tar_make()' 2>&1 | tee "$LOG"

# 4a. Zero-fault: tar_meta error+warnings all NA (tar_make's exit code ignores CAPTURED
#     warnings; this does not). Scope to the CURRENT manifest -> drops tar_source'd
#     functions/globals AND any stale dead-target rows that could false-fail the gate.
echo "== tar_meta zero-fault =="
Rscript -e '
  cur <- targets::tar_manifest(fields = "name")$name
  m <- targets::tar_meta(fields = c("error", "warnings"))
  m <- m[m$name %in% cur, , drop = FALSE]
  bad <- m[!is.na(m$error) | !is.na(m$warnings), c("name", "error", "warnings"), drop = FALSE]
  if (nrow(bad) > 0L) { print(bad); stop("tar_meta: ", nrow(bad), " target(s) with error/warning", call. = FALSE) }
  cat("ok - tar_meta clean across", nrow(m), "current targets\n")
'

# 4c. Zero-fault: no Quarto / pandoc / knitr warning in the render log. "warning" (any case)
#     catches [WARNING] / "Warning message:" / "Warning in ..."; "\bwarn\b" catches a bare WARN.
#     command grep -> the real GNU binary, not the rg-fff shell shadow (clean, file-order).
echo "== render-log zero-fault =="
if command grep -Eni 'warning|\bwarn\b' "$LOG"; then
  echo "check.sh: render log contains the warning(s) shown above" >&2
  exit 1
fi
echo "ok - render log clean"

echo "PASS - quality gate green"

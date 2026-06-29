#!/usr/bin/env bash
# Quality gate (P0-S5). Fail-loud, zero-fault. Enforces MORE than tar_make's exit code,
# which returns 0 even when targets CAPTURE warnings and never sees warnings raised inside
# the report's SEPARATE knitr/Quarto render process. Layered enforcement:
#   1  env sync (rv + uv)        -- reproducible toolchain
#   2  unit tests at warn=2      -- a stray R warning in any test -> error
#   3  tar_make, report FORCE-RENDERED each run (cheap: reads cached targets, no heavy rebuild)
#        - the report qmd sets options(warn=2) -> any R chunk warning -> render error -> tar_make fails
#        - tar_quarto(quiet=FALSE) -> Quarto/Pandoc [WARNING] lines reach the captured log for 4c
#   4a tar_meta(error,warnings) all-NA across current targets + their dynamic branches
#   4c render-log scan for ANCHORED Quarto / Pandoc / R warning lines
# Any fault -> non-zero exit.
#
# Why force-render: a cached-current report skips its render, leaving $LOG empty and the warn=2
# path un-exercised -> a stale or out-of-band warning could pass unseen (gate non-idempotent).
# Re-rendering each run keeps the gate sound and idempotent; it reads cached targets (~0.3 GB)
# and does NOT re-run the heavy load_snrnaseq build that materialises the modality targets.
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

# 3. Pipeline. Force a fresh report render every run (invalidate -> rebuild) so the report's
#    warn=2 path and the render-log scan (4c) ALWAYS exercise -- a cached report would skip both
#    and could hide a warning. Cheap: reads cached targets, does NOT re-run the heavy
#    load_snrnaseq build. Combined output (incl. the Quarto render log) is tee'd to $LOG. Wrap in
#    `if !` so a tar_make/render failure reports cleanly (pipefail would otherwise blame `tee`).
echo "== tar_make (report force-rendered) =="
if ! Rscript -e 'targets::tar_invalidate(any_of("report")); targets::tar_make()' 2>&1 | tee "$LOG"; then
  echo "check.sh: tar_make / report render failed (see output above)" >&2
  exit 1
fi

# 4a. Zero-fault: tar_meta error+warnings all NA (tar_make's exit ignores CAPTURED warnings;
#     this does not). Scope to current manifest names AND dynamic branches of current targets
#     (name in cur OR parent in cur) -> drops tar_source'd functions/globals and stale dead-target
#     rows that could false-fail the gate, while still catching a faulted branch of a live target.
#     warn=2 -> a stray warning from this check is itself an error.
echo "== tar_meta zero-fault =="
Rscript -e '
  options(warn = 2)
  cur <- targets::tar_manifest(fields = "name")$name
  m <- targets::tar_meta(fields = c("parent", "error", "warnings"))
  keep <- m$name %in% cur | (!is.na(m$parent) & m$parent %in% cur)
  m <- m[keep, , drop = FALSE]
  bad <- m[!is.na(m$error) | !is.na(m$warnings), c("name", "error", "warnings"), drop = FALSE]
  if (nrow(bad) > 0L) { print(bad); stop("tar_meta: ", nrow(bad), " target(s) with error/warning", call. = FALSE) }
  cat("ok - tar_meta clean across", nrow(m), "current targets/branches\n")
'

# 4c. Zero-fault: no Quarto / Pandoc / R warning line in the render log. quiet=FALSE (set on the
#     tar_quarto target) lets Quarto/Pandoc [WARNING] lines reach $LOG. Match ONLY anchored
#     warning forms so benign chatter containing the substring "warn" (paths, labels, progress)
#     cannot false-red. command grep -> the real GNU binary, not the rg-fff shell shadow.
#     Discriminate grep's exit: 0 = match (fault), 1 = clean, >=2 = grep infrastructure error.
echo "== render-log zero-fault =="
command -v grep >/dev/null 2>&1 || { echo "check.sh: grep not found -> cannot scan render log" >&2; exit 1; }
# grep in the if-condition is exempt from set -e AND the ERR trap (bash leaves the ERR trap armed
# under set +e), so an expected "no match" (exit 1 = clean) neither aborts nor prints a spurious
# GATE-FAILED line. Discriminate: 0 = match (fault), 1 = clean, >=2 = grep infrastructure error.
if hits="$(command grep -nE '^\[WARNING\]|^Warning:|^Warning in |^Warning messages?:|^WARN' "$LOG")"; then
  { echo "check.sh: render log contains warning(s):"; echo "$hits"; } >&2
  exit 1
else
  gstatus=$?
  if [ "$gstatus" -eq 1 ]; then
    echo "ok - render log clean"
  else
    echo "check.sh: grep failed (status $gstatus) scanning render log" >&2
    exit 1
  fi
fi

echo "PASS - quality gate green"

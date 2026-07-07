#!/usr/bin/env bash
# Lean report gate. This repo is currently optimised for fast iteration on the
# final analysis HTML, not a broad test harness. It optionally syncs the R env,
# force-renders the report target, then scans targets metadata and the render log
# for warnings/errors.
#
# Run:  scripts/check.sh                      rv sync + force-render + enforce
#       CHECK_SKIP_SYNC=1 scripts/check.sh    skip rv sync for fastest local iteration
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

LOG="$(mktemp)"
trap 'rm -f "$LOG"' EXIT
trap 'echo "check.sh: GATE FAILED -> ${BASH_COMMAND}" >&2' ERR

# 1. R environment sync. Its stdout is not captured into $LOG, so it cannot trip
#    the render-log scan.
if [ "${CHECK_SKIP_SYNC:-0}" != "1" ]; then
  echo "== rv sync =="; rv sync
fi

# 2. Pipeline. Force a fresh report render every run (invalidate -> rebuild) so the report's
#    warn=2 path and the render-log scan (4c) ALWAYS exercise -- a cached report would skip both
#    and could hide a warning. Cheap: reads cached targets, does NOT re-run the heavy
#    load_snrnaseq build. Combined output (incl. the Quarto render log) is tee'd to $LOG. Wrap in
#    `if !` so a tar_make/render failure reports cleanly (pipefail would otherwise blame `tee`).
echo "== tar_make (report force-rendered) =="
if ! Rscript -e 'targets::tar_invalidate(any_of("report")); targets::tar_make()' 2>&1 | tee "$LOG"; then
  echo "check.sh: tar_make / report render failed (see output above)" >&2
  exit 1
fi

# 3a. Zero-fault: tar_meta error+warnings all NA (tar_make's exit ignores CAPTURED warnings;
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

# 3b. Zero-fault: no Quarto / Pandoc / R warning line in the render log. quiet=FALSE (set in the
#     report render target) lets Quarto/Pandoc [WARNING] lines reach $LOG. Match ONLY anchored
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

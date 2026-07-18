#!/usr/bin/env bash
# Fast report iteration gate. This repo is currently optimised around one final
# user-facing artifact: report/tau-mutant-integration.html.
#
# Run: scripts/check.sh
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

trap 'echo "check.sh: GATE FAILED -> ${BASH_COMMAND}" >&2' ERR

echo "== render report/tau-mutant-integration.html + validate occupancy_harness_check =="
if ! Rscript -e 'targets::tar_invalidate(tidyselect::any_of(c("report", "occupancy_harness_check"))); targets::tar_make(names = c("report", "occupancy_harness_check"))'; then
  echo "check.sh: tar_make / report render or occupancy harness check failed (see output above)" >&2
  exit 1
fi

echo "PASS - report/tau-mutant-integration.html + occupancy_harness_check"

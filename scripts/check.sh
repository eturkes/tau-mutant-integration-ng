#!/usr/bin/env bash
# Fast report iteration gate. This repo is currently optimised around one final
# user-facing artifact: report/tau-mutant-integration.html.
#
# Run: scripts/check.sh
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

trap 'echo "check.sh: GATE FAILED -> ${BASH_COMMAND}" >&2' ERR

echo "== render report/tau-mutant-integration.html =="
if ! Rscript -e 'targets::tar_invalidate("report"); targets::tar_make(names = "report")'; then
  echo "check.sh: tar_make / report render failed (see output above)" >&2
  exit 1
fi

echo "PASS - report/tau-mutant-integration.html"

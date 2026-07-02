#!/usr/bin/env bash
# Read-only prompted Codex review wrapper. Stores prompt/output under .codex/runs/.
# Generic `codex exec` is intentional: CLI 0.142.5 rejects custom stdin prompts
# with `review --uncommitted`, so project-specific criteria live in reviewer.md.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

RUN_DIR=".codex/runs/review-$(date -u +%Y%m%dT%H%M%SZ)-$$"
mkdir -p "$RUN_DIR"
PROMPT="$RUN_DIR/prompt.md"
REVIEW="$RUN_DIR/review.md"

cp .codex/prompts/reviewer.md "$PROMPT"
if [ "$#" -gt 0 ]; then
  {
    printf '\nReview focus:\n'
    printf '%s\n' "$*"
  } >> "$PROMPT"
fi

timeout "${CODEX_REVIEW_TIMEOUT:-1500}" \
  codex exec --sandbox read-only -o "$REVIEW" < "$PROMPT"

cat "$REVIEW"
printf 'review_saved=%s\n' "$REVIEW"

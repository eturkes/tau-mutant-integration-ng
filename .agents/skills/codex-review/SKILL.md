---
name: codex-review
description: Project Codex review workflow for tau-mutant-integration-ng. Use when explicitly invoked as $codex-review to replace the old /codex-review flow.
---

This is the `$` equivalent of the old `/codex-review` slash command.

Read `.codex/prompts/review.md` and `.codex/prompts/reviewer.md` completely, then run the workflow.

Execution:
- Text after `$codex-review` = review focus; pass it to `scripts/codex-review.sh`.
- No trailing text = review current cohesive uncommitted session work.
- Relay findings, mark each accepted/rejected with rationale, fix accepted findings, run the smallest meaningful check, then close with one scoped commit when the change is complete.

Implementation note:
- `scripts/codex-review.sh` uses read-only generic `codex exec` intentionally because Codex CLI 0.142.5 rejects custom stdin prompts with `review --uncommitted`.

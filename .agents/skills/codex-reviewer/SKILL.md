---
name: codex-reviewer
description: Internal read-only reviewer rubric for tau-mutant-integration-ng. Use only when explicitly invoked to inspect uncommitted changes without editing.
---

Read `.codex/prompts/reviewer.md` completely and follow it exactly.

Constraints:
- Inspect `git status --short`, `git diff`, and relevant new files.
- Do not edit files.
- Return findings first, then open questions/assumptions, then a short summary only if useful.

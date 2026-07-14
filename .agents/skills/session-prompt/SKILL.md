---
name: session-prompt
description: Project session entry workflow with preflight work sizing for tau-mutant-integration-ng. Use when explicitly invoked as $session-prompt to continue from roadmap state or run an optional task.
---

Read `.codex/prompts/session.md` completely, then follow it.

Invocation contract:
- Text after `$session-prompt` = the prompt's `Task`.
- No trailing text = empty `Task`; run the roadmap-implied mode.
- Preserve the prompt's load order, preflight sizing/check posture, and commit rules.

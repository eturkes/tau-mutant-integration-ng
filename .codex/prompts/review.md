# Codex Review Workflow

Run a read-only Codex review of this session's cohesive change and act on it.
This is the human workflow note; `scripts/codex-review.sh` sends
`.codex/prompts/reviewer.md` as the runtime reviewer prompt.

Default command:

```bash
scripts/codex-review.sh
```

Focused command:

```bash
scripts/codex-review.sh "focus text"
```

After the review:
- Relay findings.
- Mark each accepted or rejected with rationale.
- Fix accepted findings.
- Re-run the smallest check that exercises the fix; run `scripts/check.sh` for code/report changes.
- Close with one scoped commit.

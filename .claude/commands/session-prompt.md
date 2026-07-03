Continue this project (fresh session). Non-empty task below ⇒ your sole task: do exactly it, editing `.agent/roadmap.md` only if it directs you to. Empty ⇒ run the MODE (below) from the roadmap's state.

Load `.agent/roadmap.md` (forward index: posture · cohesive story · active plan · phased backlog), then `.agent/memory.md` (standing contract: science facts, env, quality gate). CLAUDE.md (imports `AGENTS.md`) is auto-injected. Pull `map.md` (wiring) + `history.md` (decision digests) as a step implicates them; mine `archive_digest.md` for v1 threads (full v1 code on branch `archive`). Read only what the step implicates; navigate via tokensave or LSP (Serena) where available, else grep.

Posture (roadmap) = **BUILD**: streamlined rebuild from scratch, one cohesive story. Verify the quality gate (memory.md → `scripts/check.sh`) before AND after every change. Each mode ends in one scoped commit (convention below).

MODE ← roadmap state:
- no active plan → PLAN (confirm the next backlog phase with me first)
- active plan has an open step → EXECUTE its next open step
- active plan's steps all done, unreviewed → CLOSE-OUT

PLAN — open new multi-step work.
- Research first: read-only finders (`Explore`) + web-search for the SOTA fit (AGENTS.md); mine `archive_digest.md`; `git status`-reconcile. Spawn subagents to spare the window.
- Write `.agent/<topic>_plan.md`: scope, ordered steps each closeable in one window (acceptance per step); gate-independent prep first; flag gated / heavy-compute steps (planned, awaiting their gate).
- Decision gate: present the plan default PLUS ≥1 reasoned alternative via `AskUserQuestion`.
- Open: set roadmap `Active plan`, commit `roadmap (<topic> plan): …`.

EXECUTE — advance the active plan's next open step.
- Orient via the prior step's commit + ledger; restate the step + its acceptance in one line; implement reusing modules, matching style.
- Gated step? Its precondition must be met — confirm it FUNCTIONALLY through the pipeline/tooling; unmet ⇒ stop + report, so every result traces to real inputs.
- Verify the gate green; record durable lessons/decisions in `.agent/memory.md`; mark the step done; update `map.md` if wiring changed.
- Commit `<scope> (<topic>): …`.

CLOSE-OUT — all steps done.
- Adversarially review the plan body per AGENTS.md (correctness, scope, token-efficiency, obsolescence).
- Fold decisions into `history.md`; archive the plan → `.agent/completed/<topic>_plan_YYYY-MM-DD.md` (deny-Read); reset roadmap `Active plan` to (none).
- Commit `<scope> (<topic> close): …`.

Commits: scoped (scopedcommits.com), one per cohesive piece, LLM-parse-optimized; I handle remote. After a commit I may run `/codex-review` → fix accepted findings in a follow-up `<scope> (<topic> review): … (codex)`. Watch headroom via `.agent/context.sh`; near 80%, drive to a clean state and close out.

Task (may be empty): $ARGUMENTS

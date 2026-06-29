Continue this project (fresh session). Non-empty task below ⇒ your sole task: do exactly it, editing `.agent/roadmap.md` only if it directs you to. Empty ⇒ resume the active plan's next open step; if none, report the MAINTENANCE status + backlog and take direction (upkeep is the default; new analysis arcs need an explicit request).

Load `.agent/roadmap.md` (forward index: posture · active plan · backlog), then `.agent/memory.md` (standing contract: locked decisions, gotchas, quality gate). CLAUDE.md (imports `AGENTS.md`) is auto-injected. Pull `map.md` (codebase wiring) + `history.md` (decision digests) only as a step implicates them; archived plans in `.agent/completed/` are deny-Read → consult `history.md`. Read only what the step implicates; navigate via LSP where available, else grep.

Posture = **MAINTENANCE**: the arc programme (D->O + capstone) is sealed. Verify the project's quality gate (memory.md) before AND after every change.

MODE ← roadmap state:

PLAN — I direct new multi-step work; no active plan.
- Research first: read-only finders (`Explore`) + web-search for the SOTA fit (AGENTS.md); `git status`-reconcile. Spawn subagents to spare the 200K window.
- Write `.agent/<topic>_plan.md`: scope, ordered steps each closeable in one window (acceptance per step), gate-independent prep sequenced first, gated steps flagged.
- Decision gate: present the plan default PLUS ≥1 reasoned alternative before `AskUserQuestion`.
- Open: set roadmap `Active plan`, commit `roadmap (<topic> plan): …`.

EXECUTE — active plan has an open step.
- Restate the step + its acceptance in one line; implement reusing modules, matching surrounding style.
- GATE: a gated step needs its precondition met functionally (resolve through the pipeline/tooling); deny-listed inputs stay off-limits; unmet ⇒ stop and report.
- Verify the quality gate; record durable lessons/decisions in `.agent/memory.md`; mark the step done.
- Commit `<scope> (<topic>): …`.

CLOSE-OUT — all steps done.
- Adversarially review the plan body per AGENTS.md (correctness, scope/memory conformance, token-efficiency, obsolescence); fix what you find.
- Fold decisions into `history.md`; archive the plan → `.agent/completed/<topic>_plan_YYYY-MM-DD.md`; reset roadmap `Active plan` to (none).
- Commit `<scope> (<topic> close): …`.

Commits: scoped (scopedcommits.com), one per cohesive piece, LLM-parse-optimized; I handle remote. I may run `/codex-review` after a commit → fix accepted findings before closing. Watch headroom via `.agent/context.sh`; near 80%, drive to a clean state and close out.

Task (may be empty): $ARGUMENTS

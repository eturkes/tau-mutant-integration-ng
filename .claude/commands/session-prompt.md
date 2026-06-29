Continue this project (fresh session). Non-empty task below ⇒ your sole task: do
exactly it, editing `.agent/roadmap.md` only if it directs you to. Empty ⇒ resume
the active plan's next open step; if none, take direction on the next backlog phase.

Load `.agent/roadmap.md` (forward index: posture · cohesive story · active plan ·
phased backlog), then `.agent/memory.md` (standing contract: science facts, env,
quality gate). CLAUDE.md (imports `AGENTS.md`) is auto-injected. Pull `map.md`
(wiring) + `history.md` (new decision digests) as a step implicates them; mine
`archive_digest.md` for v1 promising threads (full v1 code on branch `archive`).
Read only what the step implicates; navigate via LSP (Serena) where available, else grep.

Posture = **BUILD**: streamlined rebuild from scratch, one cohesive story. Verify the
quality gate (memory.md) before AND after every change.

MODE ← roadmap state:

PLAN — new multi-step work; no active plan.
- Research first: read-only finders (`Explore`) + web-search for the SOTA fit
  (AGENTS.md); mine `archive_digest.md`; `git status`-reconcile. Spawn subagents to
  spare the window.
- Write `.agent/<topic>_plan.md`: scope, ordered steps each closeable in one window
  (acceptance per step), gate-independent prep first, gated / heavy-compute steps flagged.
- Decision gate: present the plan default PLUS ≥1 reasoned alternative before `AskUserQuestion`.
- Open: set roadmap `Active plan`, commit `roadmap (<topic> plan): …`.

EXECUTE — active plan has an open step.
- Restate the step + its acceptance in one line; implement reusing modules, matching style.
- Verify the quality gate; record durable lessons/decisions in `.agent/memory.md`;
  mark the step done; update `map.md` if wiring changed.
- Commit `<scope> (<topic>): …`.

CLOSE-OUT — all steps done.
- Adversarially review the plan body per AGENTS.md (correctness, scope, token-efficiency, obsolescence).
- Fold decisions into `history.md`; archive the plan → `.agent/completed/<topic>_plan_YYYY-MM-DD.md`
  (deny-Read); reset roadmap `Active plan` to (none).
- Commit `<scope> (<topic> close): …`.

Commits: scoped (scopedcommits.com), one per cohesive piece, LLM-parse-optimized; I
handle remote. I may run `/codex-review` after a commit → fix accepted findings before
closing. Watch headroom via `.agent/context.sh`; near 80%, drive to a clean state and close out.

Task (may be empty): $ARGUMENTS

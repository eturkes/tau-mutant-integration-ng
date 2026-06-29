# Roadmap (forward index)

Where we are going; read on launch BEFORE picking up work. The `.agent/` memory
system: `memory.md` = standing contract (rules, locked decisions, gotchas) ·
`map.md` = codebase wiring (where things are) · `history.md` = decision history
(where we have been) · `roadmap.md` = this file (where we are going).

Posture: **MAINTENANCE**. The planned analysis programme (arcs D->O plus the
arc-P capstone) closed 2026-06-15; see history.md. New analysis arcs require an
explicit user request — default to upkeep, not new science.

## Active plan
(none) — `report_styling_plan` closed 2026-06-16 (S0-S4: shared bslib v5 theme +
human-facing prose/figure dash sweep to a hyphen/comma register; see history.md).
Back to MAINTENANCE posture.

When no active plan, this section reads "(none)". Open one via `/session-prompt`;
multi-step work lives in `.agent/<topic>_plan.md`, is archived to
`.agent/completed/<topic>_plan_YYYY-MM-DD.md` on close (archived full plans are
deny-Read — consult history.md instead), and its decisions are folded into
history.md.

## Backlog (candidate upkeep; grounded, unordered)
- **security / dependency audit** — CLAUDE.md schedules these periodically; none
  run since close. Bump R / Bioconductor / Python deps, then re-verify a clean
  knit and the locked numbers.
- **reproducibility check** — exercise the cold clean-clone rebuild path:
  recreate `.venv` + the `.micromamba` scenic env from the gitignore recipes,
  pre-build the [S] caches, full knit of all 3 reports from scratch.
- **doc / notes compaction** — prune drift in map.md + history.md as arcs age and
  the codebase moves.
- **refactor sweep** — periodic R/ modularity pass (CLAUDE.md KISS / dedup
  standing instruction); no specific debt flagged yet.

## Parked analysis directions
(none) — the arc programme is closed; surface new analysis arcs only on an
explicit user request.

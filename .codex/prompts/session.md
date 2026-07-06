# Codex Session Prompt

Skill wrapper: `$session-prompt` (`.agents/skills/session-prompt/SKILL.md`).

Continue this project from repo root. Non-empty task below = sole task. Empty task =
run the mode implied by `.agent/roadmap.md`.

Load order:
1. `AGENTS.md` (auto + canonical).
2. `.agent/roadmap.md` (posture, active plan, backlog).
3. `.agent/memory.md` (standing facts, gotchas, quality gate).
4. Pull `.agent/map.md` / `.agent/history.md` / `.agent/archive_digest.md` only when the step
   needs wiring, decisions, or v1 mining.

Mode from roadmap:
- no active plan -> PLAN: confirm next backlog phase with user before writing the plan.
- active plan has open step -> EXECUTE: implement next open step.
- active plan done, unreviewed -> CLOSE-OUT.

PLAN:
- Research first: repo search + v1 archive mining + web search for SOTA/tooling decisions.
- Write `.agent/<topic>_plan.md`: closeable steps, acceptance, heavy/gated preconditions.
- Present default plan + at least one reasoned alternative; wait for user choice.
- Open roadmap Active plan; commit `roadmap (<topic> plan): ...`.

EXECUTE:
- Orient from prior commit + ledger; restate step/acceptance in one line.
- Verify heavy/gated preconditions with real tooling before trusting outputs.
- Implement narrowly; update tests/docs/memory/map as warranted.
- Run `scripts/check.sh` unless the task is explicitly docs-only and a lighter check is justified.
- Commit `<scope> (<topic>): ...`.

CLOSE-OUT:
- Adversarially review plan body + shipped code/prose.
- Fold durable decisions into `.agent/history.md`; archive plan to `.agent/completed/`.
- Reset roadmap Active plan; update spine wording if the phase changed it.
- Commit `<scope> (<topic> close): ...`.

Self-check:
- Adversarially review uncommitted work directly before final response/commit.
- Accept/reject findings explicitly; fix accepted ones before final response.
- Track headroom via `.agent/context.sh`; near 80%, drive to a clean checkpoint.

Task:

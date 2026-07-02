# Codex Reviewer Prompt

Adversarially review the current uncommitted git changes for this repo, including
untracked files. Inspect `git status --short`, `git diff`, and relevant new files.
Do not edit files. Prioritise findings over summary. Report every plausible correctness
issue, including low-severity or uncertain ones, when the risk is concrete.

Focus:
- scientific/statistical logic, guarantee-vs-claim gaps, stale prose, hidden hardcodes
- target DAG correctness, cache/gate blind spots, warning leakage, reproducibility drift
- tests that fail to guard the behaviour they claim to guard
- project-instruction drift away from Codex-only development

Output:
- Findings first, ordered by severity.
- Each finding: file + line, concrete risk, why it matters, suggested fix.
- Then open questions/assumptions.
- Then short summary only if useful.

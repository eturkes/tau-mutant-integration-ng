# P5 Synthesis - plan

## Scope
Close the fresh rebuild with one lean synthesis layer over the analyses already
shipped in P1-P4. P5 is read-only biological synthesis: it re-aggregates existing
compact targets, writes a concise final narrative, and adds a compact evidence
table. It must not add a new modality, a new statistical contest, or the v1
ledger/adjudication machinery.

Inputs already built:
- P1: `microglia_report`, `composition_results`, `pb_de_microglia`,
  `pb_de_substate`.
- P2: `trajectory_report`.
- P3: `mechanism_report`.
- P4: `crossmodality_report`, with `crossmodality_divergence` available only if
  the compact P4 bundle lacks a needed audit field.

Outputs:
- `synthesis_report` compact target.
- `_synthesis.qmd`, included near the top of `index.qmd` so the detailed
  chapters become the audit trail.
- A final lean-pass over chapter pointers and overview text.

Out of scope:
- No new biological inference, meta-analysis, MOFA/factor model, human
  validation, SCENIC/topology, SpatialDecon, full CCC, claims ledger, contest
  margins, or publication-figure side lane.
- No heavy report-time target reads. `_synthesis.qmd` loads only
  `synthesis_report`.

## Research Digest
Local / v1:
- Current report wiring already matches the cheap-render pattern: each analysis
  chapter loads one compact report target, and the force-rendered Quarto target
  is warning-fatal under `options(warn=2)`.
- P1-P4 closed the planned evidence spine: amyloid drives microglial
  homeostatic-to-DAM activation; mutant tau modulates that response mainly as
  extra DAM-cell composition rather than supported further trajectory
  progression; Myc has RNA support; NF-kB attenuation and Gsk3b are not
  recovered; cross-modality strengthens amyloid-response / synaptic-clearance
  axes while keeping SpatialDecon/full CCC unearned.
- v1 capstone's useful rule is read-only convergence with nulls shown
  explicitly. Its bloat to avoid: the 11-arc ledger, contest arithmetic,
  margin-neutral row bookkeeping, human layer, PNG side lane, and heavy prose.
- Existing current rebuild already has the material P5 needs. The synthesis
  should compress and adjudicate wording, not rerun or broaden the analysis.

Current-method / tooling sweep:
- Quarto includes are the right existing structure: included `.qmd` files share
  the same engine and underscore-prefixed include files are ignored as standalone
  project renders. Source: https://quarto.org/docs/authoring/includes.html
- `tar_quarto()` detects literal `tar_load()` / `tar_read()` calls in active
  report chunks and wires them as dependencies, which supports the existing
  one-compact-target-per-chapter pattern. Source:
  https://search.r-project.org/CRAN/refmans/tarchetypes/html/tar_quarto.html
- Quarto supports HTML/raw tables, tabsets, responsive figures, anchors and other
  HTML affordances without adding a new R table dependency. Sources:
  https://quarto.org/docs/authoring/tables.html and
  https://quarto.org/docs/output-formats/html-basics.html
- MOFA2 remains a current multi-omics factor-analysis option, but it is designed
  for unsupervised latent factors across matched or overlapping matrices. P5's
  question is not latent-factor discovery; it is target-derived synthesis of a
  pre-specified 2x2 story with modality-specific caveats. Source:
  https://bioc.r-universe.dev/MOFA2

## Default Design
Default = compact, read-only synthesis target + upfront synthesis chapter.

1. Build a small `synthesis_report` target.
   It consumes compact P1-P4 report targets and returns:
   - `headline`: ordered key statements.
   - `evidence_table`: about 8-14 rows, one row per claim / caveat, with stable
     fields:
     `claim_id`, `axis`, `status`, `direction`, `evidence`,
     `primary_sources`, `supporting_sources`, `caveat`, `report_anchor`.
   - `status_summary`: counts by status, not a contest score.
   - `open_questions`: explicit unsupported / unearned claims.
   - `provenance`: source target names, R version, alpha where relevant, and a
     short report contract.

2. Use status labels, not scored models.
   Suggested enum:
   `core_supported`, `corroborated`, `focused_support`, `not_supported`,
   `not_earned`, `open_caveat`.
   These are descriptive table labels. They do not imply a v1-style graded
   ledger or quantitative model contest.

3. Put synthesis before details.
   `index.qmd` keeps the overview short, includes `_synthesis.qmd` after the
   overview, then QC and the P1-P4 chapters. The reader gets the answer first;
   the existing chapters remain the audit trail.

4. Keep the synthesis honest about negatives.
   Nulls and unearned claims are first-class rows: no supported
   progression-beyond-composition signal, NF-kB attenuation not supported,
   Gsk3b not recovered, SpatialDecon/full CCC not earned, bulk run-index
   sensitivity, and bulk hippocampus not microglia-sorted.

5. Do a final lean pass.
   Remove stale "before final synthesis" / "P5 remains" wording from
   human-facing report text, tighten duplicated caveats only where it improves
   clarity, and preserve the detailed chapters' audit value.

## Alternatives
Alternative A - v1-style convergence matrix / claims ledger.
Pros: maximum audit structure and continuity with the archive. Cons: directly
revives the bloat the rebuild excluded, invites contest-margin over-reading, and
would make P5 a new adjudication system rather than a closing synthesis. Rejected
by default.

Alternative B - MOFA2 / latent factor integration.
Pros: current, formal multi-omics integration tooling. Cons: it answers an
unsupervised discovery question, needs reshaping unmatched modalities into
factor-model matrices, and would create new biological claims after the planned
analysis is already closed. Keep as a future exploratory branch, not P5.

Alternative C - Standalone short summary only, dropping detailed chapters.
Pros: shortest deliverable. Cons: loses the warning-gated audit trail and makes
the single HTML less self-verifying. Default keeps one report with a synthesis
front section plus details.

Default choice: compact target-derived synthesis. It closes the backlog, matches
the current DAG/report design, and gives the reader a concise answer without
re-opening inference.

## Steps
Each step is one closing unit. Resuming mid-plan: read this plan's Scope + your
step + `.agent/memory.md` relevant report/target sections + files named in the
step. Run `scripts/check.sh` unless explicitly docs-only.

### S1 - Compact synthesis target [DONE 2026-07-02]
Add `R/synthesis.R`, `tests/test_synthesis.R`, and target `synthesis_report`.

Contracts:
- `synthesis_report_data()` consumes compact report targets only. It may read
  `crossmodality_divergence` at target-build time only if `crossmodality_report`
  lacks a required audit slice; `_synthesis.qmd` still loads only
  `synthesis_report`.
- Guard every source field the synthesis body reads. Build-fatal anchors:
  amyloid-to-DAM support, DAM composition interaction, trajectory
  `comp_cf` / `progression_cf` rows, Myc interaction TF row, NF-kB gate status,
  Gsk3b rows, GeoMx/bulk run-index caveats, SpatialDecon status, and earned
  clearance-pair rows (including the empty/unearned case).
- Evidence rows use descriptive statuses only. No support/contradict columns,
  net scores, model margins, or ledger row IDs.
- Output is small and render-safe: no heavy Seurat object, no large harmonised
  table, no per-cell frame except a deliberately tiny optional summary.

Acceptance:
- Synthetic tests cover required-anchor extraction, status enum validation,
  missing-anchor failures, empty earned-pair handling, and no-ledger-column
  invariant.
- Fresh `synthesis_report` build warning-clean with `tar_meta` clean.
- `tar_manifest()` shows `synthesis_report` depends only on compact targets
  unless a documented `crossmodality_divergence` slice is intentionally needed.
- `scripts/check.sh` green.

Result:
- Built `synthesis_report_data()` over `microglia_report`, `trajectory_report`,
  `mechanism_report`, and `crossmodality_report` only; `crossmodality_divergence`
  was not needed.
- Live `synthesis_report` is ~4.8KB, warning-clean, and carries 10 descriptive
  evidence rows with status labels only.
- `tests/test_synthesis.R` covers required anchors, status enum validation,
  missing-anchor failures, empty earned-pair handling, and the no-ledger-column
  invariant.
- Manifest/raw dependency check: `synthesis_report_data(microglia_report,
  trajectory_report, mechanism_report, crossmodality_report)`.

### S2 - Synthesis chapter + report wiring [OPEN]
Add `_synthesis.qmd`, include it in `index.qmd` immediately after Overview, and
rewrite Overview as the entry point to the final report.

Contracts:
- `_synthesis.qmd` has `options(warn=2)`, `tar_source()`, and
  `targets::tar_load(synthesis_report)` only.
- It presents: one short answer paragraph, a compact evidence table, one small
  status/axis figure if it earns its keep, and a short "what remains unsupported"
  paragraph.
- All numbers and status labels are inline-computed from `synthesis_report`.
- Human-facing prose uses British English and hyphens for asides.

Acceptance:
- Force-rendered report is warning-clean; `_synthesis.qmd` reads only
  `synthesis_report`.
- `index.qmd` no longer says the final synthesis is still open.
- The generated HTML remains self-contained/offline and the report target stays
  cheap to force-render.
- `scripts/check.sh` green.

### S3 - Lean report pass + phase-ready close [TODO]
Review the full rendered report and source prose for stale pointers, duplicated
chapter conclusions, and overclaims introduced before P5.

Contracts:
- Preserve detailed P1-P4 chapters as audit trail, but remove stale forward
  pointers and tighten only repetitive summary text.
- Keep negative findings explicit: no supported progression beyond composition,
  NF-kB attenuation not supported, Gsk3b not recovered, SpatialDecon/full CCC not
  earned, and bulk run-index sensitivity.
- Do not alter the settled P1-P4 biological calls except to make them more
  concise or more honest.

Acceptance:
- Stale-pointer search across `index.qmd` and `_*.qmd` has no human-facing report
  hits for `P5`, `final synthesis`, `still open`, or `before the final synthesis`.
- Full `scripts/check.sh` green.
- Adversarial self-review of the synthesis table and report prose has accepted
  findings fixed before marking P5 steps complete.

After S3, enter CLOSE-OUT mode: fold durable P5 decisions into `history.md`,
archive this plan to `.agent/completed/`, reset roadmap Active plan, update the
roadmap spine if the synthesis changes wording, and commit the close-out.

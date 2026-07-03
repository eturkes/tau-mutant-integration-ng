# Figure-Caption-Only Report Plan

Opened: 2026-07-03. Task: extreme report reduction: no prose; just figures and
captions.

## Goal

Visible report path = headings + figures + captions only.

Hard target:
- `index.qmd` + `_*.qmd` render no visible paragraph/list/table/provenance prose.
- Allowed visible text = section headings, figure text, figure captions, axis labels,
  legends, and compact caption microcopy.
- Tables become figure panels or disappear. `knitr::kable()` and visible `cat()`
  provenance are out of the main path.
- `fig-alt` stays source-level accessibility text. It must not become visible body
  prose.

Current baseline, after the top-section removal commits:

| kind | blocks | words | status |
|---|---:|---:|---|
| paragraph | 33 | 594 | remove/figure-encode |
| caption | 45 | 388 | keep/compress |
| h1 | 5 | 34 | keep as navigation |
| h2 | 26 | 103 | keep as navigation |
| total inventory | 109 | 1,119 | source-count baseline |

S1 strict baseline (`uv run python scripts/prose_inventory.py --strict --html
_report/index.html --summary-only`, expected red before conversion):
- Source: 109 blocks / 1,119 words; allowed = 31 headings + 45 captions; blockers =
  33 paragraphs / 594 words.
- Rendered HTML main path: blockers = 34 paragraphs / 716 words (incl. visible author
  metadata), 4 tables / 168 words, 0 text-only `.cell-output-display`, 4 stdout
  provenance blocks / 80 words.
- Total strict blocker records = 75 (source + rendered-HTML lenses; source/rendered
  prose intentionally overlap).

Success = body prose blocks == 0 and no visible non-figure tables/provenance
outputs in rendered HTML. Captions stay short: target <=14 words each unless a
reviewed scientific caveat needs more; hard fail above 24 words.

## Research Digest

Local:
- `scripts/prose_inventory.py --summary-only` now reports 1,119 counted words
  across 109 blocks: 33 paragraph blocks / 594 words remain, mostly QC setup,
  chapter opening sentences, integration/status tails, and visible provenance.
- Prior visual reduction already proved the figure infrastructure: compact
  per-chapter figure targets, 49 rendered figures/captions, warning-fatal render,
  no external resources, lightbox embedded, no broken anchors.
- Recent direct-removal commits deleted `_synthesis.qmd`, `synthesis_report`,
  `report_visuals`, and Overview body text. The current report order is QC ->
  microglia -> trajectory -> mechanism -> cross-modality.
- Main residual blocker is not "too few figures"; it is visible explanatory text,
  tables, and provenance chunks.

Docs:
- Quarto supports executable-code figure labels/captions/alt text via `label`,
  `fig-cap`, and `fig-alt`; cross-referenceable figures need `fig-` labels.
- Quarto HTML code folding exists, but the caption-only path should prefer
  `echo: false` / `include: false` for audit code rather than relying on folded
  visible code as a prose substitute.
- WCAG 2.2 non-text-content guidance still requires text alternatives. Therefore
  `fig-alt` is compatible with "no visible prose" and should be preserved or
  added for dense figures.

Sources:
- https://quarto.org/docs/authoring/figures.html
- https://quarto.org/docs/output-formats/html-code.html
- https://www.w3.org/WAI/WCAG22/Understanding/non-text-content

## Default Design

Default/selected by task = strict visible-path conversion.

Rules:
- No markdown paragraphs in report qmds after setup chunks.
- No visible tables. Use heatmaps, compact grids, text-in-figure tiles, or
  invisible build checks.
- No visible provenance chunks. Convert essential audit facts into small status
  panels or hide as `include: false` with source comments.
- Keep chapter headings because TOC/navigation need them; keep figure captions
  because the task explicitly allows captions.
- Captions state the claim/caveat once. Figure annotations carry numbers and
  logic. Avoid paragraph-like captions.
- Null/blocked/not-earned states remain visible as plotted status cells.
- Accessibility: every dense figure gains `fig-alt` or an equivalent alt strategy
  before close-out.

## Alternatives

Alternative A - one-page figure wall.
Pros: most literal "figures only"; easiest visual scan. Cons: loses chapter
audit order and likely requires generated-gallery plumbing.

Alternative B - remove captions too, push all text inside figures.
Pros: absolute minimum surrounding text. Cons: worse accessibility/crossref
semantics; higher risk of unreadable text in panels.

Alternative C - hide prose in collapsible audit.
Pros: fastest visible cleanup. Cons: violates the spirit of "no prose" because
the output still contains prose; reserve only for non-main-path source code if
needed.

Selected: default strict visible-path conversion.

## Steps

Each implementation step runs `scripts/check.sh` unless explicitly docs-only.

### S1 - Strict Inventory Gate [DONE 2026-07-03]

Work:
- Extend or wrap `scripts/prose_inventory.py` with a strict mode that reports
  counts by `kind` and fails if any visible non-heading/non-caption block remains.
- Add rendered-HTML QA for visible paragraphs, tables, `.cell-output-display`
  text-only outputs, and `cat()` provenance.
- Record current baseline in the plan after the stricter parser runs.

Acceptance:
- Current blockers are listed by qmd/line/kind.
- Gate command exists and fails on the current report before conversion.
- No report-source edits beyond inventory/gate plumbing.

Outcome:
- `scripts/prose_inventory.py --strict --html _report/index.html --summary-only`
  is the strict gate; current report fails red with source qmd/line blockers and
  rendered HTML line/kind blockers.
- Existing `scripts/check.sh` stays green; strict caption-only gate remains separate
  until S2-S4 remove the blockers.

### S2 - QC Chapter Conversion

Work:
- Convert QC prose/tables into figures: modality-shape grid, genotype/AOI/sample
  tally panels, bounds-status tile, and one optional compact audit figure.
- Remove `knitr::kable()` from `_qc.qmd`.
- Hide structural `stopifnot()` checks; rendered output is figure/caption only.

Acceptance:
- `_qc.qmd` has headings + figure chunks/captions only.
- QC facts still visible: modality shapes, GeoMx genotype tally, 16 sample key,
  genotype x batch completeness, metric bounds pass.
- Strict gate passes for `_qc.qmd`.

### S3 - Result Chapter Body-Prose Removal

Work:
- Remove one-sentence chapter openers and status-tail paragraphs from
  `_microglia.qmd`, `_trajectory.qmd`, `_mechanism.qmd`, `_crossmodality.qmd`.
- Convert the remaining visible provenance chunks to small status figures or
  hidden audit checks.
- Keep/null states visible in existing boards/captions: interaction DE
  under-powered, progression not supported, NF-kB discordant/not supported,
  Gsk3b not recovered, SpatialDecon abundance blocked, full CCC absent,
  run-index caveat.

Acceptance:
- Each result qmd has zero paragraph/list/table blocks outside captions/headings.
- No `echo: true` chunk renders plain provenance text.
- Captions remain <=24 words; target median <=12 words.

### S4 - Caption, Alt, and HTML QA

Work:
- Compress captions to claim/caveat microcopy.
- Add/verify `fig-alt` for dense panels without reintroducing visible prose.
- Render HTML and run DOM-level QA: figures/captions count, no visible body
  paragraphs, no visible tables, no warning/error markers, no external refs, no
  duplicate IDs, lightbox intact.

Acceptance:
- Rendered HTML main path = headings + figures + captions.
- Caption length rule passes.
- Accessibility text coverage recorded.
- Full `scripts/check.sh` green.

### S5 - Close-Out

Work:
- Adversarial claim-parity review: every removed prose claim is encoded visually,
  redundant, or intentionally dropped as nonessential audit detail.
- Update `.agent/map.md`, `.agent/memory.md`, `.agent/history.md` only for
  durable wiring/gotchas.
- Archive this plan and reset roadmap Active plan.

Acceptance:
- Strict gate green.
- Rendered HTML QA green.
- Claim-parity review has no accepted blockers.
- Roadmap/history record final counts and status.

## Claim Rules

- Visual compression may remove words, not evidence standards.
- Caption-only output still carries explicit null/blocked/not-earned states.
- Captions can be terse; they cannot overclaim.
- `fig-alt` supports accessibility and may be longer than visible captions.

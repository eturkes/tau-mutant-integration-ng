# Prose-to-Figures Reduction Plan

Opened: 2026-07-03. Task: the output report still has far too much prose; plan a
route to replace as much prose as possible with figures.

## Goal

Convert the final standalone HTML from figure-rich-but-still-narrative to
visual-first. Preserve the closed biological claims, warning-fatal render gate,
offline self-contained output, and compact-target report contract.

Success = materially less human-facing prose in `index.qmd` and `_*.qmd`, with
claims moved into figures, schematic panels, matrices, caption microcopy, and
visual status encodings. The implementation phase must measure prose before and
after; the initial target is at least a 40% reduction in prose-only word count,
raised after S1 inventory if more can be removed safely.

## Research Digest

Local state:
- Current report source rough word count is 15,066 across `index.qmd` and the
  five result includes, but this includes code chunks. S1 must compute a
  prose-only baseline excluding fenced code, YAML, setup chunks, and hidden helper
  code.
- Figure expansion already delivered 42 captioned figure blocks, lightbox support,
  and compact per-chapter figure targets. This plan is not "add more plots" by
  default; it is "make figures carry the narration".
- Report qmds already obey the cheap-render contract: load compact report/figure
  targets, set `options(warn=2)`, and render through `scripts/check.sh`.
- The archive confirms v1's failure mode: long report bodies, large structured
  ledgers, many prose verdicts, and capstone narrative. Useful salvage is visual
  form only: convergence matrices, volcanoes, concordance heatmaps, pathway
  grids, decomposition panels, and standalone publication PNGs.

Current docs/tooling:
- Quarto figure and cross-reference docs support figure labels with `fig-`
  prefixes; keep hyphenated labels and no underscores.
- Quarto diagram docs support Mermaid and Graphviz directly, which is useful for
  claim-flow schematics without a new dependency.
- Quarto HTML tabsets/panels exist, but hiding prose in tabs is a fallback, not
  the main reduction. Prefer visual replacement over prose relocation.

Sources consulted:
- https://quarto.org/docs/authoring/figures.html
- https://quarto.org/docs/authoring/cross-references.html
- https://quarto.org/docs/authoring/diagrams.html
- https://quarto.org/docs/output-formats/html-basics.html

## Default Design

Default = aggressive inline visual conversion over the existing report, with no
new biological inference and no new heavyweight report-time reads.

Principles:
- Every paragraph earns a concrete disposition: delete, convert to a figure,
  convert to a schematic, compress into caption microcopy, move to a collapsed
  audit note, or keep with justification.
- Figures carry claims as data ink: small multiples, signed status grids,
  evidence matrices, schematic flows, and caveat glyphs. Avoid replacing prose
  with paragraph-like tables unless a table is truly the most compact visual
  object.
- Captions become microcopy: one factual sentence plus, only when needed, one
  caveat clause. Long explanatory captions count as failed prose reduction.
- Unsupported/null findings stay visible as plotted states, not buried in text.
- Methods and caveats remain auditable, but the main reading path becomes visual.

Likely conversions:
- Overview + synthesis: replace answer paragraphs with a visual abstract,
  claim-source matrix, unsupported/unearned status grid, and one compact spine
  schematic.
- QC: keep as compact diagnostics; convert explanatory text to axis labels,
  annotations, and figure subtitles.
- Microglia: promote composition/DE/dropout caveats into annotated panels; delete
  repeated paragraphs already encoded by existing figures.
- Trajectory: merge composition-vs-progression logic into the decomposition
  figure and a flow schematic; captions carry the "composition, not supported
  progression" result.
- Mechanism: encode Myc/NF-kB/Gsk3b as a signed mechanism status board plus
  pathway/TF/kinase panels; strip prose that restates figure cells.
- Cross-modality: use modality x contrast matrices and run-index loss encodings
  as the primary read; keep only the minimum text needed for blocked SpatialDecon
  and bulk caveats.

## Alternatives

Alternative A - Standalone visual abstract plus shorter chapters.
Pros: strongest first-screen improvement; easier to review. Cons: leaves much
chapter prose in place unless followed by the default inline pass.

Alternative B - Publication PNG gallery.
Pros: highest polish and easiest external sharing. Cons: reintroduces tracked
generated assets and script-refresh burden; weaker target-synchronisation than
in-knit compact targets.

Alternative C - Hide prose in collapsible audit sections.
Pros: fastest visible slimming. Cons: prose still exists in the output file and
does not satisfy the user's "replace with figures" direction except as a last
resort for methods details.

Default recommendation: choose the aggressive inline conversion. It directly
targets the complaint, reuses the current compact figure infrastructure, and
keeps the report as one warning-gated offline artifact.

## Steps

Each step is one closing unit. Run `scripts/check.sh` for implementation steps.
S1 may be docs/script-only but still needs at least the prose-inventory command
and a clean `git diff` review.

### S0 - Route Gate [OPEN]

Work:
- Confirm default aggressive inline conversion or choose an alternative.
- Set the provisional prose-reduction floor. Default starts at at least 40%
  prose-only reduction; S1 can raise it after the baseline inventory if more
  low-risk removals are visible.

Acceptance:
- User route selected.
- Roadmap updated from S0 to S1 with the chosen route and provisional reduction
  target.
- No report source edits before the route is selected.

### S1 - Prose Inventory and Replacement Manifest

Work:
- Add a lightweight prose-inventory script or command documented in the plan.
  Count prose-only words/paragraphs per qmd while excluding YAML, code chunks,
  generated code, and comments.
- Classify every human-facing block in `index.qmd` and `_*.qmd`:
  `delete`, `figure`, `caption`, `schematic`, `collapsed_audit`, or `keep`.
- Produce a compact replacement manifest in this plan or a generated TSV under a
  gitignored scratch path, with chapter, block id, disposition, and target figure
  slot.

Acceptance:
- Baseline prose-only counts by chapter are recorded.
- At least 80% of prose blocks have non-`keep` dispositions unless a chapter
  review justifies otherwise.
- The selected reduction target is concretely stated before S2.

### S2 - Visual Grammar and Data Contract

Work:
- Extend `R/figures.R` and compact figure targets only where existing slots cannot
  carry the replacement.
- Add visual-first slots for: report spine schematic, claim-source matrix,
  unsupported/unearned status grid, caveat/status glyph data, and chapter-level
  evidence boards.
- Prefer Mermaid/Graphviz for schematic flows that do not need data joins; prefer
  ggplot/patchwork for target-derived panels.
- Add tests for any new target slots, finite geom inputs, and manifest coverage.

Acceptance:
- New visual slots cover every S1 `figure`/`schematic` disposition.
- Compact target sizes remain qmd-safe; no raw modality or heavy Seurat report
  reads are introduced.
- `tests/test_figures.R` or a focused new test locks slot presence and finite
  plotted values.

### S3 - Synthesis and Overview Conversion

Work:
- Rewrite `index.qmd` and `_synthesis.qmd` so the first reading path is visual:
  visual abstract, claim-source matrix, status grid, and minimal connective text.
- Replace the compact evidence table if a matrix/board can carry the same audit
  content more tersely.
- Keep one short answer sentence if needed; remove repeated chapter previews.

Acceptance:
- Prose-only word count for Overview + synthesis drops by the selected target or
  better.
- Claims and caveats remain source-target-derived, not hand-waved by diagrams.
- Forced render is warning-clean.

### S4 - Result Chapter Conversion

Work:
- Apply the S1 manifest to `_microglia.qmd`, `_trajectory.qmd`,
  `_mechanism.qmd`, and `_crossmodality.qmd`.
- Replace explanatory paragraphs with annotated multi-panel figures, status
  boards, and schematic logic panels.
- Preserve explicit nulls/blocked states: no supported progression beyond
  composition, NF-kB attenuation not supported, Gsk3b not recovered,
  SpatialDecon abundance blocked, full CCC absent, and bulk run-index caveats.

Acceptance:
- Each result chapter meets or beats the selected prose-reduction target unless
  an accepted review note records why not.
- Figure captions are short and factual; no caption becomes a paragraph surrogate.
- All figure labels remain hyphenated `fig-*`; no duplicate labels.
- Full `scripts/check.sh` green.

### S5 - Visual QA and Close-Out

Work:
- Run full rendered HTML QA: figure count, captions, labels, lightbox assets,
  section anchors, no external resources, no warning/error markers.
- Run the prose-inventory command again and record before/after counts.
- Adversarially review claim parity: every removed prose claim is either encoded
  visually, redundant, or deliberately deleted as nonessential.
- Update `.agent/map.md`, `.agent/memory.md`, and `.agent/history.md` only for
  durable wiring/gotchas; archive this plan and reset roadmap Active plan.

Acceptance:
- Final before/after prose-only count shows the selected reduction target was met
  or a specific accepted blocker explains the residual prose.
- Full `scripts/check.sh` green.
- Rendered HTML QA green.
- Close-out review finds no unsupported visual overclaim.

## Claim Rules

- Visual compression may reduce words, not evidence standards.
- A figure cannot imply support that the underlying target does not support.
- Null, blocked, and not-earned states must have visible encodings.
- Captions and labels use British English, ASCII punctuation, and concise
  factual wording.

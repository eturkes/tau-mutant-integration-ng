# Roadmap - forward index (read first each session)

Posture · cohesive story · active plan · phased backlog. Companions in `.agent/`:
`memory.md` (standing contract), `map.md` (wiring, grows with the build),
`history.md` (new-project decision digests), `archive_digest.md` (v1 reference -
mine for promising threads; full v1 code on branch `archive`).

## Posture: BUILD (fresh)
v1 (23 rmd / 11 arcs D-O / 92-row ledger / 3 prose-heavy reports) is archived on
branch `archive`. `main` rebuilds from scratch: leaner pipeline, ONE cohesive
story, minimal prose, project-local env, no Docker. Re-derive value, drop bloat.

## Cohesive story (the spine - proposed, revisable)
Data = 4 mouse AD genotypes in a 2x2 (tau: MAPTKI ~ tau-KO vs P301S mutant; amyloid:
-/+ NLGF) x 4 modalities (snRNAseq, GeoMx spatial, bulk proteome, bulk phospho).
Question = how amyloid reshapes microglia under each tau background. Divergence =
interaction (NLGF_P301S - P301S) - (NLGF_MAPTKI - MAPTKI).
Durable findings mined from v1 (the headline to rebuild around):
1. Amyloid drives a microglial homeostatic->DAM activation programme (DAM / MG-M2 /
   antigen-presentation); microglia-led, corroborated across modalities.
2. Mutant tau MODULATES the amyloid response (the interaction) - clearest as a
   SYNERGISTIC acceleration along the activation trajectory (v1's one positive,
   significant orthogonal interaction).
3. Mechanism converges on Gsk3b / Myc signalling; tau attenuates amyloid-driven NF-kB.
4. Secondary axis: amyloid-linked synaptic suppression + TREM2 / APP clearance.
Streamlined OUT (v1 bloat): the 11-arc ledger + contest machinery, the margin-neutral
corroboration arcs (SCENIC, spatial-decon, celltype-specificity, gene-level dynamics),
the human-validation layer, the capstone convergence matrix, the heavy prose.

## Active plan: `p0_foundations_plan.md` - P0 Foundations (open; **S4 IN PROGRESS**). Stack: targets + rv + uv + project-local Quarto, P3M-pinned. Reports = ONE offline self-contained HTML (standalone Quarto doc + `{{< include _section.qmd >}}` + `theme.scss` + LOCAL IBM Plex), NOT a book (multi-file + sibling nav warnings under embed-resources would trip the S5 zero-warning gate). S1-S3 done. **S4 checkpoint committed** (gate-clean: offline render, 0 error/0 warning, bounds pass, tests green): book -> standalone `index.qmd` + `{{< include _qc.qmd >}}`; `theme.scss` reduced to working COLOURS (#B0344D primary/link, #3F5A6B code - verified inlined); `_qc.qmd` bounds hardened (16x16 genotype-batch bijection); `R/plot.R`; `_quarto.yml`/`_targets.R` reworked; woff2 gitignored. **NEXT SESSION (font wiring -> close S4; now SIMPLE):** codex review CORRECTED the premise -- IBM Plex DOES work via theme.scss (the prior "doesn't" was a URL-encoded-output measurement artifact). Add the 9 `@font-face` (relative `url("assets/fonts/<n>.woff2")`) + `$font-family-*` vars + body/headings/code rules to theme.scss; Quarto base64-inlines the woff2 under embed-resources (no build script needed for inlining). Open decision: commit the 9 woff2 (~200KB, KISS) vs keep gitignored + a sha256-pinned `scripts/build-fonts.sh` fetch; list them in `extra_files=` either way. Then re-render to PROVE inlining via a DECODED check (urllib.unquote the data:text/css; raw greps mislead) while offline + 0-warning hold + `tests/test_plot.R` + scrub residual book/_brand refs + map.md + mark S4 done. Detail: the plan's `S4 STATUS` block.

## Backlog - phased build (each phase = closeable increments; mine archive_digest per phase)
- P0 Foundations: project-local env (rv for R + uv .venv for Python), shared
  helpers (io / design+contrasts / plot theme), data load + QC sanity, 2x2
  factorial + 5 contrasts, concrete quality gate.
- P1 snRNAseq microglia core: QC, microglia subset + substates (homeostatic / DAM /
  IFN / proliferative), DE (pseudobulk + single-cell) across contrasts -> the
  amyloid->DAM activation programme + the gene/pathway interaction.
- P2 Interaction trajectory: activation pseudotime (homeostatic->DAM); test amyloid
  advance + tau x amyloid progression synergy, decomposed composition vs progression.
- P3 Mechanism: focused pathway/module survey; TF (decoupleR / CollecTRI) + kinase
  (decoupleR / OmniPath) -> Gsk3b / Myc; NF-kB attenuation check.
- P4 Cross-modality: GeoMx spatial DE (+ light deconvolution if it earns it),
  proteome + phospho DE, CCC for the synaptic/clearance axis, integrated divergence view.
- P5 Synthesis: ONE lean report - cohesive narrative + compact evidence table (no
  ledger / contest machinery).

## Ledger (trajectory)
- 2026-06-29 archived v1 -> branch `archive`; opened fresh orphan `main`; reset
  `.agent` docs + Claude config; reframed history as `archive_digest.md`; drafted
  this streamlined phase plan.
- 2026-06-29 S4 checkpoint committed (architecture pivot; codex-reviewed): report = ONE
  offline self-contained HTML -> Quarto book dropped for a standalone `index.qmd` +
  `{{< include _qc.qmd >}}`; `theme.scss` cut to working colours; `_qc.qmd` bounds hardened
  to a 16x16 design bijection; `R/plot.R`, `_quarto.yml`, `_targets.R` reworked; render
  offline + 0 error/0 warning, tests green. Remaining (next session, closes S4): IBM Plex
  via theme.scss @font-face (corrected by the codex-review entry below), then `tests/test_plot.R` +
  scrub + map.md + mark S4 done.
- 2026-06-29 codex review (S4 checkpoint): CORRECTED the font conclusion -- IBM Plex DOES work
  via theme.scss (@font-face relative url -> Quarto base64-inlines the woff2 under embed-resources);
  the prior "Bootstrap vars + url() both fail" was a measurement artifact (theme CSS embeds URL-
  ENCODED, so raw `IBM Plex`/`data:font/woff2` greps read ~0). Simplified the remaining S4 font
  wiring (no build script for inlining; open: commit woff2 vs sha256 fetch). Also hardened the GeoMx
  QC check (!anyNA + setequal vs nlevels==4), dropped the "balanced" overclaim (index.qmd), scrubbed
  stale "Quarto book"/`_brand.yml` plan refs.

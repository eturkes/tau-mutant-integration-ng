# Capstone synthesis plan (arc P — consolidation, not a 12th evidence arc)

**Goal:** Consolidate the eleven completed evidence arcs (D→O) + the four base
modalities + the human translational layer into ONE legible integrated synthesis:
an **evidence-layer × biological-axis convergence matrix** (the central new
artifact) + a unified figure + narrative, surfaced in the curated `summary.Rmd`.
PURE re-aggregation/presentation of existing ledger + verdict TSVs — **no new
biological computation**. Project-closing capstone. Does **NOT** modify the locked
`analysis.Rmd` §17 ledger, the 11-entity set, the grade defs, or any contest
margin (read-only aggregation throughout).

User chose "Capstone synthesis" (2026-06-15) over MOFA+ / meta-analysis / ROSMAP
(ROSMAP blocked on a RADC/Synapse DUA). Ruled out this session: a sex-modifier
arc — sex is perfectly aliased with batch (batch 01/03=male, 02/04=female; 2M/2F
per genotype), so the locked batch term already absorbs it (non-estimable; also
confirms the locked design is not sex-confounded).

## STATUS
- [DONE] C1 — Decision gate RESOLVED (2026-06-15): placement = **summary.Rmd ONLY** (sealed §17 untouched); figure form = **polished standalone PNG via the `scientific-visualization` skill**, embedded in summary.Rmd.
- [DONE] C2 (2026-06-15) — `scripts/build_capstone_synthesis.R` → 2 audit-sourced TSVs (`capstone_convergence_matrix.tsv` 39 cells×9 cols; `capstone_contest_summary.tsv` 3×10); runtime-validated, smoke-tested. Two curation extensions logged in the C2 body (added `mixed` cell-token + `band` column).
- [DONE] C3 (2026-06-15) — `scripts/build_capstone_figure.py` (`scientific-visualization` skill, matplotlib in `.venv`) → `storage/figures/capstone_convergence.png` (13×8.6in @300dpi). Two panels: (A) 12-row mouse convergence heatmap + separated Human band, (B) the 3 contest-margin bars. Committed asset path made tracked via a new `.gitignore` re-inclusion for `storage/figures/*.png`.
- [DONE] C4 (2026-06-15) — `storage/cache/summary_chunks/capstone_convergence.R` (one chunk: `include_graphics` PNG + `results='asis'` kable of `capstone_contest_summary.tsv`) + new "The integrated convergence model" section/panel/title in `build_summary_rmd.R`; regenerated `summary.Rmd` (22→23 panels). **USER-DIRECTED EXTENSION ("sync + full re-audit", 2026-06-15):** the summary was last built at H7 (Jun 4), predating arcs K–O, so its ledger-backed prose had drifted. Re-audited EVERY panel against live data — mtime proof: all `summary_*.rds` extracts are dated ≤ Jun 4 and the arc-K/L/M/N/O caches (scenic/spatial/trajectory/specificity/dynamics, Jun 5–8) are read by NO summary panel; the Gsk3b headline was empirically re-derived from `kinase_activity_decoupler.rds` (ULM +4.32 / FDR 0.00152 / consensus 1.8e-4 — exact match). Only the 2 ledger-driven panels were stale, both synced read-only: `biomodel_contest_verdicts` (interaction margin 53→55) and `crossaxis_theme_ranking` (lead flip — T-Inflammation 31 (6 Strong) now > T-Synergy 30 (9 Strong) > T-Compartment 17 > T-Tau-attenuates 12 > Hyp-0 7).
- [DONE] C5 (2026-06-15) — Knit clean: 0 error / 0 warning (Rscript-verified, since `summary.html` is deny-Read and `grep` on it is blocked); capstone section + figure (19 embedded PNGs / 21 `<img>`) + contest table all confirmed present, both syncs in the HTML. chown'd `summary.{Rmd,html}` + the chunk file + builder; single scoped commit; refreshed `map.md` + `DIGEST.md`; archived plan as `completed/capstone_synthesis_plan_2026-06-15.md`.

## Locked design principles (do not re-litigate; they ARE the synthesis value-add)
1. **Dual representation, always.** (a) **Contest view** = the ledger net-support
   + final margins (what *moved* the verdicts; from
   `biological_model_{adjudication,contest_verdicts}.tsv`). (b) **Landscape view**
   = each arc's *own* verdict + signed stance per axis (the full convergence,
   from the 8 per-arc verdict TSVs + chapters). Rationale: arcs K–O are
   deliberately margin-neutral (1–2 ledger rows each), so a pure ledger-net
   matrix would make five whole arcs look trivial. The landscape view is the
   point — it shows orthogonal corroboration the contest arithmetic discounts.
2. **Every matrix cell carries a `source`** (verdict-TSV path + value, or ledger
   `claim_id`) — mirror the §17 ledger audit-trail discipline; zero uncited
   cells.
3. **Nulls and non-convergence are first-class.** The interaction is
   static-null at matched power (N), composition-null (L), human-interaction
   null under FDR — these render as explicit null/− cells, never hidden. State
   the K–O margin-neutral pattern plainly (corroboration without verdict-movement).
4. **Human stays a SEPARATE band** (evidence-class-incommensurable per H7) —
   directional translational consistency, never counted into the mouse contests.
5. **Read-only.** No edit to the ledger / 11-entity set / grade defs / margins /
   `analysis.Rmd` §17. The capstone aggregates; it never re-adjudicates.

## Convergence-matrix spec (locks C2/C3 concretely)
**Rows (13 evidence classes) → verdict source:**
| row | arc | verdict source |
|---|---|---|
| Expression DE (4 modalities) | base 02b/04/05/06/06b/07 | `integration_table.rds` + H-phase expression ledger rows |
| Pathway/module | D (rmd/12) | `pathway_survey_unified_leaderboard.tsv` + ledger `pathway` |
| TF activity | E (rmd/13) | `tf_activity_verdict.tsv` + ledger `tf` |
| Kinase activity | F (rmd/14) | `kinase_activity_verdict.tsv` + ledger `kinase` |
| CCC / ligand-receptor | G (rmd/15) | `ccc_lr_verdict.tsv` + ledger `lr` |
| NF-κB attenuation | I (rmd/17) | ledger I-001..I-006 + `per_state_*` TSVs |
| Causal topology | J (rmd/19) | `causal_network_verdict.tsv` + ledger J-001/J-002 |
| SCENIC regulons | K (rmd/20) | `scenic_verdict.tsv` + ledger K-001/K-002 |
| Composition | L (rmd/21) | `spatial_decon_verdict.tsv` + ledger L-001/L-002 |
| Dynamics/progression | M (rmd/22) | `trajectory_verdict.tsv` + ledger M-001/M-002 |
| Specificity | N (rmd/23) | `celltype_specificity_verdict.tsv` + ledger N-001/N-002 |
| Dynamics-DE | O (rmd/22 sub) | `trajectory_dynamics_*.tsv` + ledger O-001 |
| Human cross-species | rmd/18 | `summary_human_validation.rds` (SEPARATE band) |

**Cols:** `amyloid_activation`, `synaptic_suppression`, `interaction_metabolic`
(the 3 contests). `cross_axis` carried as a note column, not a 4th display col.

**Cell vocabulary (locked):** signed direction × grade ∈
{`strong+`,`mod+`,`sugg+`,`null`,`sugg−`,`mod−`,`strong−`} OR a meta-tag
{`reframes`,`n.a.`} where the arc's verdict is not a directional axis call (e.g.
N specificity = `reframes`; an arc that never probed an axis = `n.a.`). Grade
inherits the ledger 3-tier (Strong/Moderate/Suggestive) where a ledger row backs
the cell; else the arc verdict's own confidence word.

**Output TSVs (storage/results/):**
- `capstone_convergence_matrix.tsv` — long: `evidence_class, arc, axis,
  cell_value, sign, grade, source, note`.
- `capstone_contest_summary.tsv` — `contest, hyp_A, hyp_B, winner, margin,
  n_strong, n_moderate, n_suggestive, n_layers, n_modalities` (re-surfaced from
  adjudication + contest_verdicts = the contest view).

## Step bodies

### C1 — Scope & matrix-design decision gate (GATE; needs user)
Resolve two genuine forks before any build (design principle #1 dual-view is
already LOCKED above, not a fork):
- **Placement.** Default: `summary.Rmd` ONLY (honours the chosen "extend
  summary.Rmd, don't add a 12th arc" framing + keeps sealed §17 untouched).
  Alt: also add a read-only convergence panel to `analysis.Rmd` (new §25 before
  session-info, OR appended to §17.5) so the main pipeline doc carries the
  capstone too — costs a full `analysis.Rmd` re-knit (~3 min) and edits near
  sealed §17 territory (read-only aggregation, but adjacency risk).
- **Figure form.** Default: in-knit `ggplot2` convergence heatmap (matches every
  other `summary.Rmd` panel; one shared knit session, fully reproducible). Alt:
  a polished standalone publication figure via the `scientific-visualization`
  skill (Python), embedded as a pre-built PNG — higher polish for the
  outward-facing capstone, but a non-reproducible-in-knit static asset + an
  extra build dependency.
Present both forks with defaults + the reasoned alternative, then
`AskUserQuestion`. Record answers in this STATUS block before C2.

**RESOLVED 2026-06-15:** placement = **summary.Rmd only** (default — sealed §17
left untouched, no `analysis.Rmd` re-knit). Figure form = **polished standalone
PNG via the `scientific-visualization` skill** (the alternative — publication
polish for the project's closing artifact). Consequence: C3 builds a static PNG
asset (not in-knit ggplot); C4 embeds it; the figure is NOT regenerated inside
the knit, so the committed PNG is the deliverable (re-run the skill to refresh).

### C2 — Build capstone synthesis data
Write `scripts/build_capstone_synthesis.R` (idempotent; no new R/ helper — keep
the aggregation script-local, KISS). It reads the ledger + adjudication +
contest_verdicts + the 8 per-arc verdict TSVs + the pathway leaderboard +
`integration_table.rds` + `summary_human_validation.rds`, and emits the two TSVs
in the spec above. The arc×axis cell encoding is a CURATED mapping (the verdict
TSVs are not all natively "signed per axis"), each cell anchored to a cited
`source` per principle #2; the per-arc signs come straight from the DIGEST verdict
reads (amyloid+ strong/multi-layer; synaptic− at NLGF; interaction localised to
dynamics M, null at static-matched N / composition L). Smoke-test against the
live TSVs via `Rscript -e` (shape + every cell has a non-empty source + signs
match the DIGEST) BEFORE any knit. chown outputs.

**DONE 2026-06-15:** Built `scripts/build_capstone_synthesis.R` (idempotent,
script-local, no new R/ helper). It loads the contest view
(`biological_model_{contest_verdicts,adjudication,claims_ledger}.tsv`), the 8
per-arc verdict TSVs, the pathway leaderboard, the arc-O interaction TSV, and the
two RDS inputs, then emits `storage/results/capstone_{convergence_matrix,
contest_summary}.tsv`. Smoke-tested via `Rscript scripts/build_capstone_synthesis.R`
(builds clean; the round-trip re-read confirms 39×9 + 3×10 with every `source`
non-empty).

*Convergence matrix* = 13 evidence classes × 3 axes = 39 long-format cells, a
CURATED arc×axis mapping whose signed calls come from the DIGEST verdict reads.
Amyloid column is the strongly-convergent + axis (8 positive cells) with one
`strong-` (the arc-I NF-κB attenuation = the tau-modifier behind Hyp-1B). The
interaction column is deliberately heterogeneous and the synthesis value-add:
mechanism-signed (TF `mod-` Myc-led / kinase `mod+` Gsk3b-led), recovered at
causal+dynamics+dynamics-DE (`sugg+`), but `null` at SCENIC (expression-floor) /
composition (L) / human-under-FDR, `mixed` at expression+pathway+LR, and
`reframes` at specificity (N, not microglia-specific). Synaptic column is the
expression-down + LR-clearance story (Hyp-2B): `mod-` at expression/pathway,
`null` non-findings at TF/kinase/causal, `mod+` at LR (the decisive TREM2/APP
clearance arc), `n.a.` elsewhere. Every cell carries a cited `source` (verdict-TSV
path + value or ledger claim_id) per principle #2; nulls/reframes are first-class
per principle #3.

*Contest summary* = the 3 contests re-surfaced verbatim from the adjudication +
contest_verdicts outputs (winner / margin / n_strong / n_moderate / n_suggestive),
plus two light descriptive aggregates over the ledger rows supporting each winner:
`n_layers` = distinct `layer` count, `n_modalities` = max `n_replicates_in_modalities`.
Margins re-read 18 / 12 / 55 (Hyp-1B / Hyp-2B / Hyp-3B); a build-time guard hard-
asserts those three values so the capstone fails loudly if the sealed §17 verdicts
ever drift.

*Read-only / drift guards:* the script touches nothing in the ledger / 11-entity
set / grade defs / margins / analysis.Rmd §17. It anchors the most load-bearing
cells to a live substring in their source TSV (Spi1/Myc, padj 0.002, rank 775,
8.9e-07, +2.46, Microglia 625, …) and recomputes the arc-O "110 of 11,466
FDR<0.10" from `trajectory_dynamics_interaction.tsv$adj.P.Val` — a future
verdict-TSV edit that drops a cited value fails the build. A directional sanity
block asserts the amyloid axis stays predominantly +, the NF-κB cell stays
`strong-`, and the interaction axis keeps ≥3 explicit null/reframes cells (so an
LLM "everything converges" prior cannot silently inflate the matrix).

*Two curation extensions (within C2's "curated mapping" latitude; flagged for the
record):* (1) added a `mixed` cell-token for genuinely sign-divergent cells
(modalities / LR pairs / difference-of-differences disagree) — forcing them to a
signed token over-states convergence, forcing them to `null` hides that the layer
IS significant; the ledger itself uses `direction="mixed"`. (2) added a `band`
column (mouse|human) so C3 can hold the human row in its own separated band
(principle #4 made machine-readable) — schema is the plan's 8 cols + `band` = 9.
`sign` carries two documented framings (effect-direction arcs vs
existence/recovery arcs); per-cell `note` makes each explicit so no cell relies on
the frame alone.

### C3 — Build the convergence figure (publication PNG via scientific-visualization)
Invoke the `scientific-visualization` skill (Python/matplotlib) to render a
publication-grade figure FROM the C2 TSVs: arc rows × 3 axis cols, diverging
signed palette (− blue / 0 grey / + red), grade encoded (alpha / border /
hatch), the human row offset into its own separated band; a companion panel of
the three contest margins (18 / 12 / 55) so the contest view sits beside the
landscape view. Confirm a stable committed asset path first (e.g.
`storage/results/capstone_convergence.png` or a `figures/` dir) and check the
`.gitignore` does not exclude it — the PNG is the committed deliverable (not
regenerated in-knit). Verify the PNG opens + matches the TSV signs before wiring.

**DONE 2026-06-15:** Built `scripts/build_capstone_figure.py` (idempotent;
matplotlib 3.10.9 in `.venv`; reads only the two C2 TSVs, computes nothing) →
`storage/figures/capstone_convergence.png` (13×8.6 in @300 dpi, 4569×2555 px).
Visually verified against the matrix TSV (all 39 cells match). *Asset path:* chose
`storage/figures/` over a top-level `figures/` to keep the project's one-root-dir
convention, and added a `.gitignore` re-inclusion (`!/storage/figures/` →
`/storage/figures/*` → `!/storage/figures/*.png`) mirroring the `storage/notes`
exception so the PNG is tracked while everything else under `storage/` stays
ignored; `git check-ignore` confirms the PNG is not excluded.

*Panel A (landscape):* a categorical convergence heatmap, 12 mouse evidence-class
rows × 3 axis columns, narrative-arc order (base, D…O), with the single Human row
held below a dashed divider in its own faint-purple band (italic purple label;
H7 separation made visual). Cells coloured by sign × grade on a colourblind-safe
diverging ramp (RdBu endpoints: + dark→light red, − dark→light blue), with four
meta categories outside the ramp — purple = `mixed` (sign-divergent), grey =
`null`, near-white = `n.a.`, gold + `////` hatch = `reframes`. Every cell ALSO
prints its token (Strong +/Mod −/…) so grade + sign survive grayscale /
colourblind reading (redundant encoding; the dark-red/dark-blue luminance tie is
broken by the text). *Panel B (contest):* the three sealed margins as horizontal
bars (55/18/12, largest at top), each annotated with the winning hypothesis +
its S/M/Su support counts + layer count — the contest view beside the landscape
view (design principle #1). A 10-entry legend + a 5-line audit caption (every
cell cites a source; the interaction's mechanism-signed-yet-null-in-
composition/regulons/human pattern stated plainly) close it.

*Layout iteration:* first render had the title block colliding with the panel
headers and Panel B's y-labels bleeding into Panel A; fixed by lowering the axes
(top=0.79), spacing the suptitle/subtitle/panel-labels on distinct baselines, and
replacing Panel B's left y-tick labels with per-bar above/below annotations
(decoupling B from A). The figure is NOT regenerated in-knit — re-run the script
to refresh; C4 embeds the committed PNG via `knitr::include_graphics`.

### C4 — Narrative + wire into summary.Rmd
New chunk body `storage/cache/summary_chunks/capstone_convergence.R` that
**embeds the pre-built PNG** (`knitr::include_graphics(...)`) + renders the
contest-summary table from the C2 TSV (`results='asis'` kable) + a `panels` entry
in `scripts/build_summary_rmd.R` (section header
"The integrated model" or similar, near the end before any closing panel; British
-English prose). Narrative states the integrated 3-axis model as the synthesis of
all arcs: (a) amyloid-driven DAM activation = the strongest, most multi-layer-
convergent axis; (b) the tau×amyloid interaction = carried by progression DYNAMICS
(M, positive) NOT static gene-level expression (N, null at matched power) NOR
composition (L, null) — the load-bearing nuance; (c) the K–O margin-neutral
pattern = orthogonal corroboration that deliberately does not move decided
contests; (d) human = directional translational consistency, kept separate.
Honest about diminishing returns / project closure. Re-run `build_summary_rmd.R`
to regenerate `summary.Rmd` (note: chunk dir is gitignored — chunk body must be
present locally to regenerate; the committed deliverable is `summary.Rmd` itself).

### C5 — Knit-verify, commit, close
`Rscript -e 'rmarkdown::render("summary.Rmd", quiet = TRUE)'` (+ `analysis.Rmd`
only if C1 chose the dual-placement alt). Verify `grep -c 'class="error"'` and
`'class="warning"'` are 0 on the rendered HTML. chown new files +
root-owned knit outputs. Single scoped commit. Refresh `map.md` (add the
capstone builder + TSVs to the caches/scripts tables; note the summary panel) and
`DIGEST.md` (add a capstone arc-P digest entry); archive this plan as
`completed/capstone_synthesis_plan_<YYYY-MM-DD>.md`.

## Execution model
C1 is the gate — resolve with the user, record answers here, then STOP for a
fresh session to implement (C2 is design-heavy: new builder + new TSV schema).
C2→C3 may chain in one session (data + figure are a tight pair). C4→C5 chain in a
closing session (narrative + wire + knit + commit). Per-step loop per
`/session-prompt`: smoke-test cache reads before knitting; mark DONE + completion
note; chown; knit-verify 0 error/0 warning; single scoped commit. End cleanly at
the gate or on low context.

## Anti-anchoring guardrails
- **Read-only.** Never edit the ledger / 11-entity set / grade defs / contest
  margins / `analysis.Rmd` §17. The capstone aggregates; it never re-adjudicates.
- **Cite every cell** (principle #2). A matrix cell with no `source` is a bug.
- **Nulls are first-class** (principle #3) — show the interaction's static/
  composition/human nulls explicitly; do NOT let an LLM "everything converges"
  prior inflate the matrix. Over-stating convergence is the primary failure mode.
- **Human separate** (principle #4); no mouse-contest counting of human evidence.
- **No pre-privileged axis** in the narrative; alphabetise hdWGCNA modules
  (never MG-M3 first); state thresholds before applying.
- **Standing bans:** never reintroduce MSigDB Hallmark, the OXPHOS confirmatory
  arc, or a specificity-null / expression-matched-random framework.
- British English throughout; chown new files `rstudio:rstudio`.

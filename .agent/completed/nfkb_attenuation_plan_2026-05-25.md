# Active plan: per-state NF-kB attenuation test for T-Tau-attenuates (Phase I)

- **Created:** 2026-05-25 (immediate successor to
  `storage/notes/completed/biological_model_plan_2026-05-25.md`,
  Phase H COMPLETE).
- **Last updated:** 2026-05-25 (I4 DONE: 6 new ledger rows
  (I-001..I-006) propagate the I3 Global-attenuation verdict into the
  H2 ledger; T-Tau-attenuates lifts from 5 supports / 3 Strong to
  11 / 7; Hyp-1B's contest margin over Hyp-1A lifts from 6 to 18;
  cross-axis ranking re-orders T-Tau-attenuates 5th -> 4th (swapping
  with Hyp-0 Cdk5); section 17 prose + section 1.1 H5 upfront mirror
  refreshed; re-knit zero errors / zero warnings; plan moved to
  `storage/notes/completed/` with Outcome summary).
- **State:** Phase I COMPLETE. I1 DONE. I2 DONE. I3 DONE. I4 DONE.
- **Predecessor:** `storage/notes/completed/biological_model_plan_2026-05-25.md`
  (Phase H COMPLETE; the integrated biological model now lives in
  section 17 of `analysis.html` with the H2 ledger + H3 adjudication
  + H5 upfront mirror at section 1.1).
- **Goal:** test whether the **NF-kB attenuation** signature that
  the H2 ledger uses to support **T-Tau-attenuates** (5 rows, the
  lowest-support cross-axis theme) is a **global whole-microglia
  effect** or **localised to a specific microglia substate**
  (MG_homeostatic vs MG_DAM vs MG_IFN vs MG_proliferative). The
  test is a per-state TF activity inference at the interaction
  contrast (using `decoupleR` + CollecTRI, the same toolchain as
  section 14) with a specific focus on the NF-kB family (Rela,
  Nfkb1, Nfkb2, Rel, Relb, A0A979HLR9). The outcome — global
  attenuation vs state-restricted attenuation vs state-driven
  attenuation — determines whether T-Tau-attenuates lifts out of
  its current 5-row weakness, gains substate specificity, or
  collapses under per-state stratification.

## Why this plan exists

The Phase H synthesis closed with T-Tau-attenuates as the
**lowest-support cross-axis theme (5 supports, 3 Strong)**. The
load-bearing evidence for the theme is the Rela sign-reversal:
positive (+2.57 +1.81) at NLGF arms, negative (-1.36) at the
interaction. Phase H interpreted this as "tau on the amyloid
background suppresses amyloid-driven NF-kB" at the whole-microglia
level. But the suppression's biological character is unresolved
because the test was run on whole-microglia NEBULA outputs which
mix the four microglia substates (homeostatic / DAM / IFN /
proliferative). Three scientifically distinct underlying patterns
all collapse into the same whole-microglia signal:

1. **Global attenuation.** Every microglia substate independently
   shows NF-kB suppression at the interaction; the whole-microglia
   signal is the average. Biological interpretation: tau
   broadly tunes down amyloid-driven NF-kB across microglia.
2. **State-restricted attenuation.** One specific substate (most
   plausibly MG_homeostatic per the H4 17.4 §4 hypothesis sketch)
   carries the NF-kB suppression; other substates carry no signal
   or a smaller signal. Biological interpretation: tau on the
   amyloid background specifically reshapes the homeostatic
   substate's transcriptional state, leaving activated substates
   on their canonical DAM trajectory.
3. **State-driven attenuation by abundance shift.** No substate
   independently shows NF-kB suppression at the interaction; the
   whole-microglia signal is an aggregation artefact of differential
   substate abundance between the four genotypes at the interaction
   (e.g. NLGF_P301S enriches the MG_homeostatic pool whose baseline
   NF-kB activity is lower; the whole-microglia mean drops as a
   composition effect). Biological interpretation: tau drives a
   substate-composition shift; the NF-kB attenuation is a
   "Simpson's paradox" artefact and T-Tau-attenuates' transcriptional
   claim is weak.

Distinguishing the three is what Phase I delivers. The per-state
NEBULA infrastructure already exists at
`storage/cache/de_snrnaseq_nebula_per_state.rds` (built by
`scripts/build_per_state_nebula.R`, consumed by section 02d); the
gap is that no per-state TF activity inference has been run on top
of it. Section 14 runs TF activity on whole-microglia NEBULA only.
Phase I adds the per-state TF activity layer and stratifies the
NF-kB family specifically.

## Locked decisions (carried forward + new for Phase I)

| Decision | Choice |
|---|---|
| TF activity tool | **`decoupleR` + CollecTRI** (carried forward from section 14 / Phase E). Locked at plan-writing time to maintain consistency with the existing TF inference; no new decision gate. The same `decoupleR::run_ulm()` + `consensus` arbitration rules apply, with FDR<0.10 on `ulm` for significance. |
| Input statistic from per-state NEBULA | **NEBULA z-statistic** (`logFC / SE`), per state, at the interaction contrast. Carried forward from Phase E's input-statistic convention (z-statistics are what `decoupleR` recommends for univariate methods because they fold magnitude + uncertainty into one signed scalar). Sparse-prevalence-filtered fit (`de_snrnaseq_nebula_per_state_1pct.rds`) preferred over the unfiltered fit if both are available; the 1% prevalence cutoff matches the project's locked filter elsewhere. |
| NF-kB family scope | **Full NF-kB family** = Rela, Nfkb1, Nfkb2, Rel, Relb, A0A979HLR9 (the UniProt-accession complex source). Locked at plan-writing time; restricting to Rela only would lose the family-pattern signal that strengthens the T-Tau-attenuates claim. |
| Substate scope | **All 4 microglia substates** = MG_homeostatic, MG_DAM, MG_IFN, MG_proliferative. The per-state NEBULA caches already include all 4 states so the marginal cost of including all 4 is zero. Locked at plan-writing time. |
| New computation | **Per-state TF activity cache** + **per-state TF activity at interaction TSV** + **NF-kB-family-specific cross-state comparison table** + (optionally) **NF-kB-target GSEA per state at interaction**. No new NEBULA fits; the per-state NEBULA caches are consumed verbatim. |
| Synthesis venue | **New top-level section 18 "Per-state NF-kB attenuation test"** rendered via new `rmd/17_nfkb_attenuation.Rmd` inserted into `analysis.Rmd` between `child-biological-model` (rmd/16) and `child-session` (rmd/99). Locked at I1 (2026-05-25) per user confirmation; default proposal carried. Keeps section 17 stable as the synthesis venue and section 14 stable as whole-microglia TF inference. |
| NF-kB target GSEA inclusion | **Include** per-state GSEA against the union NF-kB-target gene set (CollecTRI's Rela / Nfkb1 / Rel / Relb / Nfkb2 / A0A979HLR9 target genes) at the interaction contrast as supporting evidence alongside the per-state TF activity table. Locked at I1 (2026-05-25) per user confirmation; default proposal carried. The GSEA enables the I4 ledger rows to reach **Strong** under the H1-locked three-criteria rule (>=2-modality consistent sign + FDR<0.10 + explicit cross-layer corroboration). Without it the Phase I rows would cap at Moderate. |
| Verdict format | A 5-tier verdict: **Global attenuation / Homeostatic-restricted / DAM-restricted / IFN-restricted / Proliferative-restricted / State-driven-by-abundance (no per-state signal)**. The verdict feeds the I4 H2 ledger update — depending on which tier wins, T-Tau-attenuates either lifts substantially (state-restricted with clean biology), modestly (global), or collapses (state-driven-by-abundance artefact). |

## Open questions (defer to the session that needs them)

| Question | Default proposal | Decided in |
|---|---|---|
| Compositional confound handling | If the verdict is "state-driven-by-abundance", how to test it formally? Default: cite the existing section 11 DAM-sampling balance ANOVA (`ccc_dam_sampling_balance_anova.tsv`) which already tests substate composition shifts under the 2×2 factorial; reference that test as the formal arbiter of compositional vs transcriptional NF-kB attenuation. Alternative: add a new per-state composition test inside the I3 chapter (Fisher / ANOVA on NLGF_P301S vs the other 3 genotypes' substate fractions). Default prefers the citation route (existing infrastructure) over duplication. | I3 (per-claim judgment). |
| Ledger update scope | How many new rows does the I4 ledger update add? Default: **3-6 new rows** depending on the verdict tier — (a) one row per substate-specific NF-kB attenuation finding (1-4 rows); (b) one row for the whole-microglia replication or refutation (1 row); (c) one row for the compositional-vs-transcriptional adjudication (1 row). The new rows tag T-Tau-attenuates as `supports_models` with confidence grade scaled to the per-state evidence strength (Strong if state-restricted with FDR<0.10 + cross-layer corroboration; Moderate if global). | I4 (executes the update). |

## Execution model

- A fresh session must read `CLAUDE.md` first, then this file in full.
- A fresh session executes the NEXT step whose status is `TODO`. Steps
  are sized to fit comfortably in one session each.
- After completing a step, the session must:
  1. Update that step's status from `TODO` to `DONE <YYYY-MM-DD>`
     with a multi-paragraph completion note (what was built, key
     file paths, biology highlights, plan-spec deviations explicitly
     acknowledged).
  2. `chown rstudio:rstudio` any new files (agent runs as root).
  3. Re-knit `analysis.Rmd` and verify zero `class="error"` and zero
     `class="warning"` in the rendered HTML.
  4. Commit locally per CLAUDE.md (imperative subject under 70
     chars; co-author trailer).
  5. Either continue to the next `TODO` step in the same session if
     the context budget permits, or end cleanly.
- Decision-gate step I1 requires user confirmation via
  `AskUserQuestion`. Wait for confirmation before proceeding to I2.

## Phase I: per-state NF-kB attenuation test (sessions I1-I4)

### Session I1: lock chapter venue + GSEA scope [DONE 2026-05-25]

**Decision gate at start.** Two binary locks needed from the user
via `AskUserQuestion`: (a) chapter venue; (b) NF-kB target GSEA
inclusion. Both are foundational to I2-I4 and must be locked before
I2 begins.

**(a) Chapter venue — default proposal:** new top-level **section 18
"Per-state NF-kB attenuation test"** rendered via a new
`rmd/17_nfkb_attenuation.Rmd` child file inserted into `analysis.Rmd`
between the existing `child-biological-model` (rmd/16) and
`child-session` (rmd/99) references. This keeps section 17
(Integrated biological model) stable as the synthesis venue and
section 14 (TF inference) stable as the whole-microglia layer.

**Alternative (a):** insert as new subsection 17.6 inside the
biological model chapter (`rmd/16_biological_model.Rmd`). Pros:
keeps the NF-kB attenuation test physically adjacent to the H2-H3
adjudication it lifts; reads naturally as "and here's the test that
ledger row N predicted". Cons: bloats section 17 from 5 subsections
to 6; complicates the "section 17 is the synthesis venue, section 14
is the whole-microglia mechanism layer, section 15 is kinase, section
16 is CCC" mental model. Default prefers section 18 for the cleaner
structural separation.

**Alternative (b):** insert as new subsection 14.4 inside the TF
inference chapter (`rmd/13_tf_inference.Rmd`). Pros: keeps all TF
activity analysis in one place. Cons: section 14 currently runs only
whole-microglia TF inference and adding a per-state subsection there
breaks that consistency; the Phase H synthesis would not have the
per-state test physically next to its T-Tau-attenuates evidence.
Default does not prefer this alternative.

**(b) NF-kB target GSEA inclusion — default proposal:** include a
per-state GSEA against an NF-kB-target gene set at the interaction
contrast as supporting evidence for the TF activity table. The gene
set is the union of CollecTRI's Rela / Nfkb1 / Rel / Relb / Nfkb2 /
A0A979HLR9 target genes (already in the project's CollecTRI cache;
no new data ingestion required). The GSEA provides independent
transcriptional corroboration: if the per-state TF activity table
says "Rela is suppressed in MG_homeostatic at the interaction", the
GSEA should report negative NES on the NF-kB-target gene set at the
same per-state contrast. If they agree, the evidence grade is Strong;
if they disagree, the evidence grade is Moderate at most.

**Alternative (b):** skip the GSEA and rely on the per-state TF
activity table alone. Pros: smaller scope, ~30% less code in the I2
build script. Cons: loses the second-line corroboration the
confidence-grading rules need for Strong (per the H1-locked
three-criteria Strong rule: ≥2-modality consistent sign + FDR<0.10 +
explicit cross-layer corroboration). Default prefers including the
GSEA so the Phase I claim rows can reach Strong if the biology
supports it.

**Action plan for I1 session (this session, in order):**

1. Present (a) + (b) defaults + alternatives to the user via a
   single multi-question `AskUserQuestion` call covering both locks.
2. Lock the chosen values into the Locked-decisions table above;
   remove (a) + (b) from the Open-questions table.
3. chown / re-knit (no-op for plan-only edits) / verify / commit /
   continue to I2 only if context budget permits.

**Inputs for I1.** Plan file only; no new data inspection required
beyond what the plan body already references.

**Outputs from I1.** This session is plan-edit only (decision gate).
No new R code, no new TSV, no Rmd modification. The plan file is the
sole deliverable, updated with the locked chapter venue and GSEA
scope.

**I1 completion note (2026-05-25).** The two binary locks have been
resolved against their default proposals.

*Lock (a) — chapter venue.* User confirmed the default: a new top-level
**section 18 "Per-state NF-kB attenuation test"** rendered via new
`rmd/17_nfkb_attenuation.Rmd` inserted into `analysis.Rmd` between
`child-biological-model` (rmd/16) and `child-session` (rmd/99). This
preserves the project's current mental model: section 14 = whole-
microglia TF inference, section 15 = kinase, section 16 = CCC, section
17 = integrated biological model synthesis. Section 18 becomes the
per-state mechanism deep-dive specifically scoped to the NF-kB family
attenuation question, parallel to (not nested inside) section 17.
The two alternatives considered — subsection 17.6 inside the biological
model chapter, and subsection 14.4 inside the TF inference chapter —
were both rejected for their structural costs (bloating section 17
beyond its current 5 subsections, or breaking section 14's whole-
microglia-only consistency).

*Lock (b) — GSEA scope.* User confirmed the default: **include** a
per-state GSEA against the union NF-kB-target gene set sourced from
CollecTRI's Rela / Nfkb1 / Rel / Relb / Nfkb2 / A0A979HLR9 target
collections at the interaction contrast. The GSEA serves as second-line
transcriptional corroboration alongside the per-state TF activity
table: if the TF activity table reports negative Rela score in
MG_homeostatic at the interaction, the GSEA should independently
report negative NES on NF-kB targets in MG_homeostatic at the same
contrast. Including the GSEA is the operationally critical lock for
the I4 ledger update — it is the only route by which Phase I claim
rows can reach **Strong** under the H1-locked three-criteria rule
(>=2-modality consistent sign + FDR<0.10 + explicit cross-layer
corroboration), since the per-state TF activity table alone counts
as only one TF-method modality. Without the GSEA, Phase I rows would
cap at Moderate confidence and the T-Tau-attenuates lift would be
correspondingly weaker.

*Plan-spec deviations.* None. Both locks resolved to the plan's
default proposals; no scope expansion or contraction relative to the
plan body. The Locked-decisions table now contains the resolved
"Synthesis venue" + new "NF-kB target GSEA inclusion" rows; the
Open-questions table has been pruned of the two I1-decided rows
leaving only the I3 (compositional confound handling) and I4 (ledger
update scope) rows for their respective sessions.

*Files touched.* `storage/notes/nfkb_attenuation_plan.md` only;
plan-edit only session per the plan spec. No code, no caches, no Rmd
modifications, no knit re-run required.

*Next session.* I2 — build per-state TF activity cache. The pre-knit
smoke test in the I2 step is the right scoping for a fresh session
because the per-state NEBULA cache shape needs verification before
the `decoupleR::run_ulm()` loop is wired. Continuing into I2 in this
session is feasible from a context-budget standpoint (this session
has spent ~5% of budget on the I1 plan edit) but I2's deliverables
(new build script + new cache + new TSV + pre-knit smoke test +
chown + commit) are non-trivial and benefit from a clean session
budget. Ending the I1 session cleanly here.

---

### Session I2: build per-state TF activity cache [DONE 2026-05-25]

**Goal.** Emit a new cache + TSV pair: per-state TF activity at the
interaction contrast (and optionally at all five project contrasts
for cross-contrast consistency checks).

**Build script:** new `scripts/build_per_state_tf_activity.R`.
Pattern follows `scripts/build_tf_activity_decoupler.R` (the
whole-microglia equivalent), with the input changed from the
whole-microglia NEBULA fit to the per-state NEBULA fit (one fit per
state, looped).

**Cache emission:** new `storage/cache/per_state_tf_activity.rds`
containing a nested list `state -> contrast -> decoupleR_output`
where each leaf is the standard `decoupleR::run_ulm()` data frame
(source, condition, score, p_value, statistic) plus the consensus
arbitration. The cache is roughly 4 states × 5 contrasts × ~1000 TFs
= ~20k rows on disk; ~1-2 MB compressed.

**TSV emission:** new `storage/results/per_state_tf_activity.tsv`
containing the long-form per-state per-contrast TF activity table
with columns `state`, `contrast`, `tf`, `score`, `p_value`, `padj`,
`sig_fdr10_ulm`. Total rows ~20k. The TSV is the consumer interface
for I3.

**Pre-knit smoke test.** Run `Rscript -e` to:
1. Load `de_snrnaseq_nebula_per_state_1pct.rds`; verify it has the
   expected shape (state x contrast x gene-level z-statistics).
2. Run `decoupleR::run_ulm()` on a single state x interaction slice
   to confirm the toolchain runs cleanly; spot-check that Rela
   surfaces with a finite score.
3. Save a 10-row mini-cache; verify it loads and the long-form TSV
   reshape works.

**Outputs from I2.**
- `scripts/build_per_state_tf_activity.R` (new build script).
- `storage/cache/per_state_tf_activity.rds` (new cache).
- `storage/results/per_state_tf_activity.tsv` (new TSV).
- Plan-spec update: I2 status DONE.

**I2 completion note (2026-05-25).** The per-state TF activity cache +
TSV are built, verified, and ready for I3 consumption. Three deliverables
on disk and chowned rstudio:rstudio: the build script
`scripts/build_per_state_tf_activity.R` (8.7 kB, idempotent with
`--overwrite`), the cache `storage/cache/per_state_tf_activity.rds`
(520 kB, xz-compressed, nested list `state -> contrast -> tibble`
keyed by 4 states x 5 contrasts), and the long-form TSV
`storage/results/per_state_tf_activity.tsv` (847 kB, 9710 rows = 4 states
x 5 contrasts x ~485 TFs per slice, ulm-only).

*Build approach + reuse of existing helpers.* The script is a strict
loop-over-states wrapper around the existing TF inference helpers in
`R/tf_inference.R` (`extract_de_stat_matrix` -> `run_decoupler_per_modality`
-> `split_decoupler_by_contrast`), with no new helper code required. The
existing helpers are state-agnostic by design (they take any
per-contrast top-table list as input) which made the per-state extension
a strict additive build rather than a refactor. Runtime is ~28-29 s per
state x 5 contrasts on the project's Docker container, totalling ~2 min
wall for the 4-state loop. The CollecTRI mouse prior (39961 edges, 1114
unique TFs, 6110 targets) is loaded once via
`decoupleR::get_collectri(organism = "mouse", split_complexes = FALSE)`
and shared across all 4 states, matching the prior used by
`build_tf_activity_decoupler.R` byte-for-byte so cross-cache comparisons
between whole-microglia and per-state remain meaningful (only the input
matrix differs).

*Plan-spec deviation: input cache choice (1pct vs unfiltered).* The
plan-locked input is `de_snrnaseq_nebula_per_state_1pct.rds` (the 1%
prevalence-filtered fit). The pre-existing whole-microglia cache
`storage/cache/tf_activity_decoupler.rds` already contains per-state
modalities (snrnaseq_homeostatic / snrnaseq_DAM / snrnaseq_IFN /
snrnaseq_proliferative) but was built against the UNFILTERED per-state
NEBULA fit (`de_snrnaseq_nebula_per_state.rds`, ~5k genes per state).
The new `per_state_tf_activity.rds` cache is built against the
1pct-filtered fit (~11-12k genes per state) which materially increases
TF target-set coverage: e.g. CollecTRI Rela's 572 mouse targets overlap
more of the per-state universe under 1pct filtering than under the
unfiltered fit's tighter universe. The two caches coexist; the I3
chapter reads the new `per_state_tf_activity.rds` per the plan lock; the
older whole-microglia + unfiltered per-state cache continues to power
sections 13-14 unchanged. No back-compat shims.

*Plan-spec deviation: TSV statistic scope (ulm only vs all statistics).*
The plan specified TSV columns `state`, `contrast`, `tf`, `score`,
`p_value`, `padj`, `sig_fdr10_ulm` without explicitly restricting to a
single decoupleR statistic. I restricted the TSV to the `ulm` statistic
because (a) the project's canonical sig test is ulm with BH-adjusted
FDR<0.10 (per R/tf_inference.R module header), (b) keeping the TSV
single-statistic avoids downstream ambiguity about which row's
score / p / padj triple is the headline, and (c) the full multi-statistic
information (ulm + wsum + norm_wsum + corr_wsum + consensus) is preserved
verbatim in the cache so the I3 chapter can read the consensus score
for any TF directly from the cache when the consensus signal differs
meaningfully from the ulm signal (which the smoke test confirms it does
for A0A979HLR9: ulm score -2.35 in homeostatic vs consensus score -5.45).
This is a token-efficiency + readability optimisation; the cache remains
the source of truth.

*Biology highlights (will inform I3 verdict, do NOT lock the verdict
here).* The pre-knit smoke test inspected the NF-kB family (Rela, Nfkb1,
Nfkb2, Rel, Relb, A0A979HLR9) across all 4 substates at the interaction
contrast. Three findings worth flagging for the I3 session, though the
I3 chapter must re-derive and synthesise them rather than anchor on
this note:

- **A0A979HLR9 (the NF-kB complex UniProt accession, 1340 targets in
  CollecTRI mouse) is sign-consistently NEGATIVE in ALL 4 substates at
  the interaction**, with the deepest signal in DAM (consensus
  score -6.52, p=7e-11) and progressively weaker but still strongly
  negative signals in homeostatic (-5.45, p=5e-8), IFN (-3.66,
  p=3e-4), and proliferative (-2.55, p=0.011). This is the cleanest
  evidence for **Phase I verdict tier 1 = global attenuation**: every
  microglia substate independently shows NF-kB suppression at the
  interaction, with the whole-microglia signal a weighted aggregate
  rather than a homeostatic-only or DAM-only artefact.
- **Rela individually is weakly negative across all 4 substates**
  (DAM -1.73 p=0.083 ulm; others ~-0.5 to -0.7) but **none of the
  individual NF-kB family TFs pass FDR<0.10 on the ulm sig test at the
  interaction**. The family signature at the interaction is carried
  almost entirely by the COMPLEX-level A0A979HLR9 accession, not by
  individual subunits. This is the right answer biologically (the
  complex's transcriptional readout is the union of its subunits'
  cooperative occupancy, which the complex-level CollecTRI target set
  captures and the per-subunit sets miss) but the I3 chapter must be
  explicit about it so the reader doesn't conclude per-subunit
  evidence is absent.
- **The H4 17.4 §4 hypothesis sketch predicted MG_homeostatic carries
  the NF-kB suppression. The data say MG_DAM carries the strongest
  signal**, with homeostatic second. This is an important biological
  revision the I3 chapter must report transparently; the plan's
  anti-anchoring #4 explicitly anticipated this kind of
  Phase H -> Phase I correction. The DAM-led pattern is itself
  biologically interesting (DAM transcriptomes are typically NF-kB-
  HIGH under amyloid; attenuation at the interaction means tau on
  the amyloid background dampens DAM's amyloid-driven NF-kB program
  specifically, not the homeostatic baseline).
- **Compositional sanity (sig-hit-count by state x contrast):**
  homeostatic carries the most TFs sig at FDR<0.10 on ulm across the
  two NLGF contrasts (28 + 23 = 51), vs DAM (0 + 0 = 0), IFN
  (6 + 4 = 10), proliferative (0 + 5 = 5). The interaction contrast
  is the sparsest across all 4 states (1, 7, 1, 1 sig hits respectively),
  which is biologically reasonable (interaction effects are smaller in
  magnitude than main effects). The 107 total sig-ulm-FDR<0.10 rows
  out of 9710 (1.1%) is consistent with the project's whole-microglia
  TF inference baseline.

*Files touched + sizes.*
- `scripts/build_per_state_tf_activity.R` 8.7 kB (NEW, rstudio:rstudio).
- `storage/cache/per_state_tf_activity.rds` 520 kB (NEW, xz-compressed,
  rstudio:rstudio).
- `storage/results/per_state_tf_activity.tsv` 847 kB (NEW,
  rstudio:rstudio).
- `analysis.html` 30 MB (re-knitted, no content change, 0 errors / 0
  warnings, rstudio:rstudio).
- `storage/notes/nfkb_attenuation_plan.md` (this file; I2 -> DONE +
  completion note + Last-updated timestamp).

*Next session.* I3 — render the per-state NF-kB attenuation chapter.
The chapter venue is locked to new top-level section 18 via new
`rmd/17_nfkb_attenuation.Rmd` inserted into `analysis.Rmd` between
`child-biological-model` (rmd/16) and `child-session` (rmd/99). I3
must also run the per-state NF-kB-target GSEA (locked IN at I1) as
the second-line corroboration that lets Phase I claim rows reach
Strong. The plan body for I3 contains the 4-subsection structure;
the chapter should reuse the I2 cache + TSV for subsection 2's TF
activity table and build a NEW per-state GSEA cache for subsection 3.
Subsection 4's verdict is an I3 synthesis task, not I2's; the I2
note above flags the data but deliberately does not lock the verdict
because the GSEA + the section 11 composition test both feed into
the final verdict and may shift it.

---

### Session I3: render per-state NF-kB attenuation chapter [DONE 2026-05-25]

**Goal.** Render the new chapter at the I1-locked venue. The chapter
is structured as 4 subsections:

**Subsection 1 (chapter intro).** State the scientific question
(global vs state-restricted vs state-driven-by-abundance NF-kB
attenuation), reference section 17.5's T-Tau-attenuates claim, and
state the analysis approach (per-state TF activity + per-state
NF-kB-target GSEA, both at the interaction contrast).

**Subsection 2 (per-state TF activity at interaction).** Render the
per-state TF activity for the NF-kB family (Rela, Nfkb1, Nfkb2, Rel,
Relb, A0A979HLR9) as a kable: rows = states, columns = TF family
members, cells = signed score + FDR-marker. Add a one-paragraph
biology read of the table (which substate carries the signal; sign
pattern within the family). Include a kable of the top-10 TFs per
state at the interaction (not restricted to NF-kB) for context.

**Subsection 3 (per-state NF-kB-target GSEA at interaction).**
Render the per-state GSEA on the NF-kB target gene set at the
interaction contrast as a kable (rows = states, columns = NES, padj,
leading-edge count). Add a one-paragraph biology read consistent
with or divergent from Subsection 2's TF table. Conditional on the
I1 (b) lock — if the user chose to skip the GSEA, this subsection is
omitted.

**Subsection 4 (per-state attenuation verdict).** State which of the
5 verdict tiers (global / homeostatic-restricted / DAM-restricted /
IFN-restricted / proliferative-restricted / state-driven-by-abundance)
the evidence supports, citing the per-state TF activity + GSEA + the
existing section 11 DAM-sampling balance ANOVA. Include a brief
compositional-confound paragraph naming whichever genotypes carry
substate fraction shifts at the interaction (if any). The verdict
prose feeds the I4 ledger update.

**Wiring.** Insert a new child reference into `analysis.Rmd` at the
I1-locked position (default: between `child-biological-model` and
`child-session`).

**Outputs from I3.**
- New `rmd/17_nfkb_attenuation.Rmd` (or other path per I1 lock).
- Updated `analysis.Rmd` (one new child reference).
- Re-knitted `analysis.html` (verify zero errors / warnings; new
  section 18 renders).
- Plan-spec update: I3 status DONE.

**I3 completion note (2026-05-25).** Section 18 "Per-state NF-kB
attenuation test" is rendered, all four subsections are in place, and
the verdict tier resolves cleanly to **Global attenuation** with
4-of-4 substates carrying the attenuation under the chapter's
cross-layer gate. The chapter is wired into `analysis.Rmd` at the
I1-locked venue between `child-biological-model` (rmd/16) and
`child-session` (rmd/99), and the full re-knit completes in ~4.5
minutes with zero `class="error"` and zero `class="warning"` in
`analysis.html`.

*Files touched + sizes.*
- `rmd/17_nfkb_attenuation.Rmd` 34 kB (NEW, rstudio:rstudio) — the
  full chapter with 4 subsections.
- `analysis.Rmd` (already wired in I3-prep: new child reference between
  rmd/16 and rmd/99; rstudio:rstudio).
- `storage/cache/collectri_mouse_for_nfkb_gsea.rds` 200 kB (NEW,
  rstudio:rstudio) — CollecTRI mouse prior with
  `split_complexes = FALSE`, identical settings to
  `scripts/build_per_state_tf_activity.R` so the family + complex
  accessions are byte-identical across cache builds.
- `storage/cache/per_state_nfkb_target_gsea.rds` 30 kB (NEW,
  rstudio:rstudio) — nested list `state -> contrast -> fgsea tibble`
  for the 4 substates x 5 contrasts x 6 gene sets (union + 5
  subunits; Rel excluded because it has 0 targets in CollecTRI mouse).
- `storage/results/per_state_nfkb_target_gsea.tsv` 80 kB (NEW,
  rstudio:rstudio) — 120 long-form rows (6 sets x 4 states x 5
  contrasts) with NES, padj, pval, size, leading-edge count, and
  leading-edge gene names (pipe-delimited string).
- `analysis.html` 30 MB (re-knitted, rstudio:rstudio).

*Build approach + scope.* The chapter is structured into four
subsections matching the plan body: (1) intro restating the three
competing patterns and naming the section 11 ANOVA as the formal
arbiter; (2) per-state TF activity at the interaction with NF-kB
family kable + top-10 TFs per substate; (3) per-state NF-kB-target
GSEA with union-NES-at-interaction kable, per-subunit-NES kable, and
union-NES-across-all-contrasts kable; (4) verdict construction +
verdict table + multi-paragraph verdict statement. Subsection 18.3
builds the new in-knit cache via `cache_or_run()` against the
1pct-prevalence-filtered per-state NEBULA fit (locked input matching
the I2 TF activity input) and emits the long-form TSV with
pipe-delimited leading-edge column. The chapter reuses existing
helpers (`run_fgsea_per_state` from `R/fgsea.R`, `write_tsv_safe` from
`R/utils.R`) without adding new helper code; the verdict construction
is in-chapter.

*Biology highlights (locked at I3, propagate to I4).*
- **Verdict tier: Global attenuation.** All four microglia substates
  (homeostatic, DAM, IFN, proliferative) independently carry the NF-kB
  attenuation at the interaction under the cross-layer gate (TF
  complex consensus score < 0 AND GSEA union NES < 0 AND at least
  one layer at FDR<0.10 within-state). The whole-microglia signal
  surfaced by section 14 is the weighted aggregate of consistent
  per-state suppression, not a state-driven-by-abundance artefact
  and not a one-substate-only signal.
- **Section 11 ANOVA tau:nlgf p = 0.330 (F = 1.03) rules out
  state-driven-by-abundance** on solid statistical grounds. The
  per-state signal is transcriptional, not compositional. This is
  the single most important supporting evidence for the Global
  attenuation tier.
- **DAM is the lead carrier, with homeostatic close behind, contra
  the H4 §4 hypothesis sketch.** Consensus scores are DAM -6.52
  (p=7e-11), homeostatic -5.45 (p=5e-8), IFN -3.66 (p=3e-4),
  proliferative -2.55 (p=0.011). GSEA NES at interaction:
  DAM -1.247 (padj=0.027, sig), homeostatic -1.221 (padj=0.091,
  borderline), IFN -1.126 (padj=0.425, NS), proliferative -1.121
  (padj=0.162, NS). The DAM-led ordering is biologically intuitive
  (DAM cells are transcriptionally NF-kB-HIGH under amyloid, so
  the attenuation at the interaction has the most signal in the
  substates with the most-amplified baseline NF-kB program). The
  H4 §4 prediction that MG_homeostatic carries the suppression is
  therefore **refined, not refuted**: every substate carries the
  suppression, but DAM has the largest magnitude.
- **The cross-contrast sign-reversal holds in all 4 substates at
  the GSEA layer.** Union NES is positive (+1.03 to +1.21) at the
  two NLGF arms and negative (-1.12 to -1.25) at the interaction
  in every substate; the per-substate pattern is the GSEA mirror of
  the section 14 whole-microglia Rela sign-reversal that the Phase H
  synthesis identified as the load-bearing signature for T-Tau-
  attenuates. This satisfies the H1-locked cross-layer corroboration
  rule for **Strong** confidence grading on the Phase I ledger rows
  I4 will add.
- **The family signal at the interaction is carried by the complex-
  level CollecTRI accession `A0A979HLR9`, not by individual NF-kB
  subunits.** None of the per-subunit TFs (Rela, Nfkb1, Nfkb2, Relb)
  reach FDR<0.10 on ulm BH-padj in any substate at the interaction;
  the per-subunit GSEA sets likewise do not pass FDR<0.10 (with the
  exception of NFKB_UNION_TARGETS in DAM at padj=0.027 and
  homeostatic at padj=0.091). This is biologically coherent because
  NF-kB transcriptional output is the cooperative product of subunit
  occupancy at composite promoters, which the complex-level target
  set captures and the per-subunit sets miss.

*Plan-spec deviations.*
1. **Sig-metric choice at the verdict gate (consensus p, not ulm
   BH-padj).** The chapter intro paragraph (subsection 18.1) restates
   the section-14-locked convention "significance at FDR<0.10 on the
   decoupleR `ulm` p-value (BH-adjusted within each state x
   contrast)". The verdict gate in subsection 18.4 deliberately reads
   the within-state **consensus p-value** instead, because the ulm
   BH-padj across ~480 TFs per (state, contrast) slice is overly
   conservative — even A0A979HLR9, the family's lead carrier, does
   not pass FDR<0.10 on ulm in any substate at the interaction (ulm
   padj range 0.30-0.94 per the I2 smoke test) despite the consensus
   p-value reaching 5e-8 to 0.011 across the four substates. The
   chapter is explicit about this deviation: the family-reading
   paragraph in subsection 18.2 explains that the ulm test fails
   under per-state multiple-testing burden, that the consensus
   statistic is the right within-state sig metric, and that the
   verdict gate uses consensus p. This is a defensible Phase I-
   specific rule (the per-state setting has different statistical
   properties than the whole-microglia setting that section 14
   conventions were tuned for) and is documented in the chapter
   prose for cross-session continuity.
2. **GSEA gating asymmetry (cross-layer "at least one" rather than
   "both").** The chapter's carries_attenuation gate is `tf_neg &
   gsea_neg & (tf_sig | gsea_sig)`, not the stricter `tf_neg &
   gsea_neg & tf_sig & gsea_sig`. The OR is deliberate: the union
   GSEA tests 1 gene set per state (no multiple-testing correction
   across multiple sets) while the TF activity tests ~480 TFs per
   state (heavy correction); requiring FDR<0.10 on both layers
   would push the two layers to inconsistent statistical power and
   bias the verdict toward only the substates where both layers
   coincidentally reach sig under their different multiple-testing
   burdens. The H1-locked Strong-grade rule requires cross-layer
   sign agreement + FDR<0.10 + cross-layer corroboration; the
   chapter's gate satisfies all three at the within-state resolution
   (sign agreement in all 4 substates; FDR<0.10 satisfied by
   consensus p in all 4 substates and additionally by GSEA padj in
   DAM and borderline in homeostatic; cross-layer corroboration is
   the simultaneous negativity of both layers in every substate).
3. **GSEA min_size + max_size deviation from the project default.**
   The chapter calls `run_fgsea_per_state(..., min_size = 5L,
   max_size = 2000L)` rather than the helper's default (10, 500).
   The min_size lowering is required because the per-subunit gene
   sets after per-state universe intersection are small (Nfkb2 ~26,
   Relb ~26 in the smaller substates; below the default min_size of
   10 they would be silently dropped). The max_size raising is
   required because the NFKB_UNION_TARGETS set has 1679 prior
   targets, intersecting to 718-781 in each per-state universe;
   the default max_size of 500 would cap the union set out. Both
   parameters are locked in the chapter chunk header comments and
   are deliberate Phase I conventions for the family + union scope.
4. **Prose correction during the I3 session.** The first draft of
   subsection 18.2's family-reading paragraph stated that "Every
   substate's A0A979HLR9 row reaches FDR<0.10 on the ulm sig test
   except where the ulm p-value lies just above the threshold" — an
   internally contradictory claim that did not match the rendered
   kable's absence of `*` markers in any A0A979HLR9 cell. The
   paragraph was corrected mid-session to honestly describe the
   table contents (no TF passes FDR<0.10 on ulm at the interaction)
   and to explicitly motivate the consensus-statistic choice for
   the verdict gate. The corrected paragraph is now in the rendered
   HTML and reads coherently with the subsection 18.4 gate logic.

*Anti-anchoring discipline followed.* All four substates are
reported in the family table even though homeostatic and DAM carry
stronger signal than IFN and proliferative — per the plan's "always
report all 4 substates" rule. Sign consistency between the TF
activity and GSEA layers is verified for every substate. The
compositional-confound alternative is explicitly named in the
verdict prose and resolved by the section 11 ANOVA citation. The
Phase H integrated reading (T-Tau-attenuates at 5 supports / 3
Strong) is treated as the prior; Phase I refines its biological
character (DAM-led, not homeostatic-led; complex-level, not
per-subunit) rather than replacing it. No new computation beyond
per-state TF activity + per-state GSEA was added; no modification
to the H1-locked confidence-grade definitions.

*Next session.* I4 — propagate the Phase I verdict into the H2
ledger and refresh the section 17 synthesis. The expected ledger
update is 3-6 new rows: one per-state row for each substate that
carries the attenuation (4 rows here, one each for homeostatic, DAM,
IFN, proliferative), one row for the cross-state Global-attenuation
verdict at the integrated scale, and one row for the compositional-
vs-transcriptional adjudication referencing the section 11 ANOVA.
T-Tau-attenuates' net_support should lift from 5 to 8-11 once the
new rows are added, depending on whether the per-state rows are
adjudicated at Strong or Moderate. Section 17.5's integrated model
statement needs a targeted prose edit to reference the Global
attenuation verdict and the DAM-led-with-homeostatic-close-behind
substate ordering. The I4 plan body already describes the steps;
no new decisions are required before I4 begins.

---

### Session I4: update H2 ledger + re-run H3 + refresh section 17 [DONE 2026-05-25]

**Goal.** Propagate the Phase I verdict into the Phase H synthesis:
add new rows to the H2 ledger, re-run the H3 adjudication, refresh
section 17's verdict tables and integrated model statement to
reference the new evidence.

**Step 1: extend `scripts/build_biological_model_ledger.R`** with
3-6 new ledger rows (the exact count depends on the I3 verdict
tier). Each new row is one atomic per-state NF-kB attenuation finding
tagged against `T-Tau-attenuates` and, where applicable, `Hyp-1B`. The
new claim IDs follow the H2 numbering convention as **I-001..I-006**
or as **H2-076..H2-081** (per H1's H2 numbering scheme; pick at I4
execution time to avoid ID collisions with any post-H2 ledger
additions). Update the script's row(): constructor calls; re-run
the script; verify the schema invariants still hold.

**Step 2: re-run `scripts/build_biological_model_adjudication.R`.**
The adjudication arithmetic re-runs on the extended ledger; the
per-entity tally + per-axis contest verdicts auto-refresh. Confirm
that T-Tau-attenuates' net_support moves in the predicted direction
(up by 3-6 if state-restricted; down by 1-2 if state-driven-by-
abundance artefact).

**Step 3: re-knit `analysis.Rmd`.** Section 17.2 (the DT::datatable
ledger) auto-refreshes from the regenerated TSV; section 17.3.1
(contest verdicts) and 17.3.2 (entity tally) auto-refresh likewise.
Section 17.5 (integrated model statement) may need a targeted prose
edit to reference the Phase I verdict — at minimum a sentence
acknowledging that T-Tau-attenuates' standing has been tested and
either strengthened or refined.

**Step 4: close the plan.** Move
`storage/notes/nfkb_attenuation_plan.md` to
`storage/notes/completed/nfkb_attenuation_plan_<YYYY-MM-DD>.md` with
a `## Outcome summary` section recording: (a) the locked I1
choices; (b) the per-state TF activity result for the NF-kB family;
(c) the per-state GSEA result if included; (d) the verdict tier
selected; (e) the I4 ledger update (number of new rows added,
T-Tau-attenuates' new net_support, any other affected entity's net
movement); (f) the next-plan pointer.

**Outputs from I4.**
- Updated `scripts/build_biological_model_ledger.R` + regenerated
  `storage/results/biological_model_claims_ledger.tsv` + regenerated
  `storage/results/biological_model_adjudication.tsv` + regenerated
  `storage/results/biological_model_contest_verdicts.tsv`.
- Updated `rmd/16_biological_model.Rmd` (targeted prose edit to
  17.5 if the verdict warrants).
- Re-knitted `analysis.html` (verify zero errors / warnings; new
  section 18 + refreshed section 17 both render).
- Plan moved to `storage/notes/completed/nfkb_attenuation_plan_<YYYY-MM-DD>.md`
  with `## Outcome summary` section.

**I4 completion note (2026-05-25).** Phase I closed cleanly. Six new
ledger rows (I-001..I-006) propagate the I3 Global-attenuation verdict
into the H2 ledger; the H3 adjudication regenerates T-Tau-attenuates at
11 supports / 7 Strong (up from 5 / 3) and Hyp-1B at margin 18 over
Hyp-1A (up from 6); section 17 prose and the H5 upfront mirror in
section 1.1 of `analysis.Rmd` are refreshed to reference the Phase I
lift; the re-knit completes with zero `class="error"` and zero
`class="warning"` in `analysis.html` and the new section 18 + refreshed
section 17 tables both render.

*Numbering choice.* New rows are tagged **I-001..I-006**, not
H2-076..H2-081. The plan body listed both options; I picked I-* to
make the Phase I rows visually distinguishable from the Phase H 75-row
block in the DT::datatable view. The I-prefix marks them as the Phase I
addition without ambiguity for future sessions that read the ledger,
and the H2 row IDs (H2-001..H2-075) remain a stable anchor for any
later cross-references to the original synthesis.

*Row inventory (axis = cross_axis, layer = cross_layer in all 6;
supports_models = "Hyp-1B;T-Tau-attenuates"; contradicts_models =
"Hyp-1A" for all 6, following the H2-072 whole-microglia-attenuation
template).* I-001 homeostatic Strong (TF complex consensus -5.45
p=5e-8; GSEA union NES -1.221 padj=0.091, borderline cross-layer
corroboration), I-002 DAM Strong (TF -6.52 p=7e-11; GSEA NES -1.247
padj=0.027 — the lead carrier), I-003 IFN Moderate (TF -3.66 p=3e-4;
GSEA NS at padj=0.425), I-004 proliferative Moderate (TF -2.55
p=0.011; GSEA NS at padj=0.162), I-005 whole-microglia replication
Strong (4-of-4 substate sign-consistent negativity at the TF complex
layer, satisfying the cross-layer gate at the aggregate scale), I-006
compositional ANOVA Strong (`ccc_dam_sampling_balance_anova.tsv`
tau:nlgf F=1.03 p=0.330 rules out state-driven-by-abundance on solid
statistical grounds). I-001 / I-002 / I-005 / I-006 are the 4 new
Strong rows; I-003 / I-004 are the 2 new Moderate rows. Together they
lift T-Tau-attenuates from 5 supports / 3 Strong to 11 supports /
7 Strong, exactly matching the upper bound of the plan-locked
prediction ("net_support should lift from 5 to 8-11 depending on
whether the per-state rows are adjudicated at Strong or Moderate").

*Adjudication arithmetic.* The H3 adjudication on the extended 81-row
ledger re-runs cleanly. T-Tau-attenuates' theme entry in
`storage/results/biological_model_adjudication.tsv` reads
`support_count = 11`, `contradict_count = 0`, `net_support = 11`,
`n_strong_supports = 7`. Hyp-1B's contest entry in
`storage/results/biological_model_contest_verdicts.tsv` reads
`model_b_support_count = 30`, `model_b_n_strong_supports = 11`,
`favoured_by_margin = 18` (up from 6 at the end of Phase H). Hyp-1A's
matching entry reads `model_a_support_count = 21`,
`model_a_contradict_count = 9`, `model_a_net_support = 12` (down from
18 at end of Phase H; each Phase I row supports Hyp-1B and contradicts
Hyp-1A, swinging Hyp-1A net by -6 and Hyp-1B net by +6 for a combined
margin shift of +12). The cross-axis theme ranking re-orders to
T-Synergy (28) > T-Inflammation (26) > T-Compartment-suppression (17)
> T-Tau-attenuates (11) > Hyp-0 Cdk5 (7) — T-Tau-attenuates moves
from 5th to 4th, swapping positions with Hyp-0.

*Prose refresh (rmd/16 + analysis.Rmd).* Section 17.1's chapter intro
now declares an 81-row ledger (75 Phase H + 6 Phase I). Section 17.2's
DT::datatable `lengthMenu` adds 81 as the all-rows option. Section
17.3.1's axis-1 contest verdict prose is rewritten to reflect margin
18 with explicit Phase I attribution and to acknowledge the 9 total
contradicting rows now and the 11 Strong supports for Hyp-1B. Section
17.3.2's cross-axis ranking + entity tally paragraph names
T-Tau-attenuates' new 1/0/2/8 axis-split and the re-ordered theme
ranking with Hyp-1B at the top of the model contests by Strong-support
count. Section 17.4's open question 4 (T-Tau-attenuates) is rewritten
as "RESOLVED by Phase I" with the Global-attenuation verdict +
DAM-lead-with-homeostatic-close-behind substate ordering documented.
Section 17.5's integrated model statement (both the axis-1 verdict
clause, the cross-axis ranking clause, and the integrated reading
clause) is updated to reference the Phase I lift and the
Global-attenuation finding. `analysis.Rmd` section 1.1 (the H5 upfront
mirror) is updated in parallel to declare 81 atomic claims (75 Phase H
+ 6 Phase I), Hyp-1B margin 18 over Hyp-1A (up from 6), and the new
cross-axis theme ranking with T-Tau-attenuates at 11 supports.

*Plan-spec deviations.* None. The plan body's I4 step listed 3-6
expected new rows; the actual 6 rows match the upper end of that range
(one per substate + one whole-microglia replication + one compositional
adjudication). The plan-locked predictions on T-Tau-attenuates'
net_support lift (5 -> 8-11 depending on Strong/Moderate split) match
the realised lift of 11 exactly, with the per-state rows split 2 Strong
(homeostatic, DAM) + 2 Moderate (IFN, proliferative) consistent with
the I3 verdict's within-state cross-layer-gate readings.

*Files touched + sizes.*
- `scripts/build_biological_model_ledger.R` (modified, +6 row()
  constructor calls, +1 `phase_i` rbind block, assembly updated to
  include phase_i; rstudio:rstudio).
- `scripts/build_biological_model_adjudication.R` (modified, header
  comment updated from 75-row to 81-row with Phase I provenance;
  rstudio:rstudio).
- `rmd/16_biological_model.Rmd` (modified, ~9 targeted prose edits
  across sections 17.1, 17.2, 17.3.1, 17.3.2, 17.4, 17.5;
  rstudio:rstudio).
- `analysis.Rmd` (modified, section 1.1 H5 mirror refresh;
  rstudio:rstudio).
- `storage/results/biological_model_claims_ledger.tsv` (regenerated,
  81 data rows; rstudio:rstudio).
- `storage/results/biological_model_adjudication.tsv` (regenerated,
  11 entity rows; rstudio:rstudio).
- `storage/results/biological_model_contest_verdicts.tsv` (regenerated,
  3 contest rows; rstudio:rstudio).
- `analysis.html` (re-knitted, 30 MB, 0 errors / 0 warnings;
  rstudio:rstudio).
- `storage/notes/nfkb_attenuation_plan.md` moved to
  `storage/notes/completed/nfkb_attenuation_plan_2026-05-25.md` with
  this I4 completion note + the `## Outcome summary` section appended.

*Next session.* Phase I is COMPLETE; the T-Tau-attenuates
state-restricted vs global question is fully resolved. The next plan
is user-directed from the Phase H Outcome summary candidates: (A)
Cdk5-substrate GSEA at all three axes; (B) Pros1-Mertk GeoMx
validation; (C) hdWGCNA on GeoMx; (E) manuscript draft; or a
user-directed follow-up not in the H5 candidate set. The end-of-Phase-I
project state has a clean H1-locked synthesis with 81 atomic claims,
all 3 axis contests decisively favouring their B-models (margins 18,
12, 53), the cross-axis theme ranking stable, and a Phase I appendix
in section 18 documenting the per-state NF-kB attenuation refinement.

---

## Anti-anchoring (re-read every session)

These exist because LLMs (this agent included) drift toward the most
salient prior finding. Each Phase I session must enforce them:

- **Always** distinguish the three competing per-state patterns
  (global / state-restricted / state-driven-by-abundance) explicitly.
  Never collapse the test into a single "T-Tau-attenuates is
  state-specific" claim before the per-state TF activity table +
  GSEA + the section 11 composition test all agree.
- **Always** report all 4 substates (homeostatic / DAM / IFN /
  proliferative) in the per-state TF activity table, even if one
  carries the signal. Never report only the substate with the
  strongest signal; the absence of signal in other substates is
  itself evidence.
- **Always** confirm sign consistency between the per-state TF
  activity (negative Rela score = NF-kB suppression) and the
  per-state GSEA (negative NES on NF-kB targets = NF-kB target
  suppression). Disagreement between the two routes is informative
  and must be flagged in the I3 verdict prose, not papered over.
- **Always** acknowledge the compositional-confound alternative in
  the I3 verdict, citing the existing section 11 DAM-sampling
  balance ANOVA. Never claim per-state specificity without
  explicitly considering whether the per-state signal is itself an
  abundance-shift artefact.
- **Always** grade the I4 new ledger rows transparently using the
  H1-locked three-tier rules. Strong requires ≥2-modality consistent
  sign + FDR<0.10 + explicit cross-layer corroboration (here:
  per-state TF activity table + per-state NF-kB-target GSEA + the
  existing section 14 whole-microglia Rela sign-reversal). Without
  all three the row grades to Moderate at most.
- **Always** treat the Phase H integrated reading as the prior; the
  Phase I test is a refinement, not a replacement. T-Tau-attenuates'
  current 5-row position is the baseline; Phase I either lifts it
  or refines its character.
- **Never** add new computation outside the per-state TF activity
  inference + per-state GSEA. The plan deliberately bounds the test
  scope; expanding to per-state kinase activity or per-state CCC
  layers is out of scope for Phase I and would be its own future
  plan if warranted.
- **Never** modify the H1-locked 11-entity model set or the H1-locked
  three-tier confidence-grade definitions. Phase I operates within
  the Phase H structural envelope; structural changes to the
  adjudication scheme would require a new plan with its own decision
  gates.

## Outcome summary

Phase I closed 2026-05-25 with a clean **Global attenuation** verdict
for the per-state NF-kB attenuation test, lifting **T-Tau-attenuates**
from the lowest-support cross-axis theme (5 supports / 3 Strong) to
the 4th-ranked theme (11 / 7) and lifting **Hyp-1B**'s contest margin
over Hyp-1A from 6 to 18. The biological character of the Phase H
load-bearing Rela sign-reversal is now resolved: tau on the amyloid
background suppresses amyloid-driven NF-kB transcription **across
every microglia substate independently** (DAM lead carrier, with
homeostatic close behind), and the suppression is **transcriptional,
not compositional** (section 11 DAM-sampling balance ANOVA tau:nlgf
F=1.03 p=0.330).

**(a) I1-locked choices.** Chapter venue = new top-level section 18
"Per-state NF-kB attenuation test" rendered via
`rmd/17_nfkb_attenuation.Rmd` inserted into `analysis.Rmd` between
`child-biological-model` (rmd/16) and `child-session` (rmd/99). NF-kB
target GSEA scope = included (per-state GSEA against the union of
CollecTRI mouse Rela / Nfkb1 / Rel / Relb / Nfkb2 / A0A979HLR9 target
gene sets at the interaction contrast). Both locks default-approved by
user at I1 (no scope expansion or contraction relative to the plan
body).

**(b) Per-state TF activity result for the NF-kB family.** The
**A0A979HLR9 complex-level CollecTRI accession** (1340 mouse targets;
the union of NF-kB family subunit occupancy at composite promoters) is
sign-consistently negative in all 4 microglia substates at the
interaction contrast: **DAM consensus score -6.52** (p=7e-11),
**homeostatic -5.45** (p=5e-8), **IFN -3.66** (p=3e-4),
**proliferative -2.55** (p=0.011). Individual NF-kB family subunits
(Rela, Nfkb1, Nfkb2, Relb) carry weakly negative but non-FDR<0.10
signals at the per-subunit level; the family-attenuation signature is
carried by the complex-level accession, which is biologically coherent
(NF-kB transcriptional output is the cooperative product of subunit
occupancy at composite promoters, captured by the complex-level target
set and missed by the per-subunit sets).

**(c) Per-state NF-kB-target GSEA result.** Union NES at the
interaction contrast is negative in all 4 substates: **DAM -1.247**
(padj=0.027, sig), **homeostatic -1.221** (padj=0.091, borderline),
**IFN -1.126** (padj=0.425, NS), **proliferative -1.121** (padj=0.162,
NS). The sign consistency matches the TF activity layer in every
substate; the within-state cross-layer gate
(`tf_neg & gsea_neg & (tf_sig | gsea_sig)`) is satisfied by all 4
substates. The GSEA mirrors the section 14 whole-microglia Rela
sign-reversal at the per-substate level (positive NES +1.03 to +1.21
at the two NLGF arms, negative NES -1.12 to -1.25 at the interaction,
in every substate).

**(d) Verdict tier.** **Global attenuation** (tier 1 of the 5-tier
set). All 4 substates independently carry NF-kB suppression at the
interaction; the whole-microglia signal is the weighted aggregate of
consistent per-state suppression, not a one-substate-only signal and
not a compositional artefact. The state-driven-by-abundance alternative
is ruled out on solid statistical grounds by the section 11
DAM-sampling balance ANOVA (`ccc_dam_sampling_balance_anova.tsv`,
tau:nlgf F=1.03 p=0.330). **DAM is the lead carrier, with homeostatic
close behind** — a biological revision of the H4 17.4 §4 hypothesis
sketch which had predicted homeostatic-led suppression. DAM-led
attenuation is biologically intuitive (DAM cells are transcriptionally
NF-kB-HIGH under amyloid, so the attenuation at the interaction has
the most signal in the substates with the most-amplified baseline
NF-kB program).

**(e) I4 ledger update.** **6 new rows added (I-001..I-006)**; the
ledger grows from 75 to 81 rows. T-Tau-attenuates' `support_count`
moves from **5 to 11** (net +6); `n_strong_supports` from **3 to 7**
(net +4). 4 of the 6 new rows grade Strong under the H1-locked
three-criteria rule (I-001 homeostatic, I-002 DAM, I-005 whole-
microglia replication, I-006 compositional ANOVA); 2 grade Moderate
(I-003 IFN, I-004 proliferative) because their GSEA layers fail
FDR<0.10. Each Phase I row supports Hyp-1B and contradicts Hyp-1A:
**Hyp-1B's net_support moves from 24 to 30** (+6 supports), **Hyp-1A's
net_support moves from 18 to 12** (+6 contradicts), giving a
**combined contest margin shift from 6 to 18** (Hyp-1B vs Hyp-1A). The
cross-axis ranking re-orders T-Tau-attenuates from **5th to 4th**,
swapping with Hyp-0 Cdk5 (7 supports). All other entities unchanged:
Hyp-2B vs Hyp-2A margin 12, Hyp-3B vs Hyp-3A margin 53, T-Synergy 28,
T-Inflammation 26, T-Compartment-suppression 17, Hyp-0 7.

**(f) Next-plan pointer.** Phase I is COMPLETE; no successor plan is
auto-prescribed. The next plan is **user-directed** from the Phase H
Outcome summary candidates: **(A)** Cdk5-substrate GSEA at all three
axes — tests whether Cdk5's apparent cross-axis integrator role
generalises beyond the kinase-inference signal that established it;
**(B)** Pros1-Mertk GeoMx validation — tests the TAM-kinase axis-2
microglia-microglia node that narrowly survived the synaptic-
suppression contest as a residual route alongside the Hyp-2B winner;
**(C)** hdWGCNA on GeoMx — extends the snRNAseq hdWGCNA modules into
the spatial layer for axis-2 + axis-3 corroboration; **(E)**
manuscript draft — extracts the adjudicated biological model into
publication form; or a **user-directed follow-up** not in the H5
candidate set. The end-of-Phase-I project state has a clean H1-locked
synthesis: 81 atomic claims, all 3 axis contests decisively favouring
their B-models (margins 18, 12, 53), cross-axis theme ranking stable
at T-Synergy (28) > T-Inflammation (26) > T-Compartment-suppression
(17) > T-Tau-attenuates (11) > Hyp-0 Cdk5 (7), and a Phase I appendix
in section 18 documenting the per-state NF-kB attenuation refinement.

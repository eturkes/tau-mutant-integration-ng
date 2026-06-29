# Cell-type specificity of the tau×amyloid interaction (arc N)

Active plan. Read the STATUS block, then execute ONLY the next TODO step's body
(read other step bodies only if the current step depends on them). House loop:
mark the step DONE + write a completion note, chown new files, knit-verify
(0 `class="error"` / 0 `class="warning"`), commit locally.

## Premise (why this arc)

Every prior arc (D→M) reads ONLY microglia. The project's entire edifice —
DAM programme, the locked 2×2 contests (amyloid 18 / synaptic 12 / interaction
55), the §23 progression interaction (M-002, the first positive+significant one)
— is microglia-scoped by construction. **Untested question: is the tau×amyloid
interaction (and the amyloid/tau main effects) microglia-SPECIFIC, or do
astrocytes / oligodendrocytes / OPCs / neurons / vascular cells show comparable
factorial responses?** A non-microglial interaction would reframe the project's
microglia focus as a *choice*, not a data-forced conclusion; microglia-restriction
would *validate* it. Either outcome is a finding. This is the project's first
CROSS-CELL-TYPE readout. Display-only chapter; does NOT touch the locked mouse
2×2 verdicts except via the gated N5 ledger-feed.

## Grounded facts (verified against live caches this session — do not re-derive)

- Substrate = `storage/cache/seurat_full_processed.rds` (286,285 cells; loads in
  ~33 s, peak RSS well under the 62 G box). Metadata carries `genotype`, `batch`,
  `sex`, `genotype_batch`, `cell_type` (9-level), `broad_annotations`.
- `cell_type` = `broad_annotations` for non-microglia, `Microglia_<state>` for
  microglia. The **5 broad non-microglial types**: `Astrocyte`, `Neuronal`,
  `Oligodendrocyte`, `OPC`, `Vascular`. Pool the 4 `Microglia_*` labels →
  `Microglia` (6th comparison unit, microglia-as-whole).
- Cells per (cell_type × genotype_batch replicate) — **all 6 units have all 16
  replicates present**:
  | unit | total | min/rep | median/rep |
  |---|---|---|---|
  | Neuronal | 162,954 | 5,674 | 9,758 |
  | Oligodendrocyte | 47,988 | 557 | 3,134 |
  | Astrocyte | 24,894 | 372 | 1,422 |
  | Vascular | 16,192 | 364 | 987 |
  | OPC | 8,153 | 294 | 540 |
  | Microglia (pooled) | 26,104 | 289 | ~1,631 |
  → Pseudobulk is well-powered for ALL six. **THE CONFOUND: cell count varies
  ~25× (Neuronal 163k vs IFN-thin microglia), so pseudobulk precision ⇒ power is
  asymmetric. A specificity claim from NATIVE counts alone is confounded.** The
  load-bearing control is a cell-count-MATCHED (downsampled) refit (locked below).

- **N2-verified (do not re-derive):** (1) **K = 289** for the MATCHED regime —
  `min` over the 6 units of (min cells-per-replicate); the binding unit is
  pooled-Microglia (289 < OPC's 294), resolving the table's "verify in N3" cell.
  (2) `seurat_full_processed.rds` is **SYMBOL-keyed** (RNA rownames are MGI
  symbols, e.g. `Mob3b`/`Zbtb16`), NOT Ensembl like the microglia caches; RNA
  has a single joined `counts` layer (GetAssayData works on subsets). Arc-N DE
  therefore runs with `symbol_map = NULL` (the `gene` column already IS the
  symbol); `microglia_crosscheck` bridges both sides into SYMBOL space before
  comparing to the Ensembl-keyed `de_snrnaseq`/`de_snrnaseq_nebula`.
- **Cross-check preview (N2, pseudobulk):** seurat_full pooled-Microglia pb vs
  `de_snrnaseq` gives ρ(logFC) = tau_alone 0.962, nlgf_in_maptki 0.958,
  nlgf_in_p301s 0.984, tau_in_nlgf 0.958, **interaction 0.926**. The interaction
  sits just under the plan's blanket 0.95 — a GENUINE property (0 ambiguous
  symbols, so not a collapse artifact; it is the noisiest difference-of-
  differences contrast across two intentionally non-identical microglia analyses
  per the locked microglia-QC re-cluster). 0.926 still confirms reproduction
  (random ≈ 0). **N3 should report all 5 ρ with a gross-mismatch floor (~0.90),
  not a hard 0.95-per-contrast assert** (open for user steer before N3 locks it).

## Locked design (anti-anchoring: fixed BEFORE inspecting which types light up)

Mirror `rmd/02b_snrnaseq_de.Rmd` (which runs BOTH estimators on microglia) so all
six units are method-comparable. **TWO estimators × TWO power regimes, all 6 units
(= the locked N1 answer: pseudobulk + NEBULA, full stack):**
- **Estimator A — pseudobulk** (`genotype_batch`, 16 ids): `build_pseudobulk(sub,
  sample_col="genotype_batch", covariate_cols=c("genotype","batch","sex"))` →
  `fit_limma_voom(min_count=10)`, `design <- model.matrix(~ 0 + genotype + batch)`
  (strip `genotype` prefix, sanitise names), `make_contrast_matrix` → the 5
  canonical contrasts incl `interaction`. Annotate with `symbol_map`.
- **Estimator B — single-cell NEBULA** (full power; the microglia-HEADLINE method):
  reuse the GENERIC `fit_nebula_microglia(sub, id_col="genotype_batch",
  genotype_col="genotype", min_cell_frac=0.01, ncore=detectCores()-2,
  symbol_map=symbol_map)` — VERIFIED generic despite the name (builds the 2×2
  factorial internally, returns the same `$top[[contrast]]` shape as
  `de_snrnaseq_nebula`). Optionally alias `fit_nebula_celltype`.
- **Regime 1 — NATIVE** (all cells): pseudobulk fast; NATIVE NEBULA full-power but
  power-ASYMMETRIC (Neuronal 163k ≫ microglia) — report cells/unit + convergence.
  **NATIVE NEBULA on Neuronal is the compute bottleneck (≈ hours, chunked over
  workers) → run the build backgrounded/resumable.**
- **Regime 2 — MATCHED** (downsample each unit×genotype_batch to `K` cells, where
  `K = min over the 6 units of (min cells-per-replicate)`, ≈290–300; fixed RNG
  seed; sampling without replacement; RECORD K+seed): power-EQUALISED for BOTH
  estimators. MATCHED-NEBULA is CHEAP (≈K×16 ≈ 4.6k cells/unit). **The specificity
  verdict's headline rests on the MATCHED regime (both estimators agree ⇒ robust);
  NATIVE is shown but flagged power-confounded.**
- **Microglia cross-check (sanity gate):** NATIVE pseudobulk `Microglia` (pooled)
  vs existing `de_snrnaseq`, AND NATIVE NEBULA `Microglia` vs `de_snrnaseq_nebula`
  — assert Spearman(logFC) > 0.95 on shared genes per contrast (report ρ; not
  identical — seurat_full microglia may differ slightly from
  `microglia_seurat_processed`).

## Readouts (defined before inspection; report all 6 units SYMMETRICALLY)

R1. **Significant-gene tally** per (unit × contrast × estimator{pseudobulk,NEBULA}):
    n at FDR<0.05 and <0.10, up/down. Headline = the `interaction` column. NATIVE.
R2. **Power-matched tally** = R1 recomputed on the MATCHED fits (both estimators).
    **Load-bearing**: the specificity verdict reads from R2 (R1 shown but flagged
    confounded).
R3. **Interaction-effect concordance**: (a) genome-wide Spearman of `interaction`
    logFC between `Microglia` and each other unit (low ρ ⇒ distinct interaction
    response), per estimator; (b) a curated cross-type `interaction`-logFC heatmap
    on the microglia headline gene set (DAM programme + the M-002 progression
    drivers + the H-phase headline interaction genes — pull symbols from existing
    results/`trajectory_*`, `de_snrnaseq_nebula$top$interaction` top-N);
    (c) **cross-ESTIMATOR** pseudobulk-vs-NEBULA `interaction` concordance per unit
    (the 02b "pseudobulk vs NEBULA power gain" panel, generalised to 6 units).
R4. **Per-gene specificity class** for genes sig at `interaction` (FDR<0.10) in
    ANY unit: classify microglia-unique / shared / non-microglial-unique (UpSet-
    or tally-style). Both power regimes.
R5. (secondary, optional) **Pathway-level**: GO BP fGSEA per unit at `interaction`
    (reuse `R/fgsea.R` getters + `run_fgsea_*`); is the microglial interaction
    enrichment cell-type-restricted? Defer if N1=pseudobulk-only keeps scope tight.

## Anti-anchoring guardrails (enforce every step)

- Metrics R1–R5 are fixed above BEFORE running any fit. Do not retune after seeing
  which types light up.
- **Do NOT assume microglia win.** Report all 6 units in every table/figure with
  no privileged ordering (alphabetical or by total-cells, stated). A non-microglial
  interaction hit — especially **Astrocyte** (an established AD responder) — is a
  FINDING, reported plainly, not noise to explain away.
- The specificity verdict MUST cite the **MATCHED** (power-equalised) result; a
  NATIVE-only claim is confounded by the ~25× cell-count spread and is disallowed
  as the headline.
- If a non-microglial unit shows interaction breadth ≥ microglia under MATCHED
  power, state that it reframes (not refutes) the microglia focus; never bury it.

## Decision gates

- **N1 (START gate — method depth):** present default + alternative, then
  AskUserQuestion. **Default = pseudobulk-only** (mirrors 02b primary; fast;
  the cell-count confound is directly controllable via the MATCHED downsample →
  cleanest power-equalised specificity claim). **Alternative = add per-type
  single-cell NEBULA** for full-power interaction (matches the microglia *headline*
  method, but heavy — NEBULA on 163k Neuronal cells ≈ hours — AND reintroduces
  power asymmetry that the pseudobulk MATCHED regime avoids). Recommend default.
- **N5 (END gate — ledger feed):** decided AFTER results, L5/M5-style. A
  specificity result is META (about the validity of the microglia focus), only
  weakly model-discriminating. Options to present then: (a) no feed (framing, not
  a model move — the honest default if microglia-restricted as expected);
  (b) margin-neutral Suggestive row(s) contextualising the microglia-centric
  reading; (c) if a non-microglial interaction is found, a row recording
  cross-cell-type breadth. Keep the 11-entity set + grade defs untouched; engineer
  any feed margin-neutral (contests stay 18/12/55).

## Artifacts (paths to create)

- `R/celltype_specificity.R` — helpers: `subset_pseudobulk_de(sc, unit)` (subset →
  build_pseudobulk → fit_limma_voom, mirrors 02b), `subset_nebula_de(sc, unit)`
  (subset → `fit_nebula_microglia`, the generic single-cell estimator),
  `downsample_balanced(sc, units, K, seed)` (per unit×replicate cell sampling to K),
  `assemble_specificity_tables()` (R1–R4), plot fns (tally heatmap,
  interaction-logFC cross-type heatmap, cross-estimator panel, per-gene
  specificity). Source in `R/helpers.R` AFTER `trajectory.R` (verify exact slot in
  N2; depends only on de_pb.R + design.R + de_sc.R, all early).
- `scripts/build_celltype_specificity.R` — reads `seurat_full_processed.rds`,
  loops the 6 units × {pseudobulk, NEBULA} × {NATIVE, MATCHED}, assembles R1–R5 →
  `celltype_specificity.rds` + results TSVs. Idempotent; params (K, seed,
  min_count, min_cell_frac) embedded + echoed. **NATIVE NEBULA on Neuronal (163k
  cells) is the heavy step (≈ hours) → run backgrounded; consider caching per-unit
  NEBULA fits so a re-run resumes** (the build is not auto-resumable like SCENIC —
  guard each unit's NEBULA in its own `cache_or_run`/file or run order
  heaviest-last).
- `rmd/23_celltype_specificity.Rmd` — display-only, reads the cache; renders as
  **§24** (file≠§ — verify the offset in N4). Wire `child-celltype-specificity`
  into `analysis.Rmd` BEFORE `child-session`.
- Cache `storage/cache/celltype_specificity.rds`; results
  `storage/results/celltype_specificity_{tally_native,tally_matched,interaction_concordance,specificity_class,verdict}.tsv`.

## STATUS

- [x] **N0** design + grounding (cell counts verified; 02b/de_pb.R/design.R/de_sc.R infra mapped; confound identified) — DONE 2026-06-06
- [x] **N1** GATE method depth — DONE 2026-06-06: user chose **Pseudobulk + NEBULA (all)** (full-power stack, both estimators × NATIVE+MATCHED on all 6 units; mirrors the microglia 02b treatment). NEBULA reuses the generic `fit_nebula_microglia` per unit subset (verified generic).
- [x] **N2** `R/celltype_specificity.R` helper + smoke test — DONE 2026-06-06: 11 helpers + 4 plot fns written, sourced after trajectory.R; smoke-tested on OPC (both estimators) + downsample determinism + R1–R4 tables + plots. Verified K=289 and the symbol-keyed substrate; previewed the pseudobulk cross-check (4/5 ρ>0.95, interaction 0.926). See N2 completion note + grounded-facts updates.
- [x] **N3** `scripts/build_celltype_specificity.R` → `celltype_specificity.rds` + 6 results TSVs — DONE 2026-06-06: full 6-unit × {pseudobulk,NEBULA} × {NATIVE,MATCHED} stack + R1–R5 built; NATIVE NEBULA Neuronal (163k cells) ran in ~25 min via per-unit resume caches; all 10 cross-check ρ PASS the 0.90 floor (nebula interaction 0.917 nearest). See N3 completion note.
- [x] **N4** `rmd/23_celltype_specificity.Rmd` (§24) + wire child + clean re-knit — DONE 2026-06-06: display-only chapter built (7 sections, reads cache only), `child-celltype-specificity` wired between trajectory + session; file 23 renders as **§24** (session now §25); knit clean (0 error / 0 warning); `celltype_specificity_verdict.tsv` written. See N4 completion note.
- [x] **N5** END GATE ledger-feed decision + feed (if any) + refresh map.md/DIGEST + archive plan (date suffix) — DONE 2026-06-06: user chose **option C "feed both halves"**; 2 margin-neutral Suggestive rows N-001/N-002 added (`phase_n` block, `layer="specificity"`), ledger 89→91, contests UNCHANGED 18/12/55; prose synced in rmd/16 + analysis.Rmd; re-knit clean (0/0); map.md + DIGEST refreshed; plan archived. See N5 completion note.

## Execution model

- N1 is a START gate: present the default + the NEBULA alternative in prose, then
  AskUserQuestion; record the answer in this STATUS block before N2.
- N2 (new helper) and N3 (new cache shape) are each design-heavy → ideally a fresh
  session each; smoke-test helpers via `Rscript -e` against the live cache before
  the full build (build loads 1.9 G + runs 12 DE fits — minutes, not seconds).
- N4 is a knit step: confirm `grep -c 'class="error"' analysis.html` == 0 (same for
  warning) after rendering.
- N5 is the END gate: present the 3 ledger options with the post-hoc result, then
  AskUserQuestion; only then refresh notes + archive.
- Every step: chown rstudio:rstudio new files (you run as root) incl knit outputs;
  commit locally with an imperative subject + HEREDOC body + Co-Authored-By trailer.

## N2 completion note (2026-06-06)

Built `R/celltype_specificity.R` (sourced in `R/helpers.R` immediately after
`trajectory.R`; depends only on de_pb.R / design.R / de_sc.R, all earlier).
Functions: `celltype_unit_labels` (pool the 9-level `cell_type` → 6 units),
`min_cells_per_replicate` (K finder), `subset_pseudobulk_de` (byte-for-byte the
rmd/02b limma-voom chunk on a unit subset), `subset_nebula_de` (wraps the generic
`fit_nebula_microglia` per subset), `downsample_balanced` (seeded per
unit×replicate sampling to K), the R1–R4 table fns (`specificity_tally`,
`interaction_concordance`, `cross_estimator_concordance`,
`specificity_class_table`), `microglia_crosscheck` (symbol-bridged sanity gate),
`assemble_specificity_tables` (orchestrator over the regime×estimator×unit grid),
and four ggplot display fns (`plot_specificity_tally_heatmap`,
`plot_interaction_logfc_heatmap`, `plot_cross_estimator_panel`,
`plot_specificity_class`). Both DE wrappers lift the per-contrast tables to a
uniform top-level `$top[[contrast]]` so pseudobulk and NEBULA share one accessor.

Smoke-tested via `scripts/smoke_test_celltype_specificity.R` against the live
`seurat_full_processed.rds` (286,285 cells): helpers parse + source cleanly; the
6 units and per-replicate counts reproduce the plan's grounded facts EXACTLY
(Microglia pooled = 26,104; K = 289). `subset_pseudobulk_de("OPC")` → the 5
canonical contrasts with the 02b design (MAPTKI/P301S/NLGF_MAPTKI/NLGF_P301S + 3
batch cols), 8,153 cells, 11,694 genes, 6.8 s, symbol == gene (no NAs).
`subset_nebula_de("OPC", frac=0.05)` → the 5 contrasts, 6,557/6,635 converged,
26 s. `downsample_balanced(K=50)` is reproducible (identical cells across two
seeded calls) and respects the per-stratum cap. `assemble_specificity_tables`
returns tally/concordance/specificity_class with the expected columns/row counts;
all four plot constructors return ggplot objects.

Deviation flagged for N3 (see grounded-facts "Cross-check preview"): the headline
**interaction** cross-check is ρ=0.926, just under the plan's blanket 0.95 — a
genuine noisy-contrast property, not a wiring bug. N3 should report all 5 ρ with a
~0.90 gross-mismatch floor rather than a hard 0.95-per-contrast assert; user steer
invited before N3 hard-codes the floor. NATIVE-NEBULA timing intel for N3: OPC
(8,153 cells, 6,635 genes @ frac 0.05) fit in 26 s on 6 workers, so Neuronal
(162,954 cells) at frac 0.01 is the bottleneck (~tens of minutes, background it).

## N3 completion note (2026-06-06)

**Built.** `scripts/build_celltype_specificity.R` (thin I/O driver) + 3 R5 helpers
added to `R/celltype_specificity.R` (`specificity_pathway_fgsea`,
`specificity_pathway_tally`, `microglia_pathway_cross_unit`). Two user gates locked
before the ~30-min compute: **(1)** include R5 (pathway-level fGSEA); **(2)** a
single uniform **0.90** cross-check floor (not 0.95-per-contrast). The build ran
backgrounded; NATIVE NEBULA is guarded per-unit in
`storage/cache/celltype_specificity_fits/native_nebula_<unit>.rds` (lightest-first
order so a Neuronal crash leaves 5 units cached) — Neuronal (163k cells) fit in
~25 min, total run ~37 min. Params: `SEED=1`, `K=289` (MATCHED downsample, the
binding unit is pooled-Microglia), `MIN_COUNT=10`, `MIN_CELL_FRAC=0.01`,
`PADJ_CUT=0.10`, `XCHECK_FLOOR=0.90`.

**Cache** `storage/cache/celltype_specificity.rds` (46.7 MB, xz) — 11 slots: `fits`
(24 fits = NATIVE/MATCHED × pseudobulk/NEBULA × 6 units, each slimmed via
`slim_fit` to `top`/`n_cells`/`n_genes`), `tally` (120 rows), `interaction_concordance`,
`specificity_class` (21 rows), `pathway_tally` (24 rows),
`microglia_pathway_cross_unit` (1035 rows / 139 microglia-sig pathways), `crosscheck`
(10 rows), `headline_genes` (39 symbols), `unit_cell_counts`, `per_unit_min`,
`params`. Six results TSVs written (tally_native, tally_matched,
interaction_concordance, specificity_class, pathway_tally, crosscheck).

**Cross-check — all 10 PASS 0.90** (Spearman logFC vs canonical de_snrnaseq[_nebula],
symbol-bridged). pseudobulk: tau_alone 0.962, nlgf_in_maptki 0.958, nlgf_in_p301s
0.984, tau_in_nlgf 0.958, **interaction 0.926**. nebula: 0.959, 0.964, 0.978, 0.960,
**interaction 0.917**. The nebula interaction (0.917) is nearest the floor —
**vindicating the user's 0.90 choice**: a 0.95 assert would have spuriously failed
BOTH interaction cross-checks despite genuine reproduction. (The interaction is the
noisiest difference-of-differences contrast across two intentionally non-identical
microglia analyses; random ≈ 0, so 0.92 confirms reproduction.)

**R1 NATIVE interaction (gene-level, FDR<0.10)** — sparse everywhere: Microglia
4 (nebula) / 0 (pb), Oligodendrocyte 5 (neb) / 6 (pb), Vascular 1, rest ≈0;
tau_alone ≈0 everywhere. Main effects (nlgf_in_*) are broad and multi-cellular
(e.g. nlgf_in_maptki nebula: Microglia 1529, Astrocyte 950) — the factorial DESIGN
works; only the difference-of-differences interaction is thin.

**R2 MATCHED interaction (LOAD-BEARING, power-equalised at K=289)** — the gene-level
interaction collapses to **0 in ALL six units including Microglia**, the sole
exception Oligodendrocyte 6 (nebula only). Main effects survive matching. ⇒ At equal
power there is **no microglia-unique gene-level interaction signal**; the NATIVE
microglia hits are a power artefact of the ~25× cell-count spread.

**R3 microglia-vs-unit interaction concordance** — low genome-wide ρ (matched
0.08–0.15, native 0.14–0.30) ⇒ the per-unit interaction logFC profiles are
**distinct, not a shared programme**. Cross-estimator pb-vs-nebula agreement is
0.91–0.98 per unit (the two estimators agree; it is the cell-types that differ).

**R4 specificity_class** (genes sig at interaction in ANY unit; tiny 21-row pool):
NATIVE nebula = 4 microglia_unique (Gm30211, Cd276, Cmah, Plac9) + 5
non_microglial_unique (mostly Oligo); NATIVE pb = 6 non_microglial_unique (all Oligo);
MATCHED nebula = 6 non_microglial_unique (all Oligo). **No microglia-unique genes
survive matching.**

**R5 pathway-level (GO BP fGSEA at interaction, FDR<0.10)** — **Neuronal dominates
every regime×estimator.** MATCHED/nebula: Neuro 210, Micro 2, all others 0.
MATCHED/pb: Neuro 242, Vasc 123, Micro 37, OPC 38, Astro 16, Oligo 13. NATIVE/nebula:
Neuro 872, Oligo 164, Vasc 130, Astro 94, Micro 55, OPC 72. NATIVE/pb: Neuro 1426,
Vasc 427, Oligo 274, Astro 208, Micro 83, OPC 53. Pathway enrichment is overwhelmingly
DOWN (sign-consistent with the interaction shrinking neuronal programmes).

**HEADLINE.** The tau×amyloid interaction is **NOT microglia-specific.** At
load-bearing MATCHED power the gene-level interaction is essentially absent
everywhere; at pathway level the broadest response is **Neuronal**, not Microglia.
This **REFRAMES (does not refute)** the project's microglia focus: microglia were a
*choice*, and the strongest factorial pathway response actually sits in neurons.

**Caveats for N4 (must surface in the §24 chapter).** (1) **Depth confound**: MATCHED
equalises CELL COUNT, not transcriptome depth — neurons carry ~2–3× UMIs/cell, so the
neuronal fGSEA dominance is partially depth-driven; state this beside every "Neuronal
broadest" claim. (2) The MATCHED gene-level collapse to ~0 is the honest power-equalised
verdict; NATIVE tallies are shown but flagged confounded. (3) R5 was added per the N1-era
user gate (originally "optional"); it is now a headline readout. (4) fgsea emits benign
"P-values likely overestimated" multilevel warnings — out-of-knit here, but if any fGSEA
runs inside rmd/23 they must be `suppressWarnings`/`message`-routed to keep the §24
warning gate at 0.

**Next (N4).** Build `rmd/23_celltype_specificity.Rmd` (display-only, reads the cache;
renders as §24 — verify the file≠§ offset), wire `child-celltype-specificity` into
`analysis.Rmd` before `child-session`, clean re-knit (0 error / 0 warning), write
`celltype_specificity_verdict.tsv`. Design-heavy → ideally a fresh session. N5 (END
gate) decides the ledger feed AFTER the chapter, with the post-hoc result in hand.

## N4 completion note (2026-06-06)

**Built.** `rmd/23_celltype_specificity.Rmd` — display-only chapter, reads
`celltype_specificity.rds` only (no fits or fGSEA run inside the knit, so the
warning gate stays 0). Wired `child-celltype-specificity` into `analysis.Rmd`
between `child-trajectory` and `child-session`. **File 23 renders as §24** (the
file≠§ offset confirmed in the rendered HTML: §24 = "Cell-type specificity", §25 =
"Session info"); map.md's "rmd/22→§23" anchor extends cleanly. Knit clean:
`grep -c 'class="error"'` = 0 and `class="warning"` = 0. `chown rstudio:rstudio`
applied to the new rmd, the edited `analysis.Rmd`, `analysis.html`, and the new
verdict TSV (the html / results / cache are all gitignored, so only `analysis.Rmd`
+ `rmd/23` are tracked changes).

**Chapter structure (7 sections, mirrors the rmd/22 trajectory house style).**
Intro prose (premise / the ~20× cell-count confound / locked 2-estimator ×
2-regime × 6-unit method / thresholds / caveats / anti-anchoring); a hidden
top-of-file load chunk (`celltype-specificity-load`, `include = FALSE`) holding
every prose scalar via cache accessors `tget/xget/pget`; then **Sanity gate**
(10-ρ cross-check kable, all PASS 0.90), **Significant-gene tally** (R1/R2 — the
load-bearing MATCHED tally heatmap + the power-confounded NATIVE one),
**Interaction concordance** (R3a low microglia-vs-unit ρ in prose, R3c
cross-estimator scatter, R3b curated headline-gene interaction heatmap),
**Per-gene specificity class** (R4 native + matched bars), **Pathway-level**
(R5 inline unit × regime×estimator pathway-count heatmap), and **Verdict**
(writes `celltype_specificity_verdict.tsv` — 5 readout rows × matched/native/
verdict columns — + a kable). The four plot helpers from `R/celltype_specificity.R`
are called directly; the only inline ggplot is the R5 pathway heatmap.

**Headline rendered (all numbers cache-read, no drift).** The tau×amyloid
**interaction is NOT microglia-specific.** At load-bearing MATCHED power
(K=289/replicate) the gene-level interaction collapses to 0 in all six units
including Microglia — sole survivor Oligodendrocyte (6 genes, NEBULA), i.e. the
only power-equalised gene-level interaction is *non*-microglial. The amyloid MAIN
effect, by contrast, survives matching and stays microglia-led (nlgf_in_maptki
NEBULA: Microglia 625 vs Astrocyte 218), so the chapter draws the key
distinction: the *amyloid response* validates the microglia focus, but the
*divergence readout* (interaction) does not. R3: per-unit interaction profiles
are mutually distinct (matched ρ 0.08–0.15) while estimators agree (pb-vs-NEBULA
0.96–0.98) — cell types differ, not methods. R4: 0 microglia-unique interaction
genes survive matching (the 4 NATIVE ones — Gm30211, Cd276, Cmah, Plac9 — are the
power artefact the MATCHED control was built to expose). R5: Neuronal dominates
the pathway interaction in every regime×estimator (matched NEBULA Neuro 210 vs
Micro 2; matched pb Neuro 242), enrichment overwhelmingly down. **Reframes, does
not refute.**

**Deviations / decisions.** (1) The load chunk was moved to the very top of the
file (`include = FALSE`) after the first knit failed: the intro prose reads
`par`/`ucc` inline *before* the body, so the chunk must run first — pure ordering
fix, no content change. (2) The depth confound was verified live rather than
trusted from the N3 note: Neuronal median UMI 2,809 vs Microglia 1,064 = **2.64×**
(the recorded "~2–3×" confirmed exactly); this concrete figure is stated beside
every neuronal-dominance claim, and the neuronal pathway dominance is explicitly
read only as "broadest outside microglia", never a quantitative neuronal-magnitude
claim. (3) All four N3 caveats are surfaced in the chapter (depth confound;
MATCHED load-bearing vs NATIVE flagged confounded; R5 promoted from optional to a
headline readout; fGSEA kept out-of-knit). (4) The uniform 0.90 cross-check floor
is narrated as the calibrated choice (a 0.95-per-contrast assert would have
spuriously failed both interaction checks at ρ 0.917/0.926).

**Next (N5 — END GATE).** Decide the ledger feed with the post-hoc result in
hand. Present the three plan options — (a) no feed (framing-only; the honest
default given the result is META about the microglia focus, not model-
discriminating); (b) margin-neutral Suggestive row(s) contextualising the
microglia-centric reading; (c) a row recording the cross-cell-type breadth of the
interaction (neuronal pathway dominance / non-microglial gene-level survivor) —
each engineered margin-neutral so the locked contests stay 18/12/55 and the
11-entity set + grade defs are untouched. This is a user-gated decision
(AskUserQuestion). After the gate: apply any feed, refresh map.md (add the rmd/23
row + the §24 anchor + the verdict TSV) and completed/DIGEST.md, then archive this
plan with a YYYY-MM-DD suffix.

## N5 completion note (2026-06-06)

**Gate (AskUserQuestion).** Presented the three options against the post-hoc
result (OVERALL: the interaction is NOT microglia-specific, but the amyloid main
effect IS microglia-led). User chose **option C — "feed both halves"**: record
BOTH the reassuring scope-validation AND the central reframing, avoiding the
anti-anchoring trap of feeding only the comfortable half.

**Feed (2 margin-neutral Suggestive rows; `phase_n` block in
`scripts/build_biological_model_ledger.R`, new `layer="specificity"`).**
- **N-001** `amyloid_activation` — the amyloid MAIN effect is microglia-led at
  MATCHED K=289 in both tau backgrounds and both estimators (NEBULA Microglia 625
  vs Astrocyte 218 in `nlgf_in_maptki`; 324 vs 210 in `nlgf_in_p301s`) despite the
  2.64× depth disadvantage. Supports `Hyp-1A;Hyp-1B;T-Inflammation`, empty
  contradicts → **mirrors K-001/L-001/M-001**: amyloid margin holds at **18**; it
  lifts T-Inflammation **30→31** (a raw-count tie-break over T-Synergy 30).
- **N-002** `interaction_metabolic` — the interaction is NOT microglia-specific
  (MATCHED gene-level collapse to 0 in all 6 units incl Microglia, sole survivor
  the non-microglial Oligodendrocyte; per-unit ρ 0.08–0.15 distinct; broadest
  pathway response Neuronal). **Empty supports + empty contradicts ("feeds no
  model", the cross-cell-type counterpart to L-002)** → interaction margin holds
  at **55**, themes untouched. Deliberately NOT tagged T-Synergy: the per-unit
  profiles are mutually distinct, so a T-Synergy tag would overclaim a SHARED
  cross-cell-type mechanism (unlike M-002's genuine within-microglia progression
  synergy).

**Verified post-feed:** ledger **89→91** rows; `specificity` layer = 2; contests
**amyloid 18 / synaptic 12 / interaction 55** ALL UNCHANGED (rebuilt
`biological_model_contest_verdicts.tsv`); themes T-Inflammation **31** (6 Strong)
> T-Synergy 30 (9 Strong) > T-Compartment 17 > T-Tau-attenuates 12; Hyp-0 7.

**Prose sync (hand-edited, hardcoded in the Rmd):** `rmd/16_biological_model.Rmd`
(7 edits: 89→91 counts, `lengthMenu` 91, axis-1 net 34−16, Hyp-1A 25 supports,
T-Inflammation 28/31 anchoring, cross-entity ranking flip with the "rests on one
Suggestive row / T-Synergy keeps more Strong" caveat, a new Phase N paragraph) and
`analysis.Rmd` (91 rows + Phase N clause; theme ranking flip with the same
caveat). Swept clean for residual "89"/old theme values (HTML: 0 stale "89", 5
"91-row/atomic", margins 18/12/55 each present, N-001 ×8 / N-002 ×4).

**Knit gate:** full re-render via `scripts/knit_analysis.sh` → **0 `class="error"`
/ 0 `class="warning"`**; §17 reflects 91 rows + 18/12/55 + the theme flip.

**Notes refreshed:** `map.md` (rmd/23 pipeline row, §24 anchor, cache row, 91-row
+ Phase N ledger paragraph) and `completed/DIGEST.md` (arc-N entry + top narrative
arc extended D→…→M→N + "no active plan as of 2026-06-06"). Plan archived to
`completed/celltype_specificity_plan_2026-06-06.md`. Arc N CLOSED.

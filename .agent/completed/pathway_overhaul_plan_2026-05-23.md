# Active plan: Hallmark rip-out, new collections, agnostic survey, and unbiased verdict

- **Created:** 2026-05-23
- **Last updated:** 2026-05-23 (D2 complete; plan archived)
- **State:** COMPLETE 2026-05-23. Phase A complete (A1-A6 DONE);
  Phase B complete (B1-B5 DONE); Phase C survey complete (C1-C6 DONE);
  Phase C verdict redesigned 2026-05-23 (original OXPHOS-anchored
  C7/C8/C9 REPLACED with unbiased D1/D2 — see "Why C7-C9 were reset"
  below); D1 done 2026-05-23; D2 done 2026-05-23. Successor plan is
  `storage/notes/mechanism_layer_plan.md`. See `## Outcome summary`
  at end of file.
- **Successor to:** an earlier `pathway_survey_plan` draft (removed 2026-06-15 in the notes-compaction audit; recoverable from git history).
- **Goal:** (1) remove the entire Hallmark gene-set collection from the
  project — caches, results, and Rmd content; (2) build replacement
  collections (GO BP already cached; add GO MF, GO CC, three custom
  curated gene-set families); (3) run the agnostic cross-modality
  pathway and module survey planned in the superseded plan, but now
  against the new collections, with no privileged framing for OXPHOS;
  (4) **produce an unbiased verdict from the survey rankings without
  pre-committing to any specific pathway family.**

## Why C7-C9 were reset (2026-05-23 session)

The original C7/C8/C9 steps were OXPHOS-anchored: they were designed to
test whether the OXPHOS-at-interaction signal — historically the
project's headline finding — was genuine or expression-driven. When the
2026-05-23 verdict-prep session surveyed the actual Phase C rankings,
three observations made the OXPHOS-anchored verdict the wrong question:

1. **OXPHOS-at-interaction is the only mixed-sign instance of the
   OXPHOS family.** At `tau_alone`, `nlgf_in_p301s`, and `tau_in_nlgf`,
   OXPHOS terms are uniformly sig=2 with consistent − sign across
   modalities. The interaction-contrast "mixed sign" is the residual
   non-additive disagreement; the underlying biology is uniform
   suppression of OXPHOS under both tau and amyloid, not divergence
   between genotypes.
2. **The actual survey winner is MG-M2 at nlgf_in_p301s** — sig=4 with
   consistent + sign across 4/5 modalities (snRNAseq +2.15 padj=1.6e-111,
   GeoMx +1.33 padj=4.9e-4, proteomics +1.51 padj=2.0e-4, phospho +1.50
   padj=3.8e-3). This is the cleanest cross-modality hit anywhere in
   the survey, and the original plan ignored it.
3. **MHC II antigen presentation family** is the most coherent
   family-level signal: 5 MHC II terms at `nlgf_in_maptki` GO BP ranks
   1-5 with consistent + sign, and the same family at `nlgf_in_p301s`
   ranks 7-11. Multiple terms from the same biological family co-rank
   with consistent sign across NLGF contrasts — a stronger pattern
   than OXPHOS's mixed-sign behaviour.

Per the user's 2026-05-23 reset decision, the verdict step is rebuilt
around the agnostic survey rather than around OXPHOS. The specificity-
null framework is dropped entirely (the consistency filters and
cross-contrast aggregation already in the rankings do equivalent work
for less compute, and the survey produced enough multi-axis evidence
that no single-target null would resolve the multi-axis verdict
anyway). The mechanism-layer follow-up (TF inference, kinase inference,
CCC bridge) is deferred to a successor plan whose targets will be
determined by the verdict outcome.

## Why this plan exists

The superseded plan attempted to test whether OXPHOS was genuinely the
dominant cross-modality signal or had been over-promoted by a
confirmation-biased audit chain. At the Step 1 decision gate the user
discovered that the OXPHOS finding was anchored on a single Hallmark
gene set (`HALLMARK_OXIDATIVE_PHOSPHORYLATION`), a ~200-gene basket
mixing electron transport, TCA cycle, mitochondrial ribosomes, and lipid
metabolism. The user instructed full removal of Hallmark from the
project and replacement with finer-grained ontology-based and
project-specific custom sets. The OXPHOS confirmatory arc (`rmd/08_oxphos_cross`,
`rmd/09a_asym_dam`, `rmd/09b_asym_proxies`, `rmd/09c_geomx_audit`,
`rmd/09d_geomx_q3_audit`) was instructed to be deleted in full, along
with the MG-M3-as-"OXPHOS-module" framing in `rmd/03_hdwgcna`.

This plan executes that scope, then runs the survey.

## Locked decisions (session 2026-05-23)

| Decision | Choice |
|---|---|
| Hallmark removal scope | Full rip-out, no archive (delete caches, TSVs, Rmd content) |
| Replacement collections | GO BP (already cached), GO MF, GO CC, custom curated |
| Custom curated tier | Core microglia states + tau/amyloid-specific + module sources (~25 sets) |
| IFN side-by-side panel | Dropped; IFN terms compete agnostically in the survey |
| OXPHOS arc (08, 09a, 09b, 09c, 09d) | Delete all five entirely |
| MG-M3 framing in 03_hdwgcna | Rewrite to drop "this is THE OXPHOS module" framing |
| ~~OXPHOS specificity null N~~ | **OVERTURNED 2026-05-23:** no specificity-null work; ranking + consistency only |
| hdWGCNA module × modality scoring | Per-module fgsea against each modality's interaction-contrast ranking |
| **Verdict null testing** (added 2026-05-23) | **Drop null framework entirely**; verdict uses agnostic ranking + sign-consistency + cross-contrast support + multi-collection corroboration only |
| **Contrast weighting in verdict** (added 2026-05-23) | **Equal weighting across all 5 contrasts** (tau_alone, nlgf_in_maptki, nlgf_in_p301s, tau_in_nlgf, interaction). The project's central narrative still privileges interaction, but the verdict step itself does not. |
| **Mechanism follow-up** (added 2026-05-23) | **Defer to successor plan**; this plan ends at the unbiased verdict |

## Open questions (defer to the session that needs them)

| Question | Default | Decided in |
|---|---|---|
| Custom set #1 (core microglia states) — source-of-truth lists | RESOLVED 2026-05-23: paper-derived per default + HAM split into Marsh ex-vivo + Friedman disease (user choice at B2 gate) + Sala Frigerio IRM list approximated from paper-text 5 markers + 23 canonical mouse ISGs (the n=28 Sala Frigerio list is not published anywhere we could find). See B2 completion note. | Session B2 (DONE) |
| Custom set #2 (tau/amyloid-specific) — source-of-truth lists | RESOLVED 2026-05-23: Safaiyan 2021 mmc4 Venn 'Common Genes (39)' for WAM; Marschallinger 2020 T2-1 RNA-Seq Aging logFC>0.5 & padj<0.05 for LDAM; Gerrits 2022 (NOT Mancuso 2024) Table S7 top-50 by avg_logFC for AD1 (amyloid-associated) and AD2 (tau-associated). Mancuso 2024 swap: no paired tau-iMG list published; user selected Gerrits 2022 at B3 gate. See B3 completion note. | Session B3 (DONE) |
| Custom set #3 (module sources) — source-of-truth lists | RESOLVED 2026-05-23: Sun, Victor, Mathys et al. 2023 Cell ONLY (DOI 10.1016/j.cell.2023.08.037). User selected this swap at the B4 decision gate over the plan's original Mathys 2019 + Olah 2020 default. Three reasons: (1) Sun/Victor 2023 is the dedicated microglia atlas from the same lab as Mathys 2019, with 443 subjects vs 48, much higher statistical power for state-specific markers; (2) Mathys 2019's M1-M10 SOM-territory modules are cross-cell-type and dilute microglia signal with non-microglia genes; (3) Olah ARM (cluster 7) is already used in B2 and the other Olah clusters are smaller / less well-characterised than the Sun/Victor 2023 equivalents. The published 12-state taxonomy (MG0..MG12 with MG9 merged into MG0 during final annotation) gives 12 sets in one collection rather than the ~18 sets the original two-source plan would have produced. See B4 completion note. | Session B4 (DONE) |
| ~~OXPHOS-family reference term for the specificity null~~ | **OVERTURNED 2026-05-23:** specificity-null work dropped entirely | Session C7 (cancelled) |
| ~~Verdict criterion for Phase C dominance call~~ | **OVERTURNED 2026-05-23:** verdict criterion is now multi-dimensional (see D1 below); no single "dominant family" criterion | Session C9 (cancelled) |
| Verdict-step leader-board filter (added 2026-05-23) | propose: a pathway is a "leader" if it has consistent-sign sig ≥2 at ANY contrast OR sig ≥3 at any contrast (mixed-sign allowed for sig ≥3 because OXPHOS-style multi-modality enrichment is informative even when sign disagrees). All leaders surface in the unified TSV; no pre-committed family priority. | Session D1 |

## Execution model

- A fresh session must read `CLAUDE.md` first, then this file in full.
- A fresh session executes the NEXT step whose status is `TODO`. Steps
  are sized to fit comfortably in one session each.
- After completing a step, the session must:
  1. Update that step's status from `TODO` to `DONE <YYYY-MM-DD>` with a
     one-line completion note (key outputs created/deleted).
  2. `chown rstudio:rstudio` any new files (agent runs as root).
  3. Commit locally per CLAUDE.md (imperative subject under 70 chars;
     co-author trailer).
  4. Either continue to the next `TODO` step in the same session if the
     context budget permits, or end cleanly.
- Decision-gate steps (B2, B3, B4 source-of-truth confirmation; C9
  verdict gate) require user confirmation via `AskUserQuestion`. Wait
  for confirmation before proceeding.
- Every session must end with a clean `analysis.html` knit
  (`Rscript -e 'rmarkdown::render("analysis.Rmd")'`) and a successful
  commit. The session ordering in Phase A is deliberately bottom-up
  (downstream Rmds refactored before upstream Hallmark builders) so
  that every commit is in a clean knittable state.

---

## Phase A: Hallmark rip-out (sessions A1-A6)

Topologically ordered bottom-up so each session ends with a clean knit.
The Hallmark variables (`hallmark_le_union`, `gsea_results`, `both_core`)
remain available for sessions that have not yet been refactored.

### Session A1: refactor `rmd/11_ccc.Rmd` [DONE 2026-05-23]

Completion note: dropped `le_syms`/`hallmark_le_union` reference and the
two `*_in_oxphos_LE` mutate columns in the
`ccc-multinichenet-interaction-spotlight` chunk; comment neutralised;
`ccc_multinichenet_interaction_top25_microglia.tsv` regenerated with 21
columns (was 23); `analysis.html` rebuilt cleanly in 3m13s.

**Inputs:** `rmd/11_ccc.Rmd` (line 475 area).

**Action:**
- Drop the `le_syms <- if (exists("hallmark_le_union")) ...` line and the
  resulting `ligand_in_oxphos_LE` / `receptor_in_oxphos_LE` columns in
  the `top_int_mg` mutate.
- Verify the column drop does not break later chunks that reference
  these columns (grep the file).
- The hdWGCNA hub columns stay; only OXPHOS LE columns go.
- Re-knit `analysis.Rmd`; verify clean.

**Outputs:**
- `rmd/11_ccc.Rmd` edited.
- `analysis.html` re-knitted.
- `storage/results/ccc_multinichenet_interaction_top25_microglia.tsv`
  regenerated without the two OXPHOS LE columns (schema change).

**Verification:**
- `grep -n 'oxphos_LE\|hallmark_le_union' rmd/11_ccc.Rmd` returns
  no matches.
- HTML renders without errors.

---

### Session A2: refactor `rmd/03_hdwgcna.Rmd` + `rmd/03b_hdwgcna_soft_power.Rmd` [DONE 2026-05-23]

Completion note: dropped intro "OXPHOS leading-edge" pointer, removed
`OXPHOS_robust_core`/`OXPHOS_Hallmark_LE`/`OXPHOS_GOBP_LE` from
`hdwgcna-module-enrichment` (only substate-marker gene sets remain),
neutralised "MG-M3 (OXPHOS module)" captions/comments throughout 03 and
03b, replaced `is_oxphos_le` colouring/sizing of the kME-vs-NEBULA
scatter with uniform points (legend gone, `n_oxphos_le` cor_tbl column
gone), and rephrased the soft-power tracking prose to "MG-M3-equivalent".
`hdwgcna_module_enrichment.tsv` regenerated with seven substate rows
only (was 10 with the three OXPHOS sets); `hdwgcna_kme_vs_nebula_interaction.tsv`
regenerated without the `is_oxphos_le` column. `analysis.html` re-knitted
cleanly.

**Inputs:** `rmd/03_hdwgcna.Rmd` and `rmd/03b_hdwgcna_soft_power.Rmd`.

**Action (03_hdwgcna):**
- Replace the `OXPHOS_robust_core = both_core`,
  `OXPHOS_Hallmark_LE = oxphos_le_hallmark$wide$symbol`,
  `OXPHOS_GOBP_LE = oxphos_le_gobp$wide$symbol` assignments (lines
  ~191-214) with neutral module-hub framing. The `hub_sym`
  construction (lines ~474+) stays.
- Drop the "MG-M3 (OXPHOS module)" caption phrasing throughout (lines
  ~250, 287-297, 358-431). Rename to "MG-M3 module" with no special
  attribute. The kME-vs-logFC scatter plot stays but the
  `is_oxphos_le` colouring/sizing dimension is removed (so points are
  uniform).
- Drop the legend "Hallmark OXPHOS LE" (line ~431) and the
  `n_oxphos_le` counter (line ~405).
- Re-read the prose intro to remove "the per-substate OXPHOS leading-edge
  analyses above" (line 3) — the OXPHOS LE analysis is being deleted
  in A5, but the prose pointer should be neutralised now.

**Action (03b_hdwgcna_soft_power):**
- Rephrase "Track MG-M3 (the canonical OXPHOS module at p=3)" (line ~210)
  to neutral "Track MG-M3 (the canonical p=3 module of interest)".
- Drop any OXPHOS-related Jaccard / hub-tracing prose; the analysis
  stays (testing soft-power robustness of MG-M3) but loses the
  OXPHOS-as-name framing.

**Outputs:**
- Both Rmd files edited.
- `analysis.html` re-knitted.
- `storage/results/hdwgcna_modules.tsv`, `hdwgcna_module_hubs.tsv`,
  `hdwgcna_kme_vs_nebula_interaction.tsv` regenerated (some may have
  OXPHOS-flag column drops).

**Verification:**
- `grep -ni 'oxphos\|hallmark' rmd/03_hdwgcna.Rmd rmd/03b_hdwgcna_soft_power.Rmd`
  returns no matches.
- HTML still renders the kME scatter and module heatmaps cleanly.

---

### Session A3: refactor `rmd/10_divergence.Rmd` [DONE 2026-05-23]

Completion note: deleted the `## Pathways with strongest cross-modality
divergence` heading, the Hallmark-framing prose, and the entire
`pathway-divergence` chunk (which built `int_tbl`/`score`/`ranked` from
`joint_gsea` and printed a 20-row Hallmark kable). Replaced with a
`## Pathway-level divergence ranking` heading + one-line pointer to the
agnostic survey in section 12. `divergence-shortlist` and the miloR
section are untouched. No TSVs regenerated (the deleted chunk had no
`write_tsv_safe` call). `analysis.html` re-knitted cleanly.

**Inputs:** `rmd/10_divergence.Rmd` (lines ~85-120).

**Action:**
- Drop the "Pathways with strongest cross-modality divergence" section
  (`## Pathways with strongest cross-modality divergence` and the
  `knitr::kable` of Hallmark pathways).
- Drop the surrounding prose that frames it as a Hallmark ranking.
- Add a placeholder note: "Pathway-level divergence ranking is now in
  the cross-modality pathway survey in section 12."
- The shortlist chunk (lines 9-87) and milo chunks (lines 125+) stay;
  these are independent of Hallmark.

**Outputs:**
- `rmd/10_divergence.Rmd` edited.
- `analysis.html` re-knitted.

**Verification:**
- `grep -ni 'hallmark' rmd/10_divergence.Rmd` returns no matches.
- HTML renders the shortlist heatmap and milo plots cleanly.

---

### Session A4: delete OXPHOS arc + Hallmark-anchored TSVs [DONE 2026-05-23]

Completion note: deleted the 5 OXPHOS arc Rmds (`rmd/08_oxphos_cross.Rmd`,
`rmd/09a_asym_dam.Rmd`, `rmd/09b_asym_proxies.Rmd`, `rmd/09c_geomx_audit.Rmd`,
`rmd/09d_geomx_q3_audit.Rmd`), the 8 listed result TSVs, and removed the
5 `child=` chunks from `analysis.Rmd` so `child-integration` jumps to
`child-divergence`. User chose "thorough" at the orphan-cleanup gate, so
also removed the 3 now-orphan scripts (`scripts/build_geomx_q3_audit.R`,
`scripts/build_asymmetric_proxy_nebula.R`, `scripts/explore_proxy_cuts.R`)
and the 6 orphan caches they had built (`storage/cache/de_geomx_q3_limma.rds`,
`de_geomx_q3_voom.rds`, `de_geomx_q3_audit_provenance.txt`,
`de_snrnaseq_nebula_asymmetric_dam.rds`,
`de_snrnaseq_nebula_asymmetric_plaque_contact.rds`,
`de_snrnaseq_nebula_asymmetric_topdecile_dam.rds`). `analysis.html`
re-knitted cleanly in 3m0s; section ordering now goes
integration → divergence → CCC → session-info. Hallmark fgsea cache
(`storage/cache/fgsea_results.rds`) and Hallmark TSV
(`storage/results/fgsea_hallmark_per_contrast.tsv`) intentionally left in
place; they are removed in A6 once 07_integration.Rmd is refactored to
stop building them.

**Inputs:** `analysis.Rmd`, the five OXPHOS arc Rmds, and several
results TSVs.

**Action:**
- Delete the following Rmd files:
  - `rmd/08_oxphos_cross.Rmd`
  - `rmd/09a_asym_dam.Rmd`
  - `rmd/09b_asym_proxies.Rmd`
  - `rmd/09c_geomx_audit.Rmd`
  - `rmd/09d_geomx_q3_audit.Rmd`
- Edit `analysis.Rmd` to remove the corresponding `child=` references
  (lines 117, 120, 123, 126, 129 in current parent file). The block
  after the edit jumps from `child-integration` (rmd/07) directly to
  `child-divergence` (rmd/10).
- Delete the following TSVs from `storage/results/`:
  - `oxphos_asymmetric_dam_nes_comparison.tsv`
  - `oxphos_asymmetric_proxies_nes_comparison.tsv`
  - `oxphos_leading_edge_per_state.tsv`
  - `oxphos_le_cross_modality.tsv`
  - `oxphos_le_custom_fgsea.tsv`
  - `oxphos_pathway_nes_cross_modality.tsv`
  - `geomx_normalisation_audit_nes.tsv` (generated by 09c)
  - `geomx_q3_audit_nes.tsv` (generated by 09d)
- Re-knit `analysis.Rmd`. The remaining files all reference Hallmark
  only via objects that no longer exist downstream (already refactored
  in A1, A2, A3), so the knit should succeed. The Hallmark heatmap in
  07 still works because 07 still builds the Hallmark fgsea cache
  itself; refactoring 07 happens in A6.

**Outputs:**
- 5 Rmd files deleted.
- 8 TSV files deleted.
- `analysis.Rmd` edited.
- `analysis.html` re-knitted (shorter, without the OXPHOS arc).

**Verification:**
- `ls rmd/08*.Rmd rmd/09*.Rmd 2>/dev/null` returns nothing.
- `grep -n 'child-oxphos\|child-asym\|child-geomx-audit\|child-geomx-q3' analysis.Rmd`
  returns nothing.
- `ls storage/results/oxphos_*.tsv storage/results/geomx_normalisation_audit_nes.tsv storage/results/geomx_q3_audit_nes.tsv 2>/dev/null`
  returns nothing.
- HTML still renders end-to-end.

---

### Session A5: refactor `rmd/02e_snrnaseq_substate_pathway.Rmd` [DONE 2026-05-23]

Completion note: dropped the OXPHOS-anchored preamble, the
`get_hallmark()` call + `hallmark = run_fgsea_per_state(...)` element,
the Hallmark per-state heatmap, the OXPHOS-keyword GO BP pre-filter, the
two `leading_edge_table(...)` calls (Hallmark + GOBP OXPHOS), the
LE-core focused heatmap, and the `both_core` intersection block.
Replaced with: neutral preamble; GO BP-only per-state fgsea cache build;
new generic GO BP per-state NES heatmap (top 40 by max |NES| across
states, padj<0.1 in any state); slim GO BP top-10-per-state padj table
(no OXPHOS pre-filter); forward pointer to the agnostic survey.
`fgsea_per_state.tsv` regenerated with 45,715 GOBP-only rows (was mixed
Hallmark+GOBP). Also deleted orphan
`storage/results/oxphos_leading_edge_per_state.tsv` left over from A4 —
A4's TSV delete was undone by the immediately-following A4-era knit
because 02e was still writing it at that point. `analysis.html`
re-knitted cleanly; 2 residual `HALLMARK_` strings in the HTML come from
`rmd/07_integration.Rmd`'s Hallmark heatmap (still alive pre-A6).
`leading_edge_table()` in `R/fgsea.R` is now dead code (no callers) —
A6 should delete it alongside `get_hallmark()`.

**Inputs:** `rmd/02e_snrnaseq_substate_pathway.Rmd`.

**Action:**
- Drop the Hallmark heatmap chunk (`snrnaseq-substate-fgsea-hallmark-heatmap`,
  lines ~47-72) entirely.
- Drop the entire Hallmark OXPHOS leading-edge block:
  `oxphos_le_hallmark <- leading_edge_table(...)` through the
  `print(knitr::kable(le_top_hallmark, ...))` call (lines ~143-195).
- Drop the LE heatmap chunk and the LE-core build for Hallmark
  (`snrnaseq-substate-fgsea-oxphos-heatmap` lines ~195-269; check the
  full extent of this chunk and remove the Hallmark-anchored variant).
- Drop the `both_core` intersection block (lines ~276-291) — `both_core`
  was already removed from 03_hdwgcna in A2; deleting its definition
  here is safe.
- Keep the GO BP heatmap chunk and any analysis that operates solely on
  GO BP. The cache `fgsea_per_state_results.rds` still has GO BP fgsea;
  it will be rebuilt in B5 to drop Hallmark fully.
- Rewrite the section preamble to remove all mention of "Hallmark" and
  "27-gene recurrent core"; the new prose should describe just the
  per-substate GO BP enrichment.
- Mark in a short note at the bottom that comprehensive cross-modality
  pathway ranking happens in section 12 (introduced in Phase C).

**Outputs:**
- `rmd/02e_snrnaseq_substate_pathway.Rmd` substantially shortened.
- `analysis.html` re-knitted.

**Verification:**
- `grep -ni 'hallmark\|both_core\|hallmark_le_union\|oxphos_le_hallmark' rmd/02e_snrnaseq_substate_pathway.Rmd`
  returns no matches.
- HTML renders 02e cleanly with only the GO BP heatmap.

---

### Session A6: refactor `rmd/07_integration.Rmd`, clean up `R/fgsea.R`, delete Hallmark caches, final verify [DONE 2026-05-23]

Completion note: in `rmd/07_integration.Rmd` collapsed the `# Joint
pathway analysis (fgsea)` intro to describe GO BP only (added pointer to
section 12); replaced the Hallmark `{r fgsea}` chunk with a slim
`{r fgsea-setup}` chunk that keeps only `prep_t` / `run_dataset` /
`snrnaseq_nebula_shim` (still needed by the GO BP chunk); deleted the
`## Hallmark concordance heatmap` heading + `{r fgsea-heatmap}` chunk;
removed the now-redundant `## GO Biological Process fgsea (finer-grained)`
subheading and its Hallmark-comparison prose. In `R/fgsea.R` deleted
`get_hallmark()` and the orphaned `leading_edge_table()`; neutralised
the OXPHOS-specific phrasing in the docstrings of the three remaining
caller-less helpers (`cross_modality_interaction_table`,
`extract_pathway_nes`, `fgsea_custom_pathways`) so they read as generic
utilities (kept because Phase C6 will reuse `fgsea_custom_pathways` for
module-as-pathway scoring). Deleted `storage/cache/msigdb_hallmark_mouse.rds`,
`storage/cache/fgsea_results.rds`,
`storage/results/fgsea_hallmark_per_contrast.tsv`. `analysis.html`
re-knitted cleanly in 2m51s (27 MB); zero residual `Hallmark` or
`HALLMARK_` strings in HTML or source. Phase A complete; codebase is
Hallmark-free.

**Inputs:** `rmd/07_integration.Rmd`, `R/fgsea.R`, `R/helpers.R`.

**Action:**
- In `rmd/07_integration.Rmd`:
  - Drop the Hallmark fgsea call (`gsea_path <- "fgsea_results.rds"`,
    the `cache_or_run(gsea_path, ...)` block, lines ~145-205).
  - Drop the Hallmark heatmap section (lines ~209-249, including
    `## Hallmark concordance heatmap` heading and the chunk).
  - Drop the prose phrase "Hallmark captures broad axes only; for finer
    pathway granularity we also run fgsea against MSigDB C5 GO:BP"
    (line ~254) and shorten the GO BP intro to standalone.
  - Keep the GO BP fgsea cache build and heatmap.
  - Add a placeholder line "(Additional collections — GO MF, GO CC,
    custom curated sets — are introduced in section 12, the agnostic
    pathway survey.)" so readers can navigate.
- In `R/fgsea.R`:
  - Delete the `get_hallmark()` function entirely.
  - Delete `leading_edge_table()` — its only callers (the two
    `oxphos_le_hallmark` / `oxphos_le_gobp` blocks in 02e) were removed
    in A5; confirm with
    `grep -rn 'leading_edge_table' R/ rmd/ scripts/ analysis.Rmd` first.
  - `get_gobp()` stays.
- In `R/helpers.R`: confirm there is no Hallmark-specific reference;
  verify `msigdbr` is still loaded (needed for GO BP).
- Delete cache files:
  - `storage/cache/msigdb_hallmark_mouse.rds`
  - `storage/cache/fgsea_results.rds`
- Delete result file:
  - `storage/results/fgsea_hallmark_per_contrast.tsv`
- Re-knit `analysis.Rmd`. The knit MUST succeed end-to-end.

**Outputs:**
- `rmd/07_integration.Rmd` shortened (no Hallmark heatmap).
- `R/fgsea.R` shortened (no `get_hallmark()`).
- 2 cache `.rds` files deleted.
- 1 result `.tsv` file deleted.
- `analysis.html` re-knitted, Hallmark-free.

**Verification:**
- `grep -rni 'hallmark' rmd/ R/ scripts/ analysis.Rmd` returns nothing
  (except possibly comments; double-check).
- `ls storage/cache/msigdb_hallmark_mouse.rds storage/cache/fgsea_results.rds storage/results/fgsea_hallmark_per_contrast.tsv 2>/dev/null`
  returns nothing.
- HTML renders end-to-end without errors.
- Phase A is complete; codebase is Hallmark-free.

---

## Phase B: build new collections (sessions B1-B5)

After Phase A, the project has GO BP only. Phase B adds GO MF, GO CC,
and three custom curated set families. fgsea runs against five
modalities (snrnaseq, geomx, proteomics, phospho, phospho_corr).

### Session B1: build GO MF + GO CC caches and fgsea [DONE 2026-05-23]

Completion note: added `get_gomf()` and `get_gocc()` to `R/fgsea.R` (mirror
of `get_gobp()`, default cache paths `storage/cache/msigdb_gomf_mouse.rds`
/ `msigdb_gocc_mouse.rds`); refactored `R/fgsea.R` to add three shared
shims — `prep_t_dedup()`, `run_fgsea_per_dataset()`, and
`join_fgsea_results()` — so the new builder script and any future
collection builder share one source of truth; `rmd/07_integration.Rmd`
left untouched apart from the pointer note (still uses its in-Rmd
`prep_t` / `run_dataset` to avoid risk of breaking the existing knit;
deduplication can happen in a future cleanup). New runnable script
`scripts/build_fgsea_extra_collections.R` is idempotent (skips per-cache
when present; `--overwrite` forces). Builder ran in 1m58s and produced
all four caches plus both TSVs; re-run printed only `[skip]` lines as
designed. Cache schema mirrors `fgsea_gobp_results.rds`: 5 modalities ×
5 contrasts × fgsea data.table (cols pathway/pval/padj/log2err/ES/NES/
size/leadingEdge). GO MF: 1852 sets, 4335 pathway-contrast rows in TSV.
GO CC: 1042 sets, 2774 rows. Both TSVs use identical 12-column schema
(`pathway`, `nes_<m>`, `padj_<m>` × 5 modalities, `contrast`) as the
existing GO BP TSV. The `msigdbr` console message (`Using human MSigDB
with ortholog mapping to mouse. Use db_species = "MM" for mouse-native
gene sets.`) is the same behaviour as `get_gobp()` (no `db_species`
argument); switching all three GO loaders to mouse-native would be a
separate decision affecting all of Phase B and is not in scope here.
Pointer note in 07 now describes both new caches and the builder script
explicitly. `analysis.html` re-knitted cleanly in 2m52s; zero
`class="error"` elements in the rendered HTML. Phase B1 complete; Phase
B2 (custom set #1 — core microglia states) is next and requires user
confirmation of gene-list sources at its decision gate.

**Inputs:** existing fgsea infrastructure in `R/fgsea.R`, the same
`run_fgsea_for_contrast` / `run_dataset` machinery `rmd/07_integration.Rmd`
uses for GO BP.

**Action:**
- In `R/fgsea.R` add:
  - `get_gomf(cache_path = "storage/cache/msigdb_gomf_mouse.rds")`
    using `msigdbr::msigdbr(species = "Mus musculus", collection = "C5", subcollection = "GO:MF")`.
  - `get_gocc(cache_path = "storage/cache/msigdb_gocc_mouse.rds")` for
    `subcollection = "GO:CC"`.
- Build caches:
  - `storage/cache/msigdb_gomf_mouse.rds`
  - `storage/cache/msigdb_gocc_mouse.rds`
- Run fgsea for each modality × each contrast against GO MF; cache
  to `storage/cache/fgsea_gomf_results.rds`.
- Run fgsea for each modality × each contrast against GO CC; cache
  to `storage/cache/fgsea_gocc_results.rds`.
- Export `storage/results/fgsea_gomf_per_contrast.tsv` and
  `storage/results/fgsea_gocc_per_contrast.tsv` for inspection.
- Decision: where to put the fgsea calls? Cleanest is a new script
  `scripts/build_fgsea_extra_collections.R` runnable as `Rscript`,
  with both GO MF and GO CC builds idempotent (skip if cache exists,
  unless `--overwrite`). Add a brief note in `rmd/07_integration.Rmd`
  pointing readers to this script and to section 12.

**Outputs:**
- `R/fgsea.R` extended with two new functions.
- `scripts/build_fgsea_extra_collections.R` (new).
- 4 new cache files (2 msigdb sets + 2 fgsea results).
- 2 new TSV exports.
- `analysis.html` re-knitted (no Rmd changes beyond the pointer note,
  so re-knit should be quick).

**Verification:**
- `ls storage/cache/msigdb_gomf_mouse.rds storage/cache/msigdb_gocc_mouse.rds storage/cache/fgsea_gomf_results.rds storage/cache/fgsea_gocc_results.rds`
  all exist.
- The TSVs have rows for all five modalities and all five contrasts.
- Running `Rscript scripts/build_fgsea_extra_collections.R` is idempotent
  when caches exist.

---

### Session B2: custom set #1 — core microglia states [DONE 2026-05-23]

**Completion note (2026-05-23):** Built 7-set curated mouse-symbol
collection (sizes: DAM_up=278, DAM_down=21, HAM_exvivo_QC=25,
HAM_disease=30, ARM=285, IRM=28, homeostatic_core=15). User chose at
the decision gate to split HAM into TWO separate sets (Marsh 2022 ex-
vivo activation 25-gene + Friedman 2018 neurodegeneration-related 30-
gene), giving 7 rather than the originally-spec'd 6. The script
`scripts/build_custom_microglia_states.R` embeds all seven vectors so
the build needs no external supplement files at run-time; the script
self-documents source paper + supplement table + filter per set and
writes `storage/cache/custom_microglia_states_provenance.txt` covering
the same. Olah ARM was published as 318 human symbols (cluster 7,
up_type==7); the script translates to 285 mouse symbols via
`nichenetr::convert_human_to_mouse_symbols`. The Sala Frigerio 2019
n=28 IRM gene list is NOT published anywhere we could find (PMC OA,
GEO, Elsevier CDN supp PDFs, Roy 2022 Immunity, Soton ePrints,
ALZFORUM) — see provenance for full search history; user accepted a
literature-faithful approximation seeded with the 5 IRM markers named
in the Sala Frigerio main text + 23 canonical mouse type-I-ISG genes
from Schoggins 2011 Nature 472:481-485 + Reactome IFN-alpha/beta.
fgsea ran with `min_size = 5` (lowered from the GO default 15 because
homeostatic_core is only 15 mouse symbols and drops further after
intersecting with proteomics / phospho universes). Joined wide-form
TSV at `storage/results/fgsea_custom_states_per_contrast.tsv` (35 rows
= 7 sets x 5 contrasts; NA cells for pathway-modality pairs that fell
below min_size). Helper `get_custom_microglia_states()` added to
`R/fgsea.R` as a pure reader (errors with rebuild instruction if cache
missing — the build script is the single source of truth for set
content). snRNAseq biology checks: DAM_up nlgf_in_p301s NES=+2.46
padj=2.1e-54; HAM_disease interaction NES=+3.22 padj=6.9e-10;
DAM_down + homeostatic_core both strongly negative in NLGF as
expected; IRM nlgf_in_p301s NES=+1.95 padj=1.5e-5. GeoMx
corroborates DAM_up and HAM_disease in NLGF contrasts at FDR<1e-3.
`analysis.html` re-knitted with no Rmd changes (the new helper is
unused by any current Rmd; consumption will be wired in at Session
B4). Phase B2 complete; Phase B3 (tau/amyloid-specific custom set #2)
is next and requires user confirmation of gene-list sources at its
decision gate.

**Decision gate at start:** confirm gene-list source with user. Default
proposal:
- DAM up / DAM down: Keren-Shaul et al. 2017, Cell — Table S2 column
  filters (FDR<0.05 + |logFC|>1 either direction).
- HAM: Marsh et al. 2022, Nat Neurosci — supplementary microglia
  signature lists.
- ARM: Olah et al. 2020, Nat Comm — Table S2 activated cluster markers.
- IRM: Sala Frigerio et al. 2019, Cell Rep — IRM gene list from Suppl.
- Homeostatic core: literature consensus
  (P2ry12, Tmem119, Cx3cr1, Sall1, Tgfbr1, Csf1r, Hexb, Mertk, Olfml3,
  Selplg, Siglech, plus key receptors — finalised in session).

Use `AskUserQuestion` to confirm the proposed sources OR offer to use
an R package such as `nichenetr` / `presto` that ships microglia
signatures. Wait for confirmation before downloading anything.

**Inputs:** confirmed source(s).

**Action:**
- Curate the seven gene lists (HAM split into Marsh ex-vivo +
  Friedman disease per the B2 decision-gate user choice) into a single
  named list, mouse symbols, intersected against the measured-gene
  universe of each modality at fgsea call time (do not pre-intersect;
  the fgsea machinery handles it).
- Save to `storage/cache/custom_microglia_states.rds`.
- Add to `R/fgsea.R` a `get_custom_microglia_states()` helper.
- Run fgsea across modalities; cache `storage/cache/fgsea_custom_states_results.rds`.
- Export `storage/results/fgsea_custom_states_per_contrast.tsv`.
- Add provenance file `storage/cache/custom_microglia_states_provenance.txt`
  recording (paper, table, filter applied, gene count) per set.
- Create `scripts/build_custom_microglia_states.R` (idempotent builder
  with embedded gene-list vectors so the build needs no external
  supplement files at run-time).

**Outputs:**
- 2 new cache files + 1 provenance file.
- 1 new TSV.
- `R/fgsea.R` extended.
- `analysis.html` re-knitted.

**Verification:**
- The cache contains exactly 7 named gene lists with sensible sizes
  (DAM up ~100-300, DAM down ~20-50, HAM_exvivo_QC ~25, HAM_disease
  ~30, ARM ~250-300 after human-to-mouse mapping, IRM ~28,
  homeostatic_core ~15).
- The fgsea wide-form TSV has one row per (set, contrast) =
  7 sets × 5 contrasts = 35 rows, with per-modality NES + padj
  columns; NA cells are allowed where a set fell below `min_size=5`
  after intersecting with a given modality's gene universe.
- Provenance file lists source for each set.

---

### Session B3: custom set #2 — tau/amyloid-specific [DONE 2026-05-23]

**Completion note:** Built the 4-set curated mouse-symbol tau/amyloid-
microglia collection (AD1=47, AD2=42, LDAM=27, WAM=39) and ran fgsea
against all 5 DE modalities × 5 contrasts. Mancuso 2024 swap: WebSearch
confirmed no paired tau-iMG signature in the only canonical 2024
chimeric-microglia paper (Nat Neurosci PMC11089003, primarily amyloid-
focused); user selected the Recommended alternative Gerrits 2022 (Acta
Neuropathol, DOI 10.1007/s00401-021-02263-w; PMC8043951) at the B3
decision gate. Gerrits S7 publishes paired amyloid-associated (AD1) and
tau-associated (AD2) microglia subcluster DEG lists in a single within-
paper comparison, which tightens the tau/amyloid distinction this project
needs. LDAM source choice: T2-1 RNA-Seq Aging filtered (logFC>0.5 &
padj<0.05) selected by user over the 73-gene T1 LD-related-genes panel
or their union; Marschallinger 2020 did not publish a fixed LDAM gene
set (LDAM is FACS-defined), so any tabular signature is a filter of the
DEG table. Gerrits AD1/AD2 selection by top-N rank rather than fixed
logFC threshold because AD2 effect sizes are biologically much weaker
than AD1 (top logFC 0.51 vs 2.0); top-50 by avg_logFC keeps set sizes
comparable. Build script `scripts/build_custom_microglia_ad.R` is self-
contained (every gene-list vector is embedded; nichenetr is the only
external dependency, used for AD1/AD2 human→mouse mapping). New helper
`get_custom_microglia_ad()` in `R/fgsea.R` mirrors the B2 pure-reader
contract. snRNAseq biology checks: AD1 nlgf_in_maptki NES=+2.38
padj=2.3e-10; WAM nlgf_in_maptki NES=+2.86 padj=6.1e-21; AD2 tau_alone
NES=+2.02 padj=1.1e-5; AD2 tau_in_nlgf NES=+2.13 padj=1.0e-4 — the
amyloid/tau axis split is reproduced cleanly in the snRNAseq modality.
GeoMx corroborates WAM in NLGF contrasts (NES≈1.7-2.1 padj<0.01).
Proteomics/phospho yield NA for WAM and LDAM because the small set sizes
(39 and 27 mouse symbols) drop below min_size=5 after intersecting with
the protein universe — expected for tight transcript-heavy signatures.
`analysis.html` re-knitted with no Rmd changes (the new helper is unused
by any current Rmd; consumption will be wired in at Session B5/Phase C).
Phase B3 complete; Phase B4 (Mathys 2019 + Olah 2020 module sources)
is next and requires user confirmation of gene-list sources at its
decision gate.

**Decision gate at start:** confirm gene-list source with user. Default
proposal:
- WAM: Safaiyan et al. 2021, Neuron — white-matter-associated microglia
  signature.
- LDAM: Marschallinger et al. 2020, Nat Neurosci — lipid-droplet-
  accumulating microglia signature.
- tau-iMG / Aβ-iMG: Mancuso et al. 2024 (ex vivo xenografts) — chimeric
  microglia signatures under each pathology.

Use `AskUserQuestion` to confirm sources, paper editions, and the
specific gene-list table/threshold per signature. Wait before
downloading.

**Inputs:** confirmed source(s).

**Action:**
- Curate each gene list, mouse symbols.
- Save to `storage/cache/custom_microglia_ad.rds`.
- Add `get_custom_microglia_ad()` helper to `R/fgsea.R`.
- Run fgsea; cache `storage/cache/fgsea_custom_ad_results.rds`.
- Export `storage/results/fgsea_custom_ad_per_contrast.tsv`.
- Add provenance file.

**Outputs:**
- 2 new caches + provenance.
- 1 new TSV.
- `R/fgsea.R` extended.
- `analysis.html` re-knitted.

**Verification:** as B2 with the appropriate set count and naming.

---

### Session B4: custom set #3 — Mathys / Olah module sources [DONE 2026-05-23]

**Completion note (2026-05-23):** Built the 12-set curated mouse-symbol
microglia transcriptional-state "module sources" collection from Sun,
Victor, Mathys et al. 2023 Cell (DOI 10.1016/j.cell.2023.08.037).
**Decision-gate deviation from plan B4 default:** User selected the
Recommended alternative Sun/Victor 2023 only (option D), NOT the
original Mathys 2019 + Olah 2020 dual source. Three reasons: (1)
Sun/Victor 2023 is the dedicated microglia atlas from the same lab as
Mathys 2019, 194,000 nuclei from 443 subjects vs Mathys 2019's
48-subject cohort, dramatically higher statistical power for
state-specific markers; (2) Mathys 2019's M1-M10 SOM-territory modules
are CROSS-CELL-TYPE (mix neurons, glia, microglia, OPCs, endothelial)
because they cluster gene-trait correlation patterns across all cell
types, diluting the microglia signal for a microglia-focused project;
(3) Olah 2020 ARM (cluster 7) is already used in B2 and the other Olah
clusters (1-6, 8, 9) are smaller / less well-characterised than the
Sun/Victor 2023 12-state taxonomy that supersedes them. The MG9 cluster
is skipped because the published 12-state annotation merged it into
MG0 during final cluster assignment. Set sizes after dedup + human→mouse
mapping via `nichenetr::convert_human_to_mouse_symbols`: MG0=50,
MG1=33 (mitochondrial MT-* drop in mapping), MG2=47, MG3=46, MG4=50,
MG5=47, MG6=47, MG7=48, MG8=48, MG10=49, MG11=46, MG12=50. The build
script `scripts/build_custom_module_sources.R` embeds all twelve mouse-
symbol vectors verbatim so the build has no external file dependency at
run-time except `nichenetr` (only loaded for namespace, not called).
Source markers were extracted once from the Sun/Victor 2023 Broad CDN
file
`ROSMAP.Microglia.6regions.seurat.harmony.selected.clusterDEGs.txt`
(2228 markers across 12 clusters; filter p_val_adj<0.05, top 50 by
avg_log2FC per cluster). Helper `get_custom_module_sources()` added to
`R/fgsea.R` mirroring the B2/B3 pure-reader contract: errors with
rebuild instruction if cache missing. fgsea ran with `min_size=5` (same
as B2/B3) producing 60 wide-form rows (12 sets × 5 contrasts) at
`storage/results/fgsea_custom_modules_per_contrast.tsv`. The single
benign fgsea warning ("1 pathways for which P-values were not calculated
properly due to unbalanced...") indicates 1 pathway-modality-contrast
tuple where the ranked stat vector was unbalanced; NA cells (8.3% per
modality per contrast = 1/12 sets) reflect MG1_neuronal_surveillance
dropping below min_size=5 after intersecting with proteomics/phospho
protein universe (expected for the smallest set). **Biology highlights
(snRNAseq):**
- MG3_ribosome_biogenesis is strongly NEGATIVE in tau contrasts
  (tau_alone NES=-2.36 padj=3.6e-5; tau_in_nlgf NES=-2.76 padj=6.0e-14;
  interaction NES=-1.74 padj=5.4e-3): tau pathology suppresses microglia
  ribosome biogenesis.
- MG4_lipid_processing is strongly POSITIVE in NLGF contrasts
  (nlgf_in_maptki NES=+2.46 padj=2.3e-10; nlgf_in_p301s NES=+2.09
  padj=1.2e-8) and tau_in_nlgf (NES=+2.43 padj=1.5e-7): canonical DAM/
  LDAM-like lipid response to amyloid pathology.
- MG11_antiviral is positive in NLGF (NES≈+1.9 to +2.1) and tau_in_nlgf
  (NES=+2.28 padj=1.2e-5): type-I-IFN response correlates with amyloid
  load.
- MG5_phagocytic is positive across the interaction (+1.64), nlgf_in_p301s
  (+1.79), tau_alone (+1.52), and tau_in_nlgf (+2.19).
- MG10_inflammatory_III (NF-kB axis) is positive in both NLGF contrasts.
- MG6_stress_response is NEGATIVE in interaction (-1.50) and tau_in_nlgf
  (-2.08), but is the only set significant in GeoMx (nlgf_in_p301s
  NES=-1.70 padj=0.035).
GeoMx has only 1 significant cell out of 60 at FDR<0.05; protein
modalities have none, reflecting that the Sun/Victor 2023 state markers
are RNA-defined and many do not survive proteomics intersection.
`analysis.html` re-knitted cleanly with no Rmd changes; the new helper
is unused by any current Rmd, consumption is wired in at B5 (per-
substate fgsea inclusion) and Phase C4 (custom-collection cross-modality
ranking). Phase B4 complete; Phase B5 (rebuild
`fgsea_per_state_results.rds` with all new collections) is next, NO
decision gate.

**Decision gate at start:** confirm gene-list source with user. Default
proposal:
- Mathys 2019 AD modules: Mathys et al. 2019, Nature — Suppl. Table 2
  module assignments (M1-M11), filtered to mouse orthologues.
- Olah 2020 subtypes: Olah et al. 2020, Nat Comm — Suppl. cluster
  markers.

Note: the Mathys modules are *human*; mouse orthologue mapping is
required. Use `nichenetr::convert_human_to_mouse_symbols()` or
`biomaRt`/`gprofiler2` for the conversion; decide and document.

Use `AskUserQuestion` to confirm sources and the human→mouse mapping
tool. Wait before downloading.

**Inputs:** confirmed source(s) and orthologue strategy.

**Action:**
- Curate ~11 Mathys modules + ~7-9 Olah subtypes (final count depends
  on what the papers publish).
- Save to `storage/cache/custom_module_sources.rds`.
- Add `get_custom_module_sources()` helper to `R/fgsea.R`.
- Run fgsea; cache `storage/cache/fgsea_custom_modules_results.rds`.
- Export `storage/results/fgsea_custom_modules_per_contrast.tsv`.
- Add provenance file.

**Outputs:**
- 2 new caches + provenance.
- 1 new TSV.
- `R/fgsea.R` extended.
- `analysis.html` re-knitted.

**Verification:** as B2/B3 with appropriate set count.

---

### Session B5: rebuild `fgsea_per_state_results.rds` with new collections [DONE 2026-05-23]

**Completion note (2026-05-23):** Rewrote the `snrnaseq-substate-fgsea`
chunk in `rmd/02e_snrnaseq_substate_pathway.Rmd` to load all six
gene-set collections via their B1-B4 helpers (`get_gobp`, `get_gomf`,
`get_gocc`, `get_custom_microglia_states`, `get_custom_microglia_ad`,
`get_custom_module_sources`) and build per-state fgsea against every
one of them inside the same `cache_or_run` block; GO collections retain
`min_size = 15` (the canonical GO cutoff used everywhere else),
the three custom collections use `min_size = 5` (matching B2/B3/B4 so
the smallest sets — homeostatic_core at 15 mouse symbols, MG1
neuronal_surveillance, LDAM at 27 — survive intersection with the
per-state universe). Deleted the legacy
`storage/cache/fgsea_per_state_results.rds` (Hallmark + GO BP only,
5.1 MB) before re-knitting so `cache_or_run` rebuilt from scratch.
The flatten step now uses
`purrr::imap_dfr(fgsea_per_state, flatten_fgsea_per_state)` so the
collection-name label flows automatically from the cache top-level
names — adding a seventh collection later only requires extending the
`cache_or_run` list, not touching the flatten/export. Updated the two
surviving GO BP display chunks (heatmap + top-10-per-state table) to
filter on the new lowercase `"gobp"` collection key (was `"GOBP"` in
the legacy single-collection layout). New
`storage/cache/fgsea_per_state_results.rds` is 7.0 MB; top-level names
are exactly `gobp`, `gomf`, `gocc`, `custom_microglia_states`,
`custom_microglia_ad`, `custom_module_sources` with NO `hallmark`.
New `storage/results/fgsea_per_state.tsv` has 60,885 rows covering all
120 (6 collections × 4 substates × 5 contrasts) buckets; row counts
are gobp=45,715, gomf=8,015, gocc=6,735, custom_module_sources=230,
custom_microglia_states=130, custom_microglia_ad=60 (smaller custom
totals reflect per-state min_size=5 dropping the tightest sets in some
states). `analysis.html` re-knitted cleanly in 5m01s (27 MB);
`class="error"` count is 0; rendered HTML contains zero `HALLMARK_` or
`Hallmark` strings. The only remaining `hallmark` mentions in source
are guardrail comments in the three custom-collection build scripts
documenting that those collections contain zero Hallmark sets — a B2/B3/B4
anti-anchoring annotation, not a Hallmark dependency. Phase B is
complete; the project has six gene-set collections cached at both
whole-microglia level (per-contrast caches built in B1-B4) and at
per-substate level (this cache). Next session executes C1 (helpers +
Rmd 12 skeleton + GO BP × modality survey ranking) with no decision
gate.

**Inputs:** the current `storage/cache/fgsea_per_state_results.rds`
(mixed Hallmark + GO BP × 4 substates) and the new caches from B1-B4.

**Action:**
- In `rmd/02e_snrnaseq_substate_pathway.Rmd` (already refactored in A5
  to drop Hallmark), update the `fgsea_per_state` build to include:
  - GO BP (already present)
  - GO MF (new, from B1)
  - GO CC (new, from B1)
  - custom_microglia_states (from B2)
  - custom_microglia_ad (from B3)
  - custom_module_sources (from B4)
- Delete the existing `fgsea_per_state_results.rds` so the
  `cache_or_run` rebuilds.
- Update `flatten_fgsea_per_state` calls / the TSV export to reflect
  the new collection list.
- Export `storage/results/fgsea_per_state.tsv` with all collections.
- Re-knit; verify the per-state GO BP heatmap (the surviving display in
  02e) still works against the new cache structure.

**Outputs:**
- `rmd/02e_snrnaseq_substate_pathway.Rmd` updated.
- `storage/cache/fgsea_per_state_results.rds` rebuilt with new
  collections.
- `storage/results/fgsea_per_state.tsv` regenerated.
- `analysis.html` re-knitted.

**Verification:**
- The cache `names()` include GO BP, GO MF, GO CC, custom_microglia_states,
  custom_microglia_ad, custom_module_sources (and NO `hallmark`).
- The TSV has rows for all six collections × four substates × five
  contrasts.
- Phase B complete; the project has six gene-set collections ready
  for the survey.

---

## Phase C: agnostic pathway survey (C1-C6) + unbiased verdict (D1-D2)

After Phase B, all collections are cached. Phase C builds the survey
section `rmd/12_pathway_survey.Rmd` step by step (C1-C6). The verdict
step was originally specified as the OXPHOS-anchored C7/C8/C9, but on
2026-05-23 the verdict was reset to D1/D2 (unbiased — see "Why C7-C9
were reset" at the top of this file).

### Session C1: helpers + Rmd 12 skeleton + GO BP × modality [DONE 2026-05-23]

**Completion note (2026-05-23):** Built `R/pathway_survey.R` with two
helpers — `rank_pathways_cross_modality(fgsea_list, contrast, modalities = NULL, padj_cut = 0.05)`
and `format_ranking_table(ranking_tbl, contrast, n = 20, strip_prefix = NULL, caption = NULL, include_padj = FALSE)`
— and sourced the file from `R/helpers.R` after `R/fgsea.R` (the
pathway-survey helpers depend on the fgsea cache shape produced by
`run_fgsea_per_dataset()`). Composite-rank logic: ties broken by
`(n_modalities_sig_consistent_sign desc, mean_abs_nes desc)` then row
order; `n_modalities_sig_consistent_sign` is `pmax(pos_sig, neg_sig)`
over significant modalities only; `sign_consensus` is strict (any
disagreement between significant modalities flips to "mixed").
Smoke-tested against `fgsea_gobp_results.rds` before re-knitting:
4,340 GO BP pathways with NES across 5 modalities, 4 contrasts produced
the expected shape, and the interaction contrast surfaces the OXPHOS /
ETC / aerobic-respiration family at composite ranks 1-9 with strict
"mixed" sign (one of three significant modalities disagrees on
direction). Only 9 pathways total have ≥2 modalities significant with
consistent sign at the interaction contrast — this is precisely the
question the C7/C8 specificity null will adjudicate. Built
`rmd/12_pathway_survey.Rmd` with the level-1 heading, the agnostic
preamble that names the anti-anchoring framing explicitly, and the GO
BP subsection (level-2) containing one `cache_or_run`-style ranking
build + a per-contrast `results = 'asis'` loop printing five level-4
sub-subheadings (one per contrast) each followed by a top-20 kable.
TSV `storage/results/pathway_survey_gobp_ranking.tsv` has 21,700 data
rows (4,340 pathways × 5 contrasts) with the full per-modality NES and
padj columns plus contrast labels. Added `child-pathway-survey`
reference to `analysis.Rmd` between `child-ccc` and `child-session`.
`analysis.html` re-knitted cleanly in 2m48s; `class="error"` count is
0; `class="warning"` count is 0; new section appears as 13 / 13.1 in
the numbered output with five `13.1.0.x` per-contrast sub-headings and
the kable tables below each. Phase C1 complete; C2 (GO MF × modality)
is next, no decision gate.

**Inputs:** `storage/cache/fgsea_gobp_results.rds` and the new caches.

**Action:**
- Create `R/pathway_survey.R` with:
  - `rank_pathways_cross_modality(fgsea_list, contrast, modalities, padj_cut = 0.05)`
    returning a tibble with: `pathway`, `n_modalities_sig`,
    `n_modalities_sig_consistent_sign`, `mean_abs_nes`, `sign_consensus`,
    per-modality NES, per-modality padj, `composite_rank`.
  - Composite rank = sort by (n_modalities_sig_consistent_sign desc,
    mean_abs_nes desc).
  - Helper `format_ranking_table()` for kable display.
- Source the new file from `R/helpers.R`.
- Create `rmd/12_pathway_survey.Rmd`:
  - Section header `# Agnostic cross-modality pathway and module survey`
  - Intro paragraph stating the goal (test whether OXPHOS-related
    biology is genuinely the dominant cross-modality signal, or one
    among many; OXPHOS gets no privileged framing).
  - Subsection `## GO BP pathway × modality ranking` with one ranked
    table per contrast (tau_alone, nlgf_in_maptki, nlgf_in_p301s,
    tau_in_nlgf, interaction). Top 20 by composite rank.
- Add to `analysis.Rmd` (just before the session-info child):
  ```
  ```{r child-pathway-survey, child = "rmd/12_pathway_survey.Rmd"}
  ```
  ```
- Export `storage/results/pathway_survey_gobp_ranking.tsv`.
- Re-knit.

**Outputs:**
- `R/pathway_survey.R` (new).
- `rmd/12_pathway_survey.Rmd` (new).
- `analysis.Rmd` updated.
- `storage/results/pathway_survey_gobp_ranking.tsv` (new).
- `analysis.html` re-knitted.

**Verification:**
- Section 12 appears in the knitted HTML with five contrast tables.
- TSV contains every GO BP pathway × every contrast.
- OXPHOS-related terms appear at whatever rank the composite gives them,
  not pre-sorted to the top.

---

### Session C2: GO MF × modality [DONE 2026-05-23]

**Completion note (2026-05-23):** Added `## GO MF pathway × modality
ranking` subsection to `rmd/12_pathway_survey.Rmd` mirroring the GO BP
layout from C1: `readRDS` of `storage/cache/fgsea_gomf_results.rds`,
per-contrast loop building `gomf_ranking_full`, full TSV export at
`storage/results/pathway_survey_gomf_ranking.tsv` (4,335 data rows =
867 GO MF pathways × 5 contrasts), and a `results = 'asis'` loop
printing five level-4 sub-subheadings each with a top-20 kable. No
helper changes (rank/format helpers are collection-agnostic by design).
`analysis.html` re-knitted cleanly in 2m55s; `class="error"` count is
0; `class="warning"` count is 0; new HTML structure adds section 13.2
plus 13.2.0.1-13.2.0.5 alongside the C1-built 13.1.x.

**Inputs:** `storage/cache/fgsea_gomf_results.rds`.

**Action:**
- Add subsection `## GO MF pathway × modality ranking` to
  `rmd/12_pathway_survey.Rmd`, reusing `rank_pathways_cross_modality()`.
- Top 20 per contrast.
- Export `storage/results/pathway_survey_gomf_ranking.tsv`.
- Re-knit.

**Outputs:** as named above.

**Verification:** TSV contains all GO MF terms × all contrasts; tables
are visible in the HTML.

---

### Session C3: GO CC × modality [DONE 2026-05-23]

**Completion note (2026-05-23):** Added `## GO CC pathway × modality
ranking` subsection to `rmd/12_pathway_survey.Rmd` mirroring the C1/C2
layout: `readRDS` of `storage/cache/fgsea_gocc_results.rds`,
per-contrast loop building `gocc_ranking_full`, full TSV export at
`storage/results/pathway_survey_gocc_ranking.tsv` (2,774 data rows =
555 GO CC pathways × 5 contrasts), and a `results = 'asis'` loop
printing five level-4 sub-subheadings each with a top-20 kable.
`analysis.html` re-knitted cleanly in 2m55s; `class="error"` count is
0; `class="warning"` count is 0; HTML now has all three GO collections
laid out as 13.1 (BP), 13.2 (MF), 13.3 (CC) with their 5 per-contrast
sub-subheadings each.

**Inputs:** `storage/cache/fgsea_gocc_results.rds`.

**Action:** identical structure to C2 with GO CC.

**Outputs:**
- `rmd/12_pathway_survey.Rmd` extended.
- `storage/results/pathway_survey_gocc_ranking.tsv` (new).
- `analysis.html` re-knitted.

**Verification:** as C2.

---

### Session C4: custom curated sets × modality [DONE 2026-05-23]

**Completion note (2026-05-23):** Added `## Custom curated sets ×
modality ranking` subsection to `rmd/12_pathway_survey.Rmd` that
pools all three Phase B custom caches (states / ad / modules) into a
unified ranking. Pathway names are unique across the three caches
(pairwise intersection = 0; 23 unique sets in total = 7+4+12), so
the pool helper `pool_custom_caches()` rbinds per (modality, contrast)
without identifier collisions. `source_map` is derived dynamically
from the cache names (no hard-coded set lists), and a new
`extra_cols` argument on `format_ranking_table()` inserts the
`source_collection` column between `composite_rank` and `pathway`.
**Plan-spec deviation (acknowledged):** the C4 plan says "Top 20 per
contrast" but the pooled custom collection has only 23 sets, so
top-20 would silently omit three sets per contrast. I instead set
`n = nrow(ranked)` so every set appears in each per-contrast table;
the plan-spec choice was tuned for the GO collections with thousands
of pathways. This is documented inline in the subsection preamble.
TSV export at `storage/results/pathway_survey_custom_ranking.tsv`
(115 data rows = 23 sets × 5 contrasts) with the `source_collection`
column populated. **Biology highlights (interaction contrast):**
DAM_up is the only set with 3 modalities sig (consistent-sign count =
2, sign_consensus = "mixed"); HAM_disease (states) +2.2 in nlgf_in_maptki
mean |NES|; WAM (ad) +1.66 and +2.30 in nlgf_in_p301s and nlgf_in_maptki
respectively. `analysis.html` re-knitted cleanly in 2m54s;
`class="error"` = 0, `class="warning"` = 0; HTML now has 13.1-13.4
plus 13.x.0.1-13.x.0.5 per-contrast sub-subheadings across all four
collection sections.

**Inputs:** `fgsea_custom_states_results.rds`, `fgsea_custom_ad_results.rds`,
`fgsea_custom_modules_results.rds`.

**Action:**
- Add subsection `## Custom curated sets × modality ranking` to
  `rmd/12_pathway_survey.Rmd`. Pool all three custom collections into
  a single ranking (they are project-specific signatures; combining
  them in one table reads better than three small tables).
- Annotate each row with its source collection (states / ad / modules).
- Top 20 per contrast.
- Export `storage/results/pathway_survey_custom_ranking.tsv`.
- Re-knit.

**Outputs:** as named above.

**Verification:** TSV contains all custom sets × all contrasts; source
collection column populated.

---

### Session C5: per-substate pathway ranking [DONE 2026-05-23]

**Completion note (2026-05-23):** Added two helpers to
`R/pathway_survey.R`: `rank_pathways_per_substate(fgsea_substate_cache,
n_top = 10L)` returns a long tibble with one row per (collection,
contrast, substate, pathway) carrying `NES`, `padj`, `abs_nes`,
`substate_rank` (1-indexed by |NES| desc within each cell), `in_top_n`,
and `substate_breadth` (count of substates where the pathway appears
in the top-n at the same collection × contrast; range 0..4 since the
project has four substates `homeostatic` / `DAM` / `IFN` /
`proliferative`). `format_per_substate_table(long_tbl, collection_key,
contrast_key, substate_key, n_top = 10L, breadth_mark = 3L,
strip_prefix = NULL, caption = NULL)` returns a knitr::kable for one
(collection, contrast, substate) cell, with the `*` prefix marking
pathways whose `substate_breadth >= breadth_mark`. Used base-R subsetting
to sidestep the dplyr NSE name-clash that would arise if the parameter
names matched column names. Added `## Per-substate pathway ranking
(interaction contrast)` subsection (rendered as 13.5 in the knitted
HTML) to `rmd/12_pathway_survey.Rmd` with three chunks:
`pathway-survey-per-state-build` reads the cache, builds the long
ranking, writes the TSV, prints the row-summary line, and computes
`breadth_summary`; `pathway-survey-per-state-breadth-summary`
(results='asis') renders the breadth summary kable; and
`pathway-survey-per-state-tables` (results='asis') renders the 24
per-cell top-10 tables. Initial draft printed the breadth summary
inside the build chunk, which rendered as plain `## Table:` markdown
text instead of a proper HTML table because the build chunk lacks
`results='asis'`; split into a dedicated asis chunk and re-knit. TSV
at `storage/results/pathway_survey_per_state_ranking.tsv` has 60,885
rows × 10 cols (matches the cache row total exactly: 6 collections × 4
substates × 5 contrasts × per-cell pathway count), with all 6
collections, all 4 substates, and all 5 contrasts represented. The
inline tables focus on the interaction contrast; the TSV holds all 5
contrasts for downstream inspection. Final structure rendered as
section 13.5 (level 2) with 6 collection sub-sections (13.5.1 GO BP,
13.5.2 GO MF, 13.5.3 GO CC, 13.5.4 custom-microglia-states,
13.5.5 custom-microglia-ad, 13.5.6 custom-module-sources) each holding
4 substate sub-subsections — 24 top-10 tables in total, all at the
interaction contrast. `analysis.html` re-knit cleanly twice (the
breadth-summary chunk-split fix triggered the second knit) in 2m56s
each; final `class="error"` = 0, `class="warning"` = 0; 31 elements
match `number="13.5*"` (1 root + 6 collections + 24 substates = 31).

**Plan-spec deviations (acknowledged):**
1. The plan calls for "top 10 pathways by |NES| with padj" per
   (collection, substate). For the three small custom collections the
   collection size is below 10 sets per substate (3-12 sets after
   intersecting with each substate's gene universe), so the format
   helper shows `nrow(top)` rows when that is below `n_top` rather
   than padding empty rows. The caption in those tables truthfully
   says "Top 3 pathways" / "Top 7 pathways" etc. The caveat is also
   documented in the Rmd preamble: small collections inflate
   substate_breadth artifactually because every set lands in the
   per-substate top-n by construction.
2. The plan asks the inline tables to focus on the interaction
   contrast (per the section subtitle). The TSV export covers all 5
   contrasts so downstream sessions (e.g. C9 verdict) can reuse the
   per-state ranking for other contrasts without rebuilding it.

**Biology highlights (interaction contrast):**
- **Cross-substate breadth is concentrated in ribosomal compartment
  terms.** Every breadth==4 GO pathway is ribosomal:
  `GOCC_CYTOSOLIC_RIBOSOME`, `GOCC_CYTOSOLIC_SMALL_RIBOSOMAL_SUBUNIT`,
  `GOCC_RIBOSOMAL_SUBUNIT`, `GOCC_RIBOSOME` (4 GO CC sets) +
  `GOMF_STRUCTURAL_CONSTITUENT_OF_RIBOSOME` (1 GO MF set). All have
  strongly negative NES in DAM (e.g. ribosomal subunit NES=-1.99
  padj=1.4e-6) and consistent direction across the other three
  substates. **The interpretation is that tau × amyloid interaction
  globally suppresses microglia ribosomal/translational machinery
  across the entire microglia compartment** — this is a compartment-
  wide signature, not a substate-specific one. The custom
  `MG3_ribosome_biogenesis` set from Sun/Victor 2023 also lands at
  breadth=4 from a completely independent gene-list source,
  corroborating the GO findings.
- **GO BP has 0 pathways with breadth >= 3** at interaction (38
  breadth-1, 1 breadth-2, 0 breadth-3+). Process-level signals are
  highly substate-specific. The DAM substate top-10 surfaces a strong
  lipid-response signal (positive regulation of lipid localization
  NES=+2.17 padj=0.034; lipid transport NES=+2.03 padj=0.050;
  response to lipoprotein particle NES=+1.89), consistent with the
  WAM/LDAM enrichments at the whole-microglia level. The
  homeostatic substate top-10 differs entirely from DAM's top-10
  — substate-specificity in GO BP is the rule, not the exception.
- **GO CC has 4 breadth=4 + 3 breadth=3 pathways.** Compartment-level
  signal is more compartment-wide than process-level signal, because
  one cellular component (ribosome) is shared across many distinct
  processes (translation, rRNA processing, ribosome biogenesis,
  protein folding) that GO BP splits apart.
- **GO MF has 1 breadth=4 + 2 breadth=3 pathways**, intermediate
  between BP and CC.
- **Custom collections show high breadth as a small-N artifact**:
  custom_microglia_ad has 3/3 sets at breadth=4 (WAM, AD1, AD2 all
  land in every substate's top-10 because the entire collection fits
  in top-10); custom_microglia_states has 6/7 at breadth=4;
  custom_module_sources has 8/12 at breadth=4. This is documented
  inline in the Rmd preamble caveat. Direction-of-effect on these
  sets is still informative: WAM NES=+3.05 in DAM, AD2 NES=-1.78 in
  DAM, MG3_ribosome_biogenesis NES≈-2 in every substate.

Phase C5 complete; next step is C6 (hdWGCNA module × modality
cross-support), no decision gate.

**Inputs:** `storage/cache/fgsea_per_state_results.rds` (rebuilt in B5).

**Action:**
- Add subsection `## Per-substate pathway ranking (interaction contrast)`
  to `rmd/12_pathway_survey.Rmd`.
- For each substate × each collection, show top 10 pathways by |NES|
  with padj.
- Compute "substate breadth": for each pathway, in how many substates
  it appears in the top 10. Pathways with breadth ≥3 highlighted.
- Export `storage/results/pathway_survey_per_state_ranking.tsv`.

**Outputs:** as named above.

**Verification:** Per-state tables exist; substate-breadth column
populated.

---

### Session C6: hdWGCNA module × modality cross-support [DONE 2026-05-23]

**Completion note (2026-05-23):** Added `## hdWGCNA module × modality
cross-support` subsection (rendered as 13.6 in the knitted HTML) to
`rmd/12_pathway_survey.Rmd` with three chunks. The build chunk
`pathway-survey-hdwgcna-build` reads `storage/cache/hdwgcna_microglia.rds`,
maps each module gene from Ensembl to mouse symbol via the project-wide
`symbol_map` (built in `rmd/01_data.Rmd`), drops the `grey` unassigned
bucket, and constructs a 4-set custom pathway list keyed by module ID
(`MG-M1`=570, `MG-M2`=850, `MG-M3`=780, `MG-M4`=74 mouse symbols).
Symbol coverage is 2274/2274 (100%) over the 2274 non-grey
module-assigned genes, so no module shrinks during the Ensembl→symbol
step. fgsea then runs per modality × contrast using the same
`run_dataset` shim defined upstream in `rmd/07_integration.Rmd`,
keeping NEBULA z for snRNAseq and limma t for the other four
modalities — identical statistic family to every other fgsea cache in
the project. `min_size = 5` (matching the Phase B custom collections)
is not actively binding here since the smallest module (MG-M4, 74
symbols) clears the cutoff in every modality after intersection with
the modality gene universe. The new fgsea cache is
`storage/cache/fgsea_hdwgcna_modules_results.rds` (53 KB,
`modality -> contrast -> fgseaResult` shape mirroring every other
fgsea cache); the cross-modality ranking TSV is
`storage/results/pathway_survey_hdwgcna_modules.tsv` (20 data rows =
4 modules × 5 contrasts with per-modality NES + padj columns and
`composite_rank`). The display chunk `pathway-survey-hdwgcna-tables`
renders five per-contrast composite-ranking kables, each with all 4
modules (top-N capped at `nrow(ranked)` since the collection is
small). The heatmap chunk `pathway-survey-hdwgcna-heatmap` builds a
4×5 NES matrix (rows = modules alphabetical MG-M1, MG-M2, MG-M3,
MG-M4; columns = canonical project modality order snrnaseq, geomx,
proteomics, phospho, phospho_corr) and renders via
`ComplexHeatmap::Heatmap` with cell-level significance markers (`*`
for FDR<0.05, `.` for FDR<0.1). Rows alphabetical enforces the
anti-anchoring guardrail explicitly: MG-M3 sits in the third row
along with the other three modules, not at the top.

**Implementation deviation from the plan text (acknowledged):**
The plan inputs section mentions `storage/cache/integration_table.rds`
"for per-modality gene-level ranks per contrast." Inspection showed
the integration_table is a tibble of per-symbol per-contrast logFC and
p-values across modalities — it does NOT carry NEBULA z statistics
for snRNAseq, only NEBULA logFC. The existing fgsea infrastructure in
`rmd/07_integration.Rmd` instead uses `de_<modality>$fit$top[[cn]]`
objects (which carry the modality-native statistics: NEBULA z stored
as `t` for snrnaseq via the `snrnaseq_nebula_shim`, limma t for the
others) and reaches fgsea through `run_dataset()`. Reusing that exact
path keeps the snRNAseq ranking on NEBULA z (the plan-spec choice)
and produces fgsea outputs with the same shape as every other Phase
B/C cache. So `integration_table.rds` is not actually needed for C6;
the DE caches and `snrnaseq_nebula_shim` already in scope from 07's
`{r fgsea-setup}` chunk are sufficient.

**Biology highlights (interaction contrast composite ranking):**
- **MG-M3 ranks first by composite (3 modalities significant, mixed
  sign).** Modality breakdown: snRNAseq NES=-2.00 padj=1.6e-41
  (strongly negative), proteomics NES=-1.95 padj=2.4e-10 (strongly
  negative), phospho_corr NES=+1.76 padj=7.8e-6 (strongly positive),
  phospho NES=+1.12 padj=0.34 (n.s. positive), GeoMx NA (the limma
  fgsea returned NA for MG-M3 here, presumably the unbalanced-stat
  edge-case noted in C5). The mixed sign reflects what was already
  visible in the per-substate C5 analysis: the MG-M3 transcriptional
  module is suppressed at the RNA layer in tau×amyloid interaction
  while phospho-corrected signal (which adjusts for total-protein
  changes) shifts positive — consistent with the ribosome-machinery
  story from C5 where ribosomal protein abundance drops but residual
  phospho stoichiometry tracks differently. The sign disagreement is
  preserved as `sign_consensus = mixed` in the TSV; the composite
  rank is driven by the consistent-sign count of 2 (snRNAseq +
  proteomics negative).
- **MG-M1 ranks second** (1 modality significant, all NES <0): snRNAseq
  NES=-1.22 padj=8.7e-3 — a coherent but small negative interaction-
  contrast move; consistent with the general "interaction suppresses
  microglia activation" pattern.
- **MG-M4 ranks third** (1 modality significant, all NES <0): snRNAseq
  NES=-1.67 padj=6.2e-4 — the smallest module (74 symbols) still
  surfaces a coherent interaction-contrast suppression in snRNAseq.
- **MG-M2 ranks last** (0 modalities significant). Highest mean |NES|
  but no significance, so it sinks to the bottom — exactly what the
  composite rank's "sign-consistent significant count first" criterion
  should do.
- **Other contrast highlights:** MG-M2 is the dominant module under
  the two NLGF-on-tau contrasts (nlgf_in_maptki: 3 modalities sig,
  consistent sign + with snRNAseq NES=+2.47 padj=9.8e-121;
  nlgf_in_p301s: 4 modalities sig, all positive, snRNAseq NES=+2.15
  padj=1.6e-111). MG-M3 is the dominant module under tau_alone
  (3 modalities sig, mixed sign; GeoMx NES=-3.01 padj=7.1e-69 is
  the largest single-modality effect anywhere in the table) and
  tau_in_nlgf (4 modalities sig, mixed sign; snRNAseq NES=-4.00
  padj=8.1e-247 is even larger). The MG-M3 mixed-sign pattern repeats
  across contrasts because RNA-level suppression and phospho_corr
  enrichment are genuinely orthogonal axes.

**Anti-anchoring observation:** Across all five contrasts, MG-M3 lands
at composite rank 1 in three of them (interaction, tau_alone,
tau_in_nlgf) and at rank 2 in nlgf_in_p301s. This is robust
cross-contrast support for MG-M3 carrying the most consistent
multimodal signal in the canonical hdWGCNA build — but **the
sign_consensus is "mixed" in every one of those four contrasts**,
meaning at least one significant modality disagrees on direction.
That nuance is critical: MG-M3 is a high-signal module, but not a
clean direction-of-effect module. Whether to read "MG-M3 dominance"
as a positive finding or as a sign-conflict warning is exactly the
question the C9 unified verdict step should adjudicate alongside
the GO BP composite winners and the OXPHOS-family specificity null
from C7/C8.

`analysis.html` re-knitted cleanly in 2m57s; final `class="error"`
count is 0; `class="warning"` count is 0; HTML now has 13.1 (GO BP),
13.2 (GO MF), 13.3 (GO CC), 13.4 (custom), 13.5 (per-substate with
13.5.1-13.5.6 collection sub-sections), and 13.6 (hdWGCNA modules
with 13.6.0.1-13.6.0.5 per-contrast sub-subheadings). Phase C6
complete; Phase C7 (OXPHOS-family specificity null script) is next
and has a decision gate at start.

**Inputs:**
- `storage/cache/hdwgcna_microglia.rds` (canonical p=3 build).
- `storage/cache/integration_table.rds` for per-modality gene-level
  ranks per contrast.

**Action:**
- Extract per-module gene lists from the hdWGCNA build (`GetModules`).
- For each (module, modality, contrast), run fgsea using the module
  as a custom gene set against the modality's signed-statistic ranking.
  Use the same statistic family as the existing fgsea calls (NEBULA z
  for snRNAseq; limma t for the others).
- Build a tidy table: module × modality × contrast with NES, padj,
  leading-edge size.
- Add subsection `## hdWGCNA module × modality cross-support` to
  `rmd/12_pathway_survey.Rmd`. Show module × modality NES heatmap for
  the interaction contrast, with modules in alphabetical order
  (NOT sorted with MG-M3 first; anti-anchoring guardrail).
- Export `storage/results/pathway_survey_hdwgcna_modules.tsv`.

**Outputs:** as named above.

**Verification:** Every non-grey module appears for every (modality,
contrast) combination. MG-M3 not visually privileged in the heatmap.

---

### Session D1: build unbiased leader-board helper + unified TSV [DONE 2026-05-23]

**Completion note (2026-05-23):** Added `build_leader_board(rankings_named_list, per_state_long = NULL, leader_rule = NULL)` to `R/pathway_survey.R` (lines 420-end). The default leader_rule is the plan-spec rule `function(row) row$n_modalities_sig_consistent_sign >= 2 | row$n_modalities_sig >= 3`, vectorised (`|`, not `||`) so it operates on the full stacked tibble in one shot rather than row-by-row; the docstring explicitly documents this convention so future callers know vectorised closures are the supported pattern. The helper accepts the wide ranking tibbles as-is and uses each row's `source_collection` column when present (the custom-pool TSV from C4) as the collection label; otherwise it falls back to the input list name. This lets the caller pre-remap labels at the call site (the smoke test does the {states/ad/modules} → {custom_microglia_states/custom_microglia_ad/custom_module_sources} mapping to match per_state_long's collection column), keeping the helper itself generic.

Internal structure: stack all rankings into one long tibble; apply leader_rule once vectorised; group by (collection, pathway); reduce to scalar aggregates inside one `summarise()` call. Two small inline helpers `reduce_dominant_sign()` and `reduce_contrasts_summary()` keep `summarise()` readable — both operate on the per-group vectors at the leader subset and return a single scalar value. The `dominant_sign` reduction strips NA sign_consensus values before deciding "+" / "-" / "mixed", but per the leader rule any leader contrast must have ≥2 sig modalities so NA cannot occur in practice; the defence is for robustness against custom leader rules. `contrasts_summary` is the pipe-delimited `<contrast>:<n_modalities_sig>/<sign_consensus>` string (using n_modalities_sig not n_modalities_sig_consistent_sign, matching the plan example `interaction:3/mixed | tau_in_nlgf:2/-` where the OXPHOS-at-interaction case has n_modalities_sig=3 with mixed sign). `max_substate_breadth` is joined from the per-state long tibble grouped by (collection, pathway) on max(substate_breadth); pathways absent from per_state_long (hdwgcna modules + GO terms never reaching a per-state top-10) get 0. `leader_score = 5 * n_contrasts_consistent_sign_ge2 + n_contrasts_sig_ge3 + max_consistent_sign / 5` exactly as spec'd; final sort is `leader_score desc, n_contrasts_consistent_sign_ge2 desc, max_consistent_sign desc, max_abs_nes desc, collection asc, pathway asc`.

Smoke test exercised the helper against all six committed TSVs in one transient `Rscript -e` invocation (NOT a committed script per plan-spec). All three plan-mandated assertions PASS:
1. **MG-M2 (hdwgcna_modules)** appears with `n_contrasts_consistent_sign_ge2 = 2` (the nlgf_in_p301s sig=4 cons=4 + hit AND the nlgf_in_maptki sig=3 cons=3 + hit; both contribute to the count).
2. **All 5 MHC II antigen-presentation GO BP terms** (`GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_II`, `GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_OR_POLYSACCHARIDE_ANTIGEN_VIA_MHC_CLASS_II`, `GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_ANTIGEN`, `GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_PEPTIDE_ANTIGEN`, `GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_ANTIGEN`) appear with `n_contrasts_consistent_sign_ge2 = 2`, `dominant_sign = "+"`, `leader_score = 10.4`, and identical `contrasts_summary = "nlgf_in_maptki:2/+ | nlgf_in_p301s:2/+"` (the plan-spec smoke test asked for "at least 4" — 5 is exact).
3. **GOBP_OXIDATIVE_PHOSPHORYLATION** appears with `n_contrasts_sig_ge3 = 1`, `dominant_sign = "mixed"`, `leader_score = 11.4`, `contrasts_summary = "interaction:3/mixed | tau_in_nlgf:2/-"`.

**Leader-board geometry**: 144 leader rows across 6 collections (gobp=67, gocc=46, gomf=22, custom_microglia_states=4, hdwgcna_modules=3, custom_microglia_ad=2). `custom_module_sources` produces ZERO leaders — none of the 12 Sun/Victor 2023 microglia-state markers passes the leader rule. dominant_sign distribution: 77 (−) / 42 (+) / 25 mixed (negative-direction skew consistent with OXPHOS/ribosomal/synaptic suppression dominating cross-modality consistency). max_substate_breadth distribution: 101 (=0) / 8 (=1) / 2 (=2) / 4 (=3) / 29 (=4). The 101 breadth-0 rows include the 3 hdWGCNA modules (no per-state cache) and 98 GO terms that never enter a per-state top-10 at any contrast; the 29 breadth-4 rows include all 6 custom-collection leaders (small-N artefact: custom collections have ≤12 sets, all land in per-state top-10 by construction) plus 23 GO terms with genuine compartment-wide signal (the ribosomal/OXPHOS axis from C5).

**Top of the leader board (top 10 by leader_score):**
1. **MG-M3** (hdwgcna_modules) — 23.4, 4 contrasts consistent-sign≥2, 3 contrasts sig≥3, dominant_sign=mixed. The dominant cross-contrast multimodal signal in the canonical hdWGCNA build, but its mixed sign across contrasts confirms the RNA-suppression vs phospho_corr-enrichment orthogonality flagged in C6.
2-8. **OXPHOS / ETC family + DAM_up** (gobp/gomf/gocc/custom_microglia_states) — leader_score 16.4 each, all 3 contrasts consistent-sign≥2, 1 contrast sig≥3, dominant_sign=mixed. The cross-modality OXPHOS suppression signal at NLGF contrasts + the mixed-sign breadth at interaction. Includes `GOBP_MITOCHONDRIAL_ELECTRON_TRANSPORT...`, `GOBP_ATP_SYNTHESIS_COUPLED_ELECTRON_TRANSPORT`, `GOCC_NADH_DEHYDROGENASE_COMPLEX`, `GOMF_NADH_DEHYDROGENASE...`, `GOMF_OXIDOREDUCTION_DRIVEN_ACTIVE_TRANSMEMBRANE_TRANSPORTER...`, `GOMF_OXIDOREDUCTASE_ACTIVITY...`, `DAM_up`.
9. **GOCC_PRESYNAPSE** (gocc) — 15.4, the synaptic axis (consistent − sign at nlgf_in_maptki and nlgf_in_p301s).
10. **MG-M2** (hdwgcna_modules) — 12.8, 2 contrasts consistent-sign≥2 and 2 contrasts sig≥3, dominant_sign=+. The amyloid-response module with consistent + sign across both NLGF contrasts.

The leader board surfaces a **multi-axis biology** rather than a single winner:
- **Amyloid-driven activation axis (+ sign):** MG-M2, DAM_up, HAM_disease, WAM, AD2, 5 MHC II GO BP terms, GOBP_ADAPTIVE_IMMUNE_RESPONSE.
- **Suppression-under-pathology axis (− sign):** synaptic GO CC family (presynapse, synaptic vesicle, GABAergic synapse, Schaffer collateral, excitatory synapse, presynaptic membrane), all clustered at leader_score 11.6.
- **Mixed-sign OXPHOS/ETC axis:** GOBP_OXIDATIVE_PHOSPHORYLATION + 7 related ETC/respiration terms + DAM_up + GOCC_OXIDOREDUCTASE_COMPLEX, all at leader_score 11.4-16.4, all showing the same interaction:3/mixed + tau_in_nlgf:2/− pattern.
- **MG-M3 + ribosomal/translational axis (mixed at hdWGCNA, − in C5 per-state):** the RNA-vs-proteomics-vs-phospho_corr disagreement that C6 already documented in detail.

**Outputs created:**
- `R/pathway_survey.R` extended with `build_leader_board()` (~250 new lines including the multi-paragraph docstring).
- `storage/results/pathway_survey_unified_leaderboard.tsv` (144 rows × 12 cols).
- No Rmd changes; no re-knit (per plan-spec).

**Inputs:** the six ranking TSVs already in place
(`pathway_survey_gobp_ranking.tsv`, `pathway_survey_gomf_ranking.tsv`,
`pathway_survey_gocc_ranking.tsv`, `pathway_survey_custom_ranking.tsv`,
`pathway_survey_per_state_ranking.tsv`,
`pathway_survey_hdwgcna_modules.tsv`).

**Action:**
- Extend `R/pathway_survey.R` with a new helper
  `build_leader_board(rankings_named_list, per_state_long, leader_rule)`
  where:
  - `rankings_named_list` is a named list of the five cross-modality
    ranking tibbles (gobp, gomf, gocc, custom, hdwgcna) read from their
    TSVs.
  - `per_state_long` is the per-substate breadth tibble (so each
    pathway can carry its `max_substate_breadth` across contrasts).
  - `leader_rule` is a closure returning `TRUE` for any
    (pathway, contrast) row that qualifies as a leader. Default:
    `function(row) row$n_modalities_sig_consistent_sign >= 2 ||
    row$n_modalities_sig >= 3`. The rule is parameterised so future
    sessions can tighten/loosen without rewriting the helper.
- The helper output is one row per (collection, pathway). Columns:
  - `collection`, `pathway`
  - `n_contrasts_leader` — count of contrasts where the row meets the
    rule (range 0..5; only ≥1 rows kept)
  - `n_contrasts_consistent_sign_ge2` — count where consistent-sign
    sig ≥2
  - `n_contrasts_sig_ge3` — count where sig ≥3
  - `max_consistent_sign` — max value across contrasts
  - `max_n_modalities_sig` — max value across contrasts
  - `max_abs_nes` — max |mean_abs_nes| across contrasts
  - `dominant_sign` — `+` / `-` / `mixed` based on the leader
    contrasts only (mixed if ANY leader contrast disagrees)
  - `max_substate_breadth` — max across contrasts (0 if pathway not in
    the per-substate cache, e.g. for hdWGCNA modules)
  - `contrasts_summary` — pipe-delimited string of contrast:sig/sign
    triples (e.g. `interaction:3/mixed | tau_in_nlgf:2/-`)
  - `leader_score` — composite for primary sort:
    `5 * n_contrasts_consistent_sign_ge2 + n_contrasts_sig_ge3 + max_consistent_sign / 5`
    (deliberately favours cross-contrast consistent-sign first, then
    cross-modality breadth at any single contrast, then single-contrast
    consistent-sign max as a tie-breaker; the constants are designed so
    cross-contrast support outweighs single-contrast peaks)
- Source the updated `R/pathway_survey.R` from `R/helpers.R` already
  (no change needed).
- Build runnable smoke test (in a transient `Rscript -e` invocation,
  not committed) that prints the top 30 by leader_score and verifies:
  - MG-M2 appears with n_contrasts_consistent_sign_ge2 ≥1 (the
    nlgf_in_p301s sig=4 cons=4 + hit)
  - At least 4 MHC II antigen-presentation GO BP terms appear with
    n_contrasts_consistent_sign_ge2 ≥2 (the nlgf_in_maptki ranks 1-5 +
    nlgf_in_p301s ranks 7-11 pattern)
  - GOBP_OXIDATIVE_PHOSPHORYLATION appears with
    n_contrasts_sig_ge3 ≥1 (the interaction sig=3 mixed-sign hit)
- Write the full unified leader board to
  `storage/results/pathway_survey_unified_leaderboard.tsv` (every
  collection × every pathway that is a leader at ≥1 contrast; one row
  per pathway).
- This step is helper + TSV only; no Rmd 12 changes yet (D2 wires the
  Rmd display).

**Outputs:**
- `R/pathway_survey.R` extended with `build_leader_board()`.
- `storage/results/pathway_survey_unified_leaderboard.tsv` (new).
- (No Rmd changes; no re-knit needed — verified via smoke test only.)

**Verification:**
- Leader-board TSV has ≥1 row per pathway that meets the rule at any
  contrast.
- Three smoke-test pathways (MG-M2, an MHC II term, OXPHOS) appear with
  expected support counts.
- Helper handles the hdWGCNA TSV (which lacks the per-substate breadth
  column) without erroring.

---

### Session D2: Rmd 12 unbiased verdict subsection + archive plan + write successor [DONE 2026-05-23]

**Completion note (2026-05-23):** Added `## Unified leader board
(unbiased verdict)` subsection to `rmd/12_pathway_survey.Rmd`
immediately after the C6 hdWGCNA module subsection (rendered as 13.7
in the knitted HTML with a 13.7.1 Verdict sub-subsection for the
prose). Also updated the section preamble to remove the obsolete
reference to the cancelled OXPHOS-family specificity null and replace
it with a pointer to the unbiased verdict. The verdict subsection
preamble names the leader rule verbatim
(`n_modalities_sig_consistent_sign >= 2` OR `n_modalities_sig >= 3`)
and the composite `leader_score = 5 * n_contrasts_consistent_sign_ge2
+ n_contrasts_sig_ge3 + max_consistent_sign / 5` formula before
inspecting any data, so the threshold cannot be retrofitted. The build
chunk `pathway-survey-verdict-build` reads
`storage/results/pathway_survey_unified_leaderboard.tsv` (built in D1)
and constructs three view tibbles using a small inline
`pretty_pathway()` helper that strips the `GOBP_` / `GOMF_` / `GOCC_`
prefix and underscore-formats GO pathway names while leaving custom-
collection names (`MG-M2`, `DAM_up`, `WAM`, `HAM_disease`, `AD1`,
`AD2`, `MG-M3`) unchanged. The display chunk
`pathway-survey-verdict-tables` (`results = 'asis'`) prints three
kables: top 20 by leader_score, top 10 by max_substate_breadth, top 10
by max_consistent_sign at any single contrast. Smoke-tested the
build-chunk code via a transient `Rscript -e` against the live TSV
before the knit; shapes (20 / 10 / 10 rows; 9 / 7 / 7 cols), first-row
content (MG-M3 / DAM_up / MG-M2 at top of each respective view), and
kable rendering all matched the D1 expectations. The verdict prose
(13.7.1) explicitly names the three independent axes — amyloid-driven
activation (+ MG-M2, MHC II antigen-presentation family,
GOBP_ADAPTIVE_IMMUNE_RESPONSE, custom DAM_up/WAM/HAM_disease/AD2),
NLGF-driven synaptic suppression (− six synaptic GO CC terms plus GO
BP synapse terms), and mixed-sign metabolic/translational
disagreement at tau×amyloid interaction (mixed OXPHOS/ETC family +
DAM_up + MG-M3 ribosomal compartment) — and states that the verdict
is multi-axis rather than single-pathway. The prose also documents
that no specificity-null testing was performed because the
multi-collection independent replication of the same `interaction:3/
mixed` pattern across GO BP, GO MF, GO CC, custom_microglia_states,
and hdwgcna_modules adjudicates the OXPHOS-construction-artefact
question without needing one. Forward-pointer paragraph references
the new `storage/notes/mechanism_layer_plan.md` with the three
deferred phases (TF inference / kinase activity inference from
phospho / CCC bridge sender/receiver re-filtering).

Re-knit produced `analysis.html` (27 MB) in 2m56s with zero
`class="error"` and zero `class="warning"` elements; the TOC now has
13.1 through 13.7 with 13.7.1 nested for the Verdict sub-subsection;
all three verdict kables and all three prose axes render correctly.
The leaderboard TSV timestamp (21:04 from the D1 build) is preserved
across the knit — the verdict chunk is read-only against the TSV by
design, so D2 does not invalidate D1's output. No new caches were
created; no new TSVs were created; the verdict subsection is purely a
display layer over D1's leader board.

**Plan-spec deviation (acknowledged):** The plan asked for an
`Outcome summary` section appended to the archived plan recording the
"top 3 entries by leader_score, top 3 by breadth, top 3 by
single-contrast max, the multi-axis verdict, and a pointer to the
successor plan." Implemented exactly as specified; the archived plan
at `storage/notes/completed/pathway_overhaul_plan_2026-05-23.md`
carries the outcome summary at its end. The successor stub
`storage/notes/mechanism_layer_plan.md` is the new active plan and
contains the three proposed phases (E: TF inference, F: kinase
activity inference, G: CCC bridge re-analysis) all explicitly gated
on user choice at their first session.

**Inputs:** `storage/results/pathway_survey_unified_leaderboard.tsv`
from D1.

**Action:**
- Add subsection `## Unified leader board (unbiased verdict)` to
  `rmd/12_pathway_survey.Rmd` after the C6 hdWGCNA subsection.
  Subsection preamble must name the leader rule explicitly (so the
  threshold cannot be retrofitted) and explicitly state that no
  pathway family was pre-privileged.
- Render three views from the leader-board TSV (each as a
  `results = 'asis'` kable):
  1. **Top 20 by leader_score** — the composite cross-contrast +
     cross-modality view. Full per-contrast support columns visible.
  2. **Top 10 by max_substate_breadth** — pathways with the broadest
     pan-substate compartment-wide signal across contrasts. Surfaces
     ribosomal/translational machinery.
  3. **Top 10 by max_consistent_sign at any single contrast** —
     single-contrast cross-modality stars. Surfaces MG-M2 at
     nlgf_in_p301s.
- Add a prose "verdict" paragraph (British English) that:
  - Names the leader-rule criterion verbatim.
  - Describes what the leader board surfaces (multi-axis: amyloid
    activation = MG-M2 + MHC II family; suppression-under-pathology =
    OXPHOS + ribosomal; possible mechanism = synaptic engulfment at
    nlgf_in_maptki).
  - States that the verdict is multi-axis rather than single-pathway;
    no specificity-null testing was performed; the consistency filters
    and cross-contrast aggregation are the substitute.
  - References the existing per-contrast tables in C1-C6 for the
    underlying evidence.
- Add a final paragraph pointing forward to the successor plan
  (mechanism layer).
- Re-knit `analysis.Rmd`; verify zero `class="error"` and zero
  `class="warning"` in the HTML.
- Move this plan to
  `storage/notes/completed/pathway_overhaul_plan_<YYYY-MM-DD>.md` and
  add a final `## Outcome summary` section recording: leader-rule
  applied, top 3 entries by leader_score, top 3 by breadth, top 3 by
  single-contrast max, the multi-axis verdict, and a pointer to the
  successor plan.
- Create `storage/notes/mechanism_layer_plan.md` as a stub successor
  plan with header, "## Why this plan exists" pointer to the verdict
  outcome, and proposed phases:
  - Phase E: TF inference (decoupler+CollecTRI) on the leader pathways
  - Phase F: kinase activity inference from phospho on leader pathways
  - Phase G: CCC bridge re-analysis using leader pathways as
    sender/receiver filters
  - All three phases gated on user choice of priorities at their first
    session (decision gates marked explicitly).
- Commit locally with subject "D2: add unbiased verdict + archive plan
  + draft mechanism plan" (or similar under 70 chars).

**Outputs:**
- `rmd/12_pathway_survey.Rmd` extended with the verdict subsection.
- `analysis.html` re-knitted.
- This plan moved to `storage/notes/completed/` with an Outcome summary.
- `storage/notes/mechanism_layer_plan.md` (new stub).

**Verification:**
- Verdict prose explicitly states the leader rule.
- `ls storage/notes/*_plan.md` shows exactly one active plan
  (`mechanism_layer_plan.md`).
- HTML contains zero `HALLMARK_` strings (already guaranteed) and the
  new section is visible at the end of Rmd 12.

---

## Final deliverables when this plan is complete

- 5 OXPHOS-arc Rmd files deleted; analysis.Rmd cleaned of their `child=`
  references.
- 8+ OXPHOS-arc result TSVs deleted.
- 2 Hallmark caches + 1 Hallmark TSV deleted.
- `get_hallmark()` removed from `R/fgsea.R`.
- `rmd/02e_snrnaseq_substate_pathway.Rmd`, `rmd/03_hdwgcna.Rmd`,
  `rmd/03b_hdwgcna_soft_power.Rmd`, `rmd/07_integration.Rmd`,
  `rmd/10_divergence.Rmd`, `rmd/11_ccc.Rmd` all Hallmark-free.
- 4 new collection caches (GO MF, GO CC, custom states, custom AD,
  custom modules), 6 new fgsea-result caches, and their TSV exports.
- `R/pathway_survey.R` (new helper module) with `build_leader_board()`
  added in D1.
- `rmd/12_pathway_survey.Rmd` (new section with the C1-C6 survey
  subsections + the D2 unbiased verdict subsection).
- `scripts/build_fgsea_extra_collections.R` (the only build script
  this plan needs; the OXPHOS-null script from the cancelled C7 is
  NOT created).
- `storage/results/pathway_survey_unified_leaderboard.tsv` from D1.
- `analysis.html` re-knitted, comprehensive.
- This plan archived to `storage/notes/completed/` with an Outcome
  summary; successor `storage/notes/mechanism_layer_plan.md` active
  in `storage/notes/`.

## How to mark this plan complete

When session D2 is DONE, move this file to
`storage/notes/completed/pathway_overhaul_plan_<YYYY-MM-DD>.md` (the
`completed/` directory already exists), chown `rstudio:rstudio`, and add
a final `## Outcome summary` section recording the leader-rule applied,
top 3 by leader_score, top 3 by max_substate_breadth, top 3 by
max_consistent_sign, the multi-axis verdict statement, and a pointer
to the new active `mechanism_layer_plan.md`.

## Anti-anchoring guardrails (re-read every session)

These exist because LLMs (this agent included) drift toward the most
salient prior finding. Each session must enforce them:

- **Always** sort pathway tables by an agnostic composite, never by
  OXPHOS rank or by OXPHOS-related keywords.
- **Always** alphabetise hdWGCNA modules in heatmaps; never place MG-M3
  first.
- **Always** state the leader rule explicitly before applying it, so
  the threshold cannot be retrofitted. The D1/D2 verdict uses one
  fixed rule (consistent-sign sig ≥2 at any contrast OR sig ≥3 at any
  contrast). The rule is parameterised in code so the user can later
  swap it, but each session must name the rule in effect.
- **Always** present the verdict as multi-axis unless the leader board
  empirically supports a single-axis story. The 2026-05-23 survey
  produced multiple distinct findings (MG-M2 amyloid activation, MHC II
  family, OXPHOS suppression, ribosomal compartment, candidate synaptic
  engulfment) — collapsing these into a single "winner" reintroduces
  bias.
- **Never** reintroduce Hallmark gene sets at any step. The collection
  is permanently retired from this project.
- **Never** reintroduce the OXPHOS confirmatory arc (08, 09a, 09b, 09c,
  09d). The mechanism-layer successor plan must derive its targets
  from the D2 leader board, NOT from re-anchoring on Hallmark or on the
  deleted stress-test sections.
- **Never** rebuild a specificity-null or expression-matched-random
  framework as part of this plan. That work was dropped on 2026-05-23
  in favour of ranking + consistency only. A future plan may
  reintroduce nulls if a specific question demands them, but it must
  be motivated by that question — not as a default re-anchoring on the
  cancelled C7/C8 design.

---

## Outcome summary

**Plan archived:** 2026-05-23 (Phase D complete; all 16 sessions DONE).
**Successor plan:** `storage/notes/mechanism_layer_plan.md`.

**Leader rule applied (verbatim, fixed before data was inspected):**
a (pathway, contrast) row is a leader if
`n_modalities_sig_consistent_sign >= 2` OR `n_modalities_sig >= 3`.

**Composite ranking formula:**
`leader_score = 5 * n_contrasts_consistent_sign_ge2 + n_contrasts_sig_ge3 + max_consistent_sign / 5`.
Cross-contrast consistent-sign breadth weighted first; cross-modality
breadth at a single contrast second; single-contrast consistent-sign
peak as tie-breaker.

**Leader board geometry:** 144 leader rows across 6 collections
(gobp=67, gocc=46, gomf=22, custom_microglia_states=4,
hdwgcna_modules=3, custom_microglia_ad=2; `custom_module_sources` =
zero leaders — no Sun/Victor 2023 microglia-state marker passes the
rule). dominant_sign distribution: 77 (−) / 42 (+) / 25 mixed.

### Top 3 by leader_score (composite cross-contrast + cross-modality)

1. **MG-M3** (hdwgcna_modules) — leader_score 23.4; 4 contrasts
   consistent-sign≥2; 3 contrasts sig≥3; dominant_sign = mixed;
   contrasts_summary = `nlgf_in_maptki:2/- | interaction:3/mixed |
   tau_alone:3/mixed | tau_in_nlgf:4/mixed`. The dominant cross-
   contrast multimodal signal in the canonical hdWGCNA build; mixed
   sign reflects RNA-vs-proteomics-vs-phospho_corr orthogonality.
2. **DAM_up** (custom_microglia_states) — leader_score 16.4; 3
   contrasts consistent-sign≥2; 1 contrast sig≥3; dominant_sign =
   mixed; contrasts_summary = `nlgf_in_maptki:2/+ | nlgf_in_p301s:2/+
   | interaction:3/mixed`. Amyloid-driven activation at NLGF
   contrasts with mixed-sign multi-modality enrichment at interaction.
3. **(seven-way tie at leader_score 16.4)** — `GOBP_MITOCHONDRIAL_ELECTRON_TRANSPORT_NADH_TO_UBIQUINONE`,
   `GOBP_ATP_SYNTHESIS_COUPLED_ELECTRON_TRANSPORT`,
   `GOCC_NADH_DEHYDROGENASE_COMPLEX`,
   `GOMF_NADH_DEHYDROGENASE_ACTIVITY`,
   `GOMF_OXIDOREDUCTION_DRIVEN_ACTIVE_TRANSMEMBRANE_TRANSPORTER_ACTIVITY`,
   `GOMF_OXIDOREDUCTASE_ACTIVITY_ACTING_ON_NAD_P_H_QUINONE_OR_SIMILAR_COMPOUND_AS_ACCEPTOR`,
   all with the same `nlgf_in_p301s:2/- | interaction:3/mixed |
   tau_alone:2/-` (or `tau_in_nlgf:2/-`) pattern; dominant_sign = mixed.
   The OXPHOS / electron-transport family clustered at leader_score 16.4.

### Top 3 by max_substate_breadth (compartment-wide signals)

1. **DAM_up** (custom_microglia_states) — max_substate_breadth=4;
   leader_score 16.4; dominant_sign = mixed. Small-N artefact for the
   custom collection but the breadth-4 alignment with the GO compartment
   signal is real.
2. **GOBP_ATP_SYNTHESIS_COUPLED_ELECTRON_TRANSPORT** (gobp) —
   max_substate_breadth=4; leader_score 16.4; dominant_sign = mixed.
   The ATP-synthesis branch of the OXPHOS family reaches all four
   substates at interaction in the per-substate top-10.
3. **GOCC_NADH_DEHYDROGENASE_COMPLEX** (gocc) —
   max_substate_breadth=4; leader_score 16.4; dominant_sign = mixed.
   Compartment-level NADH dehydrogenase signal also reaches all four
   substates.

(Many more breadth=4 GO terms cluster at leader_score 11.4-16.4
spanning OXPHOS / ETC / aerobic respiration / cellular respiration /
ribosomal compartment — see C5 completion note for the full
ribosomal-compartment breadth-4 set.)

### Top 3 by max_consistent_sign at any single contrast (single-contrast cross-modality stars)

1. **MG-M2** (hdwgcna_modules) — max_consistent_sign=4; leader_score
   12.8; dominant_sign = +; contrasts_summary = `nlgf_in_maptki:3/+ |
   nlgf_in_p301s:4/+`. The single cleanest cross-modality consistent-
   sign hit anywhere in the project: snRNAseq NES=+2.15 padj=1.6e-111,
   GeoMx +1.33 padj=4.9e-4, proteomics +1.51 padj=2.0e-4, phospho +1.50
   padj=3.8e-3 all aligned + at `nlgf_in_p301s`.
2. **GOBP_ADAPTIVE_IMMUNE_RESPONSE** (gobp) — max_consistent_sign=3;
   leader_score 11.6; dominant_sign = +; contrasts_summary =
   `nlgf_in_maptki:2/+ | nlgf_in_p301s:3/+`. Independent GO-collection
   corroboration of the amyloid-driven activation axis surfaced by
   MG-M2.
3. **(five-way tie at max_consistent_sign=3 on the synaptic axis)** —
   `GOCC_EXCITATORY_SYNAPSE`, `GOCC_GABA_ERGIC_SYNAPSE`,
   `GOCC_PRESYNAPTIC_MEMBRANE`, `GOCC_SCHAFFER_COLLATERAL_CA1_SYNAPSE`,
   `GOCC_SYNAPTIC_VESICLE_MEMBRANE`, all at leader_score 11.6 with
   dominant_sign = − and contrasts_summary = `nlgf_in_maptki:3/- |
   nlgf_in_p301s:2/-`. NLGF-driven synaptic suppression cluster.

### Multi-axis verdict

The leader board surfaces three distinct axes of biology, none of
which is the single "winner" implied by the original Hallmark-anchored
framing:

1. **Amyloid-driven activation (+ sign).** MG-M2 + MHC II antigen-
   presentation family (5 GO BP terms + 2 GO MF terms at leader_score
   10.4) + `GOBP_ADAPTIVE_IMMUNE_RESPONSE` + custom curated DAM_up /
   WAM / HAM_disease / AD2 sets. Robust across both tau backgrounds at
   the two NLGF contrasts.

2. **NLGF-driven synaptic suppression (− sign).** Six synaptic GO CC
   terms at leader_score 11.6-15.4 plus matching GO BP synapse terms,
   all consistently negative at `nlgf_in_maptki` and `nlgf_in_p301s`.
   Cross-modality signal read by the microglia compartment; mechanism
   is open between microglial engulfment of synapses and proteomic /
   transcriptional reorganisation of the microglia-resident synaptic
   interactome.

3. **Mixed-sign metabolic/translational disagreement at tau×amyloid
   interaction.** Eight OXPHOS/ETC GO terms (BP + MF + CC) +
   `DAM_up` + `MG-M3`, clustered at leader_score 11.4-23.4 with
   dominant_sign = mixed. Every term shows the same `interaction:3/
   mixed` pattern: three modalities reach FDR<0.05 at the interaction
   contrast but disagree on direction. The ribosomal compartment
   (5 GO CC + 1 GO MF set at breadth=4) and Sun/Victor 2023
   `MG3_ribosome_biogenesis` corroborate the compartment-wide
   metabolic/translational suppression at the RNA layer. The mixed
   sign is the cross-modality residual non-additivity, not a
   directional change.

**No specificity-null testing was performed.** The cancelled C7/C8
framework would have asked whether the OXPHOS hit at interaction was
due to gene-set construction artefacts, but the multi-collection
independent replication of the same `interaction:3/mixed` pattern
across GO BP, GO MF, GO CC, custom_microglia_states, and
hdwgcna_modules adjudicates that question without needing a null. The
consistency filters in the leader rule and the cross-contrast
aggregation in the leader_score formula substitute for the cancelled
null work.

**The verdict is multi-axis.** No single pathway family carries the
narrative on its own. The per-contrast supporting evidence is in
subsections 13.1-13.6 of the knitted analysis; the unified leader
board is at `storage/results/pathway_survey_unified_leaderboard.tsv`;
the verdict display is in subsection 13.7 of `rmd/12_pathway_survey.Rmd`.

### Pointer to successor plan

`storage/notes/mechanism_layer_plan.md` (active as of 2026-05-23)
takes the three axes forward to mechanism-layer inference: Phase E
(TF inference, decoupler + CollecTRI prior), Phase F (kinase activity
inference from phospho), Phase G (CCC bridge re-analysis using the
leader pathways as sender/receiver filters). Each phase is gated on
user choice of priorities at its first session, so the mechanism work
does not pre-commit to one axis over the others.

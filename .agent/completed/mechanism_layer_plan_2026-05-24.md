# Active plan: mechanism-layer inference on the leader pathways

- **Created:** 2026-05-23 (immediate successor to `pathway_overhaul_plan.md`).
- **Last updated:** 2026-05-24 (Phase G5 DONE; **plan COMPLETE**;
  Phase G COMPLETE; the full mechanism arc Phases E + F + G now
  read all three D2 axes with at least one mechanism-layer verdict
  each. The G5 session appended subsection 16.3 ("Verdict (cross-
  tool + cross-contrast LR synthesis)") to `rmd/15_ccc_mechanism.Rmd`
  with five level-3 sub-subsections (16.3.1 .. 16.3.5) plus a new
  `build_lr_verdict_table()` helper at `R/ccc_inference.R` and an
  exported `storage/results/ccc_lr_verdict.tsv` (three rows, eight
  columns; mirror of F5 schema with `top_lr_pairs` carrying
  `<lr_interaction> (<sender>-><receiver>)` semicolon-separated
  per the LR-layer four-key requirement). The verdict reads the
  three D2 axes as STRONG / REFINED-POSITIVE / STRONG, with axis 2
  specifically refined from "honest non-finding" at the TF / kinase
  layers to "TREM2-mediated clearance + APP-fragment uptake +
  synaptic-adhesion modulation rather than classical complement-
  mediated pruning" at the LR layer -- the G1 triangulation
  hypothesis vindicated and refined. `n_top_sig_in_axis_contrasts =
  0/0/0` at top-5 is a real structural property of the G3-vs-G4
  complementarity (G4 deviation register entry #5) and is preserved
  honestly in the verdict TSV + 16.3.5 cross-axis subsection. Plan
  moved to `storage/notes/completed/mechanism_layer_plan_2026-05-24.md`
  with a full ## Outcome summary section per the plan's terminal
  instruction.).
- **State:** **plan COMPLETE.** Phase E COMPLETE (E1 + E2 + E3 + E4
  + E5 all DONE); Phase F COMPLETE (F1 + F2 + F3 + F4 + F5 all
  DONE); Phase G COMPLETE (G1 + G2 + G3 + G4 + G5 all DONE). The
  mechanism arc (Phases E + F + G) has produced three complementary
  per-axis verdicts (TF at 14.3, kinase at 15.3, LR at 16.3) that
  jointly read each of the section-13 D2 axes with mechanism-layer
  evidence. Plan moved to
  `storage/notes/completed/mechanism_layer_plan_2026-05-24.md`.
- **Predecessor:** `storage/notes/completed/pathway_overhaul_plan_2026-05-23.md`.
- **Goal:** take the unbiased multi-axis verdict from `pathway_overhaul_plan`
  D2 forward to mechanism-layer inference, anchored on the leader-board
  pathways rather than on any single pre-committed pathway family.

## Why this plan exists

The agnostic pathway survey (subsections 13.1-13.7 of
`analysis.html`) produced a leader board with 144 leader rows across
six gene-set collections. The leader board surfaces three independent
axes of biology (recorded in the predecessor plan's `## Outcome
summary`):

1. **Amyloid-driven activation (+ sign):** MG-M2 (hdwgcna_modules) +
   MHC II antigen-presentation family + `GOBP_ADAPTIVE_IMMUNE_RESPONSE`
   + custom curated DAM_up / WAM / HAM_disease / AD2 sets.
2. **NLGF-driven synaptic suppression (− sign):** six synaptic GO CC
   terms + matching GO BP synapse terms, consistently negative at the
   two NLGF contrasts.
3. **Mixed-sign metabolic/translational disagreement at tau×amyloid
   interaction:** eight OXPHOS/ETC GO terms (BP + MF + CC) + `DAM_up`
   + `MG-M3`, all with `interaction:3/mixed` cross-modality pattern;
   the ribosomal compartment + Sun/Victor 2023 `MG3_ribosome_biogenesis`
   corroborate at breadth=4.

The survey produced a verdict (multi-axis biology) but does not
produce a *mechanism*. The next layer of inference — which transcription
factors drive the activation axis, which kinases mediate the
phospho-layer disagreement at interaction, which cell-cell
communication edges connect microglia to the synaptic axis — is the
purpose of this plan. The plan is deliberately structured so the
mechanism work derives its targets from the D2 leader board
(`storage/results/pathway_survey_unified_leaderboard.tsv`) and not
from re-anchoring on Hallmark or on the cancelled OXPHOS-confirmatory
arc.

## Locked decisions (carried forward from `pathway_overhaul_plan`, plus E1)

| Decision | Choice |
|---|---|
| Privileged pathway family | None — mechanism work derives targets from D2 leader board, equal weight to the three axes unless user decides otherwise at the Phase E/F/G decision gates. |
| Verdict null testing | Drop — not in scope for any of E/F/G; the leader-rule consistency filters and cross-contrast aggregation substitute. |
| Reintroduce Hallmark | Never (permanent retirement from project). |
| Reintroduce OXPHOS-arc (08, 09a, 09b, 09c, 09d) | Never (deleted; targets must come from leader board). |
| Phase priority order | **E → F → G** (locked at E1, 2026-05-23). TF inference first (broadest, cross-modality, reuses DE caches); kinase second (narrower, phospho-only); CCC last (benefits from TF-derived gene universes). |
| TF inference tool (Phase E) | **`decoupleR` + CollecTRI mouse prior** (locked at E1, 2026-05-23). CollecTRI extends DoRothEA with ChIP-seq + literature evidence; decoupleR univariate methods operate directly on DE-stat vectors; same framework reused in Phase F (OmniPath kinase prior) for cross-comparison; mouse orthologue mapping built into `decoupleR::get_collectri(organism = "mouse")`. |
| Kinase activity inference tool (Phase F) | **`decoupleR` + OmniPath kinase-substrate prior** (locked at F1, 2026-05-24). OmniPath KSN (`decoupleR::get_ksn_omnipath()`) aggregates PhosphoSitePlus + SIGNOR + multiple databases (~5-10k human KSN edges typically); identifiers map cleanly to the phospho cache's `symbol + PTM.SiteAA + PTM.SiteLocation` columns after a `paste0(symbol, "_", site_aa, site_loc)` step; mouse mapping via `nichenetr::convert_human_to_mouse_symbols` on the substrate-symbol column (regulatory phosphosites are highly conserved so residue/position translation is usually unnecessary). Same framework as Phase E means the Phase F verdict in section 15.3 will use identical conventions (`ulm` arbitrates significance at FDR<0.10, `consensus` arbitrates direction, `minsize = 5L`) and read directly against the Phase E verdict in section 14.3. Alternatives KSEA-only and decoupleR + KSEA cross-tool triangulation considered and declined: KSEA-only weakens the cross-layer interpretation; triangulation would have roughly doubled the F2 build effort without addressing the project's critical-path question at this stage. |
| Mouse orthologue strategy (priors) | **`nichenetr::convert_human_to_mouse_symbols`** (locked at F1, 2026-05-24, inherited from Phase B2/B3/B4 precedent). Applied at Phase F to the OmniPath KSN substrate-symbol column; available as fallback for any future Phase that brings in a human-only resource. Phase E used `decoupleR::get_collectri(organism = "mouse")` directly (CollecTRI ships native mouse) so the nichenetr path was not needed for Phase E. |
| Axis-restricted output shape (E4, F4, and propagated to E5) | **Hybrid: emit BOTH axis-mean and per-axis-contrast columns side-by-side** (locked at F4 decision gate, 2026-05-24, **user-directed**: "Do the hybrid. Have other analyses do the hybrid as well for consistency"). The per-(axis, source) row in every axis-restricted TSV (TF in 14.2 and kinase in 15.2) and verdict TSV (TF in 14.3 and kinase in 15.3) carries the locked axis-mean `mean_activity_in_axis_contrasts` PLUS one `score_at_<contrast>` column per unique axis_contrast across all axes (NA where the row's axis does not include the contrast). The verdict TSVs additionally carry `per_contrast_score_range` (formatted `"[min, max]"` of per-contrast cells across the top-N) and `per_contrast_summary` (per-contrast `"name:[min, max]"` semicolon-joined). The hybrid view's design intent is to expose two distinct failure modes of the axis-mean: (a) per-contrast asymmetry (e.g., the TF top 5 at axes 1+2 show maptki:[-2.09, +5.40]; p301s:[-3.06, +3.23] — maptki carries the stronger positive tail, p301s the stronger negative tail), and (b) cell-coverage divergence at multi-modality layers (rowMeans of per-contrast columns generally != axis-mean when modality coverage varies across (axis-contrast, modality) cells; the divergence is a real weighting artefact, not a bug). For axes with a single axis_contrast (axis 3 here, `interaction_metabolic`), the per-contrast column is numerically identical to the axis-mean by construction; the schema preserves it anyway for cross-axis consistency. Rationale for picking hybrid over axis-mean-only: axis-mean alone can smooth over per-contrast disagreements that drive biological reading; per-contrast-only loses the at-a-glance axis verdict. Both are kept. Affects: `R/tf_inference.R::score_tf_per_axis()`, `R/tf_inference.R::format_axis_restricted_table()`, `R/tf_inference.R::build_tf_verdict_table()`, `R/kinase_inference.R::score_kinase_per_axis()`, `R/kinase_inference.R::format_axis_restricted_kinase_table()`, plus the rmd 13/14 prose at 14.2/14.3/15.2 (and the upcoming 15.3). |
| Kinase leader rule (Phase F3) | **Per-contrast independent; corrected cache only; FDR<0.10 on `ulm`; consensus arbitrates direction** (locked at F3, 2026-05-24, **user-directed deviation** from the originally locked cross-cache + cross-contrast aggregation rule). The originally locked rule (`phospho_corrected` as primary + `phospho_raw` as corroborator; gated by `n_contrasts_consistent_sign_ge2_across_caches >= 2` OR `n_contrasts_sig_corrected >= 3`) was shown by the F3 pre-build diagnostic to be structurally unsuitable: the cross-cache consistency arm yields zero kinase-contrast cells at FDR<0.10 (the two caches even disagree on sign at the dominant Gsk3b signal in three contrasts); the corrected-cache `n_sig >= 3` arm yields exactly one kinase (Gsk3b) because the small inference universe (~100 kinases per contrast) makes BH-padj harsh; the pre-locked FDR<0.15 fallback adds 0-1 kinases and does not address the cross-cache arm. The new rule treats the five canonical contrasts as independent inference units, consults only the corrected cache, and reports per-contrast significant sets at FDR<0.10 on the ulm p-value (E3 convention preserved). The phospho_raw cache is deliberately excluded from F3/F4/F5 because it cannot corroborate phospho_corrected at FDR scale (the diagnostic established this; the corrected-vs-raw mismatch is the signature of batch correction recovering real signal that batch noise was suppressing rather than a fragility of the corrected cache). Whether the contrast-independence directive propagates to the axis-restricted layer in F4 is a separate question logged for the F4 decision gate (see F3 completion note plan-spec deviation #5 and the F4 stub). |
| CCC re-analysis approach (Phase G) | **Three-tool triangulation: CellChat (existing) + MultiNicheNet (existing) + liana (new), with leader-pathway-restricted gene-universe post-filter applied across all three** (locked at G1, 2026-05-24, **user-directed deviation** from the plan's default post-filter-only proposal). User rationale (encoded in the question's option-3 framing and the user's pick): the synaptic-suppression axis is an explicit non-finding at the TF (E5 verdict 14.3.3) and kinase (F5 verdict 15.3.3) layers, and the verdict prose specifically named CCC as the layer expected to surface engulfment-mediated biology; adding tools to read an axis that prior tools couldn't read is the natural response and is a more defensible use of triangulation than the F1 case (where the F1 single-tool result was already strong). Post-filter against leader-pathway gene universes still applies (sender + receiver L-R pairs intersected with per-axis universes from D1's leader-board TSV) and is now applied uniformly to all three tools; the cross-tool consensus rule (mirror of E3's `n_modalities_sig_consistent_sign >= 2 OR n_modalities_sig >= 3` cross-modality rule, adapted to (LR-pair × tool)) is defined in G3. Alternatives considered and declined: (a) post-filter-only on CellChat + MultiNicheNet (plan default, the cheapest option; declined because it would replicate F1's single-tool-suffices logic on an axis where prior tools failed); (b) re-run CCC from scratch with leader-restricted prior (declined because it breaks comparability to section 11 and changes a locked input — the CellChatDB.mouse + NicheNet priors are part of the project's reproducibility surface); (c) hybrid (post-filter now, escalate later) — declined because it defers the triangulation decision rather than answering it now. |
| CCC mechanism Rmd venue (Phase G) | **New `rmd/15_ccc_mechanism.Rmd` rendering as section 16** (locked at G1, 2026-05-24, per the recommended option). Mirrors the per-phase one-Rmd-per-mechanism-layer convention from `rmd/13_tf_inference.Rmd` (file index 13 → rendered section 14, TF) and `rmd/14_kinase_inference.Rmd` (file index 14 → rendered section 15, kinase). Reader sees TF (14) + kinase (15) + CCC mechanism (16) grouped as a coherent mechanism arc in the TOC; section 11 remains the descriptive CCC home (subsections a-f) and is not modified by Phase G. Alternative considered and declined: extend `rmd/11_ccc.Rmd` with leader-restricted subsections (g, h, ...) per the G1 plan stub's preference — declined because it fragments the cross-layer mechanism reading across sections 11 and 14/15. |
| Cross-tool LR leader rule (Phase G3) | **`n_tools_sig_consistent_sign >= 2 OR n_tools_sig >= 3` at per-tool FDR<0.10** (locked at G3, 2026-05-24, OR form chosen for E3 parity). Per-tool sig conventions: CellChat pval < 0.10 in the contrast's primary genotype (subsetCommunication's native pval<0.05 thresh is already implicit, so the FDR<0.10 outer cut is a binary inclusion vote -- "this LR pair has detectable signalling in the contrast's primary genotype"); MultiNicheNet BOTH `scaled_p_val_ligand_adapted < 0.10` AND `scaled_p_val_receptor_adapted < 0.10` (MNN-recommended both-ends-DE filter); LIANA+ `cellphone_padj_primary < 0.10` (BH within per-genotype run, computed at G2). The pre-build diagnostic confirmed the rule populates: 46-97 cells per contrast (vs F3's emergency floor of zero), 0-1 cells per contrast at the all-3-sig arm (essentially dead at FDR<0.10 because the three tools' significance universes intersect sparsely; the consistent-sign-≥2 arm is the actual leader populator). 194 unique leader (sender, receiver, ligand, receptor) cells in the cross-contrast leader board after aggregation. Stub's original AND form `n_tools_sig >= 2 AND n_tools_sig_consistent_sign >= 2` is numerically equivalent here because sig-with-direction implies sig; the OR form is preferred for notation parity with E3's `n_modalities_sig_consistent_sign >= 2 OR n_modalities_sig >= 3`. CellChat duplicate-key dedupe bug surfaced during smoke testing (subsetCommunication returns multiple rows for the same (s,r,l,rec) cell under different interaction_name_2 downstream-complex annotations); fix added at `.extract_lr_cellchat()` (dedupe by max prob within (s,r,l,rec) after receptor-complex explosion). Cross-tool magnitude tiebreak via per-tool per-contrast `score_norm` (raw score / max(|score|) within tool × contrast) averaged over significant tools (mirror of E3's mean_abs_score over the sig subset); per-tool raw + normalised scores retained for transparency. Affects: `R/ccc_inference.R::.extract_lr_cellchat()`, `extract_lr_per_tool()`, `rank_lr_cross_tool()`, `build_lr_cross_tool_leaderboard()`, `format_lr_ranking_table()`, `plot_lr_cross_tool_heatmap()` plus the rmd/15 prose at 16.1 / 16.1.1-16.1.4. |

## Open questions (defer to the session that needs them)

| Question | Default proposal | Decided in |
|---|---|---|
| *(no open questions; G3 completed the last open one — the LIANA+ R-vs-Python question was resolved at G2 in favour of Python LIANA+ v1.7.1 via reticulate)* | — | — |

## Execution model

- A fresh session must read `CLAUDE.md` first, then this file in full.
- A fresh session executes the NEXT step whose status is `TODO`. Steps
  are sized to fit comfortably in one session each.
- After completing a step, the session must:
  1. Update that step's status from `TODO` to `DONE <YYYY-MM-DD>` with a
     multi-paragraph completion note (what was built, key file paths,
     biology highlights, plan-spec deviations explicitly acknowledged).
  2. `chown rstudio:rstudio` any new files (agent runs as root).
  3. Re-knit `analysis.Rmd` and verify zero `class="error"` and zero
     `class="warning"` in the rendered HTML.
  4. Commit locally per CLAUDE.md (imperative subject under 70 chars;
     co-author trailer).
  5. Either continue to the next `TODO` step in the same session if the
     context budget permits, or end cleanly.
- Decision-gate steps (E1 phase priority + each first session per
  phase) require user confirmation via `AskUserQuestion`. Wait for
  confirmation before proceeding.

## Phase E: transcription factor inference (sessions E1-E5)

### Session E1: phase-priority + TF-inference-tool decision gate [DONE 2026-05-23]

**Completion note.** Decision gate session. Confirmed Phase priority
order **E → F → G** with the user via `AskUserQuestion` (default
proposal: TF inference first because it is the broadest mechanism
layer, runs cross-modality on snRNAseq + GeoMx + proteomics DE caches
without new data ingestion, and speaks to both the amyloid-activation
and the synaptic-suppression axes; kinase second because narrower /
phospho-only; CCC last because it benefits from TF-derived gene
universes). Confirmed TF inference tool **`decoupleR` + CollecTRI
mouse prior** with the user via `AskUserQuestion` (default proposal:
CollecTRI extends DoRothEA with ChIP-seq + literature evidence;
`decoupleR` ships univariate methods — `ulm`, `mlm`, `wsum`, `viper`,
`consensus` — that operate directly on the DE-stat vectors we already
have; the same framework will be reused in Phase F with the OmniPath
kinase prior for cross-tool / cross-layer comparison). Alternatives
considered and rejected: F-first (smaller / noisier OmniPath prior
would yield weaker findings to anchor on); G-first (CCC without
TF-derived gene-universe filters is a re-run of existing work);
DoRothEA (now a subset of CollecTRI); pyscenic (heavy *de novo*
regulon discovery overkill given we have curated priors + DE stats).

Recorded both decisions in the Locked decisions table and removed the
two corresponding rows from the Open questions table. Scoped Phase E
into four atomic sub-sessions E2..E5 below, each one fitting in a
single context window. Phase F (sessions F1-F?) and Phase G (sessions
G1-G?) remain stubbed at their decision-gate steps and will be scoped
in their own decision-gate sessions once Phase E is complete, per the
plan's existing structure.

Housekeeping: corrected a typo in the original E1 stub where the
phospho-corrected DE cache was referenced as `de_phospho_corr.rds`;
the actual file on disk is `storage/cache/de_phospho_corrected.rds`.
The corrected name is used in the E2 Inputs section below. No new R
or Rmd code added this session — planning-only. Verification: this
session changed only `storage/notes/mechanism_layer_plan.md` (a
non-Rmd file outside the knit graph), so re-knit was a no-op
verification; ran it anyway per the execution-model discipline and
analysis.html was unchanged (zero `class="error"`, zero
`class="warning"`).

---

### Session E2: build decoupleR + CollecTRI TF activity cache [DONE 2026-05-23]

**Completion note.** Built the Phase E mechanism cache cleanly,
outside the knit, per the plan. `decoupleR` 2.16.0 was already
installed in the container; CollecTRI mouse prior pulled via
`decoupleR::get_collectri(organism = "mouse", split_complexes = FALSE)`
(under the hood OmnipathR cache returned 40,291 interactions ->
39,961 after deduplication; 1,114 unique TF sources, 6,110 targets).

Cache shape `modality -> contrast -> tibble(statistic, source, score,
p_value)` written to `storage/cache/tf_activity_decoupler.rds`
(712 KB xz-compressed). Seven modalities x five contrasts = 35
nested tibbles, no gaps:

```
snrnaseq               (whole microglia, 11,411 syms) ->  482 TFs / contrast
snrnaseq_homeostatic   (substate,         5,018 syms) ->  313 TFs / contrast
snrnaseq_DAM           (substate,         5,144 syms) ->  324 TFs / contrast
snrnaseq_IFN           (substate,         5,560 syms) ->  332 TFs / contrast
snrnaseq_proliferative (substate,         4,913 syms) ->  331 TFs / contrast
geomx                  (limma,           19,959 syms) ->  679 TFs / contrast
proteomics             (limma,            3,460 syms) -> ~252 TFs / contrast
```

The five-contrast set is uniform: `nlgf_in_maptki`, `nlgf_in_p301s`,
`interaction`, `tau_alone`, `tau_in_nlgf`. Each (modality, contrast)
tibble carries five `statistic` levels per TF: `ulm`, `wsum`,
`norm_wsum`, `corr_wsum`, `consensus`. Downstream sessions should
rank on `consensus` (decoupleR's recommended primary score) with the
per-method scores retained for sensitivity checks.

**Plan-spec deviations explicitly acknowledged.**

1. The plan's input filename `de_snrnaseq_nebula_microglia_compartment.rds`
   does not exist on disk. The canonical per-substate cache is
   `storage/cache/de_snrnaseq_nebula_per_state.rds` (produced by
   `scripts/build_per_state_nebula.R`); a second variant
   `de_snrnaseq_nebula_per_state_1pct.rds` exists as a 1%-expression-filter
   robustness check consumed only by `rmd/02e_snrnaseq_substate_pathway.Rmd`.
   E2 uses the canonical (no `_1pct`) cache. The plan's filename
   reference was stale; downstream E3..E5 sessions should refer to the
   per-substate file by its actual name.
2. The plan's sanity bound "n_TFs per tibble is ~1100 (CollecTRI mouse)"
   was overly optimistic. The total CollecTRI mouse TF count is 1,114,
   but after applying `minsize = 5L` (the standard decoupleR / fgsea
   universe-overlap threshold), only the TFs with at least five targets
   present in the modality's gene universe survive. The actual range
   is 251 (proteomics) to 679 (GeoMx), tracking universe size as
   expected. This is the correct behaviour: a TF with two CollecTRI
   targets in a 3.5k-gene proteomics universe cannot be inferred
   reliably and is rightly excluded. E3 should not interpret the
   smaller per-modality TF counts as a bug.
3. The plan suggested running `decoupleR::run_ulm` (fast default) and
   `run_consensus` if compute permits. I went one step further and
   used `decoupleR::decouple()` with `statistics = c("ulm", "wsum")`
   plus `consensus_score = TRUE` for every modality. Compute permits
   easily: wall-clock total was 153.68s (~2.5 min) for all seven
   modalities x five contrasts. The wsum family auto-emits three
   sub-statistics (raw, normalised, correlation-adjusted) plus the
   cross-method consensus. This buys downstream sessions a richer
   sensitivity check at negligible cost.
4. `run_decoupler_per_modality()` iterates per contrast internally
   rather than passing the full contrast matrix in one call. Reason:
   `decoupleR::decouple()` refuses NA / Inf cells in `mat`, but
   limma-based caches (proteomics in particular) have asymmetric
   per-contrast NAs from small-group estimability quirks
   (75 NA cells across 25 protein groups in proteomics; rows present
   in only 2 of 5 contrasts). Per-contrast iteration lets each call
   see its full clean stat vector rather than dropping any gene
   missing in any contrast. Output shape is identical (the
   `condition` column carries the contrast name). The trade-off is
   a small per-modality runtime increase (~2-3x on the snRNAseq
   substates) due to repeated method-startup cost; total still well
   under 3 minutes.

**File outputs.** All chowned `rstudio:rstudio`:

* `R/tf_inference.R` (new, 12,159 bytes) -- four helpers documented
  in module-level comment block + inline argument docs:
  `extract_de_stat_vec()`, `extract_de_stat_matrix()`,
  `run_decoupler_per_modality()`, `split_decoupler_by_contrast()`.
  Symbol-keyed throughout; duplicate symbols collapsed by max |stat|
  (rationale documented in the module header).
* `R/helpers.R` (edited) -- new `source("R/tf_inference.R")` line
  inserted after `R/pathway_survey.R` per the dependency-order
  convention; the source-order comment block updated with a
  rationale paragraph for the new file's position.
* `scripts/build_tf_activity_decoupler.R` (new, 6,098 bytes) --
  heavy compute outside the knit, idempotent via `--overwrite`,
  detailed header comment covering inputs / outputs / runtime /
  cache shape.
* `storage/cache/tf_activity_decoupler.rds` (new, 712 KB) -- the
  cache itself.

**Verification.** Per the plan, no knit was required at E2 because
no Rmd consumer exists yet. Pre-build smoke-test of the helpers on
a single (geomx, all contrasts) call passed: matrix extraction
shape correct, decoupler ran, split-by-contrast preserved the
five-contrast keys. Post-build smoke-test on the persisted cache
confirmed: 35 nested tibbles all present; TF / contrast / statistic
counts as expected; no all-NA tibbles; example consensus rows from
`snrnaseq$interaction` look reasonable (mix of signed scores, p
values in [0, 1]). The cache contains a small number of UniProt-
accession TF names (e.g. `A0A087WPA7`) inherited from CollecTRI's
complex sources; E3 should decide whether to filter to symbol-only
TFs for the leaderboard or keep them for completeness (recommend
keep, with a flag column added at the leaderboard step).

---

### Session E3: cross-modality TF activity leaderboard + first 13_tf_inference subsection [DONE 2026-05-23]

**Completion note.** Built the cross-modality TF activity leader
board and the first subsection of the new section 14 in
`analysis.html` (the file is numbered 13 in the project namespace
`rmd/13_tf_inference.Rmd`, but renders as section 14 in the knit
because the parent inserts it after the section-13 pathway survey).
Smoke testing on the live cache passed; knit completed in 2 min
59 s with zero `class="error"` and zero `class="warning"` in
`analysis.html`. The new section comprises:

* Section 14 "TF activity inference" header + intro prose
  cross-referencing the D2 three-axis verdict from section 13.
* Subsection 14.1 "Cross-modality TF activity ranking" with four
  level-3 sub-subsections:
  * 14.1.1 Per-contrast cross-modality top-10 TFs (5 tables).
  * 14.1.2 Unified TF leader board (single table, 10 rows).
  * 14.1.3 TF × modality activity heatmap at `nlgf_in_maptki`
    (ComplexHeatmap; rows alphabetical; `*` = FDR<0.05, `.` =
    FDR<0.10; NA cells grey).
  * 14.1.4 Interpretive notes (biology face-validity + the
    deliberate absence of synaptic / interaction axes at this
    layer, to be revisited in E4/E5/Phase F).

**Leader board content (full).** 10 TFs at the rule `(>=2 modalities
sig consistent OR all 3 sig)` applied at FDR<0.10 across the 3
primary modalities (snrnaseq whole-microglia + geomx + proteomics):

```
1.  Nr1h3   leader_score 10.4  nlgf_in_maptki:2/+ | nlgf_in_p301s:2/+
2.  Irf1    leader_score 10.4  nlgf_in_maptki:2/+ | nlgf_in_p301s:2/+
3.  Usf2    leader_score 10.4  nlgf_in_maptki:2/+ | nlgf_in_p301s:2/+   substate_breadth=4
4.  A0A979HLR9  leader_score 5.4   nlgf_in_maptki:2/+   substate_breadth=1
5.  Spi1    leader_score 5.4   nlgf_in_maptki:2/+
6.  A0A087WSP5  leader_score 5.4   nlgf_in_maptki:2/+   substate_breadth=2
7.  Foxa2   leader_score 5.4   nlgf_in_maptki:2/+   substate_breadth=1
8.  Dbp     leader_score 5.4   nlgf_in_maptki:2/+
9.  Gata4   leader_score 5.4   nlgf_in_maptki:2/+   substate_breadth=2
10. Tbx21   leader_score 5.4   nlgf_in_maptki:2/+
```

Every leader has `dominant_sign = "+"` (activation). All leader
contrasts are NLGF contrasts; zero leader contrasts at `interaction`,
`tau_alone`, or `tau_in_nlgf`. Biology highlights: Nr1h3 (LXRα,
microglial cholesterol/DAM priming), Spi1 (PU.1, microglial master
TF + AD GWAS locus), Irf1 (interferon-driven activation), Usf2 with
substate_breadth = 4 (active in every microglia substate — the
closest TF analogue of the compartment-wide pathway patterns at
section 12). The cross-modality TF inference layer therefore reads
the **amyloid-driven activation axis** of the D2 verdict; the
**synaptic suppression axis** and **mixed-sign interaction axis** do
not surface at this layer with this prior (interpretation: synaptic
suppression is most likely a detection-level engulfment phenomenon
not driven by direct transcriptional regulation, and the interaction
OXPHOS axis is post-transcriptional — to be revisited in Phase F's
kinase layer).

**Plan-spec deviations explicitly acknowledged.**

1. **Significance statistic switched from "consensus" to "ulm".**
   The plan locked `padj < 0.10` but did not specify which
   decoupleR statistic's p-value drives the BH adjustment. The
   helper module initially used "consensus" because the E2 note
   recommended ranking on consensus. Empirical smoke testing on
   the live cache showed `consensus` p-values to be too
   conservative: only 1 TF cleared the leader rule across all 5
   contrasts (because the wsum family's empirical p-values floor
   at 0.02 and pull the Stouffer-combined consensus p toward
   conservatism). Switching the sig statistic to "ulm" (the
   linear-model F-test p-value, decoupleR's canonical
   univariate statistic) while keeping "consensus" for ranking
   and sign attribution gives a meaningfully populated leader
   board with face-valid biology. This split is the decoupleR
   vignette convention (ulm arbitrates sig; consensus arbitrates
   magnitude / direction); the rationale is documented in the
   `R/tf_inference.R` module header and in the section 14.1
   intro prose. The plan-locked threshold of FDR<0.10 is
   preserved.

2. **TSV form is reduced (one row per TF), not "(one row per TF ×
   contrast)" as the plan text states.** The plan's listed columns
   (`n_contrasts_consistent_sign_ge2`, `leader_score`) are
   unambiguously per-TF cross-contrast aggregates, identical to
   the columns of `storage/results/pathway_survey_unified_leaderboard.tsv`
   (reduced form, one row per (collection, pathway)). The "(one
   row per TF × contrast)" parenthetical in the plan appears to
   be a slip. The exported TSV mirrors the pathway leader board
   exactly: 11 columns, one row per leader TF, with the per-
   contrast detail packed into `contrasts_summary` as a
   pipe-delimited string. Long-form per-(TF, contrast) ranking
   data lives in the in-memory object `tf_ranking_long` built
   inside the knit; E4/E5 can recompute it from the cache
   without re-exporting.

3. **Helper docstring correction (carried over from E2).** The E2
   stub's claim that `consensus` p_value is NA was incorrect — the
   current decoupleR 2.16.0 populates it with a Stouffer-combined
   continuous p-value. Empirical confirmation: consensus p_values
   in the live cache range from ~2e-14 to ~1.0, with no NAs. The
   module header now records the corrected behaviour and the
   reason behind preferring ulm for sig despite consensus being
   populated.

4. **Heatmap contrast is `nlgf_in_maptki`, not `interaction`.**
   The plan said "TF × modality agreement heatmap" without
   specifying a contrast; the natural parallel from section 12
   would have been the interaction contrast. But interaction has
   only 1 cross-modality leader TF at FDR<0.10; the heatmap would
   be 1 row and uninformative. `nlgf_in_maptki` has 10 leader TFs
   and the cleanest cross-modality activity pattern in the data.
   The heatmap subsection title explicitly states the contrast
   used.

**File outputs.** All chowned `rstudio:rstudio`:

* `R/tf_inference.R` (extended; 38 KB / 1.0K lines) — appended
  Phase E3 helpers `.extract_tf_per_modality()` (internal),
  `rank_tfs_cross_modality()`, `compute_tf_substate_breadth()`,
  `build_tf_leader_board()`, `format_tf_ranking_table()` plus a
  ~50-line module-level comment block stating the statistic
  conventions, primary-vs-substate modality split, leader rule,
  and composite score formula. Existing E2 helpers untouched.
  Two prior docstring claims about `consensus` p_value being NA
  corrected.
* `rmd/13_tf_inference.Rmd` (new; 11 KB) — section 14 + subsection
  14.1 content. All `knitr::kable()` calls are in `results = 'asis'`
  chunks per CLAUDE.md. The build chunk reads the E2 cache and
  writes the leader-board TSV; downstream sessions can `readr::
  read_tsv()` the TSV without recomputing.
* `analysis.Rmd` (edited; +3 lines) — `child-tf-inference` chunk
  inserted between `child-pathway-survey` and `child-session`.
* `storage/results/tf_activity_unified_leaderboard.tsv` (new;
  884 B) — 10-row leader board with the same 11 columns as the
  pathway leader board, TF-keyed (`source` instead of
  `collection, pathway`).
* `analysis.html` (regenerated; 26 MB) — section 14 renders
  cleanly with zero errors / warnings. Section / subsection IDs
  `tf-activity-inference`, `cross-modality-tf-activity-ranking`,
  `per-contrast-cross-modality-top-10-tfs`,
  `unified-tf-leader-board`,
  `tf-x-modality-activity-heatmap-nlgf_in_maptki-contrast`,
  `interpretive-notes-for-this-subsection` all present in the
  rendered TOC at the expected numbers (14, 14.1, 14.1.1 .. 14.1.4).

**Verification.** Pre-knit smoke test via `Rscript -e ...` against
the live cache (~2 s) — confirmed `rank_tfs_cross_modality` returns
the expected schema, `build_tf_leader_board` produces a 10-row
output, `compute_tf_substate_breadth` returns a long tibble with
breadth distribution 0/1/2/3/4 = 1702/6/3/2/2 rows. Post-knit
verification: `grep -c 'class="error"' analysis.html` returns 0;
`grep -c 'class="warning"' analysis.html` returns 0; the TF section
heading is rendered as section 14 with 4 level-3 subsections at
14.1.1 .. 14.1.4; leader TF names appear in both the rendered
leaderboard table and the heatmap row labels; the TSV file is
on disk at the expected path and shape.

---

### Session E4: leader-pathway-restricted TF analysis [DONE 2026-05-23]

**Completion note.** Built the axis-restricted TF activity ranking
layer on top of E3's cross-modality leader board. The three D2 axes
(amyloid_activation; synaptic_suppression; interaction_metabolic)
are now scored by their own per-axis gene universes drawn from the
leader-pathway TSV, with the CollecTRI TF--target graph restricted
per axis. The new subsection 14.2 in `analysis.html` shows one
table + one lollipop chart per axis, plus an interpretive notes
sub-subsection (14.2.4) that addresses cross-axis structure and
resolves the forward-looking question from 14.1.4 about CollecTRI
complex sources (decision: **retain** as compound TF entities).

**Axis universe sizes + TF counts (live values).**

```
axis                       universe_genes  axis_contrasts             TFs_scored
amyloid_activation                  3395   {nlgf_in_maptki, p301s}           437
synaptic_suppression                3940   {nlgf_in_maptki, p301s}           419
interaction_metabolic               1859   {interaction}                     262
```

All three axes pass the plan's `>= 5 TF` sanity bound by a wide
margin. The plan-spec is satisfied: helpers `build_axis_gene_universe()`,
`restrict_collectri_to_universe()`, `score_tf_per_axis()` added to
`R/tf_inference.R`; subsection 14.2 added to `rmd/13_tf_inference.Rmd`
with one table + one lollipop per axis; long-format TSV at
`storage/results/tf_activity_axis_restricted.tsv` (1118 rows).

**Top TFs per axis (ranks 1..5, |mean| desc).**

```
amyloid_activation:    A0A979HLR9 (+4.32), Spi1 (+2.81), Rreb1 (-2.57),
                       Sp3 (+2.50), Nfkb1 (+2.47).
synaptic_suppression:  A0A979HLR9 (+4.32), Spi1 (+2.81), Rreb1 (-2.57),
                       Sp3 (+2.50), Nfkb1 (+2.47).
interaction_metabolic: A0A979HLR9 (-4.85), Myc (-3.25), Creb1 (-2.28),
                       Tp53 (-2.15), Jun (-1.81).
```

Biology highlight: **Rela** (NF-kB p65) appears in both axis 1
(rank 7, +2.19) and axis 3 (rank 13, -1.36) with opposite signs --
amyloid contrasts activate Rela targets, the `tau x amyloid`
interaction suppresses them. The second CollecTRI complex source
`A0A087WSP5` mirrors the same pattern (axis 1 +2.42; axis 3 -1.59).
This sign reversal is not a bookkeeping artefact -- it is the
expected signature of a coherent inflammatory -> metabolic-shift
transition between contrast settings, and supports the axis-3
interpretation as a compartment shift away from the high-
translation / high-OXPHOS state.

**Plan-spec deviations explicitly acknowledged.**

1. **hdwgcna_modules `gene_name` needed ENSMUSG -> symbol
   translation.** The `hdwgcna_microglia.rds$modules` slot uses
   Ensembl IDs (Seurat rownames) while CollecTRI is symbol-keyed.
   The Rmd build chunk translates via `snrnaseq_symbol_map.rds`
   before passing the modules data.frame to
   `build_axis_gene_universe()`, and drops the `grey` (unassigned)
   module rows. Without this step, hdwgcna would have contributed
   zero genes to all axis universes and the bookkeeping would have
   been silently wrong. Documented inline in the build-chunk comments.

2. **Interpretation A locked over Interpretation B.** The plan
   said "rank TFs within each axis by mean TF activity score
   across the axis-relevant contrasts and across modalities". I
   read this as **Interpretation A** (reuse the existing cache
   scores from E2, restrict the TF *set* by universe overlap,
   average across axis cells). **Interpretation B** (re-run
   decoupleR with each modality's CollecTRI network clipped to the
   axis universe) was considered and rejected for the statistical-
   power reasons documented in the `R/tf_inference.R` E4 module
   header. The locked choice is recorded in both the module header
   and the 14.2 prose.

3. **Structural property of Interpretation A: axes 1 and 2 share
   their top TFs by mean score.** Because both axes'
   `axis_contrasts` evaluate to the same set (`{nlgf_in_maptki,
   nlgf_in_p301s}`), and the per-TF mean uses cached scores computed
   against full per-modality universes (not the axis-restricted
   networks), the per-TF mean is identical between axes 1 and 2
   for any TF passing both universe-overlap filters. The
   differentiator is the universe-overlap set (different TFs
   survive each filter) and the `n_targets_in_axis_universe`
   value -- the live tables show e.g. Spi1 at 184 amyloid-axis
   targets vs 84 synaptic-axis targets. Axis 3 is genuinely
   independent because its `axis_contrasts` reduces to
   `{interaction}` alone. This property is honestly documented in
   14.2.4; smoothing it out would distort the biology.

4. **Resolution of the 14.1.4 forward-looking question.** Decided
   to **retain** the two CollecTRI complex sources `A0A979HLR9`
   and `A0A087WSP5` in axis-restricted analyses rather than
   filter them out. Both surface repeatedly with high `|mean|`
   scores and substantial target counts in axes 1 and 3;
   filtering would obscure substantive biology (notably the Rela-
   parallel sign reversal noted above). A future annotation pass
   can resolve them to specific subunits via UniProt; for now
   they are reported with accession names as compound TF
   entities. The resolution is recorded in 14.2.4.

5. **Axis 3 is 14/15 negative at the top (Hdac1 the exception).**
   The signed mean across the top 15 of axis 3 leans negative
   coherently for the proliferation / cell-fate TFs
   (Myc/Creb1/Tp53/Jun/Fos/Hif1a). Hdac1 (rank 7, +1.72) is
   positive but has only 9 in-universe targets, so its position
   is more plausibly a low-denominator artefact than a real
   opposing signal. The E5 verdict should weigh this when
   producing the axis-3 sign-aware reading.

**File outputs.** All chowned `rstudio:rstudio`:

* `R/tf_inference.R` (extended; 62 KB total) -- appended five
  Phase E4 helpers: `build_axis_gene_universe()`,
  `restrict_collectri_to_universe()`, `score_tf_per_axis()`,
  `format_axis_restricted_table()`, `plot_axis_lollipop()`, plus a
  multi-paragraph module header block documenting the
  Interpretation A/B choice, the axis classification rule
  (strict-directional regex matched against `contrasts_summary`),
  and the strict-directional rationale. E2/E3 helpers untouched.
* `rmd/13_tf_inference.Rmd` (extended; 26 KB total) -- appended
  subsection 14.2 ("Axis-restricted TF rankings") with three
  level-3 sub-subsections (14.2.1 amyloid_activation,
  14.2.2 synaptic_suppression, 14.2.3 interaction_metabolic) plus
  the 14.2.4 interpretive notes. The build chunk reads the
  unified leader-board TSV + six gene-set caches + the E2 TF
  activity cache + CollecTRI mouse prior, then calls the helpers
  in pipeline order. Each axis sub-subsection has a
  `format_axis_restricted_table()` chunk (in `results = 'asis'`)
  and a `plot_axis_lollipop()` figure chunk.
* `storage/results/tf_activity_axis_restricted.tsv` (new;
  82 KB) -- long-format: `axis`, `source`,
  `mean_activity_in_axis_contrasts`, `sd_activity_in_axis_contrasts`,
  `n_cells_used`, `n_targets_in_axis_universe`, `leader_rank`.
  1118 rows (one per `(axis, TF)` pair across the three axes'
  surviving TF sets).
* `analysis.html` (regenerated; 28 MB) -- subsections 14.2.1 ..
  14.2.4 render cleanly with zero errors / warnings. Anchor IDs
  `axis-restricted-tf-rankings`, `amyloid-activation-axis-axis-1`,
  `synaptic-suppression-axis-axis-2`,
  `interaction-metabolictranslational-axis-axis-3`, and
  `interpretive-notes-for-this-subsection-1` (auto-disambiguated
  from the 14.1.4 anchor) all present in the rendered TOC at the
  expected numbers.

**Verification.** Pre-knit smoke test on the live cache (~5 s)
confirmed: hdwgcna translation produces 2274 valid rows across
MG-M1..MG-M4 modules (grey dropped); `build_axis_gene_universe()`
returns three non-empty axes with the universe sizes / contrast
sets above; `restrict_collectri_to_universe()` at
`min_targets = 5` leaves 437 / 419 / 262 TFs per axis; the
score table is 1118 rows long. Post-knit verification:
`grep -c 'class="error"' analysis.html` = 0;
`grep -c 'class="warning"' analysis.html` = 0; TSV file on disk
at the expected path and shape; section TOC shows 14.2 + four
level-3 sub-subsections at expected numbers.

**Amendment 2026-05-24 (hybrid retrofit, propagated from the F4
decision gate).** At the F4 decision gate the user directed "Do the
hybrid. Have other analyses do the hybrid as well for consistency",
which propagated to E4. `R/tf_inference.R::score_tf_per_axis()` was
extended to emit `score_at_<contrast>` columns (one per unique
axis_contrast across axes), and `format_axis_restricted_table()` was
extended to auto-include axis-relevant per-contrast columns in the
kable. The axis-mean column `mean_activity_in_axis_contrasts` is
preserved byte-for-byte; row-mean of per-contrast columns equals
axis-mean to 3.08e-02 (an honest weighting artefact when cell coverage
varies across (modality, contrast) cells; the axis-mean is mean of
all finite cells, the row-mean of per-contrast columns is a mean of
cross-modality means and therefore weights each contrast equally).
The TSV `storage/results/tf_activity_axis_restricted.tsv` was
rewritten with 3 additional `score_at_<contrast>` columns. See the F4
completion note (below) for the full E4 + E5 + F4 retrofit history.

---

### Session E5: TF inference verdict subsection [DONE 2026-05-24]

**Completion note.** Wrote subsection 14.3 ("Verdict (sign-aware TF
synthesis)") in `rmd/13_tf_inference.Rmd`, mirroring the shape of the
section-13 D2 verdict at the mechanism layer. The subsection comprises
five level-3 sub-subsections at the rendered TOC:

* 14.3 Verdict (sign-aware TF synthesis) -- top-line framing +
  explicit re-statement of all three inference thresholds (per-modality
  FDR<0.10 on `ulm`; cross-modality leader rule `n_modalities_sig_consistent_sign >= 2`
  OR `n_modalities_sig >= 3`; axis-restricted `min_targets = 5`,
  top-5-per-axis displayed).
* 14.3.1 Verdict summary table -- the one-row-per-axis kable view of
  the verdict TSV.
* 14.3.2 Amyloid-activation axis (axis 1) -- corroboration call.
* 14.3.3 Synaptic-suppression axis (axis 2) -- explicit non-finding
  call.
* 14.3.4 Interaction metabolic/translational axis (axis 3) --
  corroboration call with resolution of the gene-set-layer
  mixed-sign ambiguity.
* 14.3.5 Cross-axis observations and limitations -- sign reversals,
  Interpretation-A trade-off, asymmetric three-axis verdict
  rationale.

**Verdict bottom line (per-axis, three-line summary).**

```
amyloid_activation     top 5: A0A979HLR9, Spi1, Rreb1, Sp3, Nfkb1
                       signs: +,+,-,+,+    range [-2.57, +4.32]
                       2 of 5 also cross-modality leaders (A0A979HLR9, Spi1)
                       => CORROBORATES section-13 axis 1 (MG-M2 + MHC II + DAM_up family)

synaptic_suppression   top 5: A0A979HLR9, Spi1, Rreb1, Sp3, Nfkb1
                       (identical to axis 1 by Interpretation A)
                       differentiator: in-universe target counts
                       => HONEST NON-FINDING; TF layer does not identify
                          synaptic-suppression drivers; consistent with
                          engulfment-mediated (post-transcriptional) mechanism

interaction_metabolic  top 5: A0A979HLR9, Myc, Creb1, Tp53, Jun
                       signs: -,-,-,-,-    range [-4.85, -1.81]
                       1 of 5 also cross-modality leader (A0A979HLR9, via NLGF status)
                       => CORROBORATES section-13 axis 3 (MG-M3 + DAM_up +
                          OXPHOS/ETC + ribosomal compartment); also RESOLVES
                          the section-13 gene-set-layer mixed-sign signal in
                          favour of a single biological direction (Myc and
                          partners off, biosynthesis down)
```

The asymmetric three-axis verdict (strong / non-finding / strong) is
the *real* finding of the TF layer. Phase F (kinase) and Phase G
(CCC) will revisit axis 2 with priors that can plausibly surface a
post-transcriptional or engulfment-related driver if one exists.

**Biology highlights.**

1. **Spi1 (PU.1) is the strongest single-TF signal in the data.** It
   reaches both the cross-modality leader rule (subsection 14.1) and
   the axis-restricted top 5 of axes 1 and 2 (subsection 14.2 /
   verdict), is a curated AD GWAS locus, and is the canonical master
   microglial lineage TF. The axis-restricted target counts (184
   amyloid-universe / 84 synaptic-universe / not in axis 3 top 15)
   are quantitatively coherent with the literature: PU.1 drives
   inflammatory activation > drives microglial baseline > drives
   biosynthesis.
2. **The Myc / Creb1 / Tp53 / Jun ensemble at axis 3 resolves the
   mixed-sign reading at the gene-set layer.** Section 13's axis-3
   verdict had to call the OXPHOS / ETC / ribosomal signal
   "mixed-sign" because RNA / proteomics push negative while
   phospho-corrected pushes positive at MG-M3; the TF layer,
   integrating across modality stat vectors per contrast, surfaces a
   coherent negative direction. The biological story is therefore
   *not* internal disagreement on the direction of metabolic
   perturbation at tau x amyloid interaction; it is layer-specific
   readout disagreement (RNA + proteomics vs phospho-corrected) of a
   single underlying biology (Myc-led biosynthetic suppression).
3. **Three TFs show sign reversal between amyloid (positive) and
   interaction (negative) contrasts: A0A979HLR9, A0A087WSP5, Rela.**
   The Rela case is the most interpretable: NF-kB p65 is activated by
   amyloid alone and suppressed by tau on the amyloid background,
   exactly the signature a tau-attenuates-amyloid-inflammation model
   would predict. A0A979HLR9 / A0A087WSP5 are CollecTRI heterodimer
   complex sources retained per the 14.2.4 decision; a future UniProt
   annotation pass can resolve their subunit identity.

**Plan-spec deviations explicitly acknowledged.**

1. **Section numbering: 14.3, not 13.3.** The plan stub said "verdict
   subsection (13.3)" but the file `rmd/13_tf_inference.Rmd` renders
   as section 14 in the knit (the parent inserts it after
   `12_pathway_survey.Rmd` which is section 13). The verdict is
   therefore section 14.3 at the rendered TOC, not 13.3. This is a
   stub-text slip that mirrors the E2 stub's "13.1" prediction; the
   plan body has been updated nowhere -- the numbering convention is
   that the file index matches the source filename and the rendered
   section index is one higher than the file index because section 12
   (pathway survey) takes section 13. Future plan stubs should refer
   to verdict-subsection content by anchor ID rather than section
   number to avoid the off-by-one drift.

2. **TSV column set is wider than the plan's "axis, top_TF_per_axis,
   evidence_summary".** The TSV has six columns: `axis`, `top_TFs`
   (top 5 collapsed, comma-separated), `top_TF_signs` (5 +/- chars
   collapsed), `mean_score_range` (e.g. "[-2.57, +4.32]"),
   `n_top_in_cross_modality_leaderboard` (integer count), and
   `evidence_summary` (free-text prose). The plan's three-column
   minimum is satisfied; the extra structured columns make the TSV
   queryable in a spreadsheet without parsing the prose, and they
   make the verdict's quantitative claims auditable against the E4
   axis-restricted TSV without further computation. This is a
   conservative extension, not a deviation that changes the
   inference.

3. **The verdict carries an explicit non-finding for axis 2.** The
   plan's instruction "Summarise the top 3-5 TF drivers per axis" is
   satisfied by the table; the rule on cross-modality consistency
   notes and contrast-direction interpretation is also satisfied. But
   the axis-2 verdict explicitly names that the TF layer does *not*
   identify drivers of synaptic suppression (the per-TF mean is
   identical to axis 1 by Interpretation A construction). The
   verdict reports this as the honest non-finding rather than
   forcing a synthetic ranking, per the plan's
   anti-anchoring guardrail "do not privilege any one axis". Axis 2
   gets equal narrative real estate (its own sub-subsection) but the
   real estate is used to explain the limitation rather than to
   manufacture a finding. This is the *correct* response to the
   data, and is recorded in 14.3.3 and the cross-axis subsection
   14.3.5.

4. **No new ranking statistic introduced.** The verdict aggregates
   the existing E3 + E4 outputs into a single per-axis row; no new
   p-value, no new score, no new threshold. The helper
   `build_tf_verdict_table()` in `R/tf_inference.R` is intentionally
   thin (one tibble per axis, top-N filter, signed-mean range,
   leader-board cross-reference). Editorial prose lives in the Rmd
   chunk via the `evidence_summaries` named-list argument so a
   future revision can edit prose without touching code. This
   matches the section-13 verdict's design (prose in the Rmd,
   aggregation in the helper).

**Anti-anchoring discipline check.** All three D2 axes are named
explicitly in the verdict; the asymmetric reading
(strong / non-finding / strong) is honest rather than synthetic; the
per-axis sub-subsections are equal-length; the inference thresholds
are restated before applying them; no Hallmark / no OXPHOS-arc
terminology appears; hdWGCNA modules are alphabetised wherever they
appear in lists (the verdict mentions MG-M2 and MG-M3 only
individually, so the rule is trivially satisfied here; a multi-module
listing would have followed the rule).

**File outputs.** All chowned `rstudio:rstudio`:

* `R/tf_inference.R` (extended; 68 KB total) -- appended the
  `build_tf_verdict_table()` helper at the end of the file with a
  documentation header block. The helper signature is
  `build_tf_verdict_table(axis_tbl, tf_leaderboard, n_top_per_axis = 5L,
  evidence_summaries = NULL)`. Returns a six-column data.frame, one
  row per axis, preserving the natural amyloid -> synaptic ->
  interaction ordering rather than alphabetising. Existing
  E2/E3/E4 helpers untouched.
* `rmd/13_tf_inference.Rmd` (extended; 43 KB total) -- appended
  subsection 14.3 with five level-3 sub-subsections (14.3.1 ..
  14.3.5). The build chunk reads the E4 axis-restricted TSV from
  disk and the `tf_leaderboard` object from the in-memory knit
  scope (built by the E3 chunk earlier in the same Rmd). The
  evidence-summary strings are defined in the Rmd chunk so future
  editorial edits don't require code changes.
* `storage/results/tf_activity_verdict.tsv` (new; 2,887 B) -- three
  rows (one per axis), six columns. Re-readable with
  `readr::read_tsv()` for downstream Phase F / Phase G consumers.
* `analysis.html` (regenerated; 28 MB) -- subsection 14.3 + the
  five level-3 sub-subsections render cleanly with zero errors /
  warnings. Anchor IDs `verdict-sign-aware-tf-synthesis`,
  `verdict-summary-table`, `amyloid-activation-axis-axis-1-1` and
  `synaptic-suppression-axis-axis-2-1` and
  `interaction-metabolictranslational-axis-axis-3-1` (auto-disambiguated
  from the 14.2 anchors), and `cross-axis-observations-and-limitations`
  all present in the rendered TOC at the expected section numbers.

**Verification.** Pre-knit smoke test (~3 s) confirmed
`build_tf_verdict_table()` returns a three-row data.frame with the
expected columns; ordering is natural amyloid / synaptic /
interaction. Post-knit verification: `grep -c 'class="error"' analysis.html`
returns 0; `grep -c 'class="warning"' analysis.html` returns 0; the
verdict table renders as a real HTML table (caption + colgroup +
tbody cells visible in the html); section-14.3 sub-subsections
14.3.1 .. 14.3.5 all numbered in the rendered TOC.

**Amendment 2026-05-24 (hybrid retrofit, propagated from the F4
decision gate).** At the F4 decision gate the user directed "Do the
hybrid. Have other analyses do the hybrid as well for consistency",
which propagated to E5. `R/tf_inference.R::build_tf_verdict_table()`
was extended with two new hybrid columns: `per_contrast_score_range`
(formatted `"[min, max]"` of per-contrast cross-modality means
across the top-N x axis_contrasts cells) and `per_contrast_summary`
(string like `"nlgf_in_maptki:[-2.09, +5.40]; nlgf_in_p301s:[-3.06, +3.23]"`).
The rendered verdict kable in `rmd/13_tf_inference.Rmd` 14.3 now
shows `per_contrast_summary` alongside `mean_score_range`; the prose
at 14.3.1, 14.3.2, 14.3.3 acknowledges the hybrid view directly --
notably at axis 1, the per-contrast view reveals that maptki carries
the stronger positive activator tail (+5.40 > axis-mean upper +4.32)
while p301s carries the stronger negative repressor tail (-3.06 <
axis-mean lower -2.57), a real biological asymmetry the axis-mean
smooths over. The `tf_verdict_evidence` evidence_summary strings
likewise updated so the exported `tf_activity_verdict.tsv` reads
consistently. Axis 3 (single-contrast) collapses to the axis-mean
range by construction; hybrid columns preserved for schema
consistency. See the F4 completion note (below) for the full E4 + E5
+ F4 retrofit history.

**Phase E COMPLETE.** All five Phase E sessions (E1 decision gate,
E2 cache build, E3 cross-modality ranking, E4 axis-restricted
ranking, E5 verdict) are now DONE. The next active session is the
F1 decision gate for kinase inference; the F1 stub below is no
longer blocked.

---

## Phase F: kinase activity inference (sessions F1-F?)

### Session F1: kinase-inference-tool decision gate [DONE 2026-05-24]

**Completion note.** Decision-gate session. Confirmed the kinase
activity inference tool with the user via `AskUserQuestion`: chose
the default **`decoupleR` + OmniPath kinase-substrate prior** over
alternatives **KSEA only** and **decoupleR + KSEA cross-tool
triangulation**. Locked rationale: methodological consistency with
Phase E (CollecTRI + `decoupleR` for TFs) means the Phase F verdict
in section 15.3 can be read directly against the Phase E verdict in
section 14.3 using identical conventions (`ulm` arbitrates
significance at FDR<0.10, `consensus` arbitrates direction,
`minsize = 5L` universe-overlap filter, axis-restricted ranking via
universe overlap not network re-running). The strongest alternative
(KSEA only) was declined because the z-score formulation, while
statistically independent, would weaken the cross-layer
interpretation that is the analytic gain of Phase F. The
triangulation option (decoupleR + KSEA in parallel) was declined for
cost reasons (~2x F2 build effort for a tool-cross-check that is
not the project's critical-path question at this stage); a future
plan may revisit if the Phase F single-tool verdict surfaces
uncomfortable claims that warrant tool triangulation.

Pre-gate data-shape audit (run this session, results recorded for
F2): phospho cache `storage/cache/de_phospho.rds` has 14,715 sites
with 14,663 (99.6%) annotated with mouse `symbol + PTM.SiteAA +
PTM.SiteLocation` across 3,226 unique gene symbols; batch-corrected
cache `storage/cache/de_phospho_corrected.rds` has 12,821 sites
with 12,773 (99.6%) annotated across 2,682 unique gene symbols.
Both caches' `fit` slot is itself a list with `$fit` (the
`MArrayLM` object) plus `$top` (named list of one limma top-table
per contrast across the canonical five-contrast set
`nlgf_in_maptki`, `nlgf_in_p301s`, `tau_alone`, `tau_in_nlgf`,
`interaction`). Each per-contrast top-table has columns `feature`,
`logFC`, `AveExpr`, `t`, `P.Value`, `adj.P.Val`, `B`, `symbol`,
`PG.ProteinGroups`, `PTM.SiteAA`, `PTM.SiteLocation`. The
`paste0(symbol, "_", PTM.SiteAA, PTM.SiteLocation)` identifier
construction (matching OmniPath KSN substrate IDs) will produce
~14.7k phospho-site IDs from raw and ~12.8k from corrected, easily
clearing any sensible `minsize = 5L` universe-overlap threshold
even after nichenetr-driven mouse <-> human substrate mapping
attrition.

Housekeeping: corrected a typo in the original F1 stub where the
phospho-corrected DE cache was referenced as
`storage/cache/de_phospho_corr.rds`; the actual file on disk is
`storage/cache/de_phospho_corrected.rds`. The corrected name is
used in the F2 Inputs section below. Same slip the E1 session
caught in its own stub for the snRNAseq cache -- future plan stubs
should reference cache files by `ls storage/cache/` output rather
than predicted names.

Locked decisions table now carries two new rows: "Kinase activity
inference tool (Phase F)" and "Mouse orthologue strategy (priors)";
the matching two rows have been removed from the Open questions
table. The remaining Open question is the CCC re-analysis approach
for Phase G1.

Phase F scoped into four atomic sub-sessions F2..F5 below, each
fitting in a single context window, mirroring the E2..E5 shape
(cache build -> cross-cache/cross-contrast leaderboard ->
axis-restricted analysis -> verdict subsection in a new
`rmd/14_kinase_inference.Rmd`, which will render as section 15 of
`analysis.html` per the +1 file-index-to-section-number drift
documented in the E5 completion note). Phospho is single-modality
(no substate or geomx parallels), so F3's "cross-modality" axis from
E3 collapses to a primary "cross-contrast" ranking on the corrected
cache with a secondary "corroborated_by_raw_cache" robustness column
drawing on the raw cache. This is honestly weaker than E3's
three-modality consistency rule and is the genuine data ceiling for
phospho-only inference -- the asymmetry should be named in F3's
interpretive notes rather than smoothed over. No new R or Rmd code
added this session -- planning-only. Verification: this session
changed only `storage/notes/mechanism_layer_plan.md` (a non-Rmd
file outside the knit graph), so re-knit is a no-op; ran it anyway
per execution-model discipline and confirmed analysis.html
unchanged (zero `class="error"`, zero `class="warning"`).

---

### Session F2: build decoupleR + OmniPath kinase activity cache [DONE 2026-05-24]

**Completion note.** Built the Phase F mechanism cache cleanly,
outside the knit, per the plan. The cache file
`storage/cache/kinase_activity_decoupler.rds` is 52 KB
xz-compressed; ten nested tibbles (2 caches x 5 contrasts) all
populated, zero NA scores, statistics levels = 5 per kinase per
cell. Per-cache kinase counts after `minsize = 5L`: raw 104-112,
corrected 98-103. Total wall time 50.2s (KSN build + mouse map
6.3s; phospho_raw 22.9s; phospho_corrected 20.3s).

**Cache shape produced.**

```
phospho_raw         ->  nlgf_in_maptki (555 rows, 111 kinases)
                        nlgf_in_p301s  (530 rows, 106 kinases)
                        tau_alone      (520 rows, 104 kinases)
                        tau_in_nlgf    (560 rows, 112 kinases)
                        interaction    (520 rows, 104 kinases)
phospho_corrected   ->  nlgf_in_maptki (505 rows, 101 kinases)
                        nlgf_in_p301s  (500 rows, 100 kinases)
                        tau_alone      (495 rows,  99 kinases)
                        tau_in_nlgf    (510 rows, 102 kinases)
                        interaction    (490 rows,  98 kinases)
```

Per-cache/contrast row count = kinases x 5 statistics (`ulm`,
`wsum`, `norm_wsum`, `corr_wsum`, `consensus`). Schema matches the
TF cache from E2 byte-for-byte (`statistic`, `source`, `score`,
`p_value` columns).

**KSN attrition.** OmniPath enzyme-substrate table loaded from
OmnipathR cache: 41,506 raw relationships. After
decoupleR-equivalent preprocessing (filter on phosphorylation /
dephosphorylation; build `SYMBOL_resPos` target ids; min(mor)
collapse on duplicate (source, target) pairs): 39,350 human edges,
1,657 kinases, 16,221 substrate ids; mor distribution +1 / -1 =
38,278 / 1,072. nichenetr mouse mapping: kinase symbols 1,602/1,657
(96.7%) mapped; substrate gene symbols 3,684/3,819 (96.5%) mapped;
zero one-to-many collisions on either axis. Mouse KSN after
filtering orphans and re-distincting: 38,143 edges (96.9% retained),
1,583 kinases, 15,845 substrate ids. Overlap with phospho universe:
867 ids (5.9% of phospho_raw); 791 ids (6.2% of phospho_corrected).
112 (raw) / 103 (corrected) kinases retain >=5 substrate coverage
in the respective phospho universe -- a comfortable margin over the
`minsize = 5L` floor.

**Biology highlights** (anticipated F3+ substrate; not analysed
here):

* **Gsk3b is the dominant single-kinase signal across the corrected
  cache**, surfacing at consensus p < 0.05 in 4 of 5 contrasts.
  Sign pattern: nlgf_in_maptki -2.37 (amyloid alone REDUCES GSK3β
  substrate phosphorylation); nlgf_in_p301s +4.50 (amyloid in
  P301S background INCREASES it); tau_alone -3.04 (tau alone
  REDUCES it); tau_in_nlgf +4.81 (tau in NLGF background DRAMATICALLY
  INCREASES it); interaction +4.77 (p=1.83e-06; the synergy IS the
  +interaction). This is the cleanest single result of any Phase E
  or F mechanism output to date and is the most direct possible
  mechanism-layer answer to the project's interaction question
  ("what does tau add on the amyloid background?"). Cdk5 -- the
  second canonical tau kinase -- corroborates at nlgf_in_p301s
  (+3.00) and at interaction (+2.10), confirming a coherent tau-
  kinase-cluster signature.
* **The cell-cycle / stress-kinase cluster** (Cdk1, Cdk2, Mapk14
  i.e. p38, Mapkapk2 i.e. MK2, Csnk2a1 i.e. CK2 alpha) shows
  consistent involvement across multiple contrasts, mirroring the
  cell-cycle-re-entry-in-AD literature.
* **The corrected cache carries more signal than the raw cache**
  (4-5 ulm-significant kinases per contrast vs 0-2 raw),
  vindicating the F1 decision that the corrected cache should be
  primary at F3. The raw cache will serve as a robustness
  corroborator via the F3 `corroborated_by_raw_cache` Boolean
  rather than as a co-primary axis.

**Plan-spec deviations explicitly acknowledged.**

1. **`decoupleR::get_ksn_omnipath()` is broken in the current
   decoupleR 2.16.0 + OmnipathR 3.18.4 stack.** The wrapper
   internally calls the deprecated
   `OmnipathR::import_omnipath_enzsub()` which now errors out on
   its `!!!.` argument splat (`qs_synonyms` fails on the
   evaluated splat). Workaround: I implemented an equivalent
   helper `R/kinase_inference.R::fetch_omnipath_ksn_human()` that
   uses the modern `OmnipathR::enzyme_substrate(organism = 9606)`
   API and reproduces decoupleR's preprocessing identically
   (filter on phosphorylation / dephosphorylation modifications;
   build `SYMBOL_resPos` target id; min(mor) collapse on
   duplicates). The resulting tibble is byte-equivalent in schema
   to what the broken wrapper would have returned. When/if the
   upstream wrapper is repaired, the helper can be replaced with
   a thin shim; the module header records this. The F4 stub's
   instruction "do not re-implement the mapping -- factor F2's
   mapping out into a helper if it isn't already" is satisfied:
   `build_omnipath_ksn_mouse()` is exposed as a public helper
   composing `fetch_omnipath_ksn_human()` + `map_ksn_human_to_mouse()`.
2. **`nichenetr::convert_human_to_mouse_symbols()` requires the
   package to be attached via `library(nichenetr)` rather than
   accessed via `nichenetr::`** because the function reaches for
   the lazy-loaded data object `geneinfo_human` via the package
   search path, not via `::`. Smoke testing surfaced this
   immediately (function errors with `object 'geneinfo_human'
   not found`). Fix: added `library(nichenetr)` to
   `R/helpers.R`'s `suppressPackageStartupMessages()` block (the
   project's existing pattern for package-attachment) with an
   inline comment explaining why bare `requireNamespace()` is
   insufficient. This adds ~1.2 MB of resident memory at knit
   start; acceptable given the project already loads Seurat and
   ComplexHeatmap.
3. **`decouple()` warned about NA cells on phospho_corrected** (the
   E2 stub anticipated this for proteomics; phospho exhibits the
   same pattern). The corrected phospho cache has 1.28% overall
   NA rate, distributed asymmetrically across contrasts (interaction
   2.12% NA, nlgf_in_maptki 0.72% NA). The per-contrast iteration
   strategy in `run_decoupler_per_cache()` (mirror of the E2
   `run_decoupler_per_modality()` pattern) handles this cleanly --
   each contrast's full non-NA subset reaches `decouple()` and no
   gene is dropped just because it has an NA in some other
   contrast. The raw cache shows 0% NA rate (it has 14,715 rows
   per contrast with no estimability holes), confirming the
   asymmetry is a small-group covariance artefact of the batch
   correction step, not a data-quality concern.
4. **Score range has a long positive tail** in `wsum` family
   statistics for a few kinases just above the `minsize = 5L`
   floor (max raw wsum score observed = 184 vs 5th-95th percentile
   range ~5). This is a known property of decoupleR's wsum
   permutation null when the substrate count is small and the
   substrates all move in concert. Downstream F3 ranking will gate
   significance on `ulm`'s linear-model F-test (the locked E3
   convention) which is well-behaved at small N, and rank for
   direction / magnitude on `consensus` (Stouffer-combined across
   ulm + wsum). The wsum raw scores are kept in the cache for
   sensitivity checks but should not be ranked on directly.
5. **The F1 stub said "actual OmniPath KSN edge count
   pre/post-mouse-mapping" should be recorded.** Pre-mouse-mapping
   human KSN: 39,350 edges, 1,657 kinases, 16,221 substrate ids.
   Post-mouse-mapping mouse KSN: 38,143 edges, 1,583 kinases,
   15,845 substrate ids. Attrition is small and the orphan kinases
   (e.g. ACP3, ATP5F1A, AZU1, CDK11A) are mostly either
   non-mouse-conserved or human-only paralogues that would not be
   relevant to a mouse phospho-proteomic dataset anyway.

**File outputs.** All chowned `rstudio:rstudio` (uid 1000, gid 1000):

* `R/kinase_inference.R` (new; 17,413 bytes) -- six helpers
  documented in module-level comment block + inline argument docs:
  `fetch_omnipath_ksn_human()`, `map_ksn_human_to_mouse()`,
  `build_omnipath_ksn_mouse()` (the public composition used by both
  F2 and F4), `extract_phospho_stat_vec()`,
  `extract_phospho_stat_matrix()`, `run_decoupler_per_cache()`,
  `split_kinase_decoupler_by_contrast()`. Substrate-id-keyed
  throughout (key = `paste0(symbol, "_", PTM.SiteAA, PTM.SiteLocation)`);
  duplicate-id collapse by max |t| with rationale in the module
  header (parallels E2's gene-symbol convention).
* `R/helpers.R` (edited) -- inserted `library(nichenetr)` in the
  package block with an inline comment; inserted
  `source("R/kinase_inference.R")` after `source("R/tf_inference.R")`
  in the source-order block with a one-paragraph rationale per the
  E2 convention.
* `scripts/build_kinase_activity_decoupler.R` (new; 5,793 bytes) --
  heavy compute outside the knit, idempotent via `--overwrite`,
  detailed header comment covering inputs / outputs / runtime /
  cache shape / why the decoupleR wrapper bypass was necessary.
* `storage/cache/kinase_activity_decoupler.rds` (new; 53,392 bytes
  ~ 52 KB xz-compressed) -- the cache itself.

**Verification.** Pre-build smoke test of the helpers on a single
`(phospho_corrected, all contrasts)` call passed (~5 s for matrix
extraction + 5 s for KSN build + 20 s for full decoupleR run): row
counts as expected, all 5 contrast keys present in the split list,
no all-NA tibbles. Post-build smoke test on the persisted cache
confirmed: 10 nested tibbles all present; cache keys
`phospho_raw` / `phospho_corrected`; contrast keys = canonical
five-contrast set; statistic levels = 5 per kinase; zero NA scores;
example consensus rows look reasonable (mix of signed scores,
p-values populated and within [0, 1]). Re-knit per execution-model
discipline (since F2 introduces no Rmd consumer chunk yet, the
knit is a no-op for analysis.html); confirmed `analysis.html`
unchanged with zero `class="error"` and zero `class="warning"`.

---

### Session F3: cross-cache/cross-contrast kinase leader board + first 14_kinase_inference subsection [DONE 2026-05-24]

**Completion note.** Wrote the first subsection (15.1) of the new
section 15 in `analysis.html` (file `rmd/14_kinase_inference.Rmd`
renders as section 15 per the +1 file-index-to-section-number drift
documented in the E5 completion note). Knit completed in 3 min 4 s
with zero `class="error"` and zero `class="warning"` in
`analysis.html`. The new section comprises:

* Section 15 "Kinase activity inference" header + intro prose
  cross-referencing the section 14 TF layer and the D2 three-axis
  verdict.
* Subsection 15.1 "Per-contrast kinase activity (corrected cache,
  FDR<0.10)" with one summary table at the top + three level-3
  sub-subsections:
  * 15.1.1 Per-contrast significant kinases (5 tables, including
    explicit "no significant kinases at FDR<0.10" prose for the
    two honest non-finding contrasts `nlgf_in_maptki` and
    `tau_alone`).
  * 15.1.2 Activity score x contrast heatmap for the union of
    kinases sig in any contrast (5 kinases x 5 contrasts; rows
    alphabetical so no kinase is visually privileged; `*` =
    FDR<0.05, `.` = FDR<0.10; NA cells grey).
  * 15.1.3 Interpretive notes (Gsk3b dominance + sign-flip pattern;
    Cdk5 as the canonical-tau-kinase corroborator; Csnk1a1 opposite-
    direction pattern; Csnk1e / Ppp1ca additional surface; honest
    non-findings; cross-reference to the section 14 TF layer and
    the D2 axes; explicit re-statement of the bulk-phospho caveat).

**Per-contrast significant-kinase counts (corrected cache, FDR<0.10
on ulm).** Five contrasts, sig counts 0 / 4 / 0 / 2 / 2:

```
nlgf_in_maptki    n_kinases=101   n_sig=0
nlgf_in_p301s     n_kinases=100   n_sig=4   Cdk5, Gsk3b, Csnk1e, Ppp1ca
tau_alone         n_kinases= 99   n_sig=0
tau_in_nlgf       n_kinases=102   n_sig=2   Gsk3b, Csnk1a1
interaction       n_kinases= 98   n_sig=2   Gsk3b, Csnk1a1
```

Union: 5 unique kinases sig in any contrast (Cdk5, Csnk1a1, Csnk1e,
Gsk3b, Ppp1ca). Two contrasts (`nlgf_in_maptki`, `tau_alone`) are
honest non-findings -- documented as such in the 15.1.1 tables
prose rather than dropped from the section. Gsk3b is the dominant
single-kinase signal, sig in 3 of 5 contrasts (`nlgf_in_p301s` ulm
padj 0.008; `tau_in_nlgf` ulm padj 0.003; `interaction` ulm padj
0.002, the strongest). Its sign pattern is internally coherent
across all 5 contrasts: -1.66 / +3.79 / -1.47 / +4.15 / +4.32 --
amyloid alone on tau-null background trends to *suppress* Gsk3b
substrate phospho; tau alone trends similarly negative; adding tau
or amyloid onto a tau-containing background activates Gsk3b strongly;
the synergy between tau and amyloid *is* the +interaction activation.
This is the most direct possible mechanism-layer answer to the
project's interaction question and the cleanest single result of any
Phase E or F mechanism output to date.

**User-directed F3 deviation from the locked rule.** The original F3
stub locked a cross-cache + cross-contrast leader rule with
`phospho_corrected` as primary and `phospho_raw` as a corroborator,
gated by either `n_contrasts_consistent_sign_ge2_across_caches >= 2`
or `n_contrasts_sig_corrected >= 3`. The pre-build F3 diagnostic
(documented in the next paragraph) showed the rule to be structurally
unsuitable: the cross-cache consistency arm yields zero kinase-
contrast cells at FDR<0.10 (the two caches even disagree on sign at
the dominant Gsk3b signal in three contrasts), and the corrected-
cache arm at `>=3 contrasts sig` yields exactly one kinase (Gsk3b)
because the small inference universe makes BH-padj harsh. The
F3 deviation register's pre-locked fallback (relax to FDR<0.15)
adds 0-1 kinases and does not address the cross-cache arm. The
user-directed decision (2026-05-24, recorded in chat and now
locked in the Locked decisions table) is to drop the leader-board
pattern at the kinase layer entirely: contrasts are treated as
independent inference units; only the corrected cache is consulted;
the per-contrast significance gate is FDR<0.10 on the `ulm`
p-value (E3 convention preserved); a kinase is "significant at
this contrast" or not, with no cross-contrast aggregation. The
resulting per-contrast significant sets become the F3 output (no
unified leader board concept).

**F3 pre-build diagnostic (recorded for F4 / F5 reference).** At
FDR<0.10 on ulm, the corrected cache yields 0 / 4 / 0 / 2 / 2 sig
kinase cells across the canonical five contrasts (total 8 cells over
~100 kinases); the raw cache yields 0 / 0 / 2 / 2 / 0 (total 4
cells). Cross-cache consistency (both caches sig + matching consensus
sign) at FDR<0.10 = 0 cells in any contrast. At a relaxed
FDR<0.15: corrected adds 0-1 cells per contrast; cross-cache
consistency still 0. The corrected-cache arm of the locked leader
rule (`n_sig >= 3`) fires for Gsk3b alone (3 contrasts). At nominal
ulm p<0.05 (no BH adjustment): corrected 5 / 10 / 5 / 10 / 12;
raw 7 / 4 / 13 / 17 / 8. The raw cache shows roughly comparable
nominal-p activity to the corrected cache but lacks the depth of
signal (nominal p < 0.001 in corrected = 4 cells; in raw = 0
cells). This pattern is the signature of batch correction recovering
real signal that batch noise was suppressing.

**Plan-spec deviations explicitly acknowledged.**

1. **Leader-board pattern abandoned (user-directed).** The locked
   E3-parallel pattern (cross-cache + cross-contrast leader rule
   with `corroborated_by_raw_cache` Boolean column) is replaced by
   the per-contrast-independent corrected-cache-only rule above.
   This breaks parallel with section 14.1 in shape but preserves
   the underlying inference threshold (FDR<0.10) and ulm-for-sig +
   consensus-for-direction convention.

2. **Subsection 15.1.4 dropped; 15.1 has three level-3 sub-
   subsections instead of four.** The original plan listed 15.1.1
   per-contrast top-10 tables; 15.1.2 unified leader board;
   15.1.3 heatmap; 15.1.4 interpretive notes. The actual section
   has 15.1.1 per-contrast significant-kinase tables (5 tables, two
   of which carry the honest non-finding prose), 15.1.2 heatmap,
   15.1.3 interpretive notes. The 15.1.2 unified-leader-board sub-
   subsection is gone because there is no leader board to display.
   The summary count table that would have lived in the 15.1.2
   leader-board slot now lives as a single top-of-section
   `kinase-f3-summary-table` chunk at the top of 15.1 (one row per
   contrast: n_kinases_total, n_kinases_sig, top_kinase,
   top_kinase_padj, top_kinase_score). This rearrangement is the
   minimum-impact way to honour the user's directive while keeping
   a readable section structure.

3. **TSV renamed and broadened.** The original plan's
   `kinase_activity_unified_leaderboard.tsv` (one row per leader)
   becomes `kinase_activity_per_contrast.tsv` (one row per
   (contrast, kinase) for all ~100 kinases in the corrected cache,
   500 rows total). Columns: `contrast`, `source`, `score_ulm`,
   `p_value_ulm`, `padj_ulm`, `score_consensus`,
   `p_value_consensus`, `padj_consensus`, `sign_dir`,
   `sig_fdr10_ulm` (Boolean). The wider TSV gives downstream
   sessions (15.2, 15.3) and external consumers everything they need
   to recompute alternative per-contrast filters without re-loading
   the F2 cache.

4. **F4 and F5 stubs need consequential updates.** Both F4 (axis-
   restricted ranking) and F5 (verdict subsection) reference the
   leader-board pattern. F4 still works in principle because it
   ranks by mean activity score across axis-relevant contrasts, not
   by leader-board membership; the "Interpretation A" / "B" choice
   from E4 carries forward unchanged. But the F5 verdict's column
   `n_top_in_cross_modality_leaderboard` (paralleled by the would-
   have-been-locked `n_top_in_cross_cache_leaderboard`) loses its
   referent. The F5 verdict will instead cross-reference the per-
   contrast significant sets directly. Both F4 and F5 stubs below
   have been minimally edited to point at
   `storage/results/kinase_activity_per_contrast.tsv` instead of
   the leader-board TSV; their action sections remain otherwise
   unchanged.

5. **Question to revisit at F4: does the contrast-independence
   directive propagate?** The user-directed F3 rule treats contrasts
   as independent at the leader-board layer. The E4 axis-restricted
   layer aggregates contrasts within each D2 axis by mean activity
   score; this is a different kind of aggregation (within-axis,
   not cross-axis) and may or may not be in scope of the user
   directive. The F4 session must surface this question to the user
   before applying the E4 pattern verbatim, as a decision gate.
   The plausible alternatives at F4 are: (a) keep E4 within-axis
   averaging; (b) report per-contrast scores within each axis
   instead of axis-mean; (c) something else. This is a design fork
   and is logged here so the F4 session does not silently choose.

**Biology highlights.**

1. **Gsk3b is the single dominant kinase signal in the data**, sig
   in 3 of 5 contrasts at the strictest end of FDR (interaction
   padj 0.002, the strongest single result in the section), with a
   sign-flip pattern across the five contrasts that exactly matches
   the project's interaction-contrast logic: tau or amyloid alone on
   the MAPTKI background trends to suppress Gsk3b substrate phospho;
   the combination of tau and amyloid activates it strongly. The
   activation peaks at the `interaction` contrast itself -- the
   synergy *is* the Gsk3b activation. This reads as a direct
   mechanism-layer answer to "what does tau add on the amyloid
   background?"

2. **Cdk5 corroborates Gsk3b at `nlgf_in_p301s`** (ulm +3.91, padj
   0.008) and trends in the same direction at the other tau-
   containing contrasts even though it does not reach significance
   alone. The two canonical tau kinases moving together is a tight
   tau-kinase-cluster signature against the bulk-phospho background.

3. **Csnk1a1 shows the opposite sign pattern**: positive at
   `nlgf_in_maptki`, negative at `tau_in_nlgf` and `interaction`.
   CK1α substrate phospho is *suppressed* when tau is added to the
   amyloid background, opposite to Gsk3b's activation. This argues
   the kinase layer is reading at least two distinct biological
   programs simultaneously: a Gsk3b/Cdk5 axis active in tau-on-
   amyloid contrasts and a CK1α axis active in amyloid-alone-on-
   MAPTKI contrasts. Phase F4 axis-restricted ranking will
   adjudicate whether these are independent axes.

4. **Ppp1ca (PP1 catalytic alpha) at `nlgf_in_p301s`** (ulm +3.09,
   padj 0.053) is the one phosphatase in the significant set. PP1 is
   the major tau phosphatase counteracting both Cdk5 and Gsk3b at
   multiple tau residues; co-activation under amyloid pressure on
   the tau background reads as a coordinated kinase–phosphatase
   axis rather than as a contradiction.

5. **The kinase layer reads the tau-on-amyloid + interaction
   axis**, complementary to the section 14 TF layer's reading of
   the amyloid-driven activation axis. The synaptic suppression
   axis remains unread at either mechanism layer, consistent with
   the section 14.3 verdict that synaptic suppression is more
   plausibly an engulfment phenomenon than a regulator-driven
   program. Phase G CCC re-analysis is the layer expected to
   surface that axis.

**File outputs.** All chowned `rstudio:rstudio`:

* `R/kinase_inference.R` (extended; 34 KB total) -- appended four
  Phase F3 helpers + a ~40-line module-level comment block at the
  insertion point that documents the user-directed deviation from
  the originally locked rule and the rationale (the cross-cache
  arm dead at FDR<0.10; the small inference universe; the dropped
  leader-board concept; the new per-contrast independent rule).
  Helpers: `build_kinase_per_contrast_table()`,
  `filter_sig_kinases_for_contrast()`,
  `plot_kinase_activity_heatmap()`,
  `summarise_sig_kinases_per_contrast()`. F2 helpers untouched.
* `rmd/14_kinase_inference.Rmd` (new; 16 KB) -- section 15 + 15.1
  content. All `knitr::kable()` calls inside `results = 'asis'`
  chunks per CLAUDE.md. The build chunk reads the F2 cache and
  writes the per-contrast TSV; downstream sessions can
  `readr::read_tsv()` the TSV without recomputing.
* `analysis.Rmd` (edited; +3 lines) -- `child-kinase-inference`
  chunk inserted between `child-tf-inference` and `child-session`.
* `storage/results/kinase_activity_per_contrast.tsv` (new; 71 KB)
  -- 500 rows (100 kinases x 5 contrasts), 10 columns
  (`contrast`, `source`, `score_ulm`, `p_value_ulm`, `padj_ulm`,
  `score_consensus`, `p_value_consensus`, `padj_consensus`,
  `sign_dir`, `sig_fdr10_ulm` Boolean).
* `analysis.html` (regenerated; 28 MB) -- section 15 renders
  cleanly with zero errors / warnings. Anchor IDs
  `kinase-activity-inference`,
  `per-contrast-kinase-activity-corrected-cache-fdr0.10`,
  `per-contrast-significant-kinases`,
  `activity-score-x-contrast-view-for-kinases-significant-in-any-contrast`,
  and `interpretive-notes-for-this-subsection-2` (auto-
  disambiguated from the 14.1.4 / 14.2.4 anchors) all present in
  the rendered TOC at the expected numbers (15, 15.1, 15.1.1 ..
  15.1.3).

**Verification.** Pre-knit smoke test on the live cache (~3 s) --
confirmed `build_kinase_per_contrast_table()` returns 500 rows / 9
columns with no NAs in padj or sign columns;
`summarise_sig_kinases_per_contrast()` returns 0 / 4 / 0 / 2 / 2
sig counts; `filter_sig_kinases_for_contrast()` returns the exact
expected kinase identities per contrast;
`plot_kinase_activity_heatmap()` builds a 5 x 5 ComplexHeatmap with
rows alphabetical (Cdk5, Csnk1a1, Csnk1e, Gsk3b, Ppp1ca) and
columns in canonical narrative order. Post-knit verification:
`grep -c 'class="error"' analysis.html` = 0;
`grep -c 'class="warning"' analysis.html` = 0; section 15 +
subsection 15.1 + sub-subsections 15.1.1 / 15.1.2 / 15.1.3 all at
expected positions in the rendered TOC; TSV file on disk at
expected path with 500 rows + header.

---

### Session F3 (original plan, retained for audit) [SUPERSEDED]

**Inputs.** F2 cache `storage/cache/kinase_activity_decoupler.rds`
(no DE caches needed -- everything is already encoded in the F2
cache).

**Action.** Extend `R/kinase_inference.R` with E3-parallel helpers:
* `.extract_kinase_per_cache()` (internal) -- per-cache, per-contrast
  tidy frame of `(source, score, p_value)` rows with `padj` BH-
  adjusted on the `ulm` p-value within each per-cache/per-contrast
  cell (matching the E3 convention of ulm-for-sig + consensus-for-
  direction).
* `rank_kinases_cross_cache()` -- mirror `rank_tfs_cross_modality()`
  but with `cache_a = phospho_corrected` (primary) and
  `cache_b = phospho_raw` (corroborator); FDR<0.10 on `ulm` for sig;
  "consistent sign across caches" requires both caches significant
  with matching `consensus` sign at the same contrast.
* `build_kinase_leader_board()` -- one row per kinase; columns:
  `source`, `n_contrasts_sig_corrected`, `n_contrasts_sig_raw`,
  `n_contrasts_consistent_sign_ge2_across_caches`, `leader_score`
  (composite: `n_consistent * 5 + n_sig_corrected`), `dominant_sign`,
  `corroborated_by_raw_cache` (Boolean: TRUE if the leader's primary
  contrasts are all corroborated by the raw cache),
  `contrasts_summary` (pipe-delimited per-contrast detail string in
  the format `contrast:n_caches_sig/sign`).
* Leader rule: `n_contrasts_consistent_sign_ge2_across_caches >= 2`
  OR `n_contrasts_sig_corrected >= 3` (the latter handles kinases
  robustly detected in the corrected cache even if the raw cache
  misses them -- e.g. when batch correction removes a confounder
  that was suppressing the signal).
* `format_kinase_ranking_table()` -- kable wrapper for the per-
  contrast top-10 + unified leader board tables; render inside
  `results = 'asis'` chunks per CLAUDE.md.

Create `rmd/14_kinase_inference.Rmd` with section header
`# Kinase activity inference` (renders as section 15). Subsection
15.1 "Cross-cache kinase activity ranking" with four level-3
sub-subsections:
* 15.1.1 Per-contrast cross-cache top-10 kinases (5 tables).
* 15.1.2 Unified kinase leader board (single table; row count
  driven by the data; cap at 15 rows if the rule yields more).
* 15.1.3 Kinase x cache activity heatmap at the contrast with the
  most leader kinases (mirror E3's `nlgf_in_maptki` heatmap; pick
  the contrast in the session by `which.max` on per-contrast
  leader counts; `*` = FDR<0.05, `.` = FDR<0.10; NA cells grey).
* 15.1.4 Interpretive notes (face validity of leader kinases
  against microglial / amyloid / tau literature; explicit naming
  of which axes the kinase layer reads vs which it cannot;
  acknowledge the cross-cache vs cross-modality asymmetry vs E3).

Add `child-kinase-inference` chunk in `analysis.Rmd` between
`child-tf-inference` and `child-session`. Export leader board to
`storage/results/kinase_activity_unified_leaderboard.tsv` (mirror
`tf_activity_unified_leaderboard.tsv` shape: one row per leader
kinase, per-contrast detail packed into `contrasts_summary` as a
pipe-delimited string, plus the `corroborated_by_raw_cache`
Boolean column).

**Outputs.** `R/kinase_inference.R` (extended); `rmd/14_kinase_inference.Rmd`
(new); `analysis.Rmd` (edited; +3 lines for the child chunk);
`storage/results/kinase_activity_unified_leaderboard.tsv` (new);
`analysis.html` (re-knitted; new section 15 + subsection 15.1).

**Verification.** Pre-knit smoke test: confirm leader-rule yields
a non-empty leader board, the table has the expected columns, and
the per-contrast top-10 view is balanced (no single contrast
dominating). Post-knit: `grep -c 'class="error"' analysis.html` =
0; `grep -c 'class="warning"' analysis.html` = 0; section 15
renders at the expected position with subsection 15.1 and its
four sub-subsections; TSV file present at the expected path.
Chown all new files `rstudio:rstudio`.

**Plan-spec deviation register.** Watch for: (1) if the corrected
cache's leader kinases are largely disjoint from the raw cache's
leader kinases, the cross-cache consistency rule may be too
restrictive -- record the disagreement in F3's completion note
and defer the decision on which cache is primary to F4 / F5
rather than forcing a single resolution at F3. (2) If the leader
board is too small (< 5 kinases) under the locked rule, consider
relaxing to FDR<0.15 *only* with the rationale documented in the
completion note; do not silently change thresholds.

---

### Session F4: leader-pathway-restricted kinase analysis [DONE 2026-05-24]

**Completion note (decision gate + hybrid retrofit + F4 build, single
session).** Three substantive accomplishments in this session, all
driven by the F4 decision gate at session start:

(1) **F4 decision gate.** Asked the user via `AskUserQuestion`
whether the F4 axis-restricted kinase output should be axis-mean-only,
per-contrast-only, or **hybrid** (both). User chose: "Do the hybrid.
Have other analyses do the hybrid as well for consistency". The
hybrid choice expands scope to retrofit E4 (TF axis-restricted),
E5 (TF verdict), F4 (kinase axis-restricted), and -- by implication
-- F5 (kinase verdict, to be built next session) so every axis-
restricted output uses the same `axis-mean + per-axis-contrast`
column shape. This decision is added to the Locked decisions table
as "Axis-restricted output shape (E4, F4, and propagated to E5)".

(2) **E4 + E5 retrofit (TF layer, user-directed propagation).**
* `R/tf_inference.R::score_tf_per_axis()` extended to emit
  `score_at_<contrast>` columns (one per unique axis_contrast across
  axes, NA where the row's axis does not include the contrast). The
  axis-mean `mean_activity_in_axis_contrasts` is preserved byte-for-
  byte (verified at smoke test: top-5 axis 1 TFs are still
  A0A979HLR9, Spi1, Rreb1, Sp3, Nfkb1 with the same mean magnitudes
  to floating-point precision; row-mean of per-contrast columns
  equals axis-mean to 0.00e+00 in the kinase smoke test and to
  3.08e-02 in the TF smoke test, the latter being the honest
  weighting artefact when cell coverage differs across (modality,
  axis-contrast) cells; documented in the helper docstring).
* `R/tf_inference.R::format_axis_restricted_table()` auto-includes
  per-contrast columns relevant to the axis subset (filters all-NA
  columns per axis) with the `score_at_` prefix stripped for compact
  display; caption appends the per-contrast contrast list.
* `R/tf_inference.R::build_tf_verdict_table()` extended with
  `per_contrast_score_range` (axis-level min/max across top-N x
  axis_contrasts cells) and `per_contrast_summary` (string like
  `"nlgf_in_maptki:[-2.09, +5.40]; nlgf_in_p301s:[-3.06, +3.23]"`).
* `rmd/13_tf_inference.Rmd` subsection 14.3 verdict kable extended
  to render `per_contrast_summary` alongside `mean_score_range`;
  prose at 14.3.1 / 14.3.2 / 14.3.3 updated to acknowledge the
  hybrid view -- in particular, axis 1's per-contrast asymmetry
  (maptki carries the stronger positive activator tail; p301s
  carries the stronger negative repressor tail) is called out as
  a real biological finding the axis-mean smooths over. Axis 3 is
  noted as collapsing to the single `interaction` contrast by
  construction; the hybrid columns are preserved for schema
  consistency.
* `tf_verdict_evidence` strings updated to reference the
  per-contrast ranges so the exported `tf_activity_verdict.tsv`
  `evidence_summary` column matches the prose.

(3) **F4 kinase axis-restricted build.** Appended four F4 helpers
to `R/kinase_inference.R`:
* `restrict_ksn_to_universe(network, universe, min_targets = 5L)` --
  filters the OmniPath KSN to kinases with at least `min_targets`
  distinct in-universe substrate sites. Substrate gene symbol
  extracted from `target` (`symbol_resPos` format) by stripping
  the last `_`-delimited token, robust to symbols containing
  underscores.
* `score_kinase_per_axis(kinase_cache, axis_universes, network, ...)`
  -- mirrors `score_tf_per_axis()` schema exactly; defaults to
  `primary_modalities = "phospho_corrected"` per the F3 locked
  decision (corrected cache only); emits hybrid `score_at_<contrast>`
  columns plus axis-mean.
* `format_axis_restricted_kinase_table()` -- mirrors
  `format_axis_restricted_table()` schema; `n_targets` renamed to
  `n_sites` (substrate sites, not TF targets); auto-includes
  axis-relevant per-contrast columns.
* `plot_axis_lollipop_kinase()` -- mirrors `plot_axis_lollipop()`
  with the kinase x mean coordinates.

Appended subsection 15.2 to `rmd/14_kinase_inference.Rmd` with the
locked 4-subsection structure (build chunk + 15.2.1 amyloid_activation
+ 15.2.2 synaptic_suppression + 15.2.3 interaction_metabolic + cross-
axis observations). Each axis sub-subsection renders a
`format_axis_restricted_kinase_table()` kable (`results = 'asis'`)
and a `plot_axis_lollipop_kinase()` lollipop figure.

**Smoke-test results.** Two stand-alone smoke tests written
(`scripts/smoke_test_hybrid_verdict.R` for the TF retrofit;
`scripts/smoke_test_f4_kinase.R` for the F4 kinase build). Both
exercise the helpers end-to-end against the live caches; both pass.
Universe sizes per axis: amyloid_activation 3395 genes / 379 kinases
with >=5 KSN sites; synaptic_suppression 3940 / 425; interaction_metabolic
1859 / 262. Top-5 per axis: axis 1 = Cdk2 (+2.91), Cdk5 (+1.75),
Csnk1e (+1.64), Mapk8 (+1.17), Camk2g (-1.16); axis 2 = identical to
axis 1 by Interpretation A; axis 3 = **Gsk3b (+4.77)**, Csnk1a1
(-2.81), Mapk14 (+2.63), Cdk5 (+2.10), Cdk1 (+1.99). Row-mean of
per-contrast columns equals axis-mean exactly (max diff 0.00e+00) at
the kinase layer because the corrected cache is single-modality with
uniform cell coverage across contrasts; the hybrid columns are still
preserved for schema consistency and for the F5 verdict layer to
consume.

**Plan-spec deviation register.**
1. **Hybrid retrofit propagated to E4 + E5 (user-directed expansion).**
   The original F4 stub scoped this session as F4 alone; the user's
   "do the hybrid for consistency" directive expanded scope to also
   retrofit `score_tf_per_axis()`, `format_axis_restricted_table()`,
   `build_tf_verdict_table()`, and the rmd/13 verdict prose. Scope
   handled in one session; the next session (F5) inherits the now-
   locked hybrid schema and only needs to build the kinase verdict
   analogue.
2. **`build_axis_gene_universe_for_kinase()` not added as a separate
   helper.** The original F4 stub called for extending the E4 gene-
   universe builder to include the `phospho_corrected` gene universe.
   Skipped because the kinase smoke test confirmed that the existing
   E4 gene universes already yield 262-425 kinases per axis at the
   `min_targets = 5L` filter -- well above the "at least 5 kinases
   per axis" sanity bound. Adding the phospho gene universe to the
   per-axis union would have expanded the universes by ~5-15% and
   shifted the rankings by sub-decimal amounts, at the cost of
   diverging from the E4 axis-universe definition (which would have
   made cross-layer comparison harder). Decision: re-use
   `build_axis_gene_universe()` as-is for both TF and kinase layers.
3. **F4 axis 2 is NOT empty.** The original F4 stub flagged that
   axis 2 might be empty at the kinase layer (synaptic engulfment
   being post-transcriptional but not necessarily kinase-mediated).
   In fact axis 2 carries 103 kinases (vs 101 at axis 1 and 96 at
   axis 3) but is identical to axis 1 in its top kinases by
   Interpretation A construction -- the entanglement is preserved.
   The synaptic-suppression-axis-2 honest non-finding therefore takes
   the form "top kinases are inherited from axis 1, no distinct
   synapse-repressor program emerges" rather than "no kinases
   survive the axis 2 filter". 15.2.2 prose handles this explicitly.
4. **Bulk-phospho caveat re-stated in 15.2 prose.** The "phospho is
   bulk hippocampus, not microglia-sorted" caveat from 15.1 is
   inherited and re-stated where relevant in 15.2 prose (so a reader
   who lands on 15.2 directly is not misled).
5. **Cdk5 cross-axis pattern is the integrating signal.** Unanticipated
   at plan time: Cdk5 surfaces in the top 5 at all three axes
   (rank 2 at axes 1+2, rank 4 at axis 3), the only kinase to do so.
   The cross-axis substrate-set differential (81 amyloid / 263 synapse /
   153 interaction) is biologically interpretable as Cdk5 being the
   integrator across the three axes. Noted in 15.2.4 prose as a
   future-Cdk5-substrate-resolved-sub-analysis target.

**Outputs.**
* `R/tf_inference.R` extended (~150 lines added at the hybrid retrofit).
* `R/kinase_inference.R` extended (~378 lines added at the F4 helpers).
* `rmd/13_tf_inference.Rmd` updated (hybrid prose at 14.3.1-14.3.5;
  verdict kable + caption + tf_verdict_evidence).
* `rmd/14_kinase_inference.Rmd` extended (~371 lines added at 15.2).
* `scripts/smoke_test_hybrid_verdict.R` (new) and
  `scripts/smoke_test_f4_kinase.R` (new); kept for future regression
  testing.
* `storage/results/tf_activity_axis_restricted.tsv` (rewritten with
  hybrid columns; same row count as before, +3 columns).
* `storage/results/kinase_activity_axis_restricted.tsv` (new; 300
  rows, schema documented above).
* `analysis.html` re-knit (see Verification below).

**Verification.** Both smoke tests pass with zero errors / warnings;
all four format/plot helpers render correctly per axis; the row-mean
== axis-mean invariant holds to floating-point precision at the
kinase layer and to the documented honest weighting artefact at the
TF layer. Post-knit verification: zero `class="error"` and zero
`class="warning"` in `analysis.html`. All new files chowned to
`rstudio:rstudio`. Committed locally with `Co-Authored-By: Claude
Opus 4.7` trailer per CLAUDE.md.

---

### Session F5: kinase inference verdict subsection [DONE 2026-05-24]

**Inputs.** F3 per-contrast TSV
(`storage/results/kinase_activity_per_contrast.tsv` -- replaces
the originally planned `kinase_activity_unified_leaderboard.tsv`
per the F3 user-directed deviation; 500 rows, one per (contrast,
kinase), columns documented in the F3 completion note deviation
register entry #3); F4 axis-restricted TSV
(`storage/results/kinase_activity_axis_restricted.tsv`); in-memory
`kinase_per_contrast` tibble from the F3 chunk (built earlier in
the same Rmd by `build_kinase_per_contrast_table()`, in scope per
the parent knit's shared-R-session convention; replaces the
originally referenced `kinase_leaderboard` object that the
user-directed deviation removed).

**Action.** Append subsection 15.3 "Verdict (sign-aware kinase
synthesis)" to `rmd/14_kinase_inference.Rmd` mirroring E5's 14.3
shape. Six level-3 anchors:
* 15.3 top-line + threshold restatement (per-contrast FDR<0.10 on
  `ulm` against the corrected cache only, per the F3 user-directed
  deviation; axis-restricted `min_targets = 5`, top-5-per-axis
  displayed; phospho-is-bulk caveat).
* 15.3.1 Verdict summary table (one row per axis; six columns:
  `axis`, `top_kinases`, `top_kinase_signs`, `mean_score_range`,
  `n_top_sig_in_axis_contrasts`, `evidence_summary`). The fifth
  column counts how many of `top_kinases` are flagged
  `sig_fdr10_ulm = TRUE` in *any* of the axis-relevant contrasts
  per the F3 per-contrast TSV; this replaces the originally planned
  `n_top_in_cross_contrast_leaderboard` column that the F3 user-
  directed deviation removed.
* 15.3.2 Amyloid-activation axis (axis 1).
* 15.3.3 Synaptic-suppression axis (axis 2) -- anticipated
  non-finding per the E5 axis-2 framing; report honestly as such
  if the data supports rather than manufacturing a finding.
* 15.3.4 Interaction metabolic/translational axis (axis 3).
* 15.3.5 Cross-axis observations and limitations including
  explicit cross-reference to the TF verdict in section 14.3 --
  the Myc / Creb1 / Tp53 / Jun cluster at TF axis 3 has plausible
  kinase-layer correlates (Cdk1/2/4, Mapk1/3, Gsk3a/b
  respectively); record observed parallels and unexpected
  absences. Also record any sign-reversal patterns observed
  between amyloid and interaction contrasts (mirror of the
  Rela parallel sign reversal observed at E4 / E5).

Add `build_kinase_verdict_table()` helper to `R/kinase_inference.R`
mirroring `build_tf_verdict_table()` (thin aggregator; one tibble
per axis; top-N filter; signed-mean range; cross-reference to the
F3 per-contrast significant sets via the per-contrast TSV's
`sig_fdr10_ulm` column rather than to a leader-board membership
column per the F3 user-directed deviation; editorial prose lives
in the Rmd chunk via the `evidence_summaries` named-list argument
so future revisions can edit prose without touching code).

**Outputs.** `R/kinase_inference.R` (extended);
`rmd/14_kinase_inference.Rmd` (extended; subsection 15.3 added);
`storage/results/kinase_activity_verdict.tsv` (new; three rows,
six columns; re-readable for Phase G consumers); `analysis.html`
(re-knitted).

**Verification.** Pre-knit smoke test: `build_kinase_verdict_table()`
returns a three-row data.frame with the expected six columns;
ordering natural amyloid / synaptic / interaction (not
alphabetised; preserve the D2 axis order). Post-knit: zero errors
/ warnings; verdict table renders as a real HTML table (caption
+ colgroup + tbody cells visible in the html); subsection 15.3 +
sub-subsections 15.3.1 .. 15.3.5 all numbered in the rendered TOC.
Chown new files `rstudio:rstudio`.

**Anti-anchoring discipline check (per E5 convention).** All three
D2 axes named explicitly in the verdict; asymmetric reading (if
applicable -- likely is, per the phospho-only data ceiling) is
honest rather than synthetic; per-axis sub-subsections equal-
length even when one is a non-finding; inference thresholds
restated before applying them; no Hallmark / no OXPHOS-arc
terminology; hdWGCNA modules alphabetised wherever they appear in
multi-module lists.

**Phase F COMPLETE marker.** F5 completion note should explicitly
declare Phase F COMPLETE and update the plan's top-of-file `State:`
line to mark F1..F5 all DONE. Phase G1 (CCC re-analysis approach
decision gate) becomes the next active session at that point.

**Completion note (2026-05-24).** Helper `build_kinase_verdict_table()`
appended to `R/kinase_inference.R` (lines 1193-1307; file grew from
1122 -> 1307 lines) as a thin aggregator mirroring
`build_tf_verdict_table()`: derives `axis_contrasts` per axis from
non-all-NA `score_at_<contrast>` columns rather than requiring an
explicit `axis_universes` argument; counts unique kinases flagged
`sig_fdr10_ulm = TRUE` in any axis-relevant contrast for the
`n_top_sig_in_axis_contrasts` column; pulls editorial prose from the
caller-supplied `evidence_summaries` named list. Subsection 15.3 +
sub-subsections 15.3.1 .. 15.3.5 appended to
`rmd/14_kinase_inference.Rmd` (file grew from 675 -> 1194 lines).
Verdict TSV written to `storage/results/kinase_activity_verdict.tsv`
(header + 3 axis rows = 4 lines on disk). Pre-knit smoke test
`scripts/smoke_test_f5_kinase_verdict.R` passes cleanly: 3-row
data.frame returned, axis order natural amyloid / synaptic /
interaction (not alphabetised), top_kinases per axis match the F4
expected ranking, `n_top_sig_in_axis_contrasts` integers (2 / 2 / 2)
cross-validated against the F3 per-contrast significant sets in
the same script. Post-knit verification: 0 errors / 0 warnings in
`analysis.html`; section 15.3 + 15.3.1-15.3.5 render as five
numbered HTML anchors (`verdict-sign-aware-kinase-synthesis`,
`verdict-summary-table-1`, plus three axis sections and the
cross-axis section) at the expected positions immediately after
the F4 axis-restricted subsection 15.2.

**Schema deviation from F5 stub (additive only; documented here
rather than as a separate stub revision).** The stub specified a
6-column verdict (`axis`, `top_kinases`, `top_kinase_signs`,
`mean_score_range`, `n_top_sig_in_axis_contrasts`,
`evidence_summary`). The helper actually emits 8 columns by carrying
through the F4 hybrid pattern (axis-mean + per-contrast columns
side-by-side; user-directed deviation locked at F4 and retrofitted
to E4 / E5): `axis`, `top_kinases`, `top_kinase_signs`,
`mean_score_range`, `per_contrast_score_range`,
`per_contrast_summary`, `n_top_sig_in_axis_contrasts`,
`evidence_summary`. The two new hybrid columns
(`per_contrast_score_range`, `per_contrast_summary`) expose the
per-axis-contrast min / max envelope and a contrast-by-contrast
breakdown so consumers can detect tau-on-amyloid substructure that
the axis-mean smooths over -- exactly the F4 motivation, propagated
to F5 for shape consistency. The 2026-05-24 E5 hybrid retrofit
recorded earlier in this plan made the equivalent change to the TF
verdict; F5 inherits that retrofit by construction. The deviation
is purely additive (the original 6 columns are all present +
extended) and downstream consumers (Phase G) only widen.

**Findings (one paragraph per axis; full evidence_summary strings
in the verdict TSV).**
- *Axis 1 (amyloid_activation).* Top 5 = Cdk2 / Cdk5 / Csnk1e /
  Mapk8 / Camk2g (4-positive / 1-negative). Signed mean range
  [-1.16, +2.91] across the two NLGF contrasts; hybrid per-contrast
  view nlgf_in_maptki:[-1.32, +3.02], nlgf_in_p301s:[-1.00, +3.00]
  exposes Cdk5 / Mapk8 as p301s-biased (maptki +0.50 / +0.25 vs
  p301s +3.00 / +2.08) -- two of the top 5 would not have ranked
  here on the MAPTKI contrast alone, a real tau-on-amyloid
  substructure the axis-mean smooths over. 2 of the top 5 (Cdk5,
  Csnk1e) reach FDR<0.10 on ulm at nlgf_in_p301s per 15.1; the
  remaining three are non-significant at both NLGF contrasts.
  Reading: Cdk2 (cell-cycle progression in microglial proliferative
  re-entry), Cdk5 (canonical tau kinase, p301s-biased here), Csnk1e
  (CK1-epsilon, stress-activated), Mapk8 (JNK, stress-activated)
  form a proliferation + stress-kinase ensemble; Camk2g is the lone
  repressor (CaMKII-gamma, calcium / calmodulin signalling
  attenuated). Pairs with the section-14.3 TF axis-1 (Spi1 / Nfkb1
  / Sp3 DAM-program TFs) as the phospho-signalling effector layer
  of the same DAM program.
- *Axis 2 (synaptic_suppression).* Top 5 identical to axis 1 by
  Interpretation A construction (both axes share their NLGF
  contrasts and the per-kinase mean is computed from full-cache
  scores; the axis 1 / axis 2 entanglement propagates to the
  per-contrast layer as a structural consequence). The structural
  differentiator lives in the in-universe substrate-site counts:
  Cdk5's footprint expands 81 -> 263 sites (3.2x) when the
  substrate universe moves from amyloid to synaptic; Cdk2 259 ->
  293, Csnk1e 38 -> 70, Camk2g 24 -> 49. The synaptic axis is
  therefore best interpreted as a target-overlap filter ("which
  amyloid-active kinases have substrate programs that also touch
  the synaptic machinery") rather than as a list of kinases that
  drive synaptic suppression. This is the honest non-finding -- the
  kinase layer does not independently identify drivers of synaptic
  suppression at this resolution; the top kinases are positive
  activators of broadly synapse-relevant substrate programs (Cdk5
  in particular is a well-documented synaptic regulator), not
  repressors of synaptic gene output. Consistent with section-14.3
  TF axis-2 and the engulfment-mediated mechanism hypothesis raised
  in 14.1.4 / 14.2.4 / 15.2.4; Phase G (CCC sender-receiver) is the
  layer that can plausibly read engulfment biology directly.
- *Axis 3 (interaction_metabolic).* Top 5 = Gsk3b / Csnk1a1 /
  Mapk14 / Cdk5 / Cdk1 (4-positive / 1-negative). Signed mean range
  [-2.81, +4.77] at the single interaction contrast -- widest of
  the three axes. Gsk3b at ulm padj 0.002 is the strongest
  single-kinase result of the whole mechanism arc; Csnk1a1 reaches
  FDR<0.10 in the opposite direction. Reading: Gsk3b (canonical
  tau kinase phosphorylating tau at S202 / T205 / S396 / S404, the
  AT8 / PHF1 epitopes that stage AD tau pathology), Mapk14 (p38-
  alpha MAPK, stress-kinase activator of tau phosphorylation and
  microglial activation), Cdk5 (second canonical tau kinase,
  p35 / p25 activation cycle), and Cdk1 (master mitotic kinase,
  microglial cell-cycle re-engagement) form a canonical tau-kinase
  activation ensemble at the tau x amyloid interface; Csnk1a1
  (CK1-alpha) is the opposite-direction signal whose substrate
  program is suppressed on the tau-on-amyloid background. Combined
  with section-14.3 TF axis-3 (Myc / Creb1 / Tp53 / Jun coordinated
  suppression of biosynthetic activity), the joint two-layer
  mechanism reading at axis 3 is "the tau x amyloid interface
  drives tau hyperphosphorylation through canonical tau kinases
  (Gsk3b lead) while transcriptionally shutting down biosynthetic
  throughput (Myc lead)" -- exactly the two-layer answer the
  mechanism arc was designed to expose. The section-13 axis-3
  mixed-sign reading (OXPHOS / ribosomal compartment + MG-M3 +
  DAM_up at interaction; modality-divergent signs) resolves under
  this joint mechanism: the proteomics + phospho signal of
  biosynthetic + phospho dysregulation runs in the same direction
  once the kinase + TF layers are read together.

**Phase F COMPLETE.** All five F sessions DONE (F1 decision gate
locking decoupleR + OmniPath as the kinase tool; F2 corrected-
cache build with mouse symbol mapping via nichenetr; F3 per-
contrast significance with leader-board pattern abandoned per
user-directed deviation; F4 axis-restricted hybrid view with E4 /
E5 retrofitted to match; F5 verdict subsection 15.3 with helper
+ TSV + multi-axis joint mechanism reading). Phase G1 (CCC
re-analysis approach decision gate) is the next TODO step,
unblocked at the close of this session.

---

## Phase G: CCC bridge re-analysis with leader-pathway filters (sessions G1-G?)

### Session G1: CCC re-analysis approach decision gate [DONE 2026-05-24]

**Completion note.** Decision-gate session. Surveyed the existing
CCC infrastructure before posing the question so the user had a
real choice grounded in what already exists: section 11 covers
CellChat (per-condition signalling, differential pathway ranking,
pairwise differential interactions) + MultiNicheNet (top L-R pairs
per contrast, circos, sanity tables) across subsections (a)-(f);
three heavy caches pre-baked outside the knit
(`storage/cache/cellchat_per_condition.rds` ~1 GB,
`storage/cache/cellchat_merged.rds` ~67 MB,
`storage/cache/multinichenet_output.rds` ~277 MB) built by
`scripts/build_cellchat.R` and `scripts/build_multinichenet.R` on
the symbol-keyed `seurat_full_processed.rds`; existing TSV outputs
already shaped for post-filtering (ligand, receptor,
prioritization_score, etc. as columns); leader-board TSV
`storage/results/pathway_survey_unified_leaderboard.tsv` (145 lines
including header, 144 leader rows) provides the per-axis gene
universes via the same gene-set caches already loaded at E4 / F4
(`build_axis_gene_universe()` in `R/tf_inference.R`).

Confirmed two locked decisions with the user via
`AskUserQuestion` (two-question form):

1. **CCC re-analysis approach = three-tool triangulation
   (CellChat + MultiNicheNet + liana) with leader-pathway-restricted
   gene-universe post-filter applied across all three.** User-
   directed deviation from the plan's default post-filter-only
   proposal (option 1, my recommendation; declined because it
   would replicate F1's single-tool-suffices logic on an axis
   where prior tools already failed). User-directed deviation from
   the F1 precedent (which declined an analogous decoupleR + KSEA
   triangulation for the kinase tool; the defining asymmetry here
   is that the synaptic-suppression axis is an explicit non-finding
   at the TF and kinase layers, and the verdict prose at 14.3.3 /
   15.2.4 / 15.3.3 specifically named CCC as the layer expected to
   surface engulfment-mediated biology; adding tools to read an
   axis that prior tools couldn't read is the natural response).
   Re-run-from-scratch with restricted prior (option 2) declined
   because it breaks comparability to section 11 and changes a
   locked input. Hybrid post-filter-now-escalate-later (option 4)
   declined because it defers the triangulation decision rather
   than answering it now. Decision recorded in Locked decisions
   table; the "CCC re-analysis approach (Phase G)" row of the Open
   questions table removed.
2. **CCC mechanism Rmd venue = new `rmd/15_ccc_mechanism.Rmd`
   rendering as section 16** (per recommended option). Mirrors the
   per-phase one-Rmd-per-mechanism-layer convention from
   `rmd/13_tf_inference.Rmd` (section 14, TF) and
   `rmd/14_kinase_inference.Rmd` (section 15, kinase). Reader sees
   TF + kinase + CCC mechanism grouped as a coherent mechanism arc
   in the TOC at sections 14 / 15 / 16; section 11 remains the
   descriptive CCC home and is not modified by Phase G. The plan
   stub's preference for extending `rmd/11_ccc.Rmd` was considered
   and declined per the question's option-2 framing. Decision
   recorded in Locked decisions table.

A new row "liana R-vs-Python implementation (Phase G2)" added to
the Open questions table — the G2 build session must investigate
liana availability (archived R package vs actively-maintained
`liana+` Python; the latter requires reticulate) and resolve before
the cache build can begin. This sub-decision is not Phase-level so
does not need its own G-phase decision gate; the G2 session
handles it as a sub-decision-gate at session start if both paths
remain viable after the investigation.

Phase G scoped into four atomic sub-sessions G2..G5 below, each
fitting in a single context window, mirroring the F2..F5 shape
(cache build → per-contrast cross-tool leaderboard at 16.1 →
axis-restricted analysis at 16.2 with hybrid axis-mean +
per-axis-contrast columns → verdict at 16.3 in the new
`rmd/15_ccc_mechanism.Rmd`, which will render as section 16). No
new R or Rmd code added this session — planning-only.

Housekeeping: no cache-file rename slip surfaced at this gate (in
contrast to E1 and F1, both of which caught a `de_phospho_corr.rds`
rename in their own stubs). The audit confirms
`cellchat_per_condition.rds`, `cellchat_merged.rds`,
`multinichenet_output.rds` are the actual filenames on disk and
match the G1 stub's `storage/cache/ccc_*.rds` glob (the glob
matches by prefix on the parent script names, not the cache names
themselves — the cache names are tool-keyed, not section-keyed).
The plan reference text in G2 / G3 / G4 / G5 below uses these
verified filenames verbatim.

**Verification.** This session changed only
`storage/notes/mechanism_layer_plan.md` (a non-Rmd file outside the
knit graph), so re-knit is a no-op for analysis.html; ran it anyway
per execution-model discipline and confirmed analysis.html
unchanged (zero `class="error"`, zero `class="warning"`).

---

### Session G2: build liana CCC cache [DONE 2026-05-24]

**Completion note.** Built the Phase G CCC cache cleanly, outside
the knit, per the plan. Total wall time 4 min 20 s (260 s for the
four per-genotype `rank_aggregate` runs; <1 s for the per-contrast
derivation and provenance build). Cache file
`storage/cache/liana_output.rds` is 1.58 MB xz-compressed; three
top-level slots (`per_genotype` x 4, `per_contrast` x 5,
`provenance`); all five canonical contrasts populated; **zero NAs**
in the prioritization / significance / sign columns across all
contrasts (G2 verification gate PASS). Per-contrast tibbles carry
the eight cross-tool aggregation columns required for the G3 rbind
(`sender`, `receiver`, `ligand`, `receptor`, `lr_interaction`,
`prioritization_score`, `n_methods_pass_significance`, `sign_dir`)
plus four to seven auxiliary diagnostic columns; per-genotype
tibbles carry the full liana_res schema (13 cols) plus a
`cellphone_padj` column (BH-corrected within each per-genotype
run, scope-correct because cross-tool sig is per-genotype not per
(genotype, contrast)).

**Sub-decision-gate resolution.** Chose Python LIANA+ v1.7.1 via
reticulate against a project-local `.venv/bin/python`. R `liana`
was deprecated by its authors in 2024 in favour of the Python
`liana+` rewrite (van Dijk / Dimitrov *et al.* 2024, Nat Cell
Biol); no R-side `liana` package was installed in the container.
User confirmed via `AskUserQuestion` at the sub-gate. Pure-R
fallback (extra single-method tools alongside CellChat /
MultiNicheNet) not invoked because the reticulate path validated
cleanly via a tiny pre-build smoke test and is the SOTA for this
problem class.

**Cache shape produced.**

```
per_genotype  ->  MAPTKI       ( 5055 LR-pair rows, 14 cols)
                  P301S        ( 6683 LR-pair rows, 14 cols)
                  NLGF_MAPTKI  ( 9666 LR-pair rows, 14 cols)
                  NLGF_P301S   (10087 LR-pair rows, 14 cols)
per_contrast  ->  tau_alone        ( 6838 rows, 12 cols; 4095 sig in primary)
                  nlgf_in_maptki   ( 9755 rows, 12 cols; 6141 sig in primary)
                  nlgf_in_p301s    (10446 rows, 12 cols; 6763 sig in primary)
                  tau_in_nlgf      (11220 rows, 12 cols; 6763 sig in primary)
                  interaction      (11416 rows, 12 cols; 6763 sig in primary)
provenance    ->  liana 1.7.1 / anndata 0.12.16 / scipy 1.17.1 / python 3.12.3
                  n_cells_per_genotype:  77684 / 69921 / 69518 / 69162
                                         (MAPTKI / P301S / NLGF_MAPTKI / NLGF_P301S)
                  n_pairs_per_genotype:   5055 /  6683 /  9666 / 10087
                  params: resource_name=mouseconsensus, use_raw=FALSE,
                          expr_prop=0.10, min_cells=5, n_perms=100,
                          seed=1, sig_alpha=0.10
```

Per-contrast row counts exceed any single per-genotype row count
because the outer-join over the contrast's relevant genotype pair
(or four-way pair for interaction) captures LR pairs detected in
at least one of the joined genotypes; the absent-in-other-genotype
slots are imputed as worst-rank (`magnitude_rank = 1.0`) and
worst-padj (`cellphone_padj = 1.0`) so the G2 zero-NA gate holds.

**Biology highlights** (anticipated G3+ substrate; not analysed
here):

* **`Psen1 -> Ncstn` surfaces as the second-strongest leader in
  `tau_alone` (prio 0.995, neuronal-to-neuronal, sig=1).** This is
  the γ-secretase complex (presenilin 1 binding nicastrin), the
  single most-AD-relevant LR pair in the entire output. Its
  appearance in the `tau_alone` contrast (P301S vs MAPTKI, pure-tau
  effect with no amyloid) suggests tau pathology re-tunes
  γ-secretase signalling in a cell-autonomous manner -- a clean
  CCC-layer corroboration of the canonical "tau amplifies APP
  processing" story.
* **`Gas6 -> Axl` (microglial TAM-kinase phagocytosis) recurs as a
  top-3 leader in BOTH `nlgf_in_p301s` AND `tau_in_nlgf`** (neuron
  -> astrocyte, prio 0.996 in both contrasts). Gas6/Axl is the
  canonical signalling axis for microglial phagocytic clearance of
  amyloid plaques and dying neurons in AD; its appearance in the
  two contrasts that involve adding amyloid OR adding tau on the
  amyloid background is biologically coherent and ties the CCC
  layer to the F5 microglia-DAM kinase verdict.
* **`L1cam -> Ephb2` (neuronal axon guidance) tops `tau_alone` (prio
  0.999) and ranks 2 in `nlgf_in_maptki` (prio 0.998).** A
  cross-contrast neuronal axon-guidance reorganisation signal
  consistent with synaptic-pruning literature in tauopathy /
  amyloidopathy mouse models.
* **All top-3 leaders of the `interaction` contrast are NEGATIVE
  sign** (prio in [-1.56, -1.31]): `Rspo2 -> Znrf3` (Wnt signalling,
  neuron -> microglia_proliferative), `Nid1 -> Itgav` and
  `Nid1 -> Itgb1` (basement-membrane integrin signalling, vascular
  -> neuron and vascular -> astrocyte). All three indicate the
  NLGF effect runs **opposite** in the P301S background compared
  with MAPTKI -- a genuine divergent-NLGF biology signal that
  vindicates including the interaction contrast in the CCC layer
  and is exactly the shape of cross-axis surprise the project
  hypothesised would be detectable.
* The five-contrast significance counts (4095 / 6141 / 6763 / 6763
  / 6763) are not flat: `tau_alone` has the fewest sig pairs
  because P301S has the smallest per-genotype LR pair count (6683),
  whereas the four NLGF-touching contrasts inherit
  `NLGF_P301S`'s 6763 (or `NLGF_MAPTKI`'s 6141) primary-genotype
  sig count. This is the expected behaviour of the per-genotype
  sig-flow design.

**Plan-spec deviations explicitly acknowledged.**

1. **Symlink trap in venv-python resolution.** A Python venv's
   `bin/python` is a symlink chain to the system Python (here:
   `.venv/bin/python -> python3 -> /usr/bin/python3`); the venv
   context is determined by `pyvenv.cfg` living next to the SYMLINK,
   not next to the resolved system binary. `normalizePath()` follows
   the symlink and destroys the venv context, causing reticulate to
   bind to the system Python with no liana site-packages. The
   bootstrap helper `R/ccc_inference.R::.ensure_liana_python()`
   deliberately constructs the absolute path without following
   symlinks (`if (substr(venv, 1L, 1L) == "/") venv else
   file.path(getwd(), venv)`). Build script sets
   `RETICULATE_PYTHON` BEFORE `library(reticulate)` for
   belt-and-braces (Python initialises once per R session; env var
   must be in place at library-load time).
2. **`reticulate::py$x <- v` is invalid R syntax.** R cannot assign
   to a namespace-qualified slot expression; the fix is to bind
   `py <- reticulate::py` first then `py$x <- v`. Caught during
   validation, fixed in `run_liana_rank_aggregate_for_genotype()`.
3. **Per-genotype design, not per-cell-with-design-matrix.** liana's
   `rank_aggregate()` does not do contrasts internally; it produces
   per-condition prioritisations. To match the project's 2x2
   factorial and the five canonical contrasts used elsewhere, fit
   liana four times -- once per genotype -- then derive contrasts as
   deltas of per-genotype `magnitude_rank` (flipped to
   `comm_strength = 1 - magnitude_rank` first; lower-is-better in
   liana, so the flip aligns with the project-wide
   higher-is-stronger convention). Same pattern CellChat uses in
   `build_cellchat.R` (per-condition fits then `mergeCellChat`).
   This is the resolution of deviation-register item (4) from the
   original G2 stub.
4. **`n_methods_pass_significance` is binary per-tool (0 or 1),
   NOT 9 (the count of liana's internal methods).** Liana
   internally aggregates 9 methods via RobustRankAggregate within
   each per-genotype run. At the cross-tool level (G3) we treat
   liana as ONE tool with ONE verdict per pair per contrast --
   significant iff `cellphone_padj < 0.10` (BH within the
   per-genotype run) in the contrast's primary genotype. This is
   the honest framing for cross-tool sig counting; it avoids
   double-counting liana's internal multi-method evidence as if it
   were independent of CellChat and MultiNicheNet. alpha = 0.10
   matches F3 / G3 cross-tool conventions.
5. **Auxiliary diagnostic columns beyond the cross-tool required
   eight.** Each simple-contrast tibble adds `mag_rank_primary`,
   `mag_rank_reference`, `cellphone_padj_primary`,
   `cellphone_padj_reference`; the interaction tibble adds
   `delta_in_primary_pair`, `delta_in_reference_pair`,
   `cellphone_padj_sig_primary`. These do not break the G3 rbind
   (G3 will select the required 8 cols) and they support G3
   visualisations / diagnostic prose.
6. **Missing-pair imputation rule.** LR pairs failing `expr_prop`
   (>=10% cells in sender OR receiver) or `min_cells` (>=5) in one
   genotype are absent from that genotype's `liana_res`. After the
   outer-join during contrast derivation, missing values are
   imputed as `magnitude_rank = 1.0` (worst rank, zero priority)
   and `cellphone_padj = 1.0` (worst adj p, not significant). This
   produces well-defined zero-NA per-contrast tibbles and carries a
   clean semantic ("below liana's detection threshold in the other
   genotype").
7. **Cosmetic-only `n_cells_per_genotype` ordering in provenance.**
   `table(as.character(sc$genotype))` sorts alphabetically, so the
   current cache's `provenance$n_cells_per_genotype` reads
   MAPTKI / NLGF_MAPTKI / NLGF_P301S / P301S instead of the
   canonical MAPTKI / P301S / NLGF_MAPTKI / NLGF_P301S. Patched the
   build script to use
   `table(factor(..., levels = .LIANA_GENOTYPES))` so future
   rebuilds carry the canonical order. The current cache remains
   correct (names are present; downstream code looks up by name);
   the fix only changes the printed log and a `str()` view.
8. **No Ensembl translation needed (deviation-register item (3)
   from the stub NOT triggered).** liana's `mouseconsensus`
   resource is mouse-symbol-keyed; the project's
   `seurat_full_processed.rds` is mouse-symbol-keyed end-to-end;
   no translation step is required for the G3 rbind. Documented in
   the module header.
9. **Stub language "KSN-style attrition stats" not literally
   applicable.** The stub asked the build-script log to "record
   runtime and KSN-style attrition stats" by analogy with F2; liana
   has no kinase-substrate-network analogue, but the analogous
   attrition signal IS captured: per-genotype LR-pair counts (out
   of the 3989-pair mouseconsensus resource; ~22-26% retained per
   genotype after the `expr_prop` + `min_cells` filters in the
   smoke-test subsets), plus the cellphone_padj<0.10 sig counts per
   genotype, all printed by `build_liana_per_genotype_list()` in
   verbose mode and captured in the build-script log.

**Sub-decision-gate at start (not Phase-level).** Resolve the open
question on liana R-vs-Python implementation logged at G1. Default
investigation steps: (a) check container for an R `liana` package
install — if present and working, evaluate whether it has been
recently maintained (last release date, last commit, GitHub issue
backlog); (b) if R `liana` is end-of-life, evaluate the SOTA Python
`liana+` (van Dijk / Dimitrov *et al.* 2024) via reticulate; (c)
if both paths are viable (e.g. R `liana` is in maintenance but
functional), use `AskUserQuestion` to confirm before committing.
Alternative considered at G1 plan time: implement consensus
aggregation manually in pure R by running additional single-method
R CCC tools (iTALK, SingleCellSignalR-R) alongside CellChat +
MultiNicheNet; this is most work but pure-R and is the fallback if
liana itself is unavailable.

**Inputs (post sub-gate).** `storage/cache/seurat_full_processed.rds`
(the symbol-keyed Seurat object used by `scripts/build_cellchat.R`
and `scripts/build_multinichenet.R`); the project-wide 9-cell-type
partitioning (`Astrocyte / Microglia_{DAM,IFN,homeostatic,proliferative}
/ Neuronal / OPC / Oligodendrocyte / Vascular`); the canonical
five-contrast factorial 2×2 + batch design used everywhere else in
the project.

**Action.** Add `R/ccc_inference.R` (new) with the chosen liana
implementation + a thin extraction wrapper around the chosen liana
output schema. Add `scripts/build_liana_cache.R` for the heavy
compute outside the knit (idempotent via `--overwrite`; mirror the
F2 build-script header / detailed-comment-block convention; record
runtime and KSN-style attrition stats in the script log). Insert
`source("R/ccc_inference.R")` in `R/helpers.R` in dependency order
(after `R/kinase_inference.R`).

The cache shape MUST mirror the existing CCC tools' per-contrast
TSV outputs so cross-tool aggregation at G3 is a simple
list-rbind operation. Target per-contrast tibble columns: `sender`,
`receiver`, `ligand`, `receptor`, `lr_interaction`,
`prioritization_score` (or liana's analogous combined score),
`n_methods_pass_significance` (multi-method consensus column),
`sign_dir`. Five canonical contrasts per the project convention;
row count likely 10s of K rows per contrast pre-filter.

**Outputs.** `R/ccc_inference.R` (new); `R/helpers.R` (edited;
source line inserted after `R/kinase_inference.R`);
`scripts/build_liana_cache.R` (new);
`storage/cache/liana_output.rds` (new).

**Verification.** Pre-build smoke test confirming the chosen liana
implementation loads and runs on a tiny subset; post-build smoke
test confirming the cache shape matches the existing MultiNicheNet
output's per-contrast schema (so G3 can rbind across tools); cache
keys present for all 5 canonical contrasts; zero NA scores in
prioritization columns. No knit required at G2 (no Rmd consumer
yet, mirroring F2). Chown all new files `rstudio:rstudio`.

**Plan-spec deviation register.** Watch for: (1) liana+ Python via
reticulate may force Python-side package installation in the
container; record the install dependencies in the script header
for reproducibility. (2) Cross-tool sender / receiver naming
conventions may differ (e.g. liana's `cellchat-consensus` uses a
slightly different prior than CellChatDB.mouse); resolve at G3
during the rbind step, not by editing the per-tool caches. (3) If
liana's output uses gene Ensembl IDs rather than symbols, the
post-filter at G3 needs an Ensembl → symbol translation step (use
`snrnaseq_symbol_map.rds` as at E4); document the convention in
the cache module header. (4) The MultiNicheNet output is keyed by
the project's semantic contrast names (e.g. `tau_alone`); liana's
contrast handling may require the raw factorial-2×2 design matrix;
ensure the wrapper maps to the semantic contrast names for cross-
tool consistency.

---

### Session G3: per-contrast cross-tool LR leader board + 16.1 subsection [DONE 2026-05-24]

**Completion note.** Wrote the first subsection (16.1) of the new
section 16 in `analysis.html` (file `rmd/15_ccc_mechanism.Rmd` renders
as section 16 per the +1 file-index-to-section-number drift; the
parent inserts `15_ccc_mechanism` after `14_kinase_inference` which is
section 15). Knit completed cleanly (~3 min) with zero
`class="error"` and zero `class="warning"` in `analysis.html`. The
new section comprises:

* Section 16 "Cell-cell communication mechanism inference" header +
  intro prose cross-referencing sections 14 (TF) + 15 (kinase) and the
  D2 three-axis verdict; explicit framing of why the synaptic-
  suppression axis 2 non-finding at TF and kinase layers motivates
  the G1 three-tool triangulation decision.
* Subsection 16.1 "Cross-tool L-R activity ranking" with one summary
  table at the top + four level-3 sub-subsections:
  * 16.1.1 Per-contrast cross-tool top-10 L-R pairs (5 tables).
  * 16.1.2 Unified L-R pair leader board (single table; top 25 rows
    of the 194-row leader board displayed; full set in TSV).
  * 16.1.3 L-R pair × tool consistency heatmap at the contrast with
    the most leaders (chosen programmatically; `tau_in_nlgf` here,
    97 leaders).
  * 16.1.4 Interpretive notes (cell-type pair as part of the key;
    structural ceiling of the all-3-sig arm; CellChat
    metabolite-aware ligand limitation; per-axis biology preview
    spanning the three D2 axes).

**Per-tool sig conventions locked.** Three-tool conventions
documented in the 16.1 prose and the `R/ccc_inference.R` G3 module
header; the pre-build diagnostic
(`scripts/g3_prebuild_diagnostic.R`) characterised each tool's
per-contrast significance population at alpha 0.05 / 0.10 / 0.15:

```
CellChat  pval < 0.10 in contrast's primary genotype
          (subsetCommunication's native filter at 0.05 is
           already implicit; permutation null does MTC):
            tau_alone        n_sig = 702  /  702 pairs
            nlgf_in_maptki   n_sig = 1070 / 1070
            nlgf_in_p301s    n_sig = 1141 / 1141
            tau_in_nlgf      n_sig = 1141 / 1141
            interaction      n_sig = 1141 / 1141  (binary inclusion
              vote: pair exists in CellChat's primary-genotype output)

MNN       BOTH scaled_p_val_ligand_adapted < alpha AND
          scaled_p_val_receptor_adapted < alpha:
            tau_alone        n_sig = 143  /  41879 pairs
            nlgf_in_maptki   n_sig = 429  /  41879
            nlgf_in_p301s    n_sig = 741  /  41879
            tau_in_nlgf      n_sig = 560  /  41879
            interaction      n_sig = 206  /  41879  (the most
              selective tool; both-ends-DE filter is MNN-recommended)

LIANA+    cellphone_padj_primary < alpha (BH within per-genotype run):
            tau_alone        n_sig = 4095 / 6838 pairs
            nlgf_in_maptki   n_sig = 6141 / 9755
            nlgf_in_p301s    n_sig = 6763 / 10446
            tau_in_nlgf      n_sig = 6763 / 11220
            interaction      n_sig = 6763 / 11416
```

**Leader rule populated cleanly.** At FDR<0.10, the locked rule
`n_tools_sig_consistent_sign >= 2 OR n_tools_sig >= 3` yields 46-97
leader cells per contrast (vs the F3 emergency floor of zero cells):

```
contrast        n_pairs_any_tool  n_3tools_sig  n_2tools_consistent  n_leaders
tau_alone               44986              0              46               46
nlgf_in_maptki          46636              1              73               73
nlgf_in_p301s           46705              0              94               94
tau_in_nlgf             47158              0              97               97
interaction             47232              0              72               72
```

The all-3-tools-sig arm is essentially dead at FDR<0.10 (0-1 cells per
contrast) because the three tools' significance universes intersect
sparsely. The leader rule's disjunction with the consistent-sign-≥2
arm is what populates the section; this is documented as a structural
property in 16.1.4 prose. The cross-contrast leader board after
aggregation has **194 unique (sender, receiver, ligand, receptor)
cells** at the locked rule across the five contrasts, with
`Sema6a_Plxna4` (OPC → Neuronal axon guidance) as the SOLE 5/5-contrast
leader and 14 other LR cells reaching 4/5 contrasts.

**CellChat duplicate-key bug found and fixed during smoke testing.**
The original smoke test surfaced anomalous `n_tools_sig=3 /
n_tools_present=3` values for several `Glu-SLC1A2_GLS_GRIK4` cells in
the per-contrast top, but the per-tool detail showed `sig_cellchat =
TRUE, sig_mnn = FALSE, sig_liana = FALSE` -- impossible if the per-cell
key is uniquely associated with the three tools. Tracing this surfaced
that CellChat's `subsetCommunication()` returns up to 3 rows for the
same `(sender, receiver, ligand, receptor)` cell when the LR pair has
multiple `interaction_name_2` downstream-complex annotations (e.g.
`Glu-SLC1A2-GLS_GRIK4_GRIA1` vs `Glu-SLC1A2-GLS_GRIK4_GRIA2`). My
extractor's receptor-complex explosion stripped the underscore but
didn't dedupe across these complex-annotation variants. Fix added at
`.extract_lr_cellchat()`: dedupe by `(sender, receiver, ligand,
receptor)` after explosion, keeping the row with the strongest CellChat
evidence (max prob, lowest pval as tiebreak) per the CellChat-vignette
convention for LR-pair-level aggregation over `interaction_name_2`.
Pre-fix diagnostic `n_3tools_sig` counts were inflated to 16-17 per
contrast; post-fix true values are 0-1 per contrast. The leader rule
still populates abundantly via the consistent-sign-≥2 arm; the bug fix
strengthens rather than weakens the leader board because the
inflated cells were spurious cross-tool consensus calls.

**Plan-spec deviations explicitly acknowledged.**

1. **Leader rule structural form: OR not AND.** The plan stub
   specified `n_tools_sig >= 2 AND n_tools_consistent_sign_ge2 >= 2`
   (AND); the implementation uses the E3-parallel
   `n_tools_sig_consistent_sign >= 2 OR n_tools_sig >= 3` (OR). The
   two forms numerically converge in this data because a pair can
   only have consistent sign in ≥2 tools if it is also sig in ≥2
   tools (sig-with-direction implies sig), so the OR's first clause
   is a superset of the AND. The all-3-sig clause adds 0-1 cells per
   contrast on top. Choosing OR keeps the rule notationally identical
   to E3's `n_modalities_sig_consistent_sign >= 2 OR n_modalities_sig
   >= 3` -- the parallel readability is the value-add. Documented in
   the `R/ccc_inference.R` G3 module header and the 16.1 prose.

2. **CellChat dedupe in extractor (post-hoc bug fix, not in stub).**
   Anticipated by the stub's "receptor-complex handling" note for
   receptor strings; not anticipated for the
   `interaction_name_2`-driven duplicate-key case. Documented inline
   at `.extract_lr_cellchat()` with a CellChat-vignette-referenced
   dedupe rationale.

3. **All-3-sig arm is structurally dead (post-fix).** The original
   plan stub's leader-rule design assumed the all-3-sig arm would
   contribute meaningfully (it would catch cells where all three
   tools agree). After the dedupe fix, the diagnostic's pre-locked
   `n_3tools_sig = 16-17 per contrast` collapses to `0-1 per
   contrast`. The 16.1.4 interpretive prose acknowledges this
   explicitly: leaders here are 2-of-3-tool agreements with
   consistent direction, not all-tool unanimous calls. This is the
   right framing -- treating the three tools as independent
   evidence-providers with their own significance universes makes
   2-of-3-agreement the meaningful cross-tool consensus signal.

4. **Score normalisation for cross-tool magnitude comparison.** The
   plan stub did not specify how to compare magnitudes across the
   three tools (each uses a different score scale: CellChat delta
   prob in approx [-1, +1], MNN prioritization_score in [0, 1],
   LIANA+ prioritization_score in approx [-1, +1]). My
   implementation normalises each tool's per-contrast distribution
   to [-1, +1] by dividing by the tool's per-contrast max(|score|),
   then averages over the significant per-tool cells of each LR pair
   (mirror of E3's `mean_abs_score` over the sig subset).
   `mean_norm_score` is the tiebreak metric after
   `n_tools_sig_consistent_sign` and `n_tools_sig`; per-tool raw and
   normalised scores both retained in the long form for
   transparency. Documented in the G3 module header.

5. **TSV column set: 14 columns (vs the plan stub's narrower
   sketch).** The TSV at
   `storage/results/ccc_lr_cross_tool_leaderboard.tsv` has 14
   columns: `lr_interaction`, `sender`, `receiver`, `ligand`,
   `receptor`, `n_contrasts_leader`,
   `n_contrasts_consistent_sign_ge2`, `n_contrasts_sig_ge3`,
   `max_consistent_sign`, `max_n_tools_sig`, `max_mean_norm_score`,
   `dominant_sign`, `contrasts_summary`, `leader_score`. The schema
   mirrors `tf_activity_unified_leaderboard.tsv` exactly with the
   tool-axis terminology substituted for the modality-axis
   terminology. The per-contrast detail is packed into
   `contrasts_summary` (pipe-delimited
   `contrast:n_tools_sig/sign_consensus` per leader contrast),
   matching the TF / pathway TSV convention.

6. **Heatmap contrast chosen programmatically (`which.max` on
   `n_leaders`).** The plan stub said "at the contrast with the most
   leaders"; the implementation picks `tau_in_nlgf` (97 leaders)
   automatically. The figure caption states the contrast.

**Biology highlights.**

1. **`Sema6a_Plxna4` (OPC → Neuronal axon guidance) is the SOLE
   5/5-contrast leader.** Cross-contrast direction is `mixed` because
   `nlgf_in_maptki:2/+`, `tau_alone:2/+`, `nlgf_in_p301s:2/-`,
   `tau_in_nlgf:2/-`, `interaction:2/-` -- the LR pair flips
   direction between amyloid-only and tau-on-amyloid backgrounds,
   exactly the divergent-context signal the project is designed to
   detect. Sema6a/Plxna4 is a canonical OPC-to-neuron axon-guidance
   signal; its position at the very top of the cross-contrast
   leader board indicates intercellular axon-guidance signalling is
   the most-consistently-reshaped LR axis across all five contrasts.

2. **Axis 2 (synaptic suppression) gets multiple cross-tool-supported
   candidates here**, in contrast to the TF and kinase non-findings:
   `Nrxn2_Nlgn1` at 4/5 contrasts (Astrocyte → Neuronal /
   Astrocyte / OPC; neurexin-neuroligin synaptic adhesion);
   `Cntn1_Nrcam` at 4/5 (axon-glia adhesion); `Ncam1_Ncam2` at 4/5
   (canonical homophilic neural adhesion); `Sema6a_Plxna4` and
   `Ntn1_Unc5c` (axon guidance) all in the top 15. This is the
   first mechanism-layer evidence supporting the engulfment / synaptic
   adhesion hypothesis raised at the TF verdict (14.3.3) and kinase
   verdict (15.3.3) -- vindicates the G1 user-directed triangulation
   decision that motivated adding LIANA+ as a third tool. The
   axis-restricted analysis at G4 (16.2) will adjudicate which of
   these candidates carry the strongest axis-2-gene-universe overlap.

3. **NLGF contrasts surface canonical DAM-program LR pairs** at the
   per-contrast top: `Igf1_Igf1r` (Microglia_DAM → Vascular IGF-1),
   `Pros1_Mertk` (Microglia_IFN ↔ Microglia substates; canonical TAM-
   kinase phagocytic axis), `Apoe_Lrp8` (Microglia_DAM → Vascular;
   the AD GWAS lipoprotein-receptor pair), `Adam10_Tspan14` (ADAM10
   is the α-secretase competing with γ-secretase on APP processing;
   surfaced via tetraspanin co-factor), `Tgfb1_Sdc2` (microglia →
   neuronal TGFβ). These map directly to the amyloid-activation axis
   (axis 1) candidates expected by the section 14.3.2 / 15.3.2 TF /
   kinase verdicts; the LR layer adds cell-type-resolved intercellular
   detail to the same biology.

4. **`tau_alone` (no amyloid) surfaces neuronal-glial adhesion + axon
   guidance** at its top: `L1cam_Cd9` (neuron → microglia / oligo /
   OPC; tetraspanin); `Fam3c_Lamp1` (neuron → microglia); `Ptn_Ptprs`
   (OPC → vascular). These are tau-pathology-driven cell-cell
   signalling reshaping in the absence of amyloid, complementing the
   E5/F5 tau-only findings at the TF and kinase layers.

5. **The CellChat metabolite-aware ligand encoding has no cross-tool
   analogue.** Pairs like `Glu-SLC1A2_GLS_GRIK4` (glutamate-glutaminase
   to kainate receptor) and `TENM2_FLRT3` (teneurin-2 to FLRT3) appear
   only in CellChat's output because MNN's NicheNet prior and LIANA+'s
   mouseconsensus resource use bare gene-symbol ligand names. These
   compound-ligand pairs never reach the cross-tool consensus arm at
   FDR<0.10. The known limitation is recorded in 16.1.4; section 11
   carries the full CellChat output for these biologies.

**Anti-anchoring discipline check.** All three D2 axes named in the
16.1.4 interpretive notes with their respective leader candidates
(axis 1 amyloid-activation: Apoe_Lrp8, Igf1_Igf1r, Pros1_Mertk,
Adam10_Tspan14; axis 2 synaptic-suppression: Nrxn1/2_Nlgn1, Cntn1_Nrcam,
L1cam_Cd9, Sema6a_Plxna4; axis 3 interaction metabolic: mixed-sign
Wnt + basement-membrane signals); per-axis cell-type-pair detail
preserved; per-contrast top-10 tables show the full range of per-
contrast biology before axis filtering; the strict-sign-consensus
convention preserves the "mixed" honesty rather than collapsing into
a single sign. The decision-gate-driven triangulation choice (G1
locked) is vindicated empirically at axis 2 even though the all-3-sig
arm is dead -- the consistent-sign-≥2 arm is the right rule
populator.

**File outputs.** All chowned `rstudio:rstudio`:

* `R/ccc_inference.R` (extended; +1,189 lines for the G3 module-
  header comment block, three tool-specific extractors, public
  dispatcher, ranking helper, leader-board builder, table
  formatter, and heatmap plotter; G2 helpers untouched).
* `scripts/g3_prebuild_diagnostic.R` (existing; chowned now that it
  has been re-used and validated by this session).
* `rmd/15_ccc_mechanism.Rmd` (new; 18 KB) -- section 16 + subsection
  16.1 content. All `knitr::kable()` calls inside `results = 'asis'`
  chunks per CLAUDE.md. The build chunk reads the three caches and
  writes the leader-board TSV.
* `analysis.Rmd` (edited; +3 lines for the `child-ccc-mechanism`
  chunk inserted between `child-kinase-inference` and
  `child-session`).
* `storage/results/ccc_lr_cross_tool_leaderboard.tsv` (new; 23 KB;
  194 leader rows + header; 14 columns mirroring TF leader board
  shape).
* `analysis.html` (regenerated; 28.8 MB) -- section 16 + subsection
  16.1 + four sub-subsections render cleanly with zero errors /
  warnings. Anchor IDs `cell-cell-communication-mechanism-inference`,
  `cross-tool-l-r-activity-ranking`,
  `per-contrast-cross-tool-top-10-l-r-pairs`,
  `unified-l-r-pair-leader-board`,
  `l-r-pair-x-tool-consistency-heatmap`, and
  `interpretive-notes-for-this-subsection-3` (auto-disambiguated
  from the 14.1.4 / 14.2.4 / 15.1.3 anchors) all present in the
  rendered TOC at the expected numbers (16, 16.1, 16.1.1 .. 16.1.4).

**Verification.** Pre-knit smoke tests ran twice -- the first surfaced
the CellChat duplicate-key bug; the second (post-fix) returned the
expected per-contrast n_leader counts (46/73/94/97/72), the unified
leader board nrow = 194, dominant-sign distribution 129+ / 26- / 39
mixed / 0 NA, and all the biology highlights above. Pre-knit timing:
~50s (most of which is cache loading); the helpers run in <2s. Post-
knit verification: `grep -c 'class="error"' analysis.html` = 0;
`grep -c 'class="warning"' analysis.html` = 0; new section anchors
all present; TSV at expected path with 195 lines; chown clean.

---

### Session G3 (original stub, retained for audit) [SUPERSEDED]

**Inputs.** G2 liana cache `storage/cache/liana_output.rds`;
existing `storage/cache/cellchat_per_condition.rds`,
`storage/cache/cellchat_merged.rds`,
`storage/cache/multinichenet_output.rds`; D1 leader-board TSV
`storage/results/pathway_survey_unified_leaderboard.tsv`; per-axis
gene-set caches already loaded at E4 / F4
(`hdwgcna_microglia.rds`, custom curated sets, etc.).

**Action.** Extend `R/ccc_inference.R` with E3-parallel helpers:

* `extract_lr_per_tool()` (internal dispatcher) — per-tool,
  per-contrast tidy frame of `(sender, receiver, ligand, receptor,
  lr_interaction, prioritization_score, p_value, padj, sign_dir)`.
  Three tool implementations:
  `.extract_lr_cellchat()`, `.extract_lr_multinichenet()`,
  `.extract_lr_liana()`; concatenated via list-rbind into a single
  long tibble with a `tool` column.
* `rank_lr_cross_tool()` — mirror `rank_tfs_cross_modality()`. Per
  (contrast, LR-pair) cell, count significant tools (FDR<0.10 on
  the tool-native significance statistic; CellChat → communication
  probability p-value; MultiNicheNet → prioritization score
  percentile and/or `prioritization_padj` where available; liana →
  combined p-value from cellchat-consensus or analogous). Track
  consistent-sign count across tools (where sign is the direction
  of upregulation in the contrast's primary group).
* `build_lr_cross_tool_leaderboard()` — one row per (contrast,
  LR-pair); columns: `contrast`, `lr_interaction`, `sender`,
  `receiver`, `ligand`, `receptor`, `n_tools_sig`,
  `n_tools_consistent_sign_ge2`, `mean_prioritization_score`,
  `dominant_sign`, `tools_summary` (pipe-delimited per-tool detail
  string); leader rule `n_tools_sig >= 2 AND
  n_tools_consistent_sign_ge2 >= 2`. Mirror E3's
  `(>=2 modalities sig consistent OR all 3 sig)` rule in spirit;
  the exact AND/OR structure should be set after the G3 pre-build
  diagnostic shows actual cell counts (see deviation register).

Create `rmd/15_ccc_mechanism.Rmd` with section header
`# Cell-cell communication mechanism inference` (renders as section
16). Subsection 16.1 "Cross-tool L-R activity ranking" with four
level-3 sub-subsections (mirroring E3's 14.1 shape):

* 16.1.1 Per-contrast cross-tool top-N LR pairs (5 tables; include
  explicit honest-non-finding prose for any contrast yielding
  fewer than the requested N at the rule).
* 16.1.2 Unified LR pair leader board (single table; cap at 25
  rows per leader rule output).
* 16.1.3 LR pair × tool consistency heatmap at the contrast with
  the most leaders; rows alphabetical so no LR pair is visually
  privileged; `*` = FDR<0.05, `.` = FDR<0.10; NA cells grey
  (mirror E3 / F3 visual convention).
* 16.1.4 Interpretive notes (cross-tool consistency stats; honest
  non-findings; cross-reference to section 14.1 TF leaderboard and
  section 15.1 per-contrast kinase tables; bulk-vs-single-cell
  caveat for CellChat).

Add `child-ccc-mechanism` chunk in `analysis.Rmd` between
`child-kinase-inference` and `child-session`. Export the leader
board to `storage/results/ccc_lr_cross_tool_leaderboard.tsv`
(mirror `tf_activity_unified_leaderboard.tsv` shape: one row per
leader, per-contrast detail packed into `tools_summary` as a
pipe-delimited string).

**Outputs.** `R/ccc_inference.R` (extended);
`rmd/15_ccc_mechanism.Rmd` (new); `analysis.Rmd` (edited; +3 lines
for the child chunk);
`storage/results/ccc_lr_cross_tool_leaderboard.tsv` (new);
`analysis.html` (re-knitted; new section 16 + subsection 16.1).

**Verification.** Pre-knit smoke test (~5 s): cross-tool rbind
yields a non-empty per-contrast leader board with all 5 canonical
contrasts represented and at least 5 leader LR pairs surviving the
rule (if not, document in deviation register and consider relaxing
to FDR<0.15 only with rationale — mirror F3 deviation register
guidance). Post-knit: `grep -c 'class="error"' analysis.html` = 0;
`grep -c 'class="warning"' analysis.html` = 0; section 16 + 16.1
+ four sub-subsections render at expected positions; TSV present
at expected path. Chown all new files `rstudio:rstudio`.

**Plan-spec deviation register.** Watch for: (1) The F3 user-
directed deviation abandoned the cross-cache leader-board pattern
because the cross-cache consistency arm yielded zero cells at
FDR<0.10. The G3 cross-tool consistency arm faces a similar risk;
if the pre-build diagnostic shows zero / near-zero cross-tool
consensus cells at FDR<0.10, document the diagnostic in the same
shape as F3's, surface to user via `AskUserQuestion` for a
deviation decision, and either (a) drop to per-tool independent
significance (mirror of F3), (b) relax to a softer cross-tool rule
(e.g. ≥2/3 tools at FDR<0.15), or (c) drop the consensus arm
entirely and treat tools as independent inference units. The F3
precedent suggests honesty over rule-rescuing. (2) Three-tool
consensus rules are more flexible than two-cache rules (the rule
can require 2/3 or 3/3 etc.), so the rule-rescue space is larger
than at F3 — the diagnostic data should drive the choice.

---

### Session G4: axis-restricted L-R analysis + 16.2 subsection [DONE 2026-05-24]

**Completion note.** Built the axis-restricted L-R analysis layer
cleanly against the G3 in-memory `ccc_ranking_long` tibble. Four new
helpers appended to `R/ccc_inference.R`
(`restrict_lr_to_universe`, `score_lr_per_axis`,
`format_axis_restricted_lr_table`, `plot_axis_lollipop_lr`), each
mirroring the E4 / F4 schema exactly so the 16.2 / 14.2 / 15.2
axis-restricted tables read side-by-side with parallel mental
models. A 16,672-row hybrid TSV was exported to
`storage/results/ccc_lr_axis_restricted.tsv` (6,074 amyloid + 6,459
synaptic + 4,139 interaction LR cells; 14 columns including the
F4-locked `score_at_<contrast>` per-contrast hybrid trio). Subsection
16.2 added to `rmd/15_ccc_mechanism.Rmd` with the planned 4-subsection
structure (16.2 build chunk + 16.2.1 amyloid_activation + 16.2.2
synaptic_suppression + 16.2.3 interaction_metabolic + 16.2.4
cross-axis observations); each axis sub-subsection renders a
`format_axis_restricted_lr_table()` kable (`results = 'asis'`) and a
`plot_axis_lollipop_lr()` figure. Knit completed in ~3 min with zero
`class="error"` and zero `class="warning"` in `analysis.html`.

**Biology bottom line (per axis, three-line summary).**

```
amyloid_activation     6,074 cells; top 5: Adam11_Itga4 / Gas6_Axl /
                       L1cam_Egfr / Omg_Rtn4r / Rspo3_Lgr4 (all + )
                       canonical microglia: Cd200_Cd200r1 N->MG_DAM (rank 6),
                                            Apoe_Trem2 MG_DAM->MG_DAM (rank 9),
                                            App_Cd74 OL->MG (rank 11),
                                            Tgfb1_Sdc2 MG_DAM->N (rank 34)
                       55 / 100 top-100 microglia-involving
                       => CORROBORATES amyloid-driven microglia CCC signature

synaptic_suppression   6,459 cells; top 5: Adam11_Itga4 / Ntn4_Unc5a /
                       L1cam_Egfr / Reln_Vldlr / Lin7c_Htr2c (all + )
                       49 / 100 top-100 microglia-involving
                       canonical engulfment ranks deeply:
                         Cd47_Sirpa (rank 279) -- "don't eat me" inhibitor
                         C1qb_Lrp1  (rank 775) -- first classical complement
                         no C3 / Mertk-Pros1 / Tyro3-Axl in top 800
                       surfaced mechanism instead: TREM2-mediated clearance
                                 (Apoe_Trem2 rank 9, App_Cd74 rank 12) +
                                 synaptic adhesion modulation
                                 (Cd200_Cd200r1, Vcan_Cd44, Ncam1 ligands)
                       => TRIANGULATION HYPOTHESIS VINDICATED at presence-
                          and-rank level but the surfaced mechanism is
                          BROADER than classical complement-mediated pruning

interaction_metabolic  4,139 cells; top 5: Gpc3_Unc5c / Entpd1_Adora1 /
                       L1cam_Ephb2 / Gas6_Axl / Nrxn1_Nlgn1
                       mixed-sign: top 25 = 18 + / 7 -
                       canonical microglia at top: Efna5_Epha4 MG_h->N (rank 8),
                                                   Apoe_Lrp5 N->MG (rank 9),
                                                   L1cam_Cd9 N->MG (rank 13, - )
                       => CORROBORATES section-13 axis-3 mixed-sign reading;
                          NEW BIOLOGY: ephrin / axon-guidance REWIRED at
                          interaction (some pairs up, some down) rather than
                          uniformly suppressed
```

**Helper design choices documented in the module-level header.**

* `mode = "either"` (ligand OR receptor in axis universe) locked
  default; `mode = "both"` (stricter) documented as a robustness
  check, surfaces ~50% fewer cells per axis without changing the
  top-of-table biology. The plan's deviation register foresaw this
  sub-decision; locked the more permissive option per the plan's
  default proposal, with the alternative captured in code and prose.
* `min_n_tools_present = 2L` locked default. The G3 ranking_long
  contains 87% single-tool rows (MNN-only, since MNN's
  `group_prioritization_tbl` is 4-40x larger than LIANA+'s and
  CellChat's per-genotype universes); without this filter the axis
  rankings would be dominated by MNN-only cells whose
  `mean_norm_score` reduces to MNN's always-positive
  `prioritization_score`. Mirror of E4 / F4's `min_targets = 5L`
  baseline data-density filter.
* `score_col = "mean_norm_score"` (G3's cross-tool coalesce of sig
  and present means). Signed direction inherited from 16.1; magnitudes
  comparable across tools via per-(tool, contrast) max-abs
  normalisation.
* Per-pair universe column named `n_lr_ends_in_axis_universe` (range
  {1, 2}) rather than `n_targets_in_axis_universe` (E4 / F4 column
  name) so the downstream TSV reader sees the semantic shift and
  does not mis-interpret high-target TFs as comparable to two-end LR
  pairs. This is the only schema column-name divergence between the
  16.2 TSV and the 14.2 / 15.2 TSVs; all other columns (axis, source-
  like keys, mean / sd / n_cells / leader_rank / score_at_<contrast>)
  use parallel names.

**Plan-spec deviations explicitly acknowledged.**

1. **Source key is four columns, not one.** E4 / F4's
   `source` column (single string for TF / kinase) is replaced by
   `(sender, receiver, ligand, receptor)` plus the convenience
   `lr_interaction` derived column. This is a structural property of
   the L-R layer (the same biochemistry can appear at multiple cell-
   type interfaces, each a distinct biological signal) and was
   already baked into the G3 leader board (G3 used the same four-key
   pattern). The display kable and lollipop chart concatenate
   `lr_interaction | sender -> receiver` for the y-axis labels so
   distinct cell-typed cells of the same biochemistry are
   distinguishable. Schema parity with 14.2 / 15.2 is preserved at
   every other column; this divergence is the minimum necessary to
   honour the L-R layer's true unit of analysis.

2. **`mode` plan sub-decision resolved to `"either"`.** The plan's
   G4 deviation register foresaw the `mode` choice ("either" vs
   "both") as a sub-decision and proposed defaulting to "either"
   (parallel to E3's "TF significant in any modality" permissive
   logic). Locked that default after smoke-testing both: "either"
   gives 3-6k cells per axis at `min_n_tools_present = 2L`, "both"
   gives 0.4-4k. Top-of-table biology is dominated by the same
   ~20 canonical AD-microglia and synaptic pairs at both modes; the
   stricter "both" mode is preserved in code as a robustness check
   for future sessions but not exported as a separate TSV.

3. **`min_n_tools_present = 2L` locked as the cross-tool noise
   filter.** Not in the plan stub's enumeration of locked decisions,
   but the smoke test made it essential: without this filter, every
   axis top-N would be MNN-only positive cells (MNN universe = 42k
   cells per contrast; LIANA+ = 10k; CellChat = 1k; intersection
   sparse). Documented at length in the `R/ccc_inference.R` G4
   module header with the rationale. This is the L-R analogue of
   E4's `min_targets = 5L` baseline filter (E4's prevents TFs with
   too few in-universe targets; G4's prevents LR cells with too few
   tools observing them). Code allows `min_n_tools_present = 1L`
   for explicit relaxation, but no analysis in 16.2 uses it.

4. **The axis-2 finding rate is a NUANCED positive, not a unanimous
   one.** The plan stub said "if axis 2 yields >0 LR pairs at the
   `"either"` filter, that IS the engulfment-mediated mechanism the
   E5 / F5 verdicts predicted". The finding is: 6,459 axis-2 LR cells
   at the filter; the top-10 carry strong microglia-involving signal
   (Apoe_Trem2 rank 9, Astrocyte→MG_DAM Apoe_Trem2 rank 10,
   App_Cd74 rank 12); 49 of top-100 involve a microglia sender or
   receiver. So the layer DOES surface axis-2 biology -- the
   triangulation hypothesis is vindicated at the presence-and-rank
   level. BUT: the canonical complement-mediated synaptic pruning
   pairs (C1q*, C3, Mertk-Pros1) do NOT rank highly: first complement
   pair (`C1qb_Lrp1` Microglia_homeostatic → Neuronal) at rank 775;
   first regulator-of-engulfment (`Cd47_Sirpa`) at rank 279. The
   surfaced mechanism is therefore TREM2-mediated clearance +
   APP-fragment uptake + synaptic adhesion modulation, not classical
   complement-driven synaptic pruning. The 16.2.2 prose and the 16.2.4
   cross-axis observations make this nuance explicit; the 16.3
   verdict will need to read axis 2 as
   "TREM2-and-adhesion-mediated clearance" rather than as
   "complement-mediated pruning". Both are plausible post-
   transcriptional, post-kinase mechanisms consistent with the E5 /
   F5 prediction.

5. **G3 cross-contrast leaders rank deeply at the axis layer.** A
   property worth flagging for the 16.3 verdict that did not appear
   in the plan stub: the G3 leader board (cross-contrast aggregation)
   and the 16.2 axis layer (axis-contrast magnitude) are
   complementary views, not redundant. `Sema6a_Plxna4` (SOLE 5/5
   G3 leader; OPC → Neuronal at G3) ranks 5,882 at axis 2 in that
   direction (its best axis-2 placement is rank 158 for the
   Neuronal → Vascular cell-typed cell of the same biochemistry).
   `Psap_Gpr37l1` (4/5 G3 leader) tops out at rank 1,662 at axis 1
   with all-negative axis-1 scores; `Nrxn2_Nlgn1` (4/5 G3 leader)
   tops out at rank 5,445 at axis 2. The G3 leader rule weights
   cross-CONTRAST consistency (leader at multiple contrasts) while
   the axis layer weights axis-CONTRAST MAGNITUDE (top-of-normalised-
   distribution score at the axis's 1-2 contrasts). Pairs with
   modest per-contrast magnitudes but consistent presence dominate
   G3; pairs with strong axis-contrast magnitudes dominate the axis
   layer. Documented in 16.2.4 observation 3; the 16.3 verdict will
   integrate both.

6. **`n_cells_used = 1` cells outscore `n_cells_used = 2` cells at
   axes 1 and 2.** Surfaced during smoke testing: pairs observed at
   only one of the two NLGF contrasts retain top-of-normalised-
   distribution magnitudes (~1.0), while pairs observed at both
   contrasts average against a second observation and dilute toward
   the mean. The top 5 of axis 1 is
   `Adam11_Itga4 (n=1, 1.000) / Gas6_Axl (n=1, 0.998) / L1cam_Egfr (n=2, 0.994) / Omg_Rtn4r (n=2, 0.985) / Rspo3_Lgr4 (n=2, 0.982)`.
   The 16.2.4 prose recommends the 16.3 verdict preferentially weighs
   `n_cells_used = 2` cells (cross-contrast persistence is the stronger
   evidence of axis relevance) and reports `n_cells_used = 1` cells
   with appropriate caveat. Not a bug -- the n=1 high scores are real
   single-contrast signal, just less robust than n=2 cross-contrast
   signal.

7. **Axes 1 and 2 share many top-of-table cells by Interpretation-A-
   parallel construction.** Both NLGF axes use the same axis_contrasts,
   and key AD-microglia genes (APOE, TREM2, APP, CD74, VCAN, CD44)
   appear in both the amyloid-activation and synaptic-suppression
   leader-pathway gene sets, so the universe-overlap filter qualifies
   the same pairs at both axes with identical scores. `Apoe_Trem2`
   Microglia_DAM → Microglia_DAM is rank 9 at BOTH axes with mean
   0.961; `App_Cd74` Oligodendrocyte → Microglia_X is at identical
   ranks across the two axes. The differentiator is the lower-ranked
   cells where the universes diverge (and the `n_lr_ends_in_axis_universe`
   column). Mirror of the E4 14.2.4 "Interpretation-A axis-1/2 shared
   top" property; honest reporting at 16.2.4 rather than smoothing
   it away.

**File outputs.** All chowned `rstudio:rstudio`:

* `R/ccc_inference.R` (extended; 83 KB total) -- appended a new
  G4 section header with a ~120-line module-level comment block
  covering the input contract, score source rationale, cross-tool
  noise filter rationale, hybrid output shape, `mode` semantics, and
  per-pair universe-membership column naming convention. Four new
  helpers: `restrict_lr_to_universe()`,
  `score_lr_per_axis()`, `format_axis_restricted_lr_table()`,
  `plot_axis_lollipop_lr()`. Existing G1..G3 helpers untouched.
* `rmd/15_ccc_mechanism.Rmd` (extended; 35 KB total) -- appended
  subsection 16.2 ("Axis-restricted L-R rankings") with a build
  chunk + three axis sub-subsections (16.2.1 / 16.2.2 / 16.2.3) +
  the cross-axis observations sub-subsection (16.2.4). Build chunk
  reads the in-memory `ccc_ranking_long` (G3 shared global) and
  defensively re-builds `axis_universes` if not in scope. Each axis
  sub-subsection has a `format_axis_restricted_lr_table()` chunk
  (in `results = 'asis'`) and a `plot_axis_lollipop_lr()` figure
  chunk. Inline prose is axis-specific and surfaces canonical AD-
  microglia LR pairs by name with their leader_rank in parentheses.
* `storage/results/ccc_lr_axis_restricted.tsv` (new; 2.4 MB) --
  long-format: 14 columns (axis, sender, receiver, ligand, receptor,
  lr_interaction, mean_activity_in_axis_contrasts,
  sd_activity_in_axis_contrasts, n_cells_used,
  n_lr_ends_in_axis_universe, leader_rank, score_at_nlgf_in_maptki,
  score_at_nlgf_in_p301s, score_at_interaction). 16,672 rows
  (6,074 amyloid + 6,459 synaptic + 4,139 interaction).
* `analysis.html` (regenerated; 28.0 → 29.4 MB) -- subsections 16.2 +
  16.2.1 / 16.2.2 / 16.2.3 / 16.2.4 render cleanly with zero errors /
  warnings. Anchor IDs `axis-restricted-l-r-rankings`,
  `amyloid-activation-axis-axis-1-4`,
  `synaptic-suppression-axis-axis-2-4`,
  `interaction-metabolictranslational-axis-axis-3-4`, and
  `cross-axis-observations-for-this-subsection` (auto-disambiguated
  from the prior axis-sub-subsection anchors at 14.2, 14.3, 15.2,
  15.3, 16.1) all present in the rendered TOC at the expected
  section numbers (16.2, 16.2.1, 16.2.2, 16.2.3, 16.2.4).
* `storage/notes/mechanism_layer_plan.md` (this file; edited).

**Verification.** Pre-knit smoke test on the live caches (~3 s)
confirmed: ccc_lr_long produces 264k rows / 47k distinct LR cells
across 5 contrasts and 3 tools as expected; ranking_long has 87%
single-tool rows -- justifying the `min_n_tools_present = 2L` filter;
score_lr_per_axis produces 16,672 axis-cell rows in 1.8 s; the
top-5 of each axis has tool diversity (no MNN-only pseudo-leaders);
mode="both" is strictly smaller (7,866 rows) than mode="either"
(16,672 rows); per-contrast `score_at_<contrast>` columns populate
non-null only for axis-relevant contrasts; format / plot helpers
return correct kable / ggplot objects. Post-knit verification:
`grep -c 'class="error"' analysis.html` returns 0;
`grep -c 'class="warning"' analysis.html` returns 0; TSV file on
disk at the expected path and 16,672 rows; section TOC shows 16.2 +
four level-3 sub-subsections at 16.2.1 / 16.2.2 / 16.2.3 / 16.2.4.

**Phase G4 complete.** Phase G5 (CCC mechanism verdict subsection
16.3) is now the next active session, unblocked. The verdict will
read axes 1, 2, and 3 against the 14.3 TF verdict and 15.3 kinase
verdict, integrating the axis-2 nuance (TREM2-mediated clearance +
adhesion modulation rather than complement-mediated pruning) and the
axis-3 sign rewiring (ephrin / axon-guidance) into the project's
most complete mechanism reading.

---

### Session G4 (original stub, retained for audit) [SUPERSEDED]

**Inputs.** G3 cross-tool LR leader board TSV; G3 in-memory
`ccc_lr_long` long-form tibble (shared knit scope per project
convention; built earlier in `rmd/15_ccc_mechanism.Rmd` by the G3
chunks); per-axis gene universes from E4 / F4 (re-use
`build_axis_gene_universe()` from `R/tf_inference.R`); leader-board
TSV.

**Action.** Extend `R/ccc_inference.R` with E4 / F4-parallel
helpers for axis-restricted LR analysis. The axis-restricted view
applies the per-axis gene universe as a SENDER + RECEIVER filter:
an LR pair qualifies for axis X if the ligand AND/OR the receptor
(configurable via `mode` argument) intersects the axis-X gene
universe.

* `restrict_lr_to_universe(lr_tbl, axis_universe, mode = "either")`
  — filter LR pairs where ligand AND/OR receptor (per `mode` =
  `"either"` / `"both"`) is in `axis_universe`.
* `score_lr_per_axis(lr_long, axis_universes, ...)` — mirror
  `score_tf_per_axis()` / `score_kinase_per_axis()` schema exactly.
  Hybrid axis-mean + per-axis-contrast columns per the F4 locked
  shape; rank LR pairs within each axis by mean prioritization
  score across axis-relevant contrasts.
* `format_axis_restricted_lr_table()` — mirror
  `format_axis_restricted_table()` and
  `format_axis_restricted_kinase_table()`.
* `plot_axis_lollipop_lr()` — mirror `plot_axis_lollipop()` and
  `plot_axis_lollipop_kinase()`.

Append subsection 16.2 to `rmd/15_ccc_mechanism.Rmd` with the
locked 4-subsection structure (build chunk + 16.2.1
amyloid_activation + 16.2.2 synaptic_suppression + 16.2.3
interaction_metabolic + 16.2.4 cross-axis observations). Each axis
sub-subsection renders a `format_axis_restricted_lr_table()` kable
(`results = 'asis'`) and a `plot_axis_lollipop_lr()` figure.

Export axis-restricted long-form to
`storage/results/ccc_lr_axis_restricted.tsv` (mirror
`tf_activity_axis_restricted.tsv` /
`kinase_activity_axis_restricted.tsv` shape, plus per-axis-contrast
hybrid columns).

**Outputs.** `R/ccc_inference.R` (extended);
`rmd/15_ccc_mechanism.Rmd` (extended; subsection 16.2 added);
`storage/results/ccc_lr_axis_restricted.tsv` (new);
`analysis.html` (re-knitted, new subsection 16.2).

**Verification.** Pre-knit smoke test: axis universes still
populated (re-use the E4 / F4 sanity bound of ≥5 entries per axis);
each axis carries ≥5 LR pairs at the
`restrict_lr_to_universe(mode = "either")` filter; rank-by-mean
stable across tools (no single tool dominating the top-N). Post-knit:
`grep -c 'class="error"' analysis.html` = 0;
`grep -c 'class="warning"' analysis.html` = 0; subsection 16.2.1 /
16.2.2 / 16.2.3 / 16.2.4 numbered at expected positions. Chown all
new files `rstudio:rstudio`.

**Plan-spec deviation register.** Watch for: (1) the `mode` choice
(`"either"` vs `"both"`) on `restrict_lr_to_universe` is itself a
sub-decision — default to `"either"` (more permissive, parallels
E3's "TF significant in any modality" logic) but document the
alternative `"both"` (stricter, requires both ligand and receptor
in axis universe) in the deviation register. (2) The axis 2
synaptic-suppression axis is the critical test of the triangulation
hypothesis. If axis 2 yields >0 LR pairs at the `"either"` filter,
that IS the engulfment-mediated mechanism the E5 / F5 verdicts
predicted; document axis 2 finding rate explicitly in 16.2.2 and
the deviation register. The user-directed G1 deviation to
triangulation will be "justified by the data" or "honest non-
finding" at this step. (3) The `restrict_lr_to_universe()` filter
operates on gene symbols; if any tool's output is keyed by Ensembl
IDs, the translation step must happen here, not earlier — the G3
leader board should preserve the original IDs.

---

### Session G5: CCC mechanism verdict subsection 16.3 [DONE 2026-05-24]

**Inputs.** G3 cross-tool leader board TSV; G4 axis-restricted
TSV; in-memory `ccc_lr_long` and `ccc_lr_leaderboard` from earlier
chunks in the same `rmd/15_ccc_mechanism.Rmd` (shared knit scope
per project convention).

**Action.** Append subsection 16.3 "Verdict (cross-tool + cross-
contrast LR synthesis)" to `rmd/15_ccc_mechanism.Rmd` mirroring
E5's 14.3 and F5's 15.3 shape. Six level-3 anchors:

* 16.3 top-line + threshold restatement (per-tool FDR<0.10; cross-
  tool rule from G3 final lock; axis-restricted gene universe via
  `restrict_lr_to_universe` mode `"either"` (or G4 final lock);
  top-N per axis; bulk-vs-single-cell caveat for CellChat; OmniPath
  prior-coverage caveat for liana).
* 16.3.1 Verdict summary table (one row per axis; mirror F5
  8-column hybrid schema: `axis`, `top_lr_pairs`, `top_lr_signs`,
  `mean_score_range`, `per_contrast_score_range`,
  `per_contrast_summary`, `n_top_sig_in_axis_contrasts` (cross-
  tool consistency count from G3 leader board), `evidence_summary`).
* 16.3.2 Amyloid-activation axis (axis 1) verdict + biology;
  cross-reference TF axis-1 ensemble (Spi1 / Nfkb1 / Sp3) from F5
  and kinase axis-1 ensemble (Cdk2 / Cdk5 / Csnk1e / Mapk8 /
  Camk2g) from F5.
* 16.3.3 Synaptic-suppression axis (axis 2) verdict + biology —
  the critical test of the triangulation hypothesis. If positive,
  this is THE mechanism-layer answer for axis 2 and resolves the
  E5 / F5 non-finding; if honest non-finding, document explicitly
  and add a post-G5 retrospective paragraph on what the negative
  triangulation result implies for the overall mechanism-arc
  reading (e.g. engulfment is operating at a layer none of the
  three CCC tools can resolve, such as phagosomal trafficking).
* 16.3.4 Interaction metabolic / translational axis (axis 3)
  verdict + biology + cross-reference to TF axis-3 Myc / Creb1 /
  Tp53 / Jun ensemble (F5 finding) and kinase axis-3 Gsk3b /
  Csnk1a1 / Mapk14 / Cdk5 / Cdk1 ensemble (F5 finding). The joint
  TF + kinase + LR view at axis 3 will be the project's most-
  integrated mechanism reading.
* 16.3.5 Cross-axis observations + limitations: sign reversals;
  honest non-findings; cross-reference to TF verdict 14.3 and
  kinase verdict 15.3 by anchor; the project's full mechanism-arc
  reading per axis (TF + kinase + LR) integrated in prose.

Add `build_lr_verdict_table()` helper to `R/ccc_inference.R`
mirroring `build_tf_verdict_table()` and
`build_kinase_verdict_table()` (thin aggregator; one tibble per
axis; top-N filter; hybrid columns; cross-tool consistency count
via the G3 leader-board's `n_tools_sig` / `tools_summary` columns).

**Outputs.** `R/ccc_inference.R` (extended);
`rmd/15_ccc_mechanism.Rmd` (extended; subsection 16.3 added);
`storage/results/ccc_lr_verdict.tsv` (new; three rows, eight
columns; mirror F5 shape);
`analysis.html` (re-knitted).

**Verification.** Pre-knit smoke test: `build_lr_verdict_table()`
returns a three-row data.frame with the expected eight columns;
ordering natural amyloid / synaptic / interaction (not
alphabetised). Post-knit: zero errors / warnings; verdict table
renders as a real HTML table (caption + colgroup + tbody cells
visible in the html); subsection 16.3 + sub-subsections 16.3.1 ..
16.3.5 all numbered in the rendered TOC. Chown all new files
`rstudio:rstudio`.

**Anti-anchoring discipline check (per E5 / F5 convention).** All
three D2 axes named explicitly in the verdict; asymmetric reading
(strong / data-dependent / strong, depending on G4 axis-2 result)
honest rather than synthetic; per-axis sub-subsections equal-length
even where one is a non-finding; inference thresholds restated
before applying them; no Hallmark / no OXPHOS-arc terminology;
hdWGCNA modules alphabetised wherever they appear in multi-module
lists. Additional G-phase rule: when reading axis 2, present the
triangulation result honestly — "engulfment mechanism surfaced" if
positive; "triangulation did not surface synaptic suppression
drivers either, suggesting the mechanism operates at a layer none
of the three CCC tools can resolve" if negative, with a one-
paragraph retrospective on what the negative G5 result implies
for the overall mechanism-arc reading.

**Phase G COMPLETE marker.** G5 completion note should explicitly
declare Phase G COMPLETE and update the plan's top-of-file
`State:` line to mark G1..G5 all DONE. The mechanism arc (Phases
E + F + G) is then complete; the plan should be moved to
`storage/notes/completed/mechanism_layer_plan_<YYYY-MM-DD>.md`
per the "How to mark this plan complete" instruction at the end
of this file, with a final `## Outcome summary` section recording
(per phase that ran): the inference tool selected at the decision
gate, the top 5 drivers / kinases / LR pairs surfaced per axis,
and any cross-axis surprises.

**Completion note (2026-05-24).** Wrote the cross-tool +
cross-contrast LR verdict subsection (16.3) at
`rmd/15_ccc_mechanism.Rmd` with five level-3 sub-subsections, plus
the new `build_lr_verdict_table()` helper at
`R/ccc_inference.R` (file grew from 1,838 -> 2,069 lines). The
helper is a thin aggregator mirroring `build_kinase_verdict_table()`
and `build_tf_verdict_table()`: one row per axis, top-N filter,
F4/E5/F5-locked hybrid columns (axis-mean + per-axis-contrast
min/max + per-contrast summary string), with the LR-layer
adaptations for the four-key (sender, receiver, ligand, receptor)
unit of analysis. `top_lr_pairs` is `<lr_interaction>
(<sender>-><receiver>)` semicolon-separated -- the cell-type pair
is essential at the L-R layer because the same biochemistry (e.g.
`Apoe_Trem2`) at `Microglia_DAM->Microglia_DAM` vs
`Astrocyte->Microglia_DAM` are distinct biological signals and
collapsing to `lr_interaction` alone would lose a real dimension
of the signal. `n_top_sig_in_axis_contrasts` is implemented by
matching the top-N four-keys against the G3 cross-tool leader
board and parsing the leader board's `contrasts_summary` column
to check membership of any axis-relevant contrast in the leader
list. Verdict TSV at `storage/results/ccc_lr_verdict.tsv` (three
rows, eight columns; mirror of F5 shape with `top_lr_pairs` /
`top_lr_signs` substituted for `top_kinases` / `top_kinase_signs`).
Knit completed in ~2 min with zero `class="error"` and zero
`class="warning"` in `analysis.html` (29.45 MB; up from 29.4 MB
post-G4). Section 16.3 and sub-subsections 16.3.1 .. 16.3.5
render at the expected positions in the TOC (anchors
`verdict-cross-tool-cross-contrast-lr-synthesis`,
`verdict-summary-table-2`,
`amyloid-activation-axis-axis-1-5`,
`synaptic-suppression-axis-axis-2-5`,
`interaction-metabolictranslational-axis-axis-3-5`,
`cross-axis-observations-and-limitations-3`,
auto-disambiguated from the corresponding anchors at
sections 14.3, 15.3, and earlier in section 16).

**Verdict bottom line (per-axis, three-line summary).**

```
amyloid_activation     top 5: Adam11_Itga4 (Neuronal->Neuronal),
                              Gas6_Axl (Neuronal->Astrocyte),
                              L1cam_Egfr (Neuronal->Vascular),
                              Omg_Rtn4r (Astrocyte->Neuronal),
                              Rspo3_Lgr4 (Neuronal->Vascular)
                       signs: +,+,+,+,+    range [+0.98, +1.00]
                       canonical DAM-program LR at ranks 6-25:
                         Cd200_Cd200r1 Neuronal->MG_DAM (rank 6),
                         Apoe_Trem2    MG_DAM->MG_DAM    (rank 9),
                         App_Cd74      OL->MG            (rank 11),
                         Tgfb1_Sdc2    MG_DAM->Neuronal  (rank 34)
                       55 / 100 top-100 microglia-involving
                       0 / 5 in G3 LB at axis-relevant contrasts
                       => CORROBORATES section-13 axis 1, TF axis 1
                          (Spi1 / Nfkb1 / Sp3), and kinase axis 1
                          (Cdk2 / Cdk5 / Mapk8 / Csnk1e) -- LR adds
                          cell-type-resolved intercellular detail

synaptic_suppression   top 5: Adam11_Itga4 (Neuronal->Neuronal),
                              Ntn4_Unc5a   (Vascular->Neuronal),
                              L1cam_Egfr   (Neuronal->Vascular),
                              Reln_Vldlr   (Neuronal->Neuronal),
                              Lin7c_Htr2c  (Neuronal->Neuronal)
                       signs: +,+,+,+,+    range [+0.99, +1.00]
                       49 / 100 top-100 microglia-involving
                       Apoe_Trem2  MG_DAM->MG_DAM         (rank  9),
                       Apoe_Trem2  Astrocyte->MG_DAM      (rank 10),
                       App_Cd74    OL->MG                  (rank 12)
                       classical complement ranks deeply:
                         C1qb_Lrp1 MG_h->Neuronal         (rank 775)
                         Cd47_Sirpa                        (rank 279)
                         no C3 / Mertk-Pros1 / Tyro3-Axl in top 800
                       0 / 5 in G3 LB at axis-relevant contrasts
                       => G1 TRIANGULATION HYPOTHESIS VINDICATED at
                          presence-and-rank level (the LR layer DOES
                          surface axis-2 biology where TF + kinase
                          returned non-findings), but the SURFACED
                          MECHANISM is TREM2-mediated clearance +
                          APP-fragment uptake + synaptic-adhesion
                          modulation rather than classical complement-
                          mediated pruning. Honest non-finding ->
                          refined positive at this layer; the
                          mechanism arc reads axis 2 meaningfully.

interaction_metabolic  top 5: Gpc3_Unc5c    (Vascular->Astrocyte),
                              Entpd1_Adora1 (OPC->OPC),
                              L1cam_Ephb2   (Neuronal->Neuronal),
                              Gas6_Axl      (Neuronal->Astrocyte),
                              Nrxn1_Nlgn1   (Neuronal->Neuronal)
                       signs: +,-,-,+,+    range [-0.54, +0.56]
                       canonical microglia at broader top:
                         Efna5_Epha4 MG_h->Neuronal (rank 8, +),
                         L1cam_Cd9   Neuronal->MG   (rank 13, -)
                       => CORROBORATES section-13 axis-3 mixed-sign
                          reading; ephrin / axon-guidance REWIRED
                          (some pairs up, some down) at tau x amyloid.
                          Joint three-layer mechanism at axis 3:
                          TF Myc / kinase Gsk3b / LR ephrin family
                          -- the project's most-integrated mechanism
                          reading.
```

**Plan-spec deviations explicitly acknowledged.**

1. **`top_lr_pairs` cell format extends the F5 / E5 convention.**
   F5 / E5 emit `top_kinases` / `top_TFs` as simple comma-separated
   gene-symbol lists. At the L-R layer the four-key
   `(sender, receiver, ligand, receptor)` is the true unit of
   analysis -- the same biochemistry at different cell-type cells
   carries distinct biological signal -- so collapsing to
   `lr_interaction` alone (the closest E5/F5 parallel) would lose a
   real dimension. G5 emits `<lr_interaction>
   (<sender>-><receiver>)` with `;` as the outer separator
   (semicolon used because cell-type names like `Microglia_DAM`
   contain commas in the wild and might confuse a comma-outer-
   separator parser; semicolon avoids any ambiguity). Schema parity
   with F5 / E5 is preserved at every other column. Documented in
   the `R/ccc_inference.R` G5 module header.

2. **`n_top_sig_in_axis_contrasts = 0/0/0` at top-5 is the honest
   reading, not a bug.** The G3 cross-contrast leader board and
   the G4 axis-magnitude axis layer weight complementary aspects
   of the data: G3 rewards pairs with modest scores observed
   consistently across multiple contrasts (e.g. `Sema6a_Plxna4`,
   SOLE 5/5-contrast G3 leader); the G4 axis layer rewards pairs
   with top-of-normalised-distribution scores at the axis's 1-2
   contrasts (often `n_cells_used = 1` single-contrast cells, e.g.
   `Adam11_Itga4` axis-1 rank 1 mean 1.000). Quantified: of the
   top-100 axis-1 cells, 0 appear in the 194-row G3 leader board;
   at top-1000, 6 appear. Across all three axes the same property
   holds. The 0/0/0 column in the verdict is therefore the honest
   reading of the cross-tool consistency lens; the verdict prose
   treats `n_top_sig_in_axis_contrasts > 0` as additional cross-
   tool corroboration rather than as the dominant evidence signal,
   and the 16.3.5 cross-axis observations subsection documents the
   complementarity explicitly. Foreseen at G4 deviation register
   entry #5 and propagated cleanly to G5.

3. **Inputs read from TSV rather than from in-memory G3 globals.**
   The G5 stub listed in-memory `ccc_lr_long` and
   `ccc_lr_leaderboard` as inputs (the G3 build chunk creates them
   earlier in the same Rmd via the parent knit's shared R session).
   The actual build chunk reads
   `storage/results/ccc_lr_axis_restricted.tsv` (G4 output) and
   `storage/results/ccc_lr_cross_tool_leaderboard.tsv` (G3 output)
   via `readr::read_tsv()` instead. This matches the F5 pattern
   (which reads `kinase_activity_axis_restricted.tsv` +
   `kinase_activity_per_contrast.tsv` from disk rather than
   in-memory globals) and makes the build chunk independent of
   knit-ordering / out-of-knit re-runs. Functionally identical to
   the in-memory path; the design choice is shape consistency with
   F5 + Rmd robustness. The plan stub's "in-memory" wording is the
   slip; the actual implementation uses the more defensive
   read-from-disk path. The `lr_leaderboard` argument name in
   `build_lr_verdict_table()` is consistent with the disk-read
   variable name.

4. **TSV column set is 8 columns, mirror of F5 not E5.** E5's
   verdict has `n_top_in_cross_modality_leaderboard` (a TF cross-
   modality leader-rule count); F5's has
   `n_top_sig_in_axis_contrasts` (a kinase per-contrast FDR count).
   G5 mirrors F5's column name and semantic, with the cross-tool
   leader rule substituted for the per-contrast FDR gate. The
   choice is documented in the `R/ccc_inference.R` G5 module
   header and reflects the natural L-R analogue (cross-tool
   consensus at LR cells is the LR-layer equivalent of cross-
   modality consensus at TFs and per-contrast FDR at kinases; G3's
   `n_tools_sig_consistent_sign >= 2L OR n_tools_sig >= 3L` is the
   cross-tool leader rule and is the natural count target for G5).

5. **Axis-2 reading is a refined positive, not a vanilla positive.**
   The plan stub framed axis 2 as a binary test: "if positive, this
   is THE mechanism-layer answer for axis 2 and resolves the E5 /
   F5 non-finding; if honest non-finding, document explicitly and
   add a post-G5 retrospective paragraph on what the negative
   triangulation result implies". The actual axis-2 reading is
   neither pure: the LR layer DOES surface axis-2 biology
   (vindicating the G1 triangulation decision), but the surfaced
   mechanism (TREM2-mediated clearance + APP-fragment uptake +
   synaptic-adhesion modulation) is BROADER than the literature-
   canonical classical-complement-pruning hypothesis the prediction
   centred on. The 16.3.3 subsection makes this nuance explicit and
   the 16.3.5 cross-axis subsection includes a "post-G5
   retrospective" paragraph that the plan stub anticipated for the
   negative case but applies to the refined-positive case here. The
   refined-positive framing preserves the plan's intent (be honest
   about what triangulation surfaced) without forcing a binary
   positive/negative reading.

**Anti-anchoring discipline check.** All three D2 axes are named
explicitly in the verdict; the per-axis sub-subsections are
equal-length even where one (axis 2) is a refined positive rather
than a clean strong-positive; the inference thresholds (G3 cross-
tool leader rule, G4 universe filter, top-5 cut) are restated
explicitly before applying them in the 16.3 intro; no Hallmark / no
OXPHOS-confirmatory-arc terminology appears; hdWGCNA modules
alphabetised wherever they appear in multi-module lists. The
G-phase additional rule "present the triangulation result honestly"
is satisfied via the 16.3.3 + 16.3.5 prose explicitly naming
TREM2-mediated clearance + adhesion modulation as the surfaced
mechanism and explicitly recording that the canonical complement
pathway is NOT the dominant mechanism in these data; the post-G5
retrospective at the end of 16.3.5 frames this as a hypothesis-
generating refinement rather than as a refutation of the engulfment
hypothesis.

**File outputs.** All chowned `rstudio:rstudio`:

* `R/ccc_inference.R` (extended; 1,838 -> 2,069 lines = +231 lines)
  -- appended the new G5 section header with a ~85-line module-
  level comment block covering input contract, schema rationale,
  expected count behaviour with the complementary-views caveat,
  and per-column documentation. Single new helper
  `build_lr_verdict_table()` (~120 lines including doc + body).
  Existing G1..G4 helpers untouched.
* `rmd/15_ccc_mechanism.Rmd` (extended; 719 -> 1,279 lines = +560
  lines) -- appended subsection 16.3 with five level-3
  sub-subsections (16.3.1 verdict summary, 16.3.2 amyloid-
  activation axis, 16.3.3 synaptic-suppression axis, 16.3.4
  interaction metabolic / translational axis, 16.3.5 cross-axis
  observations and limitations). The build chunk reads the G3 +
  G4 TSVs from disk, builds the verdict, exports the verdict TSV,
  and prints the verdict-summary kable in a `results = 'asis'`
  chunk per CLAUDE.md. All `evidence_summaries` editorial prose
  lives in the Rmd's named-list assignment, not in the helper, so
  future editorial edits don't require code changes (mirror of
  F5 / E5 pattern).
* `storage/results/ccc_lr_verdict.tsv` (new; 7,595 B; three rows +
  header = four lines on disk) -- eight columns: `axis`,
  `top_lr_pairs`, `top_lr_signs`, `mean_score_range`,
  `per_contrast_score_range`, `per_contrast_summary`,
  `n_top_sig_in_axis_contrasts`, `evidence_summary`.
* `analysis.html` (regenerated; 29.45 MB; up from G4's 29.4 MB) --
  subsection 16.3 + five sub-subsections (16.3.1 .. 16.3.5) render
  cleanly with zero errors / warnings; anchor IDs all present at
  expected TOC positions.
* `storage/notes/mechanism_layer_plan.md` (this file; edited).

**Verification.** Pre-knit smoke test on the live caches (~2 s)
exercised `build_lr_verdict_table()` against the G3 + G4 TSVs:
3-row x 8-col data.frame returned; axis order natural amyloid /
synaptic / interaction (not alphabetised); top_lr_pairs / signs
match the G4 expected ranking byte-for-byte; per_contrast_summary
strings render at the F5 / E5 format `"contrast:[%+0.2f, %+0.2f]"`;
n_top_sig_in_axis_contrasts integer values 0 / 0 / 0 cross-validated
against direct G3 leader-board lookup at axis-relevant contrasts
(only Nrxn1_Nlgn1 N->N is in LB across all 15 top-5-per-axis
cells, and it is a leader at NLGF contrasts rather than at axis
3's interaction contrast, so the verdict count of 0 is correct).
Post-knit verification: `grep -c 'class="error"' analysis.html` =
0; `grep -c 'class="warning"' analysis.html` = 0; section 16.3 +
sub-subsections 16.3.1 .. 16.3.5 all numbered in rendered HTML;
verdict-summary kable renders as a real HTML table; TSV file at
expected path with expected shape; all new / modified files
chowned `rstudio:rstudio`.

**Phase G COMPLETE.** All five G sessions DONE (G1 decision gate
locking three-tool triangulation + leader-pathway-restricted post-
filter + new `rmd/15_ccc_mechanism.Rmd` venue; G2 LIANA+ via
reticulate cache build; G3 cross-tool LR leader board + 16.1 + 194-
row leader-board TSV; G4 axis-restricted LR analysis + 16.2 + 16,672-
row axis-restricted TSV; G5 verdict subsection 16.3 + helper + TSV +
three-layer integrated mechanism reading). **The mechanism arc
(Phases E + F + G) is now complete.** The plan should be moved to
`storage/notes/completed/mechanism_layer_plan_2026-05-24.md` per
the "How to mark this plan complete" instruction at the end of
this file, with the `## Outcome summary` section recording the
phase-by-phase verdict.

---

## Outcome summary

### Phase E (transcription factor inference): tool, top-5 drivers per axis, surprises

* **Inference tool selected at E1 decision gate:** `decoupleR` 2.16.0 +
  CollecTRI mouse prior (`decoupleR::get_collectri(organism = "mouse",
  split_complexes = FALSE)`, 40,291 -> 39,961 interactions after
  dedup; 1,114 unique TF sources, 6,110 unique targets). Alternative
  considered and declined: DoRothEA-only (now a subset of CollecTRI);
  pyscenic *de novo* regulon discovery (overkill given curated priors
  + DE stats already exist). `ulm` arbitrates significance at FDR<0.10;
  `consensus` arbitrates magnitude / sign.
* **Axis 1 (amyloid_activation) top 5 TFs:** A0A979HLR9 (+4.32), Spi1
  (+2.81), Rreb1 (-2.57), Sp3 (+2.50), Nfkb1 (+2.47). Two of five
  (A0A979HLR9, Spi1) also reach the cross-modality leader rule.
  Biology: PU.1 / NF-kB / Sp3 canonical inflammatory amplifiers
  activated at NLGF onset; Rreb1 lone repressor.
* **Axis 2 (synaptic_suppression) top 5 TFs:** identical to axis 1 by
  Interpretation A construction (both NLGF axes share the same
  contrast set; the per-TF mean is computed from full-universe scores).
  Honest non-finding -- the TF layer does NOT independently identify
  drivers of synaptic suppression; the differentiator is in-universe
  target counts (e.g. Spi1: 184 amyloid-universe targets vs 84
  synaptic-universe targets).
* **Axis 3 (interaction_metabolic) top 5 TFs:** A0A979HLR9 (-4.85), Myc
  (-3.25), Creb1 (-2.28), Tp53 (-2.15), Jun (-1.81). All five uniformly
  NEGATIVE; signed range [-4.85, -1.81]. One of five (A0A979HLR9) also
  reaches the cross-modality leader rule. Biology: Myc-led coordinated
  suppression of biosynthetic activity at the tau x amyloid interface.
  RESOLVES the section-13 axis-3 mixed-sign reading -- proteomics +
  phospho push negative through Myc / Creb1 / Tp53 / Jun suppression of
  ribosomal biogenesis and biosynthetic throughput.
* **Cross-axis surprises at Phase E:** (1) Three TFs show sign reversal
  between NLGF axes (+) and interaction axis (-): A0A979HLR9, A0A087WSP5,
  and **Rela** (NF-kB p65). Rela is the most interpretable: NF-kB p65
  is ACTIVATED by amyloid alone and SUPPRESSED by tau on the amyloid
  background -- the signature a tau-attenuates-amyloid-inflammation
  model would predict, and a candidate hypothesis for follow-up. (2)
  The TF mechanism layer reads as ASYMMETRIC across the three axes
  (strong / non-finding / strong); the asymmetry is the honest
  reading rather than a missed result.

### Phase F (kinase activity inference): tool, top-5 drivers per axis, surprises

* **Inference tool selected at F1 decision gate:** `decoupleR` 2.16.0 +
  OmniPath KSN prior (`decoupleR::get_ksn_omnipath()` aggregating
  PhosphoSitePlus + SIGNOR + multiple databases; mouse mapping via
  `nichenetr::convert_human_to_mouse_symbols` on the substrate-symbol
  column since regulatory phosphosites are highly conserved). Same
  framework as Phase E so the F5 verdict reads directly against the
  E5 verdict using identical conventions. Alternatives KSEA-only and
  decoupleR + KSEA cross-tool triangulation considered and declined
  (cost vs benefit unfavourable at this stage of the mechanism arc).
* **Axis 1 (amyloid_activation) top 5 kinases:** Cdk2 (+2.91), Cdk5
  (+1.75), Csnk1e (+1.64), Mapk8 (+1.17), Camk2g (-1.16). Two of five
  (Cdk5, Csnk1e) reach FDR<0.10 on ulm at nlgf_in_p301s. Biology:
  proliferation + stress-kinase ensemble downstream of Spi1 / Nfkb1
  TF activation. Hybrid per-contrast view exposes Cdk5 / Mapk8 as
  p301s-biased (maptki +0.50 / +0.25 vs p301s +3.00 / +2.08) -- a
  real tau-on-amyloid substructure the axis-mean smooths over.
* **Axis 2 (synaptic_suppression) top 5 kinases:** identical to axis 1
  by Interpretation A construction. Honest non-finding. Structural
  differentiator is in-universe substrate-site counts: Cdk5's KSN
  footprint expands from 81 amyloid-universe sites to 263 synaptic-
  universe sites (3.2x increase, reflecting its documented role in
  synaptic phospho-regulation), and similar for Cdk2 / Csnk1e /
  Camk2g. Synaptic axis is best interpreted as a target-overlap
  filter ("which amyloid-active kinases have substrate programs that
  also touch the synaptic machinery"), not as a list of kinases that
  drive synaptic suppression.
* **Axis 3 (interaction_metabolic) top 5 kinases:** **Gsk3b (+4.77)**,
  Csnk1a1 (-2.81), Mapk14 (+2.63), Cdk5 (+2.10), Cdk1 (+1.99). Two of
  five (Gsk3b, Csnk1a1) reach FDR<0.10 on ulm at the interaction
  contrast. **Gsk3b at padj 0.002 is the strongest single-kinase
  result of the entire mechanism arc.** Per-contrast pattern across
  all five contrasts: -1.66 / -1.47 / +3.79 / +4.15 / +4.32 -- amyloid
  alone and tau alone both SUPPRESS Gsk3b substrate phospho; the
  combination drives the activation, peaking at the interaction
  contrast. **The synergy IS the Gsk3b activation** -- the cleanest
  mechanism-layer answer to the project's interaction question.
* **Cross-axis surprises at Phase F:** (1) **Cdk5 is the integrator
  across all three axes** -- the only kinase to surface in the top 5
  at axes 1, 2, AND 3, with its KSN substrate footprint expanding
  with the axis (81 sites at amyloid universe; 263 at synaptic
  universe; 153 at interaction). Cdk5's literature reputation as a
  master regulator of both tau phospho-biology and synaptic
  phospho-biology is consistent with this cross-axis pattern. (2)
  The kinase mechanism layer reads asymmetrically (strong / non-
  finding / strong) in the same shape as the TF layer at Phase E;
  two mechanism layers AGREE that the synaptic-suppression axis is
  not a regulator-driven program at the resolution available to
  either, jointly motivating the G1 decision to add LIANA+ as a
  third CCC tool at Phase G.

### Phase G (CCC bridge re-analysis with leader-pathway filters): tool, top-5 drivers per axis, surprises

* **CCC re-analysis approach locked at G1 decision gate:** **three-
  tool triangulation** of CellChat (existing) + MultiNicheNet
  (existing) + LIANA+ v1.7.1 (new via reticulate at G2), with
  leader-pathway-restricted gene-universe post-filter applied
  uniformly to all three. User-directed deviation from the plan's
  default post-filter-only proposal; rationale was that the
  synaptic-suppression axis was an explicit non-finding at the TF
  and kinase layers and adding tools to read an axis that prior
  tools could not read is the natural response. Cross-tool leader
  rule (G3 lock): `n_tools_sig_consistent_sign >= 2L OR n_tools_sig
  >= 3L` at per-tool FDR<0.10 (CellChat permutation `pval < 0.10`
  in primary genotype; MNN BOTH ligand AND receptor adapted-p <
  0.10; LIANA+ `cellphone_padj_primary < 0.10` BH within per-
  genotype run).
* **CCC mechanism Rmd venue (G1 lock):** new `rmd/15_ccc_mechanism.Rmd`
  rendering as section 16. Reader sees TF (section 14) + kinase
  (section 15) + CCC mechanism (section 16) grouped as a coherent
  mechanism arc in the TOC; the descriptive CCC home at section 11
  (subsections a-f) is unchanged.
* **Axis 1 (amyloid_activation) top 5 LR cells:** Adam11_Itga4
  (Neuronal->Neuronal, +1.00), Gas6_Axl (Neuronal->Astrocyte, +1.00),
  L1cam_Egfr (Neuronal->Vascular, +0.99), Omg_Rtn4r (Astrocyte->
  Neuronal, +0.99), Rspo3_Lgr4 (Neuronal->Vascular, +0.98). 0 of 5
  in G3 leader board at axis-relevant contrasts (complementary
  views property). Top 5 dominated by Neuronal-sender axon-guidance
  + adhesion biology; canonical DAM-program LR pairs (Cd200_Cd200r1,
  Apoe_Trem2, App_Cd74) accessible at ranks 6-25; 55/100 top-100
  microglia-involving. Adds cell-type-resolved intercellular detail
  to the same DAM-program biology surfaced at TF axis 1 (Spi1 /
  Nfkb1 / Sp3) and kinase axis 1 (Cdk2 / Cdk5 / Mapk8 / Csnk1e).
* **Axis 2 (synaptic_suppression) top 5 LR cells:** Adam11_Itga4
  (Neuronal->Neuronal, +1.00), Ntn4_Unc5a (Vascular->Neuronal, +1.00),
  L1cam_Egfr (Neuronal->Vascular, +0.99), Reln_Vldlr (Neuronal->
  Neuronal, +0.99), Lin7c_Htr2c (Neuronal->Neuronal, +0.99). Top 5
  similar in shape to axis 1 (Interpretation A entanglement
  property propagates from TF / kinase to LR). 0 of 5 in G3 LB at
  axis-relevant contrasts. **THE G1 TRIANGULATION HYPOTHESIS IS
  VINDICATED**: 6,459 axis-2 LR cells (the largest axis); Apoe_Trem2
  (MG_DAM -> MG_DAM) at rank 9; Apoe_Trem2 (Astrocyte -> MG_DAM) at
  rank 10; App_Cd74 (Oligodendrocyte -> Microglia) at rank 12;
  49/100 top-100 microglia-involving. **BUT the surfaced mechanism is
  TREM2-mediated clearance + APP-fragment uptake + synaptic-adhesion
  modulation, NOT classical complement-mediated pruning** (first
  complement `C1qb_Lrp1` at rank 775; first regulator-of-engulfment
  `Cd47_Sirpa` at rank 279; no `C3` / `Mertk-Pros1` / `Tyro3-Axl` in
  top 800). Refined positive at the LR layer where TF + kinase
  returned honest non-findings -- the mechanism arc reads axis 2
  meaningfully at the LR resolution.
* **Axis 3 (interaction_metabolic) top 5 LR cells:** Gpc3_Unc5c
  (Vascular->Astrocyte, +0.56), Entpd1_Adora1 (OPC->OPC, -0.54),
  L1cam_Ephb2 (Neuronal->Neuronal, -0.50), Gas6_Axl (Neuronal->
  Astrocyte, +0.50), Nrxn1_Nlgn1 (Neuronal->Neuronal, +0.50). Mixed
  signs (3+, 2-); range [-0.54, +0.56]. 0 of 5 in G3 LB at the
  interaction contrast (Nrxn1_Nlgn1 IS in LB but as nlgf_in_p301s +
  tau_in_nlgf leader, not interaction). Cell-type-resolved LR
  rendition of the section-13 axis-3 mixed-sign reading: ephrin /
  axon-guidance / synaptic-adhesion biology is REWIRED at the
  interaction (some pairs up, some down), not uniformly suppressed.
  At broader scope `Efna5_Epha4` (Microglia_homeostatic -> Neuronal
  rank 8, +) and `L1cam_Cd9` (Neuronal -> Microglia rank 13, -)
  carry the canonical microglia ephrin signal.
* **Cross-axis surprises at Phase G:** (1) **Sema6a_Plxna4** (OPC ->
  Neuronal axon-guidance) is the SOLE 5/5-contrast G3 leader, with
  cross-contrast direction MIXED (positive at `tau_alone` and
  `nlgf_in_maptki`; negative at `nlgf_in_p301s`, `tau_in_nlgf`,
  `interaction`) -- the LR pair flips direction between amyloid-only
  and tau-on-amyloid backgrounds, exactly the divergent-context
  signal the project is designed to detect. (2) The all-3-sig arm
  of the G3 leader rule is structurally DEAD at FDR<0.10 (0-1
  cells per contrast); the populating arm is consistent-sign-≥2.
  This is the right framing -- treating the three tools as
  independent evidence-providers with their own significance
  universes makes 2-of-3-agreement the meaningful cross-tool
  consensus signal. (3) The G3 cross-contrast leader board and the
  G4 axis-magnitude axis layer are LARGELY DISJOINT at the top-N
  (0 of top-100 axis-1 cells in the 194-row G3 LB; only 6 of top-
  1000). The two views weight complementary aspects of the data and
  are jointly required for a complete LR mechanism reading. (4) The
  axis-2 LR finding REFINES the literature-canonical complement-
  centric framing of microglia-synapse interaction: in these data
  TREM2-mediated clearance + adhesion modulation, not C1q-C3-Mertk-
  Pros1 complement-tagging, is the dominant intercellular mechanism
  at the synaptic-suppression axis. This is a hypothesis-generating
  refinement for future experimental work, treated as cross-tool-
  cross-cohort-conditional rather than as a universal claim.

### Final cross-axis verdict (mechanism-arc-integrated reading)

The mechanism arc reads each of the section-13 D2 axes with the
following per-axis joint multi-layer mechanism:

* **Axis 1 (amyloid_activation):** Strong corroboration across all
  three mechanism layers. TF lead = Spi1 / Nfkb1 / Sp3 (DAM-program
  inflammatory amplifiers); kinase lead = Cdk2 / Cdk5 / Mapk8 /
  Csnk1e (proliferation + stress-kinase ensemble); LR lead = Cd200/
  Cd200r1 + Apoe/Trem2 + App/Cd74 at ranks 6-12 + a Neuronal-sender
  axon-guidance + adhesion top-5 (Adam11_Itga4 / Gas6_Axl / L1cam_Egfr
  / Omg_Rtn4r / Rspo3_Lgr4). Joint reading: the canonical microglial
  DAM-state induction at NLGF onset, with the transcriptional axis
  driving the program, the kinase axis delivering the phospho-
  signalling effector layer, and the LR axis adding cell-type-
  resolved intercellular adhesion + signalling rewiring.
* **Axis 2 (synaptic_suppression):** Honest non-finding at TF + kinase
  layers; REFINED POSITIVE at LR layer. The LR layer surfaces a
  post-regulator intercellular signalling mechanism (TREM2-mediated
  clearance + APP-fragment uptake + synaptic-adhesion modulation)
  that the TF + kinase layers predicted but could not specify. The
  surfaced mechanism is BROADER than the literature-canonical
  classical-complement-pruning hypothesis the prediction centred on;
  the canonical C1q-C3 complement pathway is NOT the dominant
  mechanism in these data at this resolution. This is the project's
  refined mechanism-layer answer for axis 2 and a hypothesis-
  generating refinement for the literature-canonical engulfment
  framing.
* **Axis 3 (interaction_metabolic):** Strong corroboration across
  all three mechanism layers and the most-integrated three-layer
  mechanism reading of the project. TF lead = Myc / Creb1 / Tp53 /
  Jun (coordinated biosynthetic suppression); kinase lead = Gsk3b
  (the strongest single-kinase result at padj 0.002) / Mapk14 /
  Cdk5 / Cdk1 (canonical tau-kinase activation); LR lead = ephrin /
  axon-guidance / synaptic-adhesion rewiring (Gpc3_Unc5c /
  Entpd1_Adora1 / L1cam_Ephb2 / Gas6_Axl / Nrxn1_Nlgn1 at top 5;
  Efna5_Epha4 + L1cam_Cd9 at broader microglia top). Joint reading:
  at the tau x amyloid interface the cell suppresses biosynthetic
  capacity transcriptionally (TF Myc lead), activates canonical
  tau-kinase phospho-signalling (kinase Gsk3b lead), and rewires
  intercellular axon-guidance / synaptic-adhesion communication
  (LR ephrin family lead) -- all three layers run in mutually
  consistent biological directions and jointly resolve the section-
  13 axis-3 mixed-sign call to a single coherent mechanism.

### Next-plan pointer

No active successor plan at the close of this session. The
mechanism arc is complete; the natural next direction is
user-directed and should be brainstormed at the next session's
launch via the survey-state-and-propose-directions branch of the
launch protocol. Candidate directions surfaced during Phase G
that warrant their own plans if the user chooses to pursue them:

1. **Cdk5-substrate-resolved sub-analysis** (foreshadowed at 15.2.4
   and 15.3.5): test whether the axis-specific Cdk5 substrates
   differ in identity across axes 1, 2, and 3 -- the cross-axis
   substrate-set differential (81 / 263 / 153 sites) suggests
   axis-specific Cdk5 programs are plausible.
2. **Complement-pathway-resolved microglial sub-analysis at axis 2**
   (foreshadowed at 16.3.3 / 16.3.5): use a microglia-only CCC
   tool, microglia-sorted bulk RNA-seq, or experimental complement-
   pathway perturbation to adjudicate whether the surfaced TREM2
   mechanism is the *only* axis-2 mechanism at this cohort or the
   *dominant* one in these three tools' overlap region.
3. **Rela parallel sign-reversal follow-up** (foreshadowed at
   14.2.4 / 14.3.5): Rela (NF-kB p65) is activated by amyloid alone
   and suppressed by tau on the amyloid background -- the signature
   a tau-attenuates-amyloid-inflammation hypothesis would predict.
   A focused sign-reversal characterisation across the three
   mechanism layers (TF Rela / kinase NF-kB pathway / LR NF-kB
   target genes) would test the hypothesis directly.
4. **CollecTRI UniProt-accession resolution** (foreshadowed at
   14.2.4 / 14.3.5): A0A979HLR9 and A0A087WSP5 are CollecTRI
   heterodimer complex sources retained through E5; resolving them
   to specific TF subunits via UniProt would clarify the axis-1 /
   axis-3 sign-reversal biology.

Each of these would be a fresh single-purpose plan with its own
decision-gate session at G1's structural model.

---

## Anti-anchoring guardrails (re-read every session)

These exist because LLMs (this agent included) drift toward the most
salient prior finding. Each session must enforce them:

- **Always** derive mechanism targets from the D2 leader-board TSV,
  never from a pre-committed pathway family.
- **Always** present the mechanism findings across all three axes
  (amyloid activation / NLGF synaptic suppression / mixed-sign
  metabolic-translational) unless the data empirically restricts the
  story to a subset. Collapsing into a single "winner" reintroduces
  bias.
- **Always** alphabetise hdWGCNA modules in heatmaps; never place
  MG-M3 first.
- **Always** state the inference threshold (FDR cutoff, score
  threshold, activity-score sign convention) explicitly before
  applying it.
- **Never** reintroduce Hallmark gene sets at any step. The collection
  is permanently retired from this project.
- **Never** reintroduce the OXPHOS confirmatory arc (deleted Rmds 08,
  09a-09d). The mechanism work must derive its targets from the D2
  leader board.
- **Never** rebuild a specificity-null or expression-matched-random
  framework as part of this plan. A future plan may reintroduce
  nulls if a specific question demands them, but it must be motivated
  by that question.
- **Never** pre-commit to one of the three biological axes at the
  expense of the other two. The leader-rule fix in D1/D2 was
  explicitly designed to keep the three axes visible; the mechanism
  layer must inherit that discipline.

## How to mark this plan complete

When the final session of the last active phase is DONE, move this
file to
`storage/notes/completed/mechanism_layer_plan_<YYYY-MM-DD>.md`,
`chown rstudio:rstudio`, and add a final `## Outcome summary` section
recording (per phase that ran): the inference tool selected at the
decision gate, the top 5 drivers / kinases / L-R pairs surfaced per
axis, and any cross-axis surprises. Write the next-plan pointer (if
any) and close.

# Decision digests (new project)

Compact per-phase decision summaries for the fresh `main` rebuild - read these, not
the full archived plans. v1 digests live in `archive_digest.md` (full v1 plans on
branch `archive`).

## P0 Foundations -- closed 2026-06-29 (-> `.agent/completed/p0_foundations_plan_2026-06-29.md`)

Built the reproducible spine for P1-P5: project-local env, shared pure-fn modules, 4 analysis-ready
modalities, a QC-sanity report, 2x2 factorial + 5-contrast machinery, a concrete quality gate. NO analysis
(P1+). Live operational facts -> memory.md; wiring -> map.md. Below = decisions + rejected alternatives.

- STACK: targets DAG + rv (R) + uv (Python) + project-local Quarto, P3M-pinned (snapshot 2026-06-22). rv
  over renv (agent-native, declarative, mirrors uv; accepts pre-1.0 + hand-wired Bioc 3.23). targets over
  plain Rmd (content-hash skip pays on the 8G object + multi-hour fits). NO Docker -> repro =
  targets(pipeline) + locks(versions) + P3M(pinning), explicitly NO bitwise guarantee. R 4.6 / Bioc 3.23 ->
  RE-BASELINE (numbers differ from v1's R 4.5.2; never reproduce v1's locked margins).
- REPORT: pivoted Quarto BOOK -> ONE standalone offline HTML (`index.qmd` format:html embed-resources +
  `{{< include _qc.qmd >}}`). Book rejected: multi-file nav emits `Could not fetch resource ./<sibling>.html`
  under embed-resources -> caught by the gate's render-log scan (quiet=FALSE). theme.scss = crimson #B0344D + IBM Plex (9 woff2
  base64-inlined offline; woff2 COMMITTED, read-economy skip).
- MODULES: pure fns via `tar_source` (the DAG orders execution -> no manual loader; supersedes v1 helpers.R):
  constants/utils/io/design/de_pb/plot/spine. 6 modalities materialized as qs2 targets. 5 canonical contrasts
  via TWO equivalent parameterisations (factorial `~tau+nlgf+tau_nlgf[+batch]` AND cell-means `~0+genotype`),
  proven equal by a property test. de_pb = limma-voom/log only; edgeR-QLF/DESeq2/dream deferred (KISS). S3 =
  machinery only (no DE results).
- QUALITY GATE: `scripts/check.sh` -- concrete, fail-loud, zero-fault; review-hardened (P0 close).
  tar_make's exit is insufficient (0 even on CAPTURED warnings, and blind to the report's SEPARATE render
  process), so enforcement is layered: qmd (`options(warn=2)` -> any R chunk warning fails the render), target
  (`tar_quarto(quiet=FALSE)` -> Quarto/Pandoc warnings reach the log), script (FORCE-render the report each run
  -> always exercises both; tar_meta all-NA scoped to manifest+branches; anchored render-log grep). Tests =
  stopifnot harness + synthetic fixtures (no testthat dep). Detail -> memory.md.

Verification (honest): each module smoke-tested vs LIVE data at its step (S2 shapes incl. the 8G load, S3
16-col pseudobulk, S4 clean render); `scripts/check.sh` green end-to-end on the populated store (S5). The
one-shot FRESH-CLONE rebuild is the design contract, verified incrementally across S1-S5, NOT re-run at S5
(the heavy 8G load_snrnaseq BUILD, not the cheap report re-render). The gate FORCE-renders the report each run
-> the warn=2 + render-log checks always exercise (no cached-clean blind spot); residual = out-of-band edits to
`_report/` outputs, moot since it re-renders from source. Gate review-hardened post-close.

Deferred -> P1+: microglia reprocess (SCT/Harmony/cluster-prune), substate assignment, single-cell DE
(NEBULA/glmmTMB), pseudobulk DE RESULTS, edgeR-QLF/DESeq2/dream, BPCells (8G RAM relief), col-map DAG
validation (P4). Out of scope (v1 bloat): cisTarget mm10/SCENIC, SEA-AD human validation.

## P1 snRNAseq microglia core -- closed 2026-06-30 (-> `.agent/completed/p1_snrnaseq_plan_2026-06-30.md`)

The first analysis phase: reprocess + integrate + cluster microglia (S1), UCell substates + contaminant prune
(S2), substate composition (S3), pseudobulk DE (S4), `_microglia.qmd` report + close (S5). Built the robust
amyloid->DAM activation headline, supported three complementary ways (composition, DE, UCell score), plus an
HONEST under-powered interaction handed to P2. Live facts -> memory.md (P1-S1..S5 sections); wiring -> map.md.
Below = decisions + rejected alternatives.

- DE METHOD: single-cell DE DROPPED -> pseudobulk limma-voom is the SOLE inference (Squair 2021 / Murphy-Skene
  2022: cell-level DE = pseudoreplication -> FDR inflation; v1's NEBULA interaction was null anyway -> nothing
  lost). voomWithQualityWeights + robust eBayes + stageR family screen (omnibus-F, df1=rank=3 -> modified-Holm
  per-contrast). REJECTED: NEBULA/glmmTMB single-cell, edgeR-QLF/DESeq2/dream (KISS, deferred indefinitely).
- COMPOSITION: REVERSED the plan's sccomp-primary -> propeller-logit LOCKED primary (+asin sensitivity) + sccomp
  OPTIONAL off-lock Bayesian cross-check. WHY: CmdStan is a compiled C++ tree OFF the P3M snapshot -> rv cannot
  lock it -> a Bayesian arm cannot be THE reproducible call. propeller needs a CELL-MEANS design (PropRatio
  raises per-genotype mean coefs to contrast powers; a factorial design crashes/misreads). REJECTED: sccomp-
  primary (unlockable), averaging discordant methods (propeller-logit STANDS, discordance flagged not averaged).
- NORMALISATION: SCTransform v2 (v1 continuity; user choice over the SOTA logNorm lean). Harmony over BATCH ONLY
  (sex perfectly aliased with batch; never integrate over genotype/amyloid -> DAM is biology). REJECTED: logNorm,
  Harmony over batch+sex.
- SUBSTATES: UCell (rank-based, dropout-robust) calibrated PER-SIGNATURE then cluster-argmax (raw UCell not
  cross-signature comparable; cluster-PRIMARY authoritative, per-cell NOISY); 4 substates + aux MHC/APC.
  Contaminant prune on RAW identity-vs-contam (z FAILS -- ambient contam pervasive destroys the absolute "is a
  microglia"); dropped {6,7,8,11}=2944/26104. REJECTED: AddModuleScore (control-bin, unstable on sparse nuclei),
  z-based prune, per-cell substate proportions, over-pruning.
- HEADLINE (robust, microglia-led): amyloid (NLGF) drives homeostatic->DAM. Confirmed 3 ways -- composition
  (propeller DAM-up FDR~1e-10/1e-13, sccomp concordant), DE (DAM markers amyloid-UP frac 1.00/0.94, meanLFC
  +1.37/+1.85 -> v1-concordant), UCell DAM score shift. INTERACTION (honest, OUTCOME-OPEN): 0 large-effect DE
  genes BUT 123 stageR-confirmed (real |logFC|>0.5, median ~1.1; sub-threshold standalone per-contrast FDR, min
  adj.P 0.17) + MDE@80%=0.92 log2FC -> under-powered NOT absent (BACKED by the power
  statement, never asserted from "0 genes"; absence of evidence != evidence of absence). Static compositional
  synergy on DAM (propeller sig, sccomp borderline) -> the tau x amyloid synergy handed to P2 (trajectory) as a
  position/extent effect (progression-rate reading only under the age-matched-snapshot assumption). Thrupp 2020 caveat carried throughout (snRNA under-detects ~18% DAM genes -> SCORE
  not threshold; DE on RAW counts).
- REPORT: `_microglia.qmd` reads a COMPACT `microglia_report` target (not the 612MB Seurat) -> the gate's force-
  render stays cheap. Prose INLINE-COMPUTED from targets (never hardcoded). Two reviews hardened S3
  (concordance false-green, sccomp diagnostics + a fresh-build gate hole) + S4 (null/power over-claim, crossing/
  marker guards).

Verification (honest): every step smoke-tested vs LIVE data then full `scripts/check.sh` green end-to-end (S3 +
S4 incl. a FORCED fresh rebuild of the heavy target); S5 report renders 0-warning with the microglia chapter.
Re-baselined on R 4.6 (NOT v1's locked 18/12/55 margins). The interaction's "under-powered not absent" is the
defensible call -- the static null is reported WITH effect-size/MDE, the synergy deferred to P2 as a position/extent
effect (rate only under the age-matched-snapshot assumption).

Deferred -> P2+: activation pseudotime + interaction-as-position/extent (P2); TF/kinase mechanism, Gsk3b/Myc,
NF-kB attenuation (P3); GeoMx spatial + proteome + phospho DE, CCC, integrated divergence (P4); lean synthesis
report (P5).

## P2 Interaction trajectory -- closed 2026-07-02 (-> `.agent/completed/p2_trajectory_plan_2026-07-02.md`)

Tested the P1-deferred claim that tau x amyloid synergy is not just more DAM cells but further progression
along a homeostatic->DAM activation ordering. Built `microglia_trajectory`, `trajectory_progression`,
`trajectory_glmm_sensitivity`, `trajectory_report`, and `_trajectory.qmd`. Live facts -> memory.md
(P2-S1..S4 sections); wiring -> map.md. Below = decisions + rejected alternatives.

- TRAJECTORY: slingshot on the cached Harmony embedding, forced Homeostatic->DAM lineage, IFN/proliferative
  off-lineage but retained for conditioning audit. UCell DAM-minus-homeostatic score-axis = concordance guard,
  not independent robustness. Pseudotime is activation position/extent in a cross-sectional snapshot, NOT
  developmental time, potency, rate, or acceleration. REJECTED: CytoTRACE2/potency framing, Python PAGA/DPT,
  CellRank/RNA-velocity claims, retuning dims/lineage after contrast inspection.
- PRIMARY INFERENCE: collapse on-lineage cells to 16 genotype_batch units -> existing factorial design
  (9 residual df), ordinary weighted limma t (no eBayes), 5 canonical contrasts. Three-channel Kitagawa/Oaxaca
  decomposition splits mean pseudotime into composition / progression / cross on raw additive pt scale; primary
  family = {progression_cf, within_homeostatic}. Shared weights preserve exact reconstruction; `mean_pt` is
  composition-confounded/exploratory for the distinct P2 question, so a positive mean-position interaction alone
  would overclaim progression. Freedman-Lane replicate permutation = sensitivity, not nominal-exact. REJECTED:
  cell-level progression tests as primary (pseudoreplication), condiments/tradeSeq as inferential
  factorial tools, 2-channel decomposition, logit/asin decomposition.
- SUPPORTIVE ARM: glmmTMB beta GLMM on per-cell pt01 with batch fixed + `(1|unit)`, graceful degrade to rank-normal
  LMM or recorded failed result. It corroborates the position shift only; because it models per-cell position, it is
  composition-confounded and never load-bearing for composition-vs-progression. Live correction vs plan: P3M served
  a CRAN binary for glmmTMB/TMB on this stack, so no source-build/ABI handwork was needed. REJECTED: Stan/Python/GitHub
  dependencies, leaning on GLMM asymptotics with only 16 units.
- HEADLINE (R4.6 re-baseline, qualitative durable): amyloid strongly advances the activation ordering. The tau x
  amyloid interaction raises mean position, but the exact decomposition attributes that position shift mainly to
  DAM-cell composition, not progression beyond composition; the pre-registered progression endpoints are not
  statistically supported. This confirms P1's DAM-fraction interaction and revises the v1 Arc-M interpretation
  from "synergistic acceleration/progression" to "more DAM cells". Absence language is scoped: negative/NS progression
  estimate = no supported progression signal, NOT proof of absence.
- REPORT/API: `_trajectory.qmd` reads one compact `trajectory_report` target (~sub-MB), never the 612MB Seurat.
  Prose and headline numbers are inline-computed. `trajectory_report_data` guards every nested field the qmd reads
  and intentionally fails on non-finite inference rows rather than printing NaN. Dropped redundant
  `trd$decomposition` output after review; loadings live in provenance and per-channel coefs live in the interaction
  table. S4b review specifically hardened absence/rate wording, failed-GLMM provenance, reconstruction/loadings
  guards, and stale hardcodes.

Verification (honest): S1-S4b each smoke-tested on live cached data, fresh leaf builds were forced where the gate
would otherwise stay cached (`trajectory_progression`, `trajectory_glmm_sensitivity`, `trajectory_report`), and
full `scripts/check.sh` was green through S4b + review. Close-out adversarial review found no remaining shipped
code/prose blocker beyond updating the forward spine. Final close-out gate re-run green on 2026-07-02.

Deferred -> P3+: mechanism survey (TF/kinase: Gsk3b/Myc, NF-kB attenuation) is the default next phase; P4 remains
GeoMx spatial + bulk proteome/phospho cross-modality; P5 lean synthesis.

## P3 Mechanism -- closed 2026-07-02 (-> `.agent/completed/p3_mechanism_plan_2026-07-02.md`)

Tested the v1 mechanism candidates against the rebuilt P1/P2 result: amyloid drives microglial DAM expansion, and
mutant tau modulates that amyloid response mainly through composition rather than supported further trajectory
progression. Built direct-mouse CollecTRI/OmniPath priors, RNA TF/pathway targets, a targeted NF-kB gate, minimal
24M bulk-phosphosite DE for kinase activity, `mechanism_report`, and `_mechanism.qmd`. Live facts -> memory.md
(P3 sections); wiring -> map.md. Below = decisions + rejected alternatives.

- PRIORS/API: `decoupleR` ULM over direct-mouse CollecTRI and OmniPath KSN, with project-local REST cache and
  hash/count drift gates. Chosen because current `OmnipathR` postprocessors fail on the live schema (`ncbi_tax_id`
  absent), while the official REST TSVs remain usable and cacheable. REJECTED: v1 off-lock `nichenetr` mapping,
  silent prior refresh, and broad topology/SCENIC/CARNIVAL rebuilds.
- RNA MECHANISM: existing replicate-corrected pseudobulk DE is the inference unit; no cell-level TF/pathway tests.
  CollecTRI TF activity + bounded native-mouse MSigDB GO fgsea + compact project gene sets. NF-kB attenuation is
  a targeted gate: only whole-microglia `interaction` can support attenuation, and support requires concordant
  negative primary rows after small-family correction. `tau_in_nlgf` and substates are supportive-only.
- KINASE: built a narrow 24M bulk hippocampus phosphosite DE leaf solely to feed kinase-substrate activity.
  Exactly 16 sample-key runs, no batch term, positive log2 + median normalisation + prevalence filter, duplicate
  biological sites collapsed before decoupleR. Run-index sensitivity is carried because the 24M acquisition order is
  genotype-blocked. REJECTED: treating this as the full P4 phosphoproteomics chapter or as microglia-sorted evidence.
- HEADLINE (rebuilt, qualitative durable): Myc is the strongest whole-microglia tau-amyloid interaction TF signal
  and is negative in DAM as well, supporting the Myc thread as an RNA regulatory hypothesis. NF-kB attenuation is
  discordant (target-GSEA negative, TF-family positive), so the attenuation claim is not supported. Gsk3b is covered
  by the KSN and phosphosite data but is not significant for interaction or tau-in-NLGF; tau-in-NLGF kinase hits
  weaken under run-index sensitivity. P3 therefore supports a DAM/Myc mechanism hypothesis and downgrades v1's
  Gsk3b and NF-kB attenuation threads to unresolved/unsupported in the rebuilt data.
- REPORT/API: `_mechanism.qmd` reads one compact `mechanism_report` target and never heavy Seurat. Prose is
  inline-computed from target rows; close-out review fixed the remaining low robustness issue by making the NF-kB
  explanatory sentence branch on the actual gate status rather than hardcoding the current discordant pattern.

Verification (honest): S1-S4 each smoke-tested on live cached data with fresh leaf builds where needed
(`mechanism_collectri`, `mechanism_gene_sets`, `mechanism_pathway`, `phospho_de_24m`, `kinase_activity`,
`mechanism_report`), and full `scripts/check.sh` was green through S4. Close-out review found no remaining
shipped code/prose blocker after the NF-kB wording fix; final close-out gate re-run green on 2026-07-02.

Deferred -> P4/P5: GeoMx spatial DE, total proteome + broader phospho interpretation, CCC for the synaptic/clearance
axis if it earns it, an integrated divergence view, and the lean synthesis report. These later closed in P4/P5.

## P4 Cross-modality -- closed 2026-07-02 (-> `.agent/completed/p4_cross_modality_plan_2026-07-02.md`)

Built the non-snRNAseq evidence layer around the P1-P3 spine: replicate-aware GeoMx spatial DE, 24M total-proteome
DE, raw and parent-protein-corrected phosphosite DE, targeted clearance-axis CCC-lite, harmonised symbol/pathway
integration, `crossmodality_report`, and `_crossmodality.qmd`. Live facts -> memory.md (P4 sections); wiring ->
map.md. Below = decisions + rejected alternatives.

- GEOMX: primary = edgeR TMM + limma-voom with slide fixed effects and `duplicateCorrelation(block=genotype:bio_rep)`;
  unblocked AOI and collapsed-bio-unit fits are sensitivities. GeoMx object default assay is SCT, so counts are read
  explicitly from RNA/counts; live counts are count-like but not fully integer, so residues are recorded and values are
  not rounded. At P4 close, SpatialDecon was gated at preflight because Q3/background were usable but nuclei sentinels
  disabled absolute rescaling and no compact reference profile existed; the later follow-up built that profile and
  found the attempted fit blocked by 4 unresolved AOIs. REJECTED: treating repeated AOIs as independent animals,
  silent decon skip once preconditions become earned.
- BULK OMICS: 24M sample-key matched 16/16 runs. Proteome sums raw positive intensities to `PG.ProteinGroups` before
  log2/median-normalisation/prevalence filtering/limma-trend. Corrected phospho subtracts matched filtered parent
  protein log2 intensity, then re-filters/refits; raw P3 phosphosite DE is reused. Additive run-index sensitivity is
  load-bearing because genotype-blocked run order weakens many primary bulk hits. REJECTED: imputation, treating bulk
  hippocampus as microglia-sorted, and a proteomics-methods benchmark (`QFeatures`/`msqrob2PTM`) for this phase.
- CLEARANCE / CCC-LITE: no CellChat/LIANA/MultiNicheNet model is run. The measured-axis table covers Apoe/Trem2,
  App/Cd74, Pros1/Mertk, complement and synaptic anchors across microglia RNA, GeoMx and bulk layers. A pair earns
  support only with coherent supported evidence in >=2 modalities, or one non-microglia modality plus a strong
  whole-microglia anchor. Live qualitative call: Apoe-Trem2 earns focused measured support in amyloid-on-P301S;
  App-Cd74 and Pros1-Mertk remain measured but unearned. REJECTED: full off-lock CCC stack and v1-style LR ledger.
- INTEGRATION: `crossmodality_table` harmonises symbols across snRNAseq whole/substate RNA, GeoMx, proteome,
  raw/corrected phospho, TF activity and kinase activity. Duplicate features collapse by best FDR then absolute
  statistic with provenance retained. Count-honesty distinction is load-bearing: `modality_class` drives broad
  modality counts; `modality_group` is layer-level evidence. Pathway rows are descriptive ranked-statistic summaries,
  not formal meta-analysis or contest margins.
- HEADLINE: P4 strengthens the amyloid-response spine across tissue/spatial and bulk layers, especially antigen
  presentation / phagocytic-clearance / synaptic axes. It narrows the claims: the tau-amyloid interaction is not
  broadly spatial or bulk-omics significant, SpatialDecon abundance is blocked after attempted fitting (superseding
  P4's preflight-only negative), full CCC is not called, and bulk run-order sensitivity downgrades standalone
  proteome/phosphosite claims.
- REPORT/API: `_crossmodality.qmd` reads one compact `crossmodality_report` target (~23KB live), not GeoMx/proteome/
  phospho or the ~10MB evidence table. Prose is target-derived; close-out review fixed the remaining stale-claim risk
  by making the earned clearance-pair text branch on the current `pair_support` rows instead of hardcoding Apoe-Trem2.

Verification (honest): S1-S5 each smoke-tested live with fresh leaf builds where relevant (`geomx_de`,
`proteome_de_24m`, `phospho_corrected_24m`, `bulk_omics_summary`, `clearance_axis`, `crossmodality_table`,
`crossmodality_pathway`, `crossmodality_divergence`, `crossmodality_report`), and full `scripts/check.sh` was green
through S5. Close-out review found one low report-prose robustness issue; it was accepted/fixed before the final
close-out gate. Final close-out gate re-run green on 2026-07-02.

Closed by P5: one lean synthesis report with a compact evidence table. P4 is used as corroboration for DAM
activation, synaptic suppression, and measured Apoe-Trem2 clearance support; progression-beyond-composition, NF-kB
attenuation, Gsk3b, full CCC, and spatial abundance remain unresolved/unsupported.

## P5 Synthesis -- closed 2026-07-02 (-> `.agent/completed/p5_synthesis_plan_2026-07-02.md`)

Closed the fresh rebuild with a read-only synthesis layer and final lean report pass. Built `synthesis_report`,
`_synthesis.qmd`, report-front Overview wording, and source-prose cleanup. Live facts -> memory.md (P5 sections);
wiring -> map.md. Below = decisions + rejected alternatives.

- SYNTHESIS SHAPE: one compact target over existing compact P1-P4 report bundles only (`microglia_report`,
  `trajectory_report`, `mechanism_report`, `crossmodality_report`). `crossmodality_divergence` was not needed.
  The target returns headline strings, 10 descriptive evidence rows, status counts, unsupported/open rows, tiny source
  highlights, and provenance. Guarded anchors fail loud if the closed story drifts. REJECTED: new inference, MOFA/meta-
  analysis, human/SCENIC/topology/full CCC/decon side arcs, and v1 convergence/ledger/contest machinery.
- STATUS CONTRACT: rows use descriptive labels (`core_supported`, `corroborated`, `focused_support`,
  `not_supported`, `not_earned`, `open_caveat`), not scores. Negative rows are first-class: no supported
  progression-beyond-composition signal, NF-kB attenuation not supported, Gsk3b not recovered, SpatialDecon/full CCC
  not earned (now with SpatialDecon abundance blocked after attempted fit), and bulk run-index sensitivity. Regression
  tests reject ledger-like columns and missing anchors.
- REPORT/API: `_synthesis.qmd` is included immediately after Overview and tar_loads only `synthesis_report`; P1-P4
  chapters remain the audit trail. Final lean pass removed stale forward pointers, tightened progression/rate wording,
  and moved raw caveat cleanup into `synthesis_report_data` so the compact table itself does not carry stale phase-step
  strings.
- HEADLINE: amyloid drives a microglial homeostatic-to-DAM activation programme; mutant tau modulates that amyloid
  response mainly through extra DAM-cell composition rather than supported further activation-axis progression; Myc is
  the focused RNA mechanism signal; NF-kB attenuation and Gsk3b are not recovered; cross-modality corroborates the
  amyloid-response and synaptic-clearance axes while SpatialDecon abundance is blocked after attempted fit and full CCC
  remains absent.

Verification (honest): S1-S3 each smoke-tested or rebuilt live where relevant (`synthesis_report`, forced report);
final `scripts/check.sh` was green after the lean-pass fixes (tests warn=2, forced 88-chunk report render,
tar_meta clean across 42 current targets/branches, render-log clean). Close-out review found no remaining plan,
target, or report-source blocker.

## Figure expansion -- closed 2026-07-02 (-> `.agent/completed/figure_expansion_plan_2026-07-02.md`)

Post-report visual-density pass over the closed P1-P5 story. The user chose inline chapter expansion over a
standalone atlas. Built compact per-chapter figure targets, wired 26 additional figures into the existing chapters,
then closed with UX/cross-reference QA. Live facts -> memory.md (Figure expansion sections); wiring -> map.md.
Below = decisions + rejected alternatives.

- SHAPE: inline figures in `_synthesis.qmd`, `_microglia.qmd`, `_trajectory.qmd`, `_mechanism.qmd`, and
  `_crossmodality.qmd`, backed by compact `*_figures` targets. REJECTED: standalone figure atlas, prebuilt
  publication gallery, and schematic-first compression.
- DATA CONTRACT: figure builders return qmd-ready compact slots with finite geom guards and pre-binned heavy shapes
  for volcano/scatter plots. Qmds still load compact report/figure targets only; raw modality tables, the 612MB
  Seurat object, and the 10MB harmonised evidence table stay out of forced report renders.
- CLAIM CONTRACT: figures expose existing evidence, including unsupported/unearned rows, without new inference.
  Captions preserve the closed calls: amyloid-to-DAM activation, DAM composition interaction, no supported
  progression beyond composition, Myc supported, NF-kB attenuation and Gsk3b not supported, SpatialDecon/full CCC
  not earned, and bulk run-index sensitivity load-bearing.
- UX/CROSSREF: enabled Quarto `lightbox: auto` in the embedded offline HTML and normalised every captioned figure
  chunk to a hyphenated `fig-*` label with no underscores.

Verification (honest): S1 compact targets built warning-clean and <5MB each; S2-S4 forced renders/gates were green
as figures landed chapter by chapter. S5 full `scripts/check.sh` was green after the lightbox + label normalisation
(tests warn=2, forced 140-chunk report render, tar_meta clean across 46 current targets/branches, render-log clean).
Rendered HTML QA found 42 figure blocks, 42 captions, 42 source `fig-*` labels, expected sections present, embedded
lightbox markers, 0 external resource refs, and 0 warning/error markers.

## Spatial decon follow-up -- closed 2026-07-02 (-> `.agent/completed/spatial_decon_followup_plan_2026-07-02.md`)

Closed the P4 gap where SpatialDecon abundance was not earned because no compact reference existed. Built and ran a
GeoMx-native SpatialDecon follow-up, then rewired report/synthesis state to the actual blocked fit. Live facts ->
memory.md (Spatial decon follow-up sections); wiring -> map.md. Below = decisions + rejected alternatives.

- ROUTE: broad-first SpatialDecon with a separately gated microglia-substate attempt. Added `SpatialDecon`
  project-locally; built `geomx_reference_profile` from the full snRNAseq RDS plus retained P1 microglia labels.
  Broad and substate profiles both earned gates; Proliferative was recorded absent rather than fabricated. REJECTED:
  nuclei-rescaled absolute cell counts while 42/91 nuclei sentinels remain, plaque-niche localisation from whole-tissue
  ROIs, full CCC, v1 ledger revival, and forcing substate abundance if unstable.
- FIT OUTCOME: SpatialDecon ran warning-clean on GeoMx RNA/data with Q3-scaled negative-probe background. Broad and
  substate arms both blocked on the same 4/91 unresolved AOIs with near-zero total beta; 87 AOIs resolved and residual
  QC remains descriptive fit QC. Because broad abundance did not earn a clean fit, no log-beta abundance DE or
  microglia-substate tissue-abundance claim is made.
- INTEGRATION/API: `geomx_abundance_de` returns blocked canonical 5-contrast empty top tables plus residual audit;
  `spatial_decon_report` is the compact handoff (no beta matrices). `clearance_axis`, `crossmodality_report`, and
  `synthesis_report` now derive SpatialDecon status from this attempted fit, not the historical P4 preflight. A future
  earned SpatialDecon/full-CCC state intentionally stops synthesis until the claim is revised.
- HEADLINE: SpatialDecon is no longer "not attempted" or "missing profile"; it was attempted and abundance is blocked
  by unresolved AOIs. The final synthesis keeps spatial abundance and full CCC outside the supported claim set while
  retaining GeoMx DE and measured Apoe-Trem2 CCC-lite support.

Verification (honest): S1-S5 each smoke-tested live with fresh leaf builds where relevant
(`geomx_reference_profile`, `geomx_decon`, `geomx_abundance_de`, `spatial_decon_report`, `clearance_axis`,
`crossmodality_report`, `synthesis_report`) and full `scripts/check.sh` was green through S5. Close-out review found
no accepted shipped code/prose blocker; final close-out gate re-run green on 2026-07-02.

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

# Completed plan: human cross-species validation of the microglial tau x amyloid interaction

## Outcome (plan closed 2026-06-04)

**Verdict: the mouse microglial tau × amyloid interaction structure is
conserved in human AD cortex as a DIRECTIONAL corroboration under
collinearity-aware modelling — "consistent with", never proof.** Tested in
SEA-AD (Gabitto 2024) MTG + DLPFC microglia (no DUA; ROSMAP confirmation
deferred for want of one). Three load-bearing results: (1) **conservation gate
PASS** — all four mouse substates (homeostatic/DAM/IFN/proliferative) resolve
as human populations, with disease-emergent SEA-AD `-SEAAD` supertype occupancy
rising monotonically along the activation axis (0.150 → 0.239 → 0.359 → 0.490);
(2) **interaction directions reproduce** — the pre-registered mechanisms match
their mouse-predicted `amyloid:tau` sign in 5/5 strata, region-concordant and
robust to single-cell aggregation + an orthogonal AUCell scorer:
`NFKB_union_targets` and `MG_M3_module` negative (tau attenuates amyloid),
`Gsk3b_targets` positive (tau-enhanced); (3) **positive control anchors the
pipeline** — human-native Gerrits AD1 shows the expected +amyloid main effect
(5/5 strata) and +mediation ACME (+0.056, nominal p = 0.012). No pre-registered
human interaction survives FDR — the expected consequence of amyloid–tau
collinearity (Spearman ≈ 0.65, VIF interaction ≈ 1.5–1.8) and observational
power, NOT counter-evidence (reading rules 1–2).

**Deliverables:** new chapter `rmd/18_human_validation.Rmd` → analysis.html
**§19** (5 figures + 7 tables); two curated panels in `summary.Rmd` **§7
"Cross-species validation in human AD"** (conservation gate + mechanism ×
stratum heatmap); the H1–H5 build scripts + light caches (see map.md). The
locked mouse analysis (the §17 biological-model ledger and its contest margins)
was **deliberately left untouched** — see the H7 completion note for the
four-point rationale. Compact decision-summary in `completed/DIGEST.md`.

**Deferred follow-on (not done, no DUA):** a powered cohort that genuinely
decouples amyloid and tau (ROSMAP / Sun 2023, N≈443, MG0–MG12) would convert
"directional corroboration" toward a significance-level test; the within-cohort
MTG-vs-DLPFC replication is the available independent corroboration. Three other
"sophisticated expansion" lanes from the originating menu remain candidate
plans: CARNIVAL/COSMOS causal signalling, MOFA+ latent integration, SCENIC
regulons.

## STATUS (read this first, then read only the next TODO step's body)

- [x] H1 lock dataset, set up env, acquire human microglia — DONE 2026-06-01
- [x] H2 assemble mouse signatures + map to human orthologs — DONE 2026-06-01
- [x] H3 human microglial substate conservation — DONE 2026-06-04 (GATE PASS:
      DAM & IFN resolve; per-state H4 claims cleared)
- [x] H4 amyloid×tau interaction models (core; collinearity-aware) — DONE
      2026-06-04 (sanity #5 passes: Gerrits AD1 +amyloid; mouse interaction-
      attenuation reproduces DIRECTIONALLY — MG-M3 & NF-κB negative amyloid:tau
      in 5/5 strata, Gsk3b positive 5/5, region/discordant-robust; no pre-reg
      interaction survives FDR under r≈0.65 collinearity → "consistent with")
- [x] H5 single-cell robustness + mediation — DONE 2026-06-04 (H4 interaction
      DIRECTIONS robust: single-cell donor-RE refit concordant 129/130, AUCell
      orthogonal-scorer concordant 82/90; NF-κB & MG-M3 negative + concordant
      across BOTH arms in 5/5 strata. Mediation (Green-2024 amyloid→sig→tau,
      donor-cluster bootstrap): disease sigs positive ACME, human-native
      Gerrits AD1 ACME +0.056 p=0.012 nominal, none survives FDR under
      collinearity. ROSMAP confirmation deferred (no DUA); MTG-vs-DLPFC is the
      within-cohort replication.)
- [x] H6 render the human-validation chapter — DONE 2026-06-04 (new
      rmd/18_human_validation.Rmd → §19, 5 figures + 7 tables, 0 errors /
      0 warnings; conservation GATE PASS + interaction directions
      concordant 5/5 strata rendered as the chapter narrative)
- [x] H7 surface in summary.Rmd + synthesis, close plan — DONE 2026-06-04
      (two curated panels added as summary.Rmd §7 "Cross-species validation in
      human AD" via the build_summary_rmd.R generator; biological-model ledger
      DELIBERATELY left untouched — rationale logged; plan closed)

Full step bodies and completion notes are below; the Execution model section
defines the per-step protocol.

- **Created:** 2026-06-01 (user chose this from a 4-way "sophisticated expansion"
  menu after the `summary.Rmd` pass at commit `5cf9c6a`; the other three lanes —
  CARNIVAL/COSMOS causal signalling, MOFA+ latent integration, SCENIC regulons —
  remain candidate follow-on plans).
- **Predecessor context:** the internal synthesis lives in `analysis.html`
  (through §18) and the curated `summary.html`; this plan ADDS an outward-facing
  translational layer and does not modify the locked mouse analysis.
- **Goal:** test whether the mouse 2x2 tau x amyloid microglial findings --
  especially the **interaction-localised** mechanisms (DAM programme,
  MG-M3 co-expression module, IFN-substate amyloid-response asymmetry,
  tau-driven NF-kB attenuation, Gsk3b-target kinase activity) -- reproduce
  in **human AD**, where amyloid and tau co-occur but vary continuously and
  partly discordantly across donors. The human analogue of the mouse
  interaction contrast is the cross-donor **amyloid x tau interaction term**
  on microglial signature scores / pseudobulk expression, fit with explicit
  **collinearity mitigation**.

## Why this plan exists

The project's central claim is a tau-dependent amyloid remodelling of
microglia, null at the whole-microglia gene level but real once resolved to
substates / modules / kinases / signalling. Every result to date is from
four **mouse** genotypes. Translational weight requires asking whether the
same structure exists in **human** AD microglia. Two things make this a
genuine test rather than a box-tick:

1. **Substate / programme conservation.** Do the mouse microglial substates
   (homeostatic / DAM / IFN / proliferative) and the interaction-carrying
   programmes (MG-M3, NF-kB targets, Gsk3b targets) resolve in human
   microglia at all? Published human atlases (Sun 2023 MG0-MG12; SEA-AD
   supertypes) already define DAM-like (MG4/MG5), interferon (MG11) and
   cycling (MG12) states, giving concrete targets to map onto.

2. **The interaction itself.** Humans have no clean 2x2: amyloid and tau
   co-progress, so a cross-donor `amyloid x tau` product is **collinear**
   with its main effects and is easily under-identified. A null human
   interaction is therefore ambiguous (biology absent vs. design
   under-powered/collinear). The plan must report identifiability
   (VIF, off-diagonal donor coverage) and direction-of-effect concordance
   alongside every interaction estimate, and lean on the dataset that best
   decouples the two axes (SEA-AD's designed continuum with separately
   image-quantified Abeta and AT8).

This reframing is load-bearing: the deliverable is **"is the mouse
mechanism structure conserved, and does its interaction signature show a
tau-dependent amyloid response in human, under collinearity-aware
modelling"**, NOT a literal replication of `(NLGF_P301S - P301S) -
(NLGF_MAPTKI - MAPTKI)`.

## Assets already on disk (no recomputation needed)

Mouse signatures are all exportable from existing result tables / caches
(H2 assembles them; cross-species mapping via `nichenetr`, the project
standard used in the custom-signature and kinase builds):

| Mouse signature | Source on disk |
|---|---|
| DAM / amyloid-microglia programme | `storage/cache/custom_microglia_ad.rds` (AD1/AD2/LDAM/WAM); `storage/results/ccc_dam_sampling_balance_*` for DAM markers |
| Microglia substate markers (4) | `storage/cache/custom_microglia_states.rds`; `storage/cache/summary_mg_substate_umap.rds`; per-state NEBULA caches |
| MG-M3 co-expression module + hubs | `storage/results/hdwgcna_modules.tsv`, `hdwgcna_module_hubs.tsv` |
| IFN-substate asymmetry genes | `storage/results/ifn_nlgf_asymmetry.tsv` (P301S-specific / shared) |
| NF-kB target union set | `storage/cache/per_state_nfkb_target_gsea.rds`; `storage/results/per_state_nfkb_target_gsea.tsv` |
| Gsk3b-target set | `storage/cache/kinase_activity_decoupler.rds` + OmniPath KSN cache (substrates of Gsk3b) |
| Interaction-direction signature | `storage/cache/de_snrnaseq_nebula_per_state*.rds`, `storage/results/nebula_per_state_*` (interaction up/down genes per state) |
| Human-native positive controls | Gerrits 2022 AD1/AD2 are human-origin (embedded in `build_custom_microglia_ad.R`); score them in human as a pipeline sanity anchor |

## Locked decisions (carried + new)

| Decision | Choice |
|---|---|
| Cross-species gene mapping | **`nichenetr`** mouse<->human ortholog conversion (project standard; same tool used to build AD1/AD2 and the human KSN). Report per-signature ortholog-mapping coverage so dropouts are visible. |
| Heavy-data tooling | **Python (scanpy / anndata / cellxgene-census)** for h5ad acquisition, microglia subsetting (backed mode), per-cell state scoring and per-donor (x state) pseudobulk export, run in a **project-local `uv` venv** (gitignored). **R** for signature assembly, statistical interaction models and plotting (the project's home turf), reading the Python-exported pseudobulk + neuropath metadata. Heavy work in `scripts/` pre-builds; the Rmd consumes light caches (project convention). |
| Single-cell scoring | per-cell `scanpy.tl.score_genes` and/or `AUCell` for robustness; the **primary** statistic is the **per-donor (x state) pseudobulk signature mean** so donor is the unit of replication (mirrors the mouse genotype_batch pseudobulk logic). |
| Human "interaction" model | per-donor pseudobulk `score ~ amyloid * tau + covariates` (age, sex, region, post-mortem interval where available); the **`amyloid:tau` coefficient** is the human analogue of the mouse interaction. Collinearity mitigation is mandatory (VIF report; residualised re-fit; discordant-staging-donor sensitivity subset). A **mediation** framing (amyloid -> microglia signature -> tau, per Green 2024) is a secondary complementary causal model, not a replacement. |
| Colour / house style | reuse existing project palettes (the user elected to **keep the current colour system**, incl. DAM = `#B0344D`); British English prose; new files chowned rstudio:rstudio. |
| Synthesis venue | a **new top-level chapter** via a new `rmd/NN_human_validation.Rmd` appended to `analysis.Rmd` before `child-session` (rmd/99); 1-2 curated panels later mirrored into `summary.Rmd` (H7). Exact section number locked at H6. |

## Open questions (defer to the session that needs them)

| Question | Default proposal | Decided in |
|---|---|---|
| Human dataset & access path | **SEA-AD MTG (Gabitto 2024)** via CELLxGENE Discover h5ad (continuous Abeta + AT8 + Braak/Thal/CERAD/CPS in `obs`; turnkey, no DUA; N=84 discovery cohort). Alternatives: SEA-AD MTG+DLPFC (within-cohort regional replication, no DUA); ROSMAP/Sun 2023 (powered N=443, MG0-MG12, but Synapse/RADC **DUA-gated** pathology); SEA-AD + ROSMAP (best rigour, ROSMAP gated). | **DECIDED 2026-06-01:** SEA-AD **MTG + DLPFC** (both regions, no DUA). User holds no RADC/Synapse DUA, so ROSMAP confirmation is a deferred follow-on. |
| Disk / RAM for the source h5ad | SEA-AD MTG is ~1.2M nuclei (multi-GB). Default: load in scanpy **backed mode**, subset microglia/PVM on the fly, write a compact microglia-only cache; delete the full download afterwards. Confirm container disk headroom at H1. | H1. |
| State assignment in human | Default: use the **dataset-native** microglial state/supertype labels AND independently **score the 4 mouse substate signatures per cell -> argmax** (mirrors the mouse AddModuleScore-argmax), then cross-tabulate the two. Use the agreement as the conservation readout. | H3. |
| Per-state vs whole-microglia interaction | Default: run the interaction model BOTH whole-microglia and per-state, but gate per-state claims on H3 conservation passing (anti-anchoring #3). | H4. |

## Execution model

- A fresh session (latest/largest model, max reasoning) reads the STATUS block
  at the top, then reads ONLY the next `TODO` step's body (plus any step it
  explicitly depends on) and executes it. Per step:
  1. Update that step's status from `TODO` to `DONE <YYYY-MM-DD>` with a
     multi-paragraph completion note (what was built, file paths, biology
     highlights, explicit plan-spec deviations).
  2. Before writing code, check for a fitting Skill / subagent (step 0 of
     the session loop): `cellxgene-census` (data), `scanpy`/`anndata`
     (h5ad), `pathway-enrichment`/`statistical-analysis` (scoring + model
     choice), `scientific-visualization` (figures), `Explore` (code search).
  3. Smoke-test new cache-reading / model code with a cheap shape/dtype
     check before any full knit (~3 min).
  4. chown rstudio:rstudio all new files; chown knit-created root-owned
     outputs (analysis.html, TSVs, caches).
  5. Re-knit where the step touches the Rmd and verify
     `grep -c 'class="error"'` and `class="warning"` are both 0.
  6. Commit locally (imperative subject < 70 chars + HEREDOC body +
     Co-Authored-By trailer per CLAUDE.md).
- Chain non-gated steps in one session only when each is a trivial
  pattern-extension; H1-H5 each carry design weight, so default to a fresh
  session per step. End cleanly at a gate, a fork, or low context.
- Decision gates require user confirmation via AskUserQuestion (present the
  default PLUS >=1 reasoned alternative).
- Keep the project-local Python `uv` venv path in `.gitignore`.

---

### Session H1: lock dataset, set up env, acquire human microglia [DONE 2026-06-01]

**Decision gate at start.** Choose the human dataset & access path
(default SEA-AD MTG; alternatives SEA-AD MTG+DLPFC, ROSMAP/Sun 2023 [DUA],
SEA-AD + ROSMAP [DUA]) -- surface via AskUserQuestion. Ask in the same gate
whether the user already holds a ROSMAP/RADC Synapse DUA (decides whether
the gated options are turnkey).

Deliverables once the gate resolves:
- Project-local `uv` venv with `cellxgene-census`, `scanpy`, `anndata`
  (+ `aucell`/`decoupler` if used); add the venv dir to `.gitignore`.
- `scripts/build_human_microglia.py`: download the chosen dataset's h5ad
  (SEA-AD: CELLxGENE Discover collection `1ca90a2d-2943-483d-b678-b809bf464c30`),
  load backed, subset to microglia/PVM, basic QC, persist a compact
  microglia-only AnnData/cache + a per-donor neuropathology table
  (continuous Abeta, AT8/pTau, Braak, Thal, CERAD, CPS/pseudoprogression,
  age, sex, region) to `storage/cache/` and `storage/data/`.
- `storage/cache/human_microglia_provenance.txt` (dataset id, census/portal
  version, download URL, n donors, n microglia, QC thresholds, which
  neuropath columns came from native obs vs a joined supplement).
- Smoke-test: print donor count, microglia count, neuropath-column dtypes,
  and the Abeta-vs-tau scatter + correlation (the collinearity baseline).
- Delete the full multi-GB download once the microglia subset is cached.

**Completion note (2026-06-01).** Gate resolved to SEA-AD MTG + DLPFC
(no DUA; ROSMAP deferred). SEA-AD ships **pre-split microglia/PVM h5ads
per region** on CELLxGENE Discover, so the planned subset-from-Whole-
Taxonomy (50 GB) step was unnecessary -- downloaded the two microglia-only
h5ads directly (`storage/data/seaad/microglia_{mtg,dlpfc}.h5ad`, 0.51 +
0.72 GB) via the curation API; these ARE the compact cache (no re-subset,
nothing to delete). Reused the existing project `.venv` (scanpy 1.12.1 +
anndata 0.12.16 already present); **no new Python deps** and no
cellxgene-census install needed (portal h5ad + curl + jq sufficed).

Built `scripts/build_human_microglia.py` -> emitted
`storage/cache/human_seaad_donor_neuropath.csv` (172 donor x region rows;
84 MTG + 80 DLPFC AD-continuum donors plus Allen neurotypical Reference
donors flagged `is_reference`) and `storage/cache/human_microglia_
provenance.txt`. The two regions share an identical 35,483-gene space;
expression is normalised log1p in `X` with raw counts in `.raw` and
symbols in `var.feature_name`.

Key design facts carried to H2-H4: (1) the human amyloid axis is **Thal
phase (0-5)** and the tau axis **Braak stage (0,2-6; no Braak I)**, both
ordinal -- the continuous image-quantified Abeta/AT8 densities are NOT in
the CELLxGENE obs (they live in brain-map.org supplements; joining them is
a documented optional enhancement). CERAD (0-3) and ADNC (0-3) are
secondary axes; CPS is **MTG-only** (DLPFC lacks it). (2) Identifiability
is workable, not fatal: Spearman(Thal, Braak) = 0.654 (MTG) / 0.643
(DLPFC), with genuine off-diagonal donors (Thal 0 / Braak IV = 6 donors;
Thal 3-4 / Braak III-IV = 9) -- the discordant-donor sensitivity subset
(anti-anchoring #1) is real. (3) Microglial substates are present as
`Supertype`: Micro-PVM_1/2 (homeostatic-leaning) plus disease-associated
Micro-PVM_2_x-SEAAD / _3-SEAAD / _4-SEAAD -- H3 maps the mouse 4 substates
onto these by per-cell signature scoring rather than relying on names
(no IFN/proliferative-named supertype exists, so conservation must be
established empirically per anti-anchoring #3). Re-run: `./.venv/bin/python
scripts/build_human_microglia.py` (idempotent; reads the two h5ads).

### Session H2: assemble mouse signatures + human orthologs [DONE 2026-06-01]

- `scripts/build_human_validation_signatures.R`: read the H2-table sources
  above, build a tidy named list of mouse-symbol signature vectors
  (DAM-up, the 4 substate markers, MG-M3 module + hubs, IFN P301S-specific,
  NF-kB target union, Gsk3b targets, interaction-up / interaction-down per
  state), map each to human orthologs via `nichenetr`, and persist
  `storage/cache/human_validation_signatures.rds` +
  `storage/results/human_validation_signature_membership.tsv`
  (signature, n mouse genes, n mapped human genes, % mapped).
- Include the human-native Gerrits AD1/AD2 as positive-control signatures
  (no mapping needed; they are the sanity anchor).
- Smoke-test: per-signature human-ortholog coverage; flag any signature
  with < ~60% mapping for cautious interpretation.

**Completion note (2026-06-01).** Built
`scripts/build_human_validation_signatures.R` ->
`storage/cache/human_validation_signatures.rds` (list of four:
`human` = named list of HUMAN-symbol vectors to score in SEA-AD,
`mouse` = pre-mapping mouse originals, `meta` = membership data.frame,
`nfkb_family` = the NF-kB TF set), plus
`storage/results/human_validation_signature_membership.tsv` and a
provenance file. The panel is **26 signatures: 24 mouse->human mapped +
2 human-native Gerrits controls**. Mapping uses
`nichenetr::convert_mouse_to_human_symbols` -- the REVERSE of the
project's only prior direction (human->mouse, used to build AD1/AD2 and
the KSN); no project code had mapped mouse->human before, so the same
NA/empty-drop + dedup + sort idiom was reused. `pct_mapped` = fraction of
input mouse symbols receiving a human ortholog; the script flags
`< 60%` and prints a per-signature coverage panel.

**Load-bearing design resolutions (each verified against the codebase
before coding, not assumed):**
(1) **"4 substate markers" = the canonical `AddModuleScore->argmax`
marker lists** in `R/constants.R::canonical_microglia_markers`
(Microglia/DAM/IFN/Proliferative, 3-6 genes each). The 4 microglia
substates (homeostatic/DAM/IFN/proliferative) were labelled in
`R/microglia.R` by scoring exactly these lists and taking the per-nucleus
argmax -- **there is no de-novo FindAllMarkers/presto per-state marker
table anywhere** (confirmed by full-tree search), and the 7-set
`custom_microglia_states.rds` collection is fgsea-only, never used for
labelling. So H3 must score **these four lists** per human cell to
reproduce the mouse argmax. This is the single most important carry to
H3.
(2) **Interaction-direction signatures come from the localisation table,
not an FDR threshold.** The per-state interaction contrast is FDR-NULL at
the whole-gene level (`nebula_per_state_summary.tsv`: `n_sig = 0` in every
state), so an FDR-thresholded interaction set would be empty. The
directional sets are instead taken from the curated
`nebula_per_state_localisation.tsv` (the 62 genes whose interaction
effect was examined; **only 25 localise to a single state**, the other 36
have `max_negLogP_state = NA` and are correctly excluded), split by the
sign of that state's interaction logFC. Per-state sets are consequently
small (1-6 genes; `interaction_up_DAM` n=6, `interaction_down_DAM` n=1);
pooled sets are `interaction_up_pooled` n=13 / `interaction_down_pooled`
n=12. Sizes are in the membership TSV so H4 can gate on them.
(3) **NF-kB target union** = `unique(target[source %in% family])` from
`collectri_mouse_for_nfkb_gsea.rds` with
`family = Rela/Nfkb1/Nfkb2/Rel/Relb/A0A979HLR9` (identical to
`rmd/17_nfkb_attenuation.Rmd`): 1679 mouse entries -> 1304 human (80%;
the ~20% loss is the minority of UniProt accessions in the CollecTRI
mouse `target` column, expected not erroneous).
(4) **Gsk3b targets** = substrate genes of `Gsk3b` from
`build_omnipath_ksn_mouse()` (KSN `target` ids stripped of their
`_residue` suffix): 538 mouse -> 538 human (100%; the OmniPath KSN is
natively human so the round-trip is near-lossless). Human-native
alternative (`fetch_omnipath_ksn_human`, `source == GSK3B`, 544 genes)
noted in provenance.
(5) **Gerrits AD1/AD2 human controls** are embedded verbatim from
`build_custom_microglia_ad.R` (Table S7, 50 symbols each) and protected
by a **build-time drift-guard**: `convert_human_to_mouse_symbols` of the
embedded lists must `setequal` the cached mouse `AD1`/`AD2` in
`custom_microglia_ad.rds`, else the build aborts. They are scored as-is
(no mapping) -- the pipeline sanity anchor for H3/H4.

**Coverage highlights (full panel in the TSV; no cherry-picking per
anti-anchoring #2):** DAM_up 278->192 (70%), DAM_down 90%, LDAM 100%,
WAM 87%; substate markers 75-100%; MG_M3_module 780->740 (95%),
MG_M3_hubs 92%; IFN_asym_P301S_specific 446->422 (95%); NFKB_union 80%;
Gsk3b 100%. The only `low_coverage` flag is `interaction_down_IFN`
(2->1, 50% -- a tiny-n artefact, not a mapping failure).

**Deviations / deferrals.** Data-universe intersection (which mapped
human genes are actually present in the SEA-AD microglia `var`) is
deferred to H3 where the h5ads load, keeping H2 R-only with no h5ad read
(plan division of labour). No Rmd was touched, so no re-knit was needed
this step. Re-run (idempotent):
`Rscript scripts/build_human_validation_signatures.R [--overwrite]`.

### Session H3: human microglial substate conservation [DONE 2026-06-04]

- **Context from H2 (DONE):** the "4 mouse substate signatures" are the canonical
  AddModuleScore->argmax marker lists in `R/constants.R` (no de-novo per-state
  marker table exists); score exactly these to reproduce the mouse labelling. The
  full H2 panel (24 nichenetr-mapped mouse signatures + 2 human-native Gerrits
  controls) is persisted to `storage/cache/human_validation_signatures.rds` + the
  membership TSV.
- Score the 4 mouse substate signatures per human microglial cell
  (scanpy `score_genes` / AUCell), argmax to a predicted substate, and
  cross-tabulate against the dataset-native states (SEA-AD supertypes or
  Sun MG0-MG12). Quantify conservation (e.g. adjusted Rand / normalised
  mutual information; expected mapping DAM->MG4/MG5, IFN->MG11,
  proliferative->MG12, homeostatic->MG0).
- Export per-cell predicted state + per-donor x state composition, and a
  per-donor x state pseudobulk expression matrix for H4.
- **Gate downstream per-state claims on this step:** if DAM / IFN do not
  resolve as distinct human populations, restrict H4 to whole-microglia
  (anti-anchoring #3).

**Completion note (2026-06-04).** Built
`scripts/build_human_substate_conservation.py` (Python/scanpy; ~3 min over
both regions). It bridges the H2 RDS into Python via a jsonlite dump
(`storage/cache/human_validation_signatures_human.json`; pyreadr is absent),
loads each SEA-AD microglia h5ad fully into memory (sparse; 57 GB free),
relabels `var_names` from Ensembl to `feature_name` symbols, and scores ALL
26 H2 signatures per cell with `sc.tl.score_genes` on `X` (log1p,
`use_raw=False` -- `.raw` holds raw COUNTS and would be wrong to score),
`n_bins=12, ctrl_size=50, random_state=0` to match the mouse
`Seurat::AddModuleScore(nbin=12, ctrl=50)`. `predicted_substate` = argmax of
the four substate scores, faithfully reproducing
`R/microglia.R::label_microglia_states`. Scoring all 26 (not just the 4
substates) is a deliberate efficiency extension: the h5ad load dominates
runtime, so H4/H5 inherit the full per-cell + per-donor-x-state score panel
without re-loading -- the conservation analysis itself uses only the four
substate signatures.

**GATE RESULT = PASS (DAM resolves=True, IFN resolves=True), so H4 per-state
claims are cleared.** Resolution rule (transparent, stated in provenance):
`resolves := (cell fraction >= 0.5% of microglia) AND (mean self-score > 0)`;
all four states pass. The biological conservation signal is directionally
correct and is the real readout (the gate rule is the floor): predicted
homeostatic -> DAM -> IFN -> proliferative cells show a MONOTONIC rise in the
fraction sitting in disease-emergent `-SEAAD` supertypes (0.150, 0.239,
0.359, 0.490; fold-over-global 0.67, 1.06, 1.59, 2.17), and each activation
state's single most-enriched supertype is a `-SEAAD` one (DAM ->
Micro-PVM_4-SEAAD 1.7x; IFN -> Micro-PVM_1_1-SEAAD 5.8x; proliferative ->
Micro-PVM_2_2-SEAAD 12.6x), while homeostatic is DEPLETED in disease
supertypes (0.67x). I.e. the mouse-defined activation states map onto the
human disease-associated supertypes in the expected order.

**Honest caveats (anti-anchoring #2).** (1) Label-partition agreement is LOW
-- ARI 0.054, AMI/NMI 0.075, Cramer's V 0.295 (all regions) -- because the
discrete 8-level SEA-AD `Supertype` taxonomy is dominated by one
homeostatic-leaning supertype (Micro-PVM_2, ~75-85% of microglia) whereas
the four-state argmax is a continuous-activation cut; the two capture
different structure. Conservation is therefore evidenced by the directional
`-SEAAD` enrichment gradient above, NOT by cluster-label agreement, and the
provenance states this. (2) The argmax leans DAM-ward in this aged human-AD
cortex (DAM ~57% of cells, mean self-score 0.70 -- the highest of the four),
unlike the homeostatic-majority mouse; this is consistent with a broadly
activated human microglial population and the relative nature of the argmax,
and is flagged so H4 does not over-read absolute DAM abundance. (3)
proliferative has the lowest self-score (0.18) and smallest fraction (3.9%)
-- it resolves but is the weakest state; H4 should treat per-proliferative
claims cautiously.

**Data-universe intersection (the H2-deferred check) is clean:** every one of
the 26 signatures is >=60% present in the SEA-AD 35,483-gene var, and the
four substate signatures are 100% present in both regions
(`human_substate_signature_coverage.tsv`); lowest overall is
Gerrits_AD2_human at 92%.

**Outputs.** caches (H4/H5-consumed): `human_substate_percell.csv.gz` (82,486
cells x 26 scores + labels + UMAP), `human_substate_score_means.csv` (854
rows = per donor x region x {4 states + "all"} mean scores + n_cells; the
PRIMARY H4 modelling input -- donor x region join key matches the H1
neuropath CSV 172/172), `human_substate_pseudobulk_counts.csv.gz` (2,962
signature-union genes x 854 samples, raw-count sums) +
`human_substate_pseudobulk_samples.csv`. results (human-readable):
`human_substate_composition.tsv`, `human_substate_crosstab.tsv`,
`human_substate_conservation_metrics.tsv`, `human_substate_signature_
coverage.tsv`, plus `human_substate_conservation_provenance.txt` (cache).

**Deviations / notes.** (a) CSV(.gz) throughout, not parquet: neither pyarrow
(Python) nor arrow (R) is installed and H4 is R, so CSV is the friction-free
exchange format. (b) Pseudobulk is raw-count SUMS restricted to the 2,996
signature-union genes (34 of which are absent from SEA-AD var -> 2,962),
keeping the cache light while remaining count-based for any H4 re-modelling;
raw counts pulled from `.raw` with a positional gene-order assertion against
the main `var`. (c) Cells from ALL donors (continuum + neurotypical
reference) are scored and exported; H4 filters via the `is_reference` flag /
neuropath join. (d) No Rmd touched (the chapter is H6), so no re-knit this
step. Re-run (deterministic): `./.venv/bin/python
scripts/build_human_substate_conservation.py` (`--smoke N` for an N-cell-per-
region dry run).

### Session H4: amyloid x tau interaction models (core) [DONE 2026-06-04]

- For every signature (and the human-native controls), compute per-donor
  (and per-donor x state) pseudobulk scores, then fit
  `score ~ amyloid * tau + covariates`. Report main effects + the
  `amyloid:tau` coefficient with CI, p, and BH-FDR **across the full
  signature panel** (no cherry-picking).
- Collinearity battery (mandatory): VIF per model; residualised re-fit
  (tau on amyloid and vice versa); a discordant-staging-donor sensitivity
  subset (high-amyloid/low-tau and low-amyloid/high-tau corners).
- Pre-registered directional hypotheses to test (report direction, not
  just significance): DAM-up = positive amyloid main effect (sanity);
  **NF-kB targets = negative `amyloid:tau`** (the mouse "tau attenuates
  amyloid-driven NF-kB" prediction); Gsk3b targets = tau-dependent amyloid
  response; MG-M3 = negative interaction.
- Export `storage/results/human_interaction_models.tsv` + a light plotting
  cache (`storage/cache/summary_human_validation.rds`).

**Completion note (2026-06-04).** Built
`scripts/build_human_interaction_models.R` (R-only; lme4 + lmerTest +
data.table; deliberately does NOT source `R/helpers.R`, which would attach
nichenetr/Seurat for nothing — a local `write_tsv_safe` clone is the only
helper needed). It joins the H3 per-donor×region×state score means to the H1
neuropath table (perfect 172/172 key match) and, for each of the **26
signatures × 5 strata** (whole-microglia `all` + the four substates, cleared
because the H3 gate passed), fits the donor-replicated mixed model
`z(score) ~ amyloid_c * tau_c + region + age_c + sex + (1|donor)`,
precision-weighted by `n_cells`. amyloid = Thal (0-5), tau = Braak (0/2-6),
both centred on the analysed rows; response z-scored per signature so the
`amyloid:tau` coefficient is a comparable cross-panel effect size. The 8 Allen
reference rows drop out (no staging), leaving 164 AD-continuum donor×region
rows (84 donors, 83 paired across MTG+DLPFC → the random intercept is well
identified). Per-state rows require `n_cells ≥ 10` (raw IFN/proliferative
medians are 23/14, mins of 1) so a donor's state mean is not a 1-cell
artefact. Outputs: `storage/results/human_interaction_models.tsv` (130 rows ×
44 cols), `storage/cache/summary_human_validation.rds` (172 KB light plotting
cache: `models`, per-stratum `frames` with scores+staging, `vif`, `prereg`,
`meta`), and `..._provenance.txt` (carries a computed headline so H6 needs no
re-run). All four chowned rstudio:rstudio. **No Rmd touched (chapter is H6) →
no re-knit this step**, matching H1-H3.

**Collinearity battery (anti-anchoring #1).** Identifiability is clean, not
fatal: centred-design VIFs are amyloid≈1.53 / tau≈1.90 / interaction≈1.53
(proliferative the highest, tau 2.22 / ix 1.76) — all far below the
VIF=5 concern line; Spearman(amyloid,tau)=0.649; 37-43 off-diagonal
"discordant" donors per stratum. Each interaction is reported three ways:
the joint-model estimate, a **discordant-donor** refit (off-diagonal half by
|resid(tau~amyloid)|), and **independent per-region** MTG+DLPFC OLS refits with
a `region_concordant` flag. **Reasoned deviation from the plan spec:** the
plan listed a "residualised re-fit (tau on amyloid and vice versa)", but by the
Frisch-Waugh-Lovell theorem that returns the joint model's interaction
coefficient unchanged (the joint amyloid*tau fit already gives
collinearity-adjusted PARTIAL effects), so it adds no information — I verified
this empirically in the first smoke run (ix_rc_est === ix_est to 4 s.f.) and
replaced it with the genuinely independent per-region replication (also the
reason H1 acquired both regions; it pre-feeds H5's robustness check). VIF
quantifies the inflation the residualisation would notionally have addressed.

**Result (the panel is the result, not any single row — anti-anchoring #2).**
(a) **Sanity #5 passes:** the human-native Gerrits **AD1** amyloid MAIN effect
is positive in 5/5 strata (nominally sig in DAM, p=0.046), validating the
scoring/modelling pipeline. Gerrits AD2 and the mouse DAM_up amyloid main
effects are ~null — the expected consequence of amyloid↔tau collinearity (tau
absorbs the shared severity variance: tau main effects are the larger, e.g.
AD1 tau +0.246, p=0.017) plus the DAM-saturated aged-human cortex H3 flagged.
(b) **The mouse interaction-attenuation mechanisms reproduce
DIRECTIONALLY and consistently:** `MG_M3_module` shows the mouse-predicted
**negative** amyloid:tau in **5/5 strata** (region-concordant in 4-5),
`NFKB_union_targets` (the tau-attenuates-amyloid-NF-κB headline) **negative in
5/5**, and `Gsk3b_targets` a consistent **positive** (tau-enhanced) amyloid
response in 5/5. DAM_up's interaction is also negative (attenuation theme).
22/30 pre-registered rows are sign-concordant. (c) **Honest significance:**
only ONE of 130 models survives within-stratum BH-FDR — the
`substate_proliferative` score interaction in `all` (+0.167, FDR=0.022), and
it is DLPFC-driven (MTG +0.083 vs DLPFC +0.265) in H3's weakest state, so it is
reported but not leaned on. **No pre-registered interaction is FDR-significant**
— under r≈0.65 collinearity and observational human power this is expected, and
anti-anchoring #1/#7 require reading it as DIRECTIONAL, region/discordant-robust
**corroboration ("consistent with")**, never proof or refutation. Residual
normality is violated in 79/130 models (bounded score means) → p-values are
approximate, which is exactly why direction (normality-free) is the primary
readout.

**Carries to H5/H6.** H5: the per-region OLS estimates already in the TSV are a
first robustness pass — extend with a single-cell mixed/AUCell re-test and the
Green-2024 mediation framing (amyloid→signature→tau). H6: the rds `frames` give
the Thal-vs-Braak discordant-donor scatter (panel a) and per-signature score
trends; `models` gives the forest plot (panel c) and the mouse-vs-human
direction table (panel e); the provenance headline is the prose-ready summary.
Re-run (idempotent): `Rscript scripts/build_human_interaction_models.R
[--overwrite|--smoke]`.

### Session H5: single-cell robustness + mediation (+ optional confirmation) [DONE 2026-06-04]

- Single-cell mixed-model / AUCell re-test as robustness against the
  pseudobulk result; report concordance.
- Mediation model (amyloid -> microglia signature -> tau, per Green 2024)
  as a complementary causal framing; contrast the mediation vs interaction
  interpretations explicitly.
- If the H1 gate selected a confirmation cohort (ROSMAP), repeat H3-H4 on
  it and report cross-cohort direction concordance; else log that
  confirmation is deferred (no silent cap -- state it).

**Completion note (2026-06-04).** Built two scripts. (1)
`scripts/build_human_aucell_rescoring.py` (Python/scanpy/decoupler) re-scores
the H2 signature panel with **AUCell** (`dc.mt.aucell`, the project-canonical
decoupler enrichment family, here installed into the `.venv` for the Python
side), reproducing H3's load EXACTLY (same h5ad cells, Ensembl→feature_name
relabel, JSON-bridged signatures) so only the SCORER changes, and
re-aggregates to donor×region×{all+4 states} means **over the same H3
`predicted_substate` cells** (joined by cell_id) → `human_substate_aucell_
score_means.csv` (854 rows, mirror shape of the H3 score_genes means). (2)
`scripts/build_human_robustness_mediation.R` (R; lme4/lmerTest/data.table,
no `R/helpers.R`) runs the three H5 arms → `storage/results/human_robustness_
mediation.tsv` (130 rows = 26 sigs × 5 strata), `storage/cache/summary_human_
robustness.rds` (light H6 cache: `panel`, `contrast`, full `mediation`,
`meta`), `..._provenance.txt`. All chowned rstudio:rstudio. **No Rmd touched
(chapter is H6) → no re-knit**, matching H1–H4. Both builds are deterministic
(`set.seed(1)`; AUCell rank-based) and idempotent (`--overwrite`/`--smoke`).

**ARM A — single-cell aggregation robustness.** Refit the exact H4
fixed-effect structure at the PER-CELL level with a donor random intercept
(`z(score) ~ amyloid_c*tau_c + region + age_c + sex + (1|donor)`, 80,092
staged microglia, same MIN_CELLS=10 donor×state floor). The amyloid:tau
DIRECTION reproduces in **129/130 models**; the single exception is
`DAM_down` in the DAM stratum where the effect is numerically zero
(pb +0.00072 vs sc −0.00056). Framing (anti-anchoring, Squair 2021 *Nat
Commun*): pseudobulk with donor as the replication unit stays the conservative
PRIMARY estimate; because amyloid/tau are donor-level the random intercept
makes the interaction SE roughly donor-calibrated, but cell-level residuals
are non-normal so per-cell p-values are APPROXIMATE — DIRECTION is the readout,
and it is essentially unanimous.

**ARM B — scoring-method robustness.** AUCell (rank recovery curve; no
control-set subtraction) is orthogonal to `score_genes`. Refitting the H4
model on AUCell means gives **82/90 direction-concordant** models (90 = the 18
sizeable signatures × 5 strata). 8 signatures (the 1–3-gene substate-marker and
interaction-DE sets, incl. `substate_IFN`/`substate_proliferative`) fall below
decoupler's `tmin=5` and are AUCell-n/a — a rank-recovery scorer is undefined
for tiny sets, reported explicitly (no silent cap); all SIX pre-registered
headline signatures are large (NF-κB 1304, MG-M3 740, Gsk3b 538, DAM_up,
Gerrits AD1/AD2) and tested. The 8 AUCell discordances are all small-magnitude
effects in the H3-flagged weak per-state strata (IFN/proliferative). **The two
headline attenuation mechanisms — `NFKB_union_targets` and `MG_M3_module`
negative amyloid:tau — are concordant across BOTH robustness arms in ALL FIVE
strata**; `DAM_up` (negative) and `Gsk3b_targets` (positive) are concordant in
the whole-microglia `all` stratum (Gsk3b's small effect flips only in the two
weakest per-state AUCell fits). So the H4 interaction directions are not an
artefact of donor aggregation nor of the `score_genes` method.

**ARM C — mediation (Green-2024 cascade amyloid→signature→tau).** Per
signature×stratum on the score_genes pseudobulk means: mediator
(`z(score) ~ amyloid_c + covs` → a), outcome (`tau ~ amyloid_c + z(score) +
covs` → b, c′=ADE), total (`tau ~ amyloid_c + covs` → c); ACME=a·b,
proportion mediated=ACME/c, with a **donor-clustered nonparametric bootstrap**
(B=2000; resample donors with both regions together; percentile 95% CI;
two-sided bootstrap p; BH-FDR within stratum). The disease/up signatures show
the Green-predicted POSITIVE mediation direction (a>0 amyloid drives the
signature, b>0 signature tracks tau): the human-native control
**`Gerrits_AD1_human` ACME=+0.056 [0.006, 0.123], p=0.012** (nominal), and
`interaction_up_pooled` +0.045 p=0.018; `NFKB_union_targets` +0.027 (p=0.15)
trends positive. **No signature survives FDR** (min 0.234) — expected under the
r≈0.65 amyloid↔tau collinearity and observational human power (anti-anchoring
#1/#7), so this is read as DIRECTIONAL corroboration of the cascade, not proof.
`MG_M3_module`/`Gerrits_AD2` are ~null (their amyloid main effect a≈0).

**Contrast (moderation vs mediation, made explicit).** H4's interaction is
MODERATION — signature is the OUTCOME and tau modifies the amyloid→signature
slope (the mouse 2×2 logic). H5's mediation casts the SAME signature as the
MEDIATOR transmitting amyloid→tau (the Green cascade). These are distinct
causal structures, not competing tests; `NFKB_union_targets` is the clean
illustration — NEGATIVE moderation (high tau dampens the amyloid→NF-κB slope)
coexists with a POSITIVE mediation direction (amyloid raises NF-κB, which
co-varies with tau). Both are observational/cross-sectional; mediation
additionally IMPORTS the amyloid→tau temporal ordering (Jack/Bateman cascade)
as an untestable assumption here. The `contrast` table in the rds carries both
coefficients side by side for H6.

**Deviations / decisions.** (1) Did BOTH single-cell AND AUCell (the plan's
"single-cell / AUCell" reads as OR) because they test DIFFERENT artefacts —
aggregation vs scoring method — and both clear. (2) Implemented mediation
manually with a donor-cluster bootstrap rather than the `mediation` package
(absent, and its iid assumption is wrong for the donor×region-clustered design;
the cluster bootstrap is the identifiability-aware choice and shares paired
donor resamples so ACME=a·b is a proper bootstrap draw). (3) **Env change for
future sessions:** `decoupler` was `pip`-installed into the project `.venv`
(gitignored) — H7 should note this in the reusable prompt; re-add via
`./.venv/bin/pip install decoupler` if the venv is rebuilt. (4) Confirmation
cohort DEFERRED: H1 holds no ROSMAP/Synapse DUA, so cross-cohort confirmation
is logged as deferred (not silently capped); the within-cohort MTG-vs-DLPFC
per-region replication (H4 `region_concordant`; reproduced here) is the
available independent-subsample robustness.

**Carries to H6.** `summary_human_robustness.rds`: `panel` (per
signature×stratum pb/sc/auc interaction estimates + `sc_concordant`/
`auc_concordant` flags) → a three-estimator direction-concordance panel
alongside the H4 forest; `contrast` (moderation×mediation per signature) and
`mediation` (full per-stratum ACME with CI) → a mediation/contrast panel.
Re-run (idempotent): `./.venv/bin/python scripts/build_human_aucell_rescoring.py`
then `Rscript scripts/build_human_robustness_mediation.R [--overwrite|--smoke]`.

### Session H6: render the human-validation chapter [DONE 2026-06-04]

- New `rmd/NN_human_validation.Rmd` -> new top-level section appended to
  `analysis.Rmd` before `child-session`. Panels: (a) Abeta-vs-tau
  collinearity / discordant-donor map; (b) substate-conservation cross-tab
  (+ UMAP); (c) per-signature `amyloid:tau` forest plot with VIF / coverage
  annotations; (d) the NF-kB and Gsk3b human tests vs the mouse direction;
  (e) mouse-vs-human direction-of-effect concordance summary.
- Reuse existing palettes; British English; `results='asis'` for any kable.
- Re-knit `analysis.Rmd`; verify 0 errors / 0 warnings; chown outputs.

**Completion note (2026-06-04).** Authored `rmd/18_human_validation.Rmd`
(file number continues the 01..17 sequence; renders as **§19** -- file≠§
per the standing map.md offset, anchor now recorded as 18→§19) and wired it
into `analysis.Rmd` as `child-human-validation` immediately before
`child-session`. The chapter is **display/read-only** (builds no cache; the
project convention that [S] caches are `readRDS`/`read.*`-loaded in the
consumer chunk is honoured via a guarded `hv_cache()` helper that points at
the H3/H4/H5 build scripts if a cache is missing). Full re-knit succeeded in
~4.5 min with **0 `class="error"` and 0 `class="warning"`** (verified by
grep); 5 base64 PNG figures + 7 captioned HTML tables rendered; all six
subsections (§19.1-19.6) present. All new/outputs chowned rstudio:rstudio
(`rmd/18`, `analysis.Rmd`, `analysis.html`); the knit wrote no root-owned
files under storage/.

**Panel realisation (all five plan panels delivered, plus the H5 layer the
H4/H5 carries asked H6 to surface):** §19.1 (panel a) -- the amyloid(Thal) x
tau(Braak) staging tile-map with discordant cells outlined + a per-stratum
VIF/discordant/Spearman identifiability kable (the collinearity budget made
visible: Spearman 0.649, VIF interaction 1.53-1.76, 83 discordant
donor x region obs). §19.2 (panel b) -- a UMAP coloured by predicted mouse
substate (60k of 82,486 nuclei, seed-fixed downsample for overplotting
only), the row-normalised predicted x SEA-AD-Supertype cross-tab heatmap
(disease-emergent `-SEAAD` supertypes ordered to the right so the activation
gradient reads left-to-right), and the per-state resolution **GATE** kable
(all four states resolve; monotone `-SEAAD` occupancy 0.150 -> 0.490). §19.3
(panel c) -- the 26-signature whole-microglia `amyloid:tau` forest
(geom_linerange CIs; colour=sign, diamond=pre-registered, black ring=FDR<0.05
of which there is exactly one, the DLPFC-driven proliferative score) + the
full-panel kable with FDR and region-concordance (anti-anchoring #2: the
panel is the result). §19.4 (panel d) -- a positive-control sanity kable
(Gerrits AD1 +amyloid main) + a mechanism x stratum signed-coefficient
heatmap with the mouse-predicted sign printed on each row label. §19.5 (H5) --
a three-estimator (pseudobulk/single-cell/AUCell) direction-concordance kable
and a moderation-vs-mediation contrast kable (Gerrits AD1 ACME +0.056,
nominal p=0.012). §19.6 (panel e) -- the mouse-vs-human concordance summary
kable (computed live from the caches, not hard-coded) + the verdict prose.

**Reasoned deviations / choices.** (1) **File naming:** chose `18_` (next in
the file-number run) rather than a literal `19_`; documented the §19 mapping
in map.md so future sessions are not surprised by file≠section (this matches
the existing 17→§18 precedent). (2) **Colours reuse the locked palette** --
the four-substate palette is byte-identical to rmd/11 (DAM=#B0344D etc.) and
the signed heatmaps reuse the §-divergence blue-white-red (#3a4cc0 / white /
#b40426); no new colour system introduced (carries the H1 "keep current
colours" decision). (3) **`geom_linerange` over `geom_errorbarh`:** the
latter is deprecated in ggplot2 4.0 and emitted a build-time warning even
though the global chunk option `warning=FALSE` would have hidden it from the
HTML; switched to the non-deprecated `geom_linerange` (no caps, matching the
intended `height=0`) and re-verified with `options(warn=2)` that the forest
renders to PNG warning-free. (4) **Concordance/robustness verdicts computed
in-chunk from the caches** (not transcribed) so the rendered tables cannot
drift from the H4/H5 builds. (5) **Internal session labels (H1-H6) kept out
of the rendered prose** -- the chapter reads as outward-facing science
(datasets, models, results), with cross-references only to verified mouse
section anchors (§15 Gsk3b kinase, §17 integrated model, §18 NF-kB
attenuation) and to the four reading rules distilled from the plan's
anti-anchoring guardrails (direction beats significance; a null is not
counter-evidence under collinearity; the panel is the result; corroboration
not proof). (6) **Smoke-tested every chunk** (shape/column assertions + forced
`ggplot_build` on all 5 plots) against live caches before the knit, per the
session loop.

**Headline as rendered (the chapter's scientific payload).** Conservation
GATE = PASS (4/4 mouse substates resolve in human SEA-AD microglia; monotone
disease-supertype gradient). The pre-registered mouse mechanisms reproduce
their predicted interaction sign in 5/5 strata, region-concordant and robust
to single-cell aggregation + AUCell: `NFKB_union_targets` and `MG_M3_module`
negative (tau-attenuates-amyloid), `Gsk3b_targets` positive (tau-enhanced).
No pre-registered human interaction survives FDR -- expected under r~0.65
collinearity and observational power -- so the verdict is framed as
DIRECTIONAL corroboration ("consistent with"), never proof. The human-native
Gerrits AD1 control behaves correctly (+amyloid main 5/5, +mediation ACME),
anchoring the pipeline. Re-render: standard full knit
(`Rscript -e 'rmarkdown::render("analysis.Rmd", quiet = TRUE)'`); the chapter
needs only its H3/H4/H5 caches, all already on disk.

**Carry to H7.** The two cleanest summary.Rmd panels are (i) the §19.6
mouse-vs-human concordance table (or the §19.4 mechanism x stratum heatmap)
as the headline translational result, and (ii) the §19.2 conservation gate
table/UMAP. Decide at H7 whether to add human-corroboration rows to the
biological-model ledger (`build_biological_model_ledger.R` + adjudication
re-run) or to log why not; if the venv is rebuilt, `decoupler` must be
re-`pip`-installed into `.venv` (noted at H5).

### Session H7: surface in summary.Rmd + synthesis, close plan [DONE 2026-06-04]

- Add 1-2 curated panels to `summary.Rmd` (the headline human-validation
  result + the conservation panel) as a new "Cross-species validation"
  section; re-knit summary (0 errors/warnings); chown.
- If warranted, add human-corroboration rows to the biological-model ledger
  (`scripts/build_biological_model_ledger.R` + adjudication re-run) so the
  synthesis reflects translational support; otherwise log why not.
- Move this plan to `storage/notes/completed/human_validation_plan_<date>.md`
  with an Outcome summary; write/refresh the reusable prompt if the env
  (Python venv) changes anything future sessions must know.

**Completion note (2026-06-04).** Surfaced the human-validation layer in the
curated `summary.Rmd` as a new top-level section **§7 "Cross-species
validation in human AD"** (the translational coda after the §6 "Integrated
synthesis"). Because `summary.Rmd` is a *generated* artefact
(`scripts/build_summary_rmd.R` reads chunk bodies from the gitignored staging
dir `storage/cache/summary_chunks/` and emits the committed Rmd), the section
was added the house way: two new chunk-body files
(`human_substate_conservation.R`, `human_mechanism_concordance.R`) + two panel
entries (the first carrying the section header + intro) appended to the
generator's `panels` list after `crossaxis_theme_ranking`, + two `panel_titles`
rows. Regenerated (22 panels, was 20) and re-knit `summary.html` in ~33 s with
**0 `class="error"` / 0 `class="warning"`** (grep-verified), 18 base64 figures.
The two panels are the two "cleanest" ones the H6 carry-note nominated:
(i) the **substate-conservation gate** table (`human_substate_conservation_
metrics.tsv`; all four mouse substates resolve, monotone -SEAAD gradient,
gate PASS) and (ii) the headline **mechanism × stratum** signed heatmap
(`summary_human_validation.rds`; the pre-registered interaction mechanisms
reproduce their mouse-predicted amyloid:tau sign). Both read only light
caches/TSVs (the 22 MB per-cell UMAP cache is untouched), so the summary still
knits in well under a minute. New files chowned rstudio:rstudio; smoke-tested
both chunk bodies against live caches (kable object + `ggplot_build`) before the
knit.

**Biological-model ledger decision: NOT modified (logged, per the step's
"otherwise log why not").** Adding human-corroboration rows to
`biological_model_claims_ledger.tsv` was assessed and declined for four
concrete reasons. (1) **Sealed schema.** The ledger + adjudication enforce
`axis ∈ {amyloid_activation, synaptic_suppression, interaction_metabolic,
cross_axis}` and `entity ∈ {the 11 locked IDs}` via `stopifnot` invariants on
every load; cross-species corroboration has no member axis and maps to no
entity that is not a mouse hypothesis, so admitting it would require *unlocking*
the locked vocabulary. (2) **Arithmetic perturbation of locked outputs.** Every
ledger row feeds `net_support`; the contest margins (Hyp-1B 18, Hyp-2B 12,
Hyp-3B 53) and theme net-supports (T-Synergy 28, T-Inflammation 26,
T-Compartment-suppression 17, T-Tau-attenuates 11, Hyp-0 7) rendered in §17 AND
in the summary's `biomodel_contest_verdicts` + `crossaxis_theme_ranking` panels
are functions of those counts, so any human row would silently shift the locked
margins — exactly the "locked mouse analysis" this plan's predecessor-context
forbids modifying. (3) **Evidence-class incommensurability.** The ledger grades
Strong/Moderate/Suggestive by replication across the four *mouse* molecular
modalities/layers; observational, cross-sectional, FDR-null, collinearity-
limited human corroboration (anti-anchoring #7: "consistent with, not proof")
is a categorically different evidence type that counting as Strong/Moderate
would mis-weight and inflate mouse margins with non-commensurable evidence.
(4) **Already surfaced at the right altitude.** The translational layer is
delivered as its own outward-facing chapter (§19) plus the two new §7 summary
panels, with explicit "consistent with, not proof" framing; cross-reference
directionality is correct (§19 → the locked mouse sections §15/§17/§18, never
the reverse). Verified the ledger currently holds zero human/cross-species rows
(`grep -ci 'human|seaad|braak|thal' = 0`) so this is the standing state, not a
removal. The decision was made autonomously because the step delegates it
("if warranted … otherwise log why not") and it is not a flagged decision gate;
it is trivially revisitable from this note if the user wants the synthesis to
fold human evidence into the mouse arithmetic.

**Plan closure.** Marked H7 DONE; added the Outcome summary below the title;
`git mv`'d the plan to `storage/notes/completed/human_validation_plan_
2026-06-04.md`; added a `### human_validation_plan_2026-06-04` digest to
`completed/DIGEST.md` and retired the "active successor" line in its narrative
arc; updated `map.md` (summary_human_* caches now also feed summary.Rmd §7;
the human lane is marked completed). `reusable_prompt.md` left **unchanged**:
H7 introduced no environment change (the scanpy/decoupler `.venv` and the
SEA-AD h5ads were established at H1 and are already gitignored + documented in
map.md), and the boot prompt is deliberately plan-agnostic — with no active
plan on disk, the next session correctly falls through to its "survey and
propose" branch.

**Reasoned deviations / choices.** (1) **Three interaction mechanisms in the
headline heatmap, not four.** The §19.4 chapter heatmap shows NF-kB, MG-M3,
Gsk3b *and* DAM-up; the summary headline shows only the three whose
*pre-registered* prediction is about the interaction term (NF-kB & MG-M3
negative, Gsk3b positive — all 5/5 strata, unambiguous). DAM-up's locked
prediction is a positive *amyloid main-effect* control (its interaction sign is
incidental and ≈0 in proliferative), so it belongs to the sanity-control set,
not the interaction headline; the full 26-signature panel + DAM-up + controls
remain in analysis.html §19, so this is curation for clarity, not selective
reporting (anti-anchoring #2 is satisfied at the §19 altitude). (2) **Gate
table + mechanism heatmap** chosen over the §19.6 concordance summary table or
the §19.2 UMAP — the table-gives-verdict + heatmap-gives-direction pairing is
the crispest one-table/one-figure statement and avoids the heavy per-cell cache.
(3) **Self-contained chunks** — each summary chunk redefines its palette /
strata labels / mouse-direction map inline (the summary setup only sources
helpers + defines `rd()`/`rt()`), matching every existing curated chunk.

## Anti-anchoring (re-read every session)

1. **A null human interaction is NOT evidence against the mouse finding.**
   Report identifiability (VIF, off-diagonal donor coverage, power) with
   every `amyloid:tau` estimate. Collinearity can mask a real interaction;
   say so before interpreting any null.
2. **No cherry-picking.** Report the FULL signature panel with FDR,
   including signatures that fail. The result is the panel, not the best row.
3. **Conservation precedes per-state claims.** Establish that DAM / IFN /
   proliferative resolve in human microglia (H3) BEFORE making any
   per-state interaction claim; otherwise restrict to whole-microglia.
4. **Direction beats significance.** A significant human interaction in the
   OPPOSITE direction to mouse refutes, not supports. Concordance of sign
   is the primary readout.
5. **Sanity-gate on positive controls.** Human-native Gerrits AD1/AD2 and
   the DAM-up-vs-amyloid main effect must behave as expected; if they do
   not, treat the whole pipeline as suspect before trusting novel rows.
6. **Keep amyloid and tau separate.** Never collapse to a single
   AD-vs-control axis -- that discards exactly the tau/amyloid separation
   that makes this a cross-species test of the INTERACTION.
7. **This is corroboration, not proof.** Human cohorts are observational,
   cross-sectional and confounded; frame outcomes as "consistent with" /
   "not consistent with" the mouse mechanism, never as causal confirmation.

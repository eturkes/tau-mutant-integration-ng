# P3 Mechanism - plan

## Scope
Explain the mechanism behind the rebuilt P1/P2 result: amyloid drives microglial
homeostatic->DAM activation, and mutant tau modulates that response mainly by
COMPOSITION (more DAM cells), not by supported further progression along the activation
ordering. P3 asks which pathway / TF / kinase signals plausibly sit under that
state expansion, with explicit checks for the v1 headline mechanisms:

- TF / pathway: Myc-family interaction signal; amyloid inflammatory/APC/DAM axis; NF-kB attenuation.
- Kinase: Gsk3b interaction support from phosphoproteomics.
- Honest interpretation: RNA-derived TF/pathway = microglia-specific static expression; kinase =
  24M bulk hippocampus phosphosite support, NOT microglia-sorted.

Inputs already built: `pb_de_microglia`, `pb_de_substate`, `composition_results`,
`trajectory_report`, `symbol_map`, `phospho`, `sample_key`. Outputs: focused
mechanism targets + `_mechanism.qmd` report chapter. Full GeoMx / proteome / phospho
cross-modality interpretation stays P4; P3 may build a minimal phosphosite-DE leaf ONLY
because kinase inference requires site-level statistics.

## Research Digest
Local / v1:
- v1 mechanism layer used `decoupleR + CollecTRI` for TF, `decoupleR + OmniPath KSN`
  for kinase, and a separate CCC layer. Top v1 findings: Gsk3b kinase at interaction,
  Myc/Creb1/Tp53/Jun TFs at interaction, Spi1/Nfkb1/Sp3 at amyloid, and a later per-state
  NF-kB attenuation check.
- v1 bloat to drop: CCC (belongs P4 synaptic/clearance axis if it earns it), claim ledger,
  CARNIVAL topology, SCENIC, human validation, cross-cell-type specificity, tradeSeq dynamics.
- v1 gotcha: `decoupleR::get_ksn_omnipath()` previously called stale OmnipathR internals.
  P3 must smoke-test the current API under this repo lock before trusting it.

Current-method sweep (2026-07-02):
- Bioconductor 3.23 matches this repo's R 4.6 stack; `decoupleR` release is 2.17.0 and
  explicitly supports signed/weighted network activity over transcriptomics and
  phosphoproteomic kinase-substrate features.
  Source: https://bioconductor.posit.co/packages/3.23/bioc/html/decoupleR.html
- `decoupleR::decouple()` remains the fit-for-purpose wrapper for multiple statistics with
  optional consensus score; `run_ulm()` returns source x condition activity scores plus p-values.
  Sources: https://saezlab.github.io/decoupleR/reference/decouple.html,
  https://saezlab.github.io/decoupleR/reference/run_ulm.html
- CollecTRI remains the best on-lock TF prior here: signed TF-target edges, expanded from DoRothEA,
  organism arg supports mouse; paper reports 1186 TFs / 43175 signed interactions and stronger
  perturbation performance than comparator curated collections.
  Sources: https://saezlab.github.io/decoupleR/reference/get_collectri.html,
  https://academic.oup.com/nar/article/51/20/10934/7318114
- OmniPathR release 4.0.0 is on Bioc 3.23; its `enzyme_substrate(organism=10090)` path says
  mouse/rat enzyme-substrate orthology is available directly. Prefer this over v1's off-lock
  `nichenetr` mapping. Smoke-test required, and the web-service prior must be treated as a
  reproducibility-sensitive input: project-local cache + hashes/counts/query args + package versions.
  Sources: https://bioconductor.org/packages/release/bioc/html/OmnipathR.html,
  https://r.omnipathdb.org/reference/enzyme_substrate.html,
  https://r.omnipathdb.org/reference/omnipath_set_cachedir.html
- Mouse MSigDB now has native mouse collections (MGI symbols); use focused GO/curated subsets, not
  the full all-sets bundle. `fgsea` is Bioc 3.23 and supports fast preranked GSEA.
  Sources: https://www.gsea-msigdb.org/gsea/msigdb/mouse/collections.jsp,
  https://bioconductor.posit.co/packages/3.23/bioc/manuals/fgsea/man/fgsea.pdf,
  https://cloud.r-project.org/web/packages/msigdbr/msigdbr.pdf

## Default Design
Lean, on-lock, mechanism-specific:

1. RNA pathway / TF activity from the existing microglia pseudobulk statistics.
   Unit = replicate-correct pseudobulk DE already built in P1. Whole microglia + fit substates
   (Homeostatic, DAM); skipped substates remain skipped. Build ranked gene-symbol t matrices from
   `pb_de_microglia` / `pb_de_substate` via `symbol_map`; duplicate symbols keep max |t|.

2. NF-kB attenuation is a targeted RNA check, not a new broad arc.
   Define NF-kB family from CollecTRI sources matching `Rela`, `Nfkb1`, `Nfkb2`, `Rel`, `Relb` plus
   any complex/accession source returned by the current mouse prior. Read both:
   - PRIMARY: `interaction` contrast = tau modifies the amyloid response.
   - SUPPORTIVE: `tau_in_nlgf` contrast = mutant tau effect on the amyloid background.
   Verdict gate = negative `interaction` in the primary whole-microglia RNA family with FDR<0.10
   after small-family correction across the predeclared TF-family and NF-kB-target GSEA tests.
   `tau_in_nlgf` and substate rows can corroborate localisation but cannot by themselves prove
   attenuation of the amyloid response. If only directionally consistent, report as suggestive.

3. Kinase activity uses a minimal 24M phosphosite-DE target.
   `phospho` raw target stores 67 runs; `sample_key` defines the first 16 24M balanced samples.
   P3 wires: match 16 intensity columns -> log2 transform -> median normalise -> prevalence filter
   -> `factorial_design(add_batch=FALSE)` -> `fit_limma_log()` -> site-level top tables. This is
   a narrow kinase substrate for P3, not P4's bulk phosphoproteomics chapter. The sample-key order
   is genotype-blocked (01-04 MAPTKI, 05-08 P301S, 09-12 NLGF_MAPTKI, 13-16 NLGF_P301S), so kinase
   claims carry a load-bearing run-order caveat and a run-index trend sensitivity where estimable.

4. Kinase inference = `decoupleR` over an OmniPath mouse KSN.
   Build KSN from the official OmniPath mouse REST `enz_sub` endpoint by default; keep
   `try_package=TRUE` as an explicit drift probe for `OmnipathR::enzyme_substrate(...)` because S1 proved
   the package postprocessor currently fails on the live schema. Filter phosphorylation /
   dephosphorylation, validate exact columns before use, drop unresolved conflicting signed pairs, and
   produce `source` kinase + `target` site id `SYMBOL_AApos` + `mor`. Because OmniPath uses UniProt IDs
   as primary identifiers and gene symbols are optional, S1 must assert mouse-symbol presence, symbol
   case, and real overlap against current phosphosite IDs before any activity score. Require smoke-test
   coverage before the real target: nonempty network, >=50 kinases with minsize>=5 on current phospho
   sites, `Gsk3b` present either as source or explicitly absent-with-count recorded. If direct mouse KSN
   fails coverage, stop before implementation-choice drift and decide fallback.

5. Report = one compact target + one chapter.
   `mechanism_report` bundles only small tables / plot frames. `_mechanism.qmd` loads that target and
   no heavy Seurat object. Prose inline-computed; no hardcoded p-values. Chapter order:
   pathway survey -> TF activity -> NF-kB attenuation -> Gsk3b kinase -> synthesis + caveats.

## Alternatives
Alternative A - RNA-only P3, defer kinase to P4.
Pros: cleaner modality boundary; no minimal bulk-phospho leaf before P4. Cons: fails the explicit P3
Gsk3b kinase objective and leaves the mechanism headline half rebuilt.

Alternative B - v1-style full mechanism arc now (TF + kinase + CCC + topology).
Pros: more orthogonal corroboration in one phase. Cons: reintroduces the exact v1 bloat the rebuild is
trying to delete; CCC belongs with P4's synaptic/clearance cross-modality story; topology/SCENIC/human
layers were mostly margin-neutral or over-scoped.

Default choice: lean P3 with minimal phosphosite kinase support. It answers the backlog's Gsk3b/Myc/NF-kB
mechanism request while preserving P4's larger cross-modality remit.

## Steps
Each step is one closing unit. Resuming mid-plan: read this plan's Scope + your step + `.agent/memory.md`
relevant sections + files named in the step. Run `scripts/check.sh` unless explicitly docs-only.

### S1 - Dependencies + API contracts [DONE 2026-07-02]
Edit `rproject.toml` and `rv.lock` only through `rv sync` after adding:
`decoupleR` (BioCsoft), `OmnipathR` (BioCsoft), `fgsea` (BioCsoft), `msigdbr` (CRAN),
`digest` (CRAN; prior hashes).
Add `R/mechanism.R` with pure API/shape helpers and `tests/test_mechanism.R` synthetic tests.

Outcome note: Bioc 3.23 `OmnipathR` loads warning-clean after `TZ=UTC`, but its `collectri()` /
`enzyme_substrate()` postprocessor currently fails on the live server schema (`ncbi_tax_id` missing).
S1 therefore uses official OmniPath REST endpoints by default (`try_package=FALSE`), with
`try_package=TRUE` left as an explicit drift probe. This is still direct mouse OmniPath, not v1's
off-lock `nichenetr` mapping.

Contracts:
- `add_symbol_to_top(top_df, symbol_map, gene_col="gene")` -> top table with `symbol`; fail loud on
  missing gene col; drop unmapped only after recording attr/count.
- `extract_rank_matrix(top_list, symbol_map, stat_col="t")` -> symbol x contrast numeric matrix;
  duplicate symbols keep max |stat|; names match canonical contrasts present.
- `run_decoupler_matrix(mat, network, minsize=5)` -> canonical long tibble
  `{statistic, source, condition, score, p_value}`; current API smoke-tested on toy data; if consensus
  output changes, wrapper records `has_consensus=FALSE` and downstream uses ULM score for direction.
- `set_mechanism_prior_cache()` -> `OmnipathR::omnipath_set_cachedir("storage/cache/omnipath")`
  (gitignored project-local cache, not user-global); records cache path in provenance.
- `prior_fingerprint(x, query)` -> stable hash of sorted prior tibble + query args + package versions.
  S1 records observed CollecTRI/KSN hashes and counts in code/test fixtures or memory; later builds call
  `assert_mechanism_prior_expectations()` to fail loudly on unexpected drift rather than silently
  refreshing priors.
- `load_collectri_mouse()` -> returns signed mouse CollecTRI with `{source,target,mor}`; records edge /
  TF / target counts, source examples, query args, package version, hash.
- `load_omnipath_ksn_mouse()` -> returns `{source,target,mor}` site network via official OmniPath mouse
  REST by default; explicitly requests/validates gene symbols (`enzyme_genesymbol`,
  `substrate_genesymbol`, `residue_type`, `residue_offset`, `modification` or documented current
  equivalents). `try_package=TRUE` probes the current OmnipathR path; no v1 `nichenetr` fallback.
- `phospho_site_ids(phospho_tbl)` -> real current phospho target IDs from `{PG.Genes, PTM.SiteAA,
  PTM.SiteLocation}` with missing/multi-gene/drop counts; used in S1 for KSN overlap before S3 DE.
- `ksn_coverage_probe(ksn, phospho_site_ids, minsize=5)` -> matched site count, kinases passing
  minsize, `Gsk3b` coverage, symbol-case diagnostics.

Acceptance:
- `rv sync` green; API smoke proves package load + output schemas under warn=2.
- Unit tests cover duplicate-symbol collapse, missing columns, decoupleR schema, CollecTRI shape, KSN
  target-id parsing on a synthetic OmnipathR-like table, prior fingerprint determinism.
- Live S1 smoke runs KSN coverage against current phospho site IDs and fails here, not in S3, if direct
  mouse OmniPath coverage is unusable.
- `scripts/check.sh` green.

### S2 - RNA pathway + TF + NF-kB targets
Add targets:
- `mechanism_gene_sets` = focused mouse gene sets: GO BP / GO CC / GO MF from native Mouse MSigDB via
  `msigdbr(db_species="MM", species="Mus musculus")`, plus compact project sets: DAM, Homeostatic,
  MHC_APC, IFN, NF-kB target union from CollecTRI family sources. Apply min_size 5 for custom, 15 for GO.
- `mechanism_tf` = decoupleR CollecTRI activity on whole microglia + fit substates across 5 contrasts.
- `mechanism_pathway` = fgsea preranked GSEA on the same RNA ranked matrices.
- `nfkb_attenuation` = targeted table for NF-kB-family TF activity + NF-kB-target GSEA in `interaction`
  and `tau_in_nlgf`, by whole microglia + fit substates.

Contracts:
- Pseudobulk is the inference unit; no cell-level TF tests.
- FDR: ULM p BH within population x contrast; fgsea padj within population x collection x contrast.
  NF-kB attenuation verdict additionally applies a small primary-family BH across the two primary
  whole-microglia tests (TF family, target GSEA) at `interaction`; substate/tau_in_nlgf rows are supportive.
- Direction: decoupleR consensus score if present; ULM score fallback recorded. fgsea direction = NES.
- Fit substates only: Homeostatic and DAM unless current target says otherwise; IFN/proliferative skipped
  status travels into the result.

Acceptance:
- Fresh target builds warning-free; `tar_meta` clean for new targets.
- Unit tests cover gene-set size filtering, NF-kB source extraction, FDR/direction convention, skipped
  substate propagation, and the rule that `tau_in_nlgf` alone cannot call attenuation supported.
- Live smoke prints top TFs/pathways for amyloid, interaction, tau_in_nlgf; no hardcoded expected winner.
- `scripts/check.sh` green.

### S3 - Minimal phosphosite DE + kinase activity
Add targets:
- `phospho_de_24m` = 16 matched 24M phosphosite limma-trend fit with 5 contrasts, no batch term.
- `kinase_activity` = decoupleR activity over the mouse OmniPath KSN, per contrast.
- `kinase_mechanism_summary` = compact table centred on significant per-contrast kinases plus explicit
  Gsk3b row whether significant or not.

Contracts:
- Column matching asserts exactly 16/16 sample-key runs found, 4/genotype balanced, no duplicate stubs.
- Transform: log2(raw intensity) after positive-value guard; nonpositive -> NA and counted.
- Feature/unit contract:
  - DE input rows stay traceable by `PTM.CollapseKey` + original row index; a non-unique collapse key fails
    loud unless disambiguated.
  - Biological kinase target ID = `paste0(PG.Genes, "_", PTM.SiteAA, PTM.SiteLocation)` after dropping blank
    symbols and multi-gene groups with recorded counts (multi-gene rows cannot be assigned to one KSN target).
  - Before kinase activity, collapse duplicate biological site IDs to ONE statistic per target by highest
    `Phosphosite probability` (typo variant accepted) and then max |t| tie-break; record rows/site distribution.
  - decoupleR never sees duplicated KSN target IDs.
- Design: `factorial_design(meta, add_batch=FALSE)` because bulk sample key has no batch. State in prose:
  bulk hippocampus, not microglia-sorted. Add run-index sensitivity: build an additive scaled run-index design
  if full-rank and the interaction is estimable; if Gsk3b support depends on the unadjusted fit, downgrade the
  claim to run-order-confounded.
- KSN coverage gate: record network source count, matched site count, kinases passing minsize, Gsk3b
  coverage. A failed coverage gate stops the target with a precise message.

Acceptance:
- Fresh `phospho_de_24m` + `kinase_activity` builds warning-free.
- Unit tests cover 16-column match, log transform guards, site-id construction, duplicate-site collapse,
  multi-gene/blank-symbol dropping, add_batch=FALSE design, run-index sensitivity rank/estimability,
  KSN coverage summariser, explicit Gsk3b carry-through.
- Live smoke prints interaction/tau_in_nlgf kinase table; `scripts/check.sh` green.

### S4 - Mechanism report + integration
Add:
- `mechanism_report_data(tf, pathway, nfkb, kinase, composition_results, trajectory_report)` -> compact
  target `mechanism_report`.
- `_mechanism.qmd`, included after `_trajectory.qmd`.
- Map/memory updates for new wiring + durable gotchas.

Report must state:
- Myc/TF and Gsk3b/kinase are mechanism hypotheses supported or not by the rebuilt data, not assumed winners.
- NF-kB attenuation is transcript-level, not composition, and only "supported" if the S2 gate clears.
- Kinase evidence is bulk hippocampus phospho, not microglia-sorted.
- P2 revised the mechanism question: explain DAM-cell expansion / amyloid-response modulation, not
  progression acceleration.

Acceptance:
- `_mechanism.qmd` renders 0-warning under `options(warn=2)`.
- `mechanism_report` is compact (<5MB expected) and qmd never loads `microglia_annotated`.
- Full `scripts/check.sh` green.
- Codex review run on uncommitted P3 work; accepted findings fixed before S4 commit.

## Close-Out
After S4:
- Adversarially review plan vs shipped code/prose.
- Fold P3 digest into `.agent/history.md`; archive plan to `.agent/completed/`.
- Reset roadmap Active plan; update spine wording if Gsk3b/Myc/NF-kB support changes.
- Commit `mechanism (p3 close): ...`.

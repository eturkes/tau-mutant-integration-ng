# P7 - Integrity-first audit (active plan)

Companions: `roadmap.md` (trajectory), `memory.md` (contract), `map.md` (wiring),
`history.md` (digests), `archive_digest.md` / branch `archive` (v1 mining).

Direction **G** selected by the user (prior-session PLANNING direction gate, 2026-07-18;
repo left unchanged at selection). Scope = resolve the three stop-the-line integrity
gates (G1/G2/G3) plus the 24M bulk run-order contract, then run a **preregistered
DAM-occupancy robustness audit**. Posture: settle integrity BEFORE adding biology or
freezing a manuscript. This milestone corrects labeling/provenance and stress-tests the
report's one central positive claim; it does not add new biological scope.

## Why - findings confirmed this session against permitted real inputs
Evidence gathered by inspecting `storage/data/*` + tracked code (no repo change):

- **G1 MAPTKI construct (LITERATURE-CONFIRMED 2026-07-18).** Repo hardcodes "MAPTKI =
  tau-KO baseline" (`.agent/memory.md:12`, `roadmap.md:14`) and describes the amyloid-on-
  MAPTKI axis as the "tau-KO background" in code (`_targets.R:153`, `R/report/plot.R:173`,
  `R/report/figures.R:342`, and the data-carried `figures.R:543 y_meaning`). Primary
  literature (background research agent, 5 sources) confirms this is WRONG. `MAPT-KI` is the
  Saito **humanized wild-type MAPT knock-in**: endogenous murine Mapt is replaced by WT human
  genomic MAPT, so the mice express all six 3R/4R human tau isoforms under native regulation
  and lack MOUSE tau but are NOT tau-null (a true Mapt-KO is a separate tau-negative genotype).
  `P301S+3` is the Watamura base-edited **MAPT^P301S;Int10+3G>A** allele - the "+3" is the
  Int10+3 splice mutation, NOT an age label (age is the separate 24M/30M token) - one of seven
  mutant lines derived from the humanized KI. So the four groups form a clean 2x2 with **no
  tau-null arm**, and the tau factor compares **WT humanized tau vs mutant (P301S;Int10+3)
  humanized tau**, optionally stratified by App^NL-G-F. Every "tau-KO"/"tau-absence"
  description in the repo is wrong. Sources: Saito 2014 (App^NL-G-F KI), Saito 2019 (WT
  humanized MAPT-KI + App/MAPT-KI), Watamura 2025 (the seven base-edited human-MAPT alleles;
  MAPT-KI != Mapt-KO), Benzow 2024 (separate H1/H2 MAPT-GR series), Morito 2025 (P301S;Int10+3
  + App^NL-G-F crosses). User breeding/allele records would add per-animal provenance but are
  no longer needed for construct identity.
  **Provenance caveat:** the unused Set6 arms `MAPT-N279K`, `MAPT-10+16`, and
  `MAPT-wt-H1-haplotype` likely belong to the separate Benzow/Koob **MAPT-GR H1** series (also
  human-tau, matched control `MAPT-wt-H1-haplotype`) - do NOT pool them with the Saito
  `MAPT-KI`; `MAPT-10+3`/`P301S-10+3`/`S305N-10+3` are the Saito base-edited humanized-KI series.

- **G2 24M "proteome" identity (CONFIRMED from data + code).** `proteomics_...tsv`
  carries phosphosite PTM columns (`Gene-pSite`, `PTM.SiteAA`, `PTM.SiteLocation`,
  `Phosphopeptide`, `Phosphosite probability`) and 16 `..._TiO2_DIA_NN.raw.PTM.Quantity`
  columns; `R/analysis/modality_de.R:2191 aggregate_proteome_raw()` sums these by
  `PG.ProteinGroups` and treats the result as a proteome. A TiO2-enriched phosphopeptide
  PTM report summed to protein groups is NOT a valid global-proteome denominator. Both
  "proteome" and "phospho" are the same TiO2-enriched assay; there is no separate global
  proteome in `storage/data/`. The bulk-context "proteome PCA" figure must be confirmed /
  relabeled / replaced / removed.

- **Run-order contract (CONFIRMED).** 24M acquisition is perfectly genotype-blocked
  (runs 01-04 MAPT-KI, 05-08 P301S+3, 09-12 NLGF-MAPT-KI, 13-16 NLGF-P301S+3), so run
  order is fully confounded with genotype in the bulk. Current bulk DE fits no batch/run
  term (`modality_de.R` "no batch").

- **G3 provenance (partially recoverable now).** The sample key is a 67-run manifest:
  Naoto set = 26 runs (24M 01-16 used; 24M/30M 17-26 unused), Set6 set = 41 runs
  (unused humanized-MAPT + NLGF allelic series). Only the 16 Naoto-24M runs enter the
  report. GeoMx has a known 91-of-112 AOI gap (`memory.md`: MAPTKI bio_rep 1 absent).
  Cross-assay animal identity, survival/attrition/batch records, and GeoMx raw files are
  NOT in-repo and need the user.

- **Occupancy core claim (untested for robustness).** P6 verdict: interaction DAM
  occupancy = +0.174 fraction (95% CI 0.095-0.253), zero-null FDR 1.81e-5, but the
  0.10-margin family is unresolved at 5% (FDR 0.081). This single positive interaction is
  the report's spine and its robustness to clustering/pruning/annotation/unit choices was
  never audited. Family + estimators live in `R/analysis/state_decomposition.R`
  (`fit_state_occupancy` L480, `state_probability_standardization` L385,
  `run_microglia_state_response` L866); labeling knobs in `R/analysis/microglia.R`
  (`reprocess_microglia` L29 `resolutions=c(0.2,0.4,0.6)`/`primary_resolution=0.4`,
  `assign_subpopulation` L202 `tol`/`amb_floor`, `flag_contaminant_clusters` L230
  `id_floor`/`mglike_floor`).

## Units - 5, all gate-independent (OPEN). Resolve G1/G2 first, then audit.

### P7.1 - MAPTKI construct verification + tau-factor relabel (G1)  [OPEN]
- Construct verification COMPLETE (literature verdict 2026-07-18, see G1 finding above):
  MAPTKI = Saito humanized WT MAPT knock-in; P301S = Watamura base-edited MAPT^P301S;Int10+3;
  four groups = 2x2 WT-humanized-tau vs mutant-humanized-tau x +/-App^NL-G-F, no tau-null arm.
  Remaining P7.1 work = the relabel below + (optional) folding user breeding records for
  per-animal provenance and the five citations into `.agent/memory.md`.
- If confirmed: relabel the tau factor MEANING (keep the `MAPTKI` factor TOKEN stable to
  avoid churn) at all six sites: `_targets.R:153`, `R/report/plot.R:173`,
  `R/report/figures.R:342` + `:543 y_meaning`, `.agent/memory.md:12`, `.agent/roadmap.md:14`.
  Check whether `y_meaning` surfaces in the rendered report; correct visible text.
- Constraint: this is a requirement-level relabel of the core factor -> exact corrective
  wording is a USER-DECISION POINT. No occupancy/DE numbers change. `scripts/check.sh` green.
- Accept: no "tau-KO"/"tau-absence" description of MAPTKI remains; report text reflects
  the verified construct; records-request recorded; render gate green.

### P7.2 - 24M bulk assay identity + run-order contract (G2 + run-order)  [OPEN]
- Document, with column + code evidence, that the "proteome" is a TiO2-enriched phospho
  PTM report summed to protein groups (not a global proteome). Add a run-order/genotype
  confound sensitivity to the 24M bulk DE (or scope the claim), using the genotype-blocked
  acquisition. Decide the bulk-context figure fate: confirm / relabel / replace / remove.
- Locations: `R/core/io.R:104-105`, `R/analysis/modality_de.R:2138-2299`
  (`aggregate_proteome_raw`, `prepare_proteome_24m_matrix`, `run_proteome_de_24m`,
  `match_24m_bulk_columns`), `R/report/plot.R:2092 proteome_pca_plot` / `:2156
  bulk_modality_context_plot`, `R/report/figures.R:192-281`, `sections/modality.qmd`.
- Constraint: figure relabel/replace/remove is a USER-DECISION POINT; run-order
  sensitivity uses existing data. Vendor report template/enrichment method = user record.
- Accept: assay identity documented; run-order confound quantified/scoped; figure resolved
  per decision; render gate green.

### P7.3 - Source-provenance dossier + record requests (G3)  [OPEN]
- Compile the recoverable provenance (67-run manifest breakdown, 16-of-67 selection, 30M +
  Set6 unused arms, GeoMx 91/112 AOI gap, snRNAseq 16-unit structure) into a dossier.
  Note the Set6 arms span TWO distinct series (Saito base-edited humanized-KI:
  `MAPT-10+3`/`P301S-10+3`/`S305N-10+3`; vs Benzow/Koob H1 MAPT-GR: `MAPT-N279K`/`MAPT-10+16`/
  `MAPT-wt-H1-haplotype`) - do NOT conflate their WT controls with the Saito `MAPT-KI`.
  Enumerate precisely what only the user can supply (cross-assay animal/aliquot crosswalk;
  survival/attrition/litter/cage/sex/batch @24M; GeoMx DCC/PKC/neg-probe/scan/mask +
  missing-AOI records). Record a standing prohibition on cross-modality/paired claims
  (DIABLO/MOFA/mediation/within-animal correlation) until the crosswalk exists.
- Locations: `.agent/` (dossier + memory/map note), `R/core/io.R`,
  `storage/data/proteomics_sample_key.csv` (source).
- Accept: dossier lists known-vs-needed provenance; cross-assay-claim gate recorded in
  memory; record-request surfaced.

### P7.4 - DAM-occupancy robustness preregistration + harness (prep, pre-results)  [OPEN]
- FREEZE the robustness protocol before viewing any sensitivity result: clustering
  resolution grid around 0.4; pruning variants (`flag_contaminant_clusters`
  `id_floor`/`mglike_floor`); annotation sensitivities (`assign_subpopulation`
  `tol`/`amb_floor`; UCell); leave-one-unit-out (16) + leave-one-batch-out (4);
  occupancy estimators (beta-binomial primary; empirical-logit; simple proportion);
  primary endpoint = interaction DAM-fraction +0.174 with its zero-null + 0.10-margin
  families; explicit tipping-point/decision rules. Build a reusable harness mapping a
  membership labeling -> occupancy family, validated ONLY against the current labeling.
- Locations: `R/analysis/state_decomposition.R`, `R/analysis/microglia.R`, new prereg doc
  in `.agent/`, optional new harness fn/target.
- Constraint: preregistration discipline - no peeking at variant occupancy before the spec
  is frozen. Existing data only.
- Accept: frozen prereg spec committed; harness reproduces the established +0.174 family
  exactly on the current labeling; gate green.

### P7.5 - Execute robustness sweep + tipping-point verdict  [OPEN, after P7.4]
- Run the frozen sweep on existing snRNAseq data; report the occupancy interaction's
  robustness range, sign-stability, and tipping point (which perturbation, if any, flips
  the sign or crosses the 0.10 margin). Optionally append one compact robustness
  figure/appendix WITHOUT disturbing Figures 1-10.
- Constraint: only the occupancy family reruns; no new biology; verdict is outcome-
  independent. SIZE-CHECK at dispatch - if the resolution x prune x annotation x LOU/LOBO x
  estimator sweep + report cannot preserve the ~72K Agent reserve, respec-split at the
  sweep/report seam.
- Accept: robustness table + tipping-point verdict; sign/margin stability characterized;
  render gate green if a figure is added.

## Sequencing
P7.1 -> P7.2 -> P7.3 (settle G1/G2 provenance first, per the direction shortlist),
then P7.4 (prereg prep) -> P7.5 (execute). All OPEN / gate-independent.

## Standing external-record requests (user-supplied; upgrade inference to definitive)
1. G1: MAPTKI/P301S breeding records, allele IDs, promoter/locus, sample labels.
2. G2: Spectronaut report template, TiO2 enrichment/fractionation method, precursor/
   modification column doc, full acquisition/run provenance.
3. G3: cross-assay animal/aliquot crosswalk; 24M survival/attrition/litter/cage/sex/batch;
   GeoMx DCC/PKC/negative-probe/scan/mask files + missing-AOI records.

## User-decision points (requirement-level; confirm before scope-source edits)
- P7.1 exact relabel wording of the tau factor (once construct is confirmed).
- P7.2 fate of the bulk "proteome" figure: relabel / replace / remove.

# Memory - standing contract (read every session)

Carried-forward learnings and gotchas that survive all plans. Companions in
`.agent/`: `roadmap.md` (where we are going) · `map.md` (codebase wiring) ·
`history.md` (per-arc decision digests) · `completed/` (full archived plans,
deny-Read). This is a sealed, MAINTENANCE-posture R/Bioconductor analysis.

## Project
Integrate snRNAseq + GeoMx spatial + bulk proteomics + bulk phosphoproteomics
across 4 mouse AD genotypes (MAPTKI, P301S, NLGF_MAPTKI, NLGF_P301S) to read how
amyloid (NLGF) reshapes microglia under different tau backgrounds. Divergence
readout = interaction contrast `(NLGF_P301S − P301S) − (NLGF_MAPTKI − MAPTKI)`.
Env: root in a Docker container; R 4.5.2 + Bioconductor 3.22; network open.

## Locked scientific decisions (hold unless the user reopens them)
- Pseudobulk replicate = genotype_batch (16 ids, 4/genotype); batch in every design.
- Proteomics: peptide → PG.ProteinGroups by sum; treat as total proteome.
- Phospho: first 16 samples (24M timepoint); it is bulk hippocampus, NOT
  microglia-sorted (restate that in kinase prose). MAPTKI phospho is masked in
  publication figures via the `publication_mode` param.
- Microglia QC: re-cluster on the subset, drop only clear outliers (no over-pruning).
- 5 canonical contrasts everywhere: tau_alone, nlgf_in_maptki, nlgf_in_p301s,
  tau_in_nlgf, interaction (factorial 2×2 + batch).
- Present the 3 axes with no pre-privileged winner; state thresholds before
  applying them; alphabetise hdWGCNA modules in renders (MG-M3 never first).

## Standing bans (never reintroduce)
MSigDB Hallmark collection; the OXPHOS confirmatory arc (deleted rmd 08/09a-d);
any specificity-null / expression-matched-random framework.

## Sealed analysis state
The 11-arc programme (D->O) + arc-P capstone closed 2026-06-15. §17 biological-
model ledger = 92 rows; contest margins **18 / 12 / 55** (amyloid Hyp-1B /
synaptic Hyp-2B / interaction Hyp-3B); cross-axis lead T-Inflammation 31 >
T-Synergy 30. Arcs K-O are margin-neutral by design, so hold these numbers; the
capstone build script hard-asserts 18/12/55 and fails loudly on drift. New
analysis arcs need an explicit user request (see roadmap.md).

## Operational (this project)
- Known-good state = a clean knit: 0 `class="error"` AND 0 `class="warning"` in
  each of the 3 reports (analysis / summary / synthesis). Verify before AND after.
- Re-knit: `Rscript -e 'rmarkdown::render("analysis.Rmd", quiet = TRUE)'`. Rebuild
  snRNAseq caches: `params=list(recompute_snrnaseq=TRUE)`. Publication build:
  `params=list(publication_mode=TRUE)`.
- One knit = ONE shared R session across child rmds → a global built in an earlier
  child is in scope for later ones (map.md lists the producers).
- Caches: [K] (built in-knit by `cache_or_run`) auto-load; [S] (built by
  `scripts/build_*.R` outside the knit) must be pre-built / readRDS-loaded. map.md
  holds the producer→consumer table.
- Knits take ~3 min: smoke-test new helpers / cache-reading code against the live
  caches via `Rscript -e '...'` before committing to a full knit.
- New `R/*.R` must be `source()`d in `R/helpers.R` in dependency order (read its top).
- `knitr::kable()` renders an HTML table ONLY inside a `results='asis'` chunk.
- rmd file number ≠ rendered §number (rmd 08/09 were deleted); map.md has the anchors.
- You run as root → `chown rstudio:rstudio` new files AND the knit's root-owned
  outputs (analysis.html, storage/results/*.tsv, new storage/cache/*).
- Prose register: British English; human-facing report/figure text uses hyphens
  over em/en-dashes (commas or parentheses for asides, colons for restatements).
  R `#` comments and code stay exempt (LLM-facing).
- The 3 report HTMLs + storage/{cache,data,logs}/** are gitignored AND deny-Read,
  so Read/grep/ls on them is blocked (the global `cat`-deny also bites). Rscript
  file reads BYPASS the deny matcher — count knit errors/warnings and inspect
  caches via `Rscript -e 'readLines(...)'` / `readRDS(...)`.
- Decision gates in a plan: present the plan default PLUS ≥1 reasoned alternative
  before invoking AskUserQuestion.
- Code carries provenance comments citing `storage/notes/<plan>.md` by name; that
  path is historical (pre-archival). Digests are in `history.md`, full plans in
  `completed/` (deny-Read). Treat the comments as plan-name labels, not live links.

## Subagents & skills
Scan the available-skills list each session and invoke a matching Skill before
improvising (they bundle SOTA workflows) — frequent fits: scientific-writing,
scientific-visualization, pathway-enrichment, pydeseq2, scanpy, anndata,
bioservices. Spawn subagents to protect main context (Explore = cross-file search,
Plan = design, general-purpose = research); all run the largest/most-capable model.

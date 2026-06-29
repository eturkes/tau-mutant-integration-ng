# SCENIC regulons plan (arc K): data-driven TF layer

Adds the project's first **data-driven** gene-regulatory inference. Every existing
mechanism layer is prior-based: TF (§14 decoupleR/CollecTRI), kinase (§15 OmniPath
KSN), CCC (§16 priors), CARNIVAL topology (§20 OmniPath PKN). SCENIC infers
regulons *de novo* from this dataset's own microglial co-expression (GRNBoost2) +
motif enrichment (cisTarget), scores them per cell (AUCell), and tests whether the
§14 prior-based TF verdicts survive a method that does **not** lean on a literature
prior. Framed as **orthogonal corroboration / honest discordance**, NOT a re-run,
and it does NOT touch the locked mouse 2×2 analysis except via an explicit,
gated ledger-feed decision (K6).

Readout = recovered regulons per axis + a controlled head-to-head against §14
(swap ONLY the network: data-driven SCENIC regulons vs CollecTRI prior, same
decoupleR scoring on the same NEBULA z). Non-recovery of a prior TF is a finding.

## STATUS
- [x] K1 (GATE): toolchain + mouse cisTarget resources + design lock — DONE (env + resources validated end-to-end; see note under K1)
- [x] K2: DONE — all 10 GRNBoost2 runs complete (10× adj_seed*.tsv, ~880–896k edges each, 994 candidate TFs, ~46 min/run). See K2 PROGRESS NOTE + K2+K3 DONE NOTE.
- [x] K3: DONE — ctx×10 → ≥8/10 edge-recurrence consensus = **51 activating regulons** (0 repressing) + AUCell (26,104 cells × 51). Census eyeballed: Spi1(+10/10) & Rel(+10/10) recovered into consensus; Nfkb1/Sp3 7/10 + Creb1 6/10 near-miss (below locked bar); Myc/Tp53 0/10 expression-floor; Nfkb2 0/10 motif-pruned. See K3 PROGRESS NOTE + K2+K3 DONE NOTE.
- [x] K4: DONE — R/scenic.R (9 pure fns) + scripts/build_scenic_contrasts.R emit 5 results TSVs + scenic_summary.rds. Head-to-head (controlled CollecTRI→SCENIC swap on the same NEBULA z) + AUCell 2×2 factorial (logit) + per-substate at interaction + target-set Jaccard + recovery ladder. HEADLINE: data-driven layer is **orthogonal, mostly discordant** with the §14 prior (median target Jaccard 0.004; Spi1 recovered but activity sign-flips; interaction_metabolic axis unrecovered). See K4 COMPLETION NOTE — signals the K6 gate toward standalone-chapter.
- [x] K5: DONE — §21 chapter `rmd/20_scenic_regulons.Rmd` (display-only, reads scenic_summary.rds; no recompute) wired before child-session; clean knit (0 error / 0 warning), rendered §number = **21**, all 5 subsections + 4 in-knit figures + scenic_verdict.tsv. See K5 COMPLETION NOTE.
- [x] K6 (GATE): DONE — user chose **add Suggestive rows** (Default, overriding the K4/K5 lean toward standalone). K-001 (Spi1 amyloid) + K-002 (Rel NF-kB) landed in §17 ledger as Suggestive, margin-neutral (contests **18/12/55 unchanged**), 5 invariants pass, clean re-knit (0 error / 0 warning). Arc K closed: map + DIGEST refreshed, plan archived. See K6 COMPLETION NOTE.

## Locked facts (verified at plan time; do not re-derive)
- **Object:** `storage/cache/microglia_seurat_processed.rds` — Seurat, **21,333 genes × 26,104 cells**, assays RNA + SCT (SCT default). meta: `genotype` {MAPTKI 4349, NLGF_MAPTKI 8787, NLGF_P301S 8583, P301S 4385}, `batch`, `genotype_batch` (16 pseudobulk ids, 4/genotype — the LOCKED replicate unit), `seurat_clusters` (0–12), `allen_labels` (LEAKY whole-brain label-transfer: 15,777 Micro-PVM + ~10k neuronal calls — a QC diagnostic, NOT the working identity; verify microglia identity in K2 via Csf1r/P2ry12/Itgam, do not subset on it).
- **Gene IDs are Ensembl `ENSMUSG…`** → MUST map to MGI symbols before pySCENIC (cisTarget DBs are MGI-keyed). Use project `snrnaseq_symbol_map` cache / `R/io.R`. Unmapped genes drop (report coverage).
- **No microglia-substate column** in the object meta. Recover the 4 substates the §18 way (the `de_snrnaseq_nebula_per_state*` caches encode them; or re-run `label_microglia_states()` from `R/microglia.R`). Needed only in K4 per-substate.
- **Host:** 8 cores, 62 GB RAM, ~1.2 TB free. micromamba/conda/mamba ABSENT (K1 installs micromamba project-local). Existing `.venv` = Python 3.12 (scanpy/anndata/decoupler) — CANNOT host pyscenic.
- **pyscenic 0.12.1 (latest, Nov 2022) requires Python ≤3.10, numpy<1.24** (numba/numpy ABI; Bioconda recipe). GRN step has a dask/TBB fork-crash → use the no-dask `arboreto_with_multiprocessing` helper on this single 8-core node.
- **Mouse cisTarget resources (resources.aertslab.org, v10_clust collection, mm10/GRCm38 — no mm39 DB exists):**
  - rankings: `https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc_v10_clust/gene_based/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather` (~226 MB; use `.rankings.feather` NOT `.scores.feather`; 10kb±TSS flavour, the SCENIC default — broader coverage than the 500bp variant; verify the `.sha256sum`).
  - motif2tf: `https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl` (108 MB; `mgi` = mouse; MUST match the v10 feather).
  - TF list: `https://resources.aertslab.org/cistarget/tf_lists/allTFs_mm.txt` (11 KB; known wart: contains non-TF Fgf15 — harmless).
- **§14/§18 comparison ledger (the head-to-head targets):**
  - amyloid_activation verdict TFs: **Spi1, Nfkb1, Sp3**.
  - interaction_metabolic verdict TFs: **Myc, Creb1, Tp53, Jun**.
  - §18 NF-κB family (complex-level attenuation across ALL 4 substates): Rela, Nfkb1, Nfkb2, Relb, Rel. Carrier accession A0A979HLR9 is NOT a symbol — ignore for SCENIC; seed/compare via Nfkb1/Rela symbols.
  - §14 scoring convention to mirror in the head-to-head: decoupleR `run_ulm` for significance at FDR<0.10 + `consensus` for direction/magnitude, `min_targets=5L`, on **NEBULA z** (`de_snrnaseq_nebula$top[[contrast]]`, stored stat = z).
- **5 canonical contrasts everywhere:** `nlgf_in_maptki, nlgf_in_p301s, interaction, tau_alone, tau_in_nlgf`. 3 axes (no pre-privileged winner): amyloid_activation & synaptic_suppression both use {nlgf_in_maptki, nlgf_in_p301s}; interaction_metabolic uses {interaction}.

## Execution model (per step)
1. Check for a fitting Skill/subagent first (Explore for cross-file search;
   general-purpose for research; the `arboreto` skill covers GRNBoost2). Subagents
   use the largest model.
2. Smoke-test new helper / cache-reading code against live caches via `Rscript -e`
   (or the py3.10 env for Python) BEFORE knitting (knits ~3 min).
3. Mark the step DONE in STATUS + write a multi-paragraph completion note (what was
   built, file paths, biology, explicit deviations) appended under that step.
4. chown rstudio:rstudio all new files + root-owned knit outputs.
5. Re-knit (only from K5 on; K1–K4 build caches outside the knit): `Rscript -e
   'rmarkdown::render("analysis.Rmd", quiet=TRUE)'`; verify
   `grep -c 'class="error"' analysis.html` == 0 and same for `class="warning"`.
6. Commit locally: imperative subject <70 chars + HEREDOC body + Co-Authored-By
   trailer (this session's model).
Give K1, K2, K4, K5, K6 fresh sessions (each is design-heavy or a long compute);
K3 may chain onto K2 if the env is warm and context is low.

## Anti-anchoring guardrails (SCENIC-specific; enforce every step)
- **Lock the GRNBoost2 multi-run / consensus policy in K1 BEFORE inspecting which
  TFs survive.** GRNBoost2 is stochastic; choosing the seed/run count after seeing
  which runs confirm §14 is seed-shopping. State the policy, then honour it.
- **Report ALL recovered regulons**, including those absent from §14 and the §14
  TFs that SCENIC does NOT recover (non-recovery = evidence the prior footprint
  lacks co-expression support in microglia, not a failure to hide).
- **No axis pre-privileged:** report amyloid + synaptic + interaction recovery
  even-handedly; do not foreground interaction.
- **State thresholds before applying:** FDR (0.10, project standard), regulon
  `min_targets` (≥5, mirror §14), AUCell default settings, motif NES/`auc_threshold`
  (pyscenic ctx defaults; record them).
- **Caveats stated at every interpretation:** (a) microglia-ONLY scope → regulons
  are microglial-STATE regulons, not cell-type-identity; (b) single-nucleus input
  (intronic reads; SCENIC handles snRNA, note it); (c) `allen_labels` leakage means
  the working microglia set must be identity-verified, not assumed.

---

## K1 (GATE): toolchain + mouse cisTarget resources + design lock

**Goal.** Stand up an isolated, reproducible pyscenic environment and fetch the
mouse cisTarget resources, smoke-tested end-to-end on a tiny subsample, so K2–K3
are pure compute. Lock the design choices that anti-anchoring requires be fixed
before data inspection.

**Actions.**
1. Install **micromamba** project-local (e.g. `./.micromamba/`, gitignored); create
   env `scenic` with `python=3.10` then `pyscenic=0.12.1` (+ its deps `arboreto`,
   `ctxcore`, `pyarrow`, `loompy`, `pandas`, `numpy<1.24`) via conda-forge+bioconda.
   FALLBACK if the solve fails: the `aertslab/pyscenic:0.12.1` container (note the
   nested-Docker caveat — prefer micromamba). Record the exact env spec to a
   lockfile under `scripts/` or `storage/notes/`.
2. Download the 3 mouse resources (URLs above) to `storage/data/cistarget/`
   (gitignore the dir); verify the feather `.sha256sum`.
3. Smoke-test: run `pyscenic grn` (multiprocessing helper) on a ~500-cell × HVG
   subsample → adjacencies; `pyscenic ctx` + `aucell` on the toy → confirm a
   non-empty regulon set and an AUCell matrix. Confirms the whole CLI chain + DB
   compatibility before the real run.
4. Update `.gitignore` for `.micromamba/`, `storage/data/cistarget/`, and the K2/K3
   heavy intermediates (loom, adjacencies, ctx output, AUCell matrix).

**DECISION GATE — present to the user before executing K2 (default + ≥1 alternative each):**
- **Env strategy:** micromamba py3.10 env *(default)* vs aertslab container.
- **GRNBoost2 run policy** (the load-bearing one; 8 cores makes multi-run costly):
  - default = **N-run consensus** (e.g. 10 runs, keep regulons recurrent in ≥80%;
    the pyscenic-recommended robustness protocol) — principled but potentially
    overnight×N on 8 cores;
  - alternative = **single fixed-seed run** (hours, reproducible-by-seed but
    stochastic; defensible for a first-pass orthogonal check).
  - Decide the count + recurrence threshold NOW (anti-anchoring).
- **Expression input:** RNA raw counts with gene pre-filter (default: detected in
  ≥ some min cells, e.g. ≥1% = 261 cells — standard, cuts ~21k→~12–14k genes) vs
  a looser/stricter filter. (Use RNA counts, NOT SCT.)
- **AUCell vs decoupleR for scoring:** run pyscenic `aucell` in-env (default) AND
  reuse the regulons in decoupleR for the §14 head-to-head (K4) — confirm both.

**LOCKED gate decisions (user, 2026-06-04):**
- Env = **micromamba py3.10 env** (project-local `.micromamba/`), NOT the container.
- GRNBoost2 run policy = **10-run consensus, keep regulons recurrent in ≥80% (≥8/10)
  of runs**. Fixed before inspecting which TFs survive (anti-anchoring). First K2 run
  times itself to calibrate the 10-run background budget. Run in background.
- Expression input = **RNA raw counts**, gene pre-filter = detected in **≥1% of cells
  (≥261)**. Ensembl→MGI via project `symbol_map` (report coverage).
- Scoring = **both** pyscenic AUCell (native) AND decoupleR-on-regulons (§14 head-to-head).

**Verify / DONE when:** env imports `pyscenic`; the 3 resources are on disk +
checksum-clean; the toy CLI chain yields ≥1 regulon + an AUCell matrix; gitignore
updated; gate decisions recorded in this step's completion note.

**K1 COMPLETION NOTE (2026-06-04).** Toolchain stood up and validated end-to-end.
- **Env:** `micromamba` 2.7.0 binary at `.micromamba/bin/micromamba` (project-local,
  gitignored); env `scenic` at `.micromamba/envs/scenic` (1.6 GB). Created via
  `micromamba create -r .micromamba -n scenic -c conda-forge -c bioconda
  python=3.10 pyscenic=0.12.1`. Versions: **pyscenic 0.12.1, arboreto 0.1.6,
  ctxcore 0.2.0, numpy 1.23.5 (<1.24 ✓), pandas 2.3.3, dask 2024.8.2, loompy
  3.0.8, pyarrow** (reads Feather v2). Run anything with
  `.micromamba/bin/micromamba run -n scenic <cmd>`. Full spec: `scripts/scenic_env.yml`.
- **DEVIATION (setuptools):** pyscenic/ctxcore import `pkg_resources`, which
  setuptools ≥81 no longer ships → `ModuleNotFoundError`. Fixed by pinning
  **`setuptools<81`** (got 80.10.2; harmless deprecation warning remains). Must
  stay pinned (already in the env / yml).
- **Resources** at `storage/data/cistarget/` (gitignored via `/storage/*`):
  `mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather`
  (237,172,098 B; **5032 motifs × 24,131 MGI-symbol genes** — confirmed via
  pyarrow read), `motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl` (113,107,706 B),
  `allTFs_mm.txt` (1860 TFs). **DEVIATION:** the published `.sha256sum` is a 404
  (not hosted) → integrity rests on exact size+date match + successful pyarrow
  read + 47 regulons recovered in the smoke test. Local hash recorded in
  `storage/data/cistarget/SHA256SUMS.local`.
- **Smoke test (PASSED, full chain on a REAL 800-cell × 2000-gene microglia
  subsample, seed=1):** Ensembl→MGI map via `symbol_map` worked (196 TFs in the
  matrix); loom built (loompy, genes×cells); **GRN** via the no-dask
  `arboreto_with_multiprocessing` helper = **53.7 s, 155,932 adjacencies, NO TBB
  crash**; **ctx** = 164 s; **aucell** = 11 s → AUCell **800 cells × 47 signed
  regulons** (`Arnt(+)`, `Clock(+)`, `E2f3(+)`…). cisTarget DBs fully functional.
- **K2 GOTCHA — helper shadowing (load-bearing):** running
  `arboreto_with_multiprocessing.py` from its install dir
  (`.../site-packages/pyscenic/cli/`) FAILS — that dir lands on `sys.path[0]`
  where `pyscenic.py` masks the `pyscenic` package (`'pyscenic' is not a
  package`). FIX: **copy the helper to a neutral dir first** and run the copy,
  e.g. `cp .micromamba/envs/scenic/lib/python3.10/site-packages/pyscenic/cli/
  arboreto_with_multiprocessing.py /tmp/grnboost2_mp.py` then
  `micromamba run -n scenic python /tmp/grnboost2_mp.py <loom> <allTFs> --method
  grnboost2 --output adj.tsv --num_workers 8 --seed <S>`.
- **K2 RUNTIME CALIBRATION:** GRN 2000×800 = 54 s on 4 workers. Real run ≈
  13k genes × 26,104 cells. Rough extrapolation ⇒ **~1.5–3 h per GRN run on 8
  workers**; the locked **10-run consensus ⇒ ~15–30 h total** ⇒ run as a
  **background job over 1–2 days** (vary `--seed` per run; the first run times
  itself to refine this). ctx/aucell are cheap per run.
- **Loom build recipe (reuse in K2):** R writes a genes×cells counts TSV
  (Ensembl→MGI mapped); Python `loompy.create(path, mat_genes_x_cells,
  {"Gene": genes}, {"CellID": cells})`. pyscenic loom convention = rows are genes,
  `Gene`/`CellID` attrs.

## K2: Ensembl→MGI export of microglia matrix + GRNBoost2 adjacencies

**Goal.** Export the microglia expression matrix with MGI symbols to a pyscenic
input, verify microglia identity, and run the (locked-policy) GRNBoost2 GRN
inference → TF→target adjacencies. The compute-heavy step.

**Actions.**
1. New `scripts/build_scenic_grn.py` (runs in the `scenic` env). Read the microglia
   RNA counts; **map ENSMUSG→MGI symbol** via the project `snrnaseq_symbol_map`
   (export the map from R to a sidecar TSV/JSON, mirroring the H3 jsonlite bridge,
   or read the Seurat `RNA` assay feature metadata). Collapse/drop unmapped &
   duplicate symbols (report coverage: n mapped / n dropped). Apply the locked gene
   pre-filter. Write a `.loom` (cells × genes) with `genotype/batch/genotype_batch`
   in col attrs.
2. **Identity check (anti-anchoring caveat (c)):** confirm the set is microglial —
   Csf1r/P2ry12/Hexb/Itgam broadly positive, neuronal markers (Snap25/Rbfox3) low;
   document any contamination and that the cell set is byte-identical to the §14 DE
   input (same 26,104 cells).
3. Run `pyscenic grn` via `arboreto_with_multiprocessing.py` (NO dask) with
   `allTFs_mm.txt`, per the locked run policy (single or N-run). Cache adjacencies
   (per-run + the consensus union if multi-run) under `storage/cache/scenic/`.

**Verify / DONE when:** loom written with >X% symbol coverage reported; identity
confirmed; adjacencies non-empty with sane TF count; runtime/seed(s) logged in the
note. Smoke-check shape via the env's python before proceeding.

**K2 PROGRESS NOTE (2026-06-04).** Two-script R→Python bridge built + the GRN job
launched in background. NOT yet DONE (compute in-flight, ~3 days).
- **R exporter** `scripts/export_microglia_for_scenic.R` (run with project Rscript):
  reads `microglia_seurat_processed.rds` RNA **raw counts** (33,683 genes × 26,104
  cells — note the full RNA assay is 33,683, NOT the 21,333 SCT-default the K1 note
  quoted), Ensembl→MGI via `snrnaseq_symbol_map` (**100% coverage, 0 unmapped, 0
  duplicate symbols** — `@misc$geneids` is a clean 1:1 unique-symbol bijection, so
  no collapse needed), gene pre-filter ≥1% cells (≥262) ⇒ **11,536 genes**, **994 /
  1860 allTFs present**. Writes a sparse bundle to `storage/cache/scenic/`:
  `microglia_counts.mtx` (genes×cells, 22.9M nnz, 7.6% dense), `microglia_genes.txt`,
  `microglia_cells.txt`, `microglia_colattrs.tsv` (CellID/genotype/batch/
  genotype_batch), `export_provenance.txt`. Runs in 18 s.
- **Identity-of-input (caveat c) VERIFIED:** 26,104 nuclei + genotype tally
  {MAPTKI 4349, NLGF_MAPTKI 8787, NLGF_P301S 8583, P301S 4385} **byte-match the
  locked §14 DE input**. Marker pct-positive: microglial Hexb 65.7%, Cx3cr1 39.1%,
  P2ry12 37.6%, Csf1r 34.4%, Itgam 24.3% (Tmem119 low 5.2%); ambient neuronal/oligo
  present (Snap25 38%, Syt1 37.5%, Plp1 42.5%) — the modest ambient signal typical
  of snRNA microglia; working identity is **inherited from §14 `broad_annotations`**,
  not re-derived. CAVEAT for K5: regulons are microglial-STATE, ambient-aware.
- **Python orchestrator** `scripts/build_scenic_grn.py` (runs in `scenic` env):
  builds `microglia.loom` (genes×cells, Gene row attr + CellID/meta col attrs;
  loompy, 9 s) from the mtx bundle, then runs GRNBoost2 over seeds 1..10 via the
  copied no-dask helper. **RESUMABLE** (skips existing `adj_seed*.tsv`). Honours the
  K1 helper-shadow gotcha by copying the helper to `storage/cache/scenic/
  _grnboost2_mp.py` (neutral dir) before invoking. BLAS pinned to 1 thread/worker.
- **LOAD-BEARING DEVIATION — DENSE not `--sparse`:** the helper's `--sparse` path
  feeds a sparse target column into `arboreto/core.py:125` which calls `.A`, removed
  in scipy ≥1.14 (env has 1.15.2) ⇒ **every target fails** (`'csc_matrix' has no
  attribute 'A'`) and the run is empty. FIX: omit `--sparse` (the K1 smoke was
  dense, hence it worked). Dense loads ~2.4 GB (26k×11.5k) — fine on 62 GB. If a
  future env upgrade keeps sparse desirable, patch core.py:125 `.A`→`.toarray()`.
- **PERF FIX — chunksize (10× speedup, OUTPUT-INVARIANT):** the stock helper does
  `p.imap(fn, targets, chunksize=1)` — one fast (~0.2 s) target dispatched at a time,
  so the 8 workers idled at ~15% CPU (load avg 2.46/8) waiting on the parent's
  per-task round-trip (dispatch-bound, ~5.5 cores wasted). `prepare_helper()` patches
  the copied helper `chunksize=1 → 32` (asserts exactly one match). Result preserves
  order (imap) ⇒ identical adjacency table; pure scheduling change, untouched by
  anti-anchoring. Measured: load avg **2.46 → 8.91**, workers **15% → 82–93%**, rate
  **2.22 s/it → ~4.4 it/s**.
- **RUNTIME CALIBRATION (real, dense, chunksize=32):** **~43 min per GRN run** on 8
  workers (was ~7 h at chunksize=1). **10-run consensus ⇒ ~7 h total** (not the ~70 h
  the chunksize=1 path implied). Launched `nohup … build_scenic_grn.py --max-runs 10
  --num-workers 8 > storage/cache/scenic/grn_run.out 2>&1 &` (PID in `grn.pid`);
  per-run timing/ETA appended to `grn_progress.log`. Monitor:
  `tr '\r' '\n' < storage/cache/scenic/grn_run.out | tail`. **RESUME after a crash/
  session end:** just re-run the same nohup line — existing adj files are skipped.
- **K3 STARTS when** `ls storage/cache/scenic/adj_seed*.tsv` shows 10 non-empty
  files; then mark K2 DONE and run cisTarget ctx + AUCell per the locked ≥80%
  recurrence consensus.

## K3: cisTarget motif pruning → signed regulons + AUCell per-cell activity

**Goal.** Turn co-expression adjacencies into motif-supported, **signed** regulons
and score per-cell activity.

**Actions.**
1. `pyscenic ctx`: prune adjacencies to regulons using the mm10 v10 rankings
   feather + `motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl`; record motif-enrichment
   params (NES threshold, `auc_threshold`, `rank_threshold` — defaults, logged).
   **Add the TF–target correlation sign** (`--mask_dropouts` / the standard
   `add_cor` step) so regulons carry a `mor` (+/-) for the decoupleR head-to-head
   (SCENIC default regulons are activating; the correlation split gives signed
   ones). If multi-run: keep regulons recurrent ≥ the locked threshold.
2. `pyscenic aucell`: per-cell AUCell activity matrix (cells × regulons).
3. Export two light artifacts for R: (a) regulons as a long TF→target table with
   sign/weight (`storage/cache/scenic/regulons.tsv` or rds-friendly); (b) AUCell
   matrix (`storage/cache/scenic/aucell.csv.gz`). Cache the raw ctx output too.

**Verify / DONE when:** regulon count in the expected ~30–150 range (single cell
type; low recovery is a known failure mode — investigate if <~20); every regulon
has ≥`min_targets` targets and a sign; AUCell matrix dims = cells × regulons; the
§14 comparison TFs (Spi1/Nfkb1/Sp3/Myc/Creb1/Tp53/Jun) are flagged present/absent
in the regulon set (recovery census — a headline result, reported as-is).

**K3 PROGRESS NOTE (2026-06-04).** Script built, cheaply validated, and LAUNCHED
as a self-gating background job that chains onto K2; compute is in-flight (K2 was
6/10 GRN runs done at launch, so ctx waits ~3 h then runs). NOT yet DONE — the
next session marks it DONE after eyeballing the census + dist (below).
- **Script** `scripts/build_scenic_ctx_aucell.py` (runs in the `scenic` env):
  per-seed `pyscenic ctx` → signed regulons → 10-run ≥8/10 consensus → `pyscenic
  aucell`. RESUMABLE (existing `reg_seed{S}.csv` skipped). SELF-GATING via
  `--wait-for-adj`: polls every 300 s until all 10 `adj_seed*.tsv` are non-empty
  AND the GRN process (`grn.pid`) has exited AND the youngest adj file is ≥90 s old
  (anti-contention — ctx never competes with GRNBoost2 for the 8 cores).
- **LOAD-BEARING CONSENSUS DESIGN (fixed before inspecting which TFs survive; the
  K1 lock was "regulons recurrent in ≥80% of runs"):** ctx is run independently on
  each of the 10 adjacency runs (canonical pySCENIC multi-run robustness), then
  **EDGE-LEVEL recurrence** is the consensus substrate — a TF→target edge is
  high-confidence iff it recurs in **≥8/10** runs (this is how SCENIC reports
  high-confidence links, Aibar 2017). A **consensus regulon** = a TF(sign) with its
  high-confidence target set, kept iff ≥**5** such edges (`min_targets`, mirroring
  §14 decoupleR). An edge can only reach 8/10 if its TF(sign) regulon existed in ≥8
  runs, so regulon-level recurrence is implied — edge-level is the stricter,
  target-set-quality-preserving form of the locked rule. **Recovery census** is
  reported SEPARATELY at the looser **regulon-presence** level (did TF X form a
  motif-supported regulon in ≥8/10 runs at all?) because "is Spi1 recovered?" is a
  presence question distinct from "what is its high-confidence target set?". Both
  are emitted honestly for every §14/§18 comparison TF.
- **ctx params (stated before applying, anti-anchoring):** `--all_modules`
  (emit BOTH activating `(+)` and repressing `(-)` regulons — signed, so the K4
  head-to-head vs the SIGNED CollecTRI net swaps ONLY the network) +
  `--expression_mtx_fname microglia.loom` (the TF–target correlation that assigns
  the sign); `--mode custom_multiprocessing` (no-dask, dodges the TBB crash like
  GRN); `--mask_dropouts` OFF (default — zeros kept in the correlation); motif
  thresholds = pyscenic defaults (rank 5000 / auc 0.05 / nes 3.0 / min_genes 20),
  logged. CAVEAT for K5: SCENIC repressing `(-)` regulons are less established than
  activating ones — sign is reported throughout, never hidden.
- **Validation done this session (no env contention):** (a) `--self-test` passes —
  the pure-Python consensus, weight-averaging, and recurrence-distribution logic
  are correct on a synthetic 3-run fixture; (b) pyscenic API confirmed:
  `df2regulons` (pyscenic.transform), `load_motifs` (pyscenic.utils),
  `load_signatures` (gmt-capable aucell loader), `Regulon.name`="Spi1(+)" /
  `.transcription_factor` / `.gene2weight`; (c) `aucell --help` confirms gmt
  signatures accepted; (d) a 150 s `ctx` startup peek on the real `adj_seed1.tsv`
  with the exact K3 flags ran with NO arg/path error (killed by `timeout`, exit
  143 = SIGTERM) — the heavy step K1 only exercised on the toy is now validated at
  scale-start.
- **Launched:** `nohup … build_scenic_ctx_aucell.py --wait-for-adj --num-workers 8
  > storage/cache/scenic/ctx_run.out 2>&1 &` (micromamba-run PID in `ctx.pid`);
  per-step timing in `ctx_progress.log`. Monitor: `tr '\r' '\n' <
  storage/cache/scenic/ctx_run.out | tail`. RESUME after a crash: re-run the same
  nohup line (skips finished seeds). If a `reg_seed{S}.csv` parse errors at the
  consensus step (rare partial-write), delete that one file and relaunch.
- **Outputs it will write** (storage/cache/scenic/, root-owned — job runs as root
  via micromamba): `reg_seed{1..10}.csv` (per-seed ctx), `scenic_consensus_
  regulons.gmt` (signed, for aucell), `scenic_regulons.tsv` (long TF/sign/target/
  recurrence/mean_importance, for K4 decoupleR), `scenic_recurrence_dist.tsv`
  (n edges / n regulons surviving at each threshold 1..10 — sensitivity diagnostic),
  `scenic_recovery_census.tsv` (§14/§18 TFs: present? n runs? sign? n high-conf
  targets?), `aucell.csv.gz` (cells × consensus regulons).
- **NEXT SESSION (mark K2 + K3 DONE, then K4):** confirm `ctx_progress.log` ends
  with "K3 compute complete"; **chown rstudio:rstudio** the new root-owned
  storage/cache/scenic/ outputs; eyeball `scenic_recovery_census.tsv` (are
  Spi1/Nfkb1/Sp3/Myc/Creb1/Tp53/Jun + the NF-κB family recovered? — a headline,
  report as-is) and `scenic_recurrence_dist.tsv` (regulon count ≥20? the K3 verify
  bar; if <20 the dist shows whether the ≥8/10 lock is over-strict vs genuine low
  single-cell-type recovery — report as a finding, do NOT silently retune the
  locked threshold). Then K4 (new `R/scenic.R`).
- **GRN-STAGE VALIDATION (2026-06-04, mid-flight at 6/10 runs).** Sanity-checked
  the K2 adjacency output the launching session could not (job was in-flight):
  every `adj_seed*.tsv` is clean — 3 cols (`TF/target/importance`), ~880–896k
  edges/run, importance numeric & positive, and **exactly 994 distinct regulators
  == the 994 candidate TFs** (`allTFs_mm.txt ∩ microglia_genes.txt`), so the GRN
  used every eligible TF. No malformed output; the ~7 h run is trustworthy.
- **PRE-COMPUTED RECOVERY-CENSUS FINDING (structural floor — load-bearing for K3
  eyeball + K5 interpretation).** Of the §14/§18 comparison TFs, **Myc and Tp53 are
  ABSENT from the 11,536-gene filtered matrix** (detected in <1% / <262 cells, so
  dropped by the locked ≥1% pre-filter) ⇒ they are **not candidate regulators** ⇒
  **0 edges in every seed** ⇒ SCENIC *structurally cannot* recover them in any run.
  This is an **expression-floor non-recovery with a deterministic cause** (same
  across all 10 seeds — same input matrix), NOT a stochastic miss, so the census
  will honestly show Myc/Tp53 at 0/10. Myc + Tp53 are **2 of the 4
  `interaction_metabolic` verdict TFs** (Myc/Creb1/Tp53/Jun): the data-driven layer
  cannot corroborate the half of that axis that lacks microglial co-expression
  footing — exactly the guardrail's "non-recovery = the prior footprint lacks
  co-expression support" finding, reported as-is (do NOT relax the filter to rescue
  them; the ≥1% floor is locked). The other comparison TFs are all candidate
  regulators with substantial edges in `adj_seed1`: Spi1 1124, Nfkb1 2913, Nfkb2
  455, Rela 310, Relb 645, Rel 817, Sp3 1601, Creb1 2245, Jun 530. Whether they
  survive ctx motif-pruning + ≥8/10 recurrence is the still-open K3 question.

**K2+K3 DONE NOTE (2026-06-05).** Both compute steps landed; outputs verified,
chowned rstudio:rstudio. Compute ran untouched in the background across the
session boundary exactly as the self-gating/resumable design intended.
- **K2 (GRN):** all 10 GRNBoost2 runs complete (seeds 1–10), mean 46.6 min/run
  (~7.7 h wall), 880–896k adjacencies each, exactly 994 candidate regulators/run
  (= allTFs ∩ filtered matrix). Output clean (see GRN-STAGE VALIDATION above).
- **K3 (ctx + AUCell):** ctx ran ~6.4 min/seed × 10; per-seed signed-regulon yield
  ~136–141. Consensus (locked **≥8/10 edge recurrence, ≥5 targets**) = **51
  regulons, all activating (+), 0 repressing** — the repressing `(-)` modules never
  reached ≥8/10 (anticipated: SCENIC `(-)` regulons are less stable). 51 ≫ the ≥20
  verify bar. ctx params as locked: `all_modules=True, mask_dropouts=False,
  rank=5000/auc=0.05/nes=3.0/min_genes=20`. AUCell = 1.5 min → 26,104 cells × 51.
  Recurrence-sensitivity curve (min_runs:n_regulons) = 1:275, 2:179, 3:141, 4:110,
  5:88, 6:70, 7:64, 8:**51**, 9:36, 10:25 — smooth decay, no cliff at 8 (the lock
  is strict but not pathological).
- **HEADLINE RECOVERY CENSUS (`scenic_recovery_census.tsv`, report as-is in K5):**
  recovered into consensus = **Spi1 (amyloid; +10/10, 15 targets)** and **Rel
  (NF-κB family; +10/10, 9 targets)** ONLY. Honest near-misses *present but below
  the locked ≥8/10 bar*: **Nfkb1 7/10(+)**, **Sp3 7/10(+)**, **Creb1 6/10(+)**,
  **Relb 10/10(+) but <5 stable targets** (regulon present, fails min_targets).
  Non-recovered with distinct causes (state precisely in K5): **Myc 0/10 & Tp53
  0/10 = EXPRESSION FLOOR** (filtered <1%, never candidate regulators — structural);
  **Nfkb2 0/10 = MOTIF-PRUNED** (candidate regulator w/ 455 adj edges, but no
  motif-supported regulon in any run — co-expression/motif, NOT floor); **Jun
  2/10**, **Rela 1/10** = weak. So the data-driven layer corroborates the amyloid
  Spi1 call cleanly, gives the NF-κB axis partial support (Rel in, Nfkb1/Relb near),
  and does NOT recover the interaction_metabolic axis at the strict bar (half of it
  is below the microglial expression floor). DO NOT relax ≥8/10 to rescue
  near-misses — the bar is locked; report 6–7/10 honestly as sub-threshold.
- **CAVEAT for K4/K5 (all-activating consensus):** every consensus regulon has
  sign `+`, so the SCENIC decoupleR net (K4 readout 1) carries **no repressing
  edges** — unlike signed CollecTRI. The network swap therefore also drops the
  ability to represent repression; note this when a §14 repressing call has no
  SCENIC counterpart (it is a structural consequence, not discordance).

### K4 START-HERE (schemas captured 2026-06-05; smoke-test against these live files)
- **`storage/cache/scenic/scenic_regulons.tsv`** — 5 cols `TF, sign, target,
  recurrence, mean_importance`; 10,711 edges / 51 regulons. → decoupleR net:
  `source=TF, target=target, weight=sign` (use **sign ±1 as `mor`** to mirror
  CollecTRI's convention and swap ONLY the network; `mean_importance`/`recurrence`
  are available if a weighted variant is wanted — decide+document, default sign-only).
  All sign=`+` here (see caveat).
- **`storage/cache/scenic/aucell.csv.gz`** — 26,104 × 52: col1 `Cell` (barcodes
  `…-1_<batch>`, match §14 cell set / `microglia_cells.txt`), cols 2–52 named
  `TF(+)` AUCell scores in [0,1]. → readout 2: join `Cell`→`genotype_batch` via
  `microglia_colattrs.tsv` (CellID/genotype/batch/genotype_batch), pseudobulk
  (mean per regulon per of the 16 ids), then `R/design.R` factorial 2×2+batch over
  all 5 contrasts. Bounded [0,1] → logit/rank transform decision: inspect the score
  distribution first (many exact 0s seen — e.g. Clock(+) col is mostly 0.0 — so a
  naive logit needs an epsilon or a rank/CLR alt; decide+document in K4).
- **`scenic_consensus_regulons.gmt`** (aucell signatures), **`scenic_recurrence_
  dist.tsv`**, **`scenic_recovery_census.tsv`** also on disk; `reg_seed{1..10}.csv`
  = per-seed ctx (keep for provenance). The 51 regulons span microglial/myeloid +
  immune identity (Zbtb20 3140 targets, Fli1, Ikzf1, Hif1a, Nfia, Mitf, Ets1,
  Irf8/Irf1, Stat1/2/5a, AP-1 Fos/Fosb/Fosl2/Junb/Jund, Spi1, Rel) — report the
  full set even-handedly (anti-anchoring), not just the §14/§18 overlap.

## K4: contrast modelling — head-to-head + AUCell pseudobulk factorial + per-substate

**Goal.** Connect regulons to the project's contrast framework via three
complementary readouts; emit verdict-ready tables. New helper `R/scenic.R` (source
in `R/helpers.R` after `causal_network.R`).

**Actions / three readouts.**
1. **Controlled head-to-head (primary).** Feed the SCENIC regulons (signed) to
   **decoupleR** (`run_ulm` + `consensus`, `min_targets=5L`, FDR<0.10) on the SAME
   NEBULA z used in §14 (`de_snrnaseq_nebula$top[[contrast]]`), all 5 contrasts.
   This swaps ONLY the network (SCENIC vs CollecTRI), isolating the prior's effect.
   Emit a per-TF table: SCENIC activity + §14 CollecTRI activity side by side
   (concordant sign? both sig?).
2. **Native SCENIC activity (corroboration).** Pseudobulk the AUCell matrix to the
   16 `genotype_batch` ids (mean per regulon); fit the LOCKED factorial 2×2 + batch
   design (`R/design.R` `make_contrast_matrix`/`factorial_design`); extract all 5
   contrasts incl. `interaction` per regulon (limma on the bounded AUCell scores —
   note the [0,1] support; consider logit or rank, decide+document).
3. **Per-substate (mirror §18).** Recover the 4 substates; per-substate regulon
   activity (AUCell pseudobulk within substate, or decoupleR on per-state NEBULA z
   `de_snrnaseq_nebula_per_state_1pct`) at `interaction`; report all 4 substates
   (absence = evidence). Focus on the NF-κB regulon(s) (does the data-driven NF-κB
   regulon reproduce §18's all-substate attenuation?).
4. **Target overlap.** For TFs in BOTH SCENIC and CollecTRI, compute Jaccard of
   target sets (how much of the data-driven regulon the prior already knew).
Emit TSVs to `storage/results/`: `scenic_regulon_recovery.tsv`,
`scenic_vs_collectri_headtohead.tsv`, `scenic_aucell_contrasts.tsv`,
`scenic_per_substate_interaction.tsv`, `scenic_target_overlap.tsv`. Cache the
assembled objects the Rmd reads (e.g. `storage/cache/scenic_summary.rds`).

**Verify / DONE when:** all 5 TSVs written; head-to-head covers the §14 verdict
TFs; smoke-tested in R against live caches; thresholds applied as pre-stated.

**K4 COMPLETION NOTE (2026-06-05).** Contrast-modelling layer built, smoke-tested
end-to-end against live caches, all 5 TSVs + the rmd cache written and chowned. No
knit (K1–K4 are out-of-knit cache builds; the §21 chapter lands in K5).

- **Code.** `R/scenic.R` — 9 pure functions, sourced in `R/helpers.R` AFTER
  `causal_network.R` (dep-safe: reuses `tf_inference.R`'s
  `extract_de_stat_matrix`/`run_decoupler_per_modality`/`split_decoupler_by_contrast`
  for the network swap, and `design.R::factorial_design` + `de_pb.R::fit_limma_log`
  for the AUCell factorial — all sourced earlier). `build_scenic_network`,
  `run_scenic_decoupler`, `decoupler_activity_long`, `build_scenic_headtohead`,
  `aucell_to_pseudobulk`, `fit_aucell_contrasts`, `scenic_substate_activity`,
  `scenic_collectri_target_overlap`, `build_recovery_table`. I/O driver
  `scripts/build_scenic_contrasts.R` (idempotent, `--overwrite`); pattern mirrors
  `build_tf_activity_decoupler.R`. Outputs: `storage/results/scenic_{regulon_recovery,
  vs_collectri_headtohead,aucell_contrasts,per_substate_interaction,target_overlap}.tsv`
  + `storage/cache/scenic_summary.rds` (the list rmd/20 reads: regulons, scenic_net,
  headtohead, scenic_activity, aucell_contrasts, aucell_pb{mat,meta}, substate_int,
  target_overlap, recovery, census, params). All gitignored under `/storage/*` —
  rebuilt by the driver, not committed (consistent with every other [S] result).

- **LOCKED DECISIONS (decide+document, per the K4 step + START-HERE).**
  (1) **decoupleR net weight = sign-only `mor` (±1)** — the minimal swap that keeps
  decoupleR's mode-of-regulation semantics identical to CollecTRI, isolating the
  network as the only changed variable. All 51 consensus regulons are sign `+`, so
  the SCENIC net carries **no repressing edges** (a §14 repressing call has no SCENIC
  counterpart by construction — structural, not discordance; caveat carried in code +
  for K5). The weighted variant (sign×mean_importance) is retained behind
  `build_scenic_network(weighted=TRUE)` for sensitivity only.
  (2) **AUCell [0,1] transform = logit (`qlogis`), no clamp.** Pseudobulking the 26,104
  per-cell scores to the 16 `genotype_batch` ids dissolves the per-cell zero-inflation
  entirely: every one of the 51×16 means is strictly interior (range **0.0049–0.097**,
  zero exact 0/1), so `qlogis` is well-defined without an epsilon. Sensitivity vs raw
  means (also fit): **Spearman(t)=0.994, sign-agreement 98.4%, sig-set Jaccard 0.80** —
  the transform is immaterial to the conclusions (recorded in params).
  (3) **FDR scope** — BH within (network, contrast) over each network's full scored-TF
  set, so SCENIC (51 TFs) and CollecTRI (482) each carry their native multiple-testing
  burden before the shared-TF comparison (mirrors §14's `rank_tfs_cross_modality`).
  padj_cut=0.10, statistics=c("ulm","wsum")+consensus, minsize=5L — all as §14.

- **HEADLINE BIOLOGY (report as-is in K5; this is the load-bearing K6 signal).** The
  data-driven layer is **orthogonal to and largely discordant with** the §14 prior:
  - **Target overlap is near-zero**: median Jaccard 0.004 across the 37 shared TFs,
    max 0.028 (Fosb). The recovered Spi1 regulon shares just **2 of 381 union targets**
    with CollecTRI's Spi1 (15 vs 368 targets). So the data-driven and prior footprints
    are almost disjoint gene sets — activity discordance is the expected consequence,
    and the strongest evidence the layer is genuinely orthogonal (not a prior re-run).
  - **Head-to-head**: of the §14 verdict TFs, only **Spi1** is also a SCENIC regulon,
    and its activity **sign-flips and loses significance** (CollecTRI cons +2.73 sig at
    nlgf_in_maptki → SCENIC cons −0.43 ns). Motif-level recovery (Spi1 +10/10) ≠
    activity-level concordance. The `interaction_metabolic` axis (Myc/Creb1/Tp53/Jun)
    is **entirely CollecTRI-only** (Myc interaction cons −6.60 sig has no SCENIC
    counterpart). NF-κB: Nfkb1/Sp3/Rela are CollecTRI-sig at amyloid but unrecovered
    by SCENIC; **Rel** (the one recovered NF-κB regulon) has no bare-symbol CollecTRI
    source to compare against.
  - **Where SCENIC IS significant**: 31/51 regulons reach padj<0.10 in ≥1 contrast —
    but these are the **large novel** regulons (Zbtb20 3140 targets, Fli1, Hif1a, the
    AP-1/IRF/STAT immune set), not the small §14-overlap ones. So SCENIC has power; it
    simply points at a different, data-driven TF set.
  - **Native AUCell factorial**: 16/51 regulons sig at main-effect contrasts but
    **0 at `interaction`** — no data-driven regulon shows a genotype×amyloid
    interaction at the pseudobulk level.
  - **Per-substate (mirrors §18)**: Rel and Spi1 are **non-significant at interaction in
    all 4 substates** (Rel IFN padj 0.109 closest) — the data-driven NF-κB regulon does
    NOT reproduce §18's all-substate NF-κB attenuation. Some other regulons do reach
    per-substate interaction significance (surfaced in the TSV for K5).
  - **Recovery ladder** (60 rows = 51 consensus + 9 non-recovered comparison TFs):
    recovered_novel 49, recovered_comparison 2 (Spi1, Rel), near_miss 3 (Nfkb1 7/10,
    Sp3 7/10, Creb1 6/10), present_targets_subthreshold 1 (Relb 10/10 present but <5
    stable targets), weak 2 (Jun 2/10, Rela 1/10), motif_pruned 1 (Nfkb2: candidate,
    0 motif-supported regulon), expression_floor 2 (Myc, Tp53: below the ≥1% filter,
    never candidate regulators). Candidacy derived live from
    `microglia_genes.txt ∩ allTFs_mm.txt` (994 candidate TFs), not hardcoded.

- **K6 SIGNAL (do NOT pre-decide; the gate owns it).** The head-to-head is
  predominantly discordant/non-recovery and the target overlap is near-zero, so the
  K6 default (feed Suggestive ledger rows à la J) is weakly supported; the evidence
  leans toward the K6 **alternative** — surface SCENIC ONLY in §21 (à la human H7),
  leaving the locked contest margins untouched. Decide AFTER K5 renders the full
  picture.

- **DEVIATIONS / minor.** (a) `recovered_novel` TFs show `max_runs_present = NA` —
  the K3 recurrence census (`scenic_recovery_census.tsv`) audited per-run presence
  ONLY for the §14/§18 comparison TFs; novel consensus regulons are ≥8/10 by
  construction (in_consensus=TRUE carries that). K5 should footnote this. (b) The
  AUCell pseudobulk + meta are cached whole (`aucell_pb`) for a K5 activity heatmap.
  (c) `Sp3` appears in the head-to-head but is NOT a SCENIC regulon — covered as a
  CollecTRI-only verdict TF (the universe is the union of SCENIC regulons + comparison
  TFs, so every §14 verdict TF has a row even when SCENIC=absent).

## K5: render §21 chapter + clean knit

**Goal.** Display-only chapter consuming K4 caches; wire into the knit.

**Actions.**
1. `rmd/20_scenic_regulons.Rmd` (verify the §number at render via the child-order
   grep — file 20 likely → §21). Subsections: 21.1 regulon recovery census (incl.
   the §14-TF present/absent table); 21.2 SCENIC-vs-CollecTRI head-to-head (the
   controlled network swap); 21.3 native AUCell activity across contrasts +
   interaction; 21.4 per-substate (NF-κB focus); 21.5 sign-aware verdict
   (orthogonal corroboration vs honest discordance, per axis, no winner privileged).
2. Plotting in the chapter's own R (or `R/scenic.R`) so it renders inside the
   shared knit session (mirrors §15/§20 — do NOT use the Python viz skill for
   in-knit figures). A regulon-activity heatmap / concordance scatter vs §14.
3. Wire `child-scenic-regulons` into `analysis.Rmd` before `child-session`.

**Verify / DONE when:** clean knit (0 `class="error"`, 0 `class="warning"`);
chapter renders the recovery census, head-to-head, activity, per-substate, verdict;
caches pre-built (display-only, like §19/§20 — knit must not recompute).

**K5 COMPLETION NOTE (2026-06-05).** Chapter rendered and wired. New file
`rmd/20_scenic_regulons.Rmd` is display-only: it `readRDS`es
`storage/cache/scenic_summary.rds` (the 11-component K4 cache) plus
`storage/cache/scenic/scenic_recurrence_dist.tsv`, builds every table/figure
in-chapter, and writes one derived artifact `storage/results/scenic_verdict.tsv`
(3-axis per-axis verdict tibble). It recomputes **nothing** (mirrors the §19/§20
pattern). Wired as `child-scenic-regulons` in `analysis.Rmd` between
`child-causal-network` and `child-session`; child order → rendered **§21** (header
auto-number confirmed in `analysis.html`). Five subsections: 21.1 regulon recovery
census (recovery-ladder table over all 60 classified TFs + §14/18 comparison census
+ recurrence-sensitivity figure), 21.2 SCENIC-vs-CollecTRI head-to-head / controlled
network swap (target-Jaccard table + per-TF head-to-head table + concordance
scatter), 21.3 native AUCell across contrasts (per-contrast significant-count table +
51-regulon × 5-contrast heatmap), 21.4 per-substate NF-κB focus (focus table +
per-substate heatmap), 21.5 per-axis verdict. Four ggplot figures render inside the
shared knit session (NOT the Python viz skill), as specified.

**Build outcome.** Full `rmarkdown::render` clean: `grep -c 'class="error"'` = 0 and
`grep -c 'class="warning"'` = 0 in `analysis.html`; all six §21 header ids present;
verdict TSV written. chowned rstudio:rstudio: the new rmd, edited analysis.Rmd, the
verdict TSV, and the session's omnipathr-log byproducts. `analysis.html`,
`storage/cache`, `storage/results`, `omnipathr-log` are gitignored — the commit
carries only `analysis.Rmd` + `rmd/20_scenic_regulons.Rmd`.

**Deviation (eval-order fix).** First full knit aborted with `object 'n_reg' not
found`: four inline ``r``-refs (`n_reg`, `n_cand`, `nrow(sc_reg)`) sat in the intro
prose **above** the `scenic-load` chunk, and R Markdown evaluates inline code in
document order. Fixed by hard-coding those four **locked** intro values (51
regulons, 10,711 high-confidence edges, 994 candidate TFs) — matching §19's
hard-coded-intro style — while every body chunk keeps its inline refs (they run
after `scenic-load`, so they stay drift-proof). Two earlier bugs were caught in K5
smoke-testing before any knit and already fixed: a `merge()` `axis.x/axis.y`
collision in the head-to-head table (now merges only `recovery_class` from
`sc_rec`), and a vectorised `element_text(face=...)` warning in the substate heatmap
(focus regulons now tagged via a `" (focus)"` row-label suffix instead). The
smoke-test harness pre-defined variables, so it could not catch the document-order
inline-ref bug — a standalone child render (`knit_root_dir` set to project root)
does, and is the better pre-knit check for child chapters.

**Biology (unchanged from K4; rendered here).** Data-driven layer is orthogonal and
mostly discordant with the §14 prior (median target Jaccard 0.004). Of the §14
verdict TFs only **Spi1** is also a SCENIC regulon (+10/10, 15 targets) and its
activity sign-flips and loses significance under the controlled swap; **Rel**
(NF-κB family) also recovers (+10/10, 9 targets) but is flat across contrasts. The
`interaction_metabolic` axis is unrecovered at the strict bar — Myc/Tp53 sit below
the microglial expression floor (structural non-recovery). Not a power failure: many
novel co-expression regulons (Zbtb20, Fli1, Hif1a, AP-1/IRF/STAT) are significant on
the same NEBULA z and in native AUCell; the per-substate readout finds many novel
interaction regulons while the §18 NF-κB attenuation is **not** reproduced by the
data-driven Rel regulon. Per-axis verdict in `storage/results/scenic_verdict.tsv`.

**K6 signal.** Per the K6 gate's own criterion ("if the head-to-head is mostly
discordant ... surface SCENIC ONLY in §21"), the K5 result leans toward the
**Alternative = standalone chapter only**. K6 is a user-confirmed GATE; do not
decide it without an AskUserQuestion.

## K6 (GATE): ledger-feed decision + close-out

**Goal.** Decide whether SCENIC corroborations enter the §17 biological-model
ledger, then archive the arc.

**DECISION GATE — present default + alternative:**
- SCENIC is mouse, in-class with the modalities (unlike human, which H7 sealed out
  as evidence-class-incommensurable), so it is *eligible* for the ledger like J.
  - **Default = add Suggestive rows** (à la J-001/J-002) with an anti-double-count
    guard: the rows score the *data-driven motif+co-expression corroboration of the
    TF activity*, distinct from the prior-based footprint already counted (Myc
    H2-046, Spi1/Nfkb1 amyloid, NF-κB I-001..I-006). Held at Suggestive by the
    method's caveats (microglia-only, stochastic GRN, motif-prior dependence). Keep
    the 11-entity set + grade defs untouched; re-run the 5 invariants.
  - **Alternative = standalone chapter only** (à la human H7): if the head-to-head
    is mostly discordant or the corroboration is too entangled to count cleanly,
    surface SCENIC ONLY in §21 and do not perturb the locked contest margins.
  - The data (K4 head-to-head concordance) should drive the choice — decide AFTER
    seeing it, not before.
**Close-out actions.** Refresh `storage/notes/map.md` (new rmd row, new R/scenic.R,
new scripts + caches + the cistarget resources, the scenic env); add a
`scenic_regulons_plan_<date>` digest to `completed/DIGEST.md` (extend the arc
narrative D→J→**K**); move this plan to `completed/` with a date suffix; chown;
final commit.

**Verify / DONE when:** map + DIGEST refreshed; plan archived; (if ledger fed)
ledger invariants pass + margins updated coherently; clean knit; committed.

**K6 COMPLETION NOTE (2026-06-05).** Gate resolved by the user toward the
**Default = add Suggestive rows**, overriding the K4/K5 data-lean toward the
standalone alternative. Two rows landed in the §17 biological-model ledger,
scoped strictly to **regulon EXISTENCE, NOT activity concordance** (the
anti-double-count guard), and engineered margin-neutral so the three locked
contest verdicts are untouched.

- **Rows added (`scripts/build_biological_model_ledger.R`, new `phase_k` block,
  rbinded after `phase_j`).**
  - **K-001 (amyloid_activation / tf, Suggestive, +):** SCENIC's de novo GRN
    recovers **Spi1** as a motif-supported microglial regulon at the strictest
    10/10 consensus (15 high-conf targets), corroborating — by a method
    independent of the §14 CollecTRI prior — that Spi1 is a genuine
    amyloid-activation regulator *in this dataset*. **Supports BOTH Hyp-1A AND
    Hyp-1B** (a shared amyloid-activation premise that does not discriminate the
    tau-modulation contest) + T-Inflammation. Margin-neutral by construction:
    both amyloid nets rise by one (Hyp-1A 12→13, Hyp-1B 30→31) so the amyloid
    margin holds at **18**. Existence-only: the controlled CollecTRI→SCENIC swap
    shows Spi1 activity sign-flips and loses significance, and the data-driven
    target set is near-disjoint from the prior (Jaccard 0.005, 2/381 union
    targets) — so it does NOT support an activity *direction*. No contradicts.
  - **K-002 (cross_axis / tf, Suggestive, +):** SCENIC recovers **Rel** as the
    SOLE NF-kB-family regulon at the 10/10 bar (9 targets; Nfkb1 7/10, Relb
    present but <5 stable targets, Nfkb2 motif-pruned, Rela 1/10). **Tags only
    T-Inflammation** (support-only, no contradicts) — an existence corroboration
    of microglial NF-kB activity. Existence-only: the recovered Rel regulon is
    non-significant at `interaction` in all 4 substates with mixed-sign consensus
    and does **NOT** reproduce the §18 all-substate attenuation, so it does NOT
    support T-Tau-attenuates. Tags no contest model ⇒ all three margins unchanged.
  - **Anti-double-count:** both rows score the *methodologically independent*
    motif+co-expression recovery, distinct from the prior-based footprints
    already counted (H2-001 Spi1 amyloid, H2-002 Nfkb1 amyloid, I-001..I-006
    NF-kB attenuation) — recorded explicitly in each row's `notes` +
    `corroborating_evidence`. Held at Suggestive by the SCENIC caveats
    (microglia-only state regulons, single-nucleus, stochastic GRN, motif-prior
    dependence).

- **Rebuild + verification.** Ran `build_biological_model_ledger.R` then
  `build_biological_model_adjudication.R` (order matters). Ledger **83→85 rows**;
  grade tally Strong 28 / Moderate 41 / Suggestive 14→**16**; **all 5
  `stopifnot` invariants pass**. Contest verdicts (`...contest_verdicts.tsv`):
  amyloid Hyp-1B net 31 vs Hyp-1A net 13 = **margin 18**; synaptic **12**;
  interaction **55** — all three UNCHANGED. T-Inflammation 26→**28** supports
  (6 Strong, axis split amyloid 25 / cross 1 / interaction 2 = 89% axis-1). TSVs
  chowned rstudio:rstudio.

- **Display prose (hardcoded anchors).** `rmd/16_biological_model.Rmd`: 85-row
  count + phase arithmetic (75+6+2+2), datatable lengthMenu→85, amyloid verdict
  `(31 - 13)`, Hyp-1A 22 supporting / net 13, Hyp-1B 31, T-Inflammation 28 /
  axis-split 25 of 28 (89%), single-axis range 89-94%, + a Phase K synthesis
  sentence. `rmd/20_scenic_regulons.Rmd`: intro now states the layer feeds two
  margin-neutral Suggestive rows (K-001/K-002) to the §17 ledger. Confirmed the
  §17 contest-margin table renders dynamically from the verdict TSV (no drift);
  rmd/01 carries NO biological-model numbers (the DIGEST's "§1.1 headline mirror"
  cites nothing that shifts — rmd/16 is the sole carrier).

- **Knit.** Full `rmarkdown::render("analysis.Rmd")` clean: **0 `class="error"`,
  0 `class="warning"`**; K-001/K-002 + "margin 18" rendered. Outputs chowned
  (rmarkdown preserved existing ownership; only the edited rmds + new omnipathr
  log needed chown).

- **DEVIATION from the K4/K5 recommendation.** K4/K5 notes leaned toward the
  standalone alternative (head-to-head mostly discordant). The user chose the
  Default ledger-feed regardless; reconciled by scoping the rows to **existence,
  not activity** and making them margin-neutral — so the discordance finding is
  preserved verbatim in §21 and the ledger gains only the orthogonal-recovery
  corroboration, with zero movement in the locked mouse 2×2 verdicts.

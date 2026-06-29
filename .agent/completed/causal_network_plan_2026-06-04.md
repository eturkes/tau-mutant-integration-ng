# Causal signalling network reconstruction (CARNIVAL / COSMOS) — plan

Narrative-arc letter **J** (D pathway → E/F/G mechanism → H biological model →
I nfkb → human → **J causal topology**).

**Goal.** The mechanism layer (E/F/G) produced *endpoint activities* — kinases
(Gsk3b strongest, axis 3), TFs (Spi1/Nfkb1/Sp3 axis 1; Myc/Creb1/Tp53/Jun axis
3), LR pairs (Apoe_Trem2/App_Cd74 axis 2) — but **no causal wiring between
them**. This plan reconstructs, per contrast, the directed signed signalling
subnetwork (OmniPath prior-knowledge network, solved by ILP via CARNIVAL or
cosmosR) that connects the **kinase/phospho upstream** layer to the
**TF/transcription downstream** layer. Headline questions: (a) is the
interaction-axis Gsk3b signal topologically wired to the interaction-axis TFs?
(b) does the tau-driven NF-κB attenuation (§18) have a recoverable upstream
path? Readout = recovered topology per axis, **honestly reporting what connects
and what does not** — an ILP that omits an expected edge is a finding, not a bug.

**Scope lock.** Reuses caches already on disk (no new DE / no new activity
inference). Mirrors the mechanism-layer house pattern: one helper
`R/causal_network.R`, one heavy build script `scripts/build_causal_network.R`
(compute OUTSIDE the knit), one display-only chapter `rmd/19_causal_network.Rmd`
rendering as **§20**, X.1 ranking / X.2 axis-restricted / X.3 sign-aware verdict.

---

## STATUS (read this first; execute only the next TODO step)

- [x] **J1**  Environment + design gate — install CARNIVAL/cosmosR + a free ILP
      solver, verify a toy `runCARNIVAL` end-to-end; **DECIDE** (gate) tool/mode,
      TF-complex symbol handling, PKN expression-filter, ledger-feed. *(Decision
      gate at start; design-heavy → own session.)* **— DONE 2026-06-04; decisions locked below.**
- [x] **J2**  `R/causal_network.R` — OmniPath signed-directed PKN assembly +
      CARNIVAL input formatters (TF `measObj`, kinase `inputObj`) + runner
      wrapper; source() into R/helpers.R; smoke-test vs live caches.
      **— DONE 2026-06-04; smoke test 0 fail. KEY FINDING: split_complexes is
      a no-op on mouse CollecTRI (split cache == non-split byte-for-byte); see
      note. CARNIVAL PKN col is `interaction` (not `sign`).**
- [x] **J3**  `scripts/build_causal_network.R` — run CARNIVAL per contrast (5);
      emit `storage/cache/causal_network.rds`. **— DONE 2026-06-04. KEY: full
      microglia-pruned PKN is ILP-INTRACTABLE → user-approved per-contrast L=3
      reachability prune (`restrict_pkn_to_reachable`); solver cbc @ threads=1
      (threads>1 segfaults cbc 2.10.3). 3 nets solved (p301s 33e/18-of-22TF,
      interaction 8e/5-of-5, tau_in_nlgf 27e/13-of-13), 2 honestly empty
      (maptki/tau_alone: 0 kinase seeds). See note.**
- [x] **J4**  `rmd/19_causal_network.Rmd` (§20) + analysis.Rmd child wiring;
      X.1 per-contrast networks / X.2 axis-restricted / X.3 verdict; emit
      `results/causal_network_{nodes,edges,verdict}.tsv`; knit-verify 0 error/warning.
      **— DONE 2026-06-04. Knit clean (0 error/0 warning, 4m44s). Headline (a)
      Gsk3b⊣Mapk14→Myc/Foxo3 recovered @interaction; (b) NF-κB path mixed
      (Gsk3b⊣Ikbkg→Nfkbib @tau_in_nlgf; Gsk3b→Nfkb1 @nlgf_in_p301s; none
      @interaction). New helper `plot_causal_network` (ggraph/stress). See note.**
- [x] **J5**  *(J1 LOCKED: ADD rows)* append Suggestive-graded topology row(s) +
      re-run schema invariants + adjudication; re-knit §17 + §20.
      **— DONE 2026-06-04. 2 Suggestive rows added (J-001 interaction bridge,
      contradicts Hyp-3A; J-002 NF-κB route, support-only). 83-row ledger clean,
      all 5 invariants pass, interaction margin 53→55, knit 0 error/0 warning.
      Also fixed a pre-existing §17 mis-citation (H2-067→H2-063 gap row). See note.**
- [x] **J6**  Close-out: update map.md + completed/DIGEST.md, final knit-verify,
      commit; move this plan to completed/ with date suffix.
      **— DONE 2026-06-04. map.md (§20 pipeline row + anchor 19→§20 + causal_network.R
      helper + causal lane in scripts + 2 cache rows) & DIGEST (J-arc entry + arc
      header) updated; chowned build_causal_network.R (was root). J5 knit stands
      (0 error/0 warning) — J6 touches only notes/, no knit inputs. Committed;
      plan archived to completed/. See note.**

---

## J1 gate decisions + de-risk (locked 2026-06-04 — J2/J3/J5 obey these)

**Toolchain (installed + verified):** CARNIVAL 2.20.0, cosmosR 1.18.1, dorothea,
lpSolve. Toy `runVanillaCarnival(solver=lpSolve)` ran end-to-end (0.26 s,
non-empty `weightedSIF`); PKN format confirmed `(source, sign∈{+1,−1}, target)`.
lpSolve is the in-R free solver (lpsymphony present but NOT CARNIVAL-native).
**J3 housekeeping:** `defaultLpSolveCarnivalOptions()` writes `lpFile`/
`parsedData` to `workdir` (defaults to **project root**) with `keepLPFiles=TRUE`
→ the build script MUST set `workdir`/`outputFolder` to a tempdir and
`keepLPFiles=FALSE` to avoid polluting the tree.

**Gotcha 2 RESOLVED — mouse PKN suffices, NO human-ortholog fallback.** Mouse
`OmnipathR::omnipath_interactions(organism=10090)` = 30,501 raw →
**18,126 signed-directed symbol edges / 5,718 nodes** after
`is_directed & xor(is_stimulation, is_inhibition)`, in `source_genesymbol`/
`target_genesymbol` space (sign = +1 stim / −1 inhib). Endpoint coverage: all
verdict kinases present except Prkca; all verdict TFs present except Stat1.
Directed reachability is rich: Gsk3b→Nfkb1/Rela/Spi1 @1 hop, →Myc/Creb1/Tp53/Jun
@2 hops; every probed kinase reaches every probed TF in 1–3 hops. The
Gsk3b→Nfkb1/Rela 1-hop edge pre-supports the §18 NF-κB-attenuation path question.

**The 4 locked decisions:**
1. **Tool/mode = standard per-contrast CARNIVAL** — kinases as signed `inputObj`
   → TF activities as signed `measObj` over the OmniPath mouse PKN. (Not inverse,
   not cosmosR.) 5 contrasts.
2. **TF complex handling = re-derive `split_complexes=TRUE`** *(NOT the plan
   default (a)).* The on-disk `tf_activity_decoupler.rds` is built
   `split_complexes=FALSE` (`build_tf_activity_decoupler.R:84`); the causal layer
   needs a `split_complexes=TRUE` re-derivation so CollecTRI complexes (e.g.
   NF-κB) become single-protein symbol nodes mapping cleanly to the PKN. Produce
   `tf_activity_decoupler_split.rds` — simplest is a `--split-complexes` flag on
   `build_tf_activity_decoupler.R` toggling line 84 (all else identical);
   `tf_meas_from_cache` reads that variant (PKN mapping then identity).
3. **PKN filter = prune to microglia-expressed** (snRNAseq microglia subset);
   caveat in X.3 that bulk-phospho kinases seed a microglia-filtered network.
4. **Ledger feed = ADD Suggestive-graded topology rows** *(NOT the no-op
   default).* J5 active: append Suggestive claim(s) for a genuinely novel,
   independently-gradeable recovered path; re-run the 5 schema invariants +
   adjudication (margins shift only as intended; never reopen the 11-entity set
   or grade defs); re-knit §17 + §20.

---

## Grounded facts (verified 2026-06-04; do NOT re-discover)

### Inputs (caches on disk — the two endpoints)
Both are nested lists `modality -> contrast -> tibble`, every leaf an identical
4-col `tbl_df`: **`statistic`** (chr: `ulm`/`wsum`/`norm_wsum`/`corr_wsum`/`consensus`
— stacked, so each node×contrast appears ~5×; FILTER to one), **`source`** (chr:
node id), **`score`** (num: signed activity — **the sign IS the direction**),
**`p_value`** (num: raw p; **no padj stored** — recompute BH within
modality×contrast, threshold **0.10** per project convention). **No sign column,
no padj column.** Convention: `consensus` for sign/magnitude, `ulm` for
significance.

- `storage/cache/tf_activity_decoupler.rds` — 7 modalities: `snrnaseq`,
  `snrnaseq_{homeostatic,DAM,IFN,proliferative}`, `geomx`, `proteomics`.
  `source` space = **UniProt accessions** (CollecTRI `split_complexes=FALSE`,
  e.g. `A0A087WPA7`) **intermixed** with symbol-named TFs (Spi1, Nfkb1, Myc…).
- `storage/cache/kinase_activity_decoupler.rds` — 2 modalities: `phospho_raw`,
  `phospho_corrected`. `source` space = **mouse gene symbols** (Akt1, Gsk3b…).
- **5 canonical contrasts (identical names in both caches):** `nlgf_in_maptki`,
  `nlgf_in_p301s`, `interaction`, `tau_alone`, `tau_in_nlgf`.
- **3 axes (locked E4 regex on `contrasts_summary`):**
  `amyloid_activation = "nlgf_in_[a-z0-9]+:[0-9]+/[+]"`,
  `synaptic_suppression = "nlgf_in_[a-z0-9]+:[0-9]+/-"`,
  `interaction_metabolic = "interaction:[0-9]+/mixed"`. Amyloid+synaptic both use
  {nlgf_in_maptki, nlgf_in_p301s}; interaction uses {interaction}.

### Tooling state (Bioconductor 3.22)
- CARNIVAL **available, not installed**; cosmosR **available, not installed**.
- Installed: `decoupleR` 2.16.0, `OmnipathR` 3.18.4, `MOFA2`, `CellChat`,
  **`lpsymphony`** (the only ILP-capable pkg present). **Absent:** `lpSolve`,
  `lpSolveAPI`, `Rcplex`; **no** cbc/glpsol/cplex binary on PATH; apt has **no**
  coinor-cbc/glpk package. → free-solver path = install `lpSolve` (CRAN,
  self-contained, CARNIVAL-native, slow on big ILPs) FIRST to de-risk; obtain a
  cbc binary (conda-forge `coincbc`, or static COIN-OR download) only if lpSolve
  is too slow on the restricted networks. CLAUDE.md grants free install latitude;
  keep installs project-local where possible.

### OmniPath PKN
`OmnipathR::import_omnipath_interactions(organism = 10090)` → **30,501 mouse
interactions** (deprecation note: newer alias `omnipath_interactions()`).
Direction/sign cols: `is_directed`, `is_stimulation`, `is_inhibition`,
`consensus_direction` (TRUE=18,935), `consensus_stimulation`,
`consensus_inhibition`. Node ids: UniProt `source`/`target` **and**
`source_genesymbol`/`target_genesymbol`. CARNIVAL PKN format =
`(source, sign∈{+1,-1}, target)`: keep directed + (stimulation XOR inhibition);
work in **mouse gene-symbol** space (both endpoints map there once TF complexes
are resolved — see gotcha 1).

### File / section conventions (match house style exactly)
- **Build script pattern** (`scripts/build_*.R`): `args <- commandArgs();
  overwrite <- "--overwrite" %in% args; source("R/helpers.R"); if
  (file.exists(out) && !overwrite) quit(status=0); …compute…;
  saveRDS(obj, out, compress="xz"); Sys.chmod(out, "0644")`. Idempotent
  `file.exists` guard (NOT `cache_or_run`). Compute lives here, not in the knit.
- **Rmd pattern:** one level-1 `# ` header (→ the §); subsections `## X.1/X.2/X.3`,
  sub `### `; tables via **`knitr::kable`** in a `results='asis'` chunk (NOT
  DT::datatable — that is reserved to §17.2); read cache via
  `readRDS(file.path(params$cache_dir, "causal_network.rds"))`; emit TSVs via
  `readr::write_tsv(x, file.path(params$results_dir, "..."))`.
- **New Rmd = `rmd/19_causal_network.Rmd` → renders as §20.** Wire as
  `child-causal-network` in analysis.Rmd **between line 179 (`child-human-validation`)
  and line 182 (`child-session`)**. (File# = §#−1 because rmd 08/09 deleted;
  rmd/18→§19, so rmd/19→§20.)
- **New helper** must be `source()`d in `R/helpers.R` in dependency order **after**
  `tf_inference.R` (it reuses `build_axis_gene_universe()` and the axis rules).
- Reusable from `R/tf_inference.R`: `build_axis_gene_universe(leaderboard,
  gene_sets, axis_rules)`, the 3 `axis_rules`, the 5-contrast vector (defined
  inline in the Rmd, not a constant).

### Gotchas that will bite
1. **TF symbol-space split.** Kinase `source` = mouse symbols; TF `source` =
   UniProt accessions (CollecTRI complexes) mixed with symbols. CARNIVAL nodes
   must be single proteins in the PKN. The verdict TFs (Spi1, Nfkb1, Sp3, Myc,
   Creb1, Tp53, Jun) ARE symbol-named singles; the NF-κB **complex**
   (`A0A979HLR9`, the §18 carrier) is the exception. **J1 gate decision:** either
   (a) map symbol-named TFs to OmniPath symbol nodes and drop/decompose
   complex-accession sources (documented), or (b) re-derive TF activity with
   `split_complexes=TRUE` for this layer only (clean singles, but diverges from
   §14's prior). Default (a). **→ J1 LOCKED (b): re-derive `split_complexes=TRUE`
   → `tf_activity_decoupler_split.rds` (build_tf_activity_decoupler.R:84 flag).**
2. **Mouse vs human PKN.** Mouse OmniPath = 30,501 edges (thinner than human). If
   coverage strands the endpoints, fallback = map mouse→human orthologs
   (`nichenetr`), solve on the 9606 PKN, map back. Try mouse first.
3. **Cross-compartment bridge (state in prose).** Kinase endpoints are **bulk
   hippocampus phospho, NOT microglia-sorted** (standing lock); TF endpoints are
   microglial snRNAseq (+geomx+proteomics). The recovered network bridges
   bulk-tissue signalling to microglial transcription — an explicit
   cross-compartment inference that MUST be caveated in X.3.
4. **PKN expression-filter (J1 sub-decision).** Default: prune PKN to genes
   expressed in the snRNAseq microglia subset (microglial readout drives the
   question); document that bulk-phospho kinases seed a microglia-filtered
   network. Alternative: leave PKN unfiltered (more intermediates, larger ILP).
5. **No solver binary** — see tooling state; start with `lpSolve`.

---

## Steps

### J1 — Environment + design gate *(Decision gate at start)*
**Goal:** stand up the solver toolchain, prove a CARNIVAL run end-to-end, and
lock the design forks before any project-specific code.

**Actions:**
1. Install (project lib): `BiocManager::install(c("CARNIVAL","cosmosR"))` +
   `install.packages("lpSolve")`. Verify `requireNamespace` for all three.
2. **Toy verification:** run `CARNIVAL`'s bundled example
   (`loadModelData()` / the package vignette toy: small `netObj`, `measObj`,
   `inputObj`) with `solver="lpSolve"` end-to-end; confirm it returns a weighted
   network without error. This de-risks the solver path before PKN assembly. If
   lpSolve errors or hangs, obtain cbc (conda-forge `coincbc` or static binary
   into a project `bin/`) and re-verify with `solver="cbc", solverPath=...`.
3. **Present the gate** (plan default + ≥1 reasoned alternative, then
   AskUserQuestion) on FOUR coupled decisions:
   - **Tool/mode.** *Default:* **standard per-contrast CARNIVAL**, kinases as
     signed `inputObj` (upstream perturbation) → TF activities as signed
     `measObj` (downstream), OmniPath PKN — the explicit phospho→TF bridge.
     *Alt 1:* **inverse CARNIVAL** (no `inputObj`; infer upstream from TFs alone,
     then test whether Gsk3b/Cdk5 emerge unprompted — stronger anti-anchoring,
     loses the explicit bridge). *Alt 2:* **cosmosR signalling→transcription**
     (wraps CARNIVAL with meta-PKN + expression preprocessing; heavier dep,
     no metabolomics so it reduces to a richer-formatted CARNIVAL). *Consider
     also:* run **both standard + inverse** and report concordance (does the
     data-nominated upstream set match the kinase-seeded set?).
   - **TF complex handling** — gotcha 1: default (a) map symbol TFs, drop/decompose
     complex accessions; alt (b) re-derive TF activity `split_complexes=TRUE`.
   - **PKN expression-filter** — gotcha 4: default prune to microglia-expressed;
     alt unfiltered.
   - **Ledger feed (J5)** — *default:* **stand-alone §20, NO ledger edit** —
     the causal net is a *derivation of already-counted* TF+kinase endpoints, so
     adding rows risks double-counting evidence; mirrors how E/F/G each stood as
     their own § before H synthesised. *Alt:* add **Suggestive**-graded topology
     rows only for a genuinely novel, independently-gradeable recovered path
     (re-runs the §17 invariant checks; never reopens the 11-entity set / grade
     defs). H7 sealed *human* rows as evidence-class-incommensurable; mouse
     topology is in-class, so eligible — but the double-counting caveat is the
     reason the default is no-edit.

**Done when:** toolchain installed, toy CARNIVAL verified, and the 4 gate
decisions are recorded (as a note appended here). Marks J2–J5 unblocked. **Own
session; stop at the gate for user confirmation.**

### J2 — `R/causal_network.R` helper + smoke test
**Goal:** project-specific PKN + CARNIVAL I/O, reusing tf_inference utilities.
**Prereq (J1 dec. 2):** add a `--split-complexes` flag to
`build_tf_activity_decoupler.R` (toggles line 84 → writes
`tf_activity_decoupler_split.rds`, `split_complexes=TRUE`, all else identical)
and run it; the smoke test + `tf_meas_from_cache` consume that split cache.
**Functions (sketch; finalise to J1 mode):**
- `build_omnipath_pkn(organism=10090, expr_filter=NULL, ...)` → CARNIVAL
  `(source, sign, target)` df: directed + signed (stim XOR inhib) → ±1, symbol
  space, optional microglia-expressed prune.
- `tf_meas_from_cache(tf_cache, contrast, modalities, fdr=0.10, ...)` → named
  signed vector (`measObj`): filter `statistic=="consensus"`, BH on `ulm`
  `p_value` within modality×contrast <0.10, read the J1-dec.-2 split cache
  (`tf_activity_decoupler_split.rds`; sources already single-protein symbols →
  PKN mapping is identity), aggregate across modalities (sign-consistent mean or
  per-modality).
- `kinase_input_from_cache(kin_cache, contrast, modality="phospho_corrected",
  fdr=0.10, ...)` → named signed vector (`inputObj`), per J1 mode (NULL if inverse).
- `run_carnival_for_contrast(pkn, measObj, inputObj=NULL, solver, ...)` → wrapper
  returning a tidy `list(nodes=<tbl: node, activity_sign, ...>, edges=<tbl:
  source, sign, target, weight>, meta=<solver/penalty/contrast>)`.
- `summarise_network(net)` + `restrict_network_to_axis(net, axis_universe)`
  (reuse `build_axis_gene_universe()` from tf_inference.R for X.2).
**Smoke test (before any knit — knits cost ~3 min):** `Rscript -e` loads both
live caches, builds PKN, formats `measObj`/`inputObj` for `interaction`, asserts
non-empty + correct ID space + sign ∈ {−1,+1}; do NOT run the full ILP here.
**Done when:** helper sourced in R/helpers.R (correct order), smoke test passes,
files chowned rstudio:rstudio.

### J3 — `scripts/build_causal_network.R` + cache
**Goal:** heavy ILP runs OUTSIDE the knit, emit one cache.
**Actions:** idempotent `file.exists`+`--overwrite` guard; loop the 5 contrasts
(and per axis if J1 chose axis-restricted) calling `run_carnival_for_contrast`;
assemble `out <- list(<contrast> = net, ..., pkn_meta=..., params=...)`;
`saveRDS(out, "storage/cache/causal_network.rds", compress="xz")`;
`Sys.chmod(0644)`. Record solver, penalty, FDR, PKN-filter, n nodes/edges in/out
per contrast for the manifest. Run it; sanity-check the cache shape via a small
`Rscript` read.
**Done when:** cache exists with all contrasts solved (or honestly-empty
networks recorded), shapes verified, chowned. *(Compute may be slow — own session.)*

### J4 — `rmd/19_causal_network.Rmd` (§20) + wiring + knit
**Goal:** display-only chapter, X.1/X.2/X.3, mirroring §14.
- `# Causal signalling network reconstruction` (→ §20). Intro: question, the two
  endpoints, the cross-compartment caveat (gotcha 3), thresholds STATED upfront.
- **20.1 Per-contrast reconstructed networks** — read cache; per contrast a
  summary kable (nodes in/out, edges, key recovered upstream→TF paths) + a focused
  subnetwork figure (igraph/ggraph node-link, signed edges coloured). Emit
  `causal_network_{nodes,edges}.tsv`.
- **20.2 Axis-restricted networks** — restrict each network to the 3 axis
  universes (`build_axis_gene_universe`); per-axis tables/figures. Interaction
  axis headline: is Gsk3b→…→{Myc/Creb1/Jun} present?
- **20.3 Verdict (sign-aware topology synthesis)** — per axis: what connects /
  what is absent (absence = evidence); NF-κB-attenuation path check; restate
  phospho-bulk caveat; no pre-privileged axis. Emit `causal_network_verdict.tsv`.
- Wire `child-causal-network` into analysis.Rmd before `child-session`.
**Done when:** `rmarkdown::render("analysis.Rmd", quiet=TRUE)` →
`grep -c 'class="error"' analysis.html` == 0 AND `class="warning"` == 0; outputs
+ new caches chowned.

### J5 — Ledger feed *(J1 LOCKED: ADD Suggestive rows)*
J1 chose **add rows** (not the no-op default). Append Suggestive-graded topology
claim_id(s) for a genuinely novel, independently-gradeable recovered path via the
existing scripts (`build_biological_model_ledger.R` →
`build_biological_model_adjudication.R`); re-run the 5 schema invariants +
adjudication, verify contest margins shift only as intended, re-knit §17 + §20.
**Anti-double-count guard:** grade **Suggestive only**, and only for a path that
is independently gradeable beyond the already-counted TF/kinase endpoints; never
reopen the 11-entity set or grade definitions. **Done when:** ledger loads clean,
invariants pass, §17 + §20 re-knit with 0 error/warning.

### J6 — Close-out
Update `map.md` (add R/causal_network.R to helpers; build_causal_network.R to
scripts; causal_network.rds to the cache table [K?/S → S]; rmd/19→§20 to the
pipeline + child tables) and `completed/DIGEST.md` (new J-entry: goal, outcome,
key decisions w/ rationale, artifacts, gotchas, verdict). Final knit-verify.
Commit. Move this file to `completed/causal_network_plan_<YYYY-MM-DD>.md`.

---

## Execution model
Per step: (0) check for a fitting Skill/subagent first — `statistical-analysis`,
`pathway-enrichment`, `scientific-visualization` (network figure), and spawn
Explore/general-purpose (largest model) for cross-file search/research to protect
main context. (1) smoke-test helper/cache-read code via `Rscript -e` before
knitting. (2) mark the step DONE here + multi-paragraph completion note (what
built, paths, biology, deviations). (3) chown rstudio:rstudio new + root-owned
knit outputs. (4) re-knit; verify 0 `class="error"`/`class="warning"`. (5) commit
locally (imperative subject <70 chars + HEREDOC body + Co-Authored-By trailer).
Give J1 (gate) and J2/J3 (new PKN/ILP shape) **fresh sessions**; J4→J6 may chain
if each is a trivial pattern-extension. End cleanly at the J1 gate or low context
(monitor `./compaction.sh`, stop ≥80%).

## Anti-anchoring guardrails
- **An omitted edge is a finding.** The verdict reports what the ILP connects AND
  what it fails to connect; never assert a path the solver did not return.
- **No pre-privileged axis/contrast** (mechanism-layer lock); present all 3 axes
  and 5 contrasts; STATE thresholds (FDR 0.10, solver, penalty, PKN filter)
  before running.
- **NO specificity-null / expression-matched-random / randomised-input control**
  as a default framework (standing ban). A null here must be motivated by a
  specific question, never as default re-anchoring. If network robustness is
  questioned, report sensitivity to solver penalty / FDR instead.
- **Cross-compartment honesty:** phospho is bulk hippocampus, not
  microglia-sorted — restate wherever kinase-seeded paths are interpreted.
- **No double-counting:** the causal net derives from already-counted TF+kinase
  endpoints; default is no ledger edit (J5).
- ILP solutions can be non-unique — if CARNIVAL returns multiple optimal/near-
  optimal solutions, report the consensus across them (CARNIVAL's built-in
  multi-solution averaging), not a single arbitrary network.

---

## Completion notes

### J1 — Environment + design gate (DONE 2026-06-04)
Stood up the CARNIVAL toolchain and locked the four design forks. Installed
(site-library) `lpSolve` (CRAN, compiled from source) and, via BiocManager,
`CARNIVAL` 2.20.0 + `cosmosR` 1.18.1 (pulling `dorothea`, `bcellViper`); most
heavy deps (decoupleR, OmnipathR) were already present, so the install was small.
`lpsymphony` was already installed but CARNIVAL does not bind it — lpSolve is the
in-R free solver, with a cbc binary as the documented escalation if lpSolve is
too slow on the real (J3) ILPs. Proved the path end-to-end:
`runVanillaCarnival(perturbations, measurements, priorKnowledgeNetwork,
defaultLpSolveCarnivalOptions())` on the package's bundled `toy_*_ex1` data
returned a non-empty `weightedSIF` in 0.26 s. Recorded the housekeeping trap that
lpSolve options default `workdir`/`outputFolder` to the project root with
`keepLPFiles=TRUE` (one `.lp` leaked on the toy run; removed, tree clean) — J3
must redirect these.

Ran an unplanned-but-high-value PKN feasibility probe to settle gotcha 2
(mouse-vs-human PKN) before asking the user to commit to a mode. The mouse
OmniPath signed-directed symbol PKN (18,126 edges / 5,718 nodes) covers all
verdict kinases bar Prkca and all verdict TFs bar Stat1, and directed
reachability from kinase seeds to TF measurements is dense (Gsk3b→Nfkb1/Rela/Spi1
at 1 hop, →Myc/Creb1/Tp53/Jun at 2 hops; every probed kinase reaches every probed
TF within 1–3 hops). Conclusion: the mouse PKN suffices — the human-ortholog
fallback is shelved — and the headline interaction-axis question (Gsk3b wired to
the interaction TFs) and the §18 NF-κB-attenuation path both have prior-network
support for CARNIVAL to test against the signed activities.

Gate outcome (user, AskUserQuestion): (1) **standard per-contrast CARNIVAL** —
plan default; kinases seed `inputObj`, TF activities are `measObj`. (2) **TF
complexes = re-derive `split_complexes=TRUE`** — a deviation from the plan default
(a); the causal layer gets clean single-protein TF nodes via a new
`tf_activity_decoupler_split.rds` (a `--split-complexes` flag on
`build_tf_activity_decoupler.R`), so PKN mapping is identity and NF-κB resolves to
Nfkb1/Rela. (3) **prune PKN to microglia-expressed** — plan default. (4) **add
Suggestive-graded ledger rows** — a deviation from the no-op default; J5 is now
active (Suggestive-only, anti-double-count guard retained). Step bodies J2/J5 and
gotcha 1 were updated to these locks; J2–J5 are unblocked. No knit or new tracked
artifact in J1 (package installs only), so no render-verify step applies. Next
session: **J2** (build the split cache + `R/causal_network.R` + smoke test).

### J2 — causal_network.R + CARNIVAL formatters (DONE 2026-06-04)
Built the causal-network helper and validated it against the live caches
without running the real ILP (that is J3). `R/causal_network.R` carries seven
functions: `build_omnipath_pkn` (assembles the signed-directed mouse PKN as a
CARNIVAL-shaped tibble), `microglia_expressed_symbols` (the J1-locked
expression filter, Ensembl→symbol via the snrnaseq symbol_map), three input
formatters (`tf_meas_from_cache` → `measObj`, `kinase_input_from_cache` →
`inputObj`), the run wrapper `run_carnival_for_contrast`, and the readouts
`summarise_network` + `restrict_network_to_axis`. It is sourced into
`R/helpers.R` after `ccc_inference.R` (capstone of the E/F/G mechanism
modules, before microglia.R) and reuses `.extract_tf_per_modality` from
tf_inference.R. The `--split-complexes` flag was added to
`scripts/build_tf_activity_decoupler.R` and run, writing
`storage/cache/tf_activity_decoupler_split.rds`. Validation is
`scripts/smoke_test_causal_network.R` (26 assertions, **0 fail**): PKN shape
+ sign + no self-loops/conflicts; 8/8 verdict TFs and 4/4 probed kinases
present; expr filter prunes 18,126→10,713 edges and is pure symbol-space;
and a tiny synthetic network (Gsk3b→Nfkb1→{Tnf,Il1b}, Gsk3b⊣Sp3) drives the
real `runVanillaCarnival` solver + result parser end-to-end (4 weight>0 edges
recovered, tidy schema verified).

**Verified CARNIVAL I/O (corrects two J1-era assumptions).** (i) The PKN
column is **`interaction`** (∈{−1,+1}), NOT `sign` — the helper emits
`source/interaction/target`. (ii) lpSolve options have **no `timelimit`
field**; the only runtime levers are `betaWeight` (sparsity penalty) and the
solver choice, so a slow J3 ILP must be bounded by the solver, not a wall.
`weightedSIF$Weight` and `nodesAttributes$AvgAct` are the **consensus across
all optimal solutions** (CARNIVAL's built-in averaging) — this is the
multi-solution honesty the plan's guardrail calls for, available for free.
The wrapper runs `runVanillaCarnival` under a tempdir with
`keepLPFiles=FALSE` (J1 housekeeping trap) and returns an **honest empty
network** (typed status `no_perturbations_in_pkn` / `no_measurements_in_pkn`
/ `solver_error_or_empty`) instead of erroring when the endpoints do not
intersect the PKN — an omitted path stays a finding, not a crash.

**DEVIATION FROM J1 DEC.2 PREMISE — `split_complexes` is a no-op on mouse
CollecTRI.** J1 dec.2 chose `split_complexes=TRUE` on the premise it would
resolve TF complexes to single-protein symbols (NF-κB → Nfkb1/Rela). It does
not: `get_collectri("mouse", split_complexes=TRUE)` is **byte-for-byte
identical** to `=FALSE` (same 1,114 sources, `identical()` source sets, the
NF-κB carrier accession `A0A979HLR9` present in BOTH). The split applies to
the human prior; the mouse build is already maximally split. **Outcome is
unchanged**, so no decision fork: the interaction-contrast `measObj` is 6 TFs
(Bmal1/Clock/Foxo3/Myc/Sfpq/Tbp), all clean single-protein symbols, **0
accession-named** — the accession concern never bites here, and any
accession-named TF would drop on the PKN intersection regardless (same end
state as the plan-default alt-a). The split cache is retained for plan
fidelity (J3 reads it) but is functionally redundant; an optional J6 cleanup
could collapse the two caches. **For J4:** if an NF-κB measurement is wanted,
seed it by the **Nfkb1/Rela symbols** (the `A0A979HLR9` carrier accession is
absent from the PKN and would silently drop).

**Cross-compartment coverage @ interaction contrast (microglia-pruned PKN).**
`measObj` 6 TFs → 5 map (Bmal1 drops); `inputObj` 2 kinases (Csnk1a1⊣, Gsk3b+)
→ 1 maps (Csnk1a1 drops, **Gsk3b retained**). Both endpoints keep ≥1 node, so
the ILP is feasible and the headline interaction-axis Gsk3b question is
testable at J3. Restate at interpretation: the kinase seeds are **bulk
hippocampus phospho**, not microglia-sorted.

**J3 reminder.** The smoke test prunes against the RAW microglia subset for
speed; the J3 build MUST prune against the PROCESSED subset
(`microglia_seurat_processed.rds`) for the real expr filter, run all 5
contrasts under tempdir/`keepLPFiles=FALSE`, and emit
`storage/cache/causal_network.rds`. Next session: **J3**.

### J3 — build_causal_network.R + causal_network.rds (DONE 2026-06-04)
Ran standard per-contrast CARNIVAL over the microglia-pruned PKN and emitted
`storage/cache/causal_network.rds` (8 entries: 5 contrast nets + `summary` +
`pkn_meta` + `params`). The build script and one new helper
(`restrict_pkn_to_reachable` in `R/causal_network.R`) implement the
tractability fix below. The cache-read smoke test confirms the J4 schema: per
contrast `list(nodes[node,node_type,activity,up,down,zero], edges[source,
interaction,target,weight], meta)`; `params` records solver/solver_path/
beta_weight/max_path/threads/timelimit/fdr; `pkn_meta` records the full PKN
size (10,713 edges / 3,779 nodes from 19,698 microglia-expressed symbols).

**CENTRAL FINDING — the J1-locked full-PKN design is ILP-INTRACTABLE; fixed by
a user-approved per-contrast reachability prune (DEVIATION from J1).** J1 locked
"solve over the full microglia-expressed PKN", but that PKN (10,713 edges /
3,779 nodes → ~4,600-binary ILP) does not solve: lpSolve ran >6 min on one
contrast (killed); cbc reported *no feasible solution* within a 600 s limit.
Gated the user, who chose **reachability prune, L=3**. `restrict_pkn_to_reachable
(pkn, seeds, targets, max_path=3)` keeps only nodes on a directed seed→target
walk of ≤ max_path hops (multi-source BFS: forward hop-distance from kinase
seeds over source→target, backward from TF targets over reversed edges; keep iff
`fwd + bwd ≤ max_path`). This is **solution-preserving up to depth max_path**
(every node on any ≤L-hop seed→target path satisfies the inequality, so the
induced subgraph retains those paths intact) and collapses the headline
`nlgf_in_p301s` ILP to 807 edges / 153 nodes, which cbc proves **Optimal in
~4.6 s**. **Completeness caveat (record at interpretation):** paths LONGER than
max_path=3 hops are excluded by construction — max_path is a tunable mechanistic
horizon set to match J1's observed 1–3 hop kinase→TF reachability, NOT an
exhaustive search. It is stored in `params$max_path` and every per-net
`meta$max_path`.

**SOLVER — cbc @ threads=1 (a hard requirement; threads>1 SEGFAULTS).** The
PuLP-bundled cbc 2.10.3 (the 2019 `linux/i64` build, auto-located in `.venv` by
`locate_cbc()`) **segfaults mid branch-and-bound in multi-threaded mode**:
observed at threads=6 on the L=3 `nlgf_in_p301s` ILP, crashing near-optimal
(obj 12.03 vs bound 11.63) so CARNIVAL found no result CSV and the wrapper
correctly logged `solver_error_or_empty` (the honest-empty handling caught the
crash — no garbage network). Single-threaded cbc is stable, deterministic, and
still fast on these tiny pruned problems. The script default is now
`--threads 1` (was 6); `--threads N` remains an override but is documented as
segfault-risky. lpSolve stays selectable via `--solver lpSolve` but is
unusable at this scale.

**Per-contrast outcome (cbc, betaWeight 0.2, FDR 0.10, L=3; all solves
`Optimal`).** `nlgf_in_p301s`: 3 kinase seeds / 22-of-24 TF targets in the
pruned PKN → **33 active edges, 18/22 TFs recovered** (frac 0.82). `interaction`
(headline): 1 seed (Gsk3b) / 5-of-5 targets → **8 edges, 5/5 TFs**, and the
recovered wiring directly answers the headline question — Gsk3b⊣{Clock,Mapk14,
Sfpq}, Mapk14→{Foxo3,Myc,Csnk2b}, Csnk2b→Tbp (a coherent Gsk3b→…→interaction-TF
module with a Gsk3b↔Mapk14 feedback). `tau_in_nlgf`: 1 seed / 13-of-13 targets →
**27 edges, 13/13 TFs**. `nlgf_in_maptki` and `tau_alone`: **honestly empty**
(`no_perturbations_in_pkn`) — 0 kinases reach FDR 0.10 in `phospho_corrected`
for these contrasts, so there is no perturbation seed and the prune returns a
0-row PKN before any solve. Per the anti-anchoring guardrail these empties are
recorded, not dropped. Cache chowned rstudio:rstudio; no `.lp` leaked (tempdir +
keepLPFiles=FALSE held); no knit in J3 (that is J4).

**J4 reminders.** (i) The networks are sparse and signed — `edges$interaction`
∈{−1,+1}, `edges$weight`/`nodes$activity` are CARNIVAL's cross-optimal-solution
consensus (multi-solution honesty for free). (ii) State the L=3 horizon and the
bulk-phospho→microglial-TF cross-compartment bridge wherever paths are read.
(iii) Two of five contrasts are empty — X.1/X.3 must present that as a finding.
(iv) For an NF-κB readout, seed Nfkb1/Rela symbols (the `A0A979HLR9` carrier is
absent from the PKN). Next session: **J4** (rmd/19 → §20).

### J4 — rmd/19_causal_network.Rmd (§20) + wiring + knit (DONE 2026-06-04)
Built the display-only causal-topology chapter (`rmd/19_causal_network.Rmd`,
renders as **§20**), wired it as `child-causal-network` in analysis.Rmd between
`child-human-validation` (rmd/18→§19) and `child-session`, and re-knit the whole
book clean: **0 `class="error"`, 0 `class="warning"`** in 4m44s. New artifacts:
`results/causal_network_{nodes,edges,verdict}.tsv` (62 active nodes / 68 active
edges across the 3 solved contrasts; 3-axis verdict) and one new helper
`plot_causal_network()` in R/causal_network.R. All chowned rstudio:rstudio
(chapter, helper, analysis.Rmd, analysis.html, 3 TSVs).

**Chapter structure (mirrors the §15 kinase analog: X.1 per-contrast / X.2
axis-restricted / X.3 verdict).** Intro states both headline questions, the two
endpoints (kinase `inputObj` ← phospho_corrected; TF `measObj` ←
tf_activity_decoupler_split), and ALL thresholds upfront (FDR<0.10 both ends;
microglia-pruned OmniPath mouse PKN 10,713e/3,779n; L=3 reachability horizon;
cbc single-thread, betaWeight 0.2; consensus across optimal solutions), plus the
bulk-phospho→microglial-TF cross-compartment caveat and the anti-anchoring
"omitted edge is a finding" lock. 20.1 = augmented per-contrast summary kable
(seeds_in→active, TF_meas_in→recovered, intermediates, edges) + one
`plot_causal_network` node-link figure per solved contrast + honest reporting of
the 2 empty contrasts. 20.2 = axis-restricted membership table. 20.3 = per-axis
verdict table + prose answering (a)/(b).

**Headline findings (solver-faithful; every path verified against the cache
edge lists before writing prose).** (a) **YES** — the `interaction` net recovers
`Gsk3b ⊣ Mapk14 → Myc` (and `→ Foxo3`): the §15 interaction-axis lead kinase is
wired to the §14.3 interaction-axis lead TF through a single inferred Mapk14
(p38α), all 5 measured TFs recovered. (b) **MIXED/honest** — `tau_in_nlgf`
recovers `Gsk3b ⊣ Ikbkg(NEMO) → Nfkbib(IκBβ)` (a candidate upstream route to the
NF-κB regulatory module, NEMO-inhibition direction-consistent with §18
attenuation); `nlgf_in_p301s` recovers only `Gsk3b → Nfkb1` (activation, the
amyloid-on-tau direction); the `interaction` contrast recovers **no** NF-κB node.
So §18 attenuation has a recoverable path on `tau_in_nlgf` but not on the
interaction contrast itself. `nlgf_in_p301s` is the best-covered net (18/22 TFs,
2 seed modules off Gsk3b + Cdk5, Spi1 recovered 3 hops down via inferred Gata2;
3rd seed Ppp1ca active but isolated). `nlgf_in_maptki`/`tau_alone` honestly empty
(0 kinase seeds at FDR<0.10, matching §15.1).

**KEY DESIGN FINDING — gene-universe axis-restriction does not partition the
wiring; it is a cross-axis BRIDGE (X.2).** The plan's "restrict each network to
the 3 axis universes" via `induced` mode (both endpoints in universe) yields **0
edges for the interaction net** because the recovered regulators scatter across
all three universes — only Gsk3b is an interaction-universe member, while
Mapk14/Myc/Foxo3 annotate to the AMYLOID universe and Gsk3b/Mapk14/Clock/Sfpq to
the SYNAPTIC universe. Resolution (kept the plan's helpers, reinterpreted the
readout): X.2 reports universe MEMBERSHIP as an annotation + the
`induced`-vs-`incident` edge-count gap as the bridging metric, rather than
forcing the degenerate induced subgraph. The axes are fundamentally CONTRAST
groupings (amyloid & synaptic both → {nlgf_in_maptki(empty), nlgf_in_p301s};
interaction → {interaction}), so amyloid/synaptic share their only solved net
(nlgf_in_p301s) — the same entanglement §15.2 reported — and synaptic-suppression
contributes no distinct wiring (honest non-finding, consistent with §14.3/§15.2
engulfment hypothesis). tau_in_nlgf/tau_alone sit outside the 3-axis framework.

**`plot_causal_network()` (R/causal_network.R, self-contained ggraph helper, smoke
-tested to PNG before knitting).** Draws ACTIVE nodes only (activity≠0); shape by
CARNIVAL role (kinase seed=triangle, intermediate=circle, TF=square), fill by
activity sign (up=red/down=blue), edge colour by interaction sign
(activation=green/inhibition=red) with arrowheads. Uses `stress` layout (tolerates
the Gsk3b↔Mapk14 feedback cycles without DAG-breaking messages) and
`max.overlaps=Inf` (so ggrepel never silently drops a label — would have on the
32-node nlgf_in_p301s net). Returns NULL on an empty net so the caller emits an
honest "no network" note. Decision: the Python `scientific-visualization` skill
was NOT used — the figure must render inside the shared R knit session, so a
ggraph helper in the chapter's own module is the correct house pattern (mirrors
§15's plot_* fns living in kinase_inference.R).

**J5 reminder (next session, FRESH — design-heavy, NOT a trivial extension).** J1
LOCKED "add Suggestive-graded ledger rows". J5 appends Suggestive topology
claim(s) for a genuinely novel, independently-gradeable recovered path (candidate:
the `Gsk3b ⊣ Mapk14 → Myc` interaction-axis bridge, or the `Gsk3b → NEMO/IκBβ`
NF-κB-module route — both join already-counted endpoints causally, so they ARE
incremental topological evidence, not a re-count of the TF/kinase activities
themselves) via build_biological_model_ledger.R → build_biological_model_
adjudication.R; then re-run the 5 schema invariants + adjudication (verify contest
margins shift only as intended; NEVER reopen the 11-entity set or grade defs) and
re-knit §17 + §20. Anti-double-count guard: Suggestive-only. Verdict TSV
`results/causal_network_verdict.tsv` is the structured input for the ledger row.

### J5 — Ledger feed (DONE 2026-06-04)
Fed the recovered causal topology into the §17 biological-model claims ledger as
**two Suggestive-graded rows**, re-ran the adjudication, fixed the dependent §17
+ executive-summary prose, and re-knit the book clean (**0 error / 0 warning,
4m50s**). The ledger is now **83 rows** (75 Phase H + 6 Phase I + 2 Phase J).

**The two rows (user-chosen scope: "both", J-001 contradicts Hyp-3A).** `J-001`
(interaction_metabolic / cross_layer) scores the headline interaction bridge
`Gsk3b ⊣ Mapk14(p38α) → Myc / Foxo3` — the §15 interaction-lead kinase wired to
the §14.3 interaction-lead TFs through one inferred intermediate; supports
`Hyp-3B;T-Synergy` and **contradicts `Hyp-3A`** (the bridge is interaction-only:
the single-insult arms returned 0-seed empty nets, so a recovered
interaction-specific path is incremental evidence against "no interaction-specific
mechanism", mirroring H2-075 at the topology layer). `J-002` (cross_axis /
cross_layer) probes the §18 NF-κB attenuation route — `Gsk3b ⊣ Ikbkg(NEMO) →
Nfkbib` @tau_in_nlgf (direction-consistent), `Gsk3b → Nfkb1` @nlgf_in_p301s
(activation), and **no NF-κB node at the interaction contrast** where the
attenuation is defined; supports `T-Tau-attenuates` **support-only (no
contradict)** because the topological support is equivocal — present off the
interaction, absent on it (the interaction-contrast absence is reported as a
finding per the anti-anchoring lock, not smoothed over).

**Anti-double-count compliance.** Both rows are **Suggestive** and score the
recovered directed *signed wiring* between endpoints, NOT the endpoint activities
already counted (Gsk3b H2-051, Myc H2-046, NF-κB attenuation I-001..I-006). Each
is held at Suggestive by four documented caveats: the cross-compartment
bulk-phospho seed, the L=3 reachability horizon, the single inferred intermediate
(Mapk14 / NEMO), and the multi-solution-consensus weighting. The 11-entity set
and grade definitions were untouched.

**Adjudication deltas (verified; margins shift only as intended).** Interaction
contest Hyp-3B vs Hyp-3A margin **53 → 55** (Hyp-3B support 27→28 / net →28;
Hyp-3A contradict 26→27 / net →−27). Amyloid (margin 18) and synaptic (margin 12)
contests **unchanged**; all three `no_tie_break_needed`. Cross-axis entities:
T-Synergy support 28→29, T-Tau-attenuates support 11→12; both single-axis shares
nudged (T-Synergy 89%→90%, T-Tau-attenuates cross_axis 73%→75%). Distribution
shifts: Suggestive grade 12→14, interaction_metabolic axis 25→26, cross_axis
11→12, cross_layer layer 8→10. All 5 schema invariants pass; 83 rows, 0 duplicate
IDs. (The phase_j ledger block + the adjudication header bump to 83 were authored
the same day; this step re-verified them and completed the consuming prose + knit.)

**Pre-existing §17 bug fixed (bonus).** §17.3's axis-3 verdict attributed the
1-row Hyp-3B-support / Hyp-3A-contradict gap to **"H2-067"**, but H2-067 in fact
*contradicts* Hyp-3A (supports `Hyp-3B;T-Synergy`, contradicts `Hyp-3A`) — it is
not a gap row. Computed the true gap row: **H2-063**, a Moderate LR-layer
microglia-involving-rate context claim (51/100 axis-3 LR pairs microglia-involving)
with an empty `contradicts_models` cell *because it is not directly
model-discriminating* — a far better-grounded explanation than the bogus
"deliberate H2 caveat" the old prose invoked. The gap stays 1 row after Phase J
(J-001 supports Hyp-3B AND contradicts Hyp-3A, so it lifts both counts in step).
Also corrected the mathematically-inconsistent "8 of 17 … carry the Strong grade"
to "8 of the 28 … (17 Moderate, 3 Suggestive)".

**Files touched.** `rmd/16_biological_model.Rmd` (§17): ledger count + Phase J
composition (intro + DT comment), DT `lengthMenu` 81→83, axis-3 margin/grade/gap
sentence, cross-axis theme tallies + ranking + Hyp-3A net, single-axis range
89-94%→90-94%, T-Tau-attenuates Phase J increment. `analysis.Rmd`: the
headline-integrated-biological-model executive summary (81→83 atomic claims +
Phase J source, margin 53→55, T-Synergy 28→29, T-Tau-attenuates 11→12 with Phase J
attribution). Grade/axis/layer distribution tables in §17.2 render dynamically
(`fmt_dist`) so they update on knit with no hardcode edit. The Cdk5 KSN
substrate-footprint "81 / 263 / 153 sites" (biology, coincidentally 81) was left
untouched. All root-owned outputs chowned rstudio:rstudio; stale omnipathr-log
files pruned to the current run.

**J6 reminder (next, FRESH session — close-out).** Update `map.md` (no new
source files in J5; note the 2 Phase J ledger rows feed §17 from the §20 causal
net — the ledger/adjudication scripts + the 3 biological_model TSVs are already
mapped) + `completed/DIGEST.md` (new J-arc entry: goal, per-step outcomes, the 4
J1 locks, the L=3-prune + cbc-single-thread + split-complexes-noop deviations,
artifacts, the H2-067→H2-063 fix, verdict). Final knit-verify. Commit. Move this
file to `completed/causal_network_plan_2026-06-04.md`.

### J6 — Close-out (DONE 2026-06-04)
Closed out the causal-network (J) arc: documentation refreshed, ownership fixed,
build verified clean, committed, plan archived. No source/cache/Rmd was touched —
J6 is pure book-keeping over the J1–J5 deliverables, which were already complete
and committed (the J5 commit `3e5f29e` carried the clean §17+§20 knit).

**`map.md` (6 edits).** Added the `19_causal_network → §20` pipeline-table row; the
`19→§20` verified anchor; the `causal_network.R` helper block (9 public fns) after
`ccc_inference.R` with a note that it reuses `.extract_tf_per_modality` +
`build_axis_gene_universe` from tf_inference.R; a "Causal lane (§20)" sentence in
the scripts section recording the two hard build constraints (cbc @ threads=1 —
threads>1 segfaults; L=3 reachability prune for ILP tractability) + the
`--split-complexes` no-op; and two cache-table rows (`tf_activity_decoupler_split.rds`
[S], `causal_network.rds` [S]). **`DIGEST.md`:** new `### causal_network_plan_2026-06-04`
entry (goal, the J1 gate decisions + the three deviations — L=3 prune,
cbc-single-thread, split-complexes-noop — verdict (a)/(b), the X.2 bridge finding,
the 2 Suggestive rows + margin deltas, the H2-067→H2-063 fix, artifacts, gotchas),
plus the narrative-arc header extended with **J** (and CARNIVAL removed from the
candidate-follow-on list since it is now done).

**Ownership.** `scripts/build_causal_network.R` was still root-owned (J3 wrote it
as root) → chowned rstudio:rstudio; every other J artifact was already
rstudio:rstudio. **Knit-verify.** The committed `analysis.html` is the J5 build and
greps **0 `class="error"` / 0 `class="warning"`**; J6 edits only `storage/notes/`
files, which are NOT knit inputs (the children are `rmd/*.Rmd`), so the J5 knit
remains valid and a re-knit would byte-reproduce it — skipped as wasteful, the
clean-knit invariant holds. **Archive.** This file `git mv`-d to
`completed/causal_network_plan_2026-06-04.md`. The causal-network plan is COMPLETE;
no active plan remains on disk — the next session surveys + proposes (see DIGEST
candidate follow-ons: MOFA+ latent integration, SCENIC regulons, or ROSMAP human
confirmation).

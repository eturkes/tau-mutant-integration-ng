# Codebase map (read this to orient; saves re-grepping)

Navigation aid for `tau-mutant-integration-ng`. Source files are truth; this
encodes the wiring (pipeline order, cross-child globals, cache producer→consumer)
that is otherwise only recoverable by grepping. Refresh cheaply with the greps in
the last section if it drifts.

Project = integrate snRNAseq + GeoMx + bulk proteomics + bulk phosphoproteomics
across 4 mouse genotypes (MAPTKI, P301S, NLGF_MAPTKI, NLGF_P301S). Readout =
interaction contrast `(NLGF_P301S − P301S) − (NLGF_MAPTKI − MAPTKI)`.

## Pipeline (analysis.Rmd child order = execution order)

One knit shares ONE R session across all children, so a global assigned in an
earlier child is in scope for every later child. **Reference sections by rmd
file, not by rendered §number**: rmd 08/09 were deleted, so file-number ≠
§number. Verified anchors: rmd/12→§13, 13→§14, 15→§16, 16→§17, 17→§18, 18→§19, 19→§20, 20→§21, 21→§22, 22→§23, 23→§24 (99_session→§25).

| rmd file | title | builds (global / cache) | reads |
|---|---|---|---|
| 01_data | Data loading | `microglia_seurat` (raw subset, K:microglia_seurat_raw), `symbol_map` (K:snrnaseq_symbol_map), `geomx`, `prot_raw`, `phos_raw`, `sample_key` | storage/data/{snrnaseq,geomx}.rds, proteomics/phospho TSVs |
| 02a_snrnaseq_qc | microglia reprocess + QC | `microglia_seurat` (processed+labelled; K:microglia_seurat_processed) | — |
| 02b_snrnaseq_de | pseudobulk + NEBULA DE | `de_snrnaseq` (K), `de_snrnaseq_nebula` (K) | microglia_seurat |
| 02c_snrnaseq_triangulation | 6-method consolidation | `sn_tri_pb`,`sn_tri_tmb` (K) | S:de_snrnaseq_triangulation_pb, de_snrnaseq_glmmtmb |
| 02d_snrnaseq_substate_de | per-substate volcano | — | de_snrnaseq_nebula_per_state |
| 02e_snrnaseq_substate_pathway | per-substate fGSEA | `fgsea_per_state` (K), collection getters | S:de_snrnaseq_nebula_per_state_1pct; fgsea collections |
| 03_hdwgcna | hdWGCNA modules | — | S:hdwgcna_microglia.rds, hdwgcna_module_de.rds |
| 03b_hdwgcna_soft_power | soft-power scan | `hdw_scan` | S:hdwgcna_microglia_{p6,p8}.rds |
| 04_geomx | GeoMx norm + DE | `de_geomx` (K) | geomx |
| 05_proteomics | proteomics DE | `de_proteomics` (K) | prot_raw, sample_key |
| 06_phospho | phospho site DE | `de_phospho` (K) | phos_raw, sample_key |
| 06b_phospho_corr | protein-corrected phospho | `de_phospho_corr` (K:de_phospho_corrected) | de_phospho, de_proteomics |
| 07_integration | cross-modality + GO BP fGSEA | `integration_tbl` (K:integration_table), `gsea_gobp` (K:fgsea_gobp_results), `gobp` | all de_* globals |
| 10_divergence | tau-context NLGF divergence | `milo_obj` (K:milo), `milo_results` (K:milo_da), heatmap tbls | de_snrnaseq_nebula, microglia_seurat, integration_tbl |
| 11_ccc | CCC: CellChat + MultiNicheNet | (display) | S:cellchat_*, multinichenet_output, hdwgcna_microglia, symbol_map |
| 12_pathway_survey | agnostic pathway/module survey | (leader board) | S:fgsea_{gomf,gocc,custom_*}_results |
| 13_tf_inference | TF activity (decoupleR+CollecTRI) | (§14 verdict) | S:tf_activity_decoupler, msigdb_*_mouse |
| 14_kinase_inference | kinase activity (decoupleR+OmniPath) | (§15 verdict) | S:kinase_activity_decoupler, msigdb_*_mouse |
| 15_ccc_mechanism | CCC mechanism (3-tool LR) | (§16 verdict) | S:liana_output, cellchat_*, multinichenet_output |
| 16_biological_model | integrated claims ledger | (§17; DT::datatable) | S:results/biological_model_*.tsv |
| 17_nfkb_attenuation | per-state NF-κB attenuation | (§18; K:per_state_nfkb_target_gsea) | S:per_state_tf_activity, de_..._per_state_1pct, collectri_mouse_for_nfkb_gsea |
| 18_human_validation | human SEA-AD cross-species validation | (§19; display-only, read-only) | S:summary_human_validation, summary_human_robustness, human_substate_{crosstab,conservation_metrics,percell}, human_validation_signature_membership |
| 19_causal_network | CARNIVAL causal signalling topology | (§20; display-only) | S:causal_network (per-contrast kinase→TF ILP nets) |
| 20_scenic_regulons | SCENIC data-driven regulons (arc K) | (§21; display-only; writes results/scenic_verdict.tsv) | S:scenic_summary.rds, scenic/scenic_recurrence_dist.tsv |
| 21_spatial_deconvolution | GeoMx tissue-composition (arc L) | (§22; display-only; writes results/spatial_decon_verdict.tsv) | S:spatial_decon.rds, results/spatial_decon_{contrasts,abundance_by_genotype,spatial_autocorr}.tsv |
| 22_trajectory | microglial activation pseudotime (arc M) + gene-level differential dynamics (arc O) | (§23; display-only; writes results/trajectory_verdict.tsv) | S:trajectory.rds (15 slots; ct_get/bg_get lookups) + trajectory_dynamics.rds (arc-O; 11 slots; td$/dyn_* lookups) |
| 23_celltype_specificity | cross-cell-type interaction specificity (arc N) | (§24; display-only; writes results/celltype_specificity_verdict.tsv) | S:celltype_specificity.rds (11 slots; tget/xget/pget lookups) |
| 99_session | session info + manifest | — | — |

## Cross-child globals (the "check the producer chunk" table)

`microglia_seurat`←01_data then re-bound in 02a · `symbol_map`←01_data ·
`de_snrnaseq`,`de_snrnaseq_nebula`←02b · `de_geomx`←04 · `de_proteomics`←05 ·
`de_phospho`←06 · `de_phospho_corr`←06b · `integration_tbl`,`gsea_gobp`,`gobp`←07 ·
`fgsea_per_state`←02e. All DE objects share shape `de$fit$top[[contrast]]`
(NEBULA: `de_snrnaseq_nebula$top[[contrast]]`).

## Caches: producer → consumer ([K]=built in-knit by cache_or_run, auto-loaded; [S]=built by scripts/*.R outside the knit, pre-build before knitting)

| cache (storage/cache/) | producer | consumed by |
|---|---|---|
| microglia_seurat_raw.rds | [K] 01_data | 02a |
| microglia_seurat_processed.rds | [K] 02a (or [S] build_seurat_full) | 02b, hdwgcna/nebula/triangulation scripts |
| snrnaseq_symbol_map.rds | [K] 01_data | 11_ccc, hdwgcna/nebula scripts |
| de_snrnaseq.rds, de_snrnaseq_nebula.rds | [K] 02b | 07, 10, fgsea/custom/tf scripts |
| de_snrnaseq_nebula_per_state.rds | [S] build_per_state_nebula | 02d, build_tf_activity_decoupler |
| de_snrnaseq_nebula_per_state_1pct.rds | [S] build_per_state_nebula_1pct | 02e, build_per_state_tf_activity (§18) |
| de_snrnaseq_triangulation_pb.rds, de_snrnaseq_glmmtmb.rds | [S] build_triangulation_caches | 02c |
| de_geomx.rds | [K] 04 | 07, fgsea/custom/tf scripts |
| de_proteomics.rds | [K] 05 | 06b, 07, fgsea/custom/tf scripts |
| de_phospho.rds | [K] 06 | 06b, 07, kinase scripts |
| de_phospho_corrected.rds | [K] 06b | 07, kinase scripts |
| integration_table.rds | [K] 07 | 10 |
| fgsea_gobp_results.rds | [K] 07 | 07 |
| fgsea_{gomf,gocc}_results.rds | [S] build_fgsea_extra_collections | 12 |
| msigdb_{gobp,gocc,gomf}_mouse.rds | [S] build_fgsea_extra_collections / fgsea.R getters | 13, 14, 15 |
| fgsea_per_state_results.rds | [K] 02e | 12 |
| custom_{microglia_states,microglia_ad,module_sources}.rds | [S] build_custom_* | fgsea.R getters |
| fgsea_custom_{states,ad,modules}_results.rds | [S] build_custom_* | 12 |
| hdwgcna_microglia{,_p6,_p8}.rds, hdwgcna_module_de*.rds | [S] build_hdwgcna | 03, 03b, 11 |
| milo.rds, milo_da.rds | [K] 10 | 10 |
| seurat_full_processed.rds (all cell types, not just microglia) | [S] build_seurat_full | build_cellchat/liana/multinichenet, build_spatial_deconvolution (§22 reference profile) |
| cellchat_{merged,per_condition}.rds | [S] build_cellchat | 11, 15 |
| liana_output.rds | [S] build_liana_cache | 15 |
| multinichenet_output.rds | [S] build_multinichenet | 11, 15 |
| tf_activity_decoupler.rds | [S] build_tf_activity_decoupler | 13 |
| tf_activity_decoupler_split.rds | [S] build_tf_activity_decoupler --split-complexes | build_causal_network.R (§20; split is a no-op on mouse CollecTRI) |
| kinase_activity_decoupler.rds | [S] build_kinase_activity_decoupler | 14 |
| per_state_tf_activity.rds | [S] build_per_state_tf_activity | 17 (§18) |
| collectri_mouse_for_nfkb_gsea.rds | [S] build_human_validation_signatures | 17 (§18) |
| per_state_nfkb_target_gsea.rds | [K] 17 | 17 |
| human_validation_signatures.rds | [S] build_human_validation_signatures | human cross-species plan (H3 reads via JSON bridge) |
| human_seaad_donor_neuropath.csv (storage/cache) | [S] build_human_microglia (H1) | H4 neuropath join (donor×region; numeric Thal=amyloid, Braak=tau) |
| human_substate_{percell.csv.gz, score_means.csv, pseudobulk_counts.csv.gz, pseudobulk_samples.csv} | [S] build_human_substate_conservation (H3) | H4/H5/H6 (score_means = primary H4 input; join key matches neuropath csv) |
| summary_human_validation.rds (+ results/human_interaction_models.tsv) | [S] build_human_interaction_models (H4) | rmd/18 §19 (forest plot, collinearity map, mechanism heatmap); summary.Rmd §7 (mechanism-concordance heatmap) |
| human_substate_aucell_score_means.csv | [S] build_human_aucell_rescoring (H5; AUCell over H3 cells, decoupler) | build_human_robustness_mediation (H5 ARM B scoring-method robustness) |
| summary_human_robustness.rds (+ results/human_robustness_mediation.tsv) | [S] build_human_robustness_mediation (H5) | rmd/18 §19 (3-estimator pb/sc/auc concordance panel; mediation ACME / moderation-vs-mediation contrast) |
| causal_network.rds | [S] build_causal_network.R | 19 (§20); 5 contrast nets (3 solved, 2 honestly empty) + summary + pkn_meta + params |
| scenic_summary.rds | [S] build_scenic_contrasts.R | 20 (§21); 11-component list (regulons, scenic_net, headtohead, scenic_activity, aucell_contrasts, aucell_pb, substate_int, target_overlap, recovery, census, params) |
| scenic/ subdir (microglia_counts.mtx + genes/cells/colattrs, microglia.loom, adj_seed{1..10}.tsv, reg_seed{1..10}.csv, scenic_regulons.tsv, aucell.csv.gz, scenic_consensus_regulons.gmt, scenic_recurrence_dist.tsv, scenic_recovery_census.tsv) | [S] SCENIC py lane: export_microglia_for_scenic.R → build_scenic_grn.py → build_scenic_ctx_aucell.py | build_scenic_contrasts.R (→ scenic_summary.rds) |
| spatial_decon.rds | [S] build_spatial_deconvolution.R (reads seurat_full_processed + geomx) | 21 (§22); two-stage SpatialDecon (stage1_broad 6-level + stage2_substate 4-microglial), 5-contrast factorial log-abundance + per-slide Moran's I |
| trajectory.rds | [S] build_trajectory.R (reads microglia_seurat_processed + symbol_map; invokes build_trajectory_python.py) | 22 (§23, arc-M); 15 slots — traj_clean/traj_all (Slingshot), per_cell pseudotime, root_validation (entropy/n_genes/CytoTRACE2 — all reject homeostatic-as-potent → activation ordering), per_replicate, progression-vs-composition decomposition, summary_matrix, 5-contrast fits, cross-method concordance, by_genotype, params. Headline: tau×amyloid pseudotime interaction carried by PROGRESSION not COMPOSITION. Also writes results/trajectory_{contrasts,pseudotime_by_genotype,progression_decomposition,method_concordance}.tsv |
| trajectory/ subdir (embedding.csv + obs.csv [R→py bundle], python_pseudotime.csv [dpt + cellrank_dam_fate], python_provenance.json) | [S] build_trajectory.R exports → build_trajectory_python.py (scanpy PAGA/DPT + CellRank2 velocity-free; `.venv`) writes back | build_trajectory.R (→ trajectory.rds) |
| trajectory_dynamics.rds | [S] build_trajectory_dynamics.R (reads microglia_seurat_processed + trajectory.rds clean lineage + symbol_map; heavy fitGAM on the ≥1% genome-wide gene set, ~11.5k genes) | 22 (§23, arc-O); 11 slots — interaction (per-gene 2×2 diff-of-diffs Wald), omnibus (conditionTest), association (associationTest), vs_static (dynamic-vs section-24 matched-power static), smoothers (per-genotype fitted curves), arc_m_synergy (M-002 link), on_lineage_counts, params (nknots=6, padj_cut=0.10, 110 sig of 11,466 fitted). Headline: 110 gene-level interaction-dynamics genes ALL static-null → non-redundant layer, but ¼ ambient+ribosomal-confounded, Gsk3b/Myc absent. Writes results/trajectory_dynamics_{interaction,omnibus,association,vs_static}.tsv |
| celltype_specificity.rds | [S] build_celltype_specificity.R (reads seurat_full_processed [SYMBOL-keyed, symbol_map=NULL]; cross-checks vs de_snrnaseq + de_snrnaseq_nebula) | rmd/23 §24 arc-N; 11 slots — fits (NATIVE/MATCHED × pseudobulk/NEBULA × 6 units {Astrocyte,Microglia,Neuronal,Oligodendrocyte,OPC,Vascular}, slimmed to top/n_cells/n_genes), tally (R1/R2), interaction_concordance (R3), specificity_class (R4), pathway_tally + microglia_pathway_cross_unit (R5 GO-BP fGSEA at interaction), crosscheck (10 ρ all PASS 0.90; nebula interaction 0.917 nearest floor), headline_genes (39), unit_cell_counts, per_unit_min (K=289), params (SEED=1, MATCHED downsample). Headline: tau×amyloid interaction is NOT microglia-specific — at MATCHED power gene-level interaction collapses to ~0 in ALL units incl Microglia; R5 pathway response is Neuronal-dominated (depth-confound caveat). Also writes results/celltype_specificity_{tally_native,tally_matched,interaction_concordance,specificity_class,pathway_tally,crosscheck}.tsv |
| celltype_specificity_fits/ subdir (native_nebula_<unit>.rds × 6) | [S] build_celltype_specificity.R per-unit resume caches (lightest-first; Neuronal 163k-cell NATIVE NEBULA is the ~30-min bottleneck) | build_celltype_specificity.R (→ celltype_specificity.rds) |

Biological-model TSVs (no cache): `build_biological_model_ledger.R` →
`results/biological_model_claims_ledger.tsv`;
`build_biological_model_adjudication.R` → `..._adjudication.tsv`,
`..._contest_verdicts.tsv`. Consumed by rmd/16. Ledger rows (92 total) are
hardcoded `row()` calls grouped by phase (H/I/J/K/L/M/N/**O**); the 2 Phase K rows
(K-001 Spi1, K-002 Rel) are the SCENIC §21 → §17 corroboration feed; the 2
Phase L rows (L-001 DAM tissue-compartment expansion supporting BOTH amyloid
models + T-Inflammation, L-002 null composition interaction with empty
supports/contradicts) are the spatial-decon §22 → §17 feed; the 2 Phase M rows
(M-001 amyloid-driven homeostatic→DAM progression supporting BOTH amyloid
models + T-Inflammation, M-002 the POSITIVE+significant tau×amyloid progression
interaction supporting **T-Synergy ONLY** — recorded thematically, NOT scored
against Hyp-3A/Hyp-3B, per the user-approved option-B M5 gate) are the
trajectory §23 → §17 feed; the 2 Phase N rows (N-001 amyloid main effect
microglia-LED under cell-count matching supporting BOTH amyloid models +
T-Inflammation, N-002 the interaction-is-NOT-microglia-specific reframing with
EMPTY supports/contradicts — feeds-no-model like L-002, NOT scored against
Hyp-3A/Hyp-3B/T-Synergy because the per-unit interaction profiles are distinct,
per the user-approved **option-C N5 gate** = feed both halves) are the
cross-cell-type specificity §24 → §17 feed; the 1 Phase O row (O-001 the gene-
level differential-dynamics interaction — 110 tradeSeq genes whose pseudotime
SHAPE is tau×amyloid-modulated, ALL static-null at matched power — with EMPTY
supports/contradicts, feeds-no-model like L-002/N-002 per the user-approved
**ALT-1 O5 gate** = always-margin-neutral; the pre-registered overlap test
locked before interpretation found Gsk3b/Myc ABSENT + axon-guidance 1/110≈chance
→ scored-IF-branch not triggered; held Suggestive by ¼ ambient+ribosomal
confound; it is the gene-level COMPANION of M-002 so scoring T-Synergy would
double-count) is the trajectory-dynamics §23 → §17 feed. All five feeds
Suggestive, existence/dynamics/specificity-not-activity, margin-neutral (contests
stay 18/12/55). N-001 lifts T-Inflammation 30→31, edging it past T-Synergy on raw
count (T-Synergy keeps more Strong support; O-001 adds to neither). `layer=
"composition"` (arc-L), `layer="dynamics"` (arc-M), `layer="specificity"`
(arc-N), and `layer="dynamics_de"` (arc-O, per-gene companion to the aggregate
`dynamics`) are the additions to the builder's layer validation; the adjudication
script ignores `layer` so only the builder needed extending.

## R/ helpers (sourced via R/helpers.R in this dependency order)

`constants.R` (load first; symbols every file uses) · `utils.R`
(`cache_or_run`, `%||%`, `write_tsv_safe`, `isTRUE_vec`) · `io.R`
(`proteomics_sample_meta`, `match_intensity_columns`, `symbols_to_ensembl`) ·
`design.R` (`make_contrast_matrix`, `factorial_design`) · `de_pb.R` (pseudobulk +
bulk DE: `build_pseudobulk`, `fit_limma_voom`, `fit_limma_log`,
`fit_edger_qlf`, `fit_deseq2_pb`, `fit_dream_pb`, `prevalence_filter`,
`median_normalise`) · `de_sc.R` (`fit_nebula_microglia`, `assemble_nebula_top`,
`fit_nebula_per_state`, `fit_glmmtmb_microglia`) · `fgsea.R` (collection getters
`get_gobp/gomf/gocc`, `get_custom_microglia_states/_ad`, `get_custom_module_sources`;
`run_fgsea_*`, `join_fgsea_results`) ·
`pathway_survey.R` (`rank_pathways_cross_modality`, `rank_pathways_per_substate`,
`build_leader_board`, `format_*_table`) · `tf_inference.R` (`extract_de_stat_*`,
`run_decoupler_per_modality`, `split_decoupler_by_contrast`, `rank_tfs_cross_modality`,
`build_axis_gene_universe`, `score_tf_per_axis`, `build_tf_verdict_table`) ·
`kinase_inference.R` (`fetch_omnipath_ksn_human`, `map_ksn_human_to_mouse`,
`run_decoupler_per_cache`, `score_kinase_per_axis`, `build_kinase_verdict_table`) ·
`ccc_inference.R` (`seurat_subset_to_anndata`, `run_liana_rank_aggregate_for_genotype`,
`derive_liana_per_contrast_list`, `rank_lr_cross_tool`, `score_lr_per_axis`,
`build_lr_verdict_table`) · `causal_network.R` (§20 CARNIVAL layer:
`build_omnipath_pkn`, `microglia_expressed_symbols`, `restrict_pkn_to_reachable`,
`tf_meas_from_cache`, `kinase_input_from_cache`, `run_carnival_for_contrast`,
`summarise_network`, `restrict_network_to_axis`, `plot_causal_network`) ·
`scenic.R` (§21 arc-K data-driven layer; sourced AFTER causal_network.R; reuses
tf_inference.R's decoupleR fns for the network swap + design.R/de_pb.R for the
AUCell factorial: `build_scenic_network`, `run_scenic_decoupler`,
`decoupler_activity_long`, `build_scenic_headtohead`, `aucell_to_pseudobulk`,
`fit_aucell_contrasts`, `scenic_substate_activity`,
`scenic_collectri_target_overlap`, `build_recovery_table`) ·
`spatial_decon.R` (§22 arc-L GeoMx tissue-composition layer; sourced AFTER
scenic.R: `build_reference_profile`, `derive_geomx_background`,
`run_spatialdecon`, `combine_two_stage`, `fit_abundance_contrasts`,
`spatial_autocorrelation`) ·
`trajectory.R` (§23 arc-M aggregate activation-trajectory + arc-O gene-level
differential-dynamics layers; sourced AFTER spatial_decon.R, before microglia.R.
arc-M: `build_microglia_trajectory` [Slingshot], `cell_potency` +
`validate_root_potency` [root check], `pseudotime_per_replicate`,
`decompose_progression_vs_composition` [guardrail-#1], `fit_trajectory_contrasts`
[reuses factorial_design + fit_limma_log], `pseudotime_concordance`,
`plot_trajectory_umap`, `plot_pseudotime_density`, `plot_interaction_forest`.
arc-O (tradeSeq NB-GAM per-gene dynamics): `fit_lineage_gam` [conditioned fitGAM,
conditions=genotype, nknots=6], `interaction_dynamics_contrast` [2×2 difference-
of-differences Wald in coefficient space, weights MAPTKI+1/P301S−1/NLGF_MAPTKI−1/
NLGF_P301S+1], `omnibus_dynamics` [conditionTest], `association_dynamics`
[associationTest], `extract_condition_smoothers` [per-genotype fitted curves for
display], `plot_condition_smoothers` [faceted per-genotype smoother ribbons]) ·
`celltype_specificity.R` (§24 arc-N cross-cell-type specificity layer; sourced
AFTER trajectory.R; depends only on de_pb.R + design.R + de_sc.R: `celltype_unit_labels`
[9-level cell_type → 6 units], `min_cells_per_replicate` [K finder → 289],
`subset_pseudobulk_de` / `subset_nebula_de` [byte-for-byte 02b limma-voom / generic
`fit_nebula_microglia` on a unit subset; both lift to uniform `$top[[contrast]]`],
`downsample_balanced` [seeded per unit×replicate cell sampling to K], `specificity_tally`,
`interaction_concordance`, `cross_estimator_concordance`, `specificity_class_table`,
`microglia_crosscheck` [symbol-bridged sanity gate vs Ensembl de_snrnaseq*],
`assemble_specificity_tables` [R1–R4 orchestrator], `specificity_pathway_fgsea` /
`specificity_pathway_tally` / `microglia_pathway_cross_unit` [R5 GO-BP fGSEA at
interaction], + 4 plot fns: `plot_specificity_tally_heatmap`,
`plot_interaction_logfc_heatmap`, `plot_cross_estimator_panel`, `plot_specificity_class`) ·
`microglia.R` (`label_microglia_states`,
`flag_contaminant_hubs`) · `hdwgcna.R` (`run_hdwgcna_pipeline`, `fit_module_de[_per_state]`,
`module_enrichment`) · `plot.R` (`concordance_plot`) · `report.R` (`report_theme()`
— shared bslib v5 HTML report theme; chrome only, no pipeline globals/caches;
also wired self-contained into each report's YAML `theme: !expr`, so stands alone
at parse time).

`build_axis_gene_universe()` (tf_inference.R) is reused by all three mechanism
layers and by `causal_network.R` (§20, which also reuses `.extract_tf_per_modality`
from tf_inference.R for its TF `measObj` builder). tf/kinase/ccc_inference.R are
mutually self-contained (cross-refs in helpers.R are documentation-only). `library(nichenetr)` is attached (not `::`)
because `convert_human_to_mouse_symbols` needs lazy-loaded `geneinfo_human`.

## scripts/ (run outside the knit; build heavy caches above)

build_*: each idempotent, gene vectors / params embedded. Outputs in the cache
table. Human cross-species lane (completed plan): `build_human_microglia.py` (H1),
`build_human_validation_signatures.R` (H2), `build_human_substate_conservation.py`
(H3) run via the project `.venv` (scanpy) and read `storage/data/seaad/microglia_
{mtg,dlpfc}.h5ad`; the Python H3 script reaches the R-built signatures through a
jsonlite JSON sidecar. Causal lane (§20): `build_causal_network.R` runs
per-contrast CARNIVAL ILPs (solver **cbc @ threads=1** — threads>1 segfaults; a
per-contrast **L=3 reachability prune** `restrict_pkn_to_reachable` makes the
microglia-PKN ILP tractable) → `causal_network.rds`; `build_tf_activity_decoupler.R
--split-complexes` writes the `tf_activity_decoupler_split.rds` variant it consumes
(split is a no-op on mouse CollecTRI — byte-identical to the non-split cache).
SCENIC lane (§21, arc K): runs in a project-local **micromamba** env `scenic`
(`.micromamba/envs/scenic`, py3.10 + pyscenic 0.12.1, setuptools<81; spec
`scripts/scenic_env.yml`; both gitignored) against mm10 v10 cisTarget resources in
`storage/data/cistarget/` (rankings feather + motif2tf tbl + allTFs_mm.txt;
gitignored). Chain: `export_microglia_for_scenic.R` (Ensembl→MGI raw-count bundle,
≥1% gene filter → 11,536 genes / 994 candidate TFs) → `build_scenic_grn.py`
(10-seed GRNBoost2 via the no-dask `arboreto_with_multiprocessing` helper —
DENSE-not-`--sparse`, chunksize-patched; resumable) → `build_scenic_ctx_aucell.py`
(per-seed `pyscenic ctx` motif-pruning → **≥8/10 edge-recurrence consensus = 51
activating regulons** + AUCell; self-gating on the GRN pid) →
`build_scenic_contrasts.R` (sources R/scenic.R; head-to-head vs §14 + AUCell 2×2
factorial + per-substate → scenic_summary.rds + 5 results TSVs). All compute is
out-of-knit; rmd/20 is display-only.
Spatial lane (§22, arc L): `build_spatial_deconvolution.R` (sources
R/spatial_decon.R) reads `seurat_full_processed.rds` + `geomx.rds`, builds an
snRNAseq reference profile matrix, runs two-stage SpatialDecon (broad 6-level
then 4-substate microglia), fits the 5-contrast factorial on log-abundance +
per-slide Moran's I → `spatial_decon.rds` + 3 results TSVs; rmd/21 is
display-only and writes `spatial_decon_verdict.tsv` at knit time.
`build_summary_rmd.R` assembles
the standalone `summary.Rmd` from curated panel chunks; its §7 cross-species
panels read `summary_human_validation.rds` (mechanism heatmap) and
`results/human_substate_conservation_metrics.tsv` (conservation gate table).
Capstone lane (arc P) — display-only in `summary.Rmd`, `analysis.Rmd` §17
untouched: `build_capstone_synthesis.R` re-aggregates READ-ONLY the ledger +
adjudication + contest_verdicts + the 8 per-arc verdict TSVs + pathway
leaderboard + `integration_table.rds`/`summary_human_validation.rds` →
`results/capstone_{convergence_matrix,contest_summary}.tsv` (39×9 long matrix +
3×10 contest summary; a build-time guard hard-asserts the sealed margins 18/12/55
so the capstone fails loudly if §17 drifts). `build_capstone_figure.py`
(matplotlib in `.venv`) renders those two TSVs → `storage/figures/
capstone_convergence.png` (committed deliverable, tracked via a `.gitignore`
re-inclusion). The gitignored `summary_chunks/capstone_convergence.R` panel ("The
integrated convergence model", final summary section) embeds the PNG
(`include_graphics`) + a kable of the contest summary. Both builders run
out-of-knit — refresh the figure/TSVs by re-running them, then `build_summary_rmd.R`.
Synthesis holistic-models lane — display-only in `synthesis.Rmd`, locked analysis
untouched: `build_synthesis_models.R` writes the authored
`results/synthesis_model_{comparison,discriminating_experiments}.tsv` (lead +
4-rival systems-model comparison + discriminating experiments; ASCII-safe cells,
read by `synthesis.Rmd` kables); `build_synthesis_model_figure.py` (matplotlib in
`.venv`) renders the lead-model wiring schematic → `storage/figures/
synthesis_model.png` (self-contained authored topology, tracked via the figures
`.gitignore` re-inclusion). Both authored/out-of-knit (curated from a dynamic
multi-agent modelling workflow); refresh by re-running, then re-knit `synthesis.Rmd`.

## Where to look next

- Top-level reports (3, each knits standalone): `analysis.Rmd` (full pipeline,
  child rmds in `rmd/`, §-numbered) · `summary.Rmd` (curated figure gallery,
  reads light `rd()` caches + `rt()` TSVs) · `synthesis.Rmd` (read-only,
  literature-grounded cross-arc conclusion in prose; reads only `results/*.tsv` +
  embeds the capstone PNG **and the holistic-models schematic `synthesis_model.png`**;
  surfaces the sealed §17 contests 18/12/55 read-only, modifies nothing; its
  **"Candidate holistic models" section** leads with one unified systems model + 4
  adversarially-vetted rivals + a discriminating-experiments table, reading
  `synthesis_model_{comparison,discriminating_experiments}.tsv`). `synthesis.{html,
  _files,_cache}` are gitignored + deny-Read.
- Standing contract: `.agent/memory.md` (rules + locked decisions + gotchas; boot-loaded).
- Forward roadmap: `.agent/roadmap.md` (posture + backlog + active-plan pointer; read first on launch).
- Active work: `.agent/*_plan.md` (read the STATUS block, then only the next TODO step).
- Past decisions: `.agent/history.md` (read this, NOT the full plans; the archived full `.agent/completed/*_plan_*.md` are deny-Read in .claude/settings.json).
- House rules: `CLAUDE.md`; session boot: `/session-prompt` command (`.claude/commands/session-prompt.md`).

## Refresh hints (if this map drifts)

- child order: `grep -n 'child *=' analysis.Rmd`
- globals/caches per rmd: `grep -nE '<- *cache_or_run|<- *readRDS|<<?- ' rmd/*.Rmd`
- R symbols: `grep -nE '^[a-zA-Z_.]+ *<- *function' R/*.R`
- script outputs: `grep -hoE '"[^"]*\.(rds|tsv)"' scripts/build_*.R`

#!/usr/bin/env Rscript
# build_biological_model_ledger.R
#
# H2 deliverable: emit storage/results/biological_model_claims_ledger.tsv
# (the structured evidence base for the H3 model adjudication and the
# centrepiece of the H4 section-17 chapter).
#
# Schema (14 columns, locked at H1 2026-05-25): claim_id, axis, layer,
# claim, direction, effect_size, primary_evidence_source,
# corroborating_evidence, n_replicates_in_modalities,
# n_replicates_in_layers, confidence_grade, supports_models,
# contradicts_models, notes.
#
# Eleven adjudicated entities (locked H1): Hyp-1A, Hyp-1B (axis 1
# amyloid_activation contest); Hyp-2A, Hyp-2B (axis 2 synaptic_suppression
# contest); Hyp-3A, Hyp-3B (axis 3 interaction_metabolic contest); Hyp-0 (Cdk5
# cross-axis integrator, stand-alone); T-Inflammation,
# T-Compartment-suppression, T-Tau-attenuates, T-Synergy (cross-axis
# themes, stand-alone).
#
# Source-priority order (locked H1): verdict TSVs -> axis-restricted
# TSVs -> cross-tool / cross-modality leader boards -> integration /
# concordance -> pathway / module support -> multi-tool LR detail ->
# completed-plan outcome summaries.

suppressPackageStartupMessages(library(utils))

# ---- helpers --------------------------------------------------------------

# Single-row constructor. All arguments are required; quote each effect_size
# as a string to preserve exact formatting (rounded numerics drift the TSV).
row <- function(claim_id, axis, layer, claim, direction, effect_size,
                primary_evidence_source, corroborating_evidence,
                n_replicates_in_modalities, n_replicates_in_layers,
                confidence_grade, supports_models, contradicts_models,
                notes) {
  data.frame(
    claim_id                    = claim_id,
    axis                        = axis,
    layer                       = layer,
    claim                       = claim,
    direction                   = direction,
    effect_size                 = effect_size,
    primary_evidence_source     = primary_evidence_source,
    corroborating_evidence      = corroborating_evidence,
    n_replicates_in_modalities  = n_replicates_in_modalities,
    n_replicates_in_layers      = n_replicates_in_layers,
    confidence_grade            = confidence_grade,
    supports_models             = supports_models,
    contradicts_models          = contradicts_models,
    notes                       = notes,
    stringsAsFactors            = FALSE
  )
}

# Locked 11-entity set; any IDs in supports_models / contradicts_models
# must be a subset of this list.
ENTITY_IDS <- c(
  "Hyp-1A", "Hyp-1B", "Hyp-2A", "Hyp-2B", "Hyp-3A", "Hyp-3B", "Hyp-0",
  "T-Inflammation", "T-Compartment-suppression",
  "T-Tau-attenuates", "T-Synergy"
)

# ---- AXIS 1 (amyloid_activation) ------------------------------------------

axis1 <- rbind(
  # --- TF layer (7 rows) ---
  row("H2-001", "amyloid_activation", "tf",
      "Spi1 (PU.1) is activated at NLGF onset in both tau backgrounds; rank 2 in axis-restricted top 5 with cross-modality leader status at nlgf_in_maptki.",
      "+", "mean +2.81; score_at_nlgf_in_maptki +2.64; score_at_nlgf_in_p301s +2.98; 184 axis-1 targets",
      "tf_activity_axis_restricted.tsv (axis=amyloid_activation, source=Spi1, leader_rank=2)",
      "tf_activity_unified_leaderboard.tsv (Spi1 leader_score 5.4 at nlgf_in_maptki); tf_activity_verdict.tsv axis 1 evidence_summary",
      2L, 4L, "Moderate",
      "Hyp-1A;Hyp-1B;T-Inflammation", "",
      "Spi1 is the master microglial lineage TF + AD GWAS locus; activation at both NLGF arms is consistent with both Hyp-1A (tau-independent) and Hyp-1B (NF-kB modulation occurs downstream); Interpretation A entanglement: this row inherits identical mean-score values to the axis-2 Spi1 row (in-universe target counts 184 vs 84 differentiate)."),
  row("H2-002", "amyloid_activation", "tf",
      "Nfkb1 (NF-kB p50 subunit) is activated at NLGF onset in both tau backgrounds; rank 5 in axis-restricted top 5.",
      "+", "mean +2.47; score_at_nlgf_in_maptki +2.70; score_at_nlgf_in_p301s +2.24; 226 axis-1 targets",
      "tf_activity_axis_restricted.tsv (axis=amyloid_activation, source=Nfkb1, leader_rank=5)",
      "tf_activity_verdict.tsv axis 1 evidence_summary; mirrors the kinase-layer Mapk14 and Mapk8 stress-kinase ensemble",
      1L, 2L, "Moderate",
      "Hyp-1A;Hyp-1B;T-Inflammation", "",
      "NF-kB family canonical inflammatory amplifier; the axis-mean POSITIVE direction here is what Hyp-1B claims is REVERSED at the interaction contrast (see cross-axis Rela sign-reversal H2-072); both Hyp-1A and Hyp-1B accept Nfkb1 activation at the NLGF arms, they differ on the interaction contrast."),
  row("H2-003", "amyloid_activation", "tf",
      "Sp3 is activated at NLGF onset in both tau backgrounds; rank 4 in axis-restricted top 5.",
      "+", "mean +2.50; score_at_nlgf_in_maptki +2.68; score_at_nlgf_in_p301s +2.31; 161 axis-1 targets",
      "tf_activity_axis_restricted.tsv (axis=amyloid_activation, source=Sp3, leader_rank=4)",
      "tf_activity_verdict.tsv axis 1 evidence_summary",
      1L, 1L, "Moderate",
      "Hyp-1A;Hyp-1B;T-Inflammation", "",
      "Sp3 is a canonical inflammatory amplifier; supports both per-axis models at axis 1 by activating consistently in both backgrounds."),
  row("H2-004", "amyloid_activation", "tf",
      "Rreb1 is the lone repressor in the axis-1 TF top 5, suppressed at NLGF onset in both tau backgrounds.",
      "-", "mean -2.57; score_at_nlgf_in_maptki -2.09; score_at_nlgf_in_p301s -3.06; 19 axis-1 targets",
      "tf_activity_axis_restricted.tsv (axis=amyloid_activation, source=Rreb1, leader_rank=3)",
      "tf_activity_verdict.tsv axis 1 evidence_summary",
      1L, 1L, "Suggestive",
      "T-Inflammation", "",
      "n_cells_used=4 (rather than the typical 6); n_targets_in_axis_universe=19 (small footprint); the strongest p301s side -3.06 is the wider negative tail noted in the verdict per-contrast asymmetry."),
  row("H2-005", "amyloid_activation", "tf",
      "A0A979HLR9 (CollecTRI complex source proxying a TF activator family) shows the strongest TF signal at axis 1 and reaches cross-modality leader status at nlgf_in_maptki.",
      "+", "mean +4.32; score_at_nlgf_in_maptki +5.40; score_at_nlgf_in_p301s +3.23; 442 axis-1 targets",
      "tf_activity_axis_restricted.tsv (axis=amyloid_activation, source=A0A979HLR9, leader_rank=1)",
      "tf_activity_unified_leaderboard.tsv (A0A979HLR9 leader_score 5.4 at nlgf_in_maptki)",
      2L, 1L, "Moderate",
      "Hyp-1A;Hyp-1B;T-Inflammation", "",
      "Cross-modality leader rule met at one NLGF arm; the maptki side carries the wider positive tail (+5.40 vs +3.23) as noted in the verdict per-contrast asymmetry; complex-source identifier hampers direct biological interpretation but contributes to the axis-1 activator weight."),
  row("H2-006", "amyloid_activation", "tf",
      "Rela (NF-kB p65 subunit) is activated at NLGF onset in both tau backgrounds; rank 7 in axis-restricted top 10.",
      "+", "mean +2.19; score_at_nlgf_in_maptki +2.57; score_at_nlgf_in_p301s +1.81; 266 axis-1 targets",
      "tf_activity_axis_restricted.tsv (axis=amyloid_activation, source=Rela, leader_rank=7)",
      "tf_activity_axis_restricted.tsv axis 3 Rela row (the sign-reversal companion at interaction); Nfkb1 H2-002 covariant",
      1L, 2L, "Moderate",
      "Hyp-1A;Hyp-1B;T-Inflammation", "",
      "At axis 1 alone, Rela activation is symmetrical with respect to tau background and supports both Hyp-1A and Hyp-1B; the Hyp-1B-specific evidence is the SIGN-REVERSAL at the interaction contrast (H2-072); this row is the positive half of the sign-reversal."),
  row("H2-007", "amyloid_activation", "tf",
      "A maptki-vs-p301s per-contrast asymmetry at the top-5 axis-1 TF view: nlgf_in_maptki [-2.09, +5.40] vs nlgf_in_p301s [-3.06, +3.23] -- maptki carries the wider positive activator tail, p301s the wider negative repressor tail.",
      "mixed", "per_contrast_score_range across top 5: maptki [-2.09, +5.40]; p301s [-3.06, +3.23]",
      "tf_activity_verdict.tsv axis 1 per_contrast_summary",
      "kinase axis-1 mirror per_contrast_summary nlgf_in_maptki [-1.32, +3.02] vs nlgf_in_p301s [-1.00, +3.00] (different in scale, similar in asymmetry direction)",
      1L, 2L, "Suggestive",
      "Hyp-1B", "Hyp-1A",
      "Asymmetric reading WEAKLY supports Hyp-1B (tau modifies amyloid program); the axis-mean smooths this asymmetry over but the per-contrast view recovers it; the asymmetry is qualitatively informative but does not exceed any FDR threshold at this layer."),

  # --- kinase layer (5 rows) ---
  row("H2-008", "amyloid_activation", "kinase",
      "Cdk2 (cell-cycle progression / microglial proliferative re-entry) is activated at NLGF onset in both tau backgrounds; rank 1 in axis-restricted top 5.",
      "+", "mean +2.91; score_at_nlgf_in_maptki +3.02; score_at_nlgf_in_p301s +2.80; 259 axis-1 substrate sites",
      "kinase_activity_axis_restricted.tsv (axis=amyloid_activation, source=Cdk2, leader_rank=1)",
      "kinase_activity_verdict.tsv axis 1 evidence_summary; TF-layer Sp3 / Spi1 co-activation",
      1L, 2L, "Moderate",
      "Hyp-1A;Hyp-1B;T-Inflammation", "",
      "Cell-cycle re-entry kinase, signature of microglial proliferation in DAM-state amplification; supports both Hyp-1A and Hyp-1B by activating in both backgrounds."),
  row("H2-009", "amyloid_activation", "kinase",
      "Cdk5 (canonical tau kinase) is activated at NLGF onset with a p301s-biased asymmetry; rank 2 in axis-restricted top 5.",
      "+", "mean +1.75; score_at_nlgf_in_maptki +0.50; score_at_nlgf_in_p301s +3.00 (padj_ulm 0.0077); 81 axis-1 substrate sites",
      "kinase_activity_axis_restricted.tsv (axis=amyloid_activation, source=Cdk5, leader_rank=2); kinase_activity_per_contrast.tsv (Cdk5 nlgf_in_p301s padj_ulm 0.0077)",
      "Cross-axis Cdk5 at axes 2 + 3 (Hyp-0 evidence); kinase verdict axis 1 evidence_summary",
      1L, 2L, "Moderate",
      "Hyp-1A;Hyp-1B;Hyp-0;T-Inflammation", "",
      "p301s-biased activation (maptki +0.50 vs p301s +3.00) consistent with the axis-mean per-contrast asymmetry; FDR<0.10 reached only at nlgf_in_p301s, not nlgf_in_maptki; this row also contributes Strong evidence to the Hyp-0 cross-axis Cdk5 integrator claim (see H2-071)."),
  row("H2-010", "amyloid_activation", "kinase",
      "Csnk1e (CK1-epsilon stress-activated kinase) is activated at NLGF onset; rank 3 in axis-restricted top 5; reaches FDR<0.10 at nlgf_in_p301s.",
      "+", "mean +1.64; score_at_nlgf_in_maptki +1.19; score_at_nlgf_in_p301s +2.08 (padj_ulm < 0.10); 38 axis-1 substrate sites",
      "kinase_activity_axis_restricted.tsv (axis=amyloid_activation, source=Csnk1e, leader_rank=3)",
      "kinase_activity_verdict.tsv axis 1 evidence_summary (one of two top-5 reaching FDR<0.10 at nlgf_in_p301s)",
      1L, 1L, "Moderate",
      "Hyp-1A;Hyp-1B;T-Inflammation", "",
      "Casein kinase 1 epsilon, stress-activated, p301s-significant on the per-contrast ulm test."),
  row("H2-011", "amyloid_activation", "kinase",
      "Mapk8 (JNK stress-activated kinase) is activated at NLGF onset with a p301s-biased asymmetry; rank 4 in axis-restricted top 5.",
      "+", "mean +1.17; score_at_nlgf_in_maptki +0.25; score_at_nlgf_in_p301s +2.08; 176 axis-1 substrate sites",
      "kinase_activity_axis_restricted.tsv (axis=amyloid_activation, source=Mapk8, leader_rank=4)",
      "kinase_activity_verdict.tsv axis 1 evidence_summary",
      1L, 1L, "Suggestive",
      "Hyp-1A;Hyp-1B;T-Inflammation", "",
      "p301s-biased; would not have ranked here on maptki alone (per the verdict per-contrast asymmetry note)."),
  row("H2-012", "amyloid_activation", "kinase",
      "Camk2g (CaMKII-gamma calcium/calmodulin signalling) is the lone repressor in the axis-1 kinase top 5, suppressed at NLGF onset.",
      "-", "mean -1.16; score_at_nlgf_in_maptki -1.32; score_at_nlgf_in_p301s -1.00; 24 axis-1 substrate sites",
      "kinase_activity_axis_restricted.tsv (axis=amyloid_activation, source=Camk2g, leader_rank=5)",
      "kinase_activity_verdict.tsv axis 1 evidence_summary",
      1L, 1L, "Suggestive",
      "T-Inflammation", "",
      "n_targets_in_axis_universe=24 small footprint; calcium signalling attenuation under amyloid; consistent in both tau backgrounds."),

  # --- LR layer (10 rows) ---
  row("H2-013", "amyloid_activation", "lr",
      "Top-5 axis-1 LR pairs are dominated by Neuronal-sender axon-guidance / adhesion biology (Adam11_Itga4, Gas6_Axl, L1cam_Egfr, Omg_Rtn4r, Rspo3_Lgr4) rather than canonical DAM-program microglia-involving pairs.",
      "+", "top-5 mean_activity 0.98-1.00; all positive sign at NLGF arms",
      "ccc_lr_axis_restricted.tsv (axis=amyloid_activation, leader_rank 1-5)",
      "ccc_lr_verdict.tsv axis 1 evidence_summary (the broader 55/100 microglia-involving figure is the joint LR reading)",
      1L, 1L, "Moderate",
      "", "",
      "n_cells_used=1 for several top-5 rows (single-contrast hits dominate by axis-magnitude ranking); the top-5 SHAPE supports neither Hyp-1A nor Hyp-1B specifically; treated as a context claim rather than a model-discriminating one."),
  row("H2-014", "amyloid_activation", "lr",
      "Cd200_Cd200r1 (Neuronal->Microglia_proliferative) is a microglia-involving DAM-program LR pair at axis 1; rank 6.",
      "+", "mean_activity 0.977; leader_rank 6; n_cells_used 1",
      "ccc_lr_axis_restricted.tsv (axis=amyloid_activation, ligand=Cd200, receptor=Cd200r1, leader_rank=6)",
      "ccc_lr_verdict.tsv axis 1 evidence_summary",
      1L, 1L, "Moderate",
      "Hyp-1A;Hyp-1B;T-Inflammation", "",
      "Canonical homeostatic-to-DAM neuronal modulation; supports both Hyp-1A and Hyp-1B by appearing in both NLGF arms."),
  row("H2-015", "amyloid_activation", "lr",
      "Apoe_Trem2 (Microglia_DAM->Microglia_DAM) appears at axis 1 rank 9 as a canonical DAM autocrine signal.",
      "+", "mean_activity 0.961; leader_rank 9; n_cells_used 2",
      "ccc_lr_axis_restricted.tsv (axis=amyloid_activation, ligand=Apoe, receptor=Trem2, sender=Microglia_DAM, receiver=Microglia_DAM, leader_rank=9)",
      "ccc_lr_verdict.tsv axis 1 evidence_summary; axis 2 Apoe_Trem2 rows H2-035/H2-036",
      1L, 1L, "Moderate",
      "Hyp-1A;Hyp-1B;T-Inflammation", "",
      "DAM autocrine signature at NLGF onset; the same pair appears at axis 2 ranks 9-10 and there contradicts Hyp-2A while supporting Hyp-2B."),
  row("H2-016", "amyloid_activation", "lr",
      "Apoe_Trem2 (Astrocyte->Microglia_DAM) appears at axis 1 rank 10 as a cross-cell-type TREM2-mediated signal.",
      "+", "mean_activity 0.959; leader_rank 10; n_cells_used 2",
      "ccc_lr_axis_restricted.tsv (axis=amyloid_activation, ligand=Apoe, receptor=Trem2, sender=Astrocyte, receiver=Microglia_DAM, leader_rank=10)",
      "ccc_lr_verdict.tsv axis 1 evidence_summary",
      1L, 1L, "Moderate",
      "Hyp-1A;Hyp-1B;T-Inflammation", "",
      "Astrocyte-derived ApoE signalling to microglial Trem2; supports both per-axis models at axis 1; foreshadows the axis-2 cross-cell-type signal."),
  row("H2-017", "amyloid_activation", "lr",
      "App_Cd74 (Oligodendrocyte->Microglia_proliferative) at axis 1 rank 11 brings APP-fragment biology to the amyloid axis.",
      "+", "mean_activity 0.956; leader_rank 11; n_cells_used 1",
      "ccc_lr_axis_restricted.tsv (axis=amyloid_activation, ligand=App, receptor=Cd74, leader_rank=11)",
      "ccc_lr_verdict.tsv axis 1 evidence_summary",
      1L, 1L, "Moderate",
      "Hyp-1A;Hyp-1B;T-Inflammation", "",
      "App fragment uptake into microglia; supports both per-axis models at axis 1; foreshadows the axis-2 App_Cd74 rows."),
  row("H2-018", "amyloid_activation", "lr",
      "Tgfb1_Sdc2 (Microglia_DAM->Neuronal) at axis 1 rank 34 is a longer-range microglia-to-neuron TGF-beta signal.",
      "+", "leader_rank 34",
      "ccc_lr_axis_restricted.tsv (axis=amyloid_activation, ligand=Tgfb1, receptor=Sdc2, leader_rank=34)",
      "ccc_lr_verdict.tsv axis 1 evidence_summary",
      1L, 1L, "Suggestive",
      "T-Inflammation", "",
      "TGF-beta is a canonical microglial-state remodeller; rank 34 is below the top-10 but cited in the verdict prose as part of the joint reading."),
  row("H2-019", "amyloid_activation", "lr",
      "55 of the top 100 axis-1 LR pairs involve a microglia sender or receiver -- the broader microglia-involving rate that the top-5 magnitude view masks.",
      "+", "55/100 microglia-involving in axis-1 top 100",
      "ccc_lr_axis_restricted.tsv (axis=amyloid_activation, top 100 by leader_rank)",
      "ccc_lr_verdict.tsv axis 1 evidence_summary",
      1L, 1L, "Strong",
      "Hyp-1A;Hyp-1B;T-Inflammation", "",
      "Microglia-involving rate is substantial across the broader top 100; this is the broader-scope reading that supports both per-axis models and is one of the strongest cross-axis context claims."),
  row("H2-020", "amyloid_activation", "lr",
      "6,074 axis-1 LR cells in the universe-filtered set; axis 1 LR universe is the second-largest of the three axes (axis 2 = 6,459; axis 3 = 4,139).",
      "+", "6,074 axis-1 LR cells",
      "ccc_lr_axis_restricted.tsv (axis=amyloid_activation, total cells)",
      "ccc_lr_verdict.tsv axis 1 evidence_summary",
      1L, 1L, "Suggestive",
      "", "",
      "Universe-size context claim; not directly model-discriminating but useful for the H4 chapter's narrative."),
  row("H2-021", "amyloid_activation", "lr",
      "Multinichenet top-100 LR pairs involve Microglia_DAM as sender in 95/100 pairs across contrasts (cumulative count) -- the predominantly DAM-state sender signature.",
      "+", "Microglia_DAM total cumulative sender count = 95 across 5 contrasts",
      "ccc_multinichenet_top100_sender_per_contrast.tsv (Microglia_DAM total = 95)",
      "ccc_lr_verdict.tsv axis 1 evidence_summary",
      4L, 1L, "Moderate",
      "Hyp-1A;Hyp-1B;T-Inflammation", "",
      "Cross-contrast aggregate (not per-axis); supports both per-axis models by establishing DAM-sender dominance regardless of tau background."),
  row("H2-022", "amyloid_activation", "lr",
      "Multinichenet top-100 LR pairs involve Microglia_DAM as receiver in 134/100 pairs across contrasts (cumulative count, multiple-counting across contrasts) -- the predominantly DAM-state receiver signature.",
      "+", "Microglia_DAM total cumulative receiver count = 134 across 5 contrasts",
      "ccc_multinichenet_top100_receiver_per_contrast.tsv (Microglia_DAM total = 134)",
      "ccc_lr_verdict.tsv axis 1 evidence_summary",
      4L, 1L, "Moderate",
      "Hyp-1A;Hyp-1B;T-Inflammation", "",
      "DAM as both top sender and top receiver across contrasts; reinforces DAM-state centrality in the microglial response."),

  # --- pathway / module layer (7 rows) ---
  row("H2-023", "amyloid_activation", "pathway",
      "MG-M2 hdWGCNA module (DAM-state marker module: Apoe;Cst7;Lpl;Spp1;Tyrobp; substate_DAM enrichment padj 1.8e-3) is up-regulated at both NLGF arms with cross-modality consistent-sign breadth 4.",
      "+", "MG-M2 logFC +5.00 at nlgf_in_maptki (padj 2.2e-4); +6.21 at nlgf_in_p301s (padj 1.8e-5); leader_score 12.8; max_consistent_sign 4",
      "hdwgcna_module_de.tsv (module=MG-M2, contrasts nlgf_in_maptki + nlgf_in_p301s)",
      "pathway_survey_unified_leaderboard.tsv (MG-M2 leader_score 12.8); hdwgcna_module_enrichment.tsv (MG-M2 substate_DAM odds_ratio 6.02 padj 1.75e-3)",
      4L, 2L, "Strong",
      "Hyp-1A;Hyp-1B;T-Inflammation", "",
      "Cross-modality breadth 4 + cross-layer corroboration via substate_DAM enrichment + canonical DAM hub genes; STRONG by H1 rules; supports both per-axis models at axis 1."),
  row("H2-024", "amyloid_activation", "pathway",
      "GOBP_ADAPTIVE_IMMUNE_RESPONSE is up-regulated at both NLGF arms with breadth-3 consistent sign at nlgf_in_p301s.",
      "+", "leader_score 11.6; max_consistent_sign 3; dominant_sign +",
      "pathway_survey_unified_leaderboard.tsv (pathway=GOBP_ADAPTIVE_IMMUNE_RESPONSE)",
      "ccc_lr_verdict.tsv axis 1 evidence_summary (broader DAM-program rationale); tf-layer Spi1 / Nfkb1 / Rela ensemble",
      3L, 3L, "Strong",
      "Hyp-1A;Hyp-1B;T-Inflammation", "",
      "Cross-modality + cross-layer corroboration with the inflammatory TF ensemble; STRONG by H1 rules."),
  row("H2-025", "amyloid_activation", "pathway",
      "MHC II antigen-presentation family (GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_OR_POLYSACCHARIDE_ANTIGEN_VIA_MHC_CLASS_II + related GOBP / GOMF terms) is up-regulated at both NLGF arms.",
      "+", "leader_score 10.4 (each); dominant_sign +; max_n_modalities_sig 2",
      "pathway_survey_unified_leaderboard.tsv (MHC II / antigen-presentation family rows)",
      "GOBP_ADAPTIVE_IMMUNE_RESPONSE H2-024 covariant; TF-layer Spi1 / Nfkb1 / Rela covariant",
      2L, 3L, "Strong",
      "Hyp-1A;Hyp-1B;T-Inflammation", "",
      "Antigen-presentation family consistently up at both NLGF arms; canonical DAM-program output."),
  row("H2-026", "amyloid_activation", "pathway",
      "Custom microglial-state gene set DAM_up is up-regulated at both NLGF arms with breadth-3 consistent sign and cross-modality breadth.",
      "+", "leader_score 16.4; max_consistent_sign 2; max_n_modalities_sig 3; dominant_sign mixed (interaction:3/mixed shifts dominant)",
      "pathway_survey_unified_leaderboard.tsv (collection=custom_microglia_states, pathway=DAM_up)",
      "MG-M2 module H2-023 covariant (DAM-state marker module overlap)",
      3L, 2L, "Strong",
      "Hyp-1A;Hyp-1B;T-Inflammation", "",
      "Direct DAM-state gene set; high leader_score; STRONG; supports both per-axis models at axis 1."),
  row("H2-027", "amyloid_activation", "pathway",
      "Custom microglial-AD gene sets WAM (microglia near amyloid plaques) and HAM_disease are up-regulated at both NLGF arms.",
      "+", "WAM + HAM_disease leader_score 10.4 each; max_consistent_sign 2; dominant_sign +",
      "pathway_survey_unified_leaderboard.tsv (collection=custom_microglia_states/ad, pathways WAM + HAM_disease)",
      "DAM_up H2-026 + MG-M2 H2-023 covariant",
      2L, 2L, "Strong",
      "Hyp-1A;Hyp-1B;T-Inflammation", "",
      "Disease-associated microglial signatures up at both NLGF arms regardless of tau background."),
  row("H2-028", "amyloid_activation", "pathway",
      "P301S NLGF response is approximately 8x more asymmetric than MAPTKI NLGF response (446 P301S-specific genes vs 56 MAPTKI-specific genes; 82 shared) on the ifn_nlgf_asymmetry analysis -- the per-contrast asymmetric reading at the gene-set layer.",
      "+", "P301S-specific 446 vs MAPTKI-specific 56 vs shared 82",
      "ifn_nlgf_asymmetry.tsv (category counts)",
      "Per-contrast asymmetry mirror at TF axis 1 (H2-007) and kinase axis 1 (Cdk5 / Mapk8 p301s-biased)",
      1L, 1L, "Moderate",
      "Hyp-1B;T-Tau-attenuates", "Hyp-1A",
      "The 8x asymmetry directly contradicts Hyp-1A's prediction of qualitatively identical amyloid response across tau backgrounds; supports Hyp-1B and the cross-axis theme T-Tau-attenuates.")
)

# ---- AXIS 2 (synaptic_suppression) -----------------------------------------

axis2 <- rbind(
  # --- TF layer (3 rows: honest non-finding + entanglement + in-universe differential) ---
  row("H2-029", "synaptic_suppression", "tf",
      "Honest non-finding: the axis-2 axis-restricted TF top 5 is IDENTICAL to axis 1 by Interpretation A construction; no axis-2-distinct TF passes FDR with cross-modality replication.",
      "mixed", "axis-2 top 5 = axis-1 top 5 (A0A979HLR9, Spi1, Rreb1, Sp3, Nfkb1); mean range [-2.57, +4.32] identical",
      "tf_activity_axis_restricted.tsv (axis=synaptic_suppression); tf_activity_verdict.tsv axis 2 evidence_summary",
      "kinase axis-2 mirror H2-032 (same non-finding shape)",
      1L, 2L, "Suggestive",
      "T-Compartment-suppression", "",
      "Interpretation A entanglement at axis 2: both axes share NLGF contrasts so per-source means are identical; the TF layer does NOT independently identify drivers of synaptic suppression at this resolution; this is the engulfment-mediated mechanism hypothesis raised in 14.1.4 / 14.2.4 (the regulator-driven layers fail at axis 2 because the mechanism is post-regulator)."),
  row("H2-030", "synaptic_suppression", "tf",
      "In-universe TF target counts differ between axes (Spi1: 184 amyloid targets vs 84 synaptic targets); axis-2 TFs read as amyloid-activated TFs whose target programmes also touch synaptic genes, not as TFs that drive synaptic gene downregulation.",
      "mixed", "Spi1 in-universe target counts: amyloid 184 vs synaptic 84 (2.2x more in amyloid); analogous reductions for Nfkb1 and Sp3",
      "tf_activity_axis_restricted.tsv (axis=synaptic_suppression, n_targets_in_axis_universe column)",
      "tf_activity_verdict.tsv axis 2 evidence_summary",
      1L, 1L, "Suggestive",
      "T-Compartment-suppression", "",
      "Target-count differential is the structural differentiator that the axis-mean smooths over; supports the engulfment-mechanism reading at axis 2 by negative evidence for regulator-driven mechanisms."),
  row("H2-031", "synaptic_suppression", "tf",
      "Axis-2 TF non-finding is consistent with the engulfment-mediated synaptic suppression hypothesis: if NLGF suppresses the synaptic compartment via intercellular signalling (TREM2 / APP / classical complement) rather than via transcriptional regulator activation, then the TF layer should NOT identify axis-2-distinct drivers.",
      "mixed", "0 axis-2-distinct TFs in top 5; LR layer axis-2 has 6,459 cells (largest axis universe) by contrast",
      "tf_activity_verdict.tsv axis 2 evidence_summary (the explicit hypothesis-aligned interpretation)",
      "kinase axis-2 H2-032 (same hypothesis-aligned non-finding); LR axis-2 H2-035 + H2-036 (Apoe_Trem2 supports the alternative mechanism)",
      0L, 2L, "Moderate",
      "T-Compartment-suppression", "",
      "Hypothesis-aligned negative-evidence claim; supports the cross-axis theme T-Compartment-suppression but does not discriminate Hyp-2A vs Hyp-2B by itself (the LR layer is where Hyp-2A vs Hyp-2B is adjudicated)."),

  # --- kinase layer (3 rows: honest non-finding + entanglement + Cdk5 KSN footprint expansion) ---
  row("H2-032", "synaptic_suppression", "kinase",
      "Honest non-finding: the axis-2 axis-restricted kinase top 5 is IDENTICAL to axis 1 by Interpretation A construction; no axis-2-distinct kinase passes FDR with cross-modality replication.",
      "mixed", "axis-2 top 5 = axis-1 top 5 (Cdk2, Cdk5, Csnk1e, Mapk8, Camk2g); mean range [-1.16, +2.91] identical",
      "kinase_activity_axis_restricted.tsv (axis=synaptic_suppression); kinase_activity_verdict.tsv axis 2 evidence_summary",
      "TF axis-2 mirror H2-029 (same non-finding shape)",
      1L, 2L, "Suggestive",
      "T-Compartment-suppression", "",
      "Same Interpretation A entanglement as TF axis 2; the kinase layer does NOT independently identify drivers of synaptic suppression."),
  row("H2-033", "synaptic_suppression", "kinase",
      "Cdk5 KSN substrate footprint expands from 81 sites in the amyloid universe to 263 sites in the synaptic universe (3.2x); analogous expansions for Cdk2 (259->293), Csnk1e (38->70), and Camk2g (24->49).",
      "+", "Cdk5: 81 amyloid sites -> 263 synaptic sites (3.2x); 81 -> 153 interaction sites (1.9x)",
      "kinase_activity_axis_restricted.tsv (n_targets_in_axis_universe across axes for Cdk5)",
      "kinase_activity_verdict.tsv axis 2 evidence_summary; Hyp-0 cross-axis Cdk5 row H2-071",
      1L, 1L, "Strong",
      "Hyp-0;T-Compartment-suppression", "",
      "Cross-axis substrate-footprint scaling claim; STRONG (cross-axis corroboration of the Hyp-0 Cdk5 integrator); the synaptic-universe footprint is the LARGEST of the three axes for Cdk5, consistent with Cdk5 as a well-documented synaptic regulator whose substrate program preferentially touches synaptic machinery."),
  row("H2-034", "synaptic_suppression", "kinase",
      "Axis-2 kinase non-finding mirrors the TF axis-2 non-finding: the regulator-driven layers AGREE that the synaptic axis is not regulator-driven at this resolution.",
      "mixed", "0 axis-2-distinct kinases in top 5; 0 axis-2-distinct TFs in top 5",
      "kinase_activity_verdict.tsv axis 2 evidence_summary + tf_activity_verdict.tsv axis 2 evidence_summary",
      "LR axis-2 H2-035 / H2-036 (Apoe_Trem2 rows surfacing the post-regulator mechanism)",
      0L, 3L, "Moderate",
      "T-Compartment-suppression", "",
      "Hypothesis-aligned cross-layer agreement; the project's reading of axis 2 as engulfment / post-regulator mechanism rests on this paired TF + kinase non-finding combined with the LR-layer positive surface."),

  # --- LR layer (14 rows: Apoe_Trem2 + App_Cd74 + complement-absence + Pros1_Mertk nuance + microglia-involving + universe size) ---
  row("H2-035", "synaptic_suppression", "lr",
      "Apoe_Trem2 (Microglia_DAM->Microglia_DAM) appears at axis 2 rank 9 -- the highest-ranked microglia-involving DAM-program LR pair at the synaptic axis.",
      "+", "mean_activity 0.961; leader_rank 9; n_cells_used 2",
      "ccc_lr_axis_restricted.tsv (axis=synaptic_suppression, ligand=Apoe, receptor=Trem2, sender=Microglia_DAM, receiver=Microglia_DAM, leader_rank=9)",
      "ccc_lr_verdict.tsv axis 2 evidence_summary (G1 triangulation hypothesis vindicated)",
      1L, 1L, "Strong",
      "Hyp-2B;T-Compartment-suppression", "Hyp-2A",
      "TREM2-mediated DAM autocrine clearance signal at axis 2; SUPPORTS Hyp-2B (TREM2/APP/synaptic-adhesion mechanism); CONTRADICTS Hyp-2A (classical complement) on the axis-restricted view; STRONG by 2-modality + cross-layer corroboration with the TF + kinase non-findings (the LR layer surfaces what the regulator layers cannot)."),
  row("H2-036", "synaptic_suppression", "lr",
      "Apoe_Trem2 (Astrocyte->Microglia_DAM) appears at axis 2 rank 10 -- cross-cell-type TREM2-mediated signal at the synaptic axis.",
      "+", "mean_activity 0.959; leader_rank 10; n_cells_used 2",
      "ccc_lr_axis_restricted.tsv (axis=synaptic_suppression, ligand=Apoe, receptor=Trem2, sender=Astrocyte, receiver=Microglia_DAM, leader_rank=10)",
      "ccc_lr_verdict.tsv axis 2 evidence_summary",
      1L, 1L, "Strong",
      "Hyp-2B;T-Compartment-suppression", "Hyp-2A",
      "Astrocyte-derived ApoE signalling to microglial Trem2; supports Hyp-2B's TREM2-mediated clearance route; STRONG."),
  row("H2-037", "synaptic_suppression", "lr",
      "App_Cd74 (Oligodendrocyte->Microglia_proliferative) appears at axis 2 rank 12 -- the highest-ranked APP-fragment microglia-involving pair at the synaptic axis.",
      "+", "mean_activity 0.956; leader_rank 12; n_cells_used 1",
      "ccc_lr_axis_restricted.tsv (axis=synaptic_suppression, ligand=App, receptor=Cd74, sender=Oligodendrocyte, receiver=Microglia_proliferative, leader_rank=12)",
      "ccc_lr_verdict.tsv axis 2 evidence_summary",
      1L, 1L, "Moderate",
      "Hyp-2B;T-Compartment-suppression", "Hyp-2A",
      "APP fragment uptake by proliferating microglia; supports Hyp-2B's APP-fragment-mediated clearance route; n_cells_used=1 caveat."),
  row("H2-038", "synaptic_suppression", "lr",
      "App_Cd74 (Oligodendrocyte->Microglia_DAM) appears at axis 2 rank 20 -- APP-fragment signal to DAM microglia at the synaptic axis.",
      "+", "leader_rank 20",
      "ccc_lr_axis_restricted.tsv (axis=synaptic_suppression, ligand=App, receptor=Cd74, sender=Oligodendrocyte, receiver=Microglia_DAM, leader_rank=20)",
      "ccc_lr_verdict.tsv axis 2 evidence_summary",
      1L, 1L, "Moderate",
      "Hyp-2B;T-Compartment-suppression", "Hyp-2A",
      "APP fragment uptake by DAM microglia; second App_Cd74 row at axis 2 reinforces Hyp-2B."),
  row("H2-039", "synaptic_suppression", "lr",
      "Canonical complement engulfment LR pairs are DEEPLY ranked at axis 2: first complement pair C1qb_Lrp1 (Microglia_homeostatic->Neuronal) at rank 775; first regulator-of-engulfment Cd47_Sirpa at rank 279; no C3 / Mertk-Pros1 / Tyro3-Axl pairs in the top 800 axis-restricted cells.",
      "mixed", "rank 775 (C1qb_Lrp1); rank 279 (Cd47_Sirpa); 0 C3/Mertk-Pros1/Tyro3-Axl in top 800",
      "ccc_lr_axis_restricted.tsv (axis=synaptic_suppression, ligands C1qb / C3 / Mertk / Pros1 / Tyro3 / Axl rank columns)",
      "ccc_lr_verdict.tsv axis 2 evidence_summary (this is the load-bearing contradiction of Hyp-2A on axis-restricted scope)",
      1L, 1L, "Strong",
      "Hyp-2B;T-Compartment-suppression", "Hyp-2A",
      "Strongest direct contradiction of Hyp-2A's prediction at axis 2; the classical complement-tagging axis is NOT the dominant mechanism in these data on the axis-restricted scope."),
  row("H2-040", "synaptic_suppression", "lr",
      "Pros1_Mertk surfaces strongly in the CROSS-TOOL LR leader board (Microglia_IFN->multiple receivers; leader_score 15.4; dominant_sign +; max_n_tools_sig 2; contrasts: nlgf_in_p301s + interaction + tau_in_nlgf), even though absent from the axis-2 axis-restricted top 800.",
      "+", "leader_score 15.4 (Microglia_IFN sender); 6 distinct cross-tool rows at leader_score 15.4; contrasts_summary nlgf_in_p301s:2/+ | interaction:2/+ | tau_in_nlgf:2/+",
      "ccc_lr_cross_tool_leaderboard.tsv (Pros1_Mertk rows, leader_score 15.4)",
      "ccc_lr_axis_restricted.tsv (axis=synaptic_suppression) -- Pros1_Mertk absent from top 800 there",
      2L, 1L, "Moderate",
      "Hyp-2A", "",
      "PROTOCOL ROW per H1 H2 Open-questions: axis-restricted view contradicts Hyp-2A (H2-039); this CROSS-TOOL view partially supports Hyp-2A's classical-complement axis as Pros1-Mertk efferocytosis specifically when the scope widens to all-contrast leader status; the H4 17.3 paragraph names this scope ambiguity honestly rather than collapsing it."),
  row("H2-041", "synaptic_suppression", "lr",
      "49 of the top 100 axis-2 LR pairs involve a microglia sender or receiver -- the broader microglia-involving rate at the synaptic axis.",
      "+", "49/100 microglia-involving in axis-2 top 100",
      "ccc_lr_axis_restricted.tsv (axis=synaptic_suppression, top 100 by leader_rank)",
      "ccc_lr_verdict.tsv axis 2 evidence_summary (49/100 figure plus the joint reading)",
      1L, 1L, "Strong",
      "Hyp-2B;T-Compartment-suppression", "",
      "Microglia-involving rate at axis 2 is comparable to axis 1 (55) and axis 3 (51); the synaptic axis is densely populated with microglia-mediated cell-cell communication; supports Hyp-2B and T-Compartment-suppression but not the specific mechanism choice."),
  row("H2-042", "synaptic_suppression", "lr",
      "6,459 axis-2 LR cells in the universe-filtered set -- the LARGEST axis universe by LR cell count; the LR layer surfaces axis-2 biology that the TF and kinase layers cannot.",
      "+", "6,459 axis-2 LR cells (vs 6,074 amyloid; 4,139 interaction)",
      "ccc_lr_axis_restricted.tsv (axis=synaptic_suppression, total cells)",
      "ccc_lr_verdict.tsv axis 2 evidence_summary; the G1 triangulation hypothesis vindication",
      0L, 1L, "Strong",
      "T-Compartment-suppression", "",
      "Universe-size context claim documenting that the LR layer DOES carry information at axis 2 (the largest universe), validating the G1 triangulation rationale for adding LIANA+; the size itself is not model-discriminating."),
  row("H2-043", "synaptic_suppression", "pathway",
      "Six synaptic GO CC terms (GOCC_PRESYNAPSE, GOCC_SYNAPTIC_VESICLE_MEMBRANE, GOCC_SCHAFFER_COLLATERAL_CA1_SYNAPSE, GOCC_GABA_ERGIC_SYNAPSE, GOCC_EXCITATORY_SYNAPSE, GOCC_PRESYNAPTIC_MEMBRANE) are SUPPRESSED at NLGF arms with leader_score 11.6-15.4 and dominant_sign negative.",
      "-", "leader_score 11.6-15.4 (GOCC_PRESYNAPSE 15.4 leader); max_consistent_sign 2 (3 for the leader); dominant_sign -",
      "pathway_survey_unified_leaderboard.tsv (synaptic GO CC family rows; collection=gocc)",
      "ccc_lr_verdict.tsv axis 2 evidence_summary (the synaptic-suppression WHAT claim that Hyp-2A and Hyp-2B both accept)",
      2L, 1L, "Strong",
      "Hyp-2A;Hyp-2B;T-Compartment-suppression", "",
      "The synaptic-suppression observation itself is SHARED between Hyp-2A and Hyp-2B; only the MECHANISM differs (the LR layer adjudicates whether complement or TREM2/APP routes drive the suppression); STRONG by cross-modality consistent-sign breadth."),
  row("H2-044", "synaptic_suppression", "pathway",
      "Sun/Victor 2023 microglial substate MG3_ribosome_biogenesis (and related ribosomal compartments at breadth-1) is present at axis 2 in a manner consistent with the synaptic-axis ribosomal-pathway reading.",
      "-", "GOCC_MITOCHONDRIAL_LARGE_RIBOSOMAL_SUBUNIT leader_score 5.4 at tau_alone",
      "pathway_survey_unified_leaderboard.tsv (ribosomal compartment rows)",
      "fgsea_per_state.tsv (Sun/Victor MG3 module hits at synaptic-axis contrasts)",
      1L, 1L, "Suggestive",
      "T-Compartment-suppression", "",
      "Ribosomal-compartment signal at axis 2 is more weakly resolved than at axis 3; primarily a context claim that ties the synaptic-axis biology to the broader biosynthetic-suppression pattern."),
  row("H2-045", "synaptic_suppression", "lr",
      "Cd200_Cd200r1 (Neuronal->Microglia) and Vcan_Cd44 (OPC->Microglia_proliferative) appear in axis-2 top 25 ranks 18 + 17, supporting Hyp-2B's synaptic-adhesion-modulation route.",
      "+", "Cd200_Cd200r1 rank 18; Vcan_Cd44 rank 17 (mean_activity 0.94)",
      "ccc_lr_axis_restricted.tsv (axis=synaptic_suppression, ligands Cd200 + Vcan)",
      "ccc_lr_verdict.tsv axis 2 evidence_summary (Cd200 + Vcan + Ncam1 ligand family naming)",
      1L, 1L, "Moderate",
      "Hyp-2B;T-Compartment-suppression", "Hyp-2A",
      "Synaptic-adhesion modulation pairs cited explicitly in the Hyp-2B prediction; supports Hyp-2B over Hyp-2A.")
)

# ---- AXIS 3 (interaction_metabolic) ----------------------------------------

axis3 <- rbind(
  # --- TF layer (5 rows) ---
  row("H2-046", "interaction_metabolic", "tf",
      "Myc is coordinately suppressed at the tau x amyloid interaction; rank 2 in axis-3 axis-restricted top 5; the canonical biosynthetic master regulator points to a ribosomal / proliferation shutdown.",
      "-", "mean -3.25 at interaction (single contrast); 176 axis-3 substrate footprint",
      "tf_activity_axis_restricted.tsv (axis=interaction_metabolic, source=Myc, leader_rank=2)",
      "tf_activity_verdict.tsv axis 3 evidence_summary; kinase axis-3 Gsk3b H2-051 (the coordinated mechanism)",
      1L, 2L, "Strong",
      "Hyp-3B;T-Synergy", "Hyp-3A",
      "Coordinated negative direction across the axis-3 TF top 5; cross-layer corroboration via the kinase axis-3 ensemble; STRONG."),
  row("H2-047", "interaction_metabolic", "tf",
      "Creb1 is coordinately suppressed at the tau x amyloid interaction; rank 3 in axis-3 axis-restricted top 5.",
      "-", "mean -2.28 at interaction; 92 axis-3 substrate footprint",
      "tf_activity_axis_restricted.tsv (axis=interaction_metabolic, source=Creb1, leader_rank=3)",
      "tf_activity_verdict.tsv axis 3 evidence_summary; Myc H2-046 + Tp53 H2-048 + Jun H2-049 ensemble",
      1L, 1L, "Moderate",
      "Hyp-3B;T-Synergy", "Hyp-3A",
      "Coordinated biosynthetic-regulator suppression at the interaction; supports Hyp-3B."),
  row("H2-048", "interaction_metabolic", "tf",
      "Tp53 is coordinately suppressed at the tau x amyloid interaction; rank 4 in axis-3 axis-restricted top 5.",
      "-", "mean -2.15 at interaction; 143 axis-3 substrate footprint",
      "tf_activity_axis_restricted.tsv (axis=interaction_metabolic, source=Tp53, leader_rank=4)",
      "tf_activity_verdict.tsv axis 3 evidence_summary",
      1L, 1L, "Moderate",
      "Hyp-3B;T-Synergy", "Hyp-3A",
      "Tumour-suppressor TF suppressed at interaction; coordinated with Myc / Creb1 / Jun."),
  row("H2-049", "interaction_metabolic", "tf",
      "Jun is coordinately suppressed at the tau x amyloid interaction; rank 5 in axis-3 axis-restricted top 5.",
      "-", "mean -1.81 at interaction; 102 axis-3 substrate footprint",
      "tf_activity_axis_restricted.tsv (axis=interaction_metabolic, source=Jun, leader_rank=5)",
      "tf_activity_verdict.tsv axis 3 evidence_summary",
      1L, 1L, "Moderate",
      "Hyp-3B;T-Synergy", "Hyp-3A",
      "AP-1 family TF suppressed at interaction; reinforces the coordinated TF-layer ensemble."),
  row("H2-050", "interaction_metabolic", "tf",
      "A0A979HLR9 (the CollecTRI complex source that is POSITIVE at axis 1 NLGF arms +4.32 mean) shows a SIGN-REVERSAL at the interaction contrast (mean -4.85, the strongest negative axis-3 TF signal).",
      "-", "mean -4.85 at interaction vs +4.32 at axis 1; sign-reversal across contrasts",
      "tf_activity_axis_restricted.tsv (axis=interaction_metabolic, source=A0A979HLR9, leader_rank=1; axis=amyloid_activation leader_rank=1)",
      "tf_activity_verdict.tsv axis 3 evidence_summary; Rela sign-reversal H2-072 (analogous shape)",
      1L, 1L, "Suggestive",
      "Hyp-3B;T-Tau-attenuates;T-Synergy", "Hyp-3A",
      "Cross-axis sign-reversal claim mirroring the Rela sign-reversal; complex-source identifier limits direct biological interpretation but the reversal pattern is unambiguous."),

  # --- kinase layer (5 rows) ---
  row("H2-051", "interaction_metabolic", "kinase",
      "Gsk3b (canonical tau kinase phosphorylating tau at S202 / T205 / S396 / S404 -- the AT8 / PHF1 epitopes that stage AD tau pathology) is the STRONGEST single-kinase signal of the project: activated specifically at the tau x amyloid interaction with padj_ulm 0.00152; per-contrast amyloid_alone -1.66, tau_alone -1.47, interaction +4.32 -- the SYNERGY IS the activation.",
      "+", "score_at_interaction +4.77 (ulm) / +4.32 (per_contrast); padj_ulm 0.00152; padj_consensus 0.00018",
      "kinase_activity_axis_restricted.tsv (axis=interaction_metabolic, source=Gsk3b, leader_rank=1); kinase_activity_per_contrast.tsv (Gsk3b interaction padj_ulm 0.00152)",
      "tf_activity_axis_restricted.tsv Myc / Creb1 / Tp53 / Jun suppression at interaction (H2-046..H2-049); ccc_lr axis-3 axon-guidance rewiring (H2-057..H2-061)",
      1L, 3L, "Strong",
      "Hyp-3B;T-Synergy", "Hyp-3A",
      "THE strongest single-mechanism result of the entire mechanism arc; the per-contrast pattern (negative at single-insult arms, strongly positive at interaction) is the textbook synergy signature -- the activation only appears when amyloid and tau are combined; cross-layer corroboration via TF Myc lead + LR axon-guidance rewiring; STRONG by all three H1 criteria."),
  row("H2-052", "interaction_metabolic", "kinase",
      "Csnk1a1 (CK1-alpha) is suppressed at the tau x amyloid interaction; rank 2 in axis-3 axis-restricted top 5; the lone opposite-direction signal in the top 5.",
      "-", "mean -2.81; score_at_interaction -2.81 (ulm) / -3.22 (per_contrast); padj_ulm 0.063",
      "kinase_activity_axis_restricted.tsv (axis=interaction_metabolic, source=Csnk1a1, leader_rank=2); kinase_activity_per_contrast.tsv (Csnk1a1 interaction padj_ulm 0.063)",
      "kinase_activity_verdict.tsv axis 3 evidence_summary",
      1L, 1L, "Moderate",
      "Hyp-3B;T-Synergy", "Hyp-3A",
      "CK1-alpha substrate program suppressed on the tau-on-amyloid background; reaches FDR<0.10 only on ulm at the interaction; consistent direction at tau_in_nlgf (padj 0.041)."),
  row("H2-053", "interaction_metabolic", "kinase",
      "Mapk14 (p38-alpha MAPK stress-kinase activator of tau phosphorylation and microglial activation) is activated at the tau x amyloid interaction; rank 3 in axis-3 axis-restricted top 5.",
      "+", "mean +2.63 at interaction; padj_ulm 0.121 (just above FDR<0.10)",
      "kinase_activity_axis_restricted.tsv (axis=interaction_metabolic, source=Mapk14, leader_rank=3); kinase_activity_per_contrast.tsv (Mapk14 interaction padj_ulm 0.121)",
      "kinase_activity_verdict.tsv axis 3 evidence_summary (canonical tau-kinase ensemble with Gsk3b lead)",
      1L, 1L, "Moderate",
      "Hyp-3B;T-Synergy", "Hyp-3A",
      "Canonical stress-activated tau-kinase ensemble member; supports Hyp-3B by activating specifically at the interaction; FDR threshold narrowly missed on ulm."),
  row("H2-054", "interaction_metabolic", "kinase",
      "Cdk5 (second canonical tau kinase) is activated at the tau x amyloid interaction; rank 4 in axis-3 axis-restricted top 5.",
      "+", "mean +2.10 at interaction; 153 axis-3 substrate sites",
      "kinase_activity_axis_restricted.tsv (axis=interaction_metabolic, source=Cdk5, leader_rank=4)",
      "kinase_activity_verdict.tsv axis 3 evidence_summary; Hyp-0 cross-axis Cdk5 row H2-071",
      1L, 2L, "Moderate",
      "Hyp-0;Hyp-3B;T-Synergy", "Hyp-3A",
      "Cdk5 at axis 3 also feeds the Hyp-0 cross-axis Cdk5 integrator (Cdk5 in top 5 at all three axes); supports both Hyp-3B and Hyp-0."),
  row("H2-055", "interaction_metabolic", "kinase",
      "Cdk1 (master mitotic kinase, microglial cell-cycle re-engagement at the tau x amyloid interface) is activated at the interaction; rank 5 in axis-3 axis-restricted top 5.",
      "+", "mean +1.99 at interaction; 182 axis-3 substrate sites",
      "kinase_activity_axis_restricted.tsv (axis=interaction_metabolic, source=Cdk1, leader_rank=5)",
      "kinase_activity_verdict.tsv axis 3 evidence_summary",
      1L, 1L, "Moderate",
      "Hyp-3B;T-Synergy", "Hyp-3A",
      "Cell-cycle re-entry kinase activated at the interaction; supports Hyp-3B as part of the canonical tau-kinase + cell-cycle activation ensemble."),

  # --- LR layer (8 rows) ---
  row("H2-056", "interaction_metabolic", "lr",
      "Top-5 axis-3 LR pairs are mixed-sign (3+, 2-) -- the cell-type-resolved LR-layer rendition of the section-13 axis-3 mixed-sign reading (OXPHOS / ETC + ribosomal compartment + MG-M3 + DAM_up at interaction, modality-divergent signs).",
      "mixed", "top-5 mean_activity range [-0.54, +0.56]; 3-positive (Gpc3_Unc5c, Gas6_Axl, Nrxn1_Nlgn1) / 2-negative (Entpd1_Adora1, L1cam_Ephb2)",
      "ccc_lr_axis_restricted.tsv (axis=interaction_metabolic, leader_rank 1-5)",
      "ccc_lr_verdict.tsv axis 3 evidence_summary",
      1L, 1L, "Moderate",
      "Hyp-3B;T-Synergy", "Hyp-3A",
      "n_cells_used=1 across all top-5 rows; the mixed sign pattern is itself the signature of interaction-specific rewiring (not uniform up or down) supporting Hyp-3B's novel-mechanism claim."),
  row("H2-057", "interaction_metabolic", "lr",
      "Sema6a_Plxna4 is the sole axis-3 LR pair appearing as a 5/5-contrast G3 leader (leader_score 25.4) AND shows a SIGN-FLIP between amyloid-only and tau-on-amyloid backgrounds (nlgf_in_maptki:2/+ but nlgf_in_p301s:2/-).",
      "mixed", "leader_score 25.4 (sole top of cross-tool leader board); contrasts_summary nlgf_in_maptki:2/+ | nlgf_in_p301s:2/- | interaction:2/- | tau_alone:2/+ | tau_in_nlgf:2/-",
      "ccc_lr_cross_tool_leaderboard.tsv (Sema6a_Plxna4 OPC->Neuronal leader_score 25.4)",
      "ccc_lr_verdict.tsv axis 3 evidence_summary; Hyp-0 cross-axis Sema6a_Plxna4 row H2-073",
      2L, 1L, "Strong",
      "Hyp-0;Hyp-3B;T-Synergy;T-Tau-attenuates", "Hyp-3A",
      "Sole 5/5-contrast leader; the sign-flip between amyloid-only and tau-on-amyloid backgrounds is the strongest single-LR cross-contrast divergence claim; supports Hyp-3B (novel interaction mechanism), Hyp-0 (cross-axis Sema6a_Plxna4 divergence), T-Synergy, and T-Tau-attenuates (tau modifies the amyloid-driven sign)."),
  row("H2-058", "interaction_metabolic", "lr",
      "Gpc3_Unc5c (Vascular->Astrocyte axon-guidance co-receptor) gains signal at the tau x amyloid interface; rank 1 in axis-3 axis-restricted top 5.",
      "+", "mean_activity +0.558 at interaction; leader_rank 1; n_cells_used 1",
      "ccc_lr_axis_restricted.tsv (axis=interaction_metabolic, ligand=Gpc3, receptor=Unc5c, leader_rank=1)",
      "ccc_lr_verdict.tsv axis 3 evidence_summary (axon-guidance rewiring)",
      1L, 1L, "Moderate",
      "Hyp-3B;T-Synergy", "Hyp-3A",
      "Axon-guidance co-receptor signal at the interaction; supports Hyp-3B's axon-guidance rewiring."),
  row("H2-059", "interaction_metabolic", "lr",
      "Entpd1_Adora1 (OPC->OPC ATPase to adenosine A1 receptor purinergic signalling) loses signal at the tau x amyloid interface; rank 2 in axis-3 axis-restricted top 5.",
      "-", "mean_activity -0.544 at interaction; leader_rank 2; n_cells_used 1",
      "ccc_lr_axis_restricted.tsv (axis=interaction_metabolic, ligand=Entpd1, receptor=Adora1, leader_rank=2)",
      "ccc_lr_verdict.tsv axis 3 evidence_summary",
      1L, 1L, "Moderate",
      "Hyp-3B;T-Synergy", "Hyp-3A",
      "Purinergic signalling attenuated at the interaction; supports the mixed-sign rewiring pattern."),
  row("H2-060", "interaction_metabolic", "lr",
      "L1cam_Ephb2 (Neuronal->Neuronal Eph axon-guidance) loses signal at the tau x amyloid interface; rank 3 in axis-3 axis-restricted top 5.",
      "-", "mean_activity -0.505 at interaction; leader_rank 3; n_cells_used 1",
      "ccc_lr_axis_restricted.tsv (axis=interaction_metabolic, ligand=L1cam, receptor=Ephb2, leader_rank=3)",
      "ccc_lr_verdict.tsv axis 3 evidence_summary",
      1L, 1L, "Moderate",
      "Hyp-3B;T-Synergy", "Hyp-3A",
      "Eph axon-guidance attenuated at the interaction; supports axon-guidance rewiring."),
  row("H2-061", "interaction_metabolic", "lr",
      "Efna5_Epha4 (Microglia_homeostatic->Neuronal ephrin axon-guidance) gains signal at the interaction at broader scope; rank 8 in axis-3 axis-restricted top 25.",
      "+", "mean_activity +0.494 at interaction; leader_rank 8; n_cells_used 1",
      "ccc_lr_axis_restricted.tsv (axis=interaction_metabolic, ligand=Efna5, receptor=Epha4, leader_rank=8)",
      "ccc_lr_verdict.tsv axis 3 evidence_summary (the canonical microglia signal at axis 3)",
      1L, 1L, "Moderate",
      "Hyp-3B;T-Synergy", "Hyp-3A",
      "Direct Microglia_homeostatic->Neuronal ephrin signal; supports Hyp-3B's axon-guidance rewiring with a direct microglia-involving pair."),
  row("H2-062", "interaction_metabolic", "lr",
      "L1cam_Cd9 (Neuronal->Microglia) loses signal at the interaction at broader scope; rank 13 in axis-3 axis-restricted top 25.",
      "-", "mean_activity -0.49 at interaction; leader_rank 13; n_cells_used 1",
      "ccc_lr_axis_restricted.tsv (axis=interaction_metabolic, ligand=L1cam, receptor=Cd9, leader_rank=13)",
      "ccc_lr_verdict.tsv axis 3 evidence_summary",
      1L, 1L, "Moderate",
      "Hyp-3B;T-Synergy", "Hyp-3A",
      "Neuronal->Microglia L1CAM signal attenuated at interaction; complements Efna5_Epha4's positive direction in the rewired-not-uniform-suppressed reading."),
  row("H2-063", "interaction_metabolic", "lr",
      "51 of the top 100 axis-3 LR pairs involve a microglia sender or receiver -- the broader microglia-involving rate at the interaction axis; comparable to axes 1 (55/100) and 2 (49/100).",
      "+", "51/100 microglia-involving in axis-3 top 100",
      "ccc_lr_axis_restricted.tsv (axis=interaction_metabolic, top 100 by leader_rank)",
      "ccc_lr_verdict.tsv axis 3 evidence_summary",
      1L, 1L, "Moderate",
      "Hyp-3B;T-Synergy", "",
      "Microglia-involving rate context claim; not directly model-discriminating but documents that the interaction axis carries substantial microglia-mediated LR signal."),

  # --- pathway / module layer (5 rows) ---
  row("H2-064", "interaction_metabolic", "pathway",
      "MG-M3 hdWGCNA module (substate_Neuron_contam enrichment; Snap25;Stmn2 hubs; previously interpreted as the synaptic-fragment-containing microglial substate) is suppressed at the interaction with cross-modality breadth-3 mixed sign and the highest pathway leader_score of the project.",
      "mixed", "MG-M3 leader_score 23.4 (highest in pathway leader board); max_consistent_sign 2; max_n_modalities_sig 4; contrasts_summary nlgf_in_maptki:2/- | interaction:3/mixed | tau_alone:3/mixed | tau_in_nlgf:4/mixed; module DE logFC -2.91 at interaction (padj 0.074)",
      "pathway_survey_unified_leaderboard.tsv (MG-M3 leader_score 23.4); hdwgcna_module_de.tsv (MG-M3 interaction logFC -2.91 padj 0.074)",
      "hdwgcna_module_enrichment.tsv (MG-M3 substate_Neuron_contam odds_ratio 4.37 padj 0.58)",
      4L, 2L, "Strong",
      "Hyp-3B;T-Synergy", "Hyp-3A",
      "Highest pathway leader_score project-wide; cross-modality breadth 4 at tau_in_nlgf; supports Hyp-3B by carrying interaction-contrast-specific direction at the module level."),
  row("H2-065", "interaction_metabolic", "pathway",
      "OXPHOS / ETC family (8 GOBP/GOMF/GOCC terms including GOBP_OXIDATIVE_PHOSPHORYLATION, GOBP_ELECTRON_TRANSPORT_CHAIN, GOCC_NADH_DEHYDROGENASE_COMPLEX, GOCC_CYTOCHROME_COMPLEX) shows mixed-sign breadth-3 at the interaction.",
      "mixed", "leader_score 11.4-16.4 across the 8 terms; max_n_modalities_sig 3; contrasts_summary interaction:3/mixed for all 8",
      "pathway_survey_unified_leaderboard.tsv (GO OXPHOS / ETC family rows)",
      "MG-M3 H2-064 covariant (mitochondrial / OXPHOS overlap)",
      3L, 1L, "Strong",
      "Hyp-3B;T-Synergy", "Hyp-3A",
      "Coordinated OXPHOS / electron-transport mixed-sign signal at the interaction (proteomics + phospho + RNA divergent directions) consistent with the Myc-led biosynthetic shutdown reading; STRONG by cross-modality breadth."),
  row("H2-066", "interaction_metabolic", "pathway",
      "Custom microglial-state gene set DAM_up shows mixed-sign breadth-3 at the interaction in addition to the +/+ NLGF arms.",
      "mixed", "DAM_up leader_score 16.4; contrasts_summary nlgf_in_maptki:2/+ | nlgf_in_p301s:2/+ | interaction:3/mixed",
      "pathway_survey_unified_leaderboard.tsv (collection=custom_microglia_states, pathway=DAM_up)",
      "Axis-1 DAM_up H2-026 covariant (the +/+ NLGF arms; this row records the interaction contrast specifically)",
      3L, 1L, "Moderate",
      "Hyp-3B;T-Synergy", "Hyp-3A",
      "The DAM-state signature acquires an interaction-specific modality-divergent component beyond the axis-1 +/+ pattern; supports Hyp-3B by interaction-specific direction."),
  row("H2-067", "interaction_metabolic", "pathway",
      "Ribosomal compartment (GOCC_MITOCHONDRIAL_LARGE_RIBOSOMAL_SUBUNIT and related GO terms) is suppressed at axis-3-relevant contrasts, consistent with the Myc-led biosynthetic shutdown.",
      "-", "ribosomal GO term leader_score 5.4 each; dominant_sign - at tau_alone",
      "pathway_survey_unified_leaderboard.tsv (ribosomal compartment rows)",
      "Myc H2-046 + MG-M3 H2-064 + OXPHOS H2-065 covariant (the biosynthetic-suppression integrated reading)",
      1L, 1L, "Suggestive",
      "Hyp-3B;T-Synergy", "Hyp-3A",
      "Ribosomal-pathway hit count at axis 3 is modest in absolute number but coherent in direction with the Myc-led biosynthetic-shutdown reading."),
  row("H2-068", "interaction_metabolic", "pathway",
      "Triangulation interaction analysis identifies Colec12 (DAM-state marker; max DAM_logFC +2.29 across 4-DE-method consensus; mean_negLogP 4.43) as the top interaction-specific gene by cross-method consensus.",
      "+", "Colec12 mean_negLogP 4.43 across 6 methods; limma_voom_pb 3.18; deseq2_pb 6.41; glmmtmb_sc 6.29",
      "triangulation_top_interaction.tsv (Colec12 row 1)",
      "nebula_per_state_localisation.tsv (Colec12 DAM_neg_log10p 3.65)",
      1L, 1L, "Moderate",
      "Hyp-3B;T-Synergy;T-Inflammation", "Hyp-3A",
      "Interaction-specific gene; DAM-localised; supports Hyp-3B by per-contrast localisation; also supports T-Inflammation as the DAM-state signature is inflammatory-coded."),
  row("H2-069", "interaction_metabolic", "pathway",
      "MG-M2 hdWGCNA module DE at tau_in_nlgf is up-regulated with logFC +2.58 padj 0.027 -- the same DAM-marker module that drives axis-1 H2-023 also carries interaction-relevant signal at tau-on-amyloid.",
      "+", "MG-M2 tau_in_nlgf logFC +2.58 padj 0.027",
      "hdwgcna_module_de.tsv (module=MG-M2, contrast=tau_in_nlgf)",
      "MG-M2 axis-1 H2-023 covariant; LR axis-3 Efna5_Epha4 + Sema6a_Plxna4 cross-axis rows",
      4L, 2L, "Moderate",
      "Hyp-3B;T-Synergy;T-Inflammation", "Hyp-3A",
      "The axis-1 DAM-marker module is responsive at tau_in_nlgf; suggests the T-Inflammation theme carries some interaction-axis support and weakly weakens T-Synergy's purity claim."),
  row("H2-070", "interaction_metabolic", "pathway",
      "MG-M3 hdWGCNA module DE at tau_in_nlgf shows logFC -3.80 padj 6.7e-4 -- the strongest single-module interaction-relevant signal of the project.",
      "-", "MG-M3 tau_in_nlgf logFC -3.80 padj 6.7e-4",
      "hdwgcna_module_de.tsv (module=MG-M3, contrast=tau_in_nlgf)",
      "MG-M3 axis-3 H2-064 covariant",
      1L, 1L, "Strong",
      "Hyp-3B;T-Synergy", "Hyp-3A",
      "MG-M3 module DE at tau_in_nlgf reaches FDR<0.001; STRONG; the strongest module-DE result at an interaction-relevant contrast.")
)

# ---- CROSS-AXIS (Hyp-0) + LAYER-LEVEL OBSERVATIONS ----------------------------

cross_axis <- rbind(
  row("H2-071", "cross_axis", "kinase",
      "Cdk5 is the only kinase in the top 5 at all three axes (ranks 2 / 2 / 4) and its KSN substrate footprint scales with the axis universe size (81 sites at amyloid, 263 at synaptic, 153 at interaction).",
      "+", "Cdk5 leader_rank: amyloid 2, synaptic 2, interaction 4; substrate footprint: 81/263/153 sites",
      "kinase_activity_axis_restricted.tsv (source=Cdk5 across three axes)",
      "Per-axis Cdk5 rows H2-009 (axis 1), H2-033 (axis 2 KSN footprint), H2-054 (axis 3); per_contrast.tsv (Cdk5 nlgf_in_p301s padj_ulm 0.0077)",
      1L, 3L, "Strong",
      "Hyp-0", "",
      "Cross-axis Cdk5 integrator evidence; STRONG by cross-axis presence at all three axes; supports Hyp-0 directly."),
  row("H2-072", "cross_axis", "tf",
      "Rela sign-reversal across contrasts: positive at axis-1 NLGF arms (+2.57 / +1.81 at maptki / p301s) but negative at axis-3 interaction (-1.36). Coordinated with Nfkb1 (axis 1 +2.47 -> axis 3 -1.04) and Spi1 (axis 1 +2.81 -> axis 3 -0.93).",
      "mixed", "Rela: +2.57/+1.81 at NLGF arms, -1.36 at interaction; Nfkb1: +2.47 -> -1.04; Spi1: +2.81 -> -0.93",
      "tf_activity_axis_restricted.tsv (Rela / Nfkb1 / Spi1 rows across axes)",
      "Per-axis Rela H2-006 (the +ve half) + H2-046 / H2-049 axis-3 TF suppression ensemble",
      1L, 2L, "Strong",
      "Hyp-0;Hyp-1B;T-Tau-attenuates;T-Synergy", "Hyp-1A",
      "Load-bearing Hyp-1B + T-Tau-attenuates evidence; the sign-reversal IS the Hyp-1B signature; STRONG by cross-axis pattern + per-contrast coordination across multiple NF-kB-family TFs."),
  row("H2-073", "cross_axis", "lr",
      "Sema6a_Plxna4 cross-axis divergence: appears at axis 1 with positive cross-tool leader signal at maptki, but the same pair sign-flips between amyloid-only and tau-on-amyloid at the interaction contrast (leader_score 25.4 sole 5/5-contrast LR leader of the project).",
      "mixed", "leader_score 25.4 (sole top of cross-tool LR leader board); 5/5-contrast leader status",
      "ccc_lr_cross_tool_leaderboard.tsv (Sema6a_Plxna4 OPC->Neuronal)",
      "Per-axis Sema6a_Plxna4 H2-057 (axis 3 sign-flip detail)",
      2L, 1L, "Strong",
      "Hyp-0;Hyp-3B;T-Synergy;T-Tau-attenuates", "Hyp-3A",
      "Cross-axis single-LR divergence; STRONG by sole 5/5-contrast leader status; the project's most-divergent LR pair across contrasts."),
  row("H2-074", "cross_axis", "cross_layer",
      "TF and kinase regulator layers AGREE on the axis-2 non-finding pattern (both top 5s = axis-1 top 5; both reach 0 FDR-significant axis-2-distinct hits) -- the regulator layers cross-validate that the synaptic axis is not driven by transcriptional or phospho-signalling regulators at this resolution.",
      "mixed", "TF axis-2 top 5 = axis-1 top 5; kinase axis-2 top 5 = axis-1 top 5; 0 distinct hits in either layer",
      "tf_activity_verdict.tsv axis 2 evidence_summary + kinase_activity_verdict.tsv axis 2 evidence_summary",
      "LR axis-2 evidence H2-035..H2-042 (the post-regulator mechanism that the LR layer surfaces)",
      0L, 2L, "Strong",
      "T-Compartment-suppression", "",
      "Cross-layer agreement on the regulator-layer non-finding; STRONG by 2-layer cross-validation; supports T-Compartment-suppression's reading of the synaptic axis as post-regulator intercellular mechanism."),
  row("H2-075", "cross_axis", "cross_layer",
      "Axis-3 TF, kinase, and LR layers AGREE on a coordinated interaction-specific mechanism: TF Myc-led biosynthetic suppression (4-of-5 axis-3 top-5 negative) + kinase Gsk3b-led tau-kinase activation (4-of-5 axis-3 top-5 positive, Gsk3b padj 0.00152) + LR axon-guidance / synaptic-adhesion rewiring (mixed-sign top 5).",
      "mixed", "TF axis-3 4-/5- top-5 negative + 1+ Hdac1 positive at +1.72; kinase axis-3 4+/5+ top-5 positive; LR axis-3 3+/2- mixed",
      "tf_activity_verdict.tsv axis 3 evidence_summary + kinase_activity_verdict.tsv axis 3 evidence_summary + ccc_lr_verdict.tsv axis 3 evidence_summary",
      "Gsk3b H2-051 + Myc H2-046 + Sema6a_Plxna4 H2-057 (the three-layer mechanism anchors)",
      0L, 3L, "Strong",
      "Hyp-3B;T-Synergy", "Hyp-3A",
      "Cross-layer coordination on the interaction-specific synergy mechanism; STRONG by 3-layer cross-validation; the most-integrated three-layer mechanism reading of the project.")
)

# ---- PHASE I (per-state NF-kB attenuation test) ----------------------------
#
# Six rows added at I4 (2026-05-25) propagating the per-state NF-kB attenuation
# test outcome from rmd/17_nfkb_attenuation.Rmd (section 18 of analysis.html)
# into the H2 claims ledger. Phase I tested whether the T-Tau-attenuates
# Rela sign-reversal that anchors H2-072 is a GLOBAL whole-microglia effect, a
# STATE-RESTRICTED effect, or a STATE-DRIVEN-BY-ABUNDANCE compositional
# artefact. The verdict tier resolved to GLOBAL ATTENUATION: 4 of 4 microglia
# substates (MG_homeostatic, MG_DAM, MG_IFN, MG_proliferative) independently
# carry NF-kB suppression at the interaction contrast under the chapter's
# cross-layer gate (TF complex consensus < 0 AND NF-kB-target-GSEA NES < 0
# AND >=1 layer at FDR<0.10 within state), with the section 11 substate-
# composition ANOVA (tau:nlgf F=1.03 p=0.330) formally ruling out the
# compositional artefact. Numbering I-001..I-006 (Phase I; visually distinct
# from H2-001..H2-075 to flag Phase I origin in the ledger view).

phase_i <- rbind(
  row("I-001", "cross_axis", "cross_layer",
      "Per-state NF-kB attenuation at the interaction in MG_homeostatic: TF complex A0A979HLR9 consensus score is negative AND NF-kB-target GSEA NES is negative, with cross-layer sign agreement.",
      "-", "TF consensus score -5.45 (p=5e-8); ulm score -2.35 (padj_ulm 0.298, NS in heavy-correction regime); GSEA NFKB_UNION_TARGETS NES -1.221 (padj 0.091, sig at FDR<0.10); 718 union targets after homeostatic-universe intersection",
      "per_state_tf_activity.rds (state=homeostatic, contrast=interaction, source=A0A979HLR9, consensus statistic); per_state_nfkb_target_gsea.tsv (state=homeostatic, contrast=interaction, pathway=NFKB_UNION_TARGETS)",
      "per_state_tf_activity.tsv (ulm layer at homeostatic interaction); section 14 whole-microglia Rela sign-reversal H2-072; per-state GSEA mirror at NLGF arms NES +1.21 / +1.18 in homeostatic (the sign-reversal companion at the NLGF arms)",
      1L, 2L, "Strong",
      "Hyp-1B;T-Tau-attenuates", "Hyp-1A",
      "Substate-specific verdict gate (tf_neg AND gsea_neg AND >=1 layer FDR<0.10) passes; homeostatic also passes FDR<0.10 on GSEA (padj 0.091) so both layers reach significance at the chapter-locked sig metrics; cross-layer sign agreement + GSEA NES sign-reversal between NLGF arms and interaction = STRONG per H1 rules."),
  row("I-002", "cross_axis", "cross_layer",
      "Per-state NF-kB attenuation at the interaction in MG_DAM is the strongest per-substate signal of the Phase I test: TF complex consensus reaches p=7e-11 and GSEA NES reaches padj 0.027.",
      "-", "TF consensus score -6.52 (p=7e-11, the deepest signal across the 4 substates); ulm score -2.34 (padj_ulm 0.606, NS in heavy-correction regime); GSEA NFKB_UNION_TARGETS NES -1.247 (padj 0.027, sig at FDR<0.10); 750 union targets after DAM-universe intersection",
      "per_state_tf_activity.rds (state=DAM, contrast=interaction, source=A0A979HLR9, consensus statistic); per_state_nfkb_target_gsea.tsv (state=DAM, contrast=interaction, pathway=NFKB_UNION_TARGETS)",
      "per_state_tf_activity.tsv (ulm layer at DAM interaction); section 14 whole-microglia Rela sign-reversal H2-072; per-state GSEA mirror at NLGF arms NES +1.10 / +0.99 in DAM",
      1L, 2L, "Strong",
      "Hyp-1B;T-Tau-attenuates", "Hyp-1A",
      "DAM is the LEAD CARRIER of the per-state attenuation (largest magnitude in both TF and GSEA layers across the 4 substates) which is biologically intuitive because DAM is transcriptionally NF-kB-HIGH under amyloid so attenuation has the most signal where the baseline is most amplified; STRONG by all three H1 criteria; the H4 17.4 section-4 hypothesis sketch predicted MG_homeostatic carries the suppression and the data REFINE (rather than refute) the prediction -- every substate carries it but DAM has the largest magnitude."),
  row("I-003", "cross_axis", "cross_layer",
      "Per-state NF-kB attenuation at the interaction in MG_IFN is present but weaker than DAM / homeostatic: TF complex consensus reaches p=3e-4 but GSEA NES is sign-consistent without reaching FDR<0.10.",
      "-", "TF consensus score -3.66 (p=3e-4); ulm score -1.33 (padj_ulm 0.698, NS); GSEA NFKB_UNION_TARGETS NES -1.126 (padj 0.425, NS; raw pval 0.112 borderline); 781 union targets after IFN-universe intersection",
      "per_state_tf_activity.rds (state=IFN, contrast=interaction, source=A0A979HLR9, consensus statistic); per_state_nfkb_target_gsea.tsv (state=IFN, contrast=interaction, pathway=NFKB_UNION_TARGETS)",
      "per_state_tf_activity.tsv (ulm layer at IFN interaction); section 14 whole-microglia Rela sign-reversal H2-072",
      1L, 2L, "Moderate",
      "Hyp-1B;T-Tau-attenuates", "Hyp-1A",
      "IFN substate passes the chapter verdict gate (tf_neg + gsea_neg + TF FDR<0.10 on consensus) but GSEA does not corroborate at FDR<0.10 (padj 0.425); the GSEA sign agrees with the TF activity but the raw pval 0.112 is just outside threshold; MODERATE per H1 rules (cross-layer sign agreement without cross-layer significance)."),
  row("I-004", "cross_axis", "cross_layer",
      "Per-state NF-kB attenuation at the interaction in MG_proliferative is the weakest per-substate signal: TF complex consensus reaches p=0.011 and GSEA NES is sign-consistent without reaching FDR<0.10.",
      "-", "TF consensus score -2.55 (p=0.011); ulm score -0.24 (padj_ulm 0.943, NS; ulm and consensus statistics diverge sharply in the smallest substate); GSEA NFKB_UNION_TARGETS NES -1.121 (padj 0.162, NS; raw pval 0.081 borderline); 763 union targets after proliferative-universe intersection",
      "per_state_tf_activity.rds (state=proliferative, contrast=interaction, source=A0A979HLR9, consensus statistic); per_state_nfkb_target_gsea.tsv (state=proliferative, contrast=interaction, pathway=NFKB_UNION_TARGETS)",
      "per_state_tf_activity.tsv (ulm layer at proliferative interaction); section 14 whole-microglia Rela sign-reversal H2-072",
      1L, 2L, "Moderate",
      "Hyp-1B;T-Tau-attenuates", "Hyp-1A",
      "Proliferative is the smallest substate pool and the noisiest; consensus statistic recovers sig (p=0.011) where ulm fails; GSEA sign-consistent but raw pval 0.081 is borderline; MODERATE per H1 rules; the substate-specific weakness is biologically reasonable because proliferative microglia are transcriptionally re-engaged on a cell-cycle program rather than carrying the NF-kB-HIGH DAM program that gets attenuated."),
  row("I-005", "cross_axis", "cross_layer",
      "Per-state replication of the whole-microglia NF-kB sign-reversal: ALL 4 microglia substates independently carry the attenuation at the interaction; the section 14 whole-microglia Rela / Nfkb1 / Spi1 sign-reversal is the weighted aggregate of consistent per-state suppression, NOT a one-substate-only signal or an aggregation artefact.",
      "-", "TF complex consensus < 0 in 4/4 substates (DAM -6.52 p=7e-11; homeostatic -5.45 p=5e-8; IFN -3.66 p=3e-4; proliferative -2.55 p=0.011); GSEA NFKB_UNION_TARGETS NES < 0 in 4/4 substates (DAM -1.247 padj 0.027; homeostatic -1.221 padj 0.091; IFN -1.126 padj 0.425; proliferative -1.121 padj 0.162); 2 of 4 substates pass FDR<0.10 on GSEA; whole-microglia replication of section-14 Rela -1.36 / Nfkb1 -1.04 / Spi1 -0.93",
      "per_state_tf_activity.rds (4 substates x interaction contrast x A0A979HLR9 consensus statistic); per_state_nfkb_target_gsea.tsv (4 substates x interaction contrast x NFKB_UNION_TARGETS)",
      "tf_activity_axis_restricted.tsv (whole-microglia Rela / Nfkb1 / Spi1 sign-reversal at axis 3 interaction); H2-072 (cross-axis Rela sign-reversal anchor)",
      1L, 2L, "Strong",
      "Hyp-1B;T-Tau-attenuates", "Hyp-1A",
      "Cross-substate 4/4 sign-consistent replication is the strongest available corroboration that the whole-microglia signal is genuinely transcriptional (not an aggregation artefact or one-substate signal); STRONG by 4/4 sign consistency + 2 layers + FDR<0.10 in at least one layer per substate; this is the row that materially LIFTS T-Tau-attenuates from its Phase H position by demonstrating cross-substate breadth that the whole-microglia signal alone could not establish."),
  row("I-006", "cross_axis", "cross_layer",
      "Compositional-vs-transcriptional adjudication: the section 11 DAM-sampling-balance ANOVA rules out 'state-driven-by-abundance' as the source of the NF-kB attenuation; the per-state signal is transcriptional, not compositional.",
      "-", "ccc_dam_sampling_balance_anova.tsv: tau:nlgf F=1.03 p=0.330 (interaction term NS, rules out compositional confound); nlgf main F=63.03 p=4.1e-6 (NLGF expands DAM substate as expected); tau main F=0.07 p=0.79 (tau alone does not shift substate composition); residual df=12",
      "ccc_dam_sampling_balance_anova.tsv (full 2x2 factorial ANOVA on DAM-fraction outcome)",
      "I-001..I-005 (per-state per-substate NF-kB attenuation rows that the ANOVA validates as transcriptional rather than compositional)",
      1L, 1L, "Strong",
      "Hyp-1B;T-Tau-attenuates", "Hyp-1A",
      "The tau:nlgf interaction term in the substate-composition ANOVA is non-significant (p=0.330) which formally rules out 'state-driven-by-abundance' as the source of the per-state NF-kB attenuation; combined with I-001..I-005's per-state transcriptional evidence, the verdict tier resolves to GLOBAL ATTENUATION; STRONG by formal statistical adjudication of the compositional-vs-transcriptional alternative.")
)

# ---- PHASE J (causal signalling-network topology) --------------------------
#
# Two rows added at J5 (2026-06-04) propagating the CARNIVAL causal-network
# reconstruction (rmd/19_causal_network.Rmd, section 20 of analysis.html) into
# the H2 claims ledger. The causal layer solved a per-contrast OmniPath ILP
# linking the kinase activity endpoints (signed inputObj <- phospho_corrected)
# to the TF activity endpoints (signed measObj <- tf_activity_decoupler_split)
# over the microglia-pruned mouse PKN. These rows score the RECOVERED DIRECTED
# WIRING between the endpoints -- connectivity that the activity layers
# (sections 14 / 15) scored independently but could not wire -- NOT a re-count
# of those endpoint activities. Both are graded SUGGESTIVE per the J5
# anti-double-count guard: the seeds are bulk-hippocampus phospho (not
# microglia-sorted; a cross-compartment bridge to the microglial TF readout),
# the path length is bounded by the L=3 reachability horizon, intermediates are
# inferred, and CARNIVAL weights are the consensus across optimal solutions.
# Numbering J-001..J-002 (Phase J; visually distinct from H2-* and I-* to flag
# Phase J origin in the ledger view).

phase_j <- rbind(
  row("J-001", "interaction_metabolic", "cross_layer",
      "Causal-topology bridge at the interaction: the OmniPath ILP (CARNIVAL) wires the section-15 interaction-lead kinase Gsk3b to the section-14.3 interaction-lead TF Myc through a single inferred Mapk14 (p38-alpha) -- Gsk3b -| Mapk14 -> Myc and -> Foxo3 -- so the recovered topology explains the measured Myc / Foxo3 suppression as downstream of Gsk3b activation.",
      "mixed", "interaction net: Gsk3b input activity +100 (up) -| Mapk14 inferred -100 (down, edge sign -1) -> Myc / Foxo3 measured -100 (down, edge sign +1); 5/5 measured TFs recovered; 8 active edges; 1 inferred intermediate (Mapk14); single-insult arms nlgf_in_maptki / tau_alone returned honestly-empty networks (0 kinase seeds)",
      "causal_network_edges.tsv (contrast=interaction: Gsk3b -1 Mapk14; Mapk14 +1 Myc; Mapk14 +1 Foxo3); causal_network_verdict.tsv (axis=interaction_metabolic, key_recovered_path)",
      "kinase Gsk3b H2-051 (interaction-lead kinase, padj_ulm 0.00152); TF Myc H2-046 (interaction-lead TF suppression); cross-layer H2-075 (3-layer interaction-mechanism agreement); causal_network_nodes.tsv (interaction node activities)",
      1L, 2L, "Suggestive",
      "Hyp-3B;T-Synergy", "Hyp-3A",
      "Suggestive and NOT a re-count of the Gsk3b / Myc activity endpoints (H2-051 / H2-046) -- it scores the newly-recovered directed signed wiring between them. The bridge is interaction-specific (the single-insult arms returned 0-seed empty networks), so a recovered interaction-only path is incremental evidence for Hyp-3B / T-Synergy and against Hyp-3A's 'no interaction-specific mechanism'; mirrors the H2-075 supports/contradicts pattern at the topology layer. Held at Suggestive by the cross-compartment bulk-phospho seed, the L=3 reachability horizon, the single inferred Mapk14 intermediate, and the multi-solution-consensus weighting."),
  row("J-002", "cross_axis", "cross_layer",
      "Causal-topology probe of the section-18 NF-kB attenuation: the ILP recovers a Gsk3b-rooted route to the NF-kB regulatory module, but it is contrast-dependent and equivocal -- an inhibitory Gsk3b -| Ikbkg(NEMO) -> Nfkbib(IkBb) arm at tau_in_nlgf (NEMO-inhibition direction-consistent with reduced NF-kB activation), an activating Gsk3b -> Nfkb1 arm at nlgf_in_p301s (NF-kB up), and NO NF-kB node recovered at the interaction contrast where the attenuation is defined.",
      "mixed", "tau_in_nlgf: Gsk3b +100 (up) -| Ikbkg inferred -100 (down, edge sign -1) -> Nfkbib measured -100 (down, edge sign +1); nlgf_in_p301s: Gsk3b +100 -> Nfkb1 measured +100 (edge sign +1, NF-kB up); interaction: 0 NF-kB nodes among 8 active edges",
      "causal_network_edges.tsv (tau_in_nlgf: Gsk3b -1 Ikbkg; Ikbkg +1 Nfkbib | nlgf_in_p301s: Gsk3b +1 Nfkb1 | interaction: no NF-kB node); causal_network_verdict.tsv (nfkb_path column across axes)",
      "Phase I per-state NF-kB attenuation I-001..I-006 (the interaction-contrast attenuation EFFECT this route probes upstream of); section-14 Rela / Nfkb1 / Spi1 sign-reversal H2-072",
      1L, 2L, "Suggestive",
      "T-Tau-attenuates", "",
      "Suggestive and SUPPORT-ONLY (no contradict): the topological support for the attenuation path is partial -- direction-consistent OFF the interaction (the tau_in_nlgf NEMO-inhibition arm) but ABSENT ON the interaction contrast itself, where I-001..I-006 localise the attenuation; per the anti-anchoring lock the interaction-contrast absence is reported as a finding, not smoothed over. It does NOT contradict Hyp-1A because the evidence is equivocal, and it is NOT a re-count of I-001..I-006 (which score the attenuation effect) -- this row scores the recovered upstream route. The IkBb (Nfkbib, itself an NF-kB inhibitor) sign-chain is more ambiguous than the NEMO arm; the route is credited for the Gsk3b -| NEMO connection only. Cross-compartment bulk-phospho seed + L=3 horizon caveats apply.")
)

# ---- PHASE K (data-driven SCENIC regulon corroboration) --------------------
#
# Two rows added at K6 (2026-06-05) propagating the SCENIC data-driven
# gene-regulatory inference (rmd/20_scenic_regulons.Rmd, section 21 of
# analysis.html) into the H2 claims ledger. SCENIC infers regulons de novo from
# this dataset's own microglial co-expression (GRNBoost2) + cisTarget motif
# enrichment -- the project's first TF-layer evidence that does NOT lean on a
# literature prior (section 14 = CollecTRI prior). These rows score the
# data-driven RECOVERY of two comparison TFs as motif-supported microglial
# regulons at the strictest 10/10-run consensus -- regulon EXISTENCE that the
# prior-based activity layer (section 14) could assert only from a curated
# network -- and explicitly NOT activity concordance: the controlled
# CollecTRI -> SCENIC network swap on the same NEBULA z is mostly discordant
# (median target Jaccard 0.004; the recovered Spi1 activity sign-flips and loses
# significance). Both are graded SUGGESTIVE per the K6 anti-double-count guard:
# the rows score the methodologically independent motif + co-expression
# recovery, NOT the prior-based footprints already counted (Spi1/Nfkb1 amyloid
# H2-001/H2-002; per-state NF-kB attenuation I-001..I-006), and are held at
# Suggestive by SCENIC's caveats (microglia-only scope, single-nucleus input,
# stochastic GRN, motif-prior dependence). Both are MARGIN-NEUTRAL by
# construction: K-001 supports BOTH amyloid contest models as a shared
# activation premise (it does not discriminate the tau-modulation contest) and
# K-002 tags only the T-Inflammation theme, so the locked contest margins
# (amyloid 18 / synaptic 12 / interaction 55) are unchanged. Numbering
# K-001..K-002 (Phase K; visually distinct from H2-*, I-*, J-* to flag Phase K
# origin in the ledger view).

phase_k <- rbind(
  row("K-001", "amyloid_activation", "tf",
      "Data-driven regulon recovery of the amyloid TF Spi1: SCENIC's de novo GRN (GRNBoost2 co-expression + cisTarget motif pruning, no literature prior) recovers Spi1 as a motif-supported microglial regulon at the strictest 10/10-run consensus, corroborating -- by a method independent of the section-14 CollecTRI prior -- that Spi1 is a genuine amyloid-activation regulator in this microglial dataset; scoped to regulon EXISTENCE, not activity concordance.",
      "+", "Spi1 in_consensus=TRUE, present +:10/10 runs, 15 high-confidence targets (>=8/10 edge recurrence, motif NES>=3.0); controlled CollecTRI -> SCENIC network swap on the same NEBULA z: activity sign-flips and loses significance (CollecTRI consensus +2.73 sig padj 0.051 -> SCENIC consensus -0.43 ns padj 0.207 at nlgf_in_maptki); data-driven target set near-disjoint from the prior (Jaccard 0.005; 2 of 381 union targets shared)",
      "scenic_regulon_recovery.tsv (Spi1: in_consensus=TRUE, max_runs_present=10, n_consensus_targets=15, recovery_class=recovered_comparison); scenic_recovery_census.tsv (axis amyloid_activation, Spi1 present +:10, consensus +:15)",
      "section-14.1 prior-based Spi1 amyloid activity H2-001 (the prior footprint this data-driven recovery corroborates, NOT re-counts); scenic_vs_collectri_headtohead.tsv (Spi1 activity discordance across contrasts); scenic_target_overlap.tsv (Spi1 Jaccard 0.005); section 21 SCENIC chapter 21.1-21.2",
      1L, 1L, "Suggestive",
      "Hyp-1A;Hyp-1B;T-Inflammation", "",
      "Suggestive and NOT a re-count of the prior-based Spi1 amyloid activity (H2-001) -- it scores the methodologically independent recovery of Spi1 as a de novo data-driven microglial regulon (no CollecTRI prior) at the strictest >=8/10 consensus bar. Scoped to regulon EXISTENCE, not activity concordance: the controlled network swap shows Spi1 activity sign-flips and loses significance, and the data-driven target set is near-disjoint from the prior (Jaccard 0.005, 2 of 381 union targets), so this corroborates Spi1's data-driven regulon footprint, NOT the activity direction. Supports BOTH Hyp-1A and Hyp-1B (a shared amyloid-activation premise that does not discriminate the tau-modulation contest), so the amyloid margin is unchanged at 18. Held at Suggestive by SCENIC caveats: microglia-only scope (state regulons, not cell-type identity), single-nucleus input, stochastic GRN, motif-prior dependence."),
  row("K-002", "cross_axis", "tf",
      "Data-driven regulon recovery of an NF-kB-family TF: SCENIC's de novo GRN recovers Rel as the SOLE motif-supported NF-kB-family regulon at the 10/10-run consensus, corroborating that microglial NF-kB signalling has a genuine data-driven regulon footprint (T-Inflammation); scoped to regulon EXISTENCE, not activity concordance -- the recovered Rel regulon does NOT reproduce the section-18 interaction attenuation and the other NF-kB members are not recovered at the strict bar, so the data-driven layer gives the NF-kB axis only thin, single-member existence support.",
      "+", "Rel in_consensus=TRUE, present +:10/10 runs, 9 high-confidence targets; sole NF-kB-family regulon at the >=8/10 bar (Nfkb1 7/10 sub-threshold; Relb present +:10/10 but <5 stable targets; Nfkb2 0/10 motif-pruned; Rela 1/10 weak); per-substate interaction non-significant in all 4 substates (homeostatic / DAM / IFN / proliferative; IFN padj 0.109 closest) with mixed-sign consensus scores -- does NOT reproduce the section-18 all-substate negative attenuation",
      "scenic_regulon_recovery.tsv (Rel: in_consensus=TRUE, max_runs_present=10, n_consensus_targets=9, recovery_class=recovered_comparison); scenic_recovery_census.tsv (nfkb_family_sec18: Rel consensus +:9; Nfkb1 7/10; Relb 10/10 present; Nfkb2 0/10; Rela 1/10)",
      "section-18 per-state NF-kB attenuation I-001..I-006 (the prior-based attenuation EFFECT this row does NOT re-count and does NOT reproduce); section-14 Nfkb1 amyloid H2-002; scenic_per_substate_interaction.tsv (Rel ns at interaction in 4/4 substates); section 21 SCENIC chapter 21.1, 21.4",
      1L, 1L, "Suggestive",
      "T-Inflammation", "",
      "Suggestive and SUPPORT-ONLY (no contradict): NOT a re-count of the section-18 NF-kB attenuation (I-001..I-006) or the prior-based Nfkb1 amyloid (H2-002) -- it scores that the de novo data-driven GRN independently recovers an NF-kB-family regulon (Rel) in microglia, an existence corroboration of microglial NF-kB activity (-> T-Inflammation). Scoped strictly to regulon existence, NOT activity concordance: the recovered Rel regulon does NOT reproduce the section-18 interaction attenuation (flat / non-significant at interaction in all 4 substates, mixed-sign consensus), so this row does NOT support T-Tau-attenuates; and the other NF-kB members are unrecovered at the strict bar (Nfkb1 7/10, Relb <5 stable targets, Nfkb2 motif-pruned, Rela 1/10) -- single-member, partial support. Tags only the T-Inflammation theme, so all three contest margins are unchanged. Held at Suggestive by the same SCENIC caveats (microglia-only, single-nucleus, stochastic GRN, motif-prior dependence).")
)

# ---- PHASE L (GeoMx spatial-deconvolution tissue-composition) --------------
#
# Two rows added at L5 (2026-06-05) propagating the GeoMx SpatialDecon
# two-stage tissue-deconvolution (rmd/21_spatial_deconvolution.Rmd, section 22
# of analysis.html) into the H2 claims ledger. This is the project's first
# CELL-COMPOSITION readout -- every prior ledger row scores expression or
# activity; these score how MUCH of each microglial substate is present in the
# unsorted GeoMx tissue, deconvolved against the snRNAseq reference. The rows
# are a NEW evidence class (layer = "composition"): L-001 scores the
# amyloid-driven DAM-up / proliferative-down substate swap that exists in BOTH
# tau backgrounds (corroborating the amyloid-activation axis from tissue the
# sorted snRNAseq cannot resolve), and L-002 records the honest NON-finding
# that the tau x amyloid interaction does NOT surface as a composition effect
# (0 of 10 cell types at FDR<0.10; the DAM interaction is direction-concordant
# with the snRNAseq but non-significant). Both are graded SUGGESTIVE per the L5
# anti-double-count guard: L-001 scores the methodologically independent
# unsorted-tissue composition estimate, NOT the sorted-nuclei DAM expression DE
# already counted (section 2a / MG-M2 H2-023 / DAM_up H2-026); and both are
# held at Suggestive by the deconvolution's model-dependence (sorting-biased
# single-nucleus reference, no histology ground truth, profile9 substate
# collinearity kappa 35 -> only DAM + proliferative resolve cleanly, whole-
# tissue geometric ROIs with no plaque targeting -> regional-bulk composition
# not plaque-niche). Both are MARGIN-NEUTRAL by construction: L-001 supports
# BOTH amyloid contest models as a shared activation premise (mirroring K-001)
# and L-002 feeds NO model (empty supports + contradicts), so the locked
# contest margins (amyloid 18 / synaptic 12 / interaction 55) are unchanged.
# Numbering L-001..L-002 (Phase L; visually distinct from H2-*, I-*, J-*, K-*
# to flag Phase L origin in the ledger view).

phase_l <- rbind(
  row("L-001", "amyloid_activation", "composition",
      "Tissue-composition corroboration of amyloid microglial activation: GeoMx SpatialDecon two-stage deconvolution of the unsorted 91-ROI tissue shows an amyloid-driven DAM-up / proliferative-down microglial substate swap in BOTH tau backgrounds -- a new evidence class (cell-type abundance, not expression or activity) corroborating the amyloid-activation axis from tissue the sorted snRNAseq cannot resolve; scoped to the EXISTENCE of the composition shift, not the interaction.",
      "+", "within-microglia DAM fraction 0.18/0.18 (MAPTKI/P301S) -> 0.79/0.93 (NLGF_MAPTKI/NLGF_P301S); DAM abundance logFC +0.40 (nlgf_in_maptki, FDR 3.5e-5) / +0.53 (nlgf_in_p301s, FDR 2.4e-7); proliferative the mirror (logFC -0.45/-0.59, FDR<1e-3); pooled-microglia tissue fraction 0.37/0.45 -> 0.72/0.75 (absolute beta ~flat = relative redistribution)",
      "spatial_decon_contrasts.tsv (layer=stage2_substate, cell_type=Microglia_DAM, contrasts nlgf_in_maptki logFC +0.40 FDR 3.5e-5 + nlgf_in_p301s logFC +0.53 FDR 2.4e-7); spatial_decon_abundance_by_genotype.tsv (Microglia_DAM / Microglia_proliferative within-microglia fractions)",
      "section-22 GeoMx spatial-deconvolution chapter (22.1-22.3); spatial_decon_verdict.tsv (axis=amyloid_activation, verdict=corroborated); section-2a snRNAseq DAM expansion (the sorted-nuclei signal this UNSORTED-tissue estimate corroborates, NOT re-counts); MG-M2 DAM module H2-023; custom DAM_up gene set H2-026",
      1L, 1L, "Suggestive",
      "Hyp-1A;Hyp-1B;T-Inflammation", "",
      "Suggestive and NOT a re-count of the snRNAseq DAM expression DE (section 2a / MG-M2 H2-023 / DAM_up H2-026) -- it scores a methodologically independent tissue-composition estimate from the UNSORTED GeoMx tissue that the sorted snRNAseq cannot give, corroborating the amyloid-activation axis from a new evidence class. Scoped to EXISTENCE of the composition shift, not the interaction (see L-002 for the null composition interaction). Supports BOTH Hyp-1A and Hyp-1B (a shared amyloid-activation premise that does not discriminate the tau-modulation contest), so the amyloid margin is unchanged at 18 (mirrors the K-001 margin-neutral construction). Held at Suggestive because composition is a MODEL ESTIMATE: sorting-biased single-nucleus reference, no histology ground truth, profile9 substate collinearity (kappa 35, max column cor 0.99 -> Microglia_IFN + Microglia_homeostatic collinearity-suppressed to ~0, only DAM + proliferative trustworthy), and whole-tissue geometric ROIs with no plaque/morphology targeting (regional-bulk composition, not plaque-niche); 'proliferative' baseline dominance is plausibly a reference-profile artefact."),
  row("L-002", "interaction_metabolic", "composition",
      "Honest non-corroboration at the composition layer: the tau x amyloid interaction does NOT surface as a microglial-abundance composition effect -- 0 of 10 cell types reach FDR<0.10 at the interaction contrast, and the DAM interaction point estimate is direction-concordant with the snRNAseq (positive) but non-significant. A new-evidence-class negative result on the divergence axis, recorded for completeness, NOT a corroboration.",
      "mixed", "interaction: 0/10 cell types FDR<0.10 (best Astrocyte adj.P 0.25); Microglia_DAM interaction logFC +0.13, t=0.98, p=0.33, adj.P 0.66 (direction-concordant with the snRNAseq divergence, NS); Microglia_IFN + Microglia_homeostatic logFC = t = 0 at interaction (collinearity-suppressed)",
      "spatial_decon_contrasts.tsv (contrast=interaction across all 10 cell types in stage1_broad + stage2_substate; 0/10 sig; Microglia_DAM logFC +0.13 t=0.98 p=0.33 adj.P 0.66); spatial_decon_verdict.tsv (axis=interaction_metabolic, verdict=not corroborated)",
      "section-22 GeoMx spatial-deconvolution chapter (22.3, 22.5); the axis-3 snRNAseq / activity-layer interaction contest (H2-046..H2-063) this composition layer does NOT reproduce -- section 22 reports the DAM interaction is direction-concordant in sign but non-significant",
      1L, 1L, "Suggestive",
      "", "",
      "Suggestive and FEEDS NO MODEL (empty supports + contradicts) by construction -- an honest existence-of-non-finding, the discordance counterpart to L-001. The composition layer does NOT corroborate the tau x amyloid interaction at FDR 0.10 (0/10 cell types; the Microglia_DAM interaction +0.13 is direction-concordant with the snRNAseq divergence but p=0.33), so this row does NOT support Hyp-3B; and it is NOT scored against Hyp-3A because the point estimate is direction-concordant (not a positive no-interaction result) and double-counting must be avoided -- it therefore leaves all three contest margins unchanged (amyloid 18 / synaptic 12 / interaction 55). Recorded for completeness so the new evidence class reports what does NOT resolve, not only what does. Same composition caveats as L-001 (sorting-biased reference, no histology ground truth, kappa-35 substate collinearity, whole-tissue non-plaque ROIs); the interaction's non-significance is partly a power / estimability limit of the deconvolution, not proof of no biological interaction.")
)

# ---- PHASE M (microglial activation-trajectory dynamics) -------------------
#
# Two rows added at M5 (2026-06-06) propagating the microglial activation-
# trajectory dynamics layer (rmd/22_trajectory.Rmd, section 23 of analysis.html)
# into the H2 claims ledger. This is the project's first DYNAMICS / PROGRESSION
# readout -- every prior ledger row scores a STATIC quantity: expression DE,
# pathway / module / TF / kinase / LR / causal / SCENIC activity, or (Phase L)
# discrete cell COMPOSITION. These score how far along the homeostatic->DAM
# activation continuum each microglial cell sits, including within-cluster
# advancement that neither expression DE nor discrete composition can see. The
# rows are a NEW evidence class (layer = "dynamics"): M-001 scores the amyloid-
# driven homeostatic->DAM pseudotime advance that exists in BOTH tau backgrounds
# (corroborating the amyloid-activation axis as progression), and M-002 records
# the POSITIVE + significant tau x amyloid interaction on that progression --
# the project's FIRST positive interaction in an orthogonal (non-expression)
# layer, carried by the within-state PROGRESSION channel (not the discrete
# composition reshuffling that was null at L-002). Both are graded SUGGESTIVE
# per the M5 anti-double-count guard: M-001 scores the methodologically
# independent within-continuum progression, NOT the sorted-nuclei DAM expression
# DE already counted (section 2a / MG-M2 H2-023 / DAM_up H2-026) nor the
# discrete composition (L-001); M-002 is the progression-channel complement of
# the null composition interaction (L-002), NOT a re-count of the axis-3
# expression/activity interaction (H2-046..H2-063). Both are held at Suggestive
# by the dynamics caveats: pseudotime is a SINGLE INFERRED ORDERING partly
# confounded with composition (the root is not the most-potent state by any
# potency proxy -> activation ordering not developmental potency/time; 1 of 3
# cross-check methods, DPT, disagrees). Both are MARGIN-NEUTRAL (user-approved
# M5 gate, option B): M-001 supports BOTH amyloid contest models as a shared
# activation premise (mirroring K-001 / L-001), and M-002 supports only the
# T-Synergy THEME (not the Hyp-3A / Hyp-3B contest models) -- recording the
# positive interaction thematically without tipping the already-decided
# interaction contest, because the result is a single inferred ordering and is
# mechanism-agnostic (synergy in progression RATE, not Hyp-3B's specific
# Gsk3b / Myc / axon-guidance mechanism). The locked contest margins (amyloid
# 18 / synaptic 12 / interaction 55) are therefore unchanged. Numbering
# M-001..M-002 (Phase M; visually distinct from H2-*, I-*, J-*, K-*, L-* to
# flag Phase M origin in the ledger view).

phase_m <- rbind(
  row("M-001", "amyloid_activation", "dynamics",
      "Activation-dynamics corroboration of amyloid microglial activation: the snRNAseq microglial activation pseudotime (Slingshot homeostatic->DAM on the harmony embedding, root = homeostatic) advances under amyloid in BOTH tau backgrounds -- a new evidence class (within-continuum progression, not expression / activity / discrete composition) corroborating the amyloid-activation axis from a continuum the discrete layers cannot resolve; scoped to the amyloid MAIN effect, not the interaction.",
      "+", "mean homeostatic->DAM pseudotime logFC +7.87 (nlgf_in_maptki, FDR 8.9e-7) / +10.32 (nlgf_in_p301s, FDR 3.7e-8); fraction past the DAM-onset threshold 4%/4% (MAPTKI/P301S) -> 40%/59% (NLGF_MAPTKI/NLGF_P301S); tau alone -1.20 ns; amyloid main effect splits across BOTH channels (progression_cf +6.69/+8.99 sig AND composition_cf +1.76/+1.90 sig)",
      "trajectory_contrasts.tsv (measure=mean_pt, contrasts nlgf_in_maptki logFC +7.87 FDR 8.9e-7 + nlgf_in_p301s logFC +10.32 FDR 3.7e-8); trajectory_pseudotime_by_genotype.tsv (mean_pt + frac_past by genotype)",
      "section-23 trajectory chapter (23.1-23.3); trajectory_verdict.tsv (axis=amyloid_activation, verdict=corroborated); section-2a snRNAseq DAM expansion + section-22 spatial-decon DAM composition (L-001) -- the static signals this DYNAMICS estimate corroborates from a continuum the discrete layers cannot resolve, NOT re-counts; trajectory_progression_decomposition.tsv",
      1L, 1L, "Suggestive",
      "Hyp-1A;Hyp-1B;T-Inflammation", "",
      "Suggestive and NOT a re-count of the snRNAseq DAM expression DE (section 2a / MG-M2 H2-023 / DAM_up H2-026) or the section-22 DAM composition (L-001) -- it scores a methodologically independent DYNAMICS estimate: how far each cell sits along the homeostatic->DAM activation continuum, including within-cluster advancement that neither expression DE nor discrete composition can see. Scoped to the amyloid MAIN effect (the interaction is M-002). Supports BOTH Hyp-1A and Hyp-1B (a shared amyloid-activation premise that does not discriminate the tau-modulation contest), so the amyloid margin is unchanged at 18 (mirrors the K-001 / L-001 margin-neutral construction). Held at Suggestive because pseudotime is a SINGLE INFERRED ORDERING partly confounded with composition (guardrail #1): the asserted homeostatic root is NOT the most-potent state by any root-free proxy (entropy 3/4, n_genes 3/4, CytoTRACE2 4/4), so the axis is an ACTIVATION ordering not developmental potency/time; and 1 of 3 cross-check methods (scanpy DPT, rho -0.09) disagrees while the two curve/fate methods agree (nested Slingshot rho 0.98, CellRank2 rho 0.57)."),
  row("M-002", "interaction_metabolic", "dynamics",
      "Activation-dynamics corroboration of the tau x amyloid divergence: the tau (P301S) background supra-additively AMPLIFIES the amyloid-driven homeostatic->DAM progression -- the project's first POSITIVE + significant interaction in an orthogonal (non-expression) layer, and the decomposition shows it is carried by genuine within-state PROGRESSION (not the discrete composition reshuffling that was null at L-002); recorded thematically (T-Synergy) without moving the interaction contest, per the user-approved M5 gate.",
      "+", "interaction on mean pseudotime logFC +2.46 (FDR 0.077), carried by progression_cf +2.30 (FDR 0.077, ~94% of the total) NOT composition_cf +0.13 (FDR 0.81, null); fraction-past-DAM-threshold interaction +0.19 (FDR 0.013, the strongest); within_homeostatic interaction +2.68 (FDR 0.077, sig -- advancement of even homeostatic-labelled cells, which cannot be a composition artefact); CellRank2 DAM-fate interaction +0.05 ns (direction-concordant)",
      "trajectory_contrasts.tsv (contrast=interaction: mean_pt logFC +2.46 FDR 0.077, progression_cf +2.30 FDR 0.077, composition_cf +0.13 FDR 0.81, frac_past +0.19 FDR 0.013, within_homeostatic +2.68 FDR 0.077); trajectory_progression_decomposition.tsv (per-replicate composition vs progression channels)",
      "section-23 trajectory chapter (23.3-23.5); trajectory_verdict.tsv (axis=interaction_metabolic, verdict=corroborated as a progression effect); the section-22 composition interaction NULL (L-002) -- this DYNAMICS interaction is the progression-channel complement the discrete-composition layer could not resolve; trajectory_method_concordance.tsv (Slingshot vs CellRank2 rho 0.57; DPT the flagged outlier)",
      1L, 1L, "Suggestive",
      "T-Synergy", "",
      "Suggestive and MARGIN-NEUTRAL by the user-approved M5 gate (option B): it SUPPORTS only the T-Synergy THEME (a positive supra-additive tau x amyloid effect on activation progression) and is NOT scored against the Hyp-3A / Hyp-3B contest models, so the interaction margin is unchanged at 55. Rationale: although the result is positive AND significant (the project's first such interaction in a non-expression layer), it is a single inferred ordering and the trajectory is mechanism-agnostic -- it shows synergy in progression RATE, not the specific Gsk3b / Myc / axon-guidance mechanism that Hyp-3B names -- so a single Suggestive row is held to the theme rather than tipping an already-decided contest. NOT a re-count of the axis-3 expression/activity interaction (H2-046..H2-063): the load-bearing decomposition shows the effect is carried by within-state PROGRESSION (progression_cf +2.30 vs composition_cf +0.13 null), a channel neither the substate DE nor the discrete composition (L-002, itself null) captures. Held at Suggestive by the same dynamics caveats as M-001 (single inferred ordering, activation-not-potency root, DPT cross-method discordance, progression partly confounded with composition -- which is exactly why the decomposition is reported)."))

# ---- PHASE N (cross-cell-type specificity of the interaction) --------------
# Section 24 takes the project's first CROSS-CELL-TYPE readout (every prior row
# scores MICROGLIA only) and folds two margin-neutral rows into the ledger. The
# tau x amyloid factorial was refit in all six broad cell-type units (Astrocyte,
# Microglia-pooled, Neuronal, Oligodendrocyte, OPC, Vascular) x two estimators
# (pseudobulk limma-voom + single-cell NEBULA) x two power regimes (NATIVE all
# cells / MATCHED downsampled to K=289 cells per unit x replicate, the load-
# bearing regime that removes the ~20x cell-count confound). The rows are a NEW
# evidence class (layer = "specificity"). N-001 records that the amyloid MAIN
# effect is MICROGLIA-LED even under cell-count matching (validating the project's
# microglia focus for the amyloid-activation axis from a comparative lens the
# microglia-only ledger never had); N-002 records that the INTERACTION (the
# project's divergence readout) is NOT microglia-specific -- at MATCHED power the
# gene-level interaction collapses to ~0 in ALL six units incl Microglia (sole
# survivor Oligodendrocyte, non-microglial), and the broadest pathway interaction
# is Neuronal, not Microglia (depth-confounded). Both are graded SUGGESTIVE and
# MARGIN-NEUTRAL (user-approved N5 gate, option C -- feed BOTH halves): N-001
# supports BOTH amyloid contest models as a shared activation premise (mirroring
# K-001 / L-001 / M-001) so amyloid stays 18; N-002 FEEDS NO MODEL (empty
# supports + contradicts, mirroring L-002) -- it is NOT scored against Hyp-3A
# (the full-power microglia interaction the axis-3 contest rests on stays real;
# MATCHED collapse is a cross-cell-type COMPARABILITY statement at equalised
# power, not a retraction), NOT for Hyp-3B (mechanism-agnostic, about other cell
# types), and NOT for T-Synergy (the per-unit interaction profiles are mutually
# DISTINCT, matched rho 0.08-0.15, so the breadth is NOT a shared synergy
# programme -- tagging T-Synergy would overclaim cross-cell-type mechanistic
# unity, unlike M-002's genuine microglia supra-additive progression). The locked
# contest margins (amyloid 18 / synaptic 12 / interaction 55) are unchanged;
# N-001 lifts the T-Inflammation theme 30 -> 31 (breaking its tie with T-Synergy).
# Numbering N-001..N-002 (Phase N; visually distinct from H2-*, I-*, J-*, K-*,
# L-*, M-* to flag Phase N origin in the ledger view).

phase_n <- rbind(
  row("N-001", "amyloid_activation", "specificity",
      "Cross-cell-type validation of the microglia focus: at power-equalised resolution (MATCHED downsample to K=289 cells per unit x replicate) the amyloid MAIN effect produces more FDR<0.10 genes in microglia than in any other broad cell type, in BOTH tau backgrounds and BOTH estimators -- the first comparative confirmation, from a six-cell-type refit, that the amyloid response the project reads in microglia is genuinely microglia-led and not an artefact of microglia having been the only cell type examined; scoped to the amyloid MAIN effect (the interaction is N-002).",
      "+", "MATCHED (K=289/replicate, seed=1) amyloid main-effect gene counts FDR<0.10 are microglia-led across both backgrounds x both estimators: nlgf_in_maptki NEBULA Microglia 625 vs next-highest Astrocyte 218 (pseudobulk 232 vs 29); nlgf_in_p301s NEBULA Microglia 324 vs Astrocyte 210 (pseudobulk 406 vs 283); microglia lead DESPITE a depth disadvantage (median UMI/cell 1064 vs Neuronal 2809 = 2.64x), so the lead is conservative w.r.t. the depth confound; NATIVE all-cell tallies show the same microglia lead but are power-confounded by the ~20x cell-count spread",
      "storage/results/celltype_specificity_tally_matched.tsv (regime=matched, contrast=nlgf_in_maptki/nlgf_in_p301s, n_sig_10 by unit x estimator); storage/results/celltype_specificity_verdict.tsv (gene_interaction_tally row: main effects survive and stay microglia-led)",
      "section-24 cross-cell-type specificity chapter; section-2a snRNAseq microglial amyloid DE + the amyloid-activation axis-1 ledger (this cross-cell-type refit confirms those microglia signals are cell-type-led, NOT a re-count of them); section-22 DAM composition (L-001) + section-23 DAM pseudotime (M-001), the microglia-scoped amyloid corroborations whose SCOPE this comparative row validates",
      1L, 1L, "Suggestive",
      "Hyp-1A;Hyp-1B;T-Inflammation", "",
      "Suggestive and MARGIN-NEUTRAL (user-approved N5 gate, option C). NOT a re-count of the section-2a microglial amyloid DE or the K/L/M amyloid corroborations -- it scores a methodologically new CROSS-CELL-TYPE comparison (six-unit factorial refit) the microglia-only ledger could never make: that the amyloid response is microglia-LED relative to astrocytes / neurons / oligodendrocytes / OPCs / vascular cells under equal power. Supports BOTH Hyp-1A and Hyp-1B (a shared amyloid-activation premise that does not discriminate the tau-modulation contest), so the amyloid margin is unchanged at 18 (mirrors the K-001 / L-001 / M-001 margin-neutral construction); it does lift the T-Inflammation theme by one (30 -> 31, breaking its tie with T-Synergy). Held at Suggestive because (1) it is a comparative SCOPE-validation, not new inflammatory mechanism, and (2) MATCHED equalises CELL COUNT not transcriptome DEPTH -- microglia leading despite a 2.64x median-UMI disadvantage makes the lead conservative, but depth is uncontrolled. The microglia lead replicates across both tau backgrounds and both estimators (NEBULA + pseudobulk), so it is not an estimator artefact."),
  row("N-002", "interaction_metabolic", "specificity",
      "The tau x amyloid INTERACTION is NOT microglia-specific: at power-equalised resolution the gene-level interaction collapses to ~0 in ALL six broad cell types including microglia (sole survivor Oligodendrocyte, a NON-microglial type), the per-unit interaction profiles are mutually distinct rather than a shared programme, and the broadest pathway-level interaction response is Neuronal, not Microglia. The project's first cross-cell-type readout of the divergence axis -- it REFRAMES (does not refute) the microglia focus and is recorded for completeness, NOT as a contest move.",
      "mixed", "MATCHED (K=289/replicate) interaction genes FDR<0.10 = 0 in all six units incl Microglia; sole survivor Oligodendrocyte 6 (NEBULA only; pseudobulk 0) -- the only power-equalised gene-level interaction is non-microglial; the 4 NATIVE microglia interaction genes (Gm30211, Cd276, Cmah, Plac9) do NOT survive matching (power artefact of the ~20x cell-count spread); Microglia-vs-unit interaction logFC Spearman 0.08-0.15 (distinct profiles) while cross-estimator pb-vs-NEBULA agreement is 0.96-0.98 (cell types differ, not methods); GO-BP fGSEA interaction FDR<0.10 broadest in Neuronal (210 NEBULA / 242 pb) vs Microglia 2, enrichment overwhelmingly down (depth-confounded, neurons 2.64x UMI/cell)",
      "storage/results/celltype_specificity_tally_matched.tsv (contrast=interaction: 0/6 units except Oligodendrocyte NEBULA 6); storage/results/celltype_specificity_interaction_concordance.tsv (Microglia-vs-unit matched Spearman 0.08-0.15); storage/results/celltype_specificity_specificity_class.tsv (0 microglia-unique genes survive matching); storage/results/celltype_specificity_pathway_tally.tsv (interaction Neuronal >> Microglia); storage/results/celltype_specificity_verdict.tsv (OVERALL: interaction NOT microglia-specific)",
      "section-24 cross-cell-type specificity chapter; the axis-3 microglia interaction contest (H2-046..H2-063) whose SCOPE this row reframes -- N-002 does NOT re-count, refute, or re-score that contest: the MATCHED collapse is a cross-cell-type COMPARABILITY statement at K=289 equalised power, NOT a retraction of the full-power (all-cell) microglia interaction the contest rests on",
      1L, 1L, "Suggestive",
      "", "",
      "Suggestive and FEEDS NO MODEL (empty supports + contradicts) by construction -- the cross-cell-type counterpart to L-002, recorded for completeness so the new specificity evidence class reports the central divergence-axis finding, not only the reassuring N-001 main-effect result (anti-anchoring guardrail: feed both halves or neither; the user-approved N5 gate option C feeds both). It does NOT support Hyp-3B: N-002 is mechanism-agnostic and about NON-microglial cell types, not the specific Gsk3b / Myc / axon-guidance microglial mechanism Hyp-3B names. It is NOT scored against Hyp-3A: the interaction is NOT absent -- it is present broadly (the full-power microglia interaction the axis-3 contest rests on stays real; the MATCHED ~0 is a comparability statement at the severe K=289 downsample, not a positive no-interaction result), so scoring Hyp-3A would be a double-counting error. It does NOT support T-Synergy despite the interaction's breadth: the per-unit interaction profiles are mutually DISTINCT (matched Spearman 0.08-0.15), so the breadth is NOT a shared supra-additive programme -- tagging T-Synergy would overclaim cross-cell-type mechanistic unity, unlike M-002's genuine microglia within-state progression synergy. It therefore leaves all three contest margins unchanged (amyloid 18 / synaptic 12 / interaction 55) and the theme tallies untouched. Held at Suggestive by (1) the DEPTH confound (MATCHED equalises cell count not transcriptome depth; the Neuronal pathway dominance is read only as 'broadest outside microglia', never a quantitative neuronal-magnitude claim, because neurons carry 2.64x the UMIs/cell), (2) single-dataset snRNAseq scope, and (3) the K=289 matched power being deliberately severe (bottlenecked by the thinnest unit) so absence-of-signal is power-limited, not proof of no interaction."))

# ---- PHASE O (gene-level differential dynamics at the interaction) ---------
# One row added at O5 (2026-06-08) propagating the GENE-LEVEL differential-
# dynamics readout (rmd/22_trajectory.Rmd section 23, "Gene-level differential
# dynamics at the interaction" subsection) into the ledger. Phase M scored the
# AGGREGATE activation-trajectory interaction (per-replicate mean-pseudotime
# rate, M-002); O-001 is its PER-GENE resolution: a tradeSeq NB-GAM difference-
# of-differences Wald along the homeostatic->DAM microglial pseudotime finds 110
# genes (FDR<0.10) whose expression-trajectory SHAPE carries the tau x amyloid
# interaction -- a new evidence class (per-gene dynamics, layer = "dynamics_de")
# that the section-24 matched-power STATIC interaction DE cannot see (0/110
# static-significant in either estimator). It FEEDS NO MODEL (empty supports +
# contradicts, mirroring L-002 / N-002) per the user-approved O5 gate (option A1,
# always-margin-neutral) for two reinforcing reasons. (1) It is the gene-level
# COMPANION of M-002 (same tradeSeq fit, same trajectory), so scoring it to
# T-Synergy would double-count the theme M-002 already carries; recorded for
# completeness as a new evidence class, not as an independent interaction. (2)
# The pre-registered O5 overlap test -- locked BEFORE the genes were interpreted
# -- found NO enrichment for the Gsk3b / Myc / axon-guidance effectors Hyp-3B
# names (Gsk3b + Myc both ABSENT from the 110; canonical axon-guidance families
# 1/110 = Epha4 only, ~chance), so the contest-scoring IF-branch was not
# triggered and the ELSE (margin-neutral) branch applies. The programme is also
# quality-confounded -- 7/110 ambient / non-microglial-lineage transcripts and
# 20/110 ribosomal / translation genes (~13x over background, a recognised
# trajectory-DE confound) -- so it is held at Suggestive. The locked contest
# margins (amyloid 18 / synaptic 12 / interaction 55) and theme tallies are
# unchanged. Numbering O-001 (Phase O; visually distinct from H2-*, I-*..N-* to
# flag Phase O origin in the ledger view).

phase_o <- rbind(
  row("O-001", "interaction_metabolic", "dynamics_de",
      "Gene-level resolution of the activation-dynamics interaction: a tradeSeq NB-GAM difference-of-differences Wald along the homeostatic->DAM microglial pseudotime finds 110 genes (FDR<0.10) whose expression-trajectory SHAPE carries a tau x amyloid interaction -- the per-gene companion to M-002's aggregate progression-rate synergy, and a new evidence class (per-gene differential dynamics, not aggregate position) that the section-24 matched-power STATIC interaction DE cannot detect (0/110 static-significant). Recorded for completeness; feeds no model because the programme is quality-confounded and names none of Hyp-3B's effectors (see notes).",
      "mixed", "110 interaction genes FDR<0.10 (difference-of-differences Wald 14-105), nested within the conditionTest omnibus (1053) and associationTest (2994): association (2994) superset of omnibus (1053) superset of interaction (110); 0/110 reach FDR<0.10 in the section-24 MATCHED-power Microglia static interaction in EITHER estimator (best static FDR 0.72 NEBULA / >0.99 pseudobulk) -> non-redundant with static DE; top by Wald Ptgds, Enpp2, Bsg, Ly6c1, Ramp2, Scg5, Mal, Aplp1, Hbb-bt, Ctla2a; quality flags 7/110 ambient + 20/110 ribosomal/translation",
      "trajectory_dynamics_interaction.tsv (110 genes adj.P.Val<0.10; waldStat/df/effect_peak/effect_l2); trajectory_dynamics_vs_static.tsv (0/110 static-significant, NEBULA + pseudobulk matched-power, best static FDR 0.72/>0.99); trajectory_dynamics_omnibus.tsv (conditionTest 1053) + trajectory_dynamics_association.tsv (associationTest 2994)",
      "section-23 trajectory chapter, gene-level differential-dynamics subsection (the per-gene companion to M-002's aggregate mean-pseudotime interaction; trajectory_contrasts.tsv interaction mean_pt logFC +2.46 FDR 0.077, progression channel +2.30); M-002 itself -- O-001 is its gene-level resolution from the SAME tradeSeq fit, NOT an independent interaction",
      1L, 1L, "Suggestive",
      "", "",
      "Suggestive and FEEDS NO MODEL (empty supports + contradicts) by the user-approved O5 gate (option A1, always-margin-neutral). It is the GENE-LEVEL RESOLUTION of M-002's aggregate progression-rate interaction (same tradeSeq fit, same homeostatic->DAM trajectory), so scoring it to T-Synergy would DOUBLE-COUNT the theme M-002 already carries; recorded as a new evidence class (per-gene differential dynamics, layer dynamics_de) for completeness, not as an independent interaction. NOT scored against the Hyp-3A / Hyp-3B contest, per the pre-registered O5 overlap test (locked before the genes were interpreted): the 110-gene programme shows NO enrichment for the Gsk3b / Myc / axon-guidance effectors Hyp-3B names -- Gsk3b and Myc are both ABSENT, and canonical axon-guidance families contribute 1/110 (Epha4 only, vs ~0.4 expected = chance) -- so the IF-branch (score the contest) was not triggered and the margin-neutral ELSE-branch applies; the interaction margin is unchanged at 55. The programme is ALSO quality-confounded, holding it at Suggestive and reinforcing feeds-no-model: (1) 7/110 are ambient / non-microglial-lineage transcripts (haemoglobins Hbb-bt / Hba-a2, endothelial-vascular Ly6c1 / Ramp2 / Enpp2 / Igfbp7, myelin Mal) -- a genotype-varying ambient gradient along pseudotime reproduces this dynamic-interaction signature with no microglial biology; (2) its ONE coherent enrichment is ribosomal / translation genes (20/110, ~13x over the fitted background), a recognised trajectory-DE confound with cell-state and RNA content, NOT a microglial signalling mechanism. NON-REDUNDANT with static DE is the layer's positive contribution: 0/110 reach FDR<0.10 in the section-24 MATCHED-power Microglia static interaction in either estimator (best static FDR 0.72 NEBULA / >0.99 pseudobulk), so differential dynamics detects interactions the matched-power static contrast cannot -- but a non-redundant signal that is roughly one-quarter ambient+ribosomal and names none of Hyp-3B's effectors is recorded adjacent to M-002 and scored against nothing. Leaves all three contest margins (amyloid 18 / synaptic 12 / interaction 55) and the theme tallies untouched (anti-anchoring: the margin-neutral rule was locked before the genes were seen)."))

# ---- ASSEMBLE + WRITE ------------------------------------------------------

ledger <- rbind(axis1, axis2, axis3, cross_axis, phase_i, phase_j, phase_k, phase_l, phase_m, phase_n, phase_o)

# Validation: claim IDs are unique
stopifnot(!anyDuplicated(ledger$claim_id))

# Validation: supports_models + contradicts_models are subsets of ENTITY_IDS
validate_models <- function(s) {
  if (is.na(s) || s == "") return(TRUE)
  ids <- strsplit(s, ";")[[1]]
  all(ids %in% ENTITY_IDS)
}
stopifnot(all(vapply(ledger$supports_models,    validate_models, logical(1))))
stopifnot(all(vapply(ledger$contradicts_models, validate_models, logical(1))))

# Validation: confidence_grade is in {Strong, Moderate, Suggestive}
stopifnot(all(ledger$confidence_grade %in% c("Strong", "Moderate", "Suggestive")))

# Validation: direction is in {+, -, mixed}
stopifnot(all(ledger$direction %in% c("+", "-", "mixed")))

# Validation: axis is in the locked set
stopifnot(all(ledger$axis %in% c(
  "amyloid_activation", "synaptic_suppression", "interaction_metabolic",
  "cross_axis"
)))

# Validation: layer is in the locked set ("composition" added at Phase L for
# the GeoMx spatial-deconvolution cell-abundance evidence class; "dynamics"
# added at Phase M for the snRNAseq aggregate activation-pseudotime evidence
# class; "specificity" added at Phase N for the cross-cell-type factorial-refit
# class; "dynamics_de" added at Phase O for the per-gene tradeSeq differential-
# dynamics class -- the gene-level companion to the aggregate "dynamics" layer).
stopifnot(all(ledger$layer %in% c("tf", "kinase", "lr", "pathway", "cross_layer",
                                  "composition", "dynamics", "specificity",
                                  "dynamics_de")))

out_path <- file.path("storage", "results", "biological_model_claims_ledger.tsv")
write.table(ledger, out_path, sep = "\t", quote = FALSE, row.names = FALSE,
            na = "", fileEncoding = "UTF-8")

# Summary to stdout for log capture
cat(sprintf("Wrote %s with %d rows x %d cols.\n",
            out_path, nrow(ledger), ncol(ledger)))
cat(sprintf("  Confidence grade: Strong=%d / Moderate=%d / Suggestive=%d\n",
            sum(ledger$confidence_grade == "Strong"),
            sum(ledger$confidence_grade == "Moderate"),
            sum(ledger$confidence_grade == "Suggestive")))
cat("  Per-axis row count:\n")
print(table(ledger$axis))
cat("  Per-layer row count:\n")
print(table(ledger$layer))
cat("  Support counts per entity:\n")
support_counts <- sapply(ENTITY_IDS, function(eid) {
  sum(vapply(ledger$supports_models, function(s) {
    if (is.na(s) || s == "") return(FALSE)
    eid %in% strsplit(s, ";")[[1]]
  }, logical(1)))
})
print(support_counts)
cat("  Contradict counts per entity:\n")
contradict_counts <- sapply(ENTITY_IDS, function(eid) {
  sum(vapply(ledger$contradicts_models, function(s) {
    if (is.na(s) || s == "") return(FALSE)
    eid %in% strsplit(s, ";")[[1]]
  }, logical(1)))
})
print(contradict_counts)
cat(sprintf("  Net support per entity (support - contradict):\n"))
print(support_counts - contradict_counts)

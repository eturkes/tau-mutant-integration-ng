#!/usr/bin/env Rscript
# build_capstone_synthesis.R
#
# Arc P (capstone) deliverable C2: re-aggregate the eleven completed evidence
# arcs (D->O) + the four base modalities + the human translational layer into
# TWO legible synthesis TSVs. PURE re-aggregation / presentation of existing
# ledger + verdict TSVs -- NO new biological computation. Read-only: it never
# edits the ledger / 11-entity set / grade defs / contest margins / analysis.Rmd
# section 17 (storage/notes/capstone_synthesis_plan.md principle #5).
#
# Inputs (all in storage/results/ unless noted):
#   biological_model_contest_verdicts.tsv  (3 per-axis contests; winner+margin)
#   biological_model_adjudication.tsv       (11 entities; winner support counts)
#   biological_model_claims_ledger.tsv      (92 rows; I/J/K/L/M/N/O feeds + layers)
#   tf_activity_verdict.tsv, kinase_activity_verdict.tsv, ccc_lr_verdict.tsv,
#   causal_network_verdict.tsv, scenic_verdict.tsv, spatial_decon_verdict.tsv,
#   trajectory_verdict.tsv, celltype_specificity_verdict.tsv   (8 per-arc verdicts)
#   pathway_survey_unified_leaderboard.tsv  (arc D; 144 rows)
#   trajectory_dynamics_interaction.tsv     (arc O; 11,478 fitted, 110 FDR<0.10)
#   storage/cache/integration_table.rds     (cross-modality DE; 107,825 x 14)
#   storage/cache/summary_human_validation.rds (human layer; SEPARATE band)
#
# Outputs (storage/results/):
#   capstone_convergence_matrix.tsv  -- long: evidence_class, arc, band, axis,
#     cell_value, sign, grade, source, note  (13 classes x 3 axes = 39 cells).
#   capstone_contest_summary.tsv     -- contest, hyp_A, hyp_B, winner, margin,
#     n_strong, n_moderate, n_suggestive, n_layers, n_modalities  (3 contests).
#
# The convergence matrix is a CURATED arc x axis mapping (the verdict TSVs are
# not all natively "signed per axis"): each cell's signed call is taken from the
# completed-plan DIGEST verdict reads, and each cell carries a cited `source`
# (verdict-TSV path + value, or ledger claim_id) -- principle #2, zero uncited
# cells. The contest summary is re-surfaced straight from the adjudication +
# contest_verdicts outputs (the contest view); n_layers / n_modalities are light
# descriptive aggregates over the ledger rows supporting each winner.
#
# Two cell-vocabulary notes (curation latitude granted to C2 by the plan):
#   - `mixed`  meta-tag added for genuinely sign-divergent cells (modalities /
#     LR pairs / difference-of-differences disagree on direction). Forcing these
#     to a signed token would over-state convergence (the guardrail's primary
#     failure mode); forcing them to `null` would hide that the layer IS
#     significant. The ledger itself uses direction="mixed", so this is in-house.
#   - `band` column added (mouse|human) so the figure (C3) can hold the human
#     row in its own separated band -- principle #4 made machine-readable.
#
# `sign` frame (documented so the figure + readers interpret cells correctly):
#   effect-direction arcs (Expression, Pathway, TF, Kinase, LR, NF-kB,
#     Composition, Dynamics): + up/activated/advanced, - down/suppressed/
#     attenuated, mixed sign-divergent, 0 tested-null, na not-probed.
#   existence/recovery arcs (Causal, SCENIC, Specificity, Dynamics-DE): + the
#     prior signal is recovered / a non-redundant layer exists, 0 not recovered,
#     reframes the arc reframes (not a directional axis call), na not-probed.
# Per-cell `note` makes each frame explicit; no cell relies on the frame alone.

suppressPackageStartupMessages({
  library(readr)
})

RES <- "storage/results"
CACHE <- "storage/cache"

# ---- 0. local writer (script-local, mirrors R/utils.R::write_tsv_safe) -------
write_tsv_local <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_tsv(x, path)
  Sys.chmod(path, mode = "0664")
}

rp <- function(...) file.path(RES, ...)
need <- function(path) {
  if (!file.exists(path)) stop("missing required input: ", path, call. = FALSE)
  path
}

# ---- 1. load the contest view (adjudication + contest_verdicts) --------------
contest <- read.delim(need(rp("biological_model_contest_verdicts.tsv")),
                      stringsAsFactors = FALSE)
adj     <- read.delim(need(rp("biological_model_adjudication.tsv")),
                      stringsAsFactors = FALSE)
ledger  <- read.delim(need(rp("biological_model_claims_ledger.tsv")),
                      stringsAsFactors = FALSE)

# n_layers / n_modalities = light descriptive aggregates over the ledger rows
# that SUPPORT each contest winner (supports_models is ";"-joined; exact match
# after split so Hyp-1B never matches on Hyp-1A). NOT a re-adjudication -- the
# winner + margin + grade counts come verbatim from the loaded TSVs.
supports_winner <- function(winner_id) {
  hits <- vapply(strsplit(ledger$supports_models, ";"), function(v)
    winner_id %in% trimws(v), logical(1))
  ledger[hits, , drop = FALSE]
}

contest_summary <- do.call(rbind, lapply(seq_len(nrow(contest)), function(i) {
  r <- contest[i, ]
  winner <- r$favoured_model
  a <- adj[adj$entity_id == winner, ]
  sup <- supports_winner(winner)
  data.frame(
    contest        = r$contest_id,
    hyp_A          = sprintf("%s (%s)", r$model_a_id, r$model_a_name),
    hyp_B          = sprintf("%s (%s)", r$model_b_id, r$model_b_name),
    winner         = sprintf("%s (%s)", winner,
                             ifelse(winner == r$model_a_id, r$model_a_name,
                                    r$model_b_name)),
    margin         = r$favoured_by_margin,
    n_strong       = a$n_strong_supports,
    n_moderate     = a$n_moderate_supports,
    n_suggestive   = a$n_suggestive_supports,
    n_layers       = length(unique(sup$layer)),
    n_modalities   = if (nrow(sup)) max(sup$n_replicates_in_modalities,
                                        na.rm = TRUE) else 0L,
    stringsAsFactors = FALSE
  )
}))

# guard: the three locked margins must be exactly 18 / 12 / 55 (sealed section 17)
locked_margins <- c(amyloid_activation = 18L, synaptic_suppression = 12L,
                    interaction_metabolic = 55L)
for (cn in names(locked_margins)) {
  got <- contest_summary$margin[contest_summary$contest == cn]
  if (length(got) != 1L || got != locked_margins[[cn]])
    stop(sprintf("contest margin drift: %s expected %d got %s -- the sealed sec-17 verdicts moved",
                 cn, locked_margins[[cn]], paste(got, collapse = ",")), call. = FALSE)
}

# ---- 2. load the per-arc verdict TSVs (provenance + drift guard) --------------
vpath <- list(
  tf       = need(rp("tf_activity_verdict.tsv")),
  kinase   = need(rp("kinase_activity_verdict.tsv")),
  ccc      = need(rp("ccc_lr_verdict.tsv")),
  causal   = need(rp("causal_network_verdict.tsv")),
  scenic   = need(rp("scenic_verdict.tsv")),
  spatial  = need(rp("spatial_decon_verdict.tsv")),
  traj     = need(rp("trajectory_verdict.tsv")),
  spec     = need(rp("celltype_specificity_verdict.tsv")),
  pathway  = need(rp("pathway_survey_unified_leaderboard.tsv")),
  trajdyn  = need(rp("trajectory_dynamics_interaction.tsv"))
)
invisible(need(file.path(CACHE, "integration_table.rds")))
invisible(need(file.path(CACHE, "summary_human_validation.rds")))

# Anchor the most load-bearing cells to a live substring in their source file,
# so a future verdict-TSV edit that drops the cited value fails this build.
anchor <- function(key, ...) {
  txt <- paste(readLines(vpath[[key]], warn = FALSE), collapse = "\n")
  for (s in c(...)) if (!grepl(s, txt, fixed = TRUE))
    stop(sprintf("anchor '%s' absent from %s -- source drift", s, vpath[[key]]),
         call. = FALSE)
}
anchor("tf",      "Spi1", "Myc")
anchor("kinase",  "padj 0.002", "Gsk3b")
anchor("ccc",     "rank 775", "Apoe_Trem2")
anchor("causal",  "Gsk3b", "Mapk14")
anchor("scenic",  "Spi1 recovered", "expression-floor")
anchor("spatial", "3.5e-05", "2.4e-07")
anchor("traj",    "8.9e-07", "+2.46")
anchor("spec",    "Microglia 625", "not microglia-specific")

# arc-O: recompute the "110 of 11,466 FDR<0.10" anchor from the live cache
trajdyn <- read.delim(vpath$trajdyn, stringsAsFactors = FALSE)
n_dyn_sig <- sum(trajdyn$adj.P.Val < 0.10, na.rm = TRUE)
if (n_dyn_sig != 110L)
  stop("arc-O interaction-dynamics gene count drifted: expected 110, got ",
       n_dyn_sig, call. = FALSE)

# ---- 3. the curated convergence matrix (13 evidence classes x 3 axes) ---------
AX <- c("amyloid_activation", "synaptic_suppression", "interaction_metabolic")
.cells <- list()
cell <- function(evidence_class, arc, band, axis, cell_value, sign, grade,
                 source, note) {
  .cells[[length(.cells) + 1L]] <<- data.frame(
    evidence_class = evidence_class, arc = arc, band = band, axis = axis,
    cell_value = cell_value, sign = sign, grade = grade,
    source = source, note = note, stringsAsFactors = FALSE)
}

# --- 1. Expression DE (4 modalities) -- base, foundation layer ---
cell("Expression DE (4 modalities)", "base", "mouse", AX[1], "strong+", "+", "Strong",
     "integration_table.rds (107,825 gene x contrast cross-modality DE); pathway_survey axis-1 verdict (DAM/MHC-II/inflammatory up at both NLGF)",
     "Foundation layer: amyloid-driven DAM / inflammatory up-regulation across modalities; the strongest, most multi-layer axis (contest Hyp-1B net +34).")
cell("Expression DE (4 modalities)", "base", "mouse", AX[2], "mod-", "-", "Moderate",
     "integration_table.rds; pathway_survey axis-2 verdict (synaptic GO CC/BP consistent - at both NLGF)",
     "Synaptic-gene down-regulation at both NLGF contrasts -- the synaptic-suppression axis at the expression layer.")
cell("Expression DE (4 modalities)", "base", "mouse", AX[3], "mixed", "mixed", "Moderate",
     "pathway_survey_unified_leaderboard.tsv: MG-M3 interaction:3/mixed (leader_score 23.4); integration_table.rds",
     "Interaction present at the expression layer but SIGN-DIVERGENT across modalities (OXPHOS/ETC + ribosomal + MG-M3 + DAM_up; 3 modalities FDR<0.05, signs disagree).")

# --- 2. Pathway / module -- arc D ---
cell("Pathway / module", "D", "mouse", AX[1], "strong+", "+", "Strong",
     "pathway_survey_unified_leaderboard.tsv (MG-M2, MHC-II antigen-presentation, GOBP_ADAPTIVE_IMMUNE_RESPONSE, DAM_up/WAM/HAM_disease all +)",
     "MG-M2 is the single cleanest cross-modality amyloid hit.")
cell("Pathway / module", "D", "mouse", AX[2], "mod-", "-", "Moderate",
     "pathway_survey_unified_leaderboard.tsv (synaptic GO CC/BP terms consistent - at both NLGF)",
     "Synaptic-compartment GO terms down at NLGF.")
cell("Pathway / module", "D", "mouse", AX[3], "mixed", "mixed", "Moderate",
     "pathway_survey_unified_leaderboard.tsv: MG-M3 tops leader_score 23.4 (interaction:3/mixed); GOCC_NADH_DEHYDROGENASE_COMPLEX interaction:3/mixed",
     "Top leader-score pathway (MG-M3) is mixed-sign at interaction; metabolic / translational families likewise mixed.")

# --- 3. TF activity -- arc E ---
cell("TF activity", "E", "mouse", AX[1], "strong+", "+", "Strong",
     "tf_activity_verdict.tsv (amyloid_activation): Spi1/Nfkb1/Sp3 (+), A0A979HLR9; 2 of top 5 reach the cross-modality leader rule",
     "DAM-programme TFs activated at NLGF onset (Spi1=PU.1, master microglial TF + AD GWAS locus).")
cell("TF activity", "E", "mouse", AX[2], "null", "0", NA_character_,
     "tf_activity_verdict.tsv (synaptic_suppression): honest non-finding; top TFs identical to axis-1 by Interpretation-A entanglement",
     "No distinct synaptic TF signature -- consistent with an engulfment / post-transcriptional mechanism (not regulator-driven).")
cell("TF activity", "E", "mouse", AX[3], "mod-", "-", "Moderate",
     "tf_activity_verdict.tsv (interaction_metabolic): Myc/Creb1/Tp53/Jun (-), signed mean [-4.85,-1.81] (Myc interaction -6.60, scenic_verdict.tsv)",
     "Coherent NEGATIVE: biosynthetic / proliferation TFs (Myc-led) suppressed at the tau x amyloid interface.")

# --- 4. Kinase activity -- arc F (bulk hippocampus phospho, NOT microglia-sorted) ---
cell("Kinase activity", "F", "mouse", AX[1], "mod+", "+", "Moderate",
     "kinase_activity_verdict.tsv (amyloid_activation): Cdk2/Cdk5/Csnk1e/Mapk8 (+); Cdk5/Csnk1e FDR<0.10 at nlgf_in_p301s",
     "Proliferation + stress-kinase ensemble. Bulk-hippocampus phospho (NOT microglia-sorted).")
cell("Kinase activity", "F", "mouse", AX[2], "null", "0", NA_character_,
     "kinase_activity_verdict.tsv (synaptic_suppression): honest non-finding; identical to axis-1 by Interpretation-A entanglement",
     "No distinct synaptic kinase signature. Bulk-hippocampus phospho.")
cell("Kinase activity", "F", "mouse", AX[3], "mod+", "+", "Moderate",
     "kinase_activity_verdict.tsv (interaction_metabolic): Gsk3b padj 0.002 (strongest single mechanism-arc result), Mapk14/Cdk5/Cdk1 (+)",
     "Canonical tau-kinase ensemble, Gsk3b lead -- the interaction's standout single result. Bulk-phospho seed caveat.")

# --- 5. CCC / ligand-receptor -- arc G ---
cell("CCC / ligand-receptor", "G", "mouse", AX[1], "mod+", "+", "Moderate",
     "ccc_lr_verdict.tsv (amyloid_activation): 55/100 top-100 microglia-involving; Cd200_Cd200r1 / Apoe_Trem2 / App_Cd74 DAM-programme LR (rank 6-25)",
     "Cell-type-resolved DAM-programme intercellular detail; top-5 by magnitude are neuronal axon-guidance.")
cell("CCC / ligand-receptor", "G", "mouse", AX[2], "mod+", "+", "Moderate",
     "ccc_lr_verdict.tsv (synaptic_suppression): Apoe_Trem2 rank 9/10, App_Cd74 rank 12 (UP); triangulation VINDICATED; classical complement C1qb_Lrp1 rank 775",
     "The DECISIVE synaptic-suppression arc (Hyp-2B wins by 12): suppression via UP-regulated TREM2 / APP-fragment CLEARANCE ligands, NOT classical complement. (+ = clearance ligands gained, the mechanism OF suppression.)")
cell("CCC / ligand-receptor", "G", "mouse", AX[3], "mixed", "mixed", "Moderate",
     "ccc_lr_verdict.tsv (interaction_metabolic): 3+/2- top-5; Efna5_Epha4 + / L1cam_Cd9 - at broader scope",
     "Axon-guidance / synaptic-adhesion intercellular signalling REWIRED at the interface (not uniformly suppressed); part of the most-integrated 3-layer interaction reading.")

# --- 6. NF-kB attenuation -- arc I (the Hyp-1B / T-Tau-attenuates mechanism) ---
cell("NF-kB attenuation", "I", "mouse", AX[1], "strong-", "-", "Strong",
     "ledger I-001..I-006 (4 Strong / 2 Moderate; supports Hyp-1B;T-Tau-attenuates, contradicts Hyp-1A); per_state_tf_activity.tsv; per_state_nfkb_target_gsea.tsv",
     "Global (all 4 substates) tau-driven suppression of amyloid-driven NF-kB TRANSCRIPTION -- the mechanism behind Hyp-1B's win (amyloid margin 6 -> 18). A suppressive modifier WITHIN the amyloid axis.")
cell("NF-kB attenuation", "I", "mouse", AX[2], "n.a.", "na", NA_character_,
     "n/a -- NF-kB attenuation does not probe the synaptic-suppression axis",
     "Axis not probed.")
cell("NF-kB attenuation", "I", "mouse", AX[3], "n.a.", "na", NA_character_,
     "ledger I-rows are axis=cross_axis, booked to the amyloid contest (sec 17)",
     "The tau-modulation of amyloid-driven NF-kB is interaction-flavoured but adjudicated under the amyloid contest, not separately scored at interaction. The data-driven SCENIC layer (K) did NOT reproduce it (scenic_verdict.tsv nfkb_family_sec18: Rel ns at interaction in all 4 substates) -- an honest non-corroboration.")

# --- 7. Causal topology -- arc J (existence/recovery frame) ---
cell("Causal topology", "J", "mouse", AX[1], "sugg+", "+", "Suggestive",
     "causal_network_verdict.tsv (amyloid_activation): 18/22 TFs recovered at nlgf_in_p301s; Gsk3b -> Nfkb1, Mapk14 -> Gata2 -| Spi1; topology_strength 'strong'",
     "Amyloid TF wiring recovered (topology 'strong'); ledger feed Suggestive (anti-double-count of the endpoint activities). Bulk-phospho seed.")
cell("Causal topology", "J", "mouse", AX[2], "null", "0", NA_character_,
     "causal_network_verdict.tsv (synaptic_suppression): shares nlgf_in_p301s net, no distinct wiring; 'non-finding (shared/entangled)'",
     "No synaptic-specific wiring (the same non-finding as the TF / kinase layers).")
cell("Causal topology", "J", "mouse", AX[3], "sugg+", "+", "Suggestive",
     "causal_network_verdict.tsv (interaction_metabolic): Gsk3b -| Mapk14 -> Myc/Foxo3 (5/5 TFs); ledger J-001 (supports Hyp-3B;T-Synergy, contradicts Hyp-3A)",
     "The headline interaction BRIDGE recovered -- wires the sec-15 lead kinase (Gsk3b) to the sec-14 lead TF (Myc); interaction margin 53 -> 55. Signed edges (-|/->) net consistent with Myc suppression.")

# --- 8. SCENIC regulons -- arc K (data-driven; existence/recovery frame) ---
cell("SCENIC regulons", "K", "mouse", AX[1], "sugg+", "+", "Suggestive",
     "scenic_verdict.tsv (amyloid_activation): Spi1 recovered +10/10 (15 targets); ledger K-001/K-002 (Hyp-1A;Hyp-1B;T-Inflammation / T-Inflammation)",
     "Data-driven (motif + co-expression) recovery of Spi1 (amyloid) & Rel (NF-kB) regulons -- EXISTENCE-level corroboration; activity-level DISCORDANT (Spi1 +2.73* -> -0.43 ns); margin-neutral.")
cell("SCENIC regulons", "K", "mouse", AX[2], "n.a.", "na", NA_character_,
     "n/a -- the SCENIC verdict has no synaptic-suppression axis row (TF-regulon focused)",
     "Axis not probed (the data-driven layer tests the prior amyloid / interaction / NF-kB TFs).")
cell("SCENIC regulons", "K", "mouse", AX[3], "null", "0", NA_character_,
     "scenic_verdict.tsv (interaction_metabolic): none recovered; Myc & Tp53 expression-floor (structural); 0/51 regulons sig at interaction",
     "Half the interaction TF axis (Myc/Tp53) is below the >=1% expression floor -> STRUCTURAL non-recovery; the data-driven interaction signal lives in NOVEL regulons (Ets1/Klf2/Hif1a/Zbtb20), not the prior TFs. The prior interaction mechanism is NOT recovered de novo.")

# --- 9. Composition -- arc L (first tissue-composition class) ---
cell("Composition", "L", "mouse", AX[1], "sugg+", "+", "Suggestive",
     "spatial_decon_verdict.tsv (amyloid_activation): DAM logFC +0.40 (FDR 3.5e-05) / +0.53 (FDR 2.4e-07); proliferative the mirror; ledger L-001",
     "First TISSUE-COMPOSITION readout; amyloid-driven DAM up / proliferative down substate swap in unsorted tissue; corroborates amyloid activation.")
cell("Composition", "L", "mouse", AX[2], "n.a.", "na", NA_character_,
     "spatial_decon_verdict.tsv (synaptic_suppression): no cell-type-abundance endpoint for the synaptic/engulfment programme",
     "Composition has no synapse-specific endpoint (neuronal abundance falls as broad dilution, not a synapse-specific claim).")
cell("Composition", "L", "mouse", AX[3], "null", "0", NA_character_,
     "spatial_decon_verdict.tsv (interaction_metabolic): 0/10 cell types sig; DAM interaction +0.13 (t=0.98, p=0.33); ledger L-002 (feeds no model)",
     "Composition does NOT corroborate the interaction: direction-concordant (+0.13) but NULL (p=0.33). A first-class null -- the interaction is not a composition effect.")

# --- 10. Dynamics / progression -- arc M (first dynamics class; THE positive interaction) ---
cell("Dynamics / progression", "M", "mouse", AX[1], "sugg+", "+", "Suggestive",
     "trajectory_verdict.tsv (amyloid_activation): pseudotime logFC +7.87 (FDR 8.9e-07) / +10.32 (FDR 3.7e-08); frac-past 4% -> 40%/59%; ledger M-001",
     "First DYNAMICS / progression readout; amyloid advances homeostatic -> DAM pseudotime in both tau backgrounds.")
cell("Dynamics / progression", "M", "mouse", AX[2], "n.a.", "na", NA_character_,
     "trajectory_verdict.tsv (synaptic_suppression): no trajectory endpoint (synaptic axis is not a microglial activation state)",
     "Axis not probed.")
cell("Dynamics / progression", "M", "mouse", AX[3], "sugg+", "+", "Suggestive",
     "trajectory_verdict.tsv (interaction_metabolic): interaction +2.46 (FDR 0.077) carried by PROGRESSION (+2.30) not composition (+0.13, FDR 0.81); frac-past +0.19 (FDR 0.013); ledger M-002 (T-Synergy)",
     "THE project's first POSITIVE + significant orthogonal interaction; carried by progression NOT composition (non-redundant with the arc-L null); supports T-Synergy (margin-neutral). The load-bearing interaction nuance.")

# --- 11. Specificity -- arc N (cross-cell-type; reframes the microglia focus) ---
cell("Specificity", "N", "mouse", AX[1], "sugg+", "+", "Suggestive",
     "celltype_specificity_verdict.tsv (R1/R2 matched): NEBULA Microglia 625 vs Astrocyte 218 (nlgf_in_maptki); ledger N-001",
     "Amyloid MAIN effect IS microglia-led at MATCHED K=289 power (survives the 2.64x depth disadvantage).")
cell("Specificity", "N", "mouse", AX[2], "n.a.", "na", NA_character_,
     "celltype_specificity probed amyloid main-effect + interaction tallies, not a synaptic-suppression axis",
     "Axis not probed.")
cell("Specificity", "N", "mouse", AX[3], "reframes", "reframes", "Suggestive",
     "celltype_specificity_verdict.tsv (OVERALL): interaction FDR<0.10 = 0 in ALL 6 units incl Microglia at matched; Microglia-vs-unit rho 0.08-0.15 (distinct); ledger N-002 (feeds no model)",
     "REFRAMES the microglia focus: the tau x amyloid interaction is NOT microglia-specific (collapses to ~0 at matched power in every cell type incl microglia; per-unit profiles distinct, not a shared programme). A specificity null, not a refutation of the interaction itself.")

# --- 12. Dynamics-DE -- arc O (gene-level companion of M-002; existence frame) ---
cell("Dynamics-DE", "O", "mouse", AX[1], "n.a.", "na", NA_character_,
     "trajectory_dynamics -- the scored deliverable is the interaction difference-of-differences contrast; amyloid main effect is not this arc's endpoint",
     "Axis not the arc's deliverable (omnibus / association cover all conditions; the scored contrast is the interaction).")
cell("Dynamics-DE", "O", "mouse", AX[2], "n.a.", "na", NA_character_,
     "trajectory_dynamics -- no synaptic endpoint",
     "Axis not probed.")
cell("Dynamics-DE", "O", "mouse", AX[3], "sugg+", "+", "Suggestive",
     sprintf("trajectory_dynamics_interaction.tsv (%d of 11,466 genes FDR<0.10); trajectory_dynamics_vs_static.tsv (0/%d static-sig at matched power); ledger O-001 (feeds no model)",
             n_dyn_sig, n_dyn_sig),
     sprintf("%d gene-level interaction-DYNAMICS genes, ALL static-null at matched power -> a NON-REDUNDANT layer + a candidate gene-level signature for M-002. Direction-AGNOSTIC (difference-of-differences omnibus); ~1/4 ambient + ribosomal-confounded; Gsk3b/Myc ABSENT (axon-guidance ~ chance) -> feeds-no-model (M-002 companion, avoids double-count).",
             n_dyn_sig))

# --- 13. Human cross-species -- SEPARATE band (evidence-class-incommensurable, H7) ---
cell("Human cross-species", "human", "human", AX[1], "sugg+", "+", "Suggestive",
     "summary_human_validation.rds: conservation gate PASS (4/4 SEA-AD substates; monotone -SEAAD occupancy 0.150 -> 0.490); DAM-up +amyloid 5/5 strata",
     "SEPARATE BAND (NEVER counted in the mouse contests, H7). Directional translational consistency: the amyloid / DAM structure conserves in human AD (SEA-AD MTG+DLPFC).")
cell("Human cross-species", "human", "human", AX[2], "n.a.", "na", NA_character_,
     "summary_human_validation.rds: the human layer tested conservation + the interaction-localised mechanisms, not a synaptic-suppression axis",
     "SEPARATE BAND. Axis not directly tested.")
cell("Human cross-species", "human", "human", AX[3], "null", "mixed", "Suggestive",
     "summary_human_validation.rds: pre-reg mechanisms match the mouse amyloid:tau SIGN in 5/5 strata (NFKB_union -, MG_M3 -, Gsk3b +) but NO interaction survives FDR (Spearman ~0.65 collinearity, VIF 1.5-1.8)",
     "SEPARATE BAND. Direction-of-effect CONCORDANT in 5/5 strata (the locked primary readout: direction beats significance) but FDR-null -- under-identification, NOT counter-evidence. A first-class human-interaction null with directional concordance.")

matrix_tbl <- do.call(rbind, .cells)
rownames(matrix_tbl) <- NULL

# ---- 4. validate the matrix --------------------------------------------------
stopifnot(
  "expect 13 evidence classes x 3 axes = 39 cells" =
    nrow(matrix_tbl) == 39L,
  "expect exactly 13 distinct evidence classes" =
    length(unique(matrix_tbl$evidence_class)) == 13L,
  "every evidence class must carry all 3 axes" =
    all(table(matrix_tbl$evidence_class) == 3L),
  "axis vocabulary" = all(matrix_tbl$axis %in% AX),
  "band vocabulary" = all(matrix_tbl$band %in% c("mouse", "human")),
  "every cell carries a non-empty source" =
    all(nzchar(trimws(matrix_tbl$source))),
  "every cell carries a non-empty note" =
    all(nzchar(trimws(matrix_tbl$note))),
  "cell_value vocabulary" = all(matrix_tbl$cell_value %in%
    c("strong+", "mod+", "sugg+", "null", "sugg-", "mod-", "strong-",
      "mixed", "reframes", "n.a.")),
  "sign vocabulary" = all(matrix_tbl$sign %in% c("+", "-", "mixed", "0", "na", "reframes")),
  "grade vocabulary" = all(is.na(matrix_tbl$grade) |
    matrix_tbl$grade %in% c("Strong", "Moderate", "Suggestive")),
  "human is the only non-mouse band" =
    all(matrix_tbl$band[matrix_tbl$evidence_class == "Human cross-species"] == "human") &&
    all(matrix_tbl$band[matrix_tbl$evidence_class != "Human cross-species"] == "mouse")
)

# directional sanity vs the DIGEST verdict reads (anti-anchoring: nulls are
# first-class, convergence must not be over-stated)
amy <- matrix_tbl[matrix_tbl$axis == "amyloid_activation", ]
intx <- matrix_tbl[matrix_tbl$axis == "interaction_metabolic", ]
stopifnot(
  "amyloid axis must be predominantly positive (the strong convergent axis)" =
    sum(amy$sign == "+") >= 8L,
  "the NF-kB attenuation cell must be the strong- amyloid modifier" =
    matrix_tbl$cell_value[matrix_tbl$evidence_class == "NF-kB attenuation" &
                          matrix_tbl$axis == "amyloid_activation"] == "strong-",
  "interaction must carry >=3 explicit null/reframes cells (nulls first-class)" =
    sum(intx$cell_value %in% c("null", "reframes")) >= 3L,
  "interaction must NOT be uniformly positive (no over-stated convergence)" =
    any(intx$cell_value %in% c("null", "reframes", "mixed"))
)

# ---- 5. emit -----------------------------------------------------------------
out_matrix  <- rp("capstone_convergence_matrix.tsv")
out_contest <- rp("capstone_contest_summary.tsv")
write_tsv_local(matrix_tbl, out_matrix)
write_tsv_local(contest_summary, out_contest)

cat(sprintf("[capstone] wrote %s (%d cells x %d cols)\n",
            out_matrix, nrow(matrix_tbl), ncol(matrix_tbl)))
cat(sprintf("[capstone] wrote %s (%d contests x %d cols)\n",
            out_contest, nrow(contest_summary), ncol(contest_summary)))
cat("\n[capstone] convergence matrix (cell_value by axis):\n")
wide <- reshape(
  matrix_tbl[, c("evidence_class", "band", "axis", "cell_value")],
  idvar = c("evidence_class", "band"), timevar = "axis", direction = "wide")
names(wide) <- sub("^cell_value\\.", "", names(wide))
print(wide[, c("evidence_class", "band", AX)], row.names = FALSE)
cat("\n[capstone] contest summary:\n")
print(contest_summary[, c("contest", "winner", "margin",
                          "n_strong", "n_moderate", "n_suggestive",
                          "n_layers", "n_modalities")], row.names = FALSE)

#!/usr/bin/env Rscript
# Assemble summary.Rmd from the curated per-panel chunk files.
#
# This is a *generator*: it reads each verified chunk body verbatim from the
# (gitignored) staging dir storage/cache/summary_chunks/, wraps it in a fenced
# knitr chunk with panel-specific options and British-English narrative prose,
# and writes a single self-contained summary.Rmd at the project root. The
# emitted summary.Rmd is the committed deliverable; this script + the chunk
# dir are the build inputs (re-run after editing any chunk or the prose here).
#
#   Rscript scripts/build_summary_rmd.R
#
# The chunk bodies assume an assembler that (a) source()s R/helpers.R and
# (b) defines two loaders, rd() for storage/cache/*.rds and rt() for
# storage/results/*.tsv; both are emitted into the setup preamble below.
# Three chunks end in knitr::kable() and are tagged results='asis' so they
# render as HTML tables rather than literal markdown (see CLAUDE.md).

chunk_dir <- "storage/cache/summary_chunks"
out_path  <- "summary.Rmd"

stopifnot(dir.exists(chunk_dir))

# --- panel order: (section header | NULL), id == chunk basename, chunk opts,
#     narrative prose emitted immediately before the fenced chunk -------------
panels <- list(
  list(
    section = "Study design",
    section_intro = paste0(
      "The study crosses tau and amyloid pathology in a 2x2 factorial and reads ",
      "each animal out across four molecular modalities. The divergence question ",
      "-- does amyloid remodel microglia *differently* depending on the tau ",
      "background -- is operationalised as a single interaction contrast."),
    id = "design_cohort_table",
    opts = "results='asis'",
    prose = paste0(
      "Four hippocampal genotypes span a 2x2 tau (wild-type humanised MAPT vs ",
      "P301S) by amyloid (no-NLGF vs NLGF knock-in) design, profiled at 24 ",
      "months. snRNAseq contributes 23,465 analysed microglia nuclei (drawn from ",
      "a 286,287-nucleus all-cell-type atlas), GeoMx 91 microglia-enriched ",
      "spatial AOIs, and bulk proteome and phosphoproteome 16 samples each (four ",
      "biological replicates per genotype). Microglial nuclei roughly double ",
      "under NLGF -- the first hint of amyloid-driven microglial expansion.")
  ),
  list(
    section = NULL,
    id = "design_factorial_contrasts",
    opts = "fig.width=8.5, fig.height=6",
    prose = paste0(
      "Five named contrasts sit on the design. The two amyloid (NLGF) effects and ",
      "two tau (P301S) effects are the grid edges; the **interaction** -- ",
      "(NLGF_P301S - P301S) - (NLGF_MAPTKI - MAPTKI) -- is the ",
      "differential-of-differentials that isolates tau-dependent amyloid ",
      "remodelling. It is the project's core readout: the single strongest ",
      "mechanism (Gsk3b kinase activation, padj 0.00152) concentrates here, even ",
      "though the whole-microglia transcriptional interaction is essentially null.")
  ),

  list(
    section = "The microglial atlas",
    section_intro = paste0(
      "Microglia resolve into four transcriptional substates, and amyloid ",
      "reshapes their proportions far more than tau does."),
    id = "mg_substate_composition",
    opts = "fig.width=7, fig.height=4.8",
    prose = paste0(
      "Across 26,104 microglial nuclei (pre-QC-pruning; substates by ",
      "AddModuleScore argmax), the disease-associated (DAM) fraction roughly ",
      "doubles under amyloid on both tau backgrounds -- MAPTKI 19.8% to 40.1%, ",
      "P301S 19.3% to 41.3% -- at the expense of the homeostatic compartment ",
      "(~65% to ~44%). Tau alone barely shifts composition. This additive amyloid ",
      "expansion is the compositional backdrop against which every downstream ",
      "divergence readout is interpreted.")
  ),

  list(
    section = "The divergence is null at the whole-microglia level",
    section_intro = paste0(
      "The natural expectation -- that tau-dependent amyloid remodelling shows up ",
      "as differentially expressed genes -- fails. At the whole-microglia level ",
      "the interaction contrast is null, and the null is robust to method, ",
      "modality and analytic shrinkage."),
    id = "sn_nebula_interaction_volcano",
    opts = "fig.width=7.5, fig.height=5.2",
    prose = paste0(
      "Single-cell NEBULA (negative-binomial GLMM, 23,465 cells x 11,411 genes) ",
      "finds no gene clearing FDR < 0.05 at the interaction; the leading hit ",
      "Plac9 sits at FDR 0.063. The cloud is centred on the origin -- there is no ",
      "whole-microglia interaction programme to speak of.")
  ),
  list(
    section = NULL,
    id = "sn_triangulation_concordance",
    opts = "fig.width=10, fig.height=5.2",
    prose = paste0(
      "The null is not a NEBULA artefact. Six DE methods spanning three shrinkage ",
      "families agree: pairwise Spearman concordance of interaction log2FC (and of ",
      "-log10 P) is uniformly high *among methods*, yet they all rank the same ",
      "near-zero effects -- method choice does not rescue a signal.")
  ),
  list(
    section = NULL,
    id = "gx_sn_interaction_concordance",
    opts = "fig.width=10, fig.height=8",
    prose = paste0(
      "Spatial transcriptomics corroborates. GeoMx reproduces the snRNAseq ",
      "amyloid main effects (Spearman rho 0.117 in MAPTKI, 0.208 in P301S) but the ",
      "interaction contrast is null cross-modally (rho -0.048, flat-to-negative ",
      "slope). The spatial modality is not globally noisy -- the main effects ",
      "validate -- so the interaction's absence is real, not a sensitivity floor.")
  ),
  list(
    section = NULL,
    id = "xmod_concordance_heatmap",
    opts = "fig.width=11, fig.height=5",
    prose = paste0(
      "Pooling all four modalities, the only appreciable cross-modality agreement ",
      "sits on the proteome-phospho axis (proteome vs phospho rho 0.21-0.35; ",
      "phospho vs corrected-phospho 0.47-0.53). Transcript-versus-protein ",
      "concordance is near-null throughout, including at the interaction. This is ",
      "the empirical case for abandoning a whole-microglia transcriptional readout ",
      "and pivoting to mechanism-resolved layers.")
  ),

  list(
    section = "Where the divergence actually lives",
    section_intro = paste0(
      "The divergence has not vanished -- it has moved off the whole-microglia ",
      "gene axis. Four orthogonal readouts localise it: to specific substates, to ",
      "coordinated pathways, to an unsupervised co-expression module, and to ",
      "neighbourhood abundance."),
    id = "sn_per_state_ifn_asymmetry",
    opts = "fig.width=7.5, fig.height=6.5",
    prose = paste0(
      "Within IFN microglia the amyloid response is sharply tau-context-dependent: ",
      "446 genes respond to NLGF only on the P301S background versus 56 only on ",
      "MAPTKI (82 shared), so the P301S-specific cloud lifts off the identity ",
      "line. This per-substate asymmetry is the mathematical substrate of the ",
      "otherwise-null whole-microglia interaction, and Pros1 -- the surviving ",
      "TAM-kinase route from the synaptic axis -- is among the top P301S-specific ",
      "responders.")
  ),
  list(
    section = NULL,
    id = "sn_per_state_de_summary",
    opts = "results='asis'",
    prose = paste0(
      "Quantifying the localisation: the interaction is null in *every* substate ",
      "(0 genes at FDR < 0.05), confirming the whole-microglia null is not a ",
      "substate-averaging artefact. The divergence reads off the asymmetry between ",
      "the two NLGF arms -- IFN is P301S-leaning (528 vs 138 significant genes) ",
      "whereas proliferative inverts to MAPTKI-leaning (139 vs 500) -- and genes ",
      "reappear under the amyloid-conditioned tau effect (tau_in_nlgf: IFN 101, ",
      "DAM 25).")
  ),
  list(
    section = NULL,
    id = "sn_per_state_pathway_breadth",
    opts = "fig.width=8.5, fig.height=8.5",
    prose = paste0(
      "Where individual genes are silent, coordinated pathways are not. GSEA on ",
      "the per-state NEBULA t-statistics recovers 43 GO BP terms at padj < 0.1 in ",
      "at least one substate, ~90% negatively signed: oxidative phosphorylation, ",
      "cellular respiration, cytoplasmic translation and ribosome biogenesis are ",
      "suppressed at the interaction, predominantly in DAM and homeostatic ",
      "microglia (e.g. DAM ribonucleoprotein-complex biogenesis NES -1.73, padj ",
      "7.1e-4).")
  ),
  list(
    section = NULL,
    id = "hdwgcna_module_logfc_heatmap",
    opts = "fig.width=8, fig.height=4",
    prose = paste0(
      "An unsupervised co-expression module carries the same echo. Of four ",
      "hdWGCNA modules, only MG-M3 has a negative interaction response (logFC ",
      "-2.91), collapsing under the amyloid-conditioned tau effect (logFC -3.80, ",
      "FDR 6.7e-4). The DAM module MG-M2 instead activates under every NLGF ",
      "contrast (logFC 5.0-6.2), the expected amyloid-driven DAM programme.")
  ),
  list(
    section = NULL,
    id = "milo_da_umap_nlgf",
    opts = "fig.width=10, fig.height=5.2",
    prose = paste0(
      "At the abundance level, miloR neighbourhood DA shows NLGF expanding the DAM ",
      "island and depleting the homeostatic island on both backgrounds. But ",
      "detectability is tau-asymmetric: 132/1,023 neighbourhoods clear SpatialFDR ",
      "< 0.1 on the P301S background versus *none* on MAPTKI, despite near-",
      "identical raw DAM-fraction shifts (0.19 to 0.37 vs 0.18 to 0.40).")
  ),

  list(
    section = "The mechanism layer",
    section_intro = paste0(
      "Resolved to kinases, intercellular signalling and pathway families, the ",
      "divergence becomes concrete -- and one result dominates."),
    id = "prot_phospho_correction_scatter",
    opts = "fig.width=9, fig.height=5",
    prose = paste0(
      "A methodological prerequisite first. Subtracting the matched parent-protein ",
      "abundance from each phosphosite isolates phosphorylation stoichiometry ",
      "before limma is refit. Raw and corrected site-level log2FC correlate only ",
      "moderately (Spearman rho 0.71/0.71 across the two NLGF contrasts, ~12.8k ",
      "sites each), so a substantial fraction of apparent phospho-regulation is ",
      "parent-protein-driven and is removed here. All kinase inference downstream ",
      "runs on the corrected estimates.")
  ),
  list(
    section = NULL,
    id = "prot_gsk3b_kinase_interaction",
    opts = "fig.width=8, fig.height=4",
    prose = paste0(
      "The headline. decoupleR kinase-activity inference over the OmniPath ",
      "kinase-substrate network puts **Gsk3b at +4.32 (ULM) at the interaction, ",
      "FDR 0.00152** (consensus FDR 1.8e-4) -- the strongest single-mechanism ",
      "signal in the project. The pattern is a tau-dependent sign-flip: Gsk3b is ",
      "non-significantly negative when amyloid lands on wild-type tau ",
      "(nlgf_in_maptki -1.66, ns) but strongly positive when it lands on a ",
      "tau-mutant background (nlgf_in_p301s +3.79, FDR 0.0077). Gsk3b activation ",
      "is therefore specifically a tau x amyloid synergy. (Bulk hippocampal ",
      "phospho, not microglia-sorted: read as kinase activity anywhere in the ",
      "tissue at this contrast.)")
  ),
  list(
    section = NULL,
    id = "ccc_interaction_lr_microglia",
    opts = "fig.width=9, fig.height=7.5",
    prose = paste0(
      "The intercellular layer rewires where the transcriptome does not. Across ",
      "the top-25 microglia-receiver ligand-receptor pairs prioritised for the ",
      "interaction (MultiNicheNet), DAM is the receiver in 18/25 and sender in ",
      "10/25, and the ligands collapse to four macrophage-trophic/clearance ",
      "molecules (Apoe, C1qb, Csf1, Icam1) -- including the Apoe->Trem2 clearance ",
      "route. None of these LR genes overlaps an hdWGCNA hub, so the signalling ",
      "signal sits off the module axis.")
  ),
  list(
    section = NULL,
    id = "ccc_axis2_lr_adjudication",
    opts = "fig.width=9.5, fig.height=6.5",
    prose = paste0(
      "Adjudicating the synaptic-suppression axis across three CCC tools (6,459 ",
      "universe-filtered LR cells), TREM2-mediated clearance (Apoe_Trem2 rank 9), ",
      "APP-fragment uptake (App_Cd74 rank 12) and glial adhesion rank shallow, ",
      "whereas classical complement ranks deep (first C1qb_Lrp1 at rank 775; no C3 ",
      "or Mertk-Pros1 in the top 800). The clearance signature, not textbook ",
      "complement pruning, is what the data surface. The residual TAM-kinase ",
      "Pros1_Mertk route survives separately at IFN-microglia nodes.")
  ),
  list(
    section = NULL,
    id = "pathway_survey_leaderboard",
    opts = "fig.width=10, fig.height=8.5",
    prose = paste0(
      "Pooling six gene-set collections under one fixed leader rule, the board ",
      "refuses to name a single winner: three axes coexist -- tau x amyloid ",
      "metabolic/translational non-additivity (mixed sign, led by MG-M3 at 23.4 ",
      "with the OXPHOS/ETC family and DAM_up at 16.4), amyloid-driven activation ",
      "(positive; MG-M2, adaptive immunity), and NLGF synaptic suppression ",
      "(negative; presynaptic GO CC terms). The mixed-sign interaction family ",
      "reaches breadth 4/4 substates.")
  ),

  list(
    section = "Integrated synthesis",
    section_intro = paste0(
      "The claims ledger adjudicates the competing models. On all three axes the ",
      "data favour a divergent model over a simple additive one, and the ",
      "cross-axis themes rank consistently."),
    id = "per_state_nfkb_attenuation",
    opts = "fig.width=9.5, fig.height=8",
    prose = paste0(
      "The axis-1 mechanism, resolved. Tau on the amyloid background suppresses ",
      "amyloid-driven NF-kB output at the interaction in all four substates: the ",
      "union NF-kB-target GSEA reverses from positive at both NLGF arms (NES ~+1.0 ",
      "to +1.2) to negative at the interaction (DAM NES -1.25, padj 0.027), ",
      "mirrored by the CollecTRI NF-kB-complex TF score (DAM -6.5, p 7e-11), DAM ",
      "leading. The substate-composition ANOVA (tau:nlgf p 0.33, ns) rules out a ",
      "Simpson's-paradox abundance artefact -- the attenuation is transcriptional ",
      "and global.")
  ),
  list(
    section = NULL,
    id = "biomodel_contest_verdicts",
    opts = "results='asis'",
    prose = paste0(
      "Head-to-head, each axis favours its divergent model: axis 1 -> tau ",
      "attenuates amyloid-driven NF-kB (margin 18), axis 2 -> TREM2/APP-fragment ",
      "clearance over classical complement (margin 12), axis 3 -> a distinct ",
      "tau x amyloid synergy mechanism (margin 55, the largest), anchored by ",
      "Gsk3b. All three resolve at the first arithmetic level with no tie-break ",
      "invoked.")
  ),
  list(
    section = NULL,
    id = "crossaxis_theme_ranking",
    opts = "fig.width=8.5, fig.height=5",
    prose = paste0(
      "Ranking the cross-axis themes by supporting ledger rows: T-Inflammation 31 ",
      "(6 Strong) > T-Synergy 30 (9 Strong) > T-Compartment-suppression 17 > ",
      "T-Tau-attenuates 12 > Hyp-0 (Cdk5) 7. Additive DAM-amplification ",
      "(T-Inflammation) now edges past the Gsk3b/Myc synergy theme on raw ",
      "supporting-row count, though T-Synergy retains the most Strong-grade ",
      "support (9 versus 6); the Cdk5 integrator remains the most evenly ",
      "cross-axis-distributed entity. This is the project's standing synthesis.")
  ),

  list(
    section = "Cross-species validation in human AD",
    section_intro = paste0(
      "Every result above is from four mouse genotypes in a clean 2x2 design. ",
      "The closing question is translational: does the same interaction ",
      "structure exist in human AD cortex, where amyloid and tau co-progress ",
      "and cannot be crossed experimentally? Using SEA-AD (Gabitto 2024) ",
      "middle-temporal-gyrus and prefrontal microglia (no data-use agreement ",
      "required), the mouse substates and their interaction-localised ",
      "mechanisms are tested for conservation and for direction-of-effect ",
      "concordance under collinearity-aware modelling. Because human amyloid ",
      "(Thal phase) and tau (Braak stage) are collinear (Spearman ~0.65), ",
      "*direction* -- not significance -- is the readout, and a null is read as ",
      "under-identification, never as evidence against the mouse finding. The ",
      "full treatment is section 19 of the analysis report."),
    id = "human_substate_conservation",
    opts = "results='asis'",
    prose = paste0(
      "Per-state human claims are only admissible if the mouse substates ",
      "resolve in human microglia at all. All four resolve as human SEA-AD ",
      "populations, and the directional signal is correct: disease-emergent ",
      "supertype occupancy rises monotonically along homeostatic -> DAM -> IFN ",
      "-> proliferative, with homeostatic the only state *depleted* in disease ",
      "supertypes. The conservation gate PASSES, clearing the per-state ",
      "interaction tests below.")
  ),
  list(
    section = NULL,
    id = "human_mechanism_concordance",
    opts = "fig.width=8, fig.height=3.2",
    prose = paste0(
      "The headline translational result. The three pre-registered ",
      "interaction mechanisms reproduce their mouse-predicted amyloid:tau sign ",
      "in human: **NF-kB-target** and **MG-M3-module** attenuation are negative ",
      "(tau dampens the amyloid response) and the **Gsk3b-target** response is ",
      "positive (tau-enhanced), in 5 of 5 strata, region-concordant and robust ",
      "to single-cell aggregation and an orthogonal AUCell scorer. No human ",
      "interaction survives FDR -- the expected consequence of amyloid-tau ",
      "collinearity and observational cross-sectional power -- so the reading ",
      "is directional corroboration (*consistent with* the mouse mechanism), ",
      "never causal proof. The human-native Gerrits positive control behaves ",
      "correctly (section 19).")
  ),

  list(
    section = "The integrated convergence model",
    section_intro = paste0(
      "Eleven evidence arcs (D--O) layered onto the four base modalities, plus ",
      "the human translational test, converge on one integrated picture, ",
      "rendered below as an evidence-layer by biological-axis matrix. Three ",
      "readings sit side by side. First, amyloid-driven DAM activation is the ",
      "strongest and most multi-layer-convergent axis -- positive and at least ",
      "moderately supported in nearly every layer, from expression through ",
      "kinase activity, composition and progression. Second, the tau x amyloid ",
      "interaction, the project's core readout, does not live on the static gene ",
      "axis: matched-power gene-level testing collapses it to null in every cell ",
      "type, which the specificity arc (N) reframes as not microglia-specific, ",
      "and it is likewise null in neighbourhood composition (spatial arc L); yet ",
      "it re-emerges as a positive progression-dynamics effect (trajectory arc ",
      "M) and is mechanism-signed at the kinase and transcription-factor layers ",
      "(Gsk3b up, Myc down). Third, the later arcs (K--O) are deliberately ",
      "margin-neutral -- each contributes only one or two Suggestive ledger rows ",
      "-- so they corroborate the decided contests orthogonally without moving ",
      "them. The human cross-species layer is held separate: directionally ",
      "consistent with the mouse mechanisms, never counted into the mouse ",
      "contests."),
    id = "capstone_convergence",
    opts = "results='asis'",
    prose = paste0(
      "Panel A is the landscape view -- each arc's own signed verdict per axis, ",
      "with nulls and reframings first-class: grey null cells, the gold reframe ",
      "cell where the specificity arc recast the interaction as not ",
      "microglia-specific, and the purple mixed cells where modalities or ",
      "ligand-receptor pairs disagree in sign. Panel B is the contest view, the ",
      "three sealed section-17 margins (55 / 18 / 12). The two views are ",
      "complementary: the contest arithmetic discounts the margin-neutral arcs ",
      "that the landscape view shows are nonetheless corroborating, which is ",
      "precisely why both are shown. The table beneath re-surfaces the sealed ",
      "margins read-only. With this synthesis the project reaches its planned ",
      "close -- further arcs would add orthogonal corroboration at diminishing ",
      "return rather than move any verdict.")
  )
)

# --- concise human-readable panel titles (used as numbered ## subheadings,
#     so the TOC reads as a narrative rather than chunk ids) -----------------
panel_titles <- c(
  design_cohort_table             = "Cohort: four genotypes across four modalities",
  design_factorial_contrasts      = "The 2x2 design and the interaction contrast",
  mg_substate_composition         = "DAM doubles under amyloid",
  sn_nebula_interaction_volcano   = "snRNAseq interaction volcano (NEBULA)",
  sn_triangulation_concordance    = "Six-method concordance of the interaction",
  gx_sn_interaction_concordance   = "Spatial validates the main effects, not the interaction",
  xmod_concordance_heatmap        = "Cross-modality concordance matrix",
  sn_per_state_ifn_asymmetry      = "IFN microglia carry a tau-context-dependent response",
  sn_per_state_de_summary         = "Per-substate differential-expression summary",
  sn_per_state_pathway_breadth    = "Pathway suppression localises to DAM and homeostatic",
  hdwgcna_module_logfc_heatmap    = "hdWGCNA module MG-M3 echoes the interaction",
  milo_da_umap_nlgf               = "Neighbourhood abundance shift, tau-asymmetric",
  prot_phospho_correction_scatter = "Phosphosite stoichiometric correction",
  prot_gsk3b_kinase_interaction   = "Gsk3b activation at the interaction (headline)",
  ccc_interaction_lr_microglia    = "DAM-dominated ligand-receptor rewiring",
  ccc_axis2_lr_adjudication       = "Axis-2 clearance: TREM2/APP over complement",
  pathway_survey_leaderboard      = "Unified pathway leader board: three axes",
  per_state_nfkb_attenuation      = "Tau attenuates NF-kB across all substates",
  biomodel_contest_verdicts       = "Three-axis biological-model contest verdicts",
  crossaxis_theme_ranking         = "Cross-axis theme support ranking",
  human_substate_conservation     = "Mouse substates resolve in human microglia",
  human_mechanism_concordance     = "Interaction mechanisms reproduce in human (headline)",
  capstone_convergence            = "Convergence matrix across all evidence layers"
)
stopifnot(setdiff(vapply(panels, `[[`, "", "id"), names(panel_titles)) |> length() == 0)

# --- static document scaffolding ------------------------------------------
yaml_header <- c(
  "---",
  'title: "Microglia divergence in tau and amyloid AD models"',
  'subtitle: "Highest-impact figures and results -- internal findings summary"',
  'author: "Claude (autonomous agent)"',
  "date: \"`r format(Sys.Date(), '%d %B %Y')`\"",
  "output:",
  "  html_document:",
  "    toc: true",
  "    toc_float: true",
  "    toc_depth: 2",
  "    number_sections: true",
  "    code_folding: hide",
  "    df_print: paged",
  "    fig_width: 8",
  "    fig_height: 6",
  "params:",
  "  publication_mode: false",
  '  cache_dir: "storage/cache"',
  '  results_dir: "storage/results"',
  "---",
  ""
)

setup_chunk <- c(
  "```{r setup, include = FALSE}",
  "knitr::opts_chunk$set(",
  "  echo = TRUE, message = FALSE, warning = FALSE,",
  '  fig.align = "center", dev = "png", dpi = 110, cache = FALSE',
  ")",
  'project_root <- normalizePath(".")',
  "knitr::opts_knit$set(root.dir = project_root)",
  "options(width = 110)",
  "",
  "# Project-wide libraries + domain helpers (genotype_levels/_colours,",
  "# contrast_definitions, proteomics_sample_meta, plot_kinase_activity_heatmap).",
  'source("R/helpers.R")',
  "",
  "# Summary-specific loaders the curated chunks assume:",
  "#   rd() -> pre-extracted minimal plotting caches in storage/cache/",
  "#   rt() -> exported result tables in storage/results/",
  "rd <- function(f) readRDS(file.path(params$cache_dir, f))",
  "rt <- function(f) read.delim(file.path(params$results_dir, f),",
  "                             check.names = FALSE, stringsAsFactors = FALSE)",
  "```",
  ""
)

overview <- c(
  "# Overview {-}",
  "",
  paste0(
    "This document curates the highest-impact figures and results from an ",
    "integrated analysis of microglia across a 2x2 tau x amyloid mouse model, ",
    "spanning single-nucleus and spatial transcriptomics, bulk proteomics and ",
    "phosphoproteomics. The central finding is a *negative result with a ",
    "constructive resolution*: tau-dependent amyloid remodelling of microglia ",
    "does not appear as differentially expressed genes at the whole-microglia ",
    "level in any modality, yet it is robustly present once the readout is ",
    "resolved to substates, co-expression modules, neighbourhood abundance, ",
    "kinase activity and intercellular signalling. The strongest single ",
    "mechanism is **Gsk3b kinase activation at the tau x amyloid interaction** ",
    "(padj 0.00152), a tau-dependent synergy."),
  "",
  paste0(
    "This is an *internal synthesis*: phosphoproteomics involving the MAPT-KI ",
    "background is shown, effect sizes and caveats are stated plainly, and ",
    "findings are framed mechanistically rather than for publication. Every ",
    "panel re-reads its data from disk (pre-extracted caches and exported result ",
    "tables) and re-renders without recomputing the underlying analysis, so this ",
    "report knits in roughly a minute and never touches the heavy Seurat or ",
    "MultiNicheNet objects."),
  ""
)

caveats <- c(
  "# Scope and caveats {-}",
  "",
  paste0(
    "- Bulk phospho/proteome are whole-hippocampus, **not** microglia-sorted, so ",
    "kinase-activity readouts report tissue-level activity at each contrast."),
  paste0(
    "- snRNAseq composition figures are pre-QC-pruning (26,104 nuclei) whereas DE ",
    "uses the pruned subset (~23,465); counts differ accordingly."),
  paste0(
    "- NEBULA produces extreme natural-log estimates for near-zero-count genes; ",
    "per-substate scatters clip *display* only (no points or statistics dropped)."),
  paste0(
    "- Cross-tool CCC ranks are universe-restricted and tool-dependent; the GeoMx ",
    "interaction volcano (not shown here) is dominated by predicted-gene spatial ",
    "confounds, which is why the cross-modal concordance panel carries the null."),
  paste0(
    "- All caches and result tables are re-read from disk; this document ",
    "re-renders curated figures without recomputing the analysis."),
  ""
)

# --- assemble -------------------------------------------------------------
emit <- character(0)
emit <- c(emit, yaml_header, setup_chunk, overview)

for (p in panels) {
  if (!is.null(p$section)) {
    emit <- c(emit, paste0("# ", p$section), "")
    if (!is.null(p$section_intro)) emit <- c(emit, p$section_intro, "")
  }
  chunk_file <- file.path(chunk_dir, paste0(p$id, ".R"))
  if (!file.exists(chunk_file)) stop("missing chunk file: ", chunk_file)
  body <- readLines(chunk_file, warn = FALSE)
  # Drop a trailing blank line so the fenced block closes tightly.
  while (length(body) > 0 && !nzchar(body[length(body)])) {
    body <- body[-length(body)]
  }
  emit <- c(emit, paste0("## ", panel_titles[[p$id]]), "")
  emit <- c(emit, p$prose, "")
  emit <- c(emit, sprintf("```{r %s, %s}", p$id, p$opts))
  emit <- c(emit, body)
  emit <- c(emit, "```", "")
}

emit <- c(emit, caveats)

writeLines(emit, out_path)
cat("wrote", out_path, "with", length(panels), "panels\n")

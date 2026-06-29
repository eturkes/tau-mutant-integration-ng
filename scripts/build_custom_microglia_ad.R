#!/usr/bin/env Rscript
# scripts/build_custom_microglia_ad.R
#
# Phase B Step B3 of the pathway-overhaul plan. Builds the curated
# tau / amyloid-specific microglia gene-set collection (4 sets, mouse
# symbols) and runs fgsea against the five DE modalities x five
# contrasts.
#
# Set inventory (alphabetical):
#   AD1         (n ~47)  -- amyloid-associated microglia from Gerrits 2022,
#                           top-50 up-regulated genes vs homeostatic
#                           subclusters (p_val_adj < 0.05, top by
#                           avg_logFC), human -> mouse via nichenetr.
#   AD2         (n ~42)  -- tau-associated microglia from Gerrits 2022,
#                           top-50 up-regulated genes vs homeostatic
#                           subclusters (p_val_adj < 0.05, top by
#                           avg_logFC), human -> mouse via nichenetr.
#                           Effect sizes much weaker than AD1; signature
#                           selected by top-N rank rather than fixed
#                           logFC threshold to keep set sizes comparable.
#   LDAM        (n =27)  -- Lipid-droplet-accumulating microglia from
#                           Marschallinger 2020 Nat Neurosci, T2-1
#                           RNA-Seq Aging sheet filtered for log2FC>0.5
#                           AND padj<0.05.
#   WAM         (n =39)  -- White-matter-associated microglia from
#                           Safaiyan 2021 Neuron, mmc4 Venn Analysis
#                           "Common Genes (39)" -- pre-curated union of
#                           (WM-21m up vs GM-21m) and (WM-21m up vs WM-4m).
#
# All gene lists are EMBEDDED in this script so the build has no
# external file dependency at run-time except `nichenetr` (used to
# translate the human Gerrits AD1/AD2 lists to mouse). Source
# supplements were parsed once by the assistant during the B3 session;
# this script captures the resulting vectors verbatim so the build is
# fully reproducible offline.
#
# See `storage/cache/custom_microglia_ad_provenance.txt` (written by
# this script) for source-paper / supplement-table / filter details
# per set.
#
# Idempotent: skips any cache that already exists unless `--overwrite`
# is passed. Re-running with no args after a successful build is a
# no-op.
#
# Outputs (under storage/cache/):
#   custom_microglia_ad.rds              (named list of mouse-symbol vectors)
#   fgsea_custom_ad_results.rds          (modality -> contrast -> fgseaResult)
#   custom_microglia_ad_provenance.txt
#
# Outputs (under storage/results/):
#   fgsea_custom_ad_per_contrast.tsv     (joined wide table)

suppressPackageStartupMessages({
  library(dplyr)
  library(nichenetr)
})

args      <- commandArgs(trailingOnly = TRUE)
overwrite <- "--overwrite" %in% args

setwd("/home/rstudio/tau-mutant-integration-ng")
source("R/helpers.R")

cache_dir   <- "storage/cache"
results_dir <- "storage/results"
dir.create(cache_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

gene_set_path <- file.path(cache_dir, "custom_microglia_ad.rds")
fgsea_path    <- file.path(cache_dir, "fgsea_custom_ad_results.rds")
tsv_path      <- file.path(results_dir, "fgsea_custom_ad_per_contrast.tsv")
prov_path     <- file.path(cache_dir, "custom_microglia_ad_provenance.txt")

# =========================================================================
# Embedded gene lists.
#
# WAM:       Safaiyan 2021 mmc4 "Venn Analysis" sheet, column "Common
#            Genes (39)". Pre-curated intersection of (WM-21m up vs GM-
#            21m) and (WM-21m up vs WM-4m) DEG sets -- the paper's own
#            white-matter-microglia signature, mouse symbols verbatim.
#
# LDAM:      Marschallinger 2020 NIHMS1544918 supplement, sheet
#            "T2 1- RNA Seq Aging". Filter: log2FoldChange > 0.5 AND
#            padj < 0.05 in the BODIPY-Hi vs BODIPY-Lo aged-microglia
#            DEG, yielding 27 up-regulated mouse symbols. Includes the
#            canonical lipid-droplet coat-protein Plin3 and the
#            Pluvinage 2019 LDAM marker Cd22, plus aged-microglia
#            transcripts.
#
# AD1_human: Gerrits 2022 Acta Neuropathol Table S7, cols 1-6 ("AD1 vs
#            homeostasis"). Filter: p_val_adj < 0.05, sorted by
#            avg_logFC descending, top 50 retained. This 50-gene human
#            list is converted to mouse here via nichenetr (~47 mouse
#            symbols after mapping NA-drop and dedup). AD1 is the
#            amyloid-associated, DAM-like microglia state in Gerrits
#            2022 (Gpnmb, Spp1, Lpl, Msr1, Olr1, Cd83 -- classical
#            DAM-overlap markers).
#
# AD2_human: Gerrits 2022 Table S7, cols 29-34 ("AD2 vs homeostasis").
#            Filter: p_val_adj < 0.05, sorted by avg_logFC descending,
#            top 50 retained. AD2 effect sizes are much weaker than AD1
#            (max logFC ~0.51 vs AD1's 2.0) so a top-N rank-based cut
#            keeps the set comparable in size to AD1 (~42 mouse symbols
#            after nichenetr mapping). AD2 is the tau-associated,
#            non-DAM-like state in Gerrits 2022 (Foxp1, Foxp2, Slit2,
#            Grid2, Npas3, Mrc1, Igf1r -- axon-guidance / neuron-
#            interaction / homeostatic-disturbance markers).
# =========================================================================

WAM <- c(
  "4930556M19Rik", "Actr3b", "AF529169", "Ahnak2", "Apoc1",
  "Apoe", "Atp6v0d2", "Axl", "Cacna1a", "Ch25h",
  "Clec7a", "Cox6a2", "Cst7", "Cxcl13", "Ddo",
  "Egln3", "Etl4", "Gas2l3", "Gm26714", "Gm44620",
  "Gm5424", "Gpnmb", "Igf1", "Il1rn", "Itgax",
  "Lox", "Lpl", "Mmp12", "Olfr110", "Olr1",
  "Pdcd1", "Pianp", "Postn", "Ptchd1", "Rab7b",
  "Spp1", "St8sia6", "Stra6l", "Tnfsf8"
)

LDAM <- c(
  "2900026A02Rik", "Ace", "Angptl7", "Arhgap19", "Bfsp1",
  "Cd22", "Cryz", "Dusp27", "Eif5a2", "Fnip2",
  "Itpr3", "Kcnj13", "Kcnj9", "Kl", "Lrp11",
  "Lypd6", "Nbl1", "Npl", "Plin3", "Psat1",
  "Rapsn", "Slc13a4", "Slc38a3", "Slc7a10", "Sulf1",
  "Tmem8b", "Zfp507"
)

AD1_human <- c(
  "PTPRG", "TPRG1", "MYO1E", "GLDN", "DPYD",
  "CD83", "CPM", "GAS7", "XYLT1", "SPP1",
  "DIRC3", "ARID5B", "EYA2", "ADAMTS17", "PPARG",
  "STARD13", "ASAP1", "MITF", "COLEC12", "GPNMB",
  "SAMD4A", "AC066613.2", "RGCC", "ADARB1", "CADM1",
  "WIPF3", "OLR1", "IQGAP2", "MSR1", "TRHDE",
  "AL158071.4", "FAM135B", "ZNF804A", "PADI2", "RASGRP3",
  "ATG7", "ALCAM", "SH3PXD2A", "ELL2", "CDK14",
  "LPL", "RASGEF1B", "NAMPT", "CYTH3", "SCIN",
  "FMNL2", "PRKCE", "NPL", "MAFB", "FLT1"
)

AD2_human <- c(
  "GRID2", "DPP10", "ADGRB3", "NPAS3", "PID1",
  "LRP1B", "NRG3", "C22orf34", "SSBP2", "C12orf42",
  "GPM6A", "SLIT2", "NRXN1", "PRDM11", "LSAMP",
  "FOXP1", "MRC1", "LRRC4C", "IGF1R", "PCDH9",
  "UNC5C", "C5orf17", "FAM155A", "TPCN1", "AFF3",
  "CDH23", "FOXP2", "THRB", "TMEM71", "RORA",
  "TVP23A", "LINC01141", "AC092691.1", "NR3C2", "PICALM",
  "CACNB4", "ZNF710", "MBNL1", "ZBTB20", "MEF2C-AS1",
  "CXADR", "GPC5", "NCAM1", "TNRC6B", "XIST",
  "C4orf19", "KCNIP4", "LMCD1-AS1", "OSBPL3", "CPED1"
)

# =========================================================================
# Build the named gene-set list, translating AD1_human / AD2_human via
# nichenetr.
# =========================================================================

build_gene_sets <- function() {
  message("Converting Gerrits AD1 human symbols to mouse via nichenetr")
  ad1_mouse <- nichenetr::convert_human_to_mouse_symbols(AD1_human)
  ad1_mouse <- unique(ad1_mouse[!is.na(ad1_mouse) & nzchar(ad1_mouse)])
  message(sprintf("  AD1 mapping: human n=%d  ->  mouse n=%d",
                  length(AD1_human), length(ad1_mouse)))

  message("Converting Gerrits AD2 human symbols to mouse via nichenetr")
  ad2_mouse <- nichenetr::convert_human_to_mouse_symbols(AD2_human)
  ad2_mouse <- unique(ad2_mouse[!is.na(ad2_mouse) & nzchar(ad2_mouse)])
  message(sprintf("  AD2 mapping: human n=%d  ->  mouse n=%d",
                  length(AD2_human), length(ad2_mouse)))

  list(
    AD1  = sort(ad1_mouse),
    AD2  = sort(ad2_mouse),
    LDAM = sort(unique(LDAM)),
    WAM  = sort(unique(WAM))
  )
}

gene_sets <- cache_or_run(gene_set_path, build_gene_sets(),
                          overwrite = overwrite)

message("\n=== gene-set sizes ===")
for (nm in names(gene_sets)) {
  message(sprintf("  %-6s n=%4d", nm, length(gene_sets[[nm]])))
}

# =========================================================================
# Provenance file.
# =========================================================================

if (overwrite || !file.exists(prov_path)) {
  writeLines(c(
    "# custom_microglia_ad.rds  --  provenance",
    sprintf("# Built: %s", format(Sys.time())),
    sprintf("# Script: scripts/build_custom_microglia_ad.R"),
    "",
    "## Collection summary",
    "# Four gene sets (mouse symbols) representing tau / amyloid-",
    "# associated and lipid / white-matter-associated microglia states.",
    "# Final counts after de-duplication and human-to-mouse mapping:",
    sprintf("#   AD1   n=%d  (Gerrits human n=%d  ->  mouse via nichenetr)",
            length(gene_sets$AD1), length(AD1_human)),
    sprintf("#   AD2   n=%d  (Gerrits human n=%d  ->  mouse via nichenetr)",
            length(gene_sets$AD2), length(AD2_human)),
    sprintf("#   LDAM  n=%d  (Marschallinger T2-1, logFC>0.5 & padj<0.05)",
            length(gene_sets$LDAM)),
    sprintf("#   WAM   n=%d  (Safaiyan mmc4 Venn 'Common Genes' verbatim)",
            length(gene_sets$WAM)),
    "",
    "## Decision-gate deviation from plan B3 default",
    "# The pathway-overhaul plan B3 default proposed Mancuso 2024 (ex vivo",
    "# xenografts) as the paired tau-iMG / amyloid-iMG source. WebSearch",
    "# verified only one canonical 2024 chimeric-microglia paper",
    "# (Nat Neurosci PMC11089003) and it is primarily amyloid-focused with",
    "# no paired tau-iMG signature published. The user selected the",
    "# Recommended alternative Gerrits 2022 (Acta Neuropathol,",
    "# DOI 10.1007/s00401-021-02263-w; PMC8043951) which DOES publish",
    "# paired amyloid-associated (AD1) and tau-associated (AD2) microglia",
    "# subcluster DEG lists in a single Table S7. This switch tightens the",
    "# tau / amyloid distinction the project needs while keeping a single",
    "# coherent within-paper comparison.",
    "",
    "## WAM  --  Safaiyan et al. 2021 Neuron",
    "# Paper:        DOI 10.1016/j.neuron.2021.01.027",
    "# Supplement:   Cell Press CDN, 1-s2.0-S0896627321000738-mmc4.xlsx",
    "# Sheet:        'Venn Analysis'",
    "# Column:       'Common Genes (39)' -- pre-curated intersection of",
    "#               (WM-21m up vs GM-21m) and (WM-21m up vs WM-4m) DEG",
    "#               sets, the paper's own white-matter-associated",
    "#               microglia signature.",
    "# Genes:        Mouse symbols verbatim, 39 entries. Includes the",
    "#               canonical lipid / phagocytic / lysosomal axis",
    "#               (Gpnmb, Spp1, Cst7, Clec7a, Itgax, Lpl, Mmp12, Igf1,",
    "#               Apoe, Apoc1, Axl, Pdcd1, Lox, Atp6v0d2, Ch25h) plus",
    "#               WM-specific markers (Cxcl13, Postn, Stra6l, Cacna1a,",
    "#               Il1rn, Tnfsf8, Pianp).",
    "# Note:         Safaiyan 2021 is NOT a PMC-deposited paper; the",
    "#               supplement was fetched directly from the Cell Press",
    "#               CDN URL (PII S0896627321000738).",
    "",
    "## LDAM  --  Marschallinger et al. 2020 Nature Neuroscience",
    "# Paper:        DOI 10.1038/s41593-019-0566-1",
    "# Supplement:   PMC OA tarball PMC7595134,",
    "#               NIHMS1544918-supplement-1.xlsx",
    "# Sheet:        'T2 1- RNA Seq Aging' (12117 rows, aged vs young",
    "#               microglia DESeq2 output).",
    "# Filter:       log2FoldChange > 0.5 AND padj < 0.05 (up in aged /",
    "#               LDAM-enriched cells). Yields 27 mouse symbols.",
    "# Genes:        Plin3 (lipid-droplet coat protein, canonical LDAM",
    "#               marker) + Cd22 (Pluvinage 2019 LDAM marker) + 25",
    "#               aged-microglia DEGs (Slc7a10, Slc38a3, Slc13a4, Ace,",
    "#               Kl, Bfsp1, Sulf1, Kcnj13, Psat1, ...). The top hits",
    "#               include unusual transcripts (transporters / receptors)",
    "#               that reflect the paper's data-driven aged-microglia",
    "#               DEG signature; the canonical lipid-droplet biology",
    "#               markers (Plin2, Apoe) are present in the wider 73-",
    "#               gene candidate panel (T1 LD-related genes) but did",
    "#               not pass the DEG threshold in the aged-microglia",
    "#               comparison.",
    "# CAVEAT:       Marschallinger 2020 did NOT publish a fixed 'LDAM",
    "#               gene set'; the LDAM phenotype is defined by FACS",
    "#               (BODIPY-Hi sorting), so any tabular signature is a",
    "#               filter of the published DEG table. User confirmed",
    "#               T2-1 filtered as the canonical choice during the B3",
    "#               curation gate (alternatives considered: 73-gene T1",
    "#               LD-panel; union of T2-1 + T1).",
    "",
    "## AD1 + AD2  --  Gerrits et al. 2022 Acta Neuropathologica",
    "# Paper:        DOI 10.1007/s00401-021-02263-w",
    "# Supplement:   PMC OA tarball PMC8043951,",
    "#               401_2021_2263_MOESM2_ESM.xlsx (Supplementary",
    "#               Information 2, all data tables).",
    "# Sheet:        'Table S7' (Figure 4c/d -- DEG between AD1 / AD2",
    "#               microglia and homeostatic subclusters 0, 1 and 5).",
    "# Layout:       8 horizontal column-blocks of 6 cols each",
    "#               (gene, p_val, avg_logFC, pct.1, pct.2, p_val_adj),",
    "#               header at row 4. AD1 is cols 1-6, AD2 is cols 29-",
    "#               34; intermediate blocks (subclusters 7, 9, 10, 2,",
    "#               3, 6) are not used.",
    "# Filter:       p_val_adj < 0.05, sort by avg_logFC descending, top",
    "#               50 retained per block. For AD1 this corresponds to",
    "#               avg_logFC range 0.52-2.00 (within the strongly-up",
    "#               regime). For AD2 the same top-50 corresponds to",
    "#               avg_logFC range 0.085-0.514 (effect sizes are much",
    "#               weaker because AD2 is biologically a heterogeneous",
    "#               state distinguished from homeostasis by many small",
    "#               shifts rather than a few strong markers).",
    "# Conversion:   Human symbols -> mouse via",
    "#               nichenetr::convert_human_to_mouse_symbols.",
    "#               Unmapped symbols (mostly LINC / AC / AL lncRNAs)",
    "#               are silently dropped.",
    "# AD1 vs AD2 biology:",
    "#               AD1 (amyloid-associated) overlaps strongly with the",
    "#               mouse DAM phenotype: Gpnmb, Spp1, Lpl, Msr1, Olr1,",
    "#               Cd83 -- lipid handling / phagocytic markers. AD2",
    "#               (tau-associated) instead carries axon-guidance and",
    "#               neuron-interaction transcripts (Foxp1, Foxp2,",
    "#               Slit2, Grid2, Npas3, Pcdh9, Unc5c, Ncam1) plus",
    "#               microglia receptors (Mrc1, Igf1r) consistent with a",
    "#               homeostasis-disturbed, non-DAM-like activation",
    "#               state.",
    "",
    "## Anti-anchoring guardrails respected",
    "# - No Hallmark gene sets are introduced (the collection has zero",
    "#   overlap with the deleted Hallmark caches).",
    "# - All gene-set names are explicit and descriptive (no rank-",
    "#   position anchoring). Sort order in the .rds is alphabetical so",
    "#   downstream heatmaps cannot anchor on any one set's position.",
    "# - AD1 and AD2 set sizes are made comparable by top-N rank cut",
    "#   rather than a unified fixed logFC threshold, to avoid an",
    "#   artificial sample-size asymmetry from the published effect-size",
    "#   differential."
  ), prov_path)
  message(sprintf("Wrote provenance: %s", prov_path))
} else {
  message(sprintf("[skip] provenance already present: %s", prov_path))
}

# =========================================================================
# fgsea against the five DE modalities.
#
# min_size lowered to 5 to match B2: the curated AD signatures are
# tight (LDAM has 27 mouse symbols, which after intersecting with the
# proteomics / phospho gene universe may drop to ~5-15). The plan's
# verification requires that every state appears in every modality x
# contrast row of the wide TSV.
# =========================================================================

shim_nebula <- function(nebula_obj) list(fit = list(top = nebula_obj$top))

message("\nLoading DE caches")
de_snrnaseq_nebula <- readRDS(file.path(cache_dir, "de_snrnaseq_nebula.rds"))
de_geomx           <- readRDS(file.path(cache_dir, "de_geomx.rds"))
de_proteomics      <- readRDS(file.path(cache_dir, "de_proteomics.rds"))
de_phospho         <- readRDS(file.path(cache_dir, "de_phospho.rds"))
de_phospho_corr    <- readRDS(file.path(cache_dir, "de_phospho_corrected.rds"))

results <- cache_or_run(fgsea_path, {
  list(
    snrnaseq     = run_fgsea_per_dataset(shim_nebula(de_snrnaseq_nebula),
                                         gene_sets, min_size = 5),
    geomx        = run_fgsea_per_dataset(de_geomx,      gene_sets, min_size = 5),
    proteomics   = run_fgsea_per_dataset(de_proteomics, gene_sets, min_size = 5),
    phospho      = run_fgsea_per_dataset(de_phospho,    gene_sets, min_size = 5),
    phospho_corr = run_fgsea_per_dataset(de_phospho_corr, gene_sets,
                                         min_size = 5)
  )
}, overwrite = overwrite)

joint <- join_fgsea_results(results)
write_tsv_safe(joint, tsv_path)
message(sprintf("\nWrote: %s  (%d pathway-contrast rows)",
                tsv_path, nrow(joint)))

# Final sanity-check: `join_fgsea_results` returns one row per
# (pathway, contrast) with per-modality columns, so the expected row
# count is `n_sets * n_contrasts`. NA cells are allowed for any (state,
# modality, contrast) tuple where the state fell below min_size after
# intersecting with the modality's gene universe.
expected_rows <- length(gene_sets) * length(unique(joint$contrast))
n_rows        <- nrow(joint)
message(sprintf("Sanity: expected %d wide-form rows (%d sets x %d contrasts); got %d.",
                expected_rows, length(gene_sets),
                length(unique(joint$contrast)), n_rows))
stopifnot(n_rows == expected_rows)

message("done")

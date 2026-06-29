#!/usr/bin/env Rscript
# scripts/build_custom_module_sources.R
#
# Phase B Step B4 of the pathway-overhaul plan. Builds the curated
# microglia transcriptional-state "module sources" gene-set collection
# (12 sets, mouse symbols) from Sun, Victor, Mathys et al. 2023 Cell
# and runs fgsea against the five DE modalities x five contrasts.
#
# Set inventory (alphabetical by state ID, MG0..MG12 with MG9 skipped):
#   MG0_homeostatic            -- canonical resting-state microglia
#                                 (P2ry12, Cx3cr1, Mrc1, Sall1-equivalents).
#   MG1_neuronal_surveillance  -- neuron-interacting microglia (axon-guidance
#                                 / synaptic-pruning transcripts).
#   MG2_inflammatory_I         -- Cpeb4 / Tmem163 / Il4ra-marked inflammatory
#                                 state (paper main-text markers).
#   MG3_ribosome_biogenesis    -- Rpl / Rps / Tpt1 / Eef1a1 protein-synthesis
#                                 axis.
#   MG4_lipid_processing       -- Pparg / Gpnmb / Lipa / Itgax / Olr1
#                                 cholesterol-efflux / lipid-laden state.
#   MG5_phagocytic             -- Apoe / Msr1 / Cd163 / F13a1 / Clec5a
#                                 phagocytic effector state.
#   MG6_stress_response        -- HSP / Dnaj / Fos / Jun unfolded-protein-
#                                 response state.
#   MG7_glycolytic             -- Hif1a / Hk2 / Pgk1 / Pfkfb3 / Slc2a3 /
#                                 Vegfa hypoxia-glycolysis axis.
#   MG8_inflammatory_II        -- Lrrk2 / Spon1 / Foxp1 inflammatory state
#                                 (paper main-text markers).
#   MG10_inflammatory_III      -- Il1b / Ccl3 / Relb / Nfkbi* canonical
#                                 NF-kB inflammatory state.
#   MG11_antiviral             -- Mx1 / Oas3 / Ifit2/3 / Stat1/2 / Rsad2
#                                 type-I-IFN-stimulated state.
#   MG12_cycling               -- Brca1/2 / Brip1 / Chek1 / Cenp* / Dnmt1
#                                 proliferating microglia.
#
# All gene lists are EMBEDDED in this script so the build has no
# external file dependency at run-time except `nichenetr` (used to
# translate the human Sun/Victor 2023 cluster markers to mouse). Source
# supplements were parsed once by the assistant during the B4 session
# from the Broad CDN-hosted marker text file
# (ROSMAP.Microglia.6regions.seurat.harmony.selected.clusterDEGs.txt);
# this script captures the resulting vectors verbatim so the build is
# fully reproducible offline.
#
# See `storage/cache/custom_module_sources_provenance.txt` (written by
# this script) for source-paper / supplement-table / filter details
# per set, including the deviation note explaining why this collection
# is sourced from Sun/Victor 2023 instead of the plan's original
# Mathys 2019 + Olah 2020 defaults.
#
# Idempotent: skips any cache that already exists unless `--overwrite`
# is passed. Re-running with no args after a successful build is a
# no-op.
#
# Outputs (under storage/cache/):
#   custom_module_sources.rds            (named list of mouse-symbol vectors)
#   fgsea_custom_modules_results.rds     (modality -> contrast -> fgseaResult)
#   custom_module_sources_provenance.txt
#
# Outputs (under storage/results/):
#   fgsea_custom_modules_per_contrast.tsv  (joined wide table)

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

gene_set_path <- file.path(cache_dir, "custom_module_sources.rds")
fgsea_path    <- file.path(cache_dir, "fgsea_custom_modules_results.rds")
tsv_path      <- file.path(results_dir, "fgsea_custom_modules_per_contrast.tsv")
prov_path     <- file.path(cache_dir, "custom_module_sources_provenance.txt")

# =========================================================================
# Embedded gene lists.
#
# Source: ROSMAP.Microglia.6regions.seurat.harmony.selected.clusterDEGs.txt
#         from https://personal.broadinstitute.org/cboix/sun_victor_et_al_data/
#         (Sun, Victor, Mathys et al. 2023 Cell). Seurat FindAllMarkers output
#         with columns p_val / avg_log2FC / pct.1 / pct.2 / p_val_adj /
#         cluster / gene, 2228 markers across 12 clusters (0-12 with 9 skipped).
#
# Filter: p_val_adj < 0.05, sort by avg_log2FC descending, top 50 per
#         cluster. After human -> mouse mapping via
#         nichenetr::convert_human_to_mouse_symbols and NA-drop / dedup,
#         final mouse-symbol counts per state are 33-50.
#
# State labels follow the biological annotations from Sun/Victor 2023
# main text (Figure 1 and Table S1, pages 1-4); MG9 is skipped because
# the published annotation merged that putative cluster into MG0
# during the final 12-state taxonomy.
# =========================================================================

MG0_homeostatic <- c(
  "0610040J01Rik", "1500009L16Rik", "Abcc4", "Adamtsl2", "Ankrd44",
  "Appl2", "Arhgap22", "Blnk", "Cacna1a", "Cped1",
  "Csgalnact1", "Cx3cr1", "Deptor", "Dip2a", "Eml1",
  "Foxp2", "Frmd4a", "Grid2", "Hs3st4", "Ifngr1",
  "Il6st", "Kbtbd12", "Kcnip1", "Khdrbs3", "Klf12",
  "Lpar1", "Mrc1", "Nav2", "Nav3", "Nlrp1b",
  "Oxr1", "P2ry12", "P3h2", "Pgm5", "Pip5k1b",
  "Prdm11", "Ptchd4", "Rasgef1c", "Rttn", "Sdk1",
  "Sesn3", "Slc15a2", "St6galnac3", "Syndig1", "Tanc1",
  "Tbc1d4", "Tenm4", "Tiam1", "Tln2", "Tmem156"
)

MG1_neuronal_surveillance <- c(
  "Abr", "Aoah", "Bod1l", "C3", "Cacna1a",
  "Ccnd3", "Csf2ra", "Frmd4a", "Gak", "Golgb1",
  "Ino80d", "Lingo1", "Maf", "Mlxipl", "Mphosph8",
  "Nktr", "Nlrp1b", "Pnisr", "Pram1", "Prpf38b",
  "Prrc2c", "Rasgef1b", "Rasgef1c", "Rbm25", "Ryr1",
  "Slc26a3", "Son", "Srrm2", "Syndig1", "Tnrc18",
  "Wnt2b", "Zfp846", "Zmat1"
)

MG2_inflammatory_I <- c(
  "Acer3", "Ank2", "Bach1", "Bcl6", "Bin1",
  "Cd14", "Cpeb4", "Dennd3", "Erc2", "Fam129a",
  "Fcgbp", "Golga4", "Gramd1a", "Hamp", "Hrh2",
  "Hsp90b1", "Hspa5", "Ikzf2", "Il1rap", "Il4ra",
  "Jak3", "Klhl2", "Lap3", "Limk2", "March1",
  "Nck2", "Pim1", "Ppa1", "Prkca", "Ptpn1",
  "Ptpn2", "Ptpre", "Pygl", "Rnf149", "Runx1",
  "Scfd2", "Sh3rf3", "Slc11a1", "Slc2a5", "Spata6",
  "Spp1", "Tbc1d14", "Tbc1d8", "Tgfbr1", "Tmem163",
  "Vsig4", "Xbp1"
)

MG3_ribosome_biogenesis <- c(
  "Aif1", "Apoo", "Eef1a1", "Fau", "Fth1",
  "Ftl1", "Gm9844", "Naca", "Ooep", "Plekha7",
  "Rpl10", "Rpl11", "Rpl13", "Rpl13a", "Rpl14",
  "Rpl19", "Rpl23", "Rpl27a", "Rpl28", "Rpl29",
  "Rpl3", "Rpl32l", "Rpl35", "Rpl37a", "Rpl41",
  "Rpl6", "Rpl8", "Rplp1", "Rplp2", "Rps11",
  "Rps15", "Rps18", "Rps19", "Rps2", "Rps20",
  "Rps23", "Rps24", "Rps25", "Rps27", "Rps27a",
  "Rps3", "Rps6", "Rps8", "Tpt1", "Uba52-ps",
  "Ybx1"
)

MG4_lipid_processing <- c(
  "Arid5b", "Asah1", "Atg7", "Cadm1", "Cd83",
  "Cpm", "Dennd4c", "Dpyd", "Dscam", "Farp1",
  "Fmn1", "Fnip2", "Frmd4b", "Gas7", "Gldn",
  "Gpnmb", "Iqgap2", "Itgax", "Kcnj5", "Kcnma1",
  "Lgmn", "Lipa", "Mafb", "Mgll", "Mitf",
  "Mtss1", "Myo1e", "Nhsl1", "Npl", "Olr1",
  "Pla2g7", "Pparg", "Prkce", "Ptprg", "Rasgrp3",
  "Rttn", "Samd4", "Scin", "Sh3pxd2a", "Socs6",
  "Spire1", "Spp1", "Stard13", "Susd1", "Tanc2",
  "Tnfaip2", "Tprg", "Wipf3", "Xylt1", "Zfp804a"
)

MG5_phagocytic <- c(
  "Acsl1", "Adgrg6", "Apoe", "Arhgap18", "Cadm1",
  "Cd14", "Cd163", "Clec5a", "Cpm", "Dpyd",
  "Dram1", "Eda", "Epb41l3", "Esr1", "Eya2",
  "F13a1", "Fam110b", "Fcgr3", "Fcho2", "Fmn1",
  "Gpr155", "Ifi44l", "Il15", "Iqgap1", "Iqgap2",
  "Itsn1", "Khdrbs2", "Lrrk2", "Mafb", "Mctp1",
  "Ms4a4a", "Ms4a6d", "Msr1", "Nckap5", "Pde4b",
  "Plxnc1", "Rhobtb3", "Slc2a9", "St7", "Stard13",
  "Tgfbi", "Thrb", "Tprg", "Trps1", "Utrn",
  "Zbtb16", "Zfp804a"
)

MG6_stress_response <- c(
  "Bag3", "Bcas2", "Cacybp", "Cd83", "Chordc1",
  "Clk1", "Ddit4", "Dnaja1", "Dnaja4", "Dnajb1",
  "Dnajb6", "Dusp1", "Elovl5", "Fcgr3", "Fkbp4",
  "Fos", "Glul", "Gna13", "Gpr183", "Hif1a",
  "Hist1h2ao", "Hmox1", "Hsp90aa1", "Hsp90ab1", "Hspa1a",
  "Hspa8", "Hspb1", "Hspd1", "Hsph1", "Ier5",
  "Jun", "Msn", "Nampt", "P4ha1", "Ppp1r15a",
  "Ptges3", "Rgs1", "Rgs16", "Ripk2", "Serpinh1",
  "Slc2a3", "Slc7a5", "Snap23", "Srgn", "Stip1",
  "Ubc", "Uspl1"
)

MG7_glycolytic <- c(
  "Acsl1", "Atp1b3", "Bcat1", "Bnip3l", "Cd83",
  "Crem", "Cxcr4", "Ddit4", "Dusp1", "Ell2",
  "Epb41l3", "Fam110b", "Fam13a", "Fos", "Fosl2",
  "Gbe1", "Glul", "Gna13", "Gykl1", "Hif1a",
  "Hk2", "Malt1", "Mb21d2", "Nampt", "P4ha1",
  "Padi2", "Papolg", "Per1", "Pfkfb3", "Pgk1",
  "Plaur", "Rab20", "Ranbp2", "Rapgef1", "Rgs1",
  "Rnf144b", "Sat1", "Sipa1l1", "Slc11a1", "Slc2a3",
  "Slc7a5", "Spp1", "Srgn", "Sytl3", "Tfrc",
  "Tnfrsf1b", "Usp36", "Vegfa"
)

MG8_inflammatory_II <- c(
  "Adgrb3", "Alox8", "Ccnd3", "Cebpd", "Cpvl",
  "Csgalnact1", "Cttnbp2", "Cyp27a1", "Dusp1", "Dusp10",
  "Fam110b", "Fam49a", "Fkbp5", "Fos", "Gab1",
  "Gpr155", "Grid2", "Hif1a", "Hist1h2ao", "Hspa1a",
  "Il15", "Insr", "Irs2", "Itpr1", "Klhl6",
  "Lrrk2", "Man2a1", "Mgat1", "Msr1", "Nhsl1",
  "Oxr1", "Padi2", "Peli2", "Pik3ip1", "Plcl1",
  "Plxnc1", "Prkag2", "Rgs1", "Spon1", "Srgn",
  "Sytl3", "Tnpo1", "Tsc22d3", "Ttc7", "Txnip",
  "Zbtb16", "Zfp36l2", "Zfp804a"
)

MG10_inflammatory_III <- c(
  "Acsl1", "Arhgap31", "B4galt1", "Bcl6", "Btg2",
  "Ccl3", "Cd83", "Cdkn1a", "Clic4", "Ebi3",
  "Ets2", "Fos", "Gch1", "Hivep1", "Hivep2",
  "Ier2", "Ier3", "Il1b", "Irak2", "Klf6",
  "Lcp2", "Limk2", "Lyn", "Mb21d2", "Nampt",
  "Ncf1", "Nfkb1", "Nfkbia", "Nfkbid", "Nfkbiz",
  "Olr1", "Padi2", "Pdgfb", "Peak1", "Plek",
  "Pstpip2", "Relb", "Rgl1", "Rgs16", "Rnf144b",
  "Sat1", "Slc2a3", "Spp1", "Sqstm1", "Srgn",
  "Tnfaip3", "Traf3", "Tymp", "Zfp36"
)

MG11_antiviral <- c(
  "Adar", "Axl", "B2m", "Cul1", "Ddx58",
  "Ddx60", "Eif2ak2", "Ephb2", "Epsti1", "Fmnl2",
  "Herc6", "Ifi44", "Ifi44l", "Ifit1bl2", "Ifit2",
  "Ifit3", "Lap3", "Mx1", "Oas1g", "Oas3",
  "Parp12", "Parp14", "Parp9", "Plscr1", "Pnpt1",
  "Ppm1k", "Rabgap1l", "Rnf213", "Rsad2", "Samd4",
  "Samd9l", "Samhd1", "Siglec1", "Sp110", "Spats2l",
  "Stat1", "Stat2", "Tnfsf13b", "Trim14", "Trim25",
  "Trim56", "Tymp", "Unc93b1", "Xaf1", "Zc3hav1",
  "Zcchc2"
)

MG12_cycling <- c(
  "2610028H24Rik", "Atad2", "Atad5", "Bard1", "Brca1",
  "Brca2", "Brip1", "Ccdc18", "Cenpk", "Cenpp",
  "Cep128", "Cep152", "Chek1", "Cit", "Clspn",
  "Diaph3", "Dnmt1", "Dtl", "Ezh2", "Fanca",
  "Fancc", "Fanci", "Gins1", "Hells", "Kif15",
  "Knl1", "Kntc1", "Lig1", "Lmnb1", "Melk",
  "Mms22l", "Nasp", "Ncapd3", "Ncapg2", "Nsd2",
  "Nup210", "Pola1", "Pole2", "Polq", "Prim2",
  "Rad18", "Rbl1", "Rfc3", "Smc4", "Tacc3",
  "Tcf19", "Uhrf1", "Wdhd1", "Zfp367", "Zgrf1"
)

# =========================================================================
# Build the named gene-set list. All embedded vectors are already mouse
# symbols (the human-to-mouse mapping was performed once during marker
# extraction; see provenance for the audit trail). The build function
# here only dedups and sorts so downstream heatmaps cannot anchor on
# any one set's position.
# =========================================================================

build_gene_sets <- function() {
  list(
    MG0_homeostatic           = sort(unique(MG0_homeostatic)),
    MG1_neuronal_surveillance = sort(unique(MG1_neuronal_surveillance)),
    MG2_inflammatory_I        = sort(unique(MG2_inflammatory_I)),
    MG3_ribosome_biogenesis   = sort(unique(MG3_ribosome_biogenesis)),
    MG4_lipid_processing      = sort(unique(MG4_lipid_processing)),
    MG5_phagocytic            = sort(unique(MG5_phagocytic)),
    MG6_stress_response       = sort(unique(MG6_stress_response)),
    MG7_glycolytic            = sort(unique(MG7_glycolytic)),
    MG8_inflammatory_II       = sort(unique(MG8_inflammatory_II)),
    MG10_inflammatory_III     = sort(unique(MG10_inflammatory_III)),
    MG11_antiviral            = sort(unique(MG11_antiviral)),
    MG12_cycling              = sort(unique(MG12_cycling))
  )
}

gene_sets <- cache_or_run(gene_set_path, build_gene_sets(),
                          overwrite = overwrite)

message("\n=== gene-set sizes ===")
for (nm in names(gene_sets)) {
  message(sprintf("  %-30s n=%4d", nm, length(gene_sets[[nm]])))
}

# =========================================================================
# Provenance file.
# =========================================================================

if (overwrite || !file.exists(prov_path)) {
  writeLines(c(
    "# custom_module_sources.rds  --  provenance",
    sprintf("# Built: %s", format(Sys.time())),
    sprintf("# Script: scripts/build_custom_module_sources.R"),
    "",
    "## Collection summary",
    "# Twelve gene sets (mouse symbols) representing the 12 published",
    "# microglia transcriptional states from Sun, Victor, Mathys et al.",
    "# 2023 Cell (\"Human microglial state dynamics in Alzheimer's disease",
    "# progression\", DOI 10.1016/j.cell.2023.08.037).",
    "# Final counts after de-duplication and human-to-mouse mapping:",
    sprintf("#   MG0_homeostatic            n=%d", length(gene_sets$MG0_homeostatic)),
    sprintf("#   MG1_neuronal_surveillance  n=%d  (mitochondrial MT-* genes",
            length(gene_sets$MG1_neuronal_surveillance)),
    "#                                       drop out of nichenetr mapping)",
    sprintf("#   MG2_inflammatory_I         n=%d", length(gene_sets$MG2_inflammatory_I)),
    sprintf("#   MG3_ribosome_biogenesis    n=%d", length(gene_sets$MG3_ribosome_biogenesis)),
    sprintf("#   MG4_lipid_processing       n=%d", length(gene_sets$MG4_lipid_processing)),
    sprintf("#   MG5_phagocytic             n=%d", length(gene_sets$MG5_phagocytic)),
    sprintf("#   MG6_stress_response        n=%d", length(gene_sets$MG6_stress_response)),
    sprintf("#   MG7_glycolytic             n=%d", length(gene_sets$MG7_glycolytic)),
    sprintf("#   MG8_inflammatory_II        n=%d", length(gene_sets$MG8_inflammatory_II)),
    sprintf("#   MG10_inflammatory_III      n=%d", length(gene_sets$MG10_inflammatory_III)),
    sprintf("#   MG11_antiviral             n=%d", length(gene_sets$MG11_antiviral)),
    sprintf("#   MG12_cycling               n=%d", length(gene_sets$MG12_cycling)),
    "",
    "## Decision-gate deviation from plan B4 default",
    "# The pathway-overhaul plan B4 default proposed two sources:",
    "#   (1) Mathys 2019 Nature Suppl. Table 4 \"M1-M10\" SOM-territory",
    "#       gene-trait correlation modules (cross-cell-type),",
    "#   (2) Olah 2020 Nat Comm Suppl Data 5 cluster markers (clusters",
    "#       1-9 except 7, which is already used as ARM in B2).",
    "# The user selected at the B4 decision gate the Recommended",
    "# alternative: Sun, Victor, Mathys et al. 2023 Cell only. Three",
    "# reasons motivated the swap:",
    "#   1. Sun/Victor 2023 is the dedicated microglia atlas from the",
    "#      same lab that produced Mathys 2019, with 194,000 nuclei from",
    "#      443 subjects (vs Mathys 2019's 48-subject cohort). Statistical",
    "#      power for state-specific markers is dramatically higher.",
    "#   2. Mathys 2019's M1-M10 SOM-territory modules are CROSS-CELL-",
    "#      TYPE (mixing neurons, glia, microglia, OPCs, endothelial",
    "#      genes) because they cluster gene-trait correlation patterns",
    "#      across all cell types. For a microglia-focused project this",
    "#      dilutes the signal with non-microglia genes that GSEA would",
    "#      drop anyway.",
    "#   3. The Olah 2020 ARM cluster (cluster 7) is already in B2; the",
    "#      remaining Olah clusters are smaller and less well-characterised",
    "#      than the corresponding Sun/Victor 2023 states. The Sun/Victor",
    "#      2023 12-state taxonomy provides a cleaner, more contemporary",
    "#      replacement that subsumes most of Olah's biology.",
    "",
    "## MG0-MG12 (skip MG9)  --  Sun, Victor, Mathys et al. 2023 Cell",
    "# Paper:        DOI 10.1016/j.cell.2023.08.037",
    "#               (Cell 186(20):4386-4403.e29, PubMed 37774678,",
    "#                PMC10644954). Not PMC-OA.",
    "# Source file:  ROSMAP.Microglia.6regions.seurat.harmony.selected.clusterDEGs.txt",
    "# URL:          https://personal.broadinstitute.org/cboix/sun_victor_et_al_data/",
    "#               (Broad Institute CDN public-share; downloaded once",
    "#                during B4 session for marker extraction).",
    "# Format:       Seurat FindAllMarkers output; tab-separated; 2228",
    "#               markers across 12 clusters (cluster IDs 0-12 with",
    "#               9 skipped because the published 12-state taxonomy",
    "#               merged the putative MG9 cluster into MG0 during",
    "#               final annotation).",
    "# Columns:      p_val, avg_log2FC, pct.1, pct.2, p_val_adj, cluster,",
    "#               gene. (The leading row.names column carries gene",
    "#               symbols with numeric suffixes for duplicates; the",
    "#               canonical symbol is in the 'gene' column.)",
    "# Filter:       p_val_adj < 0.05  AND  rank top 50 by avg_log2FC",
    "#               descending per cluster. The published p_val column",
    "#               for the strongest markers is exactly 0 due to",
    "#               numerical underflow; ranking by avg_log2FC therefore",
    "#               tiebreaks consistently. Sun/Victor 2023 used the",
    "#               same Bonferroni-corrected threshold for their main-",
    "#               text figures.",
    "# Conversion:   Human symbols -> mouse via",
    "#               nichenetr::convert_human_to_mouse_symbols.",
    "#               Unmapped symbols (mostly mitochondrial MT-* in MG1,",
    "#               X-linked RPS4X / TMSB10 / TMSB4X in MG3, and lncRNA",
    "#               C-prefixed symbols throughout) are silently dropped.",
    "#               Final mouse-symbol counts are 33-50 per state.",
    "",
    "## State annotations (Sun/Victor 2023 main-text Figure 1 + Table S1)",
    "# MG0 (Homeostatic):           P2RY12 / CX3CR1 / SALL1-equivalents.",
    "#                              Canonical resting-state microglia.",
    "# MG1 (Neuronal surveillance): Synaptic / axon-guidance transcripts;",
    "#                              biologically corresponds to surveillance-",
    "#                              competent microglia engaging neurons.",
    "# MG2 (Inflammatory I):        CPEB4 / TMEM163 / IL4RA. Paper-",
    "#                              identified as the dominant inflammatory",
    "#                              state in early AD pathology.",
    "# MG3 (Ribosome biogenesis):   RPL / RPS / EEF1A1 / TPT1 protein-",
    "#                              synthesis axis. Common across activation",
    "#                              states but defines a distinct cluster.",
    "# MG4 (Lipid-processing):      PPARG / GPNMB / LIPA / ITGAX / OLR1.",
    "#                              Cholesterol-efflux / lipid-laden state",
    "#                              overlapping mouse DAM and Marschallinger",
    "#                              LDAM biology.",
    "# MG5 (Phagocytic):            APOE / MSR1 / CD163 / F13A1 / CLEC5A.",
    "#                              Phagocytic effector state.",
    "# MG6 (Stress response):       HSP / DNAJ / FOS / JUN unfolded-protein-",
    "#                              response state. Overlaps the Marsh 2022",
    "#                              ex-vivo activation signature used in B2",
    "#                              (HAM_exvivo_QC); cross-check with that",
    "#                              set is informative.",
    "# MG7 (Glycolytic):            HIF1A / HK2 / PGK1 / PFKFB3 / SLC2A3 /",
    "#                              VEGFA hypoxia-glycolysis axis. The",
    "#                              candidate metabolic-rewiring state.",
    "# MG8 (Inflammatory II):       LRRK2 / SPON1 / FOXP1. Distinct from",
    "#                              MG2 in TF program and lipid-handling",
    "#                              transcripts.",
    "# MG10 (Inflammatory III):     IL1B / CCL3 / RELB / NFKBI*. Canonical",
    "#                              NF-kB-driven inflammatory state.",
    "# MG11 (Antiviral):            MX1 / OAS3 / IFIT2/3 / STAT1/2 / RSAD2.",
    "#                              Type-I-IFN-stimulated state. Overlaps",
    "#                              IRM signature used in B2 substantially.",
    "# MG12 (Cycling):              BRCA1/2 / BRIP1 / CHEK1 / CENP* / DNMT1.",
    "#                              Proliferating microglia (cell-cycle S/G2).",
    "",
    "## Anti-anchoring guardrails respected",
    "# - No Hallmark gene sets are introduced (this collection has zero",
    "#   overlap with the deleted Hallmark caches).",
    "# - All gene-set names are explicit and descriptive (state IDs are",
    "#   prefixed MG0..MG12 in publication order and the alphabetical",
    "#   sort in build_gene_sets places homeostatic first, cycling last;",
    "#   downstream heatmaps cannot anchor on any one set's position).",
    "# - The OXPHOS-axis biology lives in MG7 (glycolytic, the metabolic-",
    "#   rewiring state) and partially in MG3 (ribosome biogenesis). This",
    "#   is the natural Sun/Victor 2023 representation of metabolic",
    "#   microglia states and is NOT framed as the dominant signal; the",
    "#   agnostic survey in Phase C tests this against the other 11 states",
    "#   on equal footing.",
    "",
    "## Other module-source candidates considered (not selected)",
    "# - Mathys 2019 Nature M1-M10 SOM-territory modules (Ext. Data Table S4):",
    "#   cross-cell-type, dominated by non-microglia genes for several",
    "#   modules; smaller cohort (48 subjects).",
    "# - Olah 2020 Nat Comm clusters 1-9 except 7 (Suppl Data 5): older",
    "#   single-cohort taxonomy (14 donors); superseded by Sun/Victor 2023",
    "#   for AD-microglia state discovery. ARM (cluster 7) is already",
    "#   captured in B2.",
    "# - Hasselmann/Blurton-Jones 2019 iPSC-MGL transcriptomic states:",
    "#   in vitro reference, not directly comparable to in vivo state",
    "#   markers from Sun/Victor 2023.",
    "# These sources can be added in a future plan revision if the agnostic",
    "# survey (Phase C9 verdict) identifies a state-axis gap."
  ), prov_path)
  message(sprintf("Wrote provenance: %s", prov_path))
} else {
  message(sprintf("[skip] provenance already present: %s", prov_path))
}

# =========================================================================
# fgsea against the five DE modalities.
#
# min_size lowered to 5 to match B2 and B3: the smallest state set is
# MG1_neuronal_surveillance at 33 mouse symbols, which after intersecting
# with the proteomics / phospho gene universe may drop to ~5-15.
# The plan's verification requires that every state appears in every
# (modality, contrast) row of the wide TSV; min_size=5 keeps the schema
# rectangular while still preserving fgsea's statistical floor.
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
# count is `n_sets * n_contrasts` (NES / padj per modality become
# columns, not rows). NA cells are allowed for any (state, modality,
# contrast) tuple where the state fell below min_size after intersecting
# with the modality's gene universe.
expected_rows <- length(gene_sets) * length(unique(joint$contrast))
n_rows        <- nrow(joint)
message(sprintf("Sanity: expected %d wide-form rows (%d sets x %d contrasts); got %d.",
                expected_rows, length(gene_sets),
                length(unique(joint$contrast)), n_rows))
stopifnot(n_rows == expected_rows)

message("done")

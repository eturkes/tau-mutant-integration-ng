#!/usr/bin/env Rscript
# scripts/build_custom_microglia_states.R
#
# Phase B Step B2 of the pathway-overhaul plan. Builds the curated
# microglia-state gene-set collection (7 sets, mouse symbols) and runs
# fgsea against the five DE modalities × five contrasts.
#
# All gene lists are EMBEDDED in this script so the build has no external
# file dependency at run-time except `nichenetr` (used to translate the
# human Olah ARM list to mouse). Source supplements were parsed once by
# the assistant during the B2 session; this script captures the resulting
# vectors verbatim so the build is fully reproducible offline.
#
# See `storage/cache/custom_microglia_states_provenance.txt` (written by
# this script) for source-paper / supplement-table / filter details per
# set.
#
# Idempotent: skips any cache that already exists unless `--overwrite` is
# passed. Re-running with no args after a successful build is a no-op.
#
# Outputs (under storage/cache/):
#   custom_microglia_states.rds          (named list of mouse-symbol vectors)
#   fgsea_custom_states_results.rds      (modality -> contrast -> fgseaResult)
#   custom_microglia_states_provenance.txt
#
# Outputs (under storage/results/):
#   fgsea_custom_states_per_contrast.tsv  (joined wide table)

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

gene_set_path <- file.path(cache_dir, "custom_microglia_states.rds")
fgsea_path    <- file.path(cache_dir, "fgsea_custom_states_results.rds")
tsv_path      <- file.path(results_dir, "fgsea_custom_states_per_contrast.tsv")
prov_path     <- file.path(cache_dir, "custom_microglia_states_provenance.txt")

# =========================================================================
# Embedded gene lists.
#
# DAM_up / DAM_down:  Keren-Shaul 2017 Cell mmc2.xlsx (Suppl. Table S2),
#                     sheet "Diff.expression_mic3tomic1", filter:
#                       up/down == +1 AND log2((mic3+0.01)/(mic1+0.01)) > 1
#                       up/down == -1 AND log2((mic3+0.01)/(mic1+0.01)) < -1
#                     Result: 278 up (DAM markers) + 21 down (homeostatic
#                     markers repressed in DAM).
#
# HAM_exvivo_QC:      Marsh 2022 Nat Neurosci 41593_2022_1022_MOESM6_ESM.xlsx
#                     (Suppl. Table 4), column "Micro/Myeloid Shared Act.
#                     Score" — 25-gene IEG / heat-shock / chemokine module
#                     induced by enzymatic dissociation across species.
#
# HAM_disease:        Friedman 2018 Cell Rep mmc6.xlsx "By Mouse Gene" sheet,
#                     skip first 11 banner rows, filter
#                       Myeloid Activation (Coarse) == "Neurodegeneration-
#                       Related". 30-gene meta-clustered MGnD-equivalent
#                     module across many AD / neurodegeneration myeloid
#                     datasets.
#
# ARM_human:          Olah 2020 Nat Comm 41467_2020_19737_MOESM7_ESM.xls
#                     (Suppl. Data 5), sheet "degenes_pairwise_microglia_
#                     only", filter up_type == 7. 318 human-symbol genes.
#                     This script converts them to mouse symbols via
#                     `nichenetr::convert_human_to_mouse_symbols`.
#
# IRM:                Sala Frigerio 2019 Cell Rep — the n=28 IRM gene list
#                     referenced in the Methods (for `AddModuleScore` on
#                     Figure 5) is NOT published in the supplement (mmc1.pdf
#                     = figures; mmc2.pdf = paper reprint). We compose 28
#                     mouse symbols by seeding with the 5 IRM markers
#                     explicitly named in the paper's main text (Ifit2,
#                     Ifit3, Ifitm3, Irf7, Oasl2) and augmenting with 23
#                     canonical mouse type-I-interferon-stimulated genes
#                     established in the Schoggins 2011 (Nature 472:481-485)
#                     screen plus the Reactome "Interferon alpha/beta
#                     signaling" mouse-MSigDB collection. This is a
#                     literature-faithful approximation of the unpublished
#                     n=28 list; provenance flags this gap explicitly.
#
# homeostatic_core:   Literature consensus per the pathway-overhaul plan,
#                     covering canonical mouse-microglia identity markers
#                     and key receptors (15 genes). Overlaps cleanly with
#                     Marsh "Microglial Identity Score" (Friedman "Microglia"
#                     coarse module) without being literally copied from
#                     either.
# =========================================================================

DAM_up <- c(
  "5430435G22Rik", "AC151602.1", "Acaca", "Actr3b", "Adarb1",
  "Adcy3", "Adssl1", "AI607873", "Aldoa", "Ank", "Anxa3", "Anxa5",
  "Apbb2", "Aplp2", "Apoc1", "Apoc4", "Apoe", "Aqp6", "Arrdc4",
  "Atf3", "Atox1", "Atp1a3", "Atp1b3", "Atp6v0c", "AW112010", "Axl",
  "B2m", "Baiap2l2", "Bambi", "Bcl2a1d", "Bhlhe40", "Bri3bp",
  "C530028O21Rik", "Cadm1", "Capg", "Cck", "Ccl3", "Ccl4", "Ccl6",
  "Cd300lb", "Cd34", "Cd5", "Cd52", "Cd63", "Cd63-ps", "Cd68",
  "Cd72", "Cd74", "Cd83", "Cd9", "Ch25h", "Chst2", "Clec7a",
  "Colec12", "Coro2a", "Cox6a1", "Cox6a2", "Cox8a", "Cpd", "Creb3l2",
  "Creg1", "Crlf2", "Csf1", "Csf2ra", "Cst7", "Cstb", "Ctsb", "Ctsd",
  "Ctse", "Ctsl", "Ctsz", "Cxcl14", "Cxcl16", "Cxcr4", "Cybb",
  "Dkk2", "Dnajb14", "Dpp7", "Egln3", "Ephx1", "Etl4", "Fabp3",
  "Fabp5", "Fam20c", "Fam46c", "Fblim1", "Flt1", "Fth1", "Fxyd5",
  "Fxyd6", "Galk1", "Gas2l3", "Gas7", "Ghr", "Glb1", "Glipr1",
  "Gm10076", "Gm10086", "Gm10116", "Gm10154", "Gm10269", "Gm10275",
  "Gm10443", "Gm11361", "Gm11428", "Gm11953", "Gm13047", "Gm13196",
  "Gm13341", "Gm13456", "Gm13532", "Gm13570", "Gm13841", "Gm13864",
  "Gm14044", "Gm14059", "Gm14328", "Gm14456", "Gm15427", "Gm15500",
  "Gm15590", "Gm15796", "Gm16020", "Gm16238", "Gm16247", "Gm16379",
  "Gm1673", "Gm17682", "Gm2574", "Gm3511", "Gm4604", "Gm4987",
  "Gm5054", "Gm5239", "Gm5559", "Gm5963", "Gm6030", "Gm6166",
  "Gm6286", "Gm6807", "Gm6977", "Gm7336", "Gm7363", "Gm7670",
  "Gm8129", "Gm8730", "Gm9294", "Gm9396", "Gm9843", "Gnas", "Gpnmb",
  "Gpr65", "Gstm5", "Gsto1", "Gusb", "H2-D1", "H2-K1", "Hexa",
  "Hif1a", "Hint1", "Hmgn3", "Hpse", "Ifi27l2b", "Ifi44", "Igf1",
  "Il18bp", "Il3ra", "Il4i1", "Inf2", "Itgax", "Kcnj2", "Kcnma1",
  "Kif1a", "Ldha", "Leprel2", "Lgals3bp", "Lgi2", "Lilrb4", "Lox",
  "Lpl", "Lrpap1", "Lyz1", "Lyz2", "Maff", "Mamdc2", "Mif",
  "mmu-mir-703", "Mrpl54", "Myeov2", "Myo1e", "Myo5a", "Nceh1",
  "Ndufa1", "Ndufa6", "Ndufs6", "Nexn", "Ngfrap1", "Npc2", "Nrp1",
  "Pdcd1", "Pdgfrl", "Pebp1", "Perp", "Pgcp", "Pgk1", "Pkm2",
  "Plaur", "Plbd2", "Pld3", "Plekhh2", "Plin2", "Prdx4", "Prr5l",
  "Psat1", "Ptchd1", "Ramp1", "Rftn1", "Rpl10a", "Rpl10a-ps1",
  "Rpl12", "Rpl17-ps1", "Rpl18-ps1", "Rpl23", "Rpl23a-ps3",
  "Rpl36a-ps1", "Rpl37", "Rpl37-ps1", "Rpl38", "Rpl38-ps2", "Rpl5",
  "Rpl5-ps1", "Rplp2", "Rps12", "Rps14", "Rps16", "Rps16-ps2",
  "Rps18", "Rps19", "Rps19-ps11", "Rps19-ps3", "Rps2", "Rps21",
  "Rps24", "Rps24-ps2", "Rps24-ps3", "Rps25", "Rps28", "Rps5",
  "Rpsa-ps10", "Rpsa-ps4", "Rpsa-ps9", "Scpep1", "Sdf2l1",
  "Serpine2", "Sh3pxd2b", "Slamf9", "Spp1", "St14", "St8sia6",
  "Sulf2", "Syngr1", "Timp2", "Tlr2", "Tmem163", "Tmem205",
  "Tmem90a", "Tnfsf8", "Tox2", "Tpi1", "Tspo", "Ttyh2", "Tyrobp",
  "Use1", "Usp12", "Vat1", "Vps13c", "Wbp5", "Zfp618"
)

DAM_down <- c(
  "1700017B05Rik", "4632428N05Rik", "Ccr5", "Cd164", "Cmtm6",
  "Crybb1", "Cx3cr1", "Glul", "Lpcat2", "Maf", "Malat1", "Marcks",
  "P2ry12", "Rgs2", "Rhob", "Selplg", "Slco2b1", "Srgap2", "Tmem119",
  "Tmem173", "Txnip"
)

HAM_exvivo_QC <- c(
  "Atf3", "Ccl3", "Ccl4", "Dusp1", "Egr1", "Fos", "Gem", "Hist1h1c",
  "Hist1h2bc", "Hist1h4i", "Hist2h2aa1", "Hsp90aa1", "Hspa1a",
  "Hspa1b", "Ier5", "Jun", "Junb", "Jund", "Klf2", "Nfkbid",
  "Nfkbiz", "Rgs1", "Rhob", "Txnip", "Zfp36"
)

HAM_disease <- c(
  "Asb10", "Ccl6", "Cd34", "Cd68", "Cd83", "Ch25h", "Chst11",
  "Clec7a", "Cspg4", "Cst7", "Ctsz", "Fam46c", "Hcar2", "Igf1",
  "Il1r1", "Itgax", "Lag3", "Ly9", "Lyz2", "Mamdc2", "Mfsd12",
  "Pdcd1", "Plau", "Prr5l", "Rai14", "Siglec5", "St8sia6", "Tlr2",
  "Tmem202", "Tyrobp"
)

ARM_human <- c(
  "A2M", "ABCD4", "ABHD3", "ABI2", "AC074117.10", "ACSL3", "ACY3",
  "ADAP2", "ADORA3", "ADPGK", "ADRB2", "AFF4", "AHSA2", "ANKRD10",
  "ANKRD44", "AOAH", "APOC2", "ARAP1", "ARFGAP2", "ARHGAP12",
  "ARHGAP22", "ARHGAP24", "ARHGAP27", "ARHGEF40", "ARID3A", "ARID4B",
  "ARL6IP5", "ARSG", "ASAH1", "ATF6B", "ATP6AP1", "ATP6AP2",
  "ATP6V0E1", "ATXN2", "ATXN3", "AXL", "BCAP31", "BCO2", "BIN2",
  "BLNK", "C10orf54", "C14orf37", "C1orf132", "C1orf56", "C1QA",
  "C1QB", "C1QC", "C2", "C5orf15", "C5orf45", "C6orf62", "CACNB4",
  "CBR4", "CCND1", "CD164", "CD37", "CD47", "CD59", "CD63", "CD68",
  "CD74", "CD9", "CD99", "CDK5RAP3", "CEPT1", "CFD", "CH17-189H20.1",
  "CHD1", "CHD6", "CHMP1B", "CIITA", "CIRBP", "CLEC9A", "CLK4",
  "CMTM7", "CMTR2", "COMT", "CRYZL1", "CSF2RA", "CSGALNACT1", "CTSL",
  "CTSS", "CTTNBP2NL", "CX3CR1", "CXADR", "CXCL16", "CYBA", "CYBB",
  "CYTL1", "DAGLB", "DDOST", "DEGS1", "DGKZ", "DHRS7", "DIP2A",
  "DNASE2", "DNMT3A", "DOCK4", "DOCK8", "DTX2", "EBI3", "ECHDC2",
  "ELF2", "ELMSAN1", "EPM2AIP1", "ERV3-1", "EVL", "F13A1", "FAM105A",
  "FAM149A", "FAM212A", "FBXL15", "FCGR2A", "FCGR3A", "FCGRT",
  "FEZ2", "FGF20", "FGL2", "FMNL1", "FMNL3", "FNIP1", "FOLR2",
  "FRMD4A", "FTX", "GAB3", "GAL3ST4", "GGA1", "GGNBP2", "GIMAP1",
  "GIMAP4", "GIMAP7", "GLIPR1", "GPM6B", "GPR155", "GPR82", "GRID2",
  "GRINA", "GRN", "HELZ", "HEXA", "HIVEP3", "HLA-C", "HLA-DMA",
  "HLA-DMB", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQA2",
  "HLA-DQB1", "HLA-DQB2", "HLA-DRA", "HLA-DRB1", "HLA-E", "HNRNPDL",
  "HOOK2", "HPS4", "HS3ST1", "HSD17B11", "HTRA1", "IFFO1", "IGSF6",
  "IKBKB", "IL10RA", "IL16", "IL18BP", "IRF5", "IRF9", "ITGB2-AS1",
  "ITM2B", "KCNK6", "KCNMB1", "KCNQ1OT1", "KIAA0141", "KIAA0907",
  "KLHDC2", "LAMP1", "LAPTM4A", "LGMN", "LIMD2", "LINC00685",
  "LINC00865", "LINC00996", "LINC01116", "LINC01374", "LIPA",
  "LMBRD1", "LONP2", "LPAR5", "LPAR6", "LPCAT2", "LRMP", "LTC4S",
  "LUZP1", "LYZ", "MAML1", "MAN2C1", "MARCKS", "MEMO1", "MGAT4A",
  "MGEA5", "MILR1", "MS4A4A", "MS4A6A", "MS4A7", "MSLN", "MTF2",
  "NCF1", "NCK1", "NCKAP1L", "NCSTN", "NDFIP1", "NDRG2", "NHLRC3",
  "NINJ1", "NME3", "NPC2", "NT5C", "OGT", "OLFML3", "P2RX4",
  "P2RY12", "P2RY13", "PAG1", "PAN3", "PAPOLG", "PCF11", "PCMTD1",
  "PFKFB3", "PIH1D1", "PILRA", "PJA2", "PLAT", "PLD3", "PLD4", "POR",
  "PPP1R3E", "PRCP", "PRKAG2", "PRKCD", "PRPF40B", "PTAFR", "PTPRC",
  "RASSF4", "RBM23", "REEP4", "RHOT2", "RICTOR", "RIN3", "RNASET2",
  "RNF13", "RNF149", "RPS4Y1", "RSRP1", "RYR1", "SAMSN1", "SARAF",
  "SDCCAG8", "SELPLG", "SEPP1", "SERINC1", "SERPINA1", "SERPINF1",
  "SF3B1", "SIGLEC10", "SLC22A18", "SLC29A3", "SLC30A7", "SLC35F6",
  "SLC43A2", "SLC50A1", "SMARCC2", "SNAP23", "SPG11", "SPSB3",
  "SRGAP2B", "SRRM1", "SRSF7", "SUMF2", "SUSD3", "SYNGR2", "TAX1BP1",
  "TBXAS1", "THUMPD1", "TLR4", "TM6SF1", "TM9SF2", "TMEM119",
  "TMEM150A", "TMEM165", "TMEM173", "TMEM176B", "TMEM57", "TMEM86A",
  "TMX1", "TNFAIP8L2", "TNFRSF13C", "TNFSF10", "TRAF3IP3", "TREM2",
  "TRIM13", "TRIM22", "UNC50", "UVRAG", "VASH1", "VPS35", "WAC",
  "WDFY2", "WIPI1", "YIPF4", "YPEL2", "YWHAH", "ZFP36L2", "ZMYM2",
  "ZNF207", "ZNF518A", "ZNF7", "ZNF791", "ZSWIM6"
)

# IRM — Sala Frigerio paper-text (5) + canonical mouse ISGs (23) = 28.
IRM <- c(
  # 5 explicitly named in Sala Frigerio 2019 Cell Rep main text
  "Ifit2", "Ifit3", "Ifitm3", "Irf7", "Oasl2",
  # 23 canonical mouse type-I-IFN-stimulated genes (Schoggins 2011 +
  # Reactome IFN-alpha/beta signalling, mouse symbols)
  "Isg15", "Usp18", "Oas1a", "Oas2", "Oas3", "Oasl1",
  "Mx1", "Mx2", "Rsad2", "Cmpk2", "Stat1", "Stat2",
  "Ifi27l2a", "Ifi204", "Ifi207", "Iigp1", "Gbp2", "Gbp3",
  "Slfn5", "Bst2", "Phf11", "Trim30a", "Ddx60"
)

homeostatic_core <- c(
  "P2ry12", "P2ry13", "Tmem119", "Cx3cr1", "Sall1",
  "Tgfbr1", "Csf1r", "Hexb", "Mertk", "Olfml3",
  "Selplg", "Siglech", "Fcrls", "Gpr34", "Adgrg1"
)

# =========================================================================
# Build the named gene-set list, translating ARM_human via nichenetr.
# =========================================================================

build_gene_sets <- function() {
  message("Converting Olah ARM human symbols to mouse via nichenetr")
  arm_mouse <- nichenetr::convert_human_to_mouse_symbols(ARM_human)
  arm_mouse <- unique(arm_mouse[!is.na(arm_mouse) & nzchar(arm_mouse)])
  message(sprintf("  ARM mapping: human n=%d  ->  mouse n=%d",
                  length(ARM_human), length(arm_mouse)))
  list(
    DAM_up           = sort(unique(DAM_up)),
    DAM_down         = sort(unique(DAM_down)),
    HAM_exvivo_QC    = sort(unique(HAM_exvivo_QC)),
    HAM_disease      = sort(unique(HAM_disease)),
    ARM              = sort(arm_mouse),
    IRM              = sort(unique(IRM)),
    homeostatic_core = sort(unique(homeostatic_core))
  )
}

gene_sets <- cache_or_run(gene_set_path, build_gene_sets(),
                          overwrite = overwrite)

message("\n=== gene-set sizes ===")
for (nm in names(gene_sets)) {
  message(sprintf("  %-18s n=%4d", nm, length(gene_sets[[nm]])))
}

# =========================================================================
# Provenance file.
# =========================================================================

if (overwrite || !file.exists(prov_path)) {
  writeLines(c(
    "# custom_microglia_states.rds  --  provenance",
    sprintf("# Built: %s", format(Sys.time())),
    sprintf("# Script: scripts/build_custom_microglia_states.R"),
    "",
    "## Collection summary",
    sprintf("# Seven gene sets (mouse symbols) representing core microglia states."),
    sprintf("# Final counts after de-duplication and human-to-mouse mapping:"),
    sprintf("#   DAM_up           n=%d", length(gene_sets$DAM_up)),
    sprintf("#   DAM_down         n=%d", length(gene_sets$DAM_down)),
    sprintf("#   HAM_exvivo_QC    n=%d", length(gene_sets$HAM_exvivo_QC)),
    sprintf("#   HAM_disease      n=%d", length(gene_sets$HAM_disease)),
    sprintf("#   ARM              n=%d  (human n=%d  ->  mouse n=%d via nichenetr)",
            length(gene_sets$ARM), length(ARM_human), length(gene_sets$ARM)),
    sprintf("#   IRM              n=%d  (literature-faithful approximation;",
            length(gene_sets$IRM)),
    sprintf("#                                   see IRM caveat below)"),
    sprintf("#   homeostatic_core n=%d", length(gene_sets$homeostatic_core)),
    "",
    "## DAM_up / DAM_down  --  Keren-Shaul et al. 2017 Cell",
    "# Paper:        DOI 10.1016/j.cell.2017.05.018",
    "# Supplement:   Cell Press CDN, 1-s2.0-S0092867417305780-mmc2.xlsx (Suppl. Table S2)",
    "# Sheet:        Diff.expression_mic3tomic1",
    "# Columns used: gene, mic1 (homeostatic UMI), mic3 (DAM UMI), up/down flag",
    "# Filter (up):  up/down == +1  AND  log2((mic3+0.01)/(mic1+0.01)) > 1",
    "# Filter (dn):  up/down == -1  AND  log2((mic3+0.01)/(mic1+0.01)) < -1",
    "# Note:         The paper's S2 table is already pre-filtered to DEGs",
    "#               (all rows have -log10(p) >= 7); the additional log2",
    "#               fold cutoff narrows DAM_up from 471 to 278 to bring",
    "#               it closer to the plan's expected size and exclude",
    "#               very-low-magnitude entries.",
    "",
    "## HAM_exvivo_QC  --  Marsh et al. 2022 Nat Neurosci",
    "# Paper:        DOI 10.1038/s41593-022-01022-8",
    "# Supplement:   PMC OA tarball PMC11645269,",
    "#               41593_2022_1022_MOESM6_ESM.xlsx (Suppl. Table 4)",
    "# Column:       Micro/Myeloid Shared Act. Score (25 genes verbatim)",
    "# Rationale:    Marsh et al.'s curated 25-gene module captures the",
    "#               IEG / heat-shock / chemokine activation signature",
    "#               induced by enzymatic dissociation at 37C, conserved",
    "#               across microglia and peripheral myeloid cells. The",
    "#               full DEG table (MOESM7) was considered but the curated",
    "#               score module is the authors' preferred QC signature.",
    "",
    "## HAM_disease  --  Friedman et al. 2018 Cell Reports",
    "# Paper:        DOI 10.1016/j.celrep.2017.12.066",
    "# Supplement:   Cell Press CDN, 1-s2.0-S2211124717319034-mmc6.xlsx",
    "#               (Suppl. Table 5)",
    "# Sheet:        By Mouse Gene (skip first 11 banner rows)",
    "# Filter:       Myeloid Activation (Coarse) == \"Neurodegeneration-Related\"",
    "# Note:         Friedman et al. clustered DEGs across many AD /",
    "#               neurodegeneration myeloid datasets and assigned each",
    "#               gene to a coarse activation category. The",
    "#               Neurodegeneration-Related category (n=30) is the",
    "#               cross-study consensus disease/MGnD module, containing",
    "#               textbook markers (Cst7, Itgax, Clec7a, Tyrobp, Cd68,",
    "#               Igf1, Ch25h) plus less-canonical Friedman-specific",
    "#               additions (Asb10, Ly9, Lag3, Pdcd1, Plau).",
    "",
    "## ARM  --  Olah et al. 2020 Nat Comm",
    "# Paper:        DOI 10.1038/s41467-020-19737-2  (note: NOT 19227-5,",
    "#               which was an early plan typo).",
    "# Supplement:   PMC OA tarball PMC7704703,",
    "#               41467_2020_19737_MOESM7_ESM.xls (Suppl. Data 5)",
    "# Sheet:        degenes_pairwise_microglia_only",
    "# Filter:       up_type == 7  (Olah cluster 7, the MHC-II-loaded",
    "#               antigen-presenting microglia state)",
    "# Conversion:   318 human symbols -> mouse symbols via",
    "#               nichenetr::convert_human_to_mouse_symbols.",
    "#               Yield: see ARM size above (typically ~250-280).",
    "# Note:         The pathway-overhaul plan literally specified Olah",
    "#               2020 as the ARM source. The ARM concept is",
    "#               conventionally attributed to Sala Frigerio 2019",
    "#               (mouse), but Sala Frigerio's full ARM gene list is",
    "#               not published either (see IRM caveat). We honour the",
    "#               plan and use Olah cluster 7 as the cross-species",
    "#               (human-derived, ortholog-mapped) ARM signature.",
    "",
    "## IRM  --  Sala Frigerio et al. 2019 Cell Reports (literature-faithful approximation)",
    "# Paper:        DOI 10.1016/j.celrep.2019.03.099  (note: NOT 01.062;",
    "#               that was an early plan typo).",
    "# Supplement:   Cell Press CDN, 1-s2.0-S2211124719304383-mmc1.pdf",
    "#               (Suppl. figures S1-S6) + mmc2.pdf (full paper reprint).",
    "# CAVEAT:       The supplement does NOT publish the n=28 IRM gene",
    "#               list that the paper's STAR Methods references for",
    "#               their Figure-5 AddModuleScore. The 5 IRM markers",
    "#               explicitly named in the main text are:",
    "#                 Ifit2, Ifit3, Ifitm3, Irf7, Oasl2",
    "#               (see mmc2 lines 184 and 217). We seed with these 5",
    "#               and augment with 23 canonical mouse type-I-IFN-",
    "#               stimulated genes from Schoggins et al. 2011 Nature",
    "#               472:481-485 plus the Reactome \"Interferon alpha/beta",
    "#               signaling\" mouse-MSigDB entries:",
    "#                 Isg15, Usp18, Oas1a, Oas2, Oas3, Oasl1, Mx1, Mx2,",
    "#                 Rsad2, Cmpk2, Stat1, Stat2, Ifi27l2a, Ifi204, Ifi207,",
    "#                 Iigp1, Gbp2, Gbp3, Slfn5, Bst2, Phf11, Trim30a, Ddx60",
    "#               This is a literature-faithful approximation of the",
    "#               unpublished Sala Frigerio IRM signature, not a",
    "#               literal reproduction. User confirmed this strategy",
    "#               during the B2 decision gate.",
    "#               Public-source search history (none yielded the list):",
    "#                 - PMC OA NIH-MS PMC6402798: supplement is figures-only PDF",
    "#                 - Elsevier CDN mmc1.pdf / mmc2.pdf: figures + paper reprint",
    "#                 - GEO GSE127893: no supplementary data files",
    "#                 - Roy 2022 Immunity Tables S1-S3: no re-tabulation",
    "#                 - Soton ePrints 430526: only main-PDF re-host",
    "#                 - ALZFORUM, ScienceDirect: 403 / no supplement listing",
    "",
    "## homeostatic_core  --  literature consensus",
    "# Source:       Pathway-overhaul plan section B2 default + key",
    "#               receptors. 15 mouse symbols:",
    "#                 P2ry12, P2ry13, Tmem119, Cx3cr1, Sall1,",
    "#                 Tgfbr1, Csf1r, Hexb, Mertk, Olfml3,",
    "#                 Selplg, Siglech, Fcrls, Gpr34, Adgrg1",
    "# Verification: matches Marsh 2022 \"Microglial Identity Score\"",
    "#               (18-gene list) on 11 of 15 entries, and Friedman 2018",
    "#               \"Microglia\" coarse module (11-gene list) on 9 of 15.",
    "# Note:         Trem2 deliberately excluded (sits at homeostatic-to-",
    "#               DAM transition).",
    "",
    "## Anti-anchoring guardrails respected",
    "# - No Hallmark gene sets are introduced (the collection has zero",
    "#   overlap with the deleted Hallmark caches).",
    "# - All gene-set names are explicit and descriptive (no rank-position",
    "#   anchoring). Sort order in the .rds is alphabetical so downstream",
    "#   heatmaps cannot anchor on any one set's position."
  ), prov_path)
  message(sprintf("Wrote provenance: %s", prov_path))
} else {
  message(sprintf("[skip] provenance already present: %s", prov_path))
}

# =========================================================================
# fgsea against the five DE modalities.
#
# min_size lowered to 5 because the curated state signatures are tight
# (the smallest is homeostatic_core at 15 mouse symbols, which may drop
# to ~7-10 after intersecting with the proteomics / phospho gene
# universe). The plan's verification requires that every state appears in
# every modality x contrast.
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

# Project constants: factor levels, contrast definitions, canonical marker
# lists, contamination markers, and raw-data paths. Loaded by tar_source(); downstream
# R/ files + targets read these symbols. Values carried byte-exact from v1
# (genotype labels and marker sets are load-bearing -- do not let them drift).

# Canonical genotype order -> factor levels everywhere downstream (MAPTKI = reference).
genotype_levels <- c("MAPTKI", "P301S", "NLGF_MAPTKI", "NLGF_P301S")

# Pairwise contrasts (numerator, denominator) driving the divergence narrative. The 5th
# canonical contrast `interaction` = (NLGF_P301S - P301S) - (NLGF_MAPTKI - MAPTKI) is a
# difference-of-differences, not a level pair -> built in R/core/design.R (S3), keyed by name.
contrast_definitions <- list(
  nlgf_in_maptki = c("NLGF_MAPTKI", "MAPTKI"),
  nlgf_in_p301s  = c("NLGF_P301S",  "P301S"),
  tau_alone      = c("P301S",       "MAPTKI"),
  tau_in_nlgf    = c("NLGF_P301S",  "NLGF_MAPTKI")
)

# Pan-microglia identity markers -- state-INDEPENDENT core, expressed across homeostatic AND
# activated/DAM microglia (Csf1r/C1q complement/Ctss/Fcrls/Hexb/Tyrobp). The QC purity signature
# (P1-S2): a cluster scoring LOW here is non-microglial contamination. Deliberately DISTINCT from
# the Homeostatic SUBPOPULATION below (DAM legitimately downregulates P2ry12/Tmem119) -> homeostatic
# markers must NEVER be used to test "is this a microglia at all".
microglia_identity_markers <- c("Csf1r", "C1qa", "C1qb", "C1qc", "Ctss", "Fcrls", "Hexb", "Tyrobp")

# Canonical mouse-microglia subpopulation signatures (P1-S2). UCell rank-based enrichment -> z-scale
# per signature -> argmax. The 4 subpopulation sets {Homeostatic, DAM, IFN, Proliferative} drive the
# subpopulation argmax; MHC_APC is an AUXILIARY antigen-presentation axis (scored + reported, NOT
# argmax'd -- it co-varies with DAM: ARM = DAM + MHC, Sala Frigerio 2019). DAM merges the
# Trem2-independent (s1) and Trem2-dependent (s2) arms: DAM~MGnD~ARM~WAM is ONE Apoe-Trem2
# convergent programme, scored broadly so the ~18% snRNA DAM-gene dropout (Thrupp 2020) cannot
# null it. Symbols -> ensembl via symbol_map at scoring time (assay rownames are ensembl).
canonical_microglia_markers <- list(
  Homeostatic   = c("P2ry12", "P2ry13", "Cx3cr1", "Tmem119", "Hexb", "Sall1",
                    "Selplg", "Siglech", "Olfml3", "Gpr34"),
  DAM           = c("Tyrobp", "Apoe", "B2m", "Ctsb", "Ctsd", "Fth1", "Lyz2",          # DAM-s1 (Trem2-indep)
                    "Trem2", "Cst7", "Lpl", "Cd9", "Itgax", "Clec7a", "Spp1",         # DAM-s2 (Trem2-dep)
                    "Gpnmb", "Igf1", "Axl", "Cd63"),
  IFN           = c("Ifit1", "Ifit2", "Ifit3", "Irf7", "Oasl2", "Isg15", "Mx1",
                    "Ifitm3", "Usp18", "Bst2", "Rsad2", "Stat1"),
  Proliferative = c("Mki67", "Top2a", "Birc5", "Mcm5", "Stmn1", "Cenpa"),
  MHC_APC       = c("Cd74", "H2-Aa", "H2-Ab1", "H2-Eb1", "H2-K1", "Tap1", "B2m")
)

# Subpopulation signatures that drive the cluster/cell ARGMAX (subset of canonical_microglia_markers;
# MHC_APC stays auxiliary -- scored, never argmax'd).
microglia_subpopulation_levels <- c("Homeostatic", "DAM", "IFN", "Proliferative")

# Non-microglial contamination signatures (UCell QC; P1-S2 prune_contaminant_clusters). A cluster
# scoring high here AND low on microglia_identity_markers is dropped. snRNA carries pervasive
# low-level ambient signal from these lineages -> the prune compares identity-vs-contaminant at
# the CLUSTER level on absolute (raw) scores; it never thresholds the pervasive per-cell background.
contam_signatures <- list(
  Oligo  = c("Mbp", "Plp1", "Mog", "Mobp", "Mag"),
  Neuron = c("Snap25", "Stmn2", "Rbfox3", "Syt1", "Meg3"),
  Astro  = c("Aqp4", "Gfap", "Slc1a2", "Slc1a3", "Aldoc")
)

# Raw-data paths (resolve through the storage/data symlink -> external read-only copy).
# Single source of truth; _targets.R registers each as a file target, loaders take the
# path as an explicit arg.
data_paths <- list(
  snrnaseq   = "storage/data/snrnaseq.rds",
  geomx      = "storage/data/geomx.rds",
  proteomics = "storage/data/proteomics_nonfiltered_nonnormalised.tsv",
  phospho    = "storage/data/phosphoproteomics_nonfiltered_nonnormalised.tsv",
  sample_key = "storage/data/proteomics_sample_key.csv"
)

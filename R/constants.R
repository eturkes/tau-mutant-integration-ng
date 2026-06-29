# Project constants: factor levels, colours, contrast definitions, canonical marker
# lists, RBC-contamination markers, and raw-data paths. Sourced first by tar_source();
# downstream R/ files + targets read these symbols. Values carried byte-exact from v1
# (genotype labels, hex colours, marker sets are load-bearing -- do not let them drift).

# Canonical genotype order -> factor levels everywhere downstream (MAPTKI = reference).
genotype_levels <- c("MAPTKI", "P301S", "NLGF_MAPTKI", "NLGF_P301S")

# Per-genotype colours (hex; carried from the v1 bslib palette).
genotype_colours <- c(
  MAPTKI      = "#7FA8C1",
  P301S       = "#B295C1",
  NLGF_MAPTKI = "#E08754",
  NLGF_P301S  = "#B0344D"
)

# Pairwise contrasts (numerator, denominator) driving the divergence narrative. The 5th
# canonical contrast `interaction` = (NLGF_P301S - P301S) - (NLGF_MAPTKI - MAPTKI) is a
# difference-of-differences, not a level pair -> built in R/design.R (S3), keyed by name.
contrast_definitions <- list(
  nlgf_in_maptki = c("NLGF_MAPTKI", "MAPTKI"),
  nlgf_in_p301s  = c("NLGF_P301S",  "P301S"),
  tau_alone      = c("P301S",       "MAPTKI"),
  tau_in_nlgf    = c("NLGF_P301S",  "NLGF_MAPTKI")
)

# Canonical mouse-microglia marker symbols. First 4 lists (Microglia/DAM/IFN/
# Proliferative) score transcriptional substates via Seurat::AddModuleScore (P1); the 3
# *_contam lists flag non-microglial contamination for QC purity checks.
canonical_microglia_markers <- list(
  Microglia     = c("Cx3cr1", "P2ry12", "Tmem119", "Trem2", "Csf1r", "Hexb"),
  DAM           = c("Cst7", "Apoe", "Lpl", "Itgax", "Spp1", "Tyrobp"),
  IFN           = c("Ifit3", "Isg15", "Stat1", "Oasl2"),
  Proliferative = c("Mki67", "Top2a", "Cenpf"),
  Oligo_contam  = c("Mbp", "Plp1", "Mog"),
  Neuron_contam = c("Snap25", "Stmn2", "Rbfox3"),
  Astro_contam  = c("Aqp4", "Gfap", "Slc1a2")
)

# Red-blood-cell / haemoglobin markers (mouse). Microglia lack adult haemoglobin, so
# these appearing among hub/module genes signal erythroid contamination (P3+ hdWGCNA).
rbc_marker_symbols <- c(
  "Hba-a1", "Hba-a2", "Hbb-b1", "Hbb-b2", "Hbb-bs", "Hbb-bt", "Hbb-y",
  "Hbb-bh1", "Hbb-bh2", "Hbq1a", "Hbq1b",
  "Alas2",
  "Slc4a1", "Epb41", "Epb42", "Ank1",
  "Klf1", "Gata1", "Bcl11a",
  "Bpgm"
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

# Project-wide constants: factor levels, colours, contrast definitions,
# canonical marker gene lists, and contamination-marker lists. Sourced
# first; every other R/ file depends on at least one symbol from here.

# Canonical genotype order used everywhere downstream.
genotype_levels <- c("MAPTKI", "P301S", "NLGF_MAPTKI", "NLGF_P301S")

genotype_colours <- c(
  MAPTKI      = "#7FA8C1",
  P301S       = "#B295C1",
  NLGF_MAPTKI = "#E08754",
  NLGF_P301S  = "#B0344D"
)

# Contrasts that drive the divergence narrative.
contrast_definitions <- list(
  nlgf_in_maptki = c("NLGF_MAPTKI", "MAPTKI"),
  nlgf_in_p301s  = c("NLGF_P301S",  "P301S"),
  tau_alone      = c("P301S",       "MAPTKI"),
  tau_in_nlgf    = c("NLGF_P301S",  "NLGF_MAPTKI")
)

# Canonical mouse-microglia marker gene-symbol lists used to score
# transcriptional states via Seurat::AddModuleScore. Single source of
# truth: both the Rmd's substate chunk and the per-state NEBULA build
# script call `label_microglia_states()` (R/microglia.R).
canonical_microglia_markers <- list(
  Microglia     = c("Cx3cr1", "P2ry12", "Tmem119", "Trem2", "Csf1r", "Hexb"),
  DAM           = c("Cst7", "Apoe", "Lpl", "Itgax", "Spp1", "Tyrobp"),
  IFN           = c("Ifit3", "Isg15", "Stat1", "Oasl2"),
  Proliferative = c("Mki67", "Top2a", "Cenpf"),
  Oligo_contam  = c("Mbp", "Plp1", "Mog"),
  Neuron_contam = c("Snap25", "Stmn2", "Rbfox3"),
  Astro_contam  = c("Aqp4", "Gfap", "Slc1a2")
)

# Canonical red-blood-cell / haemoglobin marker symbols. Used to flag
# likely RBC-contamination genes in microglia hub lists. Source: mouse
# haematology marker reviews + Tabula Muris erythroid cluster markers.
# Genuine microglia do not express adult haemoglobin; appearance of these
# symbols among module hubs signals erythroid contamination of nuclei
# preparations or droplet doublets.
rbc_marker_symbols <- c(
  # Adult and embryonic haemoglobin subunits
  "Hba-a1", "Hba-a2", "Hbb-b1", "Hbb-b2", "Hbb-bs", "Hbb-bt", "Hbb-y",
  "Hbb-bh1", "Hbb-bh2", "Hbq1a", "Hbq1b",
  # Haem biosynthesis (erythroid-restricted ALAS-2 and friends)
  "Alas2",
  # Erythroid-restricted membrane / cytoskeleton proteins
  "Slc4a1", "Epb41", "Epb42", "Ank1",
  # Erythroid transcription factors and globin regulators
  "Klf1", "Gata1", "Bcl11a",
  # Other haem / iron handling enriched in RBC lineage
  "Bpgm"
)

#!/usr/bin/env Rscript
# scripts/build_human_validation_signatures.R
#
# Session H2 of storage/notes/human_validation_plan.md.
#
# Assembles the panel of MOUSE microglial signatures that carry the
# tau x amyloid interaction story, maps each to HUMAN orthologs via
# nichenetr (mouse -> human; the reverse of the project's usual
# human -> mouse convention), and persists a single light cache the
# downstream human-scoring sessions (H3/H4) consume against the SEA-AD
# microglia h5ads.
#
# Why these signatures (one row per H2-table source):
#   - DAM / amyloid programmes  : DAM_up/DAM_down (Keren-Shaul reference,
#       custom_microglia_states.rds), LDAM (Marschallinger 2020) and WAM
#       (Safaiyan 2021) from custom_microglia_ad.rds. The amyloid-response
#       sanity axis (DAM_up = positive amyloid main effect in H4).
#   - 4 substate-defining marker sets : the EXACT canonical marker lists
#       (R/constants.R::canonical_microglia_markers) that the mouse
#       pipeline scored with Seurat::AddModuleScore and argmax'd to label
#       every nucleus homeostatic/DAM/IFN/proliferative (R/microglia.R).
#       There is NO de-novo FindAllMarkers/presto table for the four
#       states -- these curated lists ARE the substate definitions, so
#       H3 reproduces the mouse argmax in human by scoring exactly them.
#   - Interaction-carrying programmes : MG-M3 hdWGCNA module + hubs;
#       IFN P301S-specific (and shared) amyloid-response asymmetry genes;
#       the NF-kB target union (CollecTRI mouse, NF-kB-family regulon --
#       the I1-locked headline set behind the tau-driven NF-kB
#       attenuation); Gsk3b kinase substrates (OmniPath KSN).
#   - Interaction-DIRECTION signatures : the per-state interaction
#       contrast is FDR-null at the whole-gene level (n_sig = 0 in every
#       state; see storage/results/nebula_per_state_summary.tsv), so the
#       directional sets are taken from the project's curated interaction-
#       localisation table (nebula_per_state_localisation.tsv: the 62
#       genes whose interaction effect localises to a single state),
#       split by the sign of that state's interaction logFC. Encoded
#       both per state and pooled.
#   - Human-native positive controls : Gerrits 2022 AD1/AD2 ORIGINAL
#       human symbols (no cross-species mapping) as the pipeline sanity
#       anchor (anti-anchoring #5). A drift-guard asserts that mapping
#       these embedded human lists back to mouse reproduces the cached
#       mouse AD1/AD2 in custom_microglia_ad.rds byte-for-byte.
#
# Mapping convention (mirrors scripts/build_custom_microglia_ad.R, but
# reversed): drop NA + empty, dedup, sort. Coverage (% of input mouse
# symbols receiving a human ortholog) is reported per signature so
# dropouts are visible; signatures < 60% mapped are flagged low_coverage.
#
# Idempotent: skips the cache if present unless `--overwrite` is passed.
#
# Outputs (storage/cache/):
#   human_validation_signatures.rds              list(human, mouse, meta, nfkb_family)
#   human_validation_signatures_provenance.txt
# Outputs (storage/results/):
#   human_validation_signature_membership.tsv    one row per signature

suppressPackageStartupMessages({
  library(dplyr)
})

args      <- commandArgs(trailingOnly = TRUE)
overwrite <- "--overwrite" %in% args

setwd("/home/rstudio/tau-mutant-integration-ng")
source("R/helpers.R")  # attaches nichenetr (lazy geneinfo data); defines
                       # cache_or_run, write_tsv_safe, build_omnipath_ksn_mouse,
                       # canonical_microglia_markers.

cache_dir   <- "storage/cache"
results_dir <- "storage/results"
dir.create(cache_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

sig_path  <- file.path(cache_dir,   "human_validation_signatures.rds")
tsv_path  <- file.path(results_dir, "human_validation_signature_membership.tsv")
prov_path <- file.path(cache_dir,   "human_validation_signatures_provenance.txt")

# =========================================================================
# Embedded Gerrits 2022 human symbols (positive controls).
# Copied VERBATIM from scripts/build_custom_microglia_ad.R (Table S7); a
# runtime drift-guard below ties them to the canonical cached mouse sets.
# =========================================================================
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

# NF-kB family TFs whose CollecTRI mouse regulons form the union target
# set (identical to rmd/17_nfkb_attenuation.Rmd:73). Rel has 0 mouse
# targets and contributes nothing; A0A979HLR9 is the CollecTRI
# heterodimer-resolved NF-kB complex accession.
NFKB_FAMILY <- c("Rela", "Nfkb1", "Nfkb2", "Rel", "Relb", "A0A979HLR9")

SUBSTATES <- c("homeostatic", "DAM", "IFN", "proliferative")

# =========================================================================
# mouse -> human ortholog mapping (project NA/empty-drop + dedup + sort
# convention, applied to convert_mouse_to_human_symbols).
# =========================================================================
clean_syms <- function(x) sort(unique(x[!is.na(x) & nzchar(x)]))

map_mouse_to_human <- function(mouse_syms) {
  mouse_in <- clean_syms(mouse_syms)
  if (length(mouse_in) == 0L) {
    return(list(human = character(0), n_input = 0L,
                n_mapped_inputs = 0L, n_human = 0L))
  }
  raw   <- nichenetr::convert_mouse_to_human_symbols(mouse_in)  # named by input
  keep  <- !is.na(raw) & nzchar(raw)
  human <- sort(unique(unname(raw[keep])))
  list(human = human, n_input = length(mouse_in),
       n_mapped_inputs = sum(keep), n_human = length(human))
}

# Strip the trailing _residue token from a KSN target id ("Gsk3b_Y216" ->
# "Gsk3b"); mirrors R/kinase_inference.R::.ksn_substrate_to_symbol.
ksn_target_to_symbol <- function(x) sub("_[^_]+$", "", as.character(x))

# =========================================================================
# Build (wrapped in cache_or_run; the OmniPath KSN fetch lives inside).
# =========================================================================
build_signatures <- function() {

  message("Loading H2 signature sources")
  states7   <- readRDS(file.path(cache_dir,   "custom_microglia_states.rds"))
  ad4       <- readRDS(file.path(cache_dir,   "custom_microglia_ad.rds"))
  collectri <- readRDS(file.path(cache_dir,   "collectri_mouse_for_nfkb_gsea.rds"))
  hdwgcna   <- read.delim(file.path(results_dir, "hdwgcna_modules.tsv"),
                          check.names = FALSE, stringsAsFactors = FALSE)
  hubs      <- read.delim(file.path(results_dir, "hdwgcna_module_hubs.tsv"),
                          check.names = FALSE, stringsAsFactors = FALSE)
  ifn_asym  <- read.delim(file.path(results_dir, "ifn_nlgf_asymmetry.tsv"),
                          check.names = FALSE, stringsAsFactors = FALSE)
  loc       <- read.delim(file.path(results_dir, "nebula_per_state_localisation.tsv"),
                          check.names = FALSE, stringsAsFactors = FALSE)

  # ---- drift-guard: embedded human AD1/AD2 -> mouse must equal cache ----
  remap <- function(h) clean_syms(nichenetr::convert_human_to_mouse_symbols(h))
  stopifnot(setequal(remap(AD1_human), ad4$AD1),
            setequal(remap(AD2_human), ad4$AD2))
  message("  drift-guard OK: embedded Gerrits AD1/AD2 reproduce cached mouse sets")

  # ---- registries -------------------------------------------------------
  mouse_sigs <- list()
  meta_rows  <- list()
  add_sig <- function(name, mouse_genes, group,
                      type = "mapped", direction = NA_character_,
                      state = NA_character_, source = "") {
    mouse_sigs[[name]] <<- clean_syms(mouse_genes)
    meta_rows[[name]]  <<- data.frame(
      signature = name, group = group, type = type,
      state = state, direction = direction, source = source,
      stringsAsFactors = FALSE)
  }

  # DAM / amyloid programmes ------------------------------------------------
  add_sig("DAM_up",   states7$DAM_up,   "amyloid_dam",
          source = "custom_microglia_states.rds$DAM_up (Keren-Shaul reference DAM, up)")
  add_sig("DAM_down", states7$DAM_down, "amyloid_dam", direction = "down",
          source = "custom_microglia_states.rds$DAM_down (reference DAM, down)")
  add_sig("LDAM", ad4$LDAM, "amyloid_dam",
          source = "custom_microglia_ad.rds$LDAM (Marschallinger 2020)")
  add_sig("WAM",  ad4$WAM,  "amyloid_dam",
          source = "custom_microglia_ad.rds$WAM (Safaiyan 2021)")

  # 4 substate-defining marker sets (AddModuleScore->argmax labels) --------
  add_sig("substate_homeostatic", canonical_microglia_markers$Microglia,
          "substate_marker", state = "homeostatic",
          source = "R/constants.R canonical_microglia_markers$Microglia")
  add_sig("substate_DAM", canonical_microglia_markers$DAM,
          "substate_marker", state = "DAM",
          source = "R/constants.R canonical_microglia_markers$DAM")
  add_sig("substate_IFN", canonical_microglia_markers$IFN,
          "substate_marker", state = "IFN",
          source = "R/constants.R canonical_microglia_markers$IFN")
  add_sig("substate_proliferative", canonical_microglia_markers$Proliferative,
          "substate_marker", state = "proliferative",
          source = "R/constants.R canonical_microglia_markers$Proliferative")

  # Interaction-carrying programmes ----------------------------------------
  add_sig("MG_M3_module", hdwgcna$symbol[hdwgcna$module == "MG-M3"],
          "interaction_programme",
          source = "hdwgcna_modules.tsv (module == MG-M3)")
  add_sig("MG_M3_hubs", hubs$symbol[hubs$module == "MG-M3"],
          "interaction_programme",
          source = "hdwgcna_module_hubs.tsv (module == MG-M3)")
  add_sig("IFN_asym_P301S_specific",
          ifn_asym$symbol[ifn_asym$category == "P301S-specific"],
          "interaction_programme", state = "IFN",
          source = "ifn_nlgf_asymmetry.tsv (category == P301S-specific)")
  add_sig("IFN_asym_shared",
          ifn_asym$symbol[ifn_asym$category == "Shared"],
          "interaction_programme", state = "IFN",
          source = "ifn_nlgf_asymmetry.tsv (category == Shared)")
  add_sig("NFKB_union_targets",
          unique(collectri$target[collectri$source %in% NFKB_FAMILY]),
          "interaction_programme",
          source = "collectri_mouse_for_nfkb_gsea.rds (union of NF-kB-family regulons)")

  message("Building OmniPath mouse KSN for Gsk3b substrates")
  ksn_m <- build_omnipath_ksn_mouse()
  gsk3b_substrates <- unique(ksn_target_to_symbol(
    ksn_m$target[ksn_m$source == "Gsk3b"]))
  add_sig("Gsk3b_targets", gsk3b_substrates, "interaction_programme",
          source = "OmniPath KSN build_omnipath_ksn_mouse() (substrate genes of Gsk3b)")

  # Interaction-DIRECTION signatures (localisation table, signed) ----------
  # Each localised gene gets the interaction logFC of its argmax state.
  # Most genes (36/62) do not localise to a single state
  # (max_negLogP_state == NA) and are correctly excluded from every set.
  loc$loc_logFC <- vapply(seq_len(nrow(loc)), function(i) {
    s <- loc$max_negLogP_state[i]
    if (is.na(s)) return(NA_real_)
    v <- loc[[paste0(s, "_logFC")]][i]
    if (length(v) != 1L) return(NA_real_)
    as.numeric(v)
  }, numeric(1))
  for (s in SUBSTATES) {
    rows <- loc[loc$max_negLogP_state == s & !is.na(loc$loc_logFC), ]
    add_sig(paste0("interaction_up_",   s), rows$symbol[rows$loc_logFC > 0],
            "interaction_direction", direction = "up", state = s,
            source = "nebula_per_state_localisation.tsv (localised to state, logFC>0)")
    add_sig(paste0("interaction_down_", s), rows$symbol[rows$loc_logFC < 0],
            "interaction_direction", direction = "down", state = s,
            source = "nebula_per_state_localisation.tsv (localised to state, logFC<0)")
  }
  add_sig("interaction_up_pooled",
          loc$symbol[!is.na(loc$loc_logFC) & loc$loc_logFC > 0],
          "interaction_direction", direction = "up", state = "pooled",
          source = "nebula_per_state_localisation.tsv (all states, logFC>0)")
  add_sig("interaction_down_pooled",
          loc$symbol[!is.na(loc$loc_logFC) & loc$loc_logFC < 0],
          "interaction_direction", direction = "down", state = "pooled",
          source = "nebula_per_state_localisation.tsv (all states, logFC<0)")

  # ---- map every mouse signature to human --------------------------------
  message(sprintf("Mapping %d mouse signatures to human orthologs", length(mouse_sigs)))
  mapped     <- lapply(mouse_sigs, map_mouse_to_human)
  human_list <- lapply(mapped, `[[`, "human")

  # ---- human-native positive controls (no mapping) -----------------------
  human_list[["Gerrits_AD1_human"]] <- clean_syms(AD1_human)
  human_list[["Gerrits_AD2_human"]] <- clean_syms(AD2_human)

  # ---- membership metadata ----------------------------------------------
  meta <- do.call(rbind, meta_rows)
  meta$n_mouse         <- vapply(meta$signature, function(n) mapped[[n]]$n_input,         integer(1))
  meta$n_mapped_inputs <- vapply(meta$signature, function(n) mapped[[n]]$n_mapped_inputs, integer(1))
  meta$n_human         <- vapply(meta$signature, function(n) mapped[[n]]$n_human,         integer(1))
  meta$pct_mapped      <- round(100 * meta$n_mapped_inputs / pmax(meta$n_mouse, 1L), 1)

  ctrl_meta <- data.frame(
    signature = c("Gerrits_AD1_human", "Gerrits_AD2_human"),
    group = "positive_control", type = "human_native_control",
    state = NA_character_, direction = NA_character_,
    source = "Gerrits 2022 Table S7 (original human symbols, no mapping)",
    n_mouse = NA_integer_, n_mapped_inputs = NA_integer_,
    n_human = c(length(human_list[["Gerrits_AD1_human"]]),
                length(human_list[["Gerrits_AD2_human"]])),
    pct_mapped = NA_real_, stringsAsFactors = FALSE)
  meta <- rbind(meta, ctrl_meta)
  meta$low_coverage <- !is.na(meta$pct_mapped) & meta$pct_mapped < 60
  rownames(meta) <- NULL
  meta <- meta[, c("signature", "group", "type", "state", "direction",
                   "n_mouse", "n_mapped_inputs", "n_human", "pct_mapped",
                   "low_coverage", "source")]

  # mouse-side record keeps the human-native controls as their cached
  # mouse versions for provenance symmetry.
  mouse_record <- c(mouse_sigs,
                    list(Gerrits_AD1_human = ad4$AD1,
                         Gerrits_AD2_human = ad4$AD2))

  list(human = human_list, mouse = mouse_record, meta = meta,
       nfkb_family = NFKB_FAMILY)
}

signatures <- cache_or_run(sig_path, build_signatures(), overwrite = overwrite)

# =========================================================================
# Exports + smoke summary.
# =========================================================================
meta <- signatures$meta
write_tsv_safe(meta, tsv_path)
message(sprintf("\nWrote membership table: %s  (%d signatures)",
                tsv_path, nrow(meta)))

message("\n=== signature panel (mouse n -> human n, %% mapped) ===")
for (i in seq_len(nrow(meta))) {
  r <- meta[i, ]
  flag <- if (isTRUE(r$low_coverage)) "  [LOW COVERAGE <60%]" else ""
  message(sprintf("  %-26s %-22s n_mouse=%-5s n_human=%-5d %s%s",
                  r$signature, r$group,
                  ifelse(is.na(r$n_mouse), "NA", as.character(r$n_mouse)),
                  r$n_human,
                  ifelse(is.na(r$pct_mapped), " (native)",
                         sprintf("%.0f%%", r$pct_mapped)),
                  flag))
}

low <- meta$signature[which(meta$low_coverage)]
if (length(low) > 0L) {
  message(sprintf("\n%d signature(s) below 60%% ortholog coverage (interpret with caution): %s",
                  length(low), paste(low, collapse = ", ")))
}
empties <- meta$signature[meta$n_human == 0L]
if (length(empties) > 0L) {
  message(sprintf("WARNING: %d signature(s) map to zero human genes: %s",
                  length(empties), paste(empties, collapse = ", ")))
}

# =========================================================================
# Provenance.
# =========================================================================
if (overwrite || !file.exists(prov_path)) {
  writeLines(c(
    "# human_validation_signatures.rds  --  provenance",
    sprintf("# Built: %s", format(Sys.time())),
    "# Script: scripts/build_human_validation_signatures.R",
    "# Plan:   storage/notes/human_validation_plan.md  (session H2)",
    "",
    "## Object structure",
    "#   human       : named list of HUMAN-symbol character vectors (score these in SEA-AD)",
    "#   mouse       : named list of MOUSE-symbol character vectors (pre-mapping originals)",
    "#   meta        : membership data.frame (also exported as the TSV below)",
    "#   nfkb_family : the NF-kB-family TF names behind NFKB_union_targets",
    "#",
    "# Exported TSV: storage/results/human_validation_signature_membership.tsv",
    "",
    "## Cross-species mapping",
    "# nichenetr::convert_mouse_to_human_symbols, post-processed with the",
    "# project idiom unique(x[!is.na(x) & nzchar(x)]) then sort. This is the",
    "# REVERSE of the project's usual human->mouse direction (the project",
    "# organism is mouse); no prior project code mapped mouse->human, so the",
    "# same NA/empty-drop+dedup convention is reused. pct_mapped = fraction",
    "# of input mouse symbols receiving a human ortholog; < 60% => low_coverage.",
    "",
    "## Signature groups",
    "# amyloid_dam          : DAM_up/DAM_down (Keren-Shaul reference),",
    "#                        LDAM (Marschallinger 2020), WAM (Safaiyan 2021).",
    "# substate_marker      : the four canonical marker lists",
    "#                        (R/constants.R::canonical_microglia_markers) that",
    "#                        the mouse pipeline AddModuleScore'd + argmax'd to",
    "#                        label nuclei homeostatic/DAM/IFN/proliferative",
    "#                        (R/microglia.R). No de-novo per-state marker table",
    "#                        exists; these lists ARE the state definitions, so",
    "#                        H3 reproduces the mouse argmax in human by scoring",
    "#                        exactly them. Small by design (3-6 genes each).",
    "# interaction_programme: MG-M3 hdWGCNA module + hubs; IFN P301S-specific /",
    "#                        shared asymmetry; NFKB_union_targets (CollecTRI",
    "#                        mouse, family = Rela/Nfkb1/Nfkb2/Rel/Relb/",
    "#                        A0A979HLR9; the I1-locked NF-kB-attenuation set);",
    "#                        Gsk3b_targets (OmniPath KSN substrate genes).",
    "# interaction_direction: per-state + pooled up/down sets. The per-state",
    "#                        interaction contrast is FDR-NULL at the whole-gene",
    "#                        level (nebula_per_state_summary.tsv: n_sig = 0 in",
    "#                        every state), so these are taken from the curated",
    "#                        interaction-localisation table",
    "#                        (nebula_per_state_localisation.tsv: the 62 genes",
    "#                        whose interaction effect localises to one state),",
    "#                        split by the sign of that state's interaction",
    "#                        logFC. Several per-state sets are small; sizes are",
    "#                        in the membership TSV for transparency.",
    "# positive_control     : Gerrits 2022 AD1/AD2 ORIGINAL human symbols",
    "#                        (Table S7), no cross-species mapping -- the",
    "#                        pipeline sanity anchor (anti-anchoring #5). A",
    "#                        build-time drift-guard asserts that mapping the",
    "#                        embedded human lists back to mouse reproduces the",
    "#                        cached mouse AD1/AD2 in custom_microglia_ad.rds.",
    "",
    "## Notes / caveats",
    "# - NFKB_union_targets includes a minority of UniProt accessions (e.g.",
    "#   A0A...) inherited from the CollecTRI mouse 'target' column; these are",
    "#   not mouse symbols and drop during mapping, lowering pct_mapped. This",
    "#   is expected, not an error.",
    "# - Gsk3b_targets round-trips through mouse (the OmniPath KSN is natively",
    "#   human and was human->mouse-mapped for the mouse analysis); kinase",
    "#   substrate sites are conserved so the loss is small. The human-native",
    "#   GSK3B substrate set (fetch_omnipath_ksn_human, source == GSK3B) is an",
    "#   available alternative if a future session prefers it.",
    "# - Data-universe intersection (which mapped human genes are present in",
    "#   the SEA-AD microglia var) is deferred to H3, where the h5ads load.",
    "",
    "## Per-signature membership",
    sprintf("#   %-26s %-22s %8s %8s %7s", "signature", "group",
            "n_mouse", "n_human", "%mapped"),
    paste0("#   ", apply(meta, 1, function(r) {
      sprintf("%-26s %-22s %8s %8s %7s", r["signature"], r["group"],
              ifelse(is.na(r["n_mouse"]), "NA", r["n_mouse"]),
              r["n_human"],
              ifelse(is.na(r["pct_mapped"]), "native", r["pct_mapped"]))
    }))
  ), prov_path)
  message(sprintf("Wrote provenance: %s", prov_path))
} else {
  message(sprintf("[skip] provenance already present: %s", prov_path))
}

message("\ndone")

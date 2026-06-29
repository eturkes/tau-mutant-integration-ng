# Input / identifier helpers: sample-key loader for proteomics/phospho,
# TSV column-to-key matcher, and symbol -> Ensembl map.

# Map proteomics/phospho sample columns to genotype using the sample key.
proteomics_sample_meta <- function(key_path, n_keep = 16,
                                   col_pattern = "PTM.Quantity") {
  key <- readr::read_csv(key_path, show_col_types = FALSE) |>
    dplyr::rename(file_name = `File name`, label = `Sample/Condtion`) |>
    dplyr::slice_head(n = n_keep)
  key$genotype <- dplyr::recode(key$label,
    "MAPT-KI_24M"      = "MAPTKI",
    "P301S+3_24M"      = "P301S",
    "NLGF-MAPT-KI_24M" = "NLGF_MAPTKI",
    "NLGF-P301S+3_24M" = "NLGF_P301S"
  )
  key$genotype <- factor(key$genotype, levels = genotype_levels)
  # The TSVs have ".raw.PTM.Quantity" or ".PTM.Quantity" suffixes; build matcher.
  key$col_stub <- sub("\\.PTM\\.Quantity$", "", key$file_name)
  key
}

# Match TSV intensity columns to sample key rows.
match_intensity_columns <- function(col_names, key, allow_raw = TRUE) {
  stubs <- sub("\\.PTM\\.Quantity$", "", col_names)
  stubs <- sub("\\.raw$", "", stubs)
  key_stubs <- key$col_stub
  m <- match(stubs, key_stubs)
  out <- data.frame(
    column   = col_names,
    stub     = stubs,
    key_idx  = m,
    genotype = ifelse(is.na(m), NA, as.character(key$genotype[m])),
    label    = ifelse(is.na(m), NA, key$label[m]),
    stringsAsFactors = FALSE
  )
  out
}

# Map a vector of gene symbols to the ENSEMBL row IDs used by the
# snRNAseq Seurat object. Returns a named character vector of IDs,
# preserving symbol names; missing symbols are dropped silently.
symbols_to_ensembl <- function(syms, symbol_map) {
  hit <- symbol_map$ensembl[match(syms, symbol_map$symbol)]
  setNames(hit[!is.na(hit)], syms[!is.na(hit)])
}

# Data loaders + proteomics sample-key helpers. v1 kept these inline in rmd/01_data.Rmd;
# factored here into pure, testable functions consumed by the data-load targets. Each
# loader takes an explicit path (default from data_paths in R/constants.R).

# --- snRNAseq ----------------------------------------------------------------

# Load the full snRNAseq Seurat object, subset to microglia, strip the whole-object SCT
# assay + reductions (P1 recomputes both on the subset), and free the ~8G parent. Returns
# an RNA-counts-only Seurat object (~26k microglia) carrying the original meta.data
# (genotype, batch, genotype_batch, QC cols), with genotype refactored to canonical order.
# Fails loud if the subset column or the 4x4 fully-crossed genotype_batch design breaks.
load_snrnaseq <- function(path = data_paths$snrnaseq) {
  sc <- readRDS(path)
  stopifnot("broad_annotations" %in% colnames(sc@meta.data))
  micro <- subset(sc, subset = broad_annotations == "Microglia")
  rm(sc)
  invisible(gc())
  SeuratObject::DefaultAssay(micro) <- "RNA"
  micro[["SCT"]] <- NULL
  micro@reductions <- list()
  md <- micro@meta.data
  # Design contract: 4 genotypes x 4 batches = 16 fully-crossed replicate ids, none missing.
  stopifnot(
    all(c("genotype", "batch", "genotype_batch") %in% colnames(md)),
    !anyNA(md$genotype), !anyNA(md$batch),
    dplyr::n_distinct(md$genotype_batch) == 16L,
    all(table(md$genotype, md$batch) > 0)
  )
  micro$genotype <- factor(as.character(micro$genotype), levels = genotype_levels)
  invisible(gc())
  micro
}

# Build the ensembl<->symbol map from @misc$geneids (symbols, positionally aligned to the
# RNA assay's ensembl rownames). Returns data.frame {ensembl, symbol}, one row per RNA
# feature. Map symbol->ensembl via symbols_to_ensembl(); ensembl (RNA rownames) is the key.
build_symbol_map <- function(seurat_obj) {
  ens <- rownames(seurat_obj[["RNA"]])
  syms <- as.character(seurat_obj@misc$geneids)
  stopifnot(length(ens) == length(syms))
  data.frame(ensembl = ens, symbol = syms, stringsAsFactors = FALSE)
}

# Map gene symbols -> ensembl ids via symbol_map (drops symbols with no ensembl hit).
# Returns a named character vector: names = matched input symbols, values = ensembl ids.
symbols_to_ensembl <- function(syms, symbol_map) {
  hit <- symbol_map$ensembl[match(syms, symbol_map$symbol)]
  stats::setNames(hit[!is.na(hit)], syms[!is.na(hit)])
}

# --- GeoMx spatial -----------------------------------------------------------

# Load the GeoMx WTA spatial Seurat object; refactor genotype to canonical order if present.
load_geomx <- function(path = data_paths$geomx) {
  geomx <- readRDS(path)
  if (!is.null(geomx$genotype)) {
    geomx$genotype <- factor(as.character(geomx$genotype), levels = genotype_levels)
  }
  geomx
}

# --- Proteomics / phosphoproteomics ------------------------------------------

# Read a raw Spectronaut PTM export (annotation columns + per-sample intensity columns).
# Returns the tibble verbatim; intensity columns map to samples via match_intensity_columns().
# `na` covers Spectronaut's empty/Filtered sentinels so intensities parse numeric cleanly.
read_spectronaut_tsv <- function(path) {
  readr::read_tsv(path, na = c("", "NA", "NaN", "Filtered"), show_col_types = FALSE)
}

# Parse the proteomics/phospho sample key. Keeps the first n_keep rows (= the 24M
# timepoint: 4 genotypes x 4 reps), remaps raw condition labels to canonical genotypes,
# and derives col_stub (the intensity-column key, sans the .PTM.Quantity suffix).
# Returns a tibble {file_name, label, genotype (factor), col_stub}.
proteomics_sample_meta <- function(path = data_paths$sample_key, n_keep = 16L) {
  key <- readr::read_csv(path, show_col_types = FALSE)
  key <- dplyr::rename(key, file_name = `File name`, label = `Sample/Condtion`)
  key <- dplyr::slice_head(key, n = n_keep)
  remap <- c(
    "MAPT-KI_24M"      = "MAPTKI",
    "P301S+3_24M"      = "P301S",
    "NLGF-MAPT-KI_24M" = "NLGF_MAPTKI",
    "NLGF-P301S+3_24M" = "NLGF_P301S"
  )
  key$genotype <- factor(unname(remap[key$label]), levels = genotype_levels)
  key$col_stub <- sub("\\.PTM\\.Quantity$", "", key$file_name)
  stopifnot(
    nrow(key) == n_keep,
    !anyNA(key$genotype),
    dplyr::n_distinct(key$genotype) == 4L
  )
  key
}

# Map a Spectronaut export's intensity columns to sample metadata. Strips the
# .PTM.Quantity (and an optional trailing .raw -> proteomics has .raw.PTM.Quantity, phospho
# has .PTM.Quantity) suffix to a stub, then matches to key$col_stub. Returns data.frame
# {column, stub, key_idx, genotype, label}; non-intensity columns kept with NA genotype.
match_intensity_columns <- function(col_names, key) {
  stub <- sub("\\.PTM\\.Quantity$", "", col_names)
  stub <- sub("\\.raw$", "", stub)
  idx <- match(stub, key$col_stub)
  data.frame(
    column   = col_names,
    stub     = stub,
    key_idx  = idx,
    genotype = as.character(key$genotype)[idx],
    label    = key$label[idx],
    stringsAsFactors = FALSE
  )
}

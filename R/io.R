# Data loaders + proteomics sample-key helpers. v1 kept these inline in rmd/01_data.Rmd;
# factored here into pure, testable functions consumed by the data-load targets. Each
# loader takes an explicit path (default from data_paths in R/constants.R).

# --- snRNAseq ----------------------------------------------------------------

# Load the full snRNAseq Seurat object, subset to microglia, strip the whole-object SCT
# assay + reductions + neighbor graphs (P1 recomputes them on the subset), and free the ~8G parent. Returns
# an RNA-counts-only Seurat object (~26k microglia) carrying the original meta.data
# (genotype, batch, genotype_batch, QC cols), with genotype refactored to canonical order.
# Fails loud if the subset column is missing or the 4x4 fully-crossed genotype_batch design
# breaks (4 genotypes x 4 batches = 16 cells, one replicate id per cell, nothing missing).
load_snrnaseq <- function(path = data_paths$snrnaseq) {
  sc <- readRDS(path)
  stopifnot("broad_annotations" %in% colnames(sc@meta.data))
  micro <- subset(sc, subset = broad_annotations == "Microglia")
  rm(sc)
  invisible(gc())
  SeuratObject::DefaultAssay(micro) <- "RNA"
  micro[["SCT"]] <- NULL
  micro@reductions <- list()
  micro@graphs <- list()       # drop stale SCT_nn/SCT_snn shadows (P1 recomputes); keep raw target clean
  micro@neighbors <- list()
  stopifnot(all(c("genotype", "batch", "genotype_batch") %in% colnames(micro@meta.data)))
  # Refactor genotype to canonical order BEFORE asserting: factor() coerces any value outside
  # genotype_levels to NA, so the no-NA check must see the refactored column to actually catch
  # an unexpected genotype (asserting the raw column would miss a level dropped by factor()).
  micro$genotype <- factor(as.character(micro$genotype), levels = genotype_levels)
  md <- micro@meta.data
  # Design contract. n_distinct(genotype_batch)==16 alone is insufficient (NA counts as a
  # level; it proves neither the 4x4 cross nor one-id-per-cell) -> assert the full shape:
  # no NA in any design col, 4 genotypes x 4 batches, 16 populated cells, 16 ids bijecting cells.
  stopifnot(
    !anyNA(md$genotype), !anyNA(md$batch), !anyNA(md$genotype_batch),
    dplyr::n_distinct(md$genotype) == 4L,
    dplyr::n_distinct(md$batch) == 4L,
    nrow(unique(md[, c("genotype", "batch")])) == 16L,                       # 4x4 fully crossed, every cell populated
    dplyr::n_distinct(md$genotype_batch) == 16L,
    nrow(unique(md[, c("genotype", "batch", "genotype_batch")])) == 16L      # exactly one genotype_batch id per cell
  )
  invisible(gc())
  micro
}

# Build the ensembl<->symbol map from @misc$geneids (symbols, positionally aligned to the
# RNA assay's ensembl rownames). Returns data.frame {ensembl, symbol}, one row per RNA
# feature. Map symbol->ensembl via symbols_to_ensembl(); ensembl (RNA rownames) is the key.
# Positional alignment is the v1 @misc$geneids contract; equal length is necessary but not
# sufficient to prove it (no external ground-truth pairing here to fully verify). Enforce
# what is provable: equal length, no missing symbols, unique symbols (so symbols_to_ensembl's
# match() is unambiguous, not a silent first-hit).
build_symbol_map <- function(seurat_obj) {
  ens <- rownames(seurat_obj[["RNA"]])
  syms <- as.character(seurat_obj@misc$geneids)
  stopifnot(
    length(ens) == length(syms),
    !anyNA(syms),
    anyDuplicated(syms) == 0L
  )
  data.frame(ensembl = ens, symbol = syms, stringsAsFactors = FALSE)
}

# Map gene symbols -> ensembl ids via symbol_map (drops symbols with no ensembl hit).
# symbol_map$symbol is unique (asserted in build_symbol_map) -> match() is unambiguous.
# Returns a named character vector: names = matched input symbols, values = ensembl ids.
symbols_to_ensembl <- function(syms, symbol_map) {
  hit <- symbol_map$ensembl[match(syms, symbol_map$symbol)]
  stats::setNames(hit[!is.na(hit)], syms[!is.na(hit)])
}

# --- GeoMx spatial -----------------------------------------------------------

# Load the GeoMx WTA spatial Seurat object; refactor genotype to canonical order. genotype
# is required (every AOI carries it) -> fail loud if the column is absent or any value falls
# outside the 4 canonical levels (factor() would otherwise coerce it to NA silently).
load_geomx <- function(path = data_paths$geomx) {
  geomx <- readRDS(path)
  stopifnot("genotype" %in% colnames(geomx@meta.data))
  geomx$genotype <- factor(as.character(geomx$genotype), levels = genotype_levels)
  stopifnot(!anyNA(geomx$genotype), dplyr::n_distinct(geomx$genotype) == 4L)
  geomx
}

# --- Proteomics / phosphoproteomics ------------------------------------------

# Read a raw Spectronaut PTM export (annotation columns + per-sample intensity columns).
# `na` covers Spectronaut's empty/Filtered sentinels so intensities parse numeric cleanly;
# fail loud on any parse problem (turns the one-time smoke check into an enforced contract).
# Strips readr's spec/problems attrs (external pointers that go stale -> bad_weak_ptr after
# the tibble is qs-serialized + restored from the targets store) -> returns a plain tibble.
read_spectronaut_tsv <- function(path) {
  x <- readr::read_tsv(path, na = c("", "NA", "NaN", "Filtered"), show_col_types = FALSE)
  probs <- readr::problems(x)
  if (nrow(probs) > 0L) {
    stop("read_spectronaut_tsv: ", nrow(probs), " parse problem(s) in ", path, call. = FALSE)
  }
  attr(x, "spec") <- NULL
  attr(x, "problems") <- NULL
  class(x) <- c("tbl_df", "tbl", "data.frame")
  x
}

# Normalise a Spectronaut intensity-column name (or sample-key file name) to a bare run
# stub: strip the .PTM.Quantity suffix and an optional trailing .raw. proteomics columns are
# <run>.raw.PTM.Quantity, phospho columns are <run>.PTM.Quantity, key file names carry
# neither -> all collapse to the same <run>. Shared by the producer (proteomics_sample_meta)
# and consumer (match_intensity_columns) so they normalise identically (no asymmetric strip).
normalise_ptm_stub <- function(x) sub("\\.raw$", "", sub("\\.PTM\\.Quantity$", "", x))

# Parse the proteomics/phospho sample key. Keeps the first n_keep rows (= the 24M
# timepoint: 4 genotypes x 4 reps), remaps raw condition labels to canonical genotypes,
# and derives col_stub (the intensity-column join key). Returns a tibble
# {file_name, label, genotype (factor), col_stub}. Fails loud unless the 24M design is intact:
# exactly the 4 expected labels, balanced 4 reps/genotype, unique join keys.
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
  key$col_stub <- normalise_ptm_stub(key$file_name)
  key <- key[, c("file_name", "label", "genotype", "col_stub")]
  # n_distinct(genotype)==4 alone passes a degenerate 13/1/1/1 split -> assert the balanced
  # table, the exact label set, and key uniqueness (match_intensity_columns first-hits col_stub).
  stopifnot(
    nrow(key) == n_keep,
    !anyNA(key$genotype),
    setequal(key$label, names(remap)),
    all(table(key$genotype) == n_keep %/% 4L),
    anyDuplicated(key$file_name) == 0L,
    anyDuplicated(key$col_stub) == 0L
  )
  key
}

# Map a Spectronaut export's intensity columns to sample metadata: normalise each column
# name to its run stub (normalise_ptm_stub, identical to the key's) then match to key$col_stub.
# Returns data.frame {column, stub, key_idx, genotype, label}; non-intensity (annotation)
# columns are kept with NA genotype (a meaningful "not a sample column", not a failure). The
# join is well-defined only if key$col_stub is unique -> assert it (match() first-hits).
# P4 consumes this to attach genotype to intensity columns and there asserts 16/16 matched.
match_intensity_columns <- function(col_names, key) {
  stopifnot(anyDuplicated(key$col_stub) == 0L)
  stub <- normalise_ptm_stub(col_names)
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

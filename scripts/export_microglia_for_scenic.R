#!/usr/bin/env Rscript
# export_microglia_for_scenic.R
#
# SCENIC arc (K2): export the microglia snRNAseq expression matrix from the
# Seurat object into a language-neutral sparse bundle that the pyscenic env
# (Python 3.10, project-local micromamba) consumes to build a loom and run
# GRNBoost2. This is the R half of the H3-style jsonlite/anndata bridge: R
# owns the Seurat .rds + the canonical Ensembl->MGI symbol map; Python owns
# pyscenic. See storage/notes/scenic_regulons_plan.md (K2).
#
# What it does (locked K1 gate decisions, do not re-derive):
#   1. Load microglia_seurat_processed.rds RNA *raw counts* (NOT SCT).
#   2. Identity check (anti-anchoring caveat c): report pct-positive for
#      microglial vs neuronal/astro/oligo markers; assert the cell set is the
#      same 26,104 nuclei / genotype tally as the locked section-14 DE input.
#   3. Map ENSMUSG -> MGI symbol via the project snrnaseq_symbol_map cache
#      (same convention as section 14/18: symbol_map$symbol[match(gene, ens)]).
#      Drop unmapped rows; collapse duplicate symbols by SUMMING counts
#      (sparse incidence-matrix multiply). Report coverage.
#   4. Gene pre-filter: keep symbols detected (count > 0) in >= 1% of cells.
#   5. Write the sparse bundle to storage/cache/scenic/ for the Python side.
#
# Inputs:
#   storage/cache/microglia_seurat_processed.rds  (Seurat; RNA + SCT)
#   storage/cache/snrnaseq_symbol_map.rds         (ensembl, symbol)
#
# Outputs (storage/cache/scenic/, all gitignored via /storage/*):
#   microglia_counts.mtx     genes x cells sparse integer counts (MatrixMarket)
#   microglia_genes.txt      MGI symbols, matrix row order (1 per line)
#   microglia_cells.txt      cell IDs, matrix col order (1 per line)
#   microglia_colattrs.tsv   CellID, genotype, batch, genotype_batch (per cell)
#   export_provenance.txt    coverage / identity / filter report (human + LLM)
#
# Run: Rscript scripts/export_microglia_for_scenic.R
suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(Matrix)
})
project_root <- "/home/rstudio/tau-mutant-integration-ng"
setwd(project_root)
source("R/helpers.R")

cache_dir  <- "storage/cache"
out_dir    <- file.path(cache_dir, "scenic")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Locked facts from the plan (assert, do not assume silently).
N_CELLS_EXPECTED <- 26104L
GENO_EXPECTED <- c(MAPTKI = 4349L, NLGF_MAPTKI = 8787L,
                   NLGF_P301S = 8583L, P301S = 4385L)
MIN_CELL_FRAC <- 0.01   # gene pre-filter: detected in >= 1% of cells

log_lines <- character(0)
emit <- function(...) {
  msg <- sprintf(...)
  cat(msg, "\n", sep = "")
  log_lines[[length(log_lines) + 1L]] <<- msg
  flush.console()
}

t0 <- Sys.time()
emit("== K2 export: microglia -> SCENIC sparse bundle ==")
emit("started %s", format(t0, "%Y-%m-%d %H:%M:%S"))

# ---- 1. Load object + counts ----------------------------------------------
emit("loading microglia_seurat_processed.rds ...")
micro <- readRDS(file.path(cache_dir, "microglia_seurat_processed.rds"))
counts <- GetAssayData(micro, assay = "RNA", layer = "counts")
counts <- as(counts, "CsparseMatrix")
emit("RNA counts: %d genes x %d cells (class %s)",
     nrow(counts), ncol(counts), class(counts)[1])
stopifnot(ncol(counts) == N_CELLS_EXPECTED)

meta <- micro@meta.data
stopifnot(all(c("genotype", "batch", "genotype_batch") %in% colnames(meta)))
stopifnot(identical(rownames(meta), colnames(counts)))

# ---- 2. Identity check (caveat c) -----------------------------------------
symbol_map <- readRDS(file.path(cache_dir, "snrnaseq_symbol_map.rds"))
stopifnot(all(c("ensembl", "symbol") %in% colnames(symbol_map)))

geno_tally <- table(factor(as.character(meta$genotype),
                           levels = names(GENO_EXPECTED)))
emit("genotype tally: %s",
     paste(sprintf("%s=%d", names(geno_tally), as.integer(geno_tally)),
           collapse = ", "))
if (!identical(as.integer(geno_tally), as.integer(GENO_EXPECTED))) {
  emit("WARNING: genotype tally != locked section-14 DE tally (%s)",
       paste(sprintf("%s=%d", names(GENO_EXPECTED), GENO_EXPECTED),
             collapse = ", "))
} else {
  emit("identity-of-input OK: %d nuclei + genotype tally byte-match the locked section-14 DE input",
       N_CELLS_EXPECTED)
}

markers <- list(
  microglia = c("Csf1r", "P2ry12", "Hexb", "Itgam", "Cx3cr1", "Tmem119"),
  neuronal  = c("Snap25", "Rbfox3", "Syt1"),
  astro     = c("Gfap", "Aqp4"),
  oligo     = c("Plp1", "Mbp")
)
emit("-- marker pct-positive (raw counts, Ensembl-mapped) --")
for (cls in names(markers)) {
  ens <- symbols_to_ensembl(markers[[cls]], symbol_map)
  ens <- ens[ens %in% rownames(counts)]
  if (!length(ens)) { emit("  %-9s: (no markers found in matrix)", cls); next }
  for (sym in names(ens)) {
    row <- counts[ens[[sym]], ]
    pct <- 100 * Matrix::nnzero(row) / length(row)
    emit("  %-9s %-8s pct+=%5.1f%% mean=%.3f", cls, sym, pct, mean(row))
  }
}

# ---- 3. Ensembl -> MGI symbol map + duplicate-symbol collapse -------------
sym_for_row <- symbol_map$symbol[match(rownames(counts), symbol_map$ensembl)]
n_total <- nrow(counts)
keep_mapped <- !is.na(sym_for_row) & nzchar(sym_for_row)
emit("symbol coverage: %d / %d Ensembl rows mapped (%.1f%%); %d dropped (unmapped/empty)",
     sum(keep_mapped), n_total, 100 * sum(keep_mapped) / n_total,
     sum(!keep_mapped))

counts <- counts[keep_mapped, , drop = FALSE]
sym_for_row <- sym_for_row[keep_mapped]

dup_syms <- unique(sym_for_row[duplicated(sym_for_row)])
emit("duplicate symbols (multi-Ensembl): %d symbols collapsed by summing counts",
     length(dup_syms))

# Sparse collapse: (symbol x gene incidence) %*% (gene x cell) = symbol x cell.
sym_fac <- factor(sym_for_row, levels = unique(sym_for_row))
incidence <- Matrix::fac2sparse(sym_fac, to = "d")  # n_symbol x n_gene (0/1, numeric)
collapsed <- incidence %*% counts               # n_symbol x n_cell
collapsed <- as(collapsed, "CsparseMatrix")
rownames(collapsed) <- levels(sym_fac)
emit("after collapse: %d unique symbols x %d cells", nrow(collapsed), ncol(collapsed))

# ---- 4. Gene pre-filter: detected in >= 1% of cells ------------------------
min_cells <- ceiling(MIN_CELL_FRAC * ncol(collapsed))
det <- Matrix::rowSums(collapsed > 0)
keep_gene <- det >= min_cells
emit("gene pre-filter: keep symbols detected in >= %d cells (%.0f%% of %d) -> %d / %d kept (%d dropped)",
     min_cells, 100 * MIN_CELL_FRAC, ncol(collapsed),
     sum(keep_gene), length(keep_gene), sum(!keep_gene))
mat <- collapsed[keep_gene, , drop = FALSE]
mat <- as(mat, "CsparseMatrix")
emit("final SCENIC input matrix: %d genes x %d cells (nnz=%d, %.2f%% dense)",
     nrow(mat), ncol(mat), length(mat@x),
     100 * length(mat@x) / (as.double(nrow(mat)) * ncol(mat)))

# TF-list overlap sanity (how many regressors are even present).
tf_path <- "storage/data/cistarget/allTFs_mm.txt"
if (file.exists(tf_path)) {
  tfs <- readLines(tf_path, warn = FALSE)
  emit("TFs (allTFs_mm.txt) present in matrix: %d / %d",
       sum(tfs %in% rownames(mat)), length(tfs))
}

# ---- 5. Write the sparse bundle -------------------------------------------
mtx_path   <- file.path(out_dir, "microglia_counts.mtx")
genes_path <- file.path(out_dir, "microglia_genes.txt")
cells_path <- file.path(out_dir, "microglia_cells.txt")
attrs_path <- file.path(out_dir, "microglia_colattrs.tsv")
prov_path  <- file.path(out_dir, "export_provenance.txt")

Matrix::writeMM(mat, mtx_path)
writeLines(rownames(mat), genes_path)
writeLines(colnames(mat), cells_path)
col_attrs <- data.frame(
  CellID         = colnames(mat),
  genotype       = as.character(meta$genotype),
  batch          = as.character(meta$batch),
  genotype_batch = as.character(meta$genotype_batch),
  stringsAsFactors = FALSE
)
write.table(col_attrs, attrs_path, sep = "\t", quote = FALSE, row.names = FALSE)

emit("wrote bundle:")
emit("  %s", mtx_path)
emit("  %s", genes_path)
emit("  %s", cells_path)
emit("  %s", attrs_path)
emit("total wall time: %.1f s", as.numeric(difftime(Sys.time(), t0, units = "secs")))

writeLines(log_lines, prov_path)
for (p in c(mtx_path, genes_path, cells_path, attrs_path, prov_path)) {
  Sys.chmod(p, "0664")
}
emit("provenance -> %s", prov_path)

# io contract tests (deferred from S2 review -> S3). S2 verified shapes against the LIVE data
# one-time; these commit reproducible, data-free contracts on synthetic fixtures so the loaders'
# fail-loud guards stay enforced. Covers the pure helpers (normalise_ptm_stub, symbols_to_ensembl,
# build_symbol_map, match_intensity_columns) and the loaders' assertions (proteomics_sample_meta,
# read_spectronaut_tsv, load_snrnaseq design-contract, load_geomx) via tempfiles.

source("R/constants.R")
source("R/io.R")
source("tests/helpers.R")

# --- normalise_ptm_stub: strip .PTM.Quantity then optional trailing .raw -----------------
stopifnot(
  normalise_ptm_stub("Foo.raw.PTM.Quantity") == "Foo",                 # proteomics: .raw.PTM.Quantity
  normalise_ptm_stub("Bar.PTM.Quantity") == "Bar",                     # phospho: .PTM.Quantity
  normalise_ptm_stub("Baz") == "Baz",                                  # key file name: neither
  normalise_ptm_stub("Naoto-Hippo_TiO2_DIA_NN.raw.PTM.Quantity") == "Naoto-Hippo_TiO2_DIA_NN",
  identical(normalise_ptm_stub(c("A.raw.PTM.Quantity", "B.PTM.Quantity", "C")), c("A", "B", "C"))
)

# --- symbols_to_ensembl: unambiguous match, drop misses, keep names ---------------------
sm <- data.frame(ensembl = c("E1", "E2", "E3"), symbol = c("A", "B", "C"), stringsAsFactors = FALSE)
stopifnot(identical(symbols_to_ensembl(c("A", "C", "Z"), sm), c(A = "E1", C = "E3")),
          length(symbols_to_ensembl("Z", sm)) == 0L)

# --- build_symbol_map: positional ensembl<->symbol, fail loud on broken contract ---------
obj <- make_fake_seurat(n_genes = 10L)
mapdf <- build_symbol_map(obj)
stopifnot(is.data.frame(mapdf), nrow(mapdf) == 10L,
          identical(mapdf$ensembl, rownames(obj[["RNA"]])),
          identical(mapdf$symbol, paste0("Sym", 1:10)))
o_dup <- obj;   o_dup@misc$geneids   <- c("Dup", "Dup", paste0("Sym", 3:10))   # duplicate symbol
o_short <- obj; o_short@misc$geneids <- paste0("Sym", 1:9)                      # length mismatch
o_na <- obj;    o_na@misc$geneids    <- c(NA, paste0("Sym", 2:10))             # missing symbol
expect_error(build_symbol_map(o_dup))
expect_error(build_symbol_map(o_short))
expect_error(build_symbol_map(o_na))

# --- match_intensity_columns: annotation -> NA, runs matched, dup stub fails ------------
key <- data.frame(col_stub = c("Run01", "Run02"),
                  genotype = factor(c("MAPTKI", "P301S"), levels = genotype_levels),
                  label    = c("MAPT-KI_24M", "P301S+3_24M"), stringsAsFactors = FALSE)
cols <- c("PG.Genes", "Run01.raw.PTM.Quantity", "Run02.PTM.Quantity", "Run99.PTM.Quantity")
mic  <- match_intensity_columns(cols, key)
stopifnot(nrow(mic) == 4L,
          is.na(mic$genotype[1]),                                # annotation column
          mic$genotype[2] == "MAPTKI", mic$genotype[3] == "P301S",
          is.na(mic$genotype[4]))                                # unmatched run
key_dup <- key; key_dup$col_stub <- c("Run01", "Run01")
expect_error(match_intensity_columns(cols, key_dup))

# --- proteomics_sample_meta: remap + balance/label/uniqueness asserts (synthetic CSV) ----
write_key_csv <- function(labels16, n_extra = 5L) {
  df <- rbind(
    data.frame(`File name` = paste0("Run", sprintf("%02d", seq_along(labels16))),
               `Sample/Condtion` = labels16, check.names = FALSE),
    data.frame(`File name` = paste0("Set6_", seq_len(n_extra)),
               `Sample/Condtion` = "Other_xx", check.names = FALSE)
  )
  p <- tempfile(fileext = ".csv"); readr::write_csv(df, p); p
}
good_labels <- rep(c("MAPT-KI_24M", "P301S+3_24M", "NLGF-MAPT-KI_24M", "NLGF-P301S+3_24M"), each = 4L)
key_ok <- proteomics_sample_meta(write_key_csv(good_labels), n_keep = 16L)
stopifnot(nrow(key_ok) == 16L, !anyNA(key_ok$genotype),
          identical(levels(key_ok$genotype), genotype_levels),
          all(table(key_ok$genotype) == 4L),
          all(c("file_name", "label", "genotype", "col_stub") %in% names(key_ok)),
          anyDuplicated(key_ok$col_stub) == 0L)
# unbalanced labels (5/3/4/4): all 4 labels present but not 4 each -> balance assert fails
unbalanced <- c(rep("MAPT-KI_24M", 5L), rep("P301S+3_24M", 3L),
                rep("NLGF-MAPT-KI_24M", 4L), rep("NLGF-P301S+3_24M", 4L))
expect_error(proteomics_sample_meta(write_key_csv(unbalanced), n_keep = 16L))
# unknown label -> genotype NA -> fails
bad_label <- good_labels; bad_label[1] <- "Unknown_24M"
expect_error(proteomics_sample_meta(write_key_csv(bad_label), n_keep = 16L))

# --- read_spectronaut_tsv: parses, strips readr attrs, plain tibble ----------------------
tsv <- data.frame(`PG.Genes` = c("A", "B"), `Run1.raw.PTM.Quantity` = c(1.5, 2.5), check.names = FALSE)
ptsv <- tempfile(fileext = ".tsv"); readr::write_tsv(tsv, ptsv)
x <- read_spectronaut_tsv(ptsv)
stopifnot(is.data.frame(x), nrow(x) == 2L,
          is.null(attr(x, "spec")), is.null(attr(x, "problems")),
          identical(class(x), c("tbl_df", "tbl", "data.frame")))

# --- load_snrnaseq: subset + design-contract asserts -------------------------------------
sn <- make_fake_seurat(with_broad = TRUE, with_sct = TRUE, n_other = 8L)  # 64 microglia + 8 neurons
psn <- tempfile(fileext = ".rds"); saveRDS(sn, psn)
micro <- load_snrnaseq(psn)
stopifnot(ncol(micro) == 64L,                                  # neurons subset out
          all(micro$broad_annotations == "Microglia"),
          !"SCT" %in% SeuratObject::Assays(micro),             # SCT dropped
          SeuratObject::DefaultAssay(micro) == "RNA",
          identical(levels(micro$genotype), genotype_levels),
          dplyr::n_distinct(micro$genotype_batch) == 16L)
# missing broad_annotations
sn1 <- make_fake_seurat(with_sct = TRUE)
p1 <- tempfile(fileext = ".rds"); saveRDS(sn1, p1); expect_error(load_snrnaseq(p1))
# out-of-level genotype on a microglia cell -> NA after refactor -> fails
sn2 <- sn; sn2@meta.data$genotype[1] <- "WT"
p2 <- tempfile(fileext = ".rds"); saveRDS(sn2, p2); expect_error(load_snrnaseq(p2))
# subtle: 16 distinct genotype_batch but NOT bijective with (genotype, batch) -> last assert fires
sn3 <- make_fake_seurat(with_broad = TRUE, with_sct = TRUE)
i3 <- which(sn3@meta.data$genotype_batch == "MAPTKI_batch01")[1]
sn3@meta.data$genotype_batch[i3] <- "MAPTKI_batch02"           # reuse an existing id -> still 16 distinct
p3 <- tempfile(fileext = ".rds"); saveRDS(sn3, p3); expect_error(load_snrnaseq(p3))

# --- load_geomx: genotype refactor + fail-loud -------------------------------------------
gx <- make_fake_seurat()                                       # RNA-only; genotype meta present
pg <- tempfile(fileext = ".rds"); saveRDS(gx, pg)
g_ok <- load_geomx(pg)
stopifnot(identical(levels(g_ok$genotype), genotype_levels),
          dplyr::n_distinct(g_ok$genotype) == 4L)
gx2 <- gx; gx2@meta.data$genotype <- NULL                      # missing genotype column
pgm <- tempfile(fileext = ".rds"); saveRDS(gx2, pgm); expect_error(load_geomx(pgm))
gx3 <- gx; gx3@meta.data$genotype[1] <- "WT"                   # out-of-level genotype
pgo <- tempfile(fileext = ".rds"); saveRDS(gx3, pgo); expect_error(load_geomx(pgo))

cat("ok - test_io\n")

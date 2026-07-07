options(warn = 2)
source("R/constants.R")
source("R/utils.R")
source("R/design.R")
source("R/io.R")
source("R/de_pb.R")
source("R/modality_de.R")
source("tests/helpers.R")

# --- positive_log2_matrix: nonpositive -> NA, then log2; provenance counts -------------------
m <- matrix(c(4, 0, -2, NA, 8, 16), nrow = 2)   # col-major: [,1]=4,0 [,2]=-2,NA [,3]=8,16
lg <- positive_log2_matrix(m)
cnt <- attr(lg, "log2_counts")
stopifnot(
  lg[1, 1] == 2, lg[1, 3] == 3, lg[2, 3] == 4,   # log2(4/8/16)
  is.na(lg[2, 1]), is.na(lg[1, 2]), is.na(lg[2, 2]),  # 0 and -2 -> NA; input NA stays NA
  cnt$n_nonpositive_to_na == 2L, cnt$n_missing_input == 1L, cnt$n_missing_output == 3L)
cat("ok - positive_log2_matrix sets nonpositive to NA before log2 with honest counts\n")

# --- protein_group_features + aggregate_proteome_raw: group sum, gene_first, present->NA ------
prot <- data.frame(
  PG.ProteinGroups = c("P1", "P1", "P2", ""),
  PG.Genes = c("Apoe", "Apoe;Trem2", "Mapt", ""),
  s1 = c(10, 20, 5, 99),
  s2 = c(0, 4, -1, 99),
  check.names = FALSE, stringsAsFactors = FALSE)
feat <- protein_group_features(prot)
fc <- attr(feat, "protein_group_counts")
stopifnot(
  identical(rownames(feat), c("P1", "P2")),          # empty PG dropped; radix-sorted groups
  identical(feat$gene_first, c("Apoe", "Mapt")),
  identical(feat$gene_symbols, c("Apoe;Trem2", "Mapt")),
  identical(feat$n_gene_symbols, c(2L, 1L)),
  identical(feat$n_raw_rows, c(2L, 1L)),
  fc$n_missing_protein_group == 1L, fc$n_duplicate_rows_by_group == 1L)
agg <- aggregate_proteome_raw(prot, c("s1", "s2"), feat)
stopifnot(
  identical(rownames(agg), feat$protein_group),
  agg["P1", "s1"] == 30, agg["P1", "s2"] == 4,        # P1 sums positives; s2 row1=0 -> NA, only 4 present
  agg["P2", "s1"] == 5, is.na(agg["P2", "s2"]))        # P2 s2=-1 -> NA -> no present value -> NA
cat("ok - protein_group_features + aggregate_proteome_raw sum positives per group, NA when none present\n")

# --- geomx_slide_design: cell-means + slide, full rank, bare-level names, 5 canonical contrasts
gmeta <- data.frame(
  genotype = factor(rep(genotype_levels, each = 2), levels = genotype_levels),
  slide = factor(rep(c("A", "B"), 4)),
  row.names = paste0("aoi", 1:8))
gd <- geomx_slide_design(gmeta, include_slide = TRUE)
stopifnot(
  qr(gd$design)$rank == ncol(gd$design),
  all(genotype_levels %in% colnames(gd$design)),      # genotype prefix stripped to bare levels
  identical(rownames(gd$contrasts), colnames(gd$design)),
  all(c("nlgf_in_maptki", "nlgf_in_p301s", "interaction") %in% colnames(gd$contrasts)))
gmeta1 <- gmeta; gmeta1$slide <- factor(rep("A", 8))
expect_error(geomx_slide_design(gmeta1, include_slide = TRUE), "slide")   # <2 slide levels fail loud
cat("ok - geomx_slide_design builds a full-rank cell-means design with the 5 canonical contrasts\n")

# --- geomx_spatial_descriptor: AOI coordinate score from top amyloid-response genes ----------
gmeta_sp <- transform(
  gmeta,
  bio_rep = rep(1:4, each = 2),
  roi = paste0("roi", 1:8),
  SampleID = paste0("s", 1:8),
  x = rep(c(0, 1), 4),
  y = rep(c(0, 0, 1, 1), 2),
  q3_factor = seq(1, 2, length.out = 8),
  neg_background = seq(0.1, 0.8, length.out = 8),
  nuclei = seq(100, 800, length.out = 8)
)
gcounts <- matrix(seq(10, 41), nrow = 4,
                  dimnames = list(paste0("G", 1:4), rownames(gmeta_sp)))
gtop <- list(
  nlgf_in_maptki = data.frame(symbol = paste0("G", 1:4), logFC = c(1, -2, 0.5, 0.1),
                              adj.P.Val = c(0.01, 0.02, 0.50, 0.80),
                              stringsAsFactors = FALSE),
  nlgf_in_p301s = data.frame(symbol = paste0("G", 1:4), logFC = c(2, -1, 0.4, 0.2),
                             adj.P.Val = c(0.01, 0.04, 0.40, 0.70),
                             stringsAsFactors = FALSE)
)
gs <- geomx_spatial_descriptor(gcounts, gmeta_sp, gtop, top_n = 2L)
stopifnot(
  is.data.frame(gs$aoi), is.data.frame(gs$genes),
  nrow(gs$aoi) == ncol(gcounts), nrow(gs$genes) == 2L,
  all(c("x_coord", "y_coord", "signed_response_score", "score_abs") %in% names(gs$aoi)),
  all(is.finite(gs$aoi$signed_response_score)),
  identical(as.character(gs$aoi$genotype), as.character(gmeta_sp$genotype)),
  gs$provenance$n_score_genes == 2L)
cat("ok - geomx_spatial_descriptor builds finite AOI coordinate scores\n")

# --- match_24m_bulk_columns: 16/16 matched, balanced 4/genotype, imbalance fails loud ---------
key <- data.frame(
  file_name = paste0("run", 1:16),
  label = rep(c("MAPT-KI_24M", "P301S+3_24M", "NLGF-MAPT-KI_24M", "NLGF-P301S+3_24M"), each = 4),
  genotype = factor(rep(genotype_levels, each = 4), levels = genotype_levels),
  col_stub = paste0("run", 1:16),
  stringsAsFactors = FALSE)
tbl <- as.data.frame(setNames(as.list(rep(1, 16)), paste0("run", 1:16, ".PTM.Quantity")),
                     check.names = FALSE)
tbl$Annotation <- "x"                                  # a non-intensity column -> NA match, ignored
mm <- match_24m_bulk_columns(tbl, key, modality = "test")
stopifnot(length(mm$columns) == 16L, nrow(mm$meta) == 16L,
          all(as.integer(table(mm$meta$genotype)) == 4L),
          identical(mm$meta$run_index, 1:16))
key_bad <- key
key_bad$genotype <- factor(c(rep("MAPTKI", 5), rep("P301S", 3), rep("NLGF_MAPTKI", 4),
                             rep("NLGF_P301S", 4)), levels = genotype_levels)
expect_error(match_24m_bulk_columns(tbl, key_bad, modality = "test"), "balanced")
cat("ok - match_24m_bulk_columns asserts 16/16 matched columns balanced 4/genotype\n")

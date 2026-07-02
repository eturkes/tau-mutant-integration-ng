# P4-S1 GeoMx contracts: count/meta extraction, slide + block design, replicate-aware
# limma-voom fit, sensitivity status, and deconvolution preflight provenance.

source("R/constants.R")
source("R/design.R")
source("R/crossmodality.R")
source("tests/helpers.R")

canonical <- c("tau_alone", "nlgf_in_maptki", "nlgf_in_p301s", "tau_in_nlgf", "interaction")

make_fake_geomx <- function(n_genes = 80L, n_bio = 3L, aoi_per = 2L) {
  combos <- expand.grid(genotype = genotype_levels,
                        bio_rep = seq_len(n_bio),
                        aoi_rep = seq_len(aoi_per),
                        stringsAsFactors = FALSE)
  combos$slide_rep <- ((match(combos$genotype, genotype_levels) + combos$bio_rep + combos$aoi_rep) %% 4L) + 1L
  combos$roi <- paste0("roi", seq_len(nrow(combos)))
  combos$SampleID <- paste(combos$genotype, combos$bio_rep, combos$aoi_rep, sep = "_")
  combos$`ROI Coordinate X` <- seq_len(nrow(combos)) * 10
  combos$`ROI Coordinate Y` <- seq_len(nrow(combos)) * 5
  combos$q_norm_qFactors <- 0.8 + (seq_len(nrow(combos)) %% 5L) / 10
  combos$NegGeoMean_Mm_R_NGS_WTA_v1.0 <- 4 + (seq_len(nrow(combos)) %% 7L)
  combos$nuclei <- ifelse(seq_len(nrow(combos)) %% 11L == 0L, -1, 30 + seq_len(nrow(combos)))
  rownames(combos) <- paste0("AOI", seq_len(nrow(combos)))

  gidx <- match(combos$genotype, genotype_levels)
  counts <- 20 + outer(seq_len(n_genes), seq_len(nrow(combos)),
                       function(i, j) (i * 7L + j * 11L) %% 37L)
  counts[seq_len(12L), ] <- counts[seq_len(12L), ] +
    rep(c(0, 3, 5, 9)[gidx], each = 12L)
  counts[1, ] <- 0L
  rownames(counts) <- paste0("Gene", seq_len(n_genes))
  colnames(counts) <- rownames(combos)
  obj <- SeuratObject::CreateSeuratObject(counts = Matrix::Matrix(counts, sparse = TRUE),
                                          meta.data = combos,
                                          min.cells = 0, min.features = 0)
  obj[["SCT"]] <- SeuratObject::CreateAssay5Object(counts = Matrix::Matrix(counts + 1000L, sparse = TRUE))
  SeuratObject::DefaultAssay(obj) <- "SCT"
  attr(obj, "raw_counts") <- counts
  obj
}

gx <- make_fake_geomx()
raw_counts <- attr(gx, "raw_counts")

# --- count/meta extraction: explicit RNA counts despite SCT default, zero-gene drop -----
cnt <- geomx_count_matrix(gx)
cp <- attr(cnt, "geomx_count_provenance")
stopifnot(nrow(cnt) == nrow(raw_counts) - 1L, ncol(cnt) == ncol(raw_counts),
          identical(colnames(cnt), colnames(raw_counts)),
          cnt[1, 1] == raw_counts[2, 1],                    # proves RNA counts, not SCT+1000
          cp$n_genes_dropped_empty == 1L,
          cp$n_non_integer == 0L, cp$coerced_integer,
          identical(storage.mode(cnt), "integer"))

meta <- geomx_meta(gx)
mp <- attr(meta, "geomx_meta_provenance")
stopifnot(identical(rownames(meta), colnames(cnt)),
          all(c("genotype", "slide", "bio_unit", "roi", "SampleID", "x", "y",
                "q3_factor", "neg_background", "nuclei") %in% names(meta)),
          nlevels(meta$bio_unit) == 4L * 3L,
          mp$nuclei_sentinel_count == sum(meta$nuclei < 0))

# malformed metadata fails at the named missing column, not later in model fitting
gx_bad <- gx
gx_bad@meta.data$bio_rep <- NULL
expect_error(geomx_meta(gx_bad), "bio_rep")

# --- slide design rank guard -------------------------------------------------------------
fd <- geomx_slide_design(meta)
stopifnot(identical(colnames(fd$contrasts), canonical),
          all(genotype_levels %in% colnames(fd$design)))
meta_rank_bad <- meta
meta_rank_bad$slide <- factor(as.character(meta_rank_bad$genotype))
expect_error(geomx_slide_design(meta_rank_bad), "rank-deficient")

# --- replicate-aware fit: duplicateCorrelation primary + sensitivity branches ------------
de <- fit_geomx_de(cnt, meta, min_count = 1L)
stopifnot(de$primary$status == "fit",
          de$primary$duplicate_correlation$used,
          is.finite(de$primary$duplicate_correlation$consensus_correlation),
          identical(names(de$primary$top), canonical),
          de$sensitivity$unblocked$status == "fit",
          !de$sensitivity$unblocked$duplicate_correlation$used,
          de$sensitivity$collapsed_bio_unit$status == "fit",
          de$sensitivity$collapsed_bio_unit$n_bio_units == 12L)
for (tt in de$primary$top) {
  stopifnot(is.data.frame(tt), nrow(tt) == de$primary$kept,
            all(c("symbol", "contrast", "logFC", "P.Value", "adj.P.Val", "t", "CI.L", "CI.R") %in% names(tt)),
            length(unique(tt$contrast)) == 1L)
}

miss <- meta$genotype != "NLGF_P301S"
skip <- fit_geomx_collapsed_sensitivity(cnt[, miss, drop = FALSE], meta[miss, , drop = FALSE], min_count = 1L)
stopifnot(skip$status == "skipped", grepl("rank-deficient", skip$reason, fixed = TRUE))

one_unit <- meta$bio_rep == "1"
sat <- fit_geomx_collapsed_sensitivity(cnt[, one_unit, drop = FALSE], meta[one_unit, , drop = FALSE], min_count = 1L)
stopifnot(sat$status == "skipped", grepl("no residual degrees", sat$reason, fixed = TRUE))

# --- decon preflight: specific defer/block/earned reasons --------------------------------
sp_ok <- list(package = "SpatialDecon", available = TRUE, version = "0.0.0",
              repos = "synthetic", error = NA_character_, warnings = character(), messages = character())
pf <- geomx_decon_preflight(meta, cnt, spatialdecon = sp_ok)
stopifnot(pf$status == "defer",
          any(grepl("reference profile", pf$reasons, fixed = TRUE)),
          pf$nuclei$n_sentinel > 0L,
          is.finite(pf$memory$estimated_peak_mb))

meta_bad_q3 <- meta
meta_bad_q3$q3_factor[1] <- 0
pf_bad <- geomx_decon_preflight(meta_bad_q3, cnt, spatialdecon = sp_ok)
stopifnot(pf_bad$status == "blocked",
          any(grepl("Q3 normalisation", pf_bad$reasons, fixed = TRUE)))

profile_bad <- cbind(A = 1:6, B = (1:6) * 2, C = c(6, 1, 5, 2, 4, 3))
pf_col <- geomx_decon_preflight(meta, cnt, profile = profile_bad, spatialdecon = sp_ok)
stopifnot(pf_col$status == "blocked",
          any(grepl("collinearity", pf_col$reasons, fixed = TRUE)))

profile_ok <- cbind(A = c(1, 2, 3, 4, 5, 6),
                    B = c(6, 5, 3, 1, 2, 4),
                    C = c(2, 6, 1, 5, 3, 4))
pf_ok <- geomx_decon_preflight(meta, cnt, profile = profile_ok, spatialdecon = sp_ok)
stopifnot(pf_ok$status == "earned", pf_ok$reference$profile_ok)

cat("ok - test_crossmodality\n")

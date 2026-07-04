options(warn = 2)
source("R/report.R")
source("tests/helpers.R")

tmp <- tempfile(fileext = ".html")
writeLines(
  paste0(
    '<html><body>',
    '<a href="index_files/figure-html/fig-a-1.png" class="lightbox" data-gallery="g">',
    '<img role="img" src="data:image/png;base64,AAA" class="img-fluid figure-img" alt="A"></a>',
    '<a href="index_files/figure-html/fig-b-1.png" class="lightbox" data-gallery="g">',
    '<img role="img" src="data:image/png;base64,BBB" class="img-fluid figure-img" alt="B"></a>',
    '</body></html>'
  ),
  tmp,
  useBytes = TRUE
)
res <- repair_embedded_lightbox(tmp)
txt <- paste(readLines(tmp, warn = FALSE), collapse = "\n")
stopifnot(
  identical(res$n_repaired, 2L),
  identical(res$n_local_before, 2L),
  identical(res$n_local_after, 0L),
  !grepl('href="index_files/figure-html/', txt, fixed = TRUE),
  grepl('href="data:image/png;base64,AAA"', txt, fixed = TRUE),
  grepl('href="data:image/png;base64,BBB"', txt, fixed = TRUE),
  grepl('src="data:image/png;base64,AAA"', txt, fixed = TRUE),
  grepl('src="data:image/png;base64,BBB"', txt, fixed = TRUE),
  # data-type="image" forces GLightbox to render the data: URI as an image, not a blank iframe
  grepl('<a data-type="image" href="data:image/png;base64,AAA"', txt, fixed = TRUE),
  grepl('<a data-type="image" href="data:image/png;base64,BBB"', txt, fixed = TRUE)
)
cat("ok - repair_embedded_lightbox rewrites local lightbox hrefs + forces image type\n")

# idempotent: re-running finds no local hrefs and does not double-inject data-type
res2 <- repair_embedded_lightbox(tmp)
txt2 <- paste(readLines(tmp, warn = FALSE), collapse = "\n")
stopifnot(
  identical(res2$n_repaired, 0L),
  identical(sum(gregexpr('data-type="image"', txt2, fixed = TRUE)[[1]] > 0L), 2L)
)
cat("ok - repair_embedded_lightbox idempotent (no double data-type)\n")

tmp_clean <- tempfile(fileext = ".html")
writeLines('<html><body><img src="data:image/png;base64,AAA"></body></html>', tmp_clean, useBytes = TRUE)
clean <- repair_embedded_lightbox(tmp_clean)
stopifnot(identical(clean$n_repaired, 0L), identical(clean$n_local_before, 0L))
cat("ok - repair_embedded_lightbox no-ops when no local hrefs exist\n")

tmp_bad <- tempfile(fileext = ".html")
writeLines(
  '<html><body><a href="index_files/figure-html/fig-a-1.png" class="lightbox"><img src="fig-a-1.png"></a></body></html>',
  tmp_bad,
  useBytes = TRUE
)
expect_error(repair_embedded_lightbox(tmp_bad), "matched 0 of 1 local figure href(s)")
cat("ok - repair_embedded_lightbox fails loud on unembedded lightbox shape\n")

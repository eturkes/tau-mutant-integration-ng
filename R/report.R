# Report render helpers. Keep generated-HTML repair here so the targets `report`
# file is exactly the artifact the user opens.

repair_embedded_lightbox <- function(html_file) {
  stopifnot(is.character(html_file), length(html_file) == 1L, file.exists(html_file))
  html <- paste(readLines(html_file, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  count_matches <- function(x) {
    if (length(x) == 1L && x[[1]] == -1L) 0L else length(x)
  }
  local_href_pattern <- 'href="[^"]+_files/figure-html/[^"]+\\.png"'
  local_href <- gregexpr(local_href_pattern, html, perl = TRUE)[[1]]
  n_local <- count_matches(local_href)
  if (n_local == 0L) {
    return(invisible(list(n_repaired = 0L, n_local_before = 0L, n_local_after = 0L)))
  }

  pattern <- paste0(
    '<a\\s+href="[^"]+_files/figure-html/[^"]+\\.png"',
    '[^>]*class="lightbox"[^>]*>\\s*',
    '<img[^>]*\\ssrc="data:image/png;base64,[^"]+"[^>]*>\\s*</a>'
  )
  matches <- gregexpr(pattern, html, perl = TRUE)
  n_matches <- count_matches(matches[[1]])
  if (n_matches != n_local) {
    stop("lightbox href repair matched ", n_matches, " of ", n_local,
         " local figure href(s)", call. = FALSE)
  }

  matched_anchors <- regmatches(html, matches)[[1]]
  repaired_anchors <- vapply(matched_anchors, function(anchor) {
    src <- sub('^.*\\ssrc="(data:image/png;base64,[^"]+)".*$',
               "\\1", anchor, perl = TRUE)
    if (identical(src, anchor)) {
      stop("lightbox href repair could not extract embedded src", call. = FALSE)
    }
    anchor <- sub(local_href_pattern,
                  paste0('href="', src, '"'), anchor, perl = TRUE)
    # GLightbox types each slide from the href's file extension; an embedded data: URI
    # has none, so GLightbox falls back to an <iframe> and the pop-out renders BLANK.
    # `data-type="image"` forces the image slide (empirically verified against the
    # bundled GLightbox). Idempotent guard so a re-run never double-injects.
    if (!grepl('data-type="image"', anchor, fixed = TRUE)) {
      anchor <- sub("<a\\s", '<a data-type="image" ', anchor, perl = TRUE)
    }
    anchor
  }, character(1))
  repaired <- html
  regmatches(repaired, matches) <- list(repaired_anchors)
  local_after <- gregexpr(local_href_pattern, repaired, perl = TRUE)[[1]]
  n_after <- count_matches(local_after)
  if (n_after != 0L) {
    stop("lightbox href repair left ", n_after, " local figure href(s)", call. = FALSE)
  }

  con <- file(html_file, open = "wb")
  on.exit(close(con), add = TRUE)
  writeBin(charToRaw(repaired), con)
  invisible(list(n_repaired = n_matches, n_local_before = n_local, n_local_after = n_after))
}

render_report <- function(report_sources,
                          report_extra_files,
                          microglia_report,
                          microglia_figures,
                          trajectory_figures,
                          modality_scatter_figures) {
  stopifnot(all(file.exists(c(report_sources, report_extra_files))))
  invisible(list(microglia_report, microglia_figures, trajectory_figures,
                 modality_scatter_figures))
  quarto::quarto_render(
    input = ".",
    execute = TRUE,
    execute_dir = getwd(),
    execute_daemon = 0,
    execute_daemon_restart = FALSE,
    execute_debug = FALSE,
    quiet = FALSE,
    as_job = FALSE
  )
  html_file <- file.path("report", "tau-mutant-integration.html")
  repair_embedded_lightbox(html_file)
  legacy_html <- file.path("report", "index.html")
  if (file.exists(legacy_html)) {
    unlink(legacy_html)
  }
  html_file
}

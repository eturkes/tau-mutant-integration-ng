# Report render helpers. Keep generated-HTML repair here so the targets `report`
# file is exactly the artifact the user opens.

repair_embedded_lightbox <- function(html_file) {
  stopifnot(is.character(html_file), length(html_file) == 1L, file.exists(html_file))
  html <- paste(readLines(html_file, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  count_matches <- function(x) {
    if (length(x) == 1L && x[[1]] == -1L) 0L else length(x)
  }
  local_href_pattern <- 'href="index_files/figure-html/[^"]+\\.png"'
  local_href <- gregexpr(local_href_pattern, html, perl = TRUE)[[1]]
  n_local <- count_matches(local_href)
  if (n_local == 0L) {
    return(invisible(list(n_repaired = 0L, n_local_before = 0L, n_local_after = 0L)))
  }

  pattern <- paste0(
    '<a\\s+href="index_files/figure-html/[^"]+\\.png"',
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
    sub('href="index_files/figure-html/[^"]+\\.png"',
        paste0('href="', src, '"'), anchor, perl = TRUE)
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
                          qc_figures,
                          microglia_report,
                          composition_results,
                          pb_de_microglia,
                          pb_de_substate,
                          symbol_map,
                          microglia_figures,
                          trajectory_report,
                          trajectory_figures,
                          mechanism_report,
                          mechanism_figures,
                          crossmodality_report,
                          crossmodality_figures) {
  stopifnot(all(file.exists(c(report_sources, report_extra_files))))
  invisible(list(
    qc_figures, microglia_report, composition_results, pb_de_microglia,
    pb_de_substate, symbol_map, microglia_figures, trajectory_report,
    trajectory_figures, mechanism_report, mechanism_figures,
    crossmodality_report, crossmodality_figures
  ))
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
  html_file <- file.path("_report", "index.html")
  repair_embedded_lightbox(html_file)
  html_file
}

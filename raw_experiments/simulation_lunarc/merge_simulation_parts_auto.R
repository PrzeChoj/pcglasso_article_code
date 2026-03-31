#!/usr/bin/env Rscript
#
# merge_simulation_parts_auto.R
#
# Auto-detect and merge all results_*_partXofY.rds and errors_*_partXofY.rds
# files in the current folder. Supports optional "_big" and "_n<...>" tags.
#

parse_part_file <- function(fname) {
  pattern <- "^(results|errors)_([A-Za-z0-9_]+)(?:_big)?(?:_n([A-Za-z0-9_]+))?_part([0-9]+)of([0-9]+)\\.rds$"
  m <- regexec(pattern, fname)
  reg <- regmatches(fname, m)[[1]]
  if (length(reg) == 0) {
    return(NULL)
  }
  n_tag <- if (!is.na(reg[[4]]) && nzchar(reg[[4]])) reg[[4]] else NULL
  list(
    kind = reg[[2]],
    mode = reg[[3]],
    n_tag = n_tag,
    part = as.integer(reg[[5]]),
    n_parts = as.integer(reg[[6]]),
    big = grepl("_big", fname),
    file = fname
  )
}

files <- list.files(pattern = "^(results|errors)_.*_part[0-9]+of[0-9]+\\.rds$")
if (length(files) == 0) {
  stop("No part files found in the current directory.")
}

info <- lapply(files, parse_part_file)
info <- Filter(Negate(is.null), info)
if (length(info) == 0) {
  stop("No matching part files found.")
}

keys <- vapply(info, function(x) {
  tag <- if (is.null(x$n_tag)) "none" else x$n_tag
  paste(x$kind, x$mode, if (x$big) "big" else "normal", tag, x$n_parts, sep = "|")
}, character(1))

groups <- split(info, keys)
for (key in names(groups)) {
  grp <- groups[[key]]
  kind <- grp[[1]]$kind
  mode <- grp[[1]]$mode
  big <- grp[[1]]$big
  n_tag <- grp[[1]]$n_tag
  n_parts <- grp[[1]]$n_parts
  parts_present <- vapply(grp, function(x) x$part, integer(1))
  missing <- setdiff(seq_len(n_parts), parts_present)
  if (length(missing)) {
    message("Skipping ", mode, if (big) " (big)" else "",
            if (!is.null(n_tag)) paste0(" n", n_tag) else "",
            ": missing parts ", paste(missing, collapse = ", "))
    next
  }

  ordered <- grp[order(parts_present)]
  merged <- do.call(rbind, lapply(ordered, function(x) readRDS(x$file)))
  outname <- paste0(kind, "_", mode, if (big) "_big" else "",
                    if (!is.null(n_tag)) paste0("_n", n_tag) else "", ".rds")
  saveRDS(merged, file = outname)
  message("Merged ", kind, " ", mode, if (big) " (big)" else "",
          if (!is.null(n_tag)) paste0(" n", n_tag) else "",
          " into ", outname, " from ", length(ordered), " parts")
}

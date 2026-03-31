#!/usr/bin/env Rscript
#
# merge_simulation_parts.R
#
# Merge part files from estimate_simulated.R into a single results file.
#
# Usage:
#   Rscript merge_simulation_parts.R TRUE [n_parts] [outname] [big]
#   Rscript merge_simulation_parts.R FALSE [n_parts] [outname] [big]
#   Rscript merge_simulation_parts.R pcglasso [n_parts] [outname] [big]
#   Rscript merge_simulation_parts.R glasso [n_parts] [outname] [big]
#
# Examples:
#   Rscript merge_simulation_parts.R TRUE 5
#   Rscript merge_simulation_parts.R pcglasso 5 results_pcglasso.rds
#   Rscript merge_simulation_parts.R TRUE 7 big
#   Rscript merge_simulation_parts.R pcglasso 7 results_pcglasso_big.rds --big
#

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Please provide TRUE/FALSE or pcglasso/glasso.")
}
if (length(args) > 4) {
  stop("Too many arguments. Expected: mode [n_parts] [outname] [big].")
}

is_big <- any(tolower(args) %in% c("big", "--big"))
args <- args[!tolower(args) %in% c("big", "--big")]

mode_arg <- tolower(args[[1]])
if (mode_arg %in% c("true", "false")) {
  generate.pcglasso <- as.logical(mode_arg)
  mode_tag <- if (generate.pcglasso) "pcglasso" else "glasso"
} else if (mode_arg %in% c("pcglasso", "glasso")) {
  mode_tag <- mode_arg
} else {
  stop("First argument must be TRUE/FALSE or pcglasso/glasso.")
}

n_parts <- 5L
if (length(args) >= 2) {
  n_parts <- as.integer(args[[2]])
}
if (is.na(n_parts) || n_parts < 1) {
  stop("'n_parts' must be a positive integer.")
}

outname <- if (length(args) >= 3) args[[3]] else paste0("results_", mode_tag, ".rds")
if (is_big && length(args) < 3) {
  outname <- paste0("results_", mode_tag, "_big.rds")
}

file_prefix <- if (is_big) {
  paste0("results_", mode_tag, "_big_part")
} else {
  paste0("results_", mode_tag, "_part")
}
files <- sprintf("%s%sof%s.rds", file_prefix, seq_len(n_parts), n_parts)
missing <- files[!file.exists(files)]
if (length(missing)) {
  stop("Missing part files: ", paste(missing, collapse = ", "))
}

parts <- lapply(files, readRDS)
res <- do.call(rbind, parts)
saveRDS(res, file = outname)
message("Merged ", length(parts), " part files into ", outname, "\n")

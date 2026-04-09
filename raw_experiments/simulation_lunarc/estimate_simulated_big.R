#!/usr/bin/env Rscript
#
# estimate_simulated_big.R
#
# Runs experiments with reduced threads and supports splitting reps across array jobs.
#
# Usage:
#   Rscript estimate_simulated_big.R TRUE [part] [n_parts] [n_list]
#   Rscript estimate_simulated_big.R FALSE [part] [n_parts] [n_list]
#
# Examples:
#   Rscript estimate_simulated_big.R TRUE
#   Rscript estimate_simulated_big.R TRUE 3 7
#   Rscript estimate_simulated_big.R TRUE 3 7 5000
#   Rscript estimate_simulated_big.R TRUE 1 7 200,300,500
#   Rscript estimate_simulated_big.R TRUE 1 7 all
#

suppressPackageStartupMessages({
  library(space)
  library(PCGLASSOcpp)
  library(PCGLASSO)
  library(pcglassoFast)
  source("estimation_methods.R")
  source("simulation_functions.R")
  library(glasso)
  library(parallel)
  library(Matrix)
})

quiet_eval <- function(expr) {
  res <- NULL
  suppressMessages(suppressWarnings(capture.output({
    res <- eval(expr)
  }, type = "output")))
  res
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Please provide TRUE (pcglasso) or FALSE (glasso).")
}
if (length(args) > 4) {
  stop("Too many arguments. Expected: TRUE/FALSE [part] [n_parts] [n_list].")
}

generate.pcglasso <- as.logical(args[[1]])
if (is.na(generate.pcglasso)) {
  stop("First argument must be either TRUE or FALSE.")
}
graph_type <- if (generate.pcglasso) "hub" else "non_hub"

part <- NULL
n_parts <- 7L
if (length(args) >= 2) {
  part <- as.integer(args[[2]])
}
if (length(args) >= 3) {
  n_parts <- as.integer(args[[3]])
}
if (!is.null(part) && (is.na(part) || is.na(n_parts))) {
  stop("'part' and 'n_parts' must be integers.")
}

parse_ns_arg <- function(arg, default_n, all_ns) {
  if (is.null(arg)) {
    return(list(ns = default_n, tag = paste0("n", default_n)))
  }
  arg <- gsub("\\s+", "", arg)
  if (tolower(arg) == "all") {
    return(list(ns = all_ns, tag = "nall"))
  }
  parts <- unlist(strsplit(arg, ",", fixed = TRUE))
  parts <- parts[parts != ""]
  ns <- as.integer(parts)
  if (length(ns) == 0 || any(is.na(ns))) {
    stop("n_list must be 'all' or a comma-separated list of integers.")
  }
  tag <- if (length(ns) == 1) {
    paste0("n", ns)
  } else {
    paste0("n", paste(ns, collapse = "_"))
  }
  list(ns = ns, tag = tag)
}

set.seed(2)
graphics.off()

split.train        <- 0.7
ns_all             <- c(200, 300, 500, 1000, 5000)
ns_default         <- 5000
sim                <- 20
nlambda            <- 50
mc_cores           <- 8L
alpha_grid         <- 0
lambda.min.ratio   <- 0.01

ns_arg <- NULL
if (length(args) >= 4) {
  ns_arg <- args[[4]]
}
ns_info <- parse_ns_arg(ns_arg, default_n = ns_default, all_ns = ns_all)
ns <- ns_info$ns
ns_tag <- ns_info$tag

out_dir <- Sys.getenv("OUTPUT_DIR", unset = ".")
if (is.na(out_dir) || out_dir == "") {
  out_dir <- "."
}
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

message("Running n = ", paste(ns, collapse = ", "), " with mc_cores = ", mc_cores)
message("graph_type = ", graph_type)
if (!is.null(part)) {
  message("Running rep part ", part, " of ", n_parts)
}
message("Output dir = ", out_dir)

if (!generate.pcglasso) {
  data(Q_simulated_glasso)
  Q <- Q_simulated_glasso
} else {
  data(Q_simulated_pcglasso)
  Q <- Q_simulated_pcglasso
}
Q <- (Q + t(Q)) / 2

reps_all <- seq_len(sim)
if (is.null(part)) {
  reps <- reps_all
} else {
  reps <- select_ns_part(reps_all, part = part, n_parts = n_parts)
}
message("Running reps: ", paste(reps, collapse = ", "))

run_args <- list(
  Q = Q,
  ns = ns,
  sim = sim,
  rep_ids = reps,
  mc_cores = mc_cores,
  nlambda = nlambda,
  lambda.min.ratio = lambda.min.ratio,
  alpha_grid = alpha_grid
)
if ("return_errors" %in% names(formals(run_experiments))) {
  run_args$return_errors <- TRUE
}
res <- quiet_eval(quote(do.call(run_experiments, run_args)))

if (is.list(res) && all(c("results", "errors") %in% names(res))) {
  results_df <- res$results
  errors_df <- res$errors
} else {
  results_df <- res
  errors_df <- data.frame(
    iter = integer(),
    n = integer(),
    rep = integer(),
    method = character(),
    selector = character(),
    solver = character(),
    seed = character(),
    stage = character(),
    error = character(),
    stringsAsFactors = FALSE
  )
}
results_df$graph_type <- rep(graph_type, nrow(results_df))
errors_df$graph_type <- rep(graph_type, nrow(errors_df))
if (is.null(part)) {
  outname <- paste0("results_", graph_type, "_big_", ns_tag, ".rds")
} else {
  outname <- paste0("results_", graph_type, "_big_", ns_tag, "_part", part, "of", n_parts, ".rds")
}
outpath <- file.path(out_dir, outname)
saveRDS(results_df, file = outpath)
message("Saved results to ", outpath)
errname <- sub("^results_", "errors_", outname)
errpath <- file.path(out_dir, errname)
saveRDS(errors_df, file = errpath)
message("Saved errors to ", errpath)

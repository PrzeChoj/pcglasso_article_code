#!/usr/bin/env Rscript
#
# estimate_simulated.R
#
# This script runs `run_experiments(...)` under two modes of `generate.pcglasso` (TRUE/FALSE),
# depending on the first command‐line argument.  It saves the resulting `res` object to disk
# (as an .rds file named “results_glasso.rds” or “results_pcglasso.rds”).
#
# Usage:
#   Rscript estimate_simulated.R FALSE [part] [n_parts]    # runs in the “glasso” mode
#   Rscript estimate_simulated.R TRUE [part] [n_parts]     # runs in the “pcglasso” mode
#
# Examples:
#   Rscript estimate_simulated.R TRUE          # all ns
#   Rscript estimate_simulated.R TRUE 1 5      # part 1 of 5
#   Rscript estimate_simulated.R TRUE 3        # part 3 of default 5
#

suppressPackageStartupMessages({
  library(space)
  library(PCGLASSOcpp)
  library(PCGLASSO) #devtools::install_github('JackStorrorCarter/PCGLASSO')

  library(pcglassoFast)
  source("estimation_methods.R")
  source("simulation_functions.R")  # make sure this is in the same folder or adjust the path
  library(glasso)
  library(parallel)
  library(Matrix)
})



# 1) Parse command‐line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Please provide TRUE (for pcglasso) or FALSE (for glasso).")
}
if (length(args) > 3) {
  stop("Too many arguments. Expected: TRUE/FALSE [part] [n_parts].")
}
generate.pcglasso <- as.logical(args[[1]])
if (is.na(generate.pcglasso)) {
  stop("First argument must be either TRUE or FALSE.")
}

part <- NULL
n_parts <- 5L
if (length(args) >= 2) {
  part <- as.integer(args[[2]])
}
if (length(args) >= 3) {
  n_parts <- as.integer(args[[3]])
}
if (!is.null(part) && (is.na(part) || is.na(n_parts))) {
  stop("'part' and 'n_parts' must be integers.")
}

message("Running with generate.pcglasso = ", generate.pcglasso, "\n")
if (!is.null(part)) {
  message("Running part ", part, " of ", n_parts, "\n")
}

# 2) Common settings (you can adjust these as needed)
set.seed(2)
graphics.off()

split.train        <- 0.7      # not used directly below, but could be referenced in your functions
ns_all             <- c(200, 300, 500, 1000, 5000)
ns                 <- select_ns_part(ns_all, part = part, n_parts = n_parts)
sim                <- 7
nlambda            <- 50
mc_cores           <- parallel::detectCores()
#mc_cores <- 1L
alpha_grid         <- 0
lambda.min.ratio   <- 0.01

# 3) Load the appropriate Q‐matrix depending on generate.pcglasso
if (!generate.pcglasso) {
  data(Q_simulated_glasso)      # from pcglassoFast or simulation_functions.R
  Q <- Q_simulated_glasso
} else {
  data(Q_simulated_pcglasso)
  Q <- Q_simulated_pcglasso
}
# ensure symmetry
Q <- (Q + t(Q)) / 2

# 4) Run the experiment
res <- run_experiments(
  Q = Q,
  ns          = ns,
  sim         = sim,
  mc_cores    = mc_cores,
  nlambda     = nlambda,
  lambda.min.ratio = lambda.min.ratio,
  alpha_grid  = alpha_grid
)

# 5) Save the result to a file whose name reflects the mode/part
mode_tag <- if (generate.pcglasso) "pcglasso" else "glasso"
if (is.null(part)) {
  outname <- paste0("results_", mode_tag, ".rds")
} else {
  outname <- paste0("results_", mode_tag, "_part", part, "of", n_parts, ".rds")
}
saveRDS(res, file = outname)
message("Saved results to ", outname, "\n")

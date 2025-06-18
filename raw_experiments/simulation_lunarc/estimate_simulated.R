#!/usr/bin/env Rscript
#
# estimate_simulated.R
#
# This script runs `run_experiments(...)` under two modes of `generate.pcglasso` (TRUE/FALSE),
# depending on the first command‐line argument.  It saves the resulting `res` object to disk
# (as an .rds file named “results_glasso.rds” or “results_pcglasso.rds”).
#
# Usage:
#   Rscript estimate_simulated.R FALSE    # runs in the “glasso” mode
#   Rscript estimate_simulated.R TRUE     # runs in the “pcglasso” mode
#

suppressPackageStartupMessages({
  library(space)
  library(pcglassoFast)
  source("simulation_functions.R")  # make sure this is in the same folder or adjust the path
  library(glasso)
  library(parallel)
  library(Matrix)
})



# 1) Parse the single command‐line argument to decide generate.pcglasso
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Please provide a single argument: TRUE (for pcglasso) or FALSE (for glasso).")
}
generate.pcglasso <- as.logical(args[[1]])
if (is.na(generate.pcglasso)) {
  stop("First argument must be either TRUE or FALSE.")
}

message("Running with generate.pcglasso = ", generate.pcglasso, "\n")

# 2) Common settings (you can adjust these as needed)
set.seed(2)
graphics.off()

split.train        <- 0.7      # not used directly below, but could be referenced in your functions
ns                 <- c(200,300,500,1000,5000)
sim                <- 200
nlambda            <- 50
mc_cores           <- parallel::detectCores()
#mc_cores <- 1L
alpha.grid         <- 0
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
  alpha.grid  = alpha.grid
)

# 5) Save the result to a file whose name reflects the mode
mode_tag <- if (generate.pcglasso) "pcglasso" else "glasso"
outname  <- paste0("results_", mode_tag, ".rds")
saveRDS(res, file = outname)
message("Saved results to ", outname, "\n")

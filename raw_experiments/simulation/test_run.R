library(space)
library(pcglassoFast)
library(PCGLASSO)
library(glasso)
library(parallel)
library(Matrix)
library(dplyr)
library(tidyverse)

source("estimation_methods.R")
source("simulation_functions.R")
# source("./raw_experiments/simulation/estimation_methods.R")
# source("./raw_experiments/simulation/simulation_functions.R")


# 2) Common settings (you can adjust these as needed)
set.seed(1)
graphics.off()
generate.pcglasso=T
split.train        <- 0.7      # not used directly below, but could be referenced in your functions
ns                 <- c(200)
sim                <- 7#21
nlambda            <- 50
mc_cores           <- 7
alpha_grid         <- 0
#lambda.min.ratio   <- 0.1
pcglasso_tolerance <- 0.01

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
res0_8 <- run_experiments(
  Q = Q,
  ns          = ns,
  sim         = sim,
  mc_cores    = mc_cores,
  nlambda     = nlambda,
  lambda.min.ratio = 0.8,
  alpha_grid  = alpha_grid,
  pcglasso_tolerance = pcglasso_tolerance
)

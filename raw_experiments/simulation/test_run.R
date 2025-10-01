source("estimation_methods.R")
library(space)
library(pcglassoFast)
source("simulation_functions.R")  # make sure this is in the same folder or adjust the path
library(glasso)
library(parallel)
library(Matrix)


# 2) Common settings (you can adjust these as needed)
set.seed(2)
graphics.off()
generate.pcglasso=T
split.train        <- 0.7      # not used directly below, but could be referenced in your functions
ns                 <- c(200)
sim                <- 4
nlambda            <- 50
mc_cores           <- 4
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

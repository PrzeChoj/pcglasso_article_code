# run_simulation.R
#
# Main driver for the Section 4.3 / Appendix B.2 simulation study.
# Runs the precision matrix estimation benchmark on two graph structures:
#
#   results_pcglasso.rds  -- hub graph   (Q_simulated_pcglasso)  → Section 4.3, Tables 1–2
#   results_glasso.rds    -- non-hub graph (Q_simulated_glasso)  → Appendix B.2, Tables 3–4
#
# Results are written to results/ and read later by make_tables_figures.R.
#
# Runtime note: with the default settings below (sim = 50, all five sample
# sizes including n = 5000), this script takes several hours on a standard
# laptop.  For a quick local test, reduce sim and/or drop n = 5000 from ns.

suppressPackageStartupMessages({
  library(space)
  library(pcglassoFast)   # provides pcglassoPath(), evaluate_objective_path(),
                          # compare_matrices(), cov2cor_inv(), and the Q datasets
  library(PCGLASSO)       # Carter algorithm: pcglasso()
  library(glasso)
  library(Matrix)
  library(MASS)
  library(parallel)
  library(ggplot2)
  library(patchwork)
})

source("./experiments/Section_4_3/simulation_functions.R")
source("./experiments/Section_4_3/estimation_methods.R")


# ---- Parameters ----

ns               <- c(200, 500, 1000, 5000)
sim              <- 200            # replications per sample size
nlambda          <- 50             # size of lambda grid (passed to run_experiments via ...)
lambda.min.ratio <- 0.01           # (passed to run_experiments via ...)
alpha_grid       <- 0              # PC-GLasso penalty on diagonal: 0 means standard GLasso penalty
                                   # (function default is a 10-point grid; override here to match paper)
mc_cores         <- parallel::detectCores()
pcglasso_tol     <- 0.001          # convergence tolerance (= 1e-3)

# Default estimator set used in the paper (passed to run_experiments via ...)
estimators_paper <- list(
  PCGLcpp_I  = estimator_pcglasso_I_primal,
  PCGLFor_I  = estimator_pcglasso_I_dual,
  GL         = estimator_glasso,
  SPACE      = estimator_space,
  CGL        = estimator_corglasso
)


# ---- Create output directory ------------------------------------------------

dir.create("./experiments/Section_4_3/results", showWarnings = FALSE)


start_time <- Sys.time()

# ---- Hub graph simulation  (Section 4.3, Tables 1–2) ------------------------

message("=== Hub graph (Q_simulated_pcglasso) ===")
data(Q_simulated_pcglasso, package = "pcglassoFast")
Q_hub <- (Q_simulated_pcglasso + t(Q_simulated_pcglasso)) / 2

res_pcglasso <- run_experiments(
  Q                  = Q_hub,
  ns                 = ns,
  sim                = sim,
  mc_cores           = mc_cores,
  seed               = 2,
  estimators         = estimators_paper,
  pcglasso_tolerance = pcglasso_tol,
  nlambda            = nlambda,
  lambda.min.ratio   = lambda.min.ratio,
  alpha_grid         = alpha_grid
)

saveRDS(res_pcglasso, file = "./experiments/Section_4_3/results/results_pcglasso.rds")
message("Saved results/results_pcglasso.rds\n")


# ---- Non-hub graph simulation  (Appendix B.2, Tables 3–4) ------------------

message("=== Non-hub graph (Q_simulated_glasso) ===")
data(Q_simulated_glasso, package = "pcglassoFast")
Q_nohub <- (Q_simulated_glasso + t(Q_simulated_glasso)) / 2

res_glasso <- run_experiments(
  Q                  = Q_nohub,
  ns                 = ns,
  sim                = sim,
  mc_cores           = mc_cores,
  seed               = 2,
  estimators         = estimators_paper,
  pcglasso_tolerance = pcglasso_tol,
  nlambda            = nlambda,
  lambda.min.ratio   = lambda.min.ratio,
  alpha_grid         = alpha_grid
)

saveRDS(res_glasso, file = "./experiments/Section_4_3/results/results_glasso.rds")
message("Saved results/results_glasso.rds\n")

end_time <- Sys.time()
print(end_time - start_time)

message("Done. Run make_tables_figures.R to produce tables and plots.")

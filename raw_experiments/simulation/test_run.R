library(space)
library(pcglassoFast)
library(glasso)
library(parallel)
library(Matrix)
library(dplyr)

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

res0_1 <- run_experiments(
  Q = Q,
  ns          = ns,
  sim         = sim,
  mc_cores    = mc_cores,
  nlambda     = nlambda,
  lambda.min.ratio = 0.1,
  alpha_grid  = alpha_grid,
  pcglasso_tolerance = pcglasso_tolerance
)

#####
res0_1_summarized <- res0_1 %>%
  select(method, timing, rmse) %>%
  group_by(method) %>%
  summarise(timing = mean(timing), rmse = mean(rmse))
res0_8_summarized <- res0_8 %>%
  select(method, timing, rmse) %>%
  group_by(method) %>%
  summarise(timing = mean(timing), rmse = mean(rmse))

res0_1_summarized
# # A tibble: 12 × 3
#   method       timing  rmse
#   <chr>         <dbl> <dbl>
# 1 CGL_bic       0.926 1.50
# 2 CGL_cv        1.06  1.40
# 3 GL_bic        0.233 1.63
# 4 GL_cv         0.264 1.49
# 5 PCGL_C_bic   39.1   0.621
# 6 PCGL_C_cv    46.1   1.57
# 7 PCGL_I_bic   44.3   0.562
# 8 PCGL_I_cv    46.9   1.48
# 9 PCGLcpp_bic   8.77  0.441
# 10 PCGLcpp_cv   8.31  0.976
# 11 SPACE_bic    6.19  1.06
# 12 SPACE_cv     3.86  0.687

res0_8_summarized
# A tibble: 12 × 3
#   method       timing  rmse
#   <chr>         <dbl> <dbl>
# 1 CGL_bic       0.259 1.56
# 2 CGL_cv        0.263 1.56
# 3 GL_bic        0.160 1.64
# 4 GL_cv         0.166 1.64
# 5 PCGL_C_bic   31.5   0.886
# 6 PCGL_C_cv    26.4   1.00
# 7 PCGL_I_bic   39.3   0.775
# 8 PCGL_I_cv    33.0   1.17
# 9 PCGLcpp_bic   3.46  0.479
# 10 PCGLcpp_cv   3.42  0.476
# 11 SPACE_bic    5.93  1.06
# 12 SPACE_cv     3.58  0.661

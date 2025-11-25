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
lambda.min.ratio   <- 0.1

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
  pcglasso_tolerance = 0.01
)

res0_1 <- run_experiments(
  Q = Q,
  ns          = ns,
  sim         = sim,
  mc_cores    = mc_cores,
  nlambda     = nlambda,
  lambda.min.ratio = 0.1,
  alpha_grid  = alpha_grid,
  pcglasso_tolerance = 0.01
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
# A tibble: 12 × 3
#   method        timing  rmse
#   <chr>          <dbl> <dbl>
# 1 CGL_bic        0.893 1.49
# 2 CGL_cv         1.06  1.38
# 3 GL_bic         0.220 1.63
# 4 GL_cv          0.247 1.47
# 5 PCGL_bic      31.8   0.621
# 6 PCGL_cv       31.9   1.12
# 7 PCGL_old_bic  30.2   0.514
# 8 PCGL_old_cv   28.5   1.09
# 9 PCGLcpp_bic    7.95  0.346
# 10 PCGLcpp_cv    7.81  1.00
# 11 SPACE_bic     5.64  1.01
# 12 SPACE_cv      3.40  0.562

res0_8_summarized
# A tibble: 12 × 3
#   method        timing  rmse
#   <chr>          <dbl> <dbl>
# 1 CGL_bic        0.234 1.55
# 2 CGL_cv         0.236 1.55
# 3 GL_bic         0.142 1.63
# 4 GL_cv          0.140 1.63
# 5 PCGL_bic      27.3   0.659
# 6 PCGL_cv       24.9   0.850
# 7 PCGL_old_bic  20.4   0.600
# 8 PCGL_old_cv   16.6   0.775
# 9 PCGLcpp_bic    2.98  0.359
# 10 PCGLcpp_cv    3.02  0.333
# 11 SPACE_bic     4.80  1.02
# 12 SPACE_cv      2.97  0.584

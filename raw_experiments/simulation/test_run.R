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
set.seed(2)
graphics.off()
generate.pcglasso=T
split.train        <- 0.7      # not used directly below, but could be referenced in your functions
ns                 <- c(200)
sim                <- 21
nlambda            <- 50
mc_cores           <- 7
#mc_cores <- 1L
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
res1 <- run_experiments(
  Q = Q,
  ns          = ns,
  sim         = sim,
  mc_cores    = mc_cores,
  nlambda     = nlambda,
  lambda.min.ratio = lambda.min.ratio,
  alpha_grid  = alpha_grid,
  pcglasso_tolerance = 0.1
)
res2 <- run_experiments(
  Q = Q,
  ns          = ns,
  sim         = sim,
  mc_cores    = mc_cores,
  nlambda     = nlambda,
  lambda.min.ratio = lambda.min.ratio,
  alpha_grid  = alpha_grid,
  pcglasso_tolerance = 0.01
)
res3 <- run_experiments(
  Q = Q,
  ns          = ns,
  sim         = sim,
  mc_cores    = mc_cores,
  nlambda     = nlambda,
  lambda.min.ratio = lambda.min.ratio,
  alpha_grid  = alpha_grid,
  pcglasso_tolerance = 0.001
)

#####
res1_summarized <- res1 %>%
  select(method, timing, rmse) %>%
  group_by(method) %>%
  summarise(timing = mean(timing), rmse = mean(rmse))
res2_summarized <- res2 %>%
  select(method, timing, rmse) %>%
  group_by(method) %>%
  summarise(timing = mean(timing), rmse = mean(rmse))
res3_summarized <- res3 %>%
  select(method, timing, rmse) %>%
  group_by(method) %>%
  summarise(timing = mean(timing), rmse = mean(rmse))

#####

#res1 # tolerance 0.1
#res2 # tolerance 0.01
#res3 # tolerance 0.001

times <- rbind(
  res1_summarized$timing, res2_summarized$timing, res3_summarized$timing
)
colnames(times) <- res1_summarized$method

rmse <- rbind(
  res1_summarized$rmse, res2_summarized$rmse, res3_summarized$rmse
)
colnames(rmse) <- res1_summarized$method

times
#        CGL_bic   CGL_cv    GL_bic     GL_cv  PCGL_bic   PCGL_cv SPACE_bic SPACE_cv
# [1,] 0.9066667 1.045619 0.2360000 0.2700476  5.348381  8.891429  6.486762 3.875238
# [2,] 1.1954762 1.377857 0.3181429 0.3500476 44.157048 52.304667  8.404048 5.089714
# [3,] 1.1631905 1.341048 0.3062381 0.3455238 55.524048 71.630000  8.251000 5.004381
rmse
#       CGL_bic   CGL_cv   GL_bic    GL_cv  PCGL_bic  PCGL_cv SPACE_bic  SPACE_cv
# [1,] 1.500033 1.387332 1.627034 1.489864 1.1434563 3.164260  1.035873 0.5844242
# [2,] 1.498104 1.393602 1.627513 1.489561 0.7037663 1.437055  1.034759 0.6195062
# [3,] 1.499412 1.397609 1.628953 1.494241 0.5013122 0.748949  1.034365 0.6415307

# We can sey one can use big tolerance = 0.1 for approximately the same time as SPACE and similar results
# but one can set better tolerance = 0.001 for 10 times more time, but very good results

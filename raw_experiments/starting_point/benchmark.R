library(pcglassoFast)
library(foreach)
library(doSNOW)
library(snow)

source("./raw_experiments/starting_point/get_S.R")
source("./raw_experiments/starting_point/paths_functions.R")



set.seed(42)
#p_values <- c(10, 30, 50, 100)
p_values <- c(10, 20, 30, 40, 50)
alpha_values <- c(-0.1, 0, 0.2)
experiment_values <- 1:4
replications <- 30

param_grid <- expand.grid(
  p = p_values,
  alpha = alpha_values,
  experiment = experiment_values,
  rep_iter = 1:replications,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
num_jobs <- nrow(param_grid)
message(paste("Total number of jobs to run:", num_jobs))


num_cores <- max(1, parallel::detectCores() - 1)
message(paste("Setting up cluster with", num_cores, "cores."))
cl <- makeCluster(num_cores, type = "SOCK")
registerDoSNOW(cl)

pb <- txtProgressBar(min = 0, max = num_jobs, style = 3, width = 50, char = "=")
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# --- Main parallel loop to collect timing data ---
message("Starting benchmark...")
results_df <- foreach(
  i = 1:num_jobs,
  .combine = rbind,
  .options.snow = opts
) %dopar% {
  p_val <- param_grid$p[i]
  alpha_val <- param_grid$alpha[i]
  exp_val <- param_grid$experiment[i]

  S <- get_S(p = p_val, which_experiment = exp_val)

  time_ud <- {start <- proc.time(); ud(S, alpha_val); proc.time()-start}[["elapsed"]]
  time_du <- {start <- proc.time(); du(S, alpha_val); proc.time()-start}[["elapsed"]]
  time_du2 <- {start <- proc.time(); du2(S, alpha_val); proc.time()-start}[["elapsed"]]
  time_du3 <- {start <- proc.time(); du3(S, alpha_val); proc.time()-start}[["elapsed"]]
  time_du4 <- {start <- proc.time(); du4(S, alpha_val); proc.time()-start}[["elapsed"]]
  time_udu <- {start <- proc.time(); udu(S, alpha_val); proc.time()-start}[["elapsed"]]
  time_dud <- {start <- proc.time(); dud(S, alpha_val); proc.time()-start}[["elapsed"]]

  data.frame(
    p = p_val,
    alpha = alpha_val,
    experiment = exp_val,
    replication = param_grid$rep_iter[i],
    time_ud,
    time_du,
    time_du2,
    time_du3,
    time_du4,
    time_udu,
    time_dud
  )
}
message("\nBenchmark finished.")

stopCluster(cl)
close(pb)


saveRDS(results_df, "./raw_experiments/starting_point/benchmark_results.rds")

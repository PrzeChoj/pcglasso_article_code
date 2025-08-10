library(pcglassoFast)
library(foreach)
library(doSNOW)
library(snow)

source("./raw_experiments/starting_point/get_S.R")
source("./raw_experiments/starting_point/starting_point_functions.R")

timed <- function(expr) {
  val <- NULL
  t <- system.time(val <- eval.parent(substitute(expr)))[["elapsed"]]
  list(time = as.numeric(t), res = as.numeric(val))
}

set.seed(42)
p_values <- c(10, 30, 50, 100)
alpha_values <- c(-0.1, 0, 0.2)
lambda_values <- c(0.01, 0.05, 0.1)
experiment_values <- 1:4
replications <- 30

param_grid <- expand.grid(
  p = p_values,
  alpha = alpha_values,
  lambda = lambda_values,
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
  i = seq_len(num_jobs),
  .combine = rbind,
  .options.snow = opts,
  .packages = c("pcglassoFast")
) %dopar% {
  p_val     <- param_grid$p[i]
  alpha_val <- param_grid$alpha[i]
  lambda_val<- param_grid$lambda[i]
  exp_val   <- param_grid$experiment[i]

  S <- get_S(p = p_val, which_experiment = exp_val)

  tud  <- timed(path_up_down(S, alpha_val, lambda_val));    time_path_ud      <- tud$time;  res_path_ud      <- tud$res
  tudu <- timed(path_up_down_up(S, alpha_val, lambda_val)); time_path_udu     <- tudu$time; res_path_udu     <- tudu$res
  ti   <- timed(start_I(S, alpha_val, lambda_val));         time_start_I      <- ti$time;   res_start_I      <- ti$res
  tc   <- timed(start_cor(S, alpha_val, lambda_val));       time_start_cor    <- tc$time;   res_start_cor    <- tc$res
  tg   <- timed(start_glasso(S, alpha_val, lambda_val));    time_start_glasso <- tg$time;   res_start_glasso <- tg$res
  tL2  <- timed(start_L2(S, alpha_val, lambda_val));        time_start_L2     <- tL2$time;  res_start_L2     <- tL2$res

  data.frame(
    p = p_val,
    alpha = alpha_val,
    lambda = lambda_val,
    experiment = exp_val,
    replication = param_grid$rep_iter[i],
    time_path_ud,
    time_path_udu,
    time_start_I,
    time_start_cor,
    time_start_glasso,
    time_start_L2,
    res_path_ud,
    res_path_udu,
    res_start_I,
    res_start_cor,
    res_start_glasso,
    res_start_L2,
    stringsAsFactors = FALSE
  )
}
message("\nBenchmark finished.")

stopCluster(cl)
close(pb)

# saveRDS(results_df, "./raw_experiments/starting_point/benchmark_results.rds")

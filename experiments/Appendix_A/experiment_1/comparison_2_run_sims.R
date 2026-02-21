library(parallel)
library(pbmcapply)

source("./experiments/Appendix_A/experiment_1/comparison_1_utils.R")

data_dir <- "./experiments/Appendix_A/experiment_1/res_data"
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

# hyper-parameters
M <- 50
p_vec                 <- c(10, 50, 100, 150)
cor_modifier_vec_hub  <- c(1, 0.9)
cor_modifier_vec_line <- c(0.8, 0.9)
lambda_vec            <- c(0.1, 0.2)
alpha_vec             <- c(0, 0.5)

n_cores <- max(1, min(7, detectCores(logical = FALSE) - 1))

tolerance_list <- exp(seq(log(0.01), log(0.00000001), length.out = 12))

set.seed(123)
Sys.setenv(OMP_NUM_THREADS = 1)


#####
param_grid_hub <- expand.grid(
  p            = p_vec,
  cor_modifier = cor_modifier_vec_hub,
  lambda       = lambda_vec,
  alpha        = alpha_vec,
  K_structure  = "hub",
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

param_grid_line <- expand.grid(
  p            = p_vec,
  cor_modifier = cor_modifier_vec_line,
  lambda       = lambda_vec,
  alpha        = alpha_vec,
  K_structure  = "line",
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

param_grid <- rbind(param_grid_hub, param_grid_line)

run_one <- function(row) {
  p            <- row$p
  cor_modifier <- row$cor_modifier
  lambda       <- row$lambda
  alpha        <- row$alpha
  K_structure  <- row$K_structure

  tryCatch({
    cat("Running M =", M,
        " p =", p,
        " cor =", cor_modifier,
        " lambda =", lambda,
        " alpha =", alpha,
        " K_structure =", K_structure, "\n")

    df <- simulate_pcglasso(
      M              = M,
      p              = p,
      cor_modifier   = cor_modifier,
      lambda         = lambda,
      alpha          = alpha,
      tolerance_list = tolerance_list,
      K_structure    = K_structure
    )

    cor_str    <- gsub("\\.", "_", as.character(cor_modifier))
    lambda_str <- gsub("\\.", "_", as.character(lambda))
    alpha_str  <- gsub("\\.", "_", as.character(alpha))
    K_str      <- as.character(K_structure)

    file_name <- sprintf(
      "%s/comparison_%s_M%d_p%d_cor%s_lambda%s_alpha%s.rds",
      data_dir, K_str, M, p, cor_str, lambda_str, alpha_str
    )

    saveRDS(df, file_name)

    list(
      p = p,
      cor = cor_modifier,
      lambda = lambda,
      alpha = alpha,
      K_structure = K_structure,
      status = "ok"
    )
  }, error = function(e) {
    message("FAILED: p=", p,
            " cor=", cor_modifier,
            " lambda=", lambda,
            " alpha=", alpha,
            " K_structure=", K_structure,
            " | ", conditionMessage(e))
    list(
      p = p,
      cor = cor_modifier,
      lambda = lambda,
      alpha = alpha,
      K_structure = K_structure,
      status = "fail"
    )
  })
}

start_time <- Sys.time()
status <- pbmclapply(
  seq_len(nrow(param_grid)),
  function(i) run_one(param_grid[i, ]),
  mc.cores = n_cores,
  mc.set.seed = TRUE
)
end_time <- Sys.time()

statuses <- vapply(status, `[[`, character(1), "status")

if (all(statuses == "ok")) {
  message("Simulation Successful")
  print(end_time - start_time)
} else {
  message("Simulation FAILED")
}

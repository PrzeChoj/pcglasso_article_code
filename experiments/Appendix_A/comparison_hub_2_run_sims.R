library(parallel)
library(pbmcapply)

source("./experiments/Appendix_A/comparison_hub_1_utils.R")

# hyper-parameters
M <- 50
p_vec            <- c(50, 70)
cor_modifier_vec <- c(1, 0.9)
lambda_vec       <- c(0.1, 0.2)
alpha_vec        <- c(0, 0.5)

n_cores <- min(7, detectCores(logical = FALSE) - 1)

tolerance_list <- exp(seq(log(0.01), log(0.0001), length.out = 12))

set.seed(123)
Sys.setenv(OMP_NUM_THREADS = 1)



#####
start_time <- Sys.time()
param_grid <- expand.grid(
  p            = p_vec,
  cor_modifier = cor_modifier_vec,
  lambda       = lambda_vec,
  alpha        = alpha_vec,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

run_one <- function(row) {
  p            <- row$p
  cor_modifier <- row$cor_modifier
  lambda       <- row$lambda
  alpha        <- row$alpha

  tryCatch({
    cat("Running M =", M,
        " p =", p,
        " cor =", cor_modifier,
        " lambda =", lambda,
        " alpha =", alpha, "\n")

    df <- simulate_pcglasso(
      M              = M,
      p              = p,
      cor_modifier   = cor_modifier,
      lambda         = lambda,
      alpha          = alpha,
      tolerance_list = tolerance_list
    )

    cor_str    <- gsub("\\.", "_", as.character(cor_modifier))
    lambda_str <- gsub("\\.", "_", as.character(lambda))
    alpha_str  <- gsub("\\.", "_", as.character(alpha))

    file_name <- sprintf(
      "./experiments/Appendix_A/res_data/comparison_hub_M%d_p%d_cor%s_lambda%s_alpha%s.rds",
      M, p, cor_str, lambda_str, alpha_str
    )

    saveRDS(df, file_name)

    list(
      p = p,
      cor = cor_modifier,
      lambda = lambda,
      alpha = alpha,
      status = "ok"
    )
  }, error = function(e) {
    message("FAILED: ", row$p, " ", row$cor_modifier,
            " ", row$lambda, " ", row$alpha)
    list(
      p = p,
      cor = cor_modifier,
      lambda = lambda,
      alpha = alpha,
      status = "fail"
    )
  })
}

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

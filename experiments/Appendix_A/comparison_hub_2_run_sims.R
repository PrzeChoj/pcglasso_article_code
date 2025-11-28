source("./experiments/Appendix_A/comparison_hub_1_utils.R")

# hyper-parameters
M <- 50
p_vec            <- c(50, 70)
cor_modifier_vec <- c(1, 0.9)
lambda_vec       <- c(0.1, 0.2)
alpha_vec        <- c(0, 0.5)

tolerance_list <- exp(seq(log(0.01), log(0.0001), length.out = 12))

for (p in p_vec) {
  for (cor_modifier in cor_modifier_vec) {
    for (lambda in lambda_vec) {
      for (alpha in alpha_vec) {

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
      }
    }
  }
}

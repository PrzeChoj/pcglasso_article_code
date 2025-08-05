# Load required libraries
library(huge)
library(snow)
library(doSNOW)
library(foreach)
library(Matrix)
library(ggplot2)
library(dplyr)
library(tidyr)

# remotes::install_github("PrzeChoj/pcglassoFast", ref = "DOptimizationTime")
library(pcglassoFast)

# test whether proper pcglassoFast package version is installed
if (!("full_time_D_optim" %in% names(pcglassoFast(S = matrix(c(1.1, 0, 0, 1.1), nrow = 2), 0.1, 0.1)))) {
  stop('Special version of `pcglassoFast` package is needed for this script: `remotes::install_github("PrzeChoj/pcglassoFast", ref = "DOptimizationTime")`')
}


# Load stock market data from the 'huge' package
data("stockdata", package = "huge")
my_data <- stockdata$data # Rows: time points, Columns: companies
log_returns <- log(my_data[-1, ] / my_data[-nrow(my_data), ])
log_returns <- sweep(log_returns, 1, rowMeans(log_returns))

# Simulation parameters
set.seed(42)
sim <- 200                        # Number of replications per (p, lambda) combination
n   <- 400                        # Sample size (number of time points) for each replication
p_vec <- c(50, 100, 150, 300)     # Different numbers of companies to test
alpha_vec <- c(-1, 0, 0.9)

# Create a grid of all parameter combinations: (p, lambda, replication)
param_grid <- expand.grid(
  p = p_vec,
  alpha = alpha_vec,
  rep_iter = seq_len(sim),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
num_jobs <- nrow(param_grid)

# Set up a SNOW cluster with a progress bar
num_cores <- min(parallel::detectCores(), 7)
cl <- makeCluster(num_cores, type = "SOCK")
registerDoSNOW(cl)
pb <- txtProgressBar(min = 0, max = num_jobs, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# MAIN PARALLEL LOOP: iterate over all (p, alpha, replication) combinations
results_list <- foreach(
  i = seq_len(num_jobs),
  .combine = rbind,
  .options.snow = opts,
  .packages = c("PCGLASSO")
) %dopar% {
  # Extract current simulation parameters
  p <- param_grid$p[i]
  lambda <- 0.1
  alpha <- param_grid$alpha[i]
  rep_iter <- param_grid$rep_iter[i]

  # --- Data Processing using Stock Market Data ---
  # Randomly select p companies (columns)
  selected <- sample(ncol(my_data), p)
  log_returns_selected <- log_returns[, selected]

  # Remove rows with outliers: any value deviating more than 5 SDs from its column mean
  col_means <- colMeans(log_returns_selected)
  col_sds <- apply(log_returns_selected, 2, sd)
  outlier_matrix <- abs(sweep(log_returns_selected, 2, col_means)) > 5 * rep(col_sds, each = nrow(log_returns_selected))
  keep_rows <- !apply(outlier_matrix, 1, any)
  clean_data <- log_returns_selected[keep_rows, ]

  # Sample n observations (rows) from the cleaned data (with replacement if needed)
  if (nrow(clean_data) < n) {
    X <- clean_data[sample(nrow(clean_data), n, replace = TRUE), ]
  } else {
    X <- clean_data[sample(nrow(clean_data), n), ]
  }

  # Compute the empirical covariance matrix
  S <- cov(X) + diag(0.0001, p)

  # --- Time the Two Methods ---
  res_exact_newton <- pcglassoFast::pcglassoFast(
    S,
    lambda,
    alpha = alpha,
    tolerance = 1e-4,
    max_iter = 1000,
    diagonal_Newton = FALSE
  )
  t_exact_newton <- res_exact_newton[["full_time_D_optim"]]

  res_diagonal_newton <- pcglassoFast::pcglassoFast(
    S,
    lambda,
    alpha = alpha,
    tolerance = 1e-4,
    max_iter = 1000,
    diagonal_Newton = TRUE
  )
  t_diagonal_newton <- res_diagonal_newton[["full_time_D_optim"]]

  # Compute the number of nonzero off-diagonal elements (edges)
  exact_newton_nonzero <- (sum(res_exact_newton$Sinv != 0) - p) / 2
  diagonal_newton_nonzero <- (sum(res_diagonal_newton$Sinv != 0) - p) / 2

  # Extract final loss values (if available)
  loss_exact_newton <- if (is.null(res_exact_newton$loss)) NA else tail(res_exact_newton$loss, 1)
  loss_diagonal_newton <- if (is.null(res_diagonal_newton$loss)) NA else tail(res_diagonal_newton$loss, 1)

  # Return simulation results for this iteration as a one-row data frame
  data.frame(
    p = p,
    alpha = alpha,
    replication = rep_iter,
    time_exact_newton = t_exact_newton,
    time_diagonal_newton = t_diagonal_newton,
    exact_newton_nonzero = exact_newton_nonzero,
    diagonal_newton_nonzero = diagonal_newton_nonzero,
    loss_exact_newton = loss_exact_newton,
    loss_diagonal_newton = loss_diagonal_newton
  )
}

stopCluster(cl)
close(pb)

# Reshape results to long format for plotting run times
results_long <- results_list %>%
  pivot_longer(
    cols = c("time_exact_newton", "time_diagonal_newton"),
    names_to = "Method",
    values_to = "elapsed_time"
  ) %>%
  mutate(Method = ifelse(
    Method == "time_exact_newton",
    "Exact Newton",
    "Diagonal Newton"
  ))

# Summarize average run times and confidence intervals by (p, lambda, method)
plot_data <- results_long %>%
  group_by(p, alpha, Method) %>%
  summarise(
    mean_time = mean(elapsed_time),
    sd_time = sd(elapsed_time),
    n = n(),
    se_time = sd_time / sqrt(n),
    lower_ci = mean_time - 1.96 * se_time,
    upper_ci = mean_time + 1.96 * se_time,
    .groups = "drop"
  ) %>%
  mutate(alpha = factor(
    alpha,
    levels = unique(alpha),
    labels = paste0("alpha == ", unique(alpha))
  ))

# Plot mean run times vs. dimension (p) for each method, faceted by lambda value
fig <- ggplot(plot_data, aes(x = p, y = mean_time, color = Method)) +
  geom_ribbon(
    aes(ymin = lower_ci, ymax = upper_ci, fill = Method),
    alpha = 0.4, colour = NA, show.legend = FALSE
  ) +
  geom_line(aes(linetype = Method), linewidth = 0.4) +
  geom_point(size = 0.4) +
  facet_wrap(~alpha, scales = "free_y", labeller = label_parsed) +
  labs(
    x = "Dimension (p)",
    y = "Mean Time (seconds) of D optimization",
    color = "Method"
  ) +
  scale_y_log10(
    limits = c(0.0003, 0.3),
    breaks = c(0.001, 0.01, 0.1, 1),
    labels = c("0.001", "0.01", "0.1", "1")
  ) +
  scale_x_continuous(breaks = p_vec) +
  theme_minimal()

print(fig)

# ggsave(
#   "raw_experiments/stockmarket/diag_newton.png",
#   plot = fig, width = 7, height = 4
# )

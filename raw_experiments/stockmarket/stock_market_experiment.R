# Load required libraries
library(huge)
library(snow)
library(doSNOW)
library(foreach)
library(Matrix)
library(PCGLASSO) # devtools::install_github('JackStorrorCarter/PCGLASSO')
library(ggplot2)
library(dplyr)
library(tidyr)
library(pcglassoFast)

# Load stock market data from the 'huge' package
data("stockdata", package = "huge")
my_data <- stockdata$data # Rows: time points, Columns: companies
log_returns <- log(my_data[-1, ] / my_data[-nrow(my_data), ])
log_returns <- sweep(log_returns, 1, rowMeans(log_returns))

# Simulation parameters
set.seed(42)
sim <- 100                        # Number of replications per (p, lambda) combination
n   <- 400                        # Sample size (number of time points) for each replication
p_vec <- c(10, 50, 100, 150)      # Different numbers of companies to test
lambda_vec <- c(0.001, 0.01, 0.1) # Regularization parameter values

# Create a grid of all parameter combinations: (p, lambda, replication)
param_grid <- expand.grid(
  p = p_vec,
  lambda = lambda_vec,
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

# MAIN PARALLEL LOOP: iterate over all (p, lambda, replication) combinations
results_list <- foreach(
  i = seq_len(num_jobs),
  .combine = rbind,
  .options.snow = opts,
  .packages = c("PCGLASSO")
) %dopar% {
  # Extract current simulation parameters
  p <- param_grid$p[i]
  lambda <- param_grid$lambda[i]
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
  # Time PCGLASSO:::pcglasso.carter()
  start_carter <- proc.time()
  res_carter <- PCGLASSO::pcglasso(
      as.matrix(S),
      lambda,
      c            = 1,
      Theta_start = diag(p),
      threshold    = 1e-4,
      max_iter     = 1000
    )
  t_carter <- (proc.time() - start_carter)[["elapsed"]]
  res_carter <- list(
    Sinv = res_carter,
    loss = pcglassoFast:::function_to_optimize(cov2cor(res_carter), sqrt(diag(res_carter)), S, lambda, 0)
  )

  # Time pcglassoFast::pcglassoFast()
  start_fortran <- proc.time()
  res_fortran <- pcglassoFast::pcglassoFast(
    S,
    lambda,
    alpha = 0,
    tolerance = 1e-4,
    max.iter = 1000
  )
  t_fortran <- (proc.time() - start_fortran)[["elapsed"]]

  # Compute the number of nonzero off-diagonal elements (edges)
  carter_nonzero <- (sum(res_carter$Sinv != 0) - p) / 2
  fortran_nonzero <- (sum(res_fortran$Sinv != 0) - p) / 2

  # Extract final loss values (if available)
  loss_carter <- if (is.null(res_carter$loss)) NA else tail(res_carter$loss, 1)
  loss_fortran <- if (is.null(res_fortran$loss)) NA else tail(res_fortran$loss, 1)

  # Return simulation results for this iteration as a one-row data frame
  data.frame(
    p = p,
    lambda = lambda,
    replication = rep_iter,
    time_carter = t_carter,
    time_Fortran = t_fortran,
    carter_nonzero = carter_nonzero,
    fortran_nonzero = fortran_nonzero,
    loss_carter = loss_carter,
    loss_fortran = loss_fortran
  )
}

stopCluster(cl)
close(pb)

# Reshape results to long format for plotting run times
results_long <- results_list %>%
  pivot_longer(
    cols = c("time_carter", "time_Fortran"),
    names_to = "method",
    values_to = "elapsed_time"
  )

# Summarize average run times and confidence intervals by (p, lambda, method)
plot_data <- results_long %>%
  group_by(p, lambda, method) %>%
  summarise(
    mean_time = mean(elapsed_time),
    sd_time = sd(elapsed_time),
    n = n(),
    se_time = sd_time / sqrt(n),
    lower_ci = mean_time - 1.96 * se_time,
    upper_ci = mean_time + 1.96 * se_time,
    .groups = "drop"
  ) %>%
  mutate(lambda = factor(lambda,
                         levels = unique(lambda),
                         labels = paste0("lambda == ", unique(lambda))
  ))

# Plot mean run times vs. dimension (p) for each method, faceted by lambda value
fig <- ggplot(plot_data, aes(x = p, y = mean_time, color = method)) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = method),
              alpha = 0.4, colour = NA, show.legend = FALSE
  ) +
  geom_line(aes(linetype = method), linewidth = 0.4) +
  geom_point(size = 0.4) +
  facet_wrap(~lambda, scales = "free_y", labeller = label_parsed) +
  labs(
    x = "Dimension (p)",
    y = "Mean Time (seconds)",
    color = "Method"
  ) +
  scale_y_continuous(limits = c(0, 5.5), expand = c(0,0)) +
  theme_minimal()

print(fig)
ggsave(
  "./raw_experiments/simulation_stock_market.pdf", fig,
  width = 7.5, height = 4.7
)

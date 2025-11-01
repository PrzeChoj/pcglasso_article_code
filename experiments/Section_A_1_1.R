#devtools::install_github('JackStorrorCarter/PCGLASSO')
#devtools::install_github("PrzeChoj/pcglassoFast")

library(PCGLASSO)
library(pcglassoFast)
library(MASS)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork) # or library(gridExtra)

source("./experiments/utils.R")

M <- 100
p_vector <- c(10, 50, 100, 150)
n <- 400
lambda_vector <- c(0.01, 0.05, 0.1)
alpha <- 0
tolerance <- 0.001

get_plot <- function(which_figure) {
  S_maker <- if(which_figure == 1) {
    function(p) {
      S_star <- diag(1, p); S_star[1,2:p] <- S_star[2:p,1] <- -1/sqrt(p); S_star[1,1] <- 1
      Z <- mvrnorm(n = n, mu = rep(0, p), Sigma = S_star)
      S <- t(Z) %*% Z

      S
    }
  } else if (which_figure == 2) {
    function(p) {
      stopifnot(abs(p / 5 - round(p / 5)) < 0.00001)
      K <- round(p/5)
      S_star <- diag(1, p);
      for (i in 1:5) {
        S_star[(K*(i-1)+2):(K*i),(K*(i-1)+1)] <- S_star[(K*(i-1)+1),(K*(i-1)+2):(K*i)] <- -1/sqrt(p)
      }

      Z <- mvrnorm(n = n, mu = rep(0, p), Sigma = S_star)
      S <- t(Z) %*% Z

      S
    }
  } else {
    stop()
  }

  # positive -> pcglasso is better; negative -> pcglassoFast is better
  res_val_diff <- matrix(nrow = length(p_vector), ncol = length(lambda_vector))
  res_time_pcglassoFast_mean <- matrix(nrow = length(p_vector), ncol = length(lambda_vector))
  res_time_pcglassoFast_sd <- matrix(nrow = length(p_vector), ncol = length(lambda_vector))
  res_time_pcglasso_mean <- matrix(nrow = length(p_vector), ncol = length(lambda_vector))
  res_time_pcglasso_sd <- matrix(nrow = length(p_vector), ncol = length(lambda_vector))
  counter <- 0
  pb <- txtProgressBar(min = counter, max = length(p_vector) * length(lambda_vector) * M, style = 3)
  for (p_i in 1:length(p_vector)) {
    p <- p_vector[p_i]
    for (lambda_i in 1:length(lambda_vector)) {
      lambda <- lambda_vector[lambda_i]

      res_val_diff_single <- numeric(M)
      res_time_pcglassoFast_single <- numeric(M)
      res_time_pcglasso_single <- numeric(M)
      for (m in 1:M) {
        counter <- counter + 1
        setTxtProgressBar(pb, counter)

        S <- S_maker(p)

        # pcglassoFast
        start_time <- Sys.time()
        res <- pcglassoFast(S, lambda = lambda, alpha = alpha, tolerance = tolerance, verbose = 0)
        end_time <- Sys.time()
        res_time_pcglassoFast_single[m] <- as.numeric(difftime(end_time, start_time, units = "secs"))

        Sinv <- res$Sinv
        res_val_pcglassoFast <- pcglasso_goal_function(S, lambda, alpha, cov2cor(Sinv), sqrt(diag(Sinv)))

        # pcglasso
        start_time <- Sys.time()
        res <- pcglasso(S, lambda, 1 - alpha, threshold = tolerance)
        end_time <- Sys.time()
        res_time_pcglasso_single[m] <- as.numeric(difftime(end_time, start_time, units = "secs"))

        Sinv <- res
        res_val_pcglasso <- pcglasso_goal_function(S, lambda, alpha, cov2cor(Sinv), sqrt(diag(Sinv)))

        # res val diff; positive -> pcglasso is better; negative -> pcglassoFastis better
        res_val_diff_single[m] <- res_val_pcglassoFast - res_val_pcglasso
      }
      res_val_diff[p_i, lambda_i] <- mean(res_val_diff_single)
      res_time_pcglassoFast_mean[p_i, lambda_i] <- mean(res_time_pcglassoFast_single)
      res_time_pcglassoFast_sd[p_i, lambda_i] <- sd(res_time_pcglassoFast_single)
      res_time_pcglasso_mean[p_i, lambda_i] <- mean(res_time_pcglasso_single)
      res_time_pcglasso_sd[p_i, lambda_i] <- sd(res_time_pcglasso_single)
    }
  }

  # plot
  df <- data.frame(
    p = rep(p_vector, times = 2 * length(lambda_vector)),
    lambda = rep(rep(lambda_vector, each = length(p_vector)), times = 2),
    Method = rep(c("pcglassoFast", "pcglasso"), each = length(p_vector) * length(lambda_vector)),
    Mean = c(as.vector(res_time_pcglassoFast_mean),
             as.vector(res_time_pcglasso_mean)),
    SD   = c(as.vector(res_time_pcglassoFast_sd),
             as.vector(res_time_pcglasso_sd))
  )
  ymax <- max(df$Mean + df$SD, na.rm = TRUE)

  figure <- ggplot(df, aes(x = p, y = Mean, color = Method)) +
    geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD, fill = Method),
                alpha = 0.25, colour = NA) +
    geom_line(aes(linetype = Method), linewidth = 1) +
    geom_point(size = 2) +
    facet_wrap(
      ~ lambda,
      labeller = labeller(lambda = function(x) paste0("\u03BB = ", x)),
      scales = "fixed"
    ) +
    scale_x_continuous(breaks = p_vector) +
    #scale_y_log10() +
    coord_cartesian(ylim = c(0, ymax)) +
    labs(x = "Dimension (p)", y = "Mean Time (seconds)") +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "none",
      strip.text = element_text(size = 14),
      panel.spacing.x = unit(0.8, "cm")
    ) +
    scale_linetype_manual(values = c(pcglassoFast = "dashed", pcglasso = "solid"))

  list(figure, res_val_diff)
}



# Hub
list_hub <- get_plot(1)
figure_Hub <- list_hub[[1]]
print(figure_Hub)
list_hub[[2]] # close numbers, the results of pcglasso and pcglassoFast are similar; positive -> pcglasso is better; negative -> pcglassoFastis better

# Block Hub
list_Blockhub <- get_plot(2)
figure_BlockHub <- list_Blockhub[[1]]
print(figure_BlockHub)
list_Blockhub[[2]] # close numbers, the results of pcglasso and pcglassoFast are similar; positive -> pcglasso is better; negative -> pcglassoFastis better

fig <- (figure_Hub | figure_BlockHub)
print(fig)
# ggsave("./outputs/Figure_5_simulation.png", plot = fig, width = 11, height = 4)

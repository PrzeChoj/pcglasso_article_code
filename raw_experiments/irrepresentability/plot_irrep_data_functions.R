library(MASS)
library(glasso)
library(pcglassoFast)
library(parallel)
library(ggplot2)
library(tidyr)
library(dplyr)

#' @examples
#' perform_job(Sigma_true, 80, "GLASSO", (1:10)*0.04)
perform_job <- function(Sigma_true, n, optim_mode, lambda_list) {
  p <- dim(Sigma_true)[1]
  X <- mvrnorm(n, Sigma = Sigma_true, mu = rep(0, p))
  Cmle <- t(X) %*% X / n # we know the true mean.

  if (optim_mode == "GLASSO") {
    ans <- glassopath(
      Cmle, rholist = lambda_list,
      penalize.diagonal = FALSE, trace = 0
    )

    list(
      Cmle = Cmle,
      lambdas = lambda_list,
      S_est = ans$wi
    )
  } else if (optim_mode == "PCGLASSO") {
    ans <- pcglassoPath(
      Cmle, alpha = 0, lambdas = lambda_list
    )

    list(
      Cmle = Cmle,
      lambdas = lambda_list,
      S_est = simplify2array(ans$W_path) # TODO: investigate
    )
  }
}

perform_jobs <- function(K_true, n_list, lambda_list, N, n_cores = 1) {
  Sigma_true <- solve(K_true)

  table_of_tasks <- data.frame(
    Ni = rep(1:N, times = 2 * length(n_list)),
    n = rep(n_list, times = 2, each = N),
    optim_mode = rep(c("GLASSO", "PCGLASSO"), each = length(n_list) * N)
  )

  my_f <- function(i) {
    n <- table_of_tasks$n[i]
    optim_mode <- table_of_tasks$optim_mode[i]

    perform_job(Sigma_true, n, optim_mode, lambda_list)
  }

  jobs_results <- if (n_cores == 1) {
    lapply(1:(dim(table_of_tasks)[1]), my_f)
  } else {
    mclapply(1:(dim(table_of_tasks)[1]), my_f, mc.cores = n_cores)
  }

  attr(jobs_results, "table_of_tasks") <- table_of_tasks
  attr(jobs_results, "K_true") <- K_true

  jobs_results
}

analyse_results <- function(jobs_results) {
  table_of_tasks <- attr(jobs_results, "table_of_tasks")
  K_true <- attr(jobs_results, "K_true")

  N <- max(table_of_tasks$Ni)
  n_list <- unique(table_of_tasks$n)

  table_prob <- data.frame(
    GLASSO = rep(0, length(n_list)),
    PCGLASSO = rep(0, length(n_list))
  )

  rownames(table_prob) <- as.character(n_list)

  for (i in 1:nrow(table_of_tasks)) {
    all_wi <- jobs_results[[i]]$S_est

    if (any_wi_succeeded(all_wi, K_true)) {
      this_n <- as.character(table_of_tasks[i, "n"])
      this_optim_mode <- table_of_tasks[i, "optim_mode"]

      table_prob[this_n, this_optim_mode] <- 1/N + table_prob[this_n, this_optim_mode]
    }
  }

  table_prob
}

any_wi_succeeded <- function(all_wi, K_true) {
  for (i in 1:dim(all_wi)[3]) {
    wi <- all_wi[, , i]

    if (all((K_true == 0) == (abs(wi) < 1e-16))){
      return(TRUE)
    }
  }

  return(FALSE)
}

my_plot <- function(jobs_results) {
  table_prob <- analyse_results(jobs_results)

  # Reshape the dataframe into long format
  table_prob_long <- table_prob %>%
    mutate(n = as.numeric(row.names(table_prob))) %>%
    pivot_longer(cols = c(GLASSO, PCGLASSO), names_to = "Method", values_to = "Probability")

  # Plot using ggplot2
  ggplot(table_prob_long, aes(x = n, y = Probability, color = Method)) +
    geom_line(linewidth = 1) +
    scale_x_log10() +
    labs(x = "n", y = "Probability", title = "Comparison of GLASSO and PCGLASSO") +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    scale_color_manual(values = c("blue", "red"))
}

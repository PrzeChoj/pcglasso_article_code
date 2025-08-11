library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

method_levels <- c("start_glasso", "start_cor", "start_I", "start_L2", "path_udu", "path_ud")

results_file <- "./raw_experiments/starting_point/benchmark_results.rds"
if (!file.exists(results_file)) {
  stop(paste(
    "Error: The data file was not found.",
    "Please ensure 'benchmark_results.rds' is in the working directory."
  ))
}
results_df <- readRDS(results_file)


results_long <- results_df %>%
  pivot_longer(
    cols = starts_with("res_"),
    names_to = "Method",
    values_to = "found_target",
    names_prefix = "res_"
  ) %>%
  select(!starts_with("time_"))

results_long$found_target <- round(results_long$found_target, 1)

approx_rank <- function(x, tol = 1e-3) {
  ord <- order(-x)
  ranks <- integer(length(x))

  current_rank <- 1
  ranks[ord[1]] <- current_rank

  for (i in 2:length(x)) {
    prev <- x[ord[i - 1]]
    curr <- x[ord[i]]
    if (abs(curr - prev) > tol) {
      current_rank <- current_rank + 1
    }
    ranks[ord[i]] <- current_rank
  }
  ranks
}

results_ranked <- results_long %>%
  group_by(p, alpha, lambda, experiment, replication) %>%
  mutate(rank = approx_rank(found_target)) %>%
  ungroup() %>%
  select(!found_target)

results_times_best <- results_ranked %>%
  group_by(p, alpha, lambda, experiment, Method) %>%
  mutate(times_best = sum(rank == 1)) %>%
  ungroup() %>%
  select(!c("replication", "rank")) %>%
  distinct() %>%
  mutate(
    experiment = factor(case_when(
      experiment == 1 ~ "Hub",
      experiment == 2 ~ "AR",
      experiment == 3 ~ "Sanger",
      experiment == 4 ~ "StockMarket",
      TRUE ~ as.character(experiment)
    ), levels = c("Hub", "AR", "Sanger", "StockMarket"))
  )

num_p <- length(unique(results_df$p))
num_alpha <- length(unique(results_df$alpha))
num_lambda <- length(unique(results_df$lambda))
num_experiments <- length(unique(results_df$experiment))
num_methods <- length(unique(results_long$Method))
stopifnot(
  nrow(results_times_best) == num_p * num_alpha * num_lambda * num_experiments * num_methods
)

plot_alpha <- function(a) {
  df <- results_times_best |>
    filter(alpha == a) |>
    mutate(times_best = 100 * times_best / max(times_best)) |>
    mutate(method = factor(Method, levels = method_levels),
           p = factor(p))

  ggplot(df, aes(p, times_best, fill = method)) +
    geom_col(position = position_dodge2(preserve = "single", padding = 0.2),
             width = 0.7) +
    facet_grid(experiment ~ lambda) +
    labs(title = paste("Times best per method — alpha =", a), x = "p", y = "times best (%)", fill = NULL) +
    theme_bw(11) + theme(legend.position = "bottom")
}

alphas <- sort(unique(results_times_best$alpha))

stopifnot(length(alphas) == 3)
plot_alpha(alphas[1])
plot_alpha(alphas[2])
plot_alpha(alphas[3])

# invisible(lapply(
#   seq_along(alphas), \(i)
#   ggsave(sprintf("./raw_experiments/starting_point/plots/times_best_%s.png",
#                  as.character(alphas[i])),
#          plot_alpha(alphas[i]), width = 4, height = 5)))

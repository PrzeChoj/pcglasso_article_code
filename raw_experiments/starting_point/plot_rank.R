# R SCRIPT FOR VISUALIZING BENCHMARK RESULTS

library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

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

approx_rank <- function(x, tol = 1e-1) {
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

max_ranked <- results_ranked$rank %>% max

results_times_best <- results_ranked %>%
  group_by(p, alpha, lambda, experiment, Method) %>%
  mutate(times_best = sum(rank == 1)) %>%
  ungroup() %>%
  select(!c("replication", "rank")) %>%
  distinct()

num_p <- length(unique(results_df$p))
num_alpha <- length(unique(results_df$alpha))
num_lambda <- length(unique(results_df$lambda))
num_experiments <- length(unique(results_df$experiment))
num_methods <- length(unique(results_long$Method))
stopifnot(
  nrow(results_times_best) == num_p * num_alpha * num_lambda * num_experiments * num_methods
)

plot_alpha <- function(a, df = results_times_best) {
  # consistent method order within an alpha (best overall on top)
  method_levels <- df %>%
    filter(alpha == a) %>%
    group_by(Method) %>%
    summarise(total = sum(times_best), .groups = "drop") %>%
    arrange(desc(total)) %>%
    pull(Method)

  ggplot(
    df %>% filter(alpha == a) %>% mutate(Method = factor(Method, method_levels)),
    aes(x = Method, y = times_best, fill = Method)
  ) +
    geom_col(width = 0.75) +
    coord_flip() +
    facet_grid(experiment ~ lambda, labeller = label_both) + # 4×3 = 12 panels
    labs(title = paste("Times best per method — alpha =", a),
         x = NULL, y = "times_best") +
    theme_bw() +
    theme(legend.position = "none",
          strip.background = element_rect(fill = "grey95"))
}

alphas <- sort(unique(results_times_best$alpha))

stopifnot(length(alphas) == 3)
plot_alpha(alphas[1])
plot_alpha(alphas[2])
plot_alpha(alphas[3])

invisible(lapply(
  seq_along(alphas), \(i)
  ggsave(sprintf("./raw_experiments/starting_point/plots/times_best_%s.pdf",
                 as.character(alphas[i])),
         plot_alpha(alphas[i]), width = 11, height = 8)))

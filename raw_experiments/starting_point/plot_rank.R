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

results_df <- results_df %>%
  filter(p != 10)

number_of_experiments <- results_df %>%
  group_by(p, alpha, lambda, experiment) %>%
  summarise(number_of_experiments = n(), .groups = "drop") %>%
  pull(number_of_experiments) %>%
  unique

stopifnot(length(number_of_experiments) == 1)

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

plot_alpha <- function(alpha_val, plot_label = TRUE) {
  df <- results_times_best |>
    dplyr::filter(alpha == alpha_val)

  if (!"lambda_lab" %in% names(df)) df$lambda_lab <- paste0("\u03BB = ", df$lambda)
  if (!"experiment_lab" %in% names(df)) df$experiment_lab <- df$experiment

  df <- df |>
    dplyr::mutate(
      Method = factor(Method, levels = method_levels),
      p_fac  = factor(paste0("p = ", p), levels = paste0("p = ", sort(unique(p)))),
      times_best = 100 * times_best / number_of_experiments
    )

  ggplot(df, aes(x = times_best, y = Method, fill = Method)) +
    geom_col(width = .8) +
    facet_grid(experiment_lab ~ lambda_lab + p_fac, switch = "y", scales = "fixed") +
    scale_x_continuous(
      limits = c(0, 100),
      breaks  = seq(0, 100, 25),
      labels  = scales::label_number(accuracy = 1, suffix = "%"),
      expand  = expansion(mult = c(0, .02))
    ) +
    scale_fill_manual(
      values = c("#F08A84","#C7A74A","#2AA054","#3BC7CF","#6B8FE7","#E86AD8"),
      breaks = method_levels
    ) +
    labs(
      title = paste0("Times best per method for \u03B1 = ", alpha_val),
      x = "times best (%)", y = NULL
    ) +
    theme_bw(base_size = 14) +
    theme(
      legend.position = if(plot_label){"bottom"}else{"none"},
      legend.title = element_blank(),
      strip.placement = "outside",
      plot.title = element_text(hjust = .5, face = "bold"),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
      axis.text.y = element_text(hjust = 0),
      panel.spacing.y = unit(0.1, "lines")
    )
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
#          plot_alpha(alphas[i], plot_label = (i != 1)), width = 8, height = if(i == 1){6} else {7})))

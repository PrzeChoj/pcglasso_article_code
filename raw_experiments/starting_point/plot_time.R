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
    cols = starts_with("time_"),
    names_to = "Method",
    values_to = "elapsed_time",
    names_prefix = "time_"
  )

plot_data <- results_long %>%
  group_by(p, alpha, lambda, experiment, Method) %>%
  summarise(
    mean_time = mean(elapsed_time, na.rm = TRUE),
    sd_time   = sd(elapsed_time, na.rm = TRUE),
    n         = n(),
    se_time   = sd_time / sqrt(n),
    lower_ci  = pmax(0.00001, mean_time - 1.96 * se_time),
    upper_ci  = mean_time + 1.96 * se_time,
    .groups   = "drop"
  ) %>%
  mutate(
    experiment_label = factor(case_when(
      experiment == 1 ~ "Hub",
      experiment == 2 ~ "AR",
      experiment == 3 ~ "Sanger",
      experiment == 4 ~ "StockMarket",
      TRUE ~ as.character(experiment)
    ), levels = c("Hub", "AR", "Sanger", "StockMarket")),
    alpha_label = factor(paste("alpha ==", alpha),
                         levels = paste("alpha ==", sort(unique(results_df$alpha)))),
    lambda_label = factor(paste("lambda ==", lambda),
                         levels = paste("lambda ==", sort(unique(results_df$lambda)))
    )
  )


plot_runtime <- function(a, data = plot_data) {
  df <- dplyr::filter(data, alpha == a) %>%
    mutate(Method = factor(Method, levels = method_levels))
  if (nrow(df) == 0) stop("No data for alpha = ", a)

  ggplot(df, aes(p, mean_time, color = Method, group = Method)) +
    geom_line() +
    geom_point(size = 1.6) +
    #geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.08, alpha = 0.6) +
    scale_y_log10(
      breaks = scales::log_breaks(base = 10),
      labels = function(x) {
        s <- formatC(x, format = "f", digits = 3)
        sub("\\.?0+$", "", s)
      }
    ) +
    scale_x_continuous(breaks = sort(unique(df$p))) +
    labs(
      x = "p (number of variables)",
      y = "Mean runtime [s] (log10 scale)",
      title = paste0("Benchmark runtimes ±95% CI — alpha = ", a),
      subtitle = "Faceted by experiment (rows) and lambda (columns)"
    ) +
    facet_grid(
      rows = vars(experiment_label),
      cols  = vars(lambda_label),
      labeller = label_parsed,
      scales = "free_y"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text = element_text(size = 10)
    ) +
    guides(color = guide_legend(nrow = 1, byrow = TRUE))
}
alphas <- sort(unique(plot_data$alpha))

stopifnot(length(alphas) == 3)
plot_runtime(alphas[1])
plot_runtime(alphas[2])
plot_runtime(alphas[3])

# invisible(lapply(
#   seq_along(alphas), \(i)
#   ggsave(sprintf("./raw_experiments/starting_point/plots/runtime_%s.pdf",
#                  as.character(alphas[i])),
#          plot_runtime(alphas[i]), width = 11, height = 8)))

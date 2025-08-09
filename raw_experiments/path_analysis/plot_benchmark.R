# R SCRIPT FOR VISUALIZING BENCHMARK RESULTS

library(ggplot2)
library(dplyr)
library(tidyr)

results_file <- "./raw_experiments/starting_point/benchmark_results.rds"
if (!file.exists(results_file)) {
  stop(paste(
    "Error: The data file was not found.",
    "Please ensure 'benchmark_results.rds' is in the working directory."
  ))
}
results_df <- readRDS(results_file)


results_long <- results_df %>%
  mutate(time_du_and_ud = time_du + time_ud) %>%
  pivot_longer(
    cols = starts_with("time_"),
    names_to = "Method",
    values_to = "elapsed_time",
    names_prefix = "time_"
  ) %>%
  mutate(Method = recode(Method, du_and_ud = "du and ud"))

plot_data <- results_long %>%
  group_by(p, alpha, experiment, Method) %>%
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
      levels = paste("alpha ==", sort(unique(results_df$alpha)))
    )
  )

create_benchmark_plot <- function(data, title) {
  ggplot(data, aes(x = p, y = mean_time, color = Method, group = Method)) +
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = Method),
      alpha = 0.2, linetype = 0, show.legend = FALSE
    ) +
    geom_line(aes(linetype = Method), linewidth = 0.6) +
    geom_point(size = 1.5) +
    facet_grid(experiment_label ~ alpha_label, labeller = label_parsed, scales = "free_y") +
    scale_y_log10(
      breaks = c(0.0001, 0.001, 0.01, 0.1, 0.3, 1, 3, 10, 30),
      limits = function(cur_lims) {c(0.06, cur_lims[2])}
    ) +
    labs(
      title    = title,
      subtitle = "Mean execution time across different matrix types and parameters",
      x        = "Matrix Dimension (p)",
      y        = "Mean Time (seconds, log scale)",
      color    = "Function",
      linetype = "Function"
    ) +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 10, face = "bold"),
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 10),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 9)
    )
}

# Plot 1:
plot1_data <- plot_data %>% filter(Method %in% c("ud", "du", "du2", "du3", "du4"))
plot1 <- create_benchmark_plot(plot1_data, "Performance of 'ud' and 'du' variants")

# Plot 2:
plot2_data <- plot_data %>% filter(Method %in% c("udu", "dud", "du and ud"))
plot2 <- create_benchmark_plot(plot2_data, "Performance of Combined Functions")

print(plot1)
print(plot2)

# ggsave("./raw_experiments/starting_point/path_plot1.png", plot1, width = 5, height = 5, dpi = 300)
# ggsave("./raw_experiments/starting_point/path_plot2.png",  plot2, width = 5, height = 5, dpi = 300)

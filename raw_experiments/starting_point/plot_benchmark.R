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
    scale_y_log10() +
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

# Plot 1: All 'du' variants
plot1_data <- plot_data %>% filter(Method %in% c("du", "du2", "du3", "du4"))
plot1 <- create_benchmark_plot(plot1_data, "Performance of 'du' variants")

# Plot 2: Plain 'du' versus 'ud'
plot2_data <- plot_data %>% filter(Method %in% c("du", "ud"))
plot2 <- create_benchmark_plot(plot2_data, "Performance Comparison: 'du' vs 'ud'")

# Plot 3: Combined/iterated functions
plot3_data <- plot_data %>% filter(Method %in% c("udu", "dud", "du and ud"))
plot3 <- create_benchmark_plot(plot3_data, "Performance of Combined Functions")

print(plot1)
print(plot2)
print(plot3)

# ggsave("plot1_du_variants.png", plot1, width = 12, height = 8, dpi = 300)
# ggsave("plot2_du_vs_ud.png",  plot2, width = 12, height = 8, dpi = 300)
# ggsave("plot3_combined.png",  plot3, width = 12, height = 8, dpi = 300)

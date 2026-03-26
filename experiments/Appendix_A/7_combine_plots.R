# 2 seconds

start_time <- Sys.time()

library(ggplot2)
library(patchwork)
library(dplyr)
library(readr)
library(ragg)

source("./experiments/Appendix_A/utils.R")
source("./experiments/Appendix_A/0_parameters.R")
source("./experiments/Appendix_A/4_functions_plot.R")

load("./experiments/Appendix_A/res_data/instances.RData")

data_dir <- "./experiments/Appendix_A/res_data"
out_dir <- "./outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------
# Read data
# ------------------------------------------------------------------

summary_path <- file.path(data_dir, sprintf("summary_M%d.csv", M))
df_all <- read_csv(summary_path, show_col_types = FALSE)

expected_columns <- c(
  "p", "lambda", "alpha", "K_structure", "solver", "starting_point",
  "tolerance", "time_median", "f_end"
)
stopifnot(all(expected_columns %in% names(df_all)))

raw_path <- file.path(data_dir, sprintf("raw_M%d.csv", M))
df_raw <- read_csv(raw_path, show_col_types = FALSE)

if (any(df_raw$status != "ok")) stop("Simulation FAILED: non-ok status present")
if (anyNA(df_raw$time) || anyNA(df_raw$f_end)) stop("Simulation FAILED: NA in time or f_end")
df_raw <- dplyr::select(df_raw, -c("m", "status", "error")) %>%
  filter(p >= 50)

df_all <- df_all %>%
  mutate(
    time = time_median,
    value = f_end,
    alg = solver,
    init = starting_point
  )

group_keys_path <- file.path(data_dir, sprintf("group_keys_with_best_value_M%d.csv", M))
group_keys <- read_csv(group_keys_path, show_col_types = FALSE)
stopifnot(all(c("p", "lambda", "alpha", "K_structure", "best_value") %in% names(group_keys)))

# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------

selected_graphs <- c("hub_1", "hub_09", "AR2", "random")

graph_label <- function(K_structure) {
  switch(
    K_structure,
    "hub_1"  = "Hub, part_cor = -1 / sqrt(p)",
    "hub_09" = "Hub, part_cor = -0.9 / sqrt(p)",
    "AR2"    = "AR2",
    "random" = "Random",
    K_structure
  )
}

thr_f_diff_to_best <- 1e-5
make_type_1_single_plot <- function(i) {
  x <- prepare_data_for_plot_type_1(i)

  ggplot(x$df, aes(x = time, y = value_shifted, color = alg, shape = init)) +
    geom_point(size = 4) +
    scale_shape_manual(values = c(C = 20, I = 8)) +
    scale_y_log10(
      breaks = sort(unique(c(scales::log_breaks()(range(x$df$value_shifted, na.rm = TRUE)),
                             thr_f_diff_to_best))),
      labels = scales::label_scientific()
    ) +
    geom_hline(yintercept = thr_f_diff_to_best, linetype = "dashed") +
    expand_limits(x = 0) +
    labs(
      title = graph_label(x$df$K_structure[1]),
      x = "Time [s]",
      y = "Objective difference to best (log-scale)",
      color = "Algorithm",
      shape = "Init"
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "bottom"
    )
}

make_type_2_single_plot <- function(df_sub, K_val) {
  make_median_time_vs_p_plot(df_sub, K_val) +
    labs(
      title = graph_label(K_val),
      subtitle = NULL,
      color = "Algorithm",
      linetype = "Init"
    ) +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "bottom"
    )
}

# ------------------------------------------------------------------
# Combined figure 1: median_time_vs_p
# ------------------------------------------------------------------

df_type_2_list <- prepare_data_for_plot_type_2(df_raw, group_keys, thr_f_diff_to_best, M)

df_raw_filtered <- df_type_2_list$df_raw_filtered
df_median_time <- df_type_2_list$df_median_time

plots_type_2 <- lapply(selected_graphs, function(K_val) {
  df_sub <- df_median_time %>% filter(K_structure == K_val)
  make_type_2_single_plot(df_sub, K_val)
})

combined_type_2 <- wrap_plots(plots_type_2, ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

ggsave(
  filename = file.path(out_dir, "Appendix_A_median_time_vs_p.png"),
  plot = combined_type_2,
  device = ragg::agg_png,
  width = 12,
  height = 8,
  units = "in",
  dpi = 150
)

message("Saved: ", file.path(out_dir, "Appendix_A_median_time_vs_p.png"))

# ------------------------------------------------------------------
# Combined figure 2: selected type_1 plots
# ------------------------------------------------------------------

group_keys_sub <- group_keys %>%
  filter(
    p == 200,
    lambda == 0.1,
    alpha == 0,
    K_structure %in% selected_graphs
  ) %>%
  mutate(K_structure = factor(K_structure, levels = selected_graphs)) %>%
  arrange(K_structure)

stopifnot(nrow(group_keys_sub) == 4)

make_key <- function(df) {
  paste(
    df$p,
    sprintf("%.1f", df$lambda),
    sprintf("%.1f", df$alpha),
    df$K_structure,
    sep = "|"
  )
}

idx <- match(
  make_key(group_keys_sub),
  make_key(group_keys)
)

plots_type_1 <- lapply(idx, make_type_1_single_plot)

combined_type_1 <- wrap_plots(plots_type_1, ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

ggsave(
  filename = file.path(out_dir, "Appendix_A_p200_lambda0_1_alpha0.png"),
  plot = combined_type_1,
  device = ragg::agg_png,
  width = 12,
  height = 8,
  units = "in",
  dpi = 150
)

message("Saved: ", file.path(out_dir, "Appendix_A_p200_lambda0_1_alpha0.png"))

end_time <- Sys.time()
message("Done.")
print(end_time - start_time)

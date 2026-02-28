# 1 minute of 7 cores of Apple's M2

library(ggplot2)
library(ggpattern)
library(patchwork)
library(tidyr)
library(scales)
library(dplyr)
library(readr)
library(future.apply)

source("./experiments/Appendix_A/utils.R")
source("./experiments/Appendix_A/0_parameters.R")
source("./experiments/Appendix_A/4_functions_plot.R")

load("./experiments/Appendix_A/res_data/instances.RData")

data_dir <- "./experiments/Appendix_A/res_data"
plot_dir <- "./experiments/Appendix_A/plots"
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

n_cores <- max(1, parallel::detectCores(logical = FALSE) - 1)

# read summary
summary_path <- file.path(data_dir, sprintf("summary_M%d.csv", M))
df_all <- read_csv(summary_path, show_col_types = FALSE)
# expected columns
expected_columns <- c("p", "lambda", "alpha", "K_structure", "solver", "starting_point", "tolerance", "time_trimmed_mean", "f_end")
stopifnot(all(expected_columns %in% names(df_all)))

# read raw
raw_path <- file.path(data_dir, sprintf("raw_M%d.csv", M))
df_raw <- read_csv(raw_path, show_col_types = FALSE)
# --- Validate df_raw ---
if (any(df_raw$status != "ok")) stop("Simulation FAILED: non-ok status present")
if (anyNA(df_raw$time) || anyNA(df_raw$f_end)) stop("Simulation FAILED: NA in time or f_end")
df_raw <- dplyr::select(df_raw, -c("m", "status", "error"))


# labels used by plot
df_all <- df_all %>%
  mutate(
    time = time_trimmed_mean,
    value = f_end,
    alg = solver,
    init = starting_point
  )

group_keys <- df_all %>%
  distinct(p, lambda, alpha, K_structure) %>%
  arrange(desc(p))

# add the best_value to `group_keys`
start_time <- Sys.time()
plan(multisession, workers = n_cores)
best_values <- future_sapply(
  seq_len(nrow(group_keys)),
  compute_best_value_for_i
)
group_keys$best_value <- as.numeric(best_values)
stopifnot(all(is.finite(group_keys$best_value)))


#####
# make plots type 1:
invisible(future_lapply(
  seq_len(nrow(group_keys)),
  save_plot_for_i
))


#####
# plots type 2
thr_f_diff_to_best <- 1e-5
df_type_2_list <- prepare_type2(df_raw, group_keys, thr_f_diff_to_best, M)

df_raw_filtered <- df_type_2_list$df_raw_filtered
df_mean_time <- df_type_2_list$df_mean_time

make_mean_time_vs_p_plot <- function(df_sub, K_val) {
  df_sub <- add_lam_alp_label(df_sub)

  ggplot(
    df_sub,
    aes(
      x = p,
      y = mean_time,
      color = solver,
      linetype = starting_point,
      group = interaction(solver, starting_point)
    )
  ) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.3) +
    facet_wrap(~ lam_alp, ncol = 2, scales = "free_y") +
    labs(
      title = sprintf("Mean time vs p (%s)", K_val),
      x = "p",
      y = "Mean time [s]",
      color = "solver",
      linetype = "start"
    ) +
    theme_bw(base_size = 14) +
    theme(
      legend.position = "right",
      strip.background = element_rect(fill = "grey92"),
      panel.grid.minor = element_blank()
    )
}

K_vals <- unique(df_mean_time$K_structure)

for (K_val in K_vals) {
  df_sub <- df_mean_time %>% filter(K_structure == K_val)
  plt <- make_mean_time_vs_p_plot(df_sub, K_val)

  out_file <- file.path(
    plot_dir, "type_2",
    sprintf("mean_time_vs_p_%s.png", as.character(K_val))
  )

  ggsave(out_file, plt, width = 12, height = 6, dpi = 150)
  message("Saved: ", out_file)
}

sol_order <- c(
  "pcglasso / C", "pcglasso / I",
  "pcglassoFast_Dual / C", "pcglassoFast_Dual / I",
  "pcglassoFast_Primal / C", "pcglassoFast_Primal / I"
)
sol_order_2x3 <- c(
  "pcglasso / C", "pcglassoFast_Dual / C", "pcglassoFast_Primal / C",
  "pcglasso / I", "pcglassoFast_Dual / I", "pcglassoFast_Primal / I"
)

make_violin_time_vs_p_plot <- function(df_sub, K_val) {
  df_sub <- add_lam_alp_label(df_sub) %>%
    mutate(
      sol_start = factor(interaction(solver, starting_point, sep = " / "), levels = sol_order)
    )

  sol_levels <- levels(df_sub$solver)
  base_cols <- setNames(hue_pal()(length(sol_levels)), sol_levels)
  bright_cols <- setNames(brighten(base_cols, amount = 0.45), sol_levels)

  col_map <- c(
    setNames(base_cols,  paste0(sol_levels, " / C")),
    setNames(bright_cols, paste0(sol_levels, " / I"))
  )

  plt_main <- ggplot(df_sub, aes(x = p, y = time, fill = sol_start, col = sol_start)) +
    geom_violin(
      position = position_dodge(width = 0.85),
      trim = TRUE,
      scale = "width",
      linewidth = 0.25
    ) +
    facet_wrap(~ lam_alp, ncol = 2) +
    scale_y_log10(labels = label_number()) +
    scale_fill_manual(values = col_map) +
    scale_color_manual(values = col_map) +
    labs(
      title = sprintf("Time vs p (%s)", as.character(K_val)),
      x = "p",
      y = "Time [s]",
      fill = "solver / start"
    ) +
    theme_bw(base_size = 14) +
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey92")
    )

  my_legend <- make_matrix_legend()

  plt_main + my_legend + plot_layout(widths = c(4.7, 1))
}

K_vals <- unique(df_raw_filtered$K_structure)

for (K_val in K_vals) {
  df_sub <- df_raw_filtered %>% filter(K_structure == K_val)
  plt <- make_violin_time_vs_p_plot(df_sub, K_val)

  out_file <- file.path(
    plot_dir, "type_2",
    sprintf("violin_time_vs_p_%s.png", as.character(K_val))
  )

  ggsave(out_file, plt, width = 12, height = 7, dpi = 150)
  message("Saved: ", out_file)
}


end_time <- Sys.time()
message("Done.")
print(end_time - start_time)

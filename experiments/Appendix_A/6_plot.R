# 15 seconds

start_time <- Sys.time()

library(ggplot2)
library(patchwork)
library(scales)
library(dplyr)
library(readr)
library(ragg)

source("./experiments/Appendix_A/utils.R")
source("./experiments/Appendix_A/0_parameters.R")
source("./experiments/Appendix_A/4_functions_plot.R")

load("./experiments/Appendix_A/res_data/instances.RData")

data_dir <- "./experiments/Appendix_A/res_data"
plot_dir <- "./experiments/Appendix_A/plots"
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(plot_dir, "type_1"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(plot_dir, "type_2"), showWarnings = FALSE, recursive = TRUE)

# read summary
summary_path <- file.path(data_dir, sprintf("summary_M%d.csv", M))
df_all <- read_csv(summary_path, show_col_types = FALSE)
# expected columns
expected_columns <- c("p", "lambda", "alpha", "K_structure", "solver", "starting_point", "tolerance", "time_median", "f_end")
stopifnot(all(expected_columns %in% names(df_all)))

# read raw
raw_path <- file.path(data_dir, sprintf("raw_M%d.csv", M))
df_raw <- read_csv(raw_path, show_col_types = FALSE)
# --- Validate df_raw ---
if (any(df_raw$status != "ok")) stop("Simulation FAILED: non-ok status present")
if (anyNA(df_raw$time) || anyNA(df_raw$f_end)) stop("Simulation FAILED: NA in time or f_end")
df_raw <- dplyr::select(df_raw, -c("m", "status", "error"))


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


#####
# make plots type 1:
for (i in seq_len(nrow(group_keys))) {
  x <- prepare_data_for_plot_type_1(i)

  plt <- ggplot(x$df, aes(x = time, y = value_shifted, color = alg, shape = init)) +
    geom_point(size = 4) +
    scale_shape_manual(values = c(C = 20, I = 8)) +
    scale_y_log10() +
    expand_limits(x = 0) +
    labs(
      title = "PCGLASSO vs pcglassoFast Dual vs pcglassoFast Primal",
      subtitle = x$subtitle,
      x = "Median Time [s]",
      y = "Objective difference to best (log-scale)"
    ) +
    theme_bw(base_size = 14)

  ggsave(x$file, plt, device = ragg::agg_png, width = 9, height = 6, units = "in", dpi = 150)
  message("Saved: ", x$file)
}


#####
# plots type 2
thr_f_diff_to_best <- 1e-5
df_type_2_list <- prepare_data_for_plot_type_2(df_raw, group_keys, thr_f_diff_to_best, M)

df_raw_filtered <- df_type_2_list$df_raw_filtered
df_median_time <- df_type_2_list$df_median_time

K_vals <- unique(df_median_time$K_structure)

for (K_val in K_vals) {
  df_sub <- df_median_time %>% filter(K_structure == K_val)
  plt <- make_median_time_vs_p_plot(df_sub, K_val)

  out_file <- file.path(
    plot_dir, "type_2",
    sprintf("median_time_vs_p_%s.png", as.character(K_val))
  )

  ggsave(out_file, plt, width = 12, height = 6, dpi = 150)
  message("Saved: ", out_file)
}

sol_order <- c(
  "pcglasso / C", "pcglasso / I",
  "pcglassoFast_Dual / C", "pcglassoFast_Dual / I",
  "pcglassoFast_Primal / C", "pcglassoFast_Primal / I",
  "pcglassoFast_PrimalDual / C", "pcglassoFast_PrimalDual / I",
  "path_up / C", "path_up / I",
  "path_down / C", "path_down / I"
)
sol_order_2x4 <- c(
  "pcglasso / C", "pcglassoFast_Dual / C", "pcglassoFast_Primal / C", "pcglassoFast_PrimalDual / C",
  "pcglasso / I", "pcglassoFast_Dual / I", "pcglassoFast_Primal / I", "pcglassoFast_PrimalDual / I"
)

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

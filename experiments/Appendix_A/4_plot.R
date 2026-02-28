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
source("./experiments/Appendix_A/2_functions.R")

load("./experiments/Appendix_A/res_data/instances.RData")

data_dir <- "./experiments/Appendix_A/res_data"
plot_dir <- "./experiments/Appendix_A/plots"
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

n_cores <- max(1, min(7, parallel::detectCores(logical = FALSE) - 1))

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

# make plots type 1:
invisible(future_lapply(
  seq_len(nrow(group_keys)),
  save_plot_for_i
))




# plots type 2
thr_f_diff_to_best <- 1e-5
# --- Consistency check: within each group f_end is constant ---
if (
  any(df_raw %>%
        group_by(p, lambda, alpha, K_structure, solver, starting_point, tolerance) %>%
        summarise(span = max(f_end) - min(f_end), .groups = "drop") %>%
        pull(span) != 0)) {
  stop("Some simulations arrived at different f_end within the same group; analysis compromised.")
}

df_raw_no_time <- df_raw %>%
  dplyr::select(-time) %>%
  distinct()

stopifnot(nrow(df_raw_no_time) == (
  length(unique(df_raw_no_time$p)) *
    length(unique(df_raw_no_time$lambda)) *
    length(unique(df_raw_no_time$alpha)) *
    length(unique(df_raw_no_time$K_structure)) *
    length(unique(df_raw_no_time$solver)) *
    length(unique(df_raw_no_time$starting_point)) *
    length(unique(df_raw_no_time$tolerance))
))

df_diff_to_best <- df_raw_no_time %>%
  left_join(
    group_keys,
    by = c("p", "lambda", "alpha", "K_structure")
  ) %>%
  mutate(
    f_diff_to_best = f_end - best_value,
    .keep = "unused"
  )
stopifnot(all(df_diff_to_best$f_diff_to_best >= 0))

df_tol_found <- df_diff_to_best %>%
  group_by(p, lambda, alpha, K_structure, solver, starting_point) %>%
  summarise(
    tol_found = {
      ok <- tolerance[f_diff_to_best < thr_f_diff_to_best]
      if (length(ok) == 0) NA_real_ else max(ok)
    },
    .groups = "drop"
  )
stopifnot(all(!is.na(df_tol_found$tol_found)))

grouping_rows <- c("p", "lambda", "alpha", "K_structure", "solver", "starting_point")
df_raw_filtered <- df_raw %>%
  inner_join(df_tol_found, by = grouping_rows) %>%
  filter(!is.na(tol_found), tolerance == tol_found) %>%
  dplyr::select(-c(tol_found, tolerance, f_end)) %>%
  mutate(
    p = factor(p, levels = sort(unique(p))),
    K_structure = factor(K_structure),
    solver = factor(solver),
    starting_point = factor(starting_point, levels = c("C", "I")),
    lambda = as.numeric(lambda),
    alpha = as.numeric(alpha)
  )
stopifnot(nrow(df_raw_filtered) == (
  length(unique(df_raw_no_time$p)) *
    length(unique(df_raw_no_time$lambda)) *
    length(unique(df_raw_no_time$alpha)) *
    length(unique(df_raw_no_time$K_structure)) *
    length(unique(df_raw_no_time$solver)) *
    length(unique(df_raw_no_time$starting_point)) *
    M
))

trim <- if (M >= 10) { 0.1 } else { 0 }
df_mean_time <- df_raw_filtered %>%
  group_by(p, lambda, alpha, K_structure, solver, starting_point) %>%
  summarise(
    mean_time = mean(time, trim = trim),
    .groups = "drop"
  ) %>%
  mutate(
    p = factor(p, levels = sort(unique(p))),
    K_structure = factor(K_structure, levels = sort(unique(K_structure))),
    solver = factor(solver),
    starting_point = factor(starting_point, levels = c("C", "I")),
    lambda = as.numeric(lambda),
    alpha = as.numeric(alpha)
  )

make_one_plot <- function(df_sub, K_val) {
  df_sub <- df_sub %>%
    mutate(
      lam_alp = factor(
        sprintf("lambda = %.2g, alpha = %.2g", lambda, alpha),
        levels = (df_sub %>% distinct(lambda, alpha) %>%
                    arrange(lambda, alpha) %>%
                    transmute(lab = sprintf("lambda = %.2g, alpha = %.2g", lambda, alpha)) %>%
                    pull(lab))
      )
    )

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
  plt <- make_one_plot(df_sub, K_val)

  out_file <- file.path(
    plot_dir, "type_2",
    sprintf("mean_time_vs_p_%s.png", as.character(K_val))
  )

  ggsave(out_file, plt, width = 12, height = 6, dpi = 150)
  message("Saved: ", out_file)
}

brighten <- function(cols, amount = 0.45) {
  rgb <- grDevices::col2rgb(cols)
  rgb2 <- rgb + (255 - rgb) * amount
  grDevices::rgb(rgb2[1,], rgb2[2,], rgb2[3,], maxColorValue = 255)
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

make_matrix_legend <- function(base_cols, bright_cols) {
  solvers <- c("pcglasso","pcglassoFast_Dual","pcglassoFast_Primal")

  legend_df <- tibble::tibble(
    solver = factor(solvers, levels = solvers),
    C = unname(base_cols[solvers]),
    I = unname(bright_cols[solvers])
  ) %>%
    tidyr::pivot_longer(c(C, I), names_to = "start", values_to = "col") %>%
    mutate(start = factor(start, levels = c("C","I")))

  ggplot(legend_df, aes(x = solver, y = start)) +
    geom_tile(aes(fill = col),
              width = 0.85, height = 0.85,
              color = "grey40", linewidth = 0.3) +
    scale_fill_identity() +
    scale_y_discrete(limits = c("I","C"), expand = expansion(mult = c(0, 0))) +
    scale_x_discrete(expand = expansion(mult = c(0.05, 0.05))) +
    coord_fixed(ratio = 1, clip = "off") +
    annotate("text", x = 2, y = 2.65, label = "starting point / algorithm",
             size = 4, hjust = 0.5) +
    expand_limits(y = 2.8) +
    annotate("text", x = 1, y = 0.05, label = "pcglasso",            angle = 90, size = 4) +
    annotate("text", x = 2, y = -0.5, label = "pcglassoFast_Dual",   angle = 90, size = 4) +
    annotate("text", x = 3, y = -0.6, label = "pcglassoFast_Primal", angle = 90, size = 4) +
    expand_limits(y = 0.35) +
    theme_void(base_size = 11) +
    theme(
      axis.text.x = element_blank(),
      plot.margin = margin(0, 0, 0, 0),
      axis.text.y = element_text(
        margin = margin(r = 2),
        size = 12
      )
    )
}

make_one_violin <- function(df_sub, K_val) {
  df_sub <- df_sub %>%
    mutate(
      lam_alp = factor(
        sprintf("lambda = %.2g, alpha = %.2g", lambda, alpha),
        levels = (df_sub %>% distinct(lambda, alpha) %>%
                    arrange(lambda, alpha) %>%
                    transmute(lab = sprintf("lambda = %.2g, alpha = %.2g", lambda, alpha)) %>%
                    pull(lab))
      ),
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

  my_legend <- make_matrix_legend(base_cols, bright_cols)

  plt_main + my_legend + plot_layout(widths = c(4.7, 1))
}

K_vals <- unique(df_raw_filtered$K_structure)

for (K_val in K_vals) {
  df_sub <- df_raw_filtered %>% filter(K_structure == K_val)
  plt <- make_one_violin(df_sub, K_val)

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

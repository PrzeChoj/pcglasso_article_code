# 1 minute of 7 cores of Apple's M2

library(ggplot2)
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

# make plots:
invisible(future_lapply(
  seq_len(nrow(group_keys)),
  save_plot_for_i
))
end_time <- Sys.time()

message("Done.")
print(end_time - start_time)




# Other plot
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

thr_f_diff_to_best <- 1e-5
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
  dplyr::select(-c(tol_found, tolerance, f_end))
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

make_one_plot <- function(df_sub, lambda_val, alpha_val) {
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
    facet_wrap(~ K_structure, ncol = 2, scales = "free_y") +
    labs(
      title = sprintf("Mean time vs p (lambda = %.2g, alpha = %.2g)", lambda_val, alpha_val),
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

# 4 plots: one per (lambda, alpha)
pairs <- df_mean_time %>%
  distinct(lambda, alpha) %>%
  arrange(lambda, alpha)

plots <- vector("list", nrow(pairs))

for (k in seq_len(nrow(pairs))) {
  lam <- pairs$lambda[k]
  alp <- pairs$alpha[k]

  df_sub <- df_mean_time %>% filter(lambda == lam, alpha == alp)
  plt <- make_one_plot(df_sub, lam, alp)

  out_file <- file.path(
    plot_dir, "type_2",
    sprintf("mean_time_vs_p_lambda%s_alpha%s.png",
            gsub("\\.", "_", as.character(lam)),
            gsub("\\.", "_", as.character(alp)))
  )

  ggsave(out_file, plt, width = 12, height = 6, dpi = 150)
  message("Saved: ", out_file)
}

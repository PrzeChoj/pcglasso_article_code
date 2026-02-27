library(ggplot2)
library(dplyr)
library(readr)

source("./experiments/Appendix_A/experiment_1/0_parameters.R")
source("./experiments/Appendix_A/utils.R")

load("./experiments/Appendix_A/experiment_1/res_data/instances.RData")

data_dir <- "./experiments/Appendix_A/experiment_1/res_data"
plot_dir <- "./experiments/Appendix_A/experiment_1/plots"
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# read summary
summary_path <- file.path(data_dir, sprintf("experiment_1_summary_M%d.csv", M))
df_all <- read_csv(summary_path, show_col_types = FALSE)

# expected columns
expected_columns <- c("p", "lambda", "alpha", "K_structure", "solver", "starting_point", "tolerance", "time_trimmed_mean", "f_end")
stopifnot(all(expected_columns %in% names(df_all)))

# labels used by plot
df_all <- df_all %>%
  mutate(
    time = time_trimmed_mean,
    value = f_end,
    alg = solver,
    init = starting_point
  )

fmt_part_cor <- function(K_structure) {
  switch (K_structure,
    "hub_1" = "part_cor = -1 / sqrt(p)",
    "hub_09" = "part_cor = -0.9 / sqrt(p)",
    "AR2" = "a",
    "random" = "r",
    stop("Unknown K_structure: ", K_structure)
  )
}

group_keys <- df_all %>%
  distinct(p, lambda, alpha, K_structure)

for (i in seq_len(nrow(group_keys))) {
  info <- group_keys[i, ]

  df <- df_all %>%
    filter(
      p == info$p,
      lambda == info$lambda,
      alpha == info$alpha,
      K_structure == info$K_structure
    )

  # baseline best value
  S <- instances[[info$K_structure]][[as.character(info$p)]]
  if (is.null(S)) stop("Missing S for: ", K_structure, " p=", p)
  best_value <- get_best_value(S, info$p, info$K_structure, info$lambda, info$alpha)

  df <- df %>%
    mutate(value_shifted = value - best_value + 1e-12)

  stopifnot(all(df$value_shifted > 0))

  plt <- ggplot(df, aes(x = time, y = value_shifted, color = alg, shape = init)) +
    geom_point(size = 4) +
    scale_shape_manual(values = c(C = 20, I = 8)) +
    theme_minimal(base_size = 14) +
    scale_y_log10() +
    expand_limits(x = 0) +
    labs(
      title = "PCGLASSO vs pcglassoFast Dual vs pcglassoFast Primal",
      subtitle = sprintf(
        "p = %d   |   %s graph   |   %s   |   lambda = %.1f   |   alpha = %.1f",
        info$p, info$K_structure,
        fmt_part_cor(info$K_structure),
        info$lambda, info$alpha
      ),
      x = "Time [s] (trimmed mean over M)",
      y = "Objective difference to best (log-scale)"
    )

  plot_file <- sprintf(
    "%s/plot_%s_p%d_lambda%s_alpha%s.pdf",
    plot_dir,
    info$K_structure,
    info$p,
    gsub("\\.", "_", as.character(info$lambda)),
    gsub("\\.", "_", as.character(info$alpha))
  )

  ggsave(plot_file, plt, width = 8, height = 6, dpi = 150)
  message("Saved: ", plot_file)
}

message("Done.")

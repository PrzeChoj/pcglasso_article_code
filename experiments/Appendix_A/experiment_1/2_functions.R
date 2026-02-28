source("./experiments/Appendix_A/utils.R")
source("./experiments/Appendix_A/experiment_1/0_parameters.R")

load("./experiments/Appendix_A/experiment_1/res_data/instances.RData")

# For simulation:
run_one <- function(row) {
  p <- row$p
  lambda <- row$lambda
  alpha <- row$alpha
  K_structure <- row$K_structure
  solver <- row$solver
  starting_point <- row$starting_point
  tolerance <- row$tolerance

  S <- instances[[K_structure]][[as.character(p)]]
  if (is.null(S)) stop("Missing S for: ", K_structure, " p=", p)

  if (solver == "pcglasso") {
    tolerance <- tolerance * pcglasso_tolerance_multiplier
  }

  res <- tryCatch(
    {
      t0 <- proc.time()[["elapsed"]]
      f_end <- value_after_optimization(S, solver, starting_point, tolerance, lambda, alpha)
      t1 <- proc.time()[["elapsed"]]
      list(status = "ok", time = t1 - t0, f_end = f_end, error = NA_character_)
    },
    error = function(e) {
      list(status = "fail", time = NA_real_, f_end = NA_real_, error = conditionMessage(e))
    }
  )

  c(as.list(row), res)
}

# For plotting:
compute_best_value_for_i <- function(i) {
  info <- group_keys[i, ]
  message("Optimize for p = ", info$p, ", lambda = ", info$lambda, ", alpha = ", info$alpha, ", K_structure = ", info$K_structure)
  S <- instances[[info$K_structure]][[as.character(info$p)]]
  if (is.null(S)) {
    stop("Missing S for: ", info$K_structure, " p=", info$p)
  }
  best_value_calculated <- get_best_value(S, info$p, info$K_structure, info$lambda, info$alpha)

  best_value_from_simulations <- df_raw %>%
    filter(p == info$p, lambda == info$lambda, alpha == info$alpha, K_structure == info$K_structure) %>%
    pull(f_end) %>%
    min

  min(best_value_calculated, best_value_from_simulations)
}

fmt_part_cor <- function(K_structure) {
  switch(K_structure,
    "hub_1" = "part_cor = -1 / sqrt(p)",
    "hub_09" = "part_cor = -0.9 / sqrt(p)",
    "AR2" = "a",
    "random" = "r",
    stop("Unknown K_structure: ", K_structure)
  )
}

save_plot_for_i <- function(i) {
  info <- group_keys[i, ]

  if (is.na(info$best_value) || !is.finite(info$best_value)) {
    stop("Missing/invalid best_value for row i = ", i)
  }

  df <- df_all %>%
    filter(
      p == info$p,
      lambda == info$lambda,
      alpha == info$alpha,
      K_structure == info$K_structure
    )

  # shift value based on precomputed best_value
  best_value <- info$best_value
  df <- df %>%
    mutate(value_shifted = value - best_value)
  if (any(df$value_shifted <= 1e-12)) {
    df$value_shifted <- df$value_shifted + 1e-10
  }
  stopifnot(all(df$value_shifted > 0))

  # labels
  graph_label <- if (info$K_structure %in% c("hub_1", "hub_09")) {
    "Graph hub"
  } else {
    sprintf("Graph %s", info$K_structure)
  }
  part_cor_label <- if (info$K_structure %in% c("hub_1", "hub_09")) {
    fmt_part_cor(info$K_structure)
  } else {
    NULL
  }
  subtitle_txt <- paste(
    c(
      sprintf("p = %d", info$p),
      graph_label,
      part_cor_label,
      sprintf("lambda = %.1f", info$lambda),
      sprintf("alpha = %.1f", info$alpha)
    ),
    collapse = "   |   "
  )

  # make graph
  plt <- ggplot(df, aes(x = time, y = value_shifted, color = alg, shape = init)) +
    geom_point(size = 4) +
    scale_shape_manual(values = c(C = 20, I = 8)) +
    scale_y_log10() +
    expand_limits(x = 0) +
    labs(
      title = "PCGLASSO vs pcglassoFast Dual vs pcglassoFast Primal",
      subtitle = subtitle_txt,
      x = "Time [s]",
      y = "Objective difference to best (log-scale)"
    ) +
    theme_bw(base_size = 14)

  plot_file <- sprintf(
    "%s/plot_%s_p%d_lambda%s_alpha%s.png",
    plot_dir,
    info$K_structure,
    info$p,
    gsub("\\.", "_", as.character(info$lambda)),
    gsub("\\.", "_", as.character(info$alpha))
  )

  ggsave(plot_file, plt, width = 9, height = 6, dpi = 150)
  message("Saved: ", plot_file)

  NULL
}

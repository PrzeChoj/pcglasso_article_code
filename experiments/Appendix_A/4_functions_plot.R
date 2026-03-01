source("./experiments/Appendix_A/utils.R")
source("./experiments/Appendix_A/0_parameters.R")

load("./experiments/Appendix_A/res_data/instances.RData")

brighten <- function(cols, amount = 0.45) {
  rgb <- grDevices::col2rgb(cols)
  rgb2 <- rgb + (255 - rgb) * amount
  grDevices::rgb(rgb2[1, ], rgb2[2, ], rgb2[3, ], maxColorValue = 255)
}

make_matrix_legend <- function() {
  solvers <- c("pcglasso", "pcglassoFast_Dual", "pcglassoFast_Primal")

  base_cols <- setNames(scales::hue_pal()(length(solvers)), solvers)
  bright_cols <- setNames(brighten(base_cols, amount = 0.45), solvers)

  legend_df <- tibble::tibble(
    solver = factor(solvers, levels = solvers),
    C = unname(base_cols[solvers]),
    I = unname(bright_cols[solvers])
  ) %>%
    tidyr::pivot_longer(c(C, I), names_to = "start", values_to = "col") %>%
    mutate(start = factor(start, levels = c("C", "I")))

  ggplot(legend_df, aes(x = solver, y = start)) +
    geom_tile(aes(fill = col),
      width = 0.85, height = 0.85,
      color = "grey40", linewidth = 0.3
    ) +
    scale_fill_identity() +
    scale_y_discrete(limits = c("I", "C"), expand = expansion(mult = c(0, 0))) +
    scale_x_discrete(expand = expansion(mult = c(0.05, 0.05))) +
    coord_fixed(ratio = 1, clip = "off") +
    annotate("text",
      x = 2, y = 2.65, label = "starting point / algorithm",
      size = 4, hjust = 0.5
    ) +
    expand_limits(y = 2.8) +
    annotate("text", x = 1, y = 0.05, label = "pcglasso", angle = 90, size = 4) +
    annotate("text", x = 2, y = -0.5, label = "pcglassoFast_Dual", angle = 90, size = 4) +
    annotate("text", x = 3, y = -0.6, label = "pcglassoFast_Primal", angle = 90, size = 4) +
    expand_limits(y = 0.35) +
    theme_void(base_size = 11) +
    theme(
      axis.text.x = element_blank(),
      plot.margin = margin(0, 0, 0, 0),
      axis.text.y = element_text(margin = margin(r = 2), size = 12)
    )
}

add_lam_alp_label <- function(df_sub) {
  df_sub %>%
    mutate(
      lam_alp = factor(
        sprintf("lambda = %.2g, alpha = %.2g", lambda, alpha),
        levels = (df_sub %>% distinct(lambda, alpha) %>%
          arrange(lambda, alpha) %>%
          transmute(lab = sprintf("lambda = %.2g, alpha = %.2g", lambda, alpha)) %>%
          pull(lab))
      )
    )
}

prepare_data_for_plot_type_1 <- function(i) {
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
    ) %>%
    mutate(value_shifted = pmax(value - info$best_value, 1e-10))

  graph_label <- if (info$K_structure %in% c("hub_1", "hub_09")) "Graph hub" else sprintf("Graph %s", info$K_structure)
  part_cor_label <- if (info$K_structure %in% c("hub_1", "hub_09")) fmt_part_cor(info$K_structure) else NULL

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

  plot_file <- sprintf(
    "%s/type_1/plot_%s_p%d_lambda%s_alpha%s.png",
    plot_dir,
    info$K_structure,
    info$p,
    gsub("\\.", "_", as.character(info$lambda)),
    gsub("\\.", "_", as.character(info$alpha))
  )

  list(df = df, subtitle = subtitle_txt, file = plot_file)
}

prepare_data_for_plot_type_2 <- function(df_raw, group_keys, thr, M) {
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
        ok <- tolerance[f_diff_to_best < thr]
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

  list(
    df_raw_filtered = df_raw_filtered,
    df_mean_time = df_mean_time
  )
}

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
    min()

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
    facet_wrap(~lam_alp, ncol = 2, scales = "free_y") +
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

make_violin_time_vs_p_plot <- function(df_sub, K_val) {
  df_sub <- add_lam_alp_label(df_sub) %>%
    mutate(
      sol_start = factor(interaction(solver, starting_point, sep = " / "), levels = sol_order)
    )

  sol_levels <- levels(df_sub$solver)
  base_cols <- setNames(hue_pal()(length(sol_levels)), sol_levels)
  bright_cols <- setNames(brighten(base_cols, amount = 0.45), sol_levels)

  col_map <- c(
    setNames(base_cols, paste0(sol_levels, " / C")),
    setNames(bright_cols, paste0(sol_levels, " / I"))
  )

  plt_main <- ggplot(df_sub, aes(x = p, y = time, fill = sol_start, col = sol_start)) +
    geom_violin(
      position = position_dodge(width = 0.85),
      trim = TRUE,
      scale = "width",
      linewidth = 0.25
    ) +
    facet_wrap(~lam_alp, ncol = 2) +
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

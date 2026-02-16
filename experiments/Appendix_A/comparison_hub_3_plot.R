library(ggplot2)
library(stringr)

source("./experiments/Appendix_A/comparison_hub_1_utils.R")

data_dir <- "./experiments/Appendix_A/res_data"
plot_dir <- "./experiments/Appendix_A/plots"
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# best method table
hub_methods <- c(
  "pcglasso_C",
  "pcglasso_cpp_C",
  "pcglasso_cpp_I",
  "pcglassoFast_I",
  "pcglasso_cpp_I",
  "pcglassoFast_C",
  "pcglassoFast_I",
  "pcglasso_I",
  "pcglasso_C",
  "pcglassoFast_C",
  "pcglassoFast_C",
  "pcglassoFast_C",
  "pcglassoFast_I",
  "pcglassoFast_C",
  "pcglassoFast_I",
  "pcglasso_I"
)
stopifnot(length(hub_methods) == 16)
grid_hub <- expand.grid(
  p            = c(50, 70),
  cor_modifier = c(1.0, 0.9),
  lambda       = c(0.1, 0.2),
  alpha        = c(0.0, 0.5),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

grid_line <- expand.grid(
  p            = c(50, 70),
  cor_modifier = c(0.8, 0.9),
  lambda       = c(0.1, 0.2),
  alpha        = c(0.0, 0.5),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

best_method_table <- rbind(
  transform(grid_line,
            K_structure = "line",
            best_method = "pcglassoFast_C"),
  transform(grid_hub,
            K_structure = "hub",
            best_method = hub_methods)
)

files <- list.files(
  data_dir,
  pattern = "^comparison_(hub|line)_M[0-9]+_p[0-9]+_cor.*_lambda.*_alpha.*\\.rds$",
  full.names = TRUE
)

parse_filename <- function(path) {
  fname <- basename(path)

  M <- as.numeric(str_match(fname, "_M([0-9]+)_")[, 2])
  p <- as.numeric(str_match(fname, "_p([0-9]+)_")[, 2])

  cor_str     <- str_match(fname, "_cor(.+?)_lambda")[, 2]
  lambda_str  <- str_match(fname, "_lambda(.+?)_alpha")[, 2]
  alpha_str   <- str_match(fname, "_alpha(.+?)\\.rds$")[, 2]
  K_structure <- str_match(fname, "^comparison_(hub|line)_")[, 2]

  cor_modifier <- as.numeric(gsub("_", ".", cor_str))
  lambda       <- as.numeric(gsub("_", ".", lambda_str))
  alpha        <- as.numeric(gsub("_", ".", alpha_str))

  list(
    M = M,
    p = p,
    cor_modifier = cor_modifier,
    lambda = lambda,
    alpha = alpha,
    K_structure = K_structure
  )
}

fmt_part_cor <- function(p, cor_modifier, K_structure) {
  if (K_structure == "hub") {
    sprintf("part_cor = -%.1f / sqrt(p)", cor_modifier)
  } else if (K_structure == "line") {
    sprintf("part_cor = %.1f * max", cor_modifier)
  } else {
    stop("Unknown K_structure: ", K_structure)
  }
}

for (file in files) {
  info <- parse_filename(file)
  df   <- readRDS(file)

  # overwrite metadata from filename (safer)
  df$p            <- info$p
  df$cor_modifier <- info$cor_modifier
  df$lambda       <- info$lambda
  df$alpha        <- info$alpha
  df$K_structure  <- info$K_structure

  row_best <- subset(
    best_method_table,
    p            == info$p &
    cor_modifier == info$cor_modifier &
    lambda       == info$lambda &
    alpha        == info$alpha &
    K_structure  == info$K_structure
  )
  if (nrow(row_best) != 1L) {
    stop("best_method not found or not unique for: p=",
         info$p, " cor=", info$cor_modifier,
         " lambda=", info$lambda, " alpha=", info$alpha,
         " K=", info$K_structure)
  }
  best_method <- row_best$best_method[1]
  best_value <- compute_best_value(
    p            = info$p,
    cor_modifier = info$cor_modifier,
    lambda       = info$lambda,
    alpha        = info$alpha,
    best_method  = best_method,
    K_structure  = info$K_structure,
    tolerance_best = 1e-12,
    seed         = 1234
  )

  df$value_shifted <- df$value - best_value + 1e-12
  if (any(df$value_shifted == 0)) {
    df$value_shifted <- df$value_shifted
  }
  stopifnot(all(df$value_shifted > 0)) # when error, make tolerance_best smaller

  # ---------------- plotting ----------------
  plt <- ggplot(df, aes(x = time, y = value_shifted,
                        color = alg, shape = init)) +
    geom_point(size = 4) +
    scale_shape_manual(values = c(C = 20, I = 8)) +
    theme_minimal(base_size = 14) +
    theme(
      panel.background = element_rect(fill = "white", colour = "white"),
      plot.background  = element_rect(fill = "white", colour = "white")
    ) +
    scale_y_log10() +
    expand_limits(x = 0) +
    labs(
      title = "PCGLASSO vs PCGLASSOFast vs PCGLASSOcpp",
      subtitle = sprintf(
        "p = %d   |   K = %s   |   %s   |   lambda = %.1f   |   alpha = %.1f",
        info$p, info$K_structure,
        fmt_part_cor(info$p, info$cor_modifier, info$K_structure),
        info$lambda, info$alpha
      ),
      x = "Time [s]",
      y = "Objective difference to best (log-scale)"
    )

  # -------- save plot --------
  plot_file <- sprintf(
    "%s/plot_%s_p%d_cor%s_lambda%s_alpha%s.png",
    plot_dir,
    info$K_structure,
    info$p,
    gsub("\\.", "_", as.character(info$cor_modifier)),
    gsub("\\.", "_", as.character(info$lambda)),
    gsub("\\.", "_", as.character(info$alpha))
  )

  ggsave(plot_file, plt, width = 8, height = 6, dpi = 150)
  message("Saved: ", plot_file)
}

message("Done.")

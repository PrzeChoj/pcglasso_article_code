library(ggplot2)
library(stringr)

data_dir <- "./experiments/Appendix_A/res_data"
plot_dir <- "./experiments/Appendix_A/plots"
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

files <- list.files(
  data_dir,
  pattern = "^comparison_hub_M[0-9]+_p[0-9]+_cor.*_lambda.*_alpha.*\\.rds$",
  full.names = TRUE
)

parse_filename <- function(path) {
  fname <- basename(path)

  M <- as.numeric(str_match(fname, "_M([0-9]+)_")[, 2])
  p <- as.numeric(str_match(fname, "_p([0-9]+)_")[, 2])

  cor_str    <- str_match(fname, "_cor(.+?)_lambda")[, 2]
  lambda_str <- str_match(fname, "_lambda(.+?)_alpha")[, 2]
  alpha_str  <- str_match(fname, "_alpha(.+?)\\.rds$")[, 2]

  cor_modifier <- as.numeric(gsub("_", ".", cor_str))
  lambda       <- as.numeric(gsub("_", ".", lambda_str))
  alpha        <- as.numeric(gsub("_", ".", alpha_str))

  list(
    M = M,
    p = p,
    cor_modifier = cor_modifier,
    lambda = lambda,
    alpha = alpha
  )
}

for (file in files) {

  info <- parse_filename(file)
  df   <- readRDS(file)

  # overwrite metadata from filename (safer)
  df$p            <- info$p
  df$cor_modifier <- info$cor_modifier
  df$lambda       <- info$lambda
  df$alpha        <- info$alpha

  # shift objective
  eps <- 1e-7
  df$value_shifted <- df$value - min(df$value) + eps

  # ---------------- plotting ----------------
  plt <- ggplot(df, aes(x = time, y = value_shifted,
                        color = alg, shape = init)) +
    geom_point(size = 4) +
    scale_shape_manual(values = c(C = 20, I = 8)) +
    theme_minimal(base_size = 14) +
    scale_y_log10() +
    expand_limits(x = 0) +
    labs(
      title = "PCGLASSO vs PCGLASSOFast vs PCGLASSOcpp",
      subtitle = sprintf(
        "p = %d   |   corr = -%.1f / sqrt(p)   |   lambda = %.1f   |   alpha = %.1f",
        info$p, info$cor_modifier, info$lambda, info$alpha
      ),
      x = "Time [s]",
      y = "Objective (shifted, log-scale)"
    )

  # -------- save plot --------
  plot_file <- sprintf(
    "%s/plot_p%d_cor%s_lambda%s_alpha%s.png",
    plot_dir,
    info$p,
    gsub("\\.", "_", as.character(info$cor_modifier)),
    gsub("\\.", "_", as.character(info$lambda)),
    gsub("\\.", "_", as.character(info$alpha))
  )

  ggsave(plot_file, plt, width = 8, height = 6, dpi = 150)
  message("Saved: ", plot_file)
}

message("Done.")

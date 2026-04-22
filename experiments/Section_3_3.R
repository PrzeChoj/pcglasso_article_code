# 2 minutes on 7 cores of Apple M2

library(pcglassoFast)

library(parallel)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(dplyr)


which_graph <- 1 # 1 -> "hub", 2 -> "KAR"

num_of_cores <- 7

p <- 15
b <- 1

Khub <- function(p, a, b, c) {
  K <- matrix(0, nrow = p, ncol = p)
  K[, 1] <- K[1, ] <- c
  diag(K) <- b
  K[1, 1] <- a
  K
}

is_Khub_positive_definite <- function(p, a, b, c) {
  a > 0 && a * b > (p - 1) * c * c
}

KAR <- function(p, a, b, c) {
  K <- matrix(0, nrow = p, ncol = p)
  diag(K) <- b
  K[1, 1] <- a
  K[p, p] <- a

  for (i in 2:p) {
    K[i - 1, i] <- K[i, i - 1] <- c
  }

  K
}

is_KAR_positive_definite <- function(p, a, b, c) {
  K <- KAR(p, a, b, c)
  min(eigen(K, symmetric = TRUE, only.values = TRUE)$values) > 1e-5
}

Kgenerator <- c(
  Khub, KAR
)[[which_graph]]

is_K_positive_definite <- c(
  is_Khub_positive_definite, is_KAR_positive_definite
)[[which_graph]]

Ktext <- c(
  "hub", "KAR"
)[which_graph]

single_irrep_value <- function(p, a, b, c) {
  if (!is_K_positive_definite(p, a, b, c)) {
    return(c(NA_real_, NA_real_))
  }

  K <- Kgenerator(p, a, b, c)
  c(irrepGLASSO(K), irrepPCGLASSO(K))
}

all_a <- list(
  (1:200) * 0.5,
  (1:200) * 0.05
)[[which_graph]]

all_c <- list(
  (-120:120) / 60,
  (-63:63) / 120
)[[which_graph]]

my_res_GLASSO <- matrix(nrow = length(all_c), ncol = length(all_a))
my_res_PCGLASSO <- matrix(nrow = length(all_c), ncol = length(all_a))

for (i_for_a in seq_along(all_a)) {
  print(paste0(i_for_a, " of ", length(all_a)))

  my_results <- if (num_of_cores == 1) {
    lapply(all_c, function(c_value) {
      single_irrep_value(p, all_a[i_for_a], b, c_value)
    })
  } else {
    mclapply(all_c, function(c_value) {
      single_irrep_value(p, all_a[i_for_a], b, c_value)
    }, mc.cores = num_of_cores)
  }

  for (i_for_c in seq_along(all_c)) {
    my_res <- my_results[[i_for_c]]
    my_res_GLASSO[i_for_c, i_for_a] <- my_res[1]
    my_res_PCGLASSO[i_for_c, i_for_a] <- my_res[2]
  }
}

make_irrep_df <- function(mat, all_a, all_c) {
  df <- reshape2::melt(mat, varnames = c("c_id", "a_id"), value.name = "value")
  names(df) <- c("c_id", "a_id", "value")
  df$c <- all_c[df$c_id]
  df$a <- all_a[df$a_id]
  df
}

irrep_axis_breaks <- function(all_a, all_c, a_every = 40, c_breaks = NULL) {
  # Create evenly spaced breaks at multiples of a_every, starting from 0
  y_breaks <- seq(0, max(all_a), by = 20)
  # Keep only breaks that exist in all_a or are close to it
  y_breaks <- y_breaks[y_breaks <= max(all_a)]

  x_breaks <- if (is.null(c_breaks)) pretty(all_c, n = 7) else c_breaks

  list(x = x_breaks, y = y_breaks)
}

plot_irrep_publication <- function(mat,
                                   all_a,
                                   all_c,
                                   panel_title,
                                   fill_limits = NULL,
                                   midpoint = 1,
                                   c_breaks = NULL,
                                   a_every = 40,
                                   legend_title = "IRR") {
  df <- make_irrep_df(mat, all_a, all_c)
  brks <- irrep_axis_breaks(all_a, all_c, a_every = a_every, c_breaks = c_breaks)

  if (is.null(fill_limits)) {
    fill_limits <- c(0, max(df$value, na.rm = TRUE))
  }

  ggplot(df, aes(x = c, y = a, fill = value)) +
    geom_raster(interpolate = FALSE) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = midpoint,
      limits = fill_limits,
      oob = scales::squish,
      na.value = "grey70",
      name = legend_title,
      guide = guide_colorbar(
        title.position = "top",
        barheight = grid::unit(28, "mm"),
        barwidth = grid::unit(3, "mm"),
        ticks.colour = "black",
        frame.colour = "black"
      )
    ) +
    scale_x_continuous(breaks = brks$x, expand = c(0, 0)) +
    scale_y_continuous(breaks = brks$y, expand = c(0, 0)) +
    coord_cartesian(expand = FALSE) +
    labs(
      x = "c",
      y = "a",
      title = panel_title
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey92", linewidth = 0.25),
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 9, colour = "black"),
      legend.position = "right",
      legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 8),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0),
      plot.margin = margin(3, 3, 3, 3)
    )
}

plot_irrep_pair_publication <- function(my_res_GLASSO,
                                        my_res_PCGLASSO,
                                        all_a,
                                        all_c,
                                        which_graph,
                                        glasso_limits = NULL,
                                        pcglasso_limits = NULL,
                                        midpoint = 1,
                                        out_file = NULL,
                                        width = 6.2,
                                        height = 5.2,
                                        dpi = 300,
                                        a_every = 40) {
  plot_c_breaks <- list(
    (1:9) / 2 - 2.5,
    (1:9) / 8 - 0.625
  )[[which_graph]]

  if (is.null(glasso_limits)) {
    glasso_limits <- c(0, max(my_res_GLASSO, na.rm = TRUE))
  }
  if (is.null(pcglasso_limits)) {
    pcglasso_limits <- c(0, max(my_res_PCGLASSO, na.rm = TRUE))
  }

  p1 <- plot_irrep_publication(
    mat = my_res_GLASSO,
    all_a = all_a,
    all_c = all_c,
    panel_title = "GLASSO",
    fill_limits = glasso_limits,
    midpoint = midpoint,
    c_breaks = plot_c_breaks,
    legend_title = "IRR",
    a_every = a_every
  )

  p2 <- plot_irrep_publication(
    mat = my_res_PCGLASSO,
    all_a = all_a,
    all_c = all_c,
    panel_title = "PCGLASSO",
    fill_limits = pcglasso_limits,
    midpoint = midpoint,
    c_breaks = plot_c_breaks,
    legend_title = "IRR",
    a_every = a_every
  )

  p_both <- gridExtra::arrangeGrob(
    p1, p2,
    ncol = 1,
    heights = c(1, 1)
  )

  if (!is.null(out_file)) {
    ggplot2::ggsave(
      filename = out_file,
      plot = p_both,
      width = width,
      height = height,
      units = "in",
      dpi = dpi,
      bg = "white"
    )
  }

  p_both
}

common_limits <- c(0, max(c(my_res_GLASSO, my_res_PCGLASSO), na.rm = TRUE))

p_both <- plot_irrep_pair_publication(
  my_res_GLASSO = my_res_GLASSO,
  my_res_PCGLASSO = my_res_PCGLASSO,
  all_a = all_a,
  all_c = all_c,
  which_graph = which_graph,
  glasso_limits = common_limits,
  pcglasso_limits = common_limits,
  midpoint = 1,
  a_every = 40
)

grid::grid.newpage()
grid::grid.draw(p_both)

# ggsave(
#   paste0(
#     "outputs/Figure_1_irrep_on_", Ktext, ".png"
#   ), p_both,
#   width = 6.2, height = 5.2, units = "in", dpi = 300, bg = "white"
# )

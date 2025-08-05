# 2 minutes on 7 cores of Apple M2

library(pcglassoFast)
library(parallel)

library(reshape2)
library(ggplot2)
library(gridExtra)
library(dplyr)


which_graph <- 2 # 1 -> "KAR", 2 -> "hub"

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
  (a > 0 && a * b > (p - 1) * c * c)
}

KAR <- function(p, a, b, c) {
  K <- matrix(0, nrow = p, ncol = p)
  diag(K) <- b
  K[1, 1] <- a
  K[p, p] <- a

  for (i in 2:p) {
    K[i-1, i] <- K[i, i-1] <- c
  }

  K
}
is_KAR_positive_definite <- function(p, a, b, c) {
  K <- KAR(p, a, b, c)
  min(eigen(K, TRUE, TRUE)$values) > 0.00001
}


Kgenerator <- c(
  KAR, Khub
)[[which_graph]]
is_K_positive_definite <- c(
  is_KAR_positive_definite, is_Khub_positive_definite
)[[which_graph]]
Ktext <- c(
  "KAR", "hub"
)[which_graph]


single_irrep_value <- function(p, a, b, c) {
  if (!is_K_positive_definite(p, a, b, c)) {
    return(NA)
  }

  K <- Kgenerator(p, a, b, c)

  c(irrepGLASSO(K), irrepPCGLASSO(K))
}


all_a <- list(
  (1:200) * 0.05,
  (1:200) * 0.5
)[[which_graph]]
all_c <- list(
  (-63:63) / 120,
  (-120:120) / 60
)[[which_graph]]

my_res_GLASSO <- matrix(nrow = length(all_c), ncol = length(all_a))
my_res_PCGLASSO <- matrix(nrow = length(all_c), ncol = length(all_a))
for (i_for_a in 1:length(all_a)) {
  print(paste0(i_for_a, " of ", length(all_a)))

  my_results <- if (num_of_cores == 1) {
    lapply(all_c, function(c){
      single_irrep_value(p, all_a[i_for_a], b, c)
    })
  } else {
    mclapply(all_c, function(c){
      single_irrep_value(p, all_a[i_for_a], b, c)
    }, mc.cores = num_of_cores)
  }


  for (i_for_c in 1:length(all_c)) {
    my_res <- my_results[[i_for_c]]

    my_res_GLASSO[i_for_c, i_for_a] <- my_res[1]
    my_res_PCGLASSO[i_for_c, i_for_a] <- my_res[2]
  }
}



plot_a_breaks_scalling <- 20
plot_c_breaks <- list(
  (1:9)/8 - 0.625,
  (1:9)/2 - 2.5
)[[which_graph]]

mat <- my_res_GLASSO

mat_df <- melt(mat)
colnames(mat_df) <- c("x", "y", "value")
mat_df$x <- all_c[mat_df$x]
mat_df$y <- all_a[mat_df$y]

max_IRR <- c(2, 4)[which_graph]

p1 <- ggplot(mat_df, aes(x = x, y = y, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "green",
    mid = "white",
    high = "red",
    midpoint = 1,
    limits = c(0, max_IRR)
  ) +
  labs(x = "c", y = "a", fill = "IRR GLASSO") +
  theme_minimal() +
  scale_y_continuous(breaks = c(all_a[1] - (all_a[2] - all_a[1]), all_a[plot_a_breaks_scalling * (1:floor(length(all_a) / plot_a_breaks_scalling))])) +
  scale_x_continuous(breaks = plot_c_breaks)

mat_df$value <- case_when(
  is.na(mat_df$value) ~ "NA",
  mat_df$value >= 1 ~ "Not satisfied",
  TRUE ~ "Satisfied"
)
p1_two_colors <- ggplot(mat_df, aes(x = x, y = y, fill = value)) +
  geom_tile() +
  scale_fill_manual(values = c(
    "NA" = "dimgrey", "Not satisfied" = "red", "Satisfied" = "green"
  )) +
  labs(x = "c", y = "a", fill = "IRR GLASSO") +
  theme_minimal() +
  scale_y_continuous(breaks = c(all_a[1] - (all_a[2] - all_a[1]), all_a[plot_a_breaks_scalling * (1:floor(length(all_a) / plot_a_breaks_scalling))])) +
  scale_x_continuous(breaks = plot_c_breaks)


mat <- my_res_PCGLASSO
mat_df <- melt(mat)
colnames(mat_df) <- c("x", "y", "value")
mat_df$x <- all_c[mat_df$x]
mat_df$y <- all_a[mat_df$y]

p2 <- ggplot(mat_df, aes(x = x, y = y, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "green",
    mid = "white",
    high = "red",
    midpoint = 1,
    limits = c(0, max_IRR)
  ) +
  labs(x = "c", y = "a", fill = "IRR PCGLASSO") +
  theme_minimal() +
  scale_y_continuous(breaks = c(all_a[1] - (all_a[2] - all_a[1]), all_a[plot_a_breaks_scalling * (1:floor(length(all_a) / plot_a_breaks_scalling))])) +
  scale_x_continuous(breaks = plot_c_breaks)

mat_df$value <- case_when(
  is.na(mat_df$value) ~ "NA",
  mat_df$value >= 1 ~ "Not satisfied",
  TRUE ~ "Satisfied"
)
p2_two_colors <- ggplot(mat_df, aes(x = x, y = y, fill = value)) +
  geom_tile() +
  scale_fill_manual(values = c(
    "NA" = "dimgrey", "Not satisfied" = "red", "Satisfied" = "green"
  )) +
  labs(x = "c", y = "a", fill = "IRR PCGLASSO") +
  theme_minimal() +
  scale_y_continuous(breaks = c(all_a[1] - (all_a[2] - all_a[1]), all_a[plot_a_breaks_scalling * (1:floor(length(all_a) / plot_a_breaks_scalling))])) +
  scale_x_continuous(breaks = plot_c_breaks)

p_both <- grid.arrange(p1, p2, ncol = 1)
print(p_both)

ggsave(
  paste0(
    "raw_experiments/irrepresentability/plots/irrep_on_", Ktext, ".png"
  ), p_both,
  width = 7, height = 4, units = "in"
)


p_both_two_colors <- grid.arrange(p1_two_colors, p2_two_colors, ncol = 1)
print(p_both_two_colors)

ggsave(
  paste0(
    "raw_experiments/irrepresentability/plots/irrep_on_", Ktext, "_two_colors.png"
  ), p_both_two_colors,
  width = 7, height = 4, units = "in"
)

#setwd("./raw_experiments/genetic/")
source("../estimation_function.R")

gamma <- 0.5
n.lambda <- 100
lambda_min_ratio = 0.3
data(X.Sanger)
Sigma.Sanger <- cov(X.Sanger)
n <- dim(X)[1]
PC.est <- cov2cor(solve(Sigma.Sanger))

alpha.pcest<-  get_alpha(PC.est,0.5)

lam_max <- max(abs(Sigma.Sanger - diag(diag(Sigma.Sanger))))
lam_min <- lambda_min_ratio * lam_max
lambdas <- exp(seq(log(lam_max), log(lam_min), length.out = n.lambda))

glasso.res   <- estimator_glasso(Sigma.Sanger,n, lambdas,gamma = gamma)
Optim.glasso <- get_optimal_matrix(glasso.res$path, glasso.res$path.loss)


space.res <-    estimator_space(Sigma.Sanger,n, lambdas, X,gamma = gamma)
Optim.space <- get_optimal_matrix(space.res$path, space.res$path.loss)

lam_max <- max(abs(cov2cor(Sigma.Sanger) - diag(diag(cov2cor(Sigma.Sanger)))))
lam_min <- lambda_min_ratio * lam_max
lambdas_corr <- seq(lam_max, lam_min, length.out = n.lambda)

pcglasso.res <- estimator_pcglasso(Sigma.Sanger, n, lambdas_corr, alpha_grid = 0, gamma = gamma)
Optim.pcglasso <- get_optimal_matrix(pcglasso.res$path, pcglasso.res$path.loss)

Hcorrglasso.res <- estimator_hubcorglasso(Sigma.Sanger,n, 15*lambdas_corr, gamma = gamma)
Optim.Hcorrglasso <- get_optimal_matrix(Hcorrglasso.res$path, Hcorrglasso.res$path.loss)

corglasso.res <- estimator_corglasso(Sigma.Sanger,n, lambdas_corr,gamma = gamma)
Optim.corglasso <- get_optimal_matrix(corglasso.res$path, corglasso.res$path.loss)


alpha_grid <- sort(unique(c(
  seq(-0.1, -0.01, length.out = 6),
  0)))
pcglasso.res <- estimator_pcglasso(Sigma.Sanger,n, lambdas_corr,alpha_grid,gamma,max_edge_fraction=0.7)

# Flatten all BIC values with corresponding alpha and lambda
bic_df <- do.call(rbind, lapply(names(pcglasso.res$path.loss), function(a) {
  loss <- pcglasso.res$path.loss[[a]]
  lambda <- pcglasso.res$path.all[[a]]$lambda
  data.frame(
    alpha = as.numeric(a),
    lambda = lambda,
    BIC = loss$BIC_gamma
  )
}))

# Find the row with minimum BIC
bic_min <- bic_df[which.min(bic_df$BIC), ]

# Extract the optimal alpha
optimal_alpha <- bic_min$alpha

library(ggplot2)

# Create tidy data frames for each method
df_glasso <- data.frame(
  Edges  = glasso.res$path.loss$nEdges,
  BIC    = glasso.res$path.loss$BIC_gamma,
  Method = "GLasso"
)

df_corglasso <- data.frame(
  Edges  = corglasso.res$path.loss$nEdges,
  BIC    = corglasso.res$path.loss$BIC_gamma,
  Method = "Cor-GLasso"
)

df_space <- data.frame(
  Edges  = space.res$path.loss$nEdges,
  BIC    = space.res$path.loss$BIC_gamma,
  Method = "SPACE"
)
df_pcglasso <- data.frame(
  Edges  = pcglasso.res$path.loss[[as.character(0)]]$nEdges,
  BIC    = pcglasso.res$path.loss[[as.character(0)]]$BIC_gamma,
  Method = "PC-GLasso"
)
df_pcglasso_alpha <- data.frame(
  Edges  = pcglasso.res$path.loss[[as.character(optimal_alpha)]]$nEdges,
  BIC    = pcglasso.res$path.loss[[as.character(optimal_alpha)]]$BIC_gamma,
  Method = paste("PC-GLasso"," alpha=",round(optimal_alpha,3),sep="")
)
df_hubglasso <- data.frame(
  Edges  = Hcorrglasso.res$path.loss$nEdges,
  BIC    = Hcorrglasso.res$path.loss$BIC_gamma,
  Method = "Hub-Glasso"
)
# Combine all
df_all <- rbind(df_glasso, df_corglasso, df_space, df_pcglasso)#,df_pcglasso_alpha)

# Plot
# Dynamic method label
#label_pcglasso_opt <- paste("PC-GLasso", " alpha=", round(optimal_alpha, 3), sep = "")


colors_named <- c(
  "GLasso"     = "#1b9e77",
  "Cor-GLasso" = "#7570b3",
  "SPACE"      = "#e7298a",
  "PC-GLasso"  = "#d95f02")
  #"Hub-Glasso" = "#e6ab02"
  # Use the *value* of the label as name
#  setNames("#e6ab02", label_pcglasso_opt)
#)

fig <- ggplot(df_all, aes(x = Edges, y = BIC, color = Method)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = colors_named) +
  labs(x = "#Edges", y = "BIC") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    legend.box = "horizontal"
  ) +
  coord_cartesian(xlim = c(0, 1000), ylim = c(28000, 35000))

print(fig)
ggsave(
  "BIC_fig.png",
  plot = fig, width = 7, height = 4
)


library(patchwork)  # or library(gridExtra)

p_glasso   <- make_plot_matrix(Optim.glasso$Theta_opt, "GLasso")
p_corglasso     <- make_plot_matrix(Optim.corglasso$Theta_opt, "Cor-GLasso")
p_pcglasso    <- make_plot_matrix(Optim.pcglasso$Theta_opt, "PCGLasso")
p_space   <- make_plot_matrix(Optim.space$Theta_opt, "SPACE")

# Combine all 4 plots in 2×2 layout
fig <- ((p_glasso | p_corglasso) /
  (p_pcglasso | p_space))
print(fig)
ggsave(
  "matrices.png",
  plot = fig, width = 7, height = 4
)


Optim.nEdge.space <- get_optimal_matrix(space.res$path, space.res$path.loss, max_edges=Optim.pcglasso$nEdges)
Optim.nEdge.corrglasso <- get_optimal_matrix(corglasso.res$path, corglasso.res$path.loss, max_edges=Optim.pcglasso$nEdges)
Optim.nEdge.glasso <- get_optimal_matrix(glasso.res$path, glasso.res$path.loss, max_edges=Optim.pcglasso$nEdges)

p2_glasso   <- make_plot_matrix(Optim.nEdge.glasso$Theta_opt, "GLasso")
p2_corglasso     <- make_plot_matrix(Optim.nEdge.corrglasso$Theta_opt, "Cor-GLasso")
p2_pcglasso    <- make_plot_matrix(Optim.pcglasso$Theta_opt, "PCGLasso")
p2_space   <- make_plot_matrix(Optim.nEdge.space$Theta_opt, "SPACE")

fig <- ((p2_glasso | p2_corglasso) /
          (p2_pcglasso | p2_space))
print(fig)
ggsave(
  "matrices2.png",
  plot = fig, width = 7, height = 4
)


# Helper to make a single sorted-alpha plot
make_alpha_plot <- function(alpha, title, ylims = NULL) {
  a <- alpha[is.finite(alpha)]
  stopifnot(length(a) > 0)
  a <- sort(a)
  df <- data.frame(idx = seq_along(a), alpha = a)

  p <- ggplot(df, aes(x = idx, y = alpha)) +
    geom_point(size = 0.9, alpha = 0.9) +
    geom_line(linewidth = 0.3) +
    labs(title = title, x = "Order statistic index", y = expression(alpha)) +
    theme_minimal(base_size = 11) +
    theme(plot.title.position = "plot", panel.grid.minor = element_blank())

  if (!is.null(ylims)) p <- p + scale_y_continuous(limits = ylims)
  p
}


# Build a grid of alpha plots from a named list
make_alpha_grid <- function(alpha_list, ncol = 2, common_y = TRUE) {
  stopifnot(length(alpha_list) > 0)
  if (is.null(names(alpha_list))) {
    names(alpha_list) <- paste0("Series ", seq_along(alpha_list))
  }
  ylims <- if (common_y) range(unlist(alpha_list), finite = TRUE, na.rm = TRUE) else NULL

  plots <- lapply(names(alpha_list), function(nm) {
    make_alpha_plot(alpha_list[[nm]], nm, ylims)
  })
  # Arrange as a grid with patchwork
  wrap_plots(plots, ncol = ncol)
}

# --- your original set (if you still want it) ------------------------------
alphas <- list(
  "GLasso"     = Optim.nEdge.glasso$alpha,
  "Cor-GLasso" = Optim.nEdge.corrglasso$alpha,
  "PCGLasso"   = Optim.pcglasso$alpha,
  "SPACE"      = Optim.nEdge.space$alpha
)
fig_alpha <- make_alpha_grid(alphas, ncol = 2, common_y = TRUE)
print(fig_alpha)
ggsave("alphas_grid.png", plot = fig_alpha, width = 7, height = 4, dpi = 300)

# --- your new request: alphas.optim ---------------------------------------
alphas.optim <- list(
  "GLasso"     = Optim.glasso$alpha,
  "Cor-GLasso" = Optim.corglasso$alpha,
  "PCGLasso"   = Optim.pcglasso$alpha,
  "SPACE"      = Optim.space$alpha
)

fig_alpha_opt <- make_alpha_grid(alphas.optim, ncol = 2, common_y = TRUE)
print(fig_alpha_opt)
ggsave("alphas_optim_grid.png", plot = fig_alpha_opt, width = 7, height = 4, dpi = 300)

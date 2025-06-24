source("../estimation_function.R")

gamma <- 0.3
n.lambda <- 100
lambda.min.ratio = 0.03
data(Sanger)
X <- Sanger
Sigma.Sanger <- cov(X)
n <- dim(X)[1]

lam_max <- 2*max(abs(Sigma.Sanger - diag(diag(Sigma.Sanger))))
lam_min <- lambda.min.ratio * lam_max
lambdas <- exp(seq(log(lam_max), log(lam_min), length.out = n.lambda))


lam_max <- 2*max(abs(cov2cor(Sigma.Sanger) - diag(diag(cov2cor(Sigma.Sanger)))))
lam_min <- lambda.min.ratio * lam_max
lambdas_corr <- exp(seq(log(lam_max), log(lam_min), length.out = n.lambda))

pcglasso.res <- estimator_pcglasso(Sigma.Sanger,n, lambdas,alpha.grid=0, gamma = gamma)
glasso.res   <- estimator_glasso(Sigma.Sanger,n, lambdas,gamma = gamma)
corglasso.res <- estimator_corglasso(Sigma.Sanger,n, lambdas_corr,gamma = gamma)
space.res <-    estimator_space(Sigma.Sanger,n, lambdas, X,gamma = gamma)

alpha.grid <- sort(unique(c(
  seq(-0.1, -0.01, length.out = 6),
  0)))
pcglasso.res <- estimator_pcglasso(Sigma.Sanger,n, lambdas,alpha.grid,gamma,max.edge.fraction=0.7)

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

# Combine all
df_all <- rbind(df_glasso, df_corglasso, df_space, df_pcglasso,df_pcglasso_alpha)

# Plot
# Dynamic method label
label_pcglasso_opt <- paste("PC-GLasso", " alpha=", round(optimal_alpha, 3), sep = "")


colors_named <- c(
  "GLasso"     = "#1b9e77",
  "Cor-GLasso" = "#7570b3",
  "SPACE"      = "#e7298a",
  "PC-GLasso"  = "#d95f02",
  # Use the *value* of the label as name
  setNames("#e6ab02", label_pcglasso_opt)
)

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
  coord_cartesian(xlim = c(0, 2000), ylim = c(28000, 40000))

print(fig)
ggsave(
  "BIC_fig.png",
  plot = fig, width = 7, height = 4
)


library(patchwork)  # or library(gridExtra)

# GLasso
idx_glasso <- which.min(glasso.res$path.loss$BIC_gamma)
Q_glasso   <- glasso.res$path[,,idx_glasso]
p_glasso   <- make_plot_matrix(Q_glasso, "GLasso")

# Cor-GLasso
idx_corg   <- which.min(corglasso.res$path.loss$BIC_gamma)
Q_corg     <- corglasso.res$path[,,idx_corg]
p_corg     <- make_plot_matrix(Q_corg, "Cor-GLasso")

idx_pcglasso_opt <- which.min(pcglasso.res$path.loss[[as.character(optimal_alpha)]]$BIC_gamma)
Q_pcglasso_opt   <- pcglasso.res$path[[as.character(optimal_alpha)]][,,idx_pcglasso_opt]
label_pcglasso_opt <- paste("PC-GLasso", "\u03B1=", round(optimal_alpha, 3))
p_pcglasso <- make_plot_matrix(Q_pcglasso_opt, label_pcglasso_opt)

idx_space <- which.min(space.res$path.loss$BIC_gamma)
Q_space   <- space.res$path[,,idx_space]
p_space   <- make_plot_matrix(Q_space, "SPACE")

# Combine all 4 plots in 2×2 layout
fig <- ((p_glasso | p_corg) /
  (p_pcglasso | p_space))
print(fig)
ggsave(
  "matrices.png",
  plot = fig, width = 7, height = 4
)

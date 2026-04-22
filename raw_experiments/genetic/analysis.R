#setwd("./raw_experiments/genetic/")
source("../estimation_function.R")
#source("./raw_experiments/estimation_function.R")

gamma <- 0.5
n.lambda <- 100
lambda_min_ratio = 0.3
data(X.Sanger)
Sigma.Sanger <- cov(X.Sanger)
n <- dim(X.Sanger)[1]
PC.est <- cov2cor(solve(Sigma.Sanger))

alpha.pcest<-  get_alpha(PC.est,0.5)

lam_max <- max(abs(Sigma.Sanger - diag(diag(Sigma.Sanger))))
lam_min <- lambda_min_ratio * lam_max
lambdas <- exp(seq(log(lam_max), log(lam_min), length.out = n.lambda))

glasso.res   <- estimator_glasso(Sigma.Sanger,n, lambdas,gamma = gamma)
Optim.glasso <- get_optimal_matrix(glasso.res$path, glasso.res$path.loss)


space.res <-    estimator_space(Sigma.Sanger,n, lambdas, X.Sanger,gamma = gamma)
Optim.space <- get_optimal_matrix(space.res$path, space.res$path.loss)

lam_max <- max(abs(cov2cor(Sigma.Sanger) - diag(diag(cov2cor(Sigma.Sanger)))))
lam_min <- lambda_min_ratio * lam_max
lambdas_corr <- seq(lam_max, lam_min, length.out = n.lambda)

pcglasso.res <- estimator_pcglasso(Sigma.Sanger, n, lambdas_corr, alpha_grid = 0, gamma = gamma)
Optim.pcglasso <- get_optimal_matrix(pcglasso.res$path, pcglasso.res$path.loss)

#pcglasso_cpp.res <- estimator_pcglasso_cpp(Sigma.Sanger, n, lambdas_corr, alpha_grid = 0, gamma = gamma)
#Optim.pcglasso_cpp <- get_optimal_matrix(pcglasso_cpp.res$path, pcglasso_cpp.res$path.loss)

Hcorrglasso.res <- estimator_hubcorglasso(Sigma.Sanger,n, 15*lambdas_corr, gamma = gamma)
Optim.Hcorrglasso <- get_optimal_matrix(Hcorrglasso.res$path, Hcorrglasso.res$path.loss)

corglasso.res <- estimator_corglasso(Sigma.Sanger,n, lambdas_corr,gamma = gamma)
Optim.corglasso <- get_optimal_matrix(corglasso.res$path, corglasso.res$path.loss)


#alpha_grid <- sort(unique(c(
#  seq(-0.1, -0.01, length.out = 6),
#  0)))
#pcglasso.res <- estimator_pcglasso(Sigma.Sanger,n, lambdas_corr,alpha_grid,gamma,max_edge_fraction=0.7)

# Flatten all BIC values with corresponding alpha and lambda


# Find the row with minimum BIC

library(ggplot2)

# Create tidy data frames for each method
df_glasso <- data.frame(
  Edges  = glasso.res$path.loss$nEdges,
  BIC    = glasso.res$path.loss$BIC_gamma,
  Method = "GLASSO"
)

df_corglasso <- data.frame(
  Edges  = corglasso.res$path.loss$nEdges,
  BIC    = corglasso.res$path.loss$BIC_gamma,
  Method = "Cor-GLASSO"
)

df_space <- data.frame(
  Edges  = space.res$path.loss$nEdges,
  BIC    = space.res$path.loss$BIC_gamma,
  Method = "SPACE"
)
df_pcglasso <- data.frame(
  Edges  = pcglasso.res$path.loss$nEdges,
  BIC    = pcglasso.res$path.loss$BIC_gamma,
  Method = "PCGLASSO"
)


# Combine all
df_all <- rbind(df_glasso, df_corglasso, df_space, df_pcglasso)#,df_pcglasso_alpha)

# Plot
# Dynamic method label
#label_pcglasso_opt <- paste("PC-GLasso", " alpha=", round(optimal_alpha, 3), sep = "")


colors_named <- c(
  "GLASSO"     = "#1b9e77",
  "Cor-GLASSO" = "#7570b3",
  "SPACE"      = "#e7298a",
  "PCGLASSO"  = "#d95f02")
  #"Hub-Glasso" = "#e6ab02"
  # Use the *value* of the label as name
#  setNames("#e6ab02", label_pcglasso_opt)
#)

fig <- ggplot(df_all, aes(x = Edges, y = BIC, colour = Method)) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = colors_named) +
  scale_x_continuous(
    name = "Number of edges",
    breaks = seq(0, 300, 50),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_y_continuous(
    name = expression(EBIC(gamma == 0.5)),
    breaks = seq(28000, 35000, 1000),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  coord_cartesian(xlim = c(0, 300), ylim = c(29000, 34000)) +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE)) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.key.width = grid::unit(1.3, "lines"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(linewidth = 0.4),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

print(fig)
 ggsave(
   "BIC_fig.png",
   plot = fig, width = 7, height = 4
 )


library(patchwork)  # or library(gridExtra)

p_glasso   <- make_plot_matrix(Optim.glasso$Theta_opt, "GLASSO", x_lab="", y_lab="")
p_corglasso     <- make_plot_matrix(Optim.corglasso$Theta_opt, "Cor-GLASSO", x_lab="", y_lab="")
p_pcglasso    <- make_plot_matrix(Optim.pcglasso$Theta_opt, "PCGLASSO", x_lab="", y_lab="")
p_space   <- make_plot_matrix(Optim.space$Theta_opt, "SPACE", x_lab="", y_lab="")

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

p2_glasso   <- make_plot_matrix(Optim.nEdge.glasso$Theta_opt, "GLASSO", x_lab="", y_lab="")
p2_corglasso     <- make_plot_matrix(Optim.nEdge.corrglasso$Theta_opt, "Cor-GLASSO", x_lab="", y_lab="")
p2_pcglasso    <- make_plot_matrix(Optim.pcglasso$Theta_opt, "PCGLASSO", x_lab="", y_lab="")
p2_space   <- make_plot_matrix(Optim.nEdge.space$Theta_opt, "SPACE", x_lab="", y_lab="")

fig <- ((p2_glasso | p2_corglasso) /
          (p2_pcglasso | p2_space))
print(fig)
ggsave(
  "matrices2.png",
  plot = fig, width = 7, height = 4
)



# --- your original set (if you still want it) ------------------------------
alphas <- list(
  "GLASSO"     = Optim.nEdge.glasso$alpha,
  "Cor-GLASSO" = Optim.nEdge.corrglasso$alpha,
  "PCGLASSO"   = Optim.pcglasso$alpha,
  "SPACE"      = Optim.nEdge.space$alpha
)
fig_alpha <- make_alpha_grid(alphas, ncol = 2, common_y = TRUE)
print(fig_alpha)
ggsave("alphas_grid.png", plot = fig_alpha, width = 7, height = 4, dpi = 300)

# --- your new request: alphas.optim ---------------------------------------
alphas.optim <- list(
  "GLASSO"     = Optim.glasso$alpha,
  "Cor-GLASSO" = Optim.corglasso$alpha,
  "PCGLASSO"   = Optim.pcglasso$alpha,
  "SPACE"      = Optim.space$alpha
)

fig_alpha_opt <- make_alpha_grid(alphas.optim, ncol = 2, common_y = TRUE)
print(fig_alpha_opt)
ggsave("alphas_optim_grid.png", plot = fig_alpha_opt, width = 7, height = 4, dpi = 300)


theta.pcglasso <- Optim.pcglasso$Theta_opt
theta.glasso <- Optim.corglasso$Theta_opt
library(igraph)
par(mfrow = c(1, 1))
glasso.corr.ind <- which(theta.pcglasso[124, -124] != 0)
res.pcglasso.ind <- (which(theta.glasso[124, -124] != 0))
joint.ind <- intersect(glasso.corr.ind, res.pcglasso.ind)
glasso.corr.ind <- setdiff(glasso.corr.ind, joint.ind)
res.pcglasso.ind <- setdiff(res.pcglasso.ind, joint.ind)
sub.graph <- c(joint.ind, res.pcglasso.ind, glasso.corr.ind, 124)
sub.graph.col <- c(
  "lightgreen",
  rep(
    "yellow",
    length(res.pcglasso.ind)
  ),
  rep("lightblue", length(glasso.corr.ind)), "red"
)
ppcglasso <- cov2cor(theta.pcglasso[sub.graph, sub.graph])
vertex.label <- sub.graph
vertex.label[length(vertex.label)] <- "CCT8"
diag(ppcglasso) <- 0
infered_graph.pcgl <- graph_from_adjacency_matrix(ppcglasso != 0, mode = "undirected", weighted = FALSE)
V(infered_graph.pcgl)$color <- sub.graph.col
V(infered_graph.pcgl)$size <- 22
pdf("PCGLASSO_sub.pdf", width = 7, height = 7)
lay <- layout_with_fr(infered_graph.pcgl)
plot(infered_graph.pcgl, layout = lay, vertex.label = vertex.label, main = "PCGLASSO")
dev.off()


glasso_cor <- cov2cor(theta.glasso[sub.graph, sub.graph])
diag(glasso_cor) <- 0
infered_graph.gl_cor <- graph_from_adjacency_matrix(glasso_cor != 0, mode = "undirected", weighted = FALSE)
V(infered_graph.gl_cor)$color <- sub.graph.col
V(infered_graph.gl_cor)$size <- 22
pdf("PCGLASSO_corr_sub.pdf", width = 7, height = 7)
plot(infered_graph.gl_cor, layout = lay, vertex.label = vertex.label, main = "Cor-GLASSO")
dev.off()



# Diagonal analysis
D_corr_sd <- diag(Optim.corglasso_sd$Theta_opt)
D_corr <- diag(Optim.corglasso$Theta_opt)
D_pcglasso <- diag(Optim.pcglasso$Theta_opt)


df <- data.frame(
  id          = seq_along(D_corr),        # index for optional labels
  D_corr      = as.numeric(D_corr),
  D_pcglasso  = as.numeric(D_pcglasso),
  D_corr_sd   = as.numeric(D_corr_sd)
)

# remove any non-finite rows (safer plotting)
df <- subset(df, is.finite(D_corr) & is.finite(D_pcglasso) & is.finite(D_corr_sd))



p <- ggplot(df, aes(D_corr, D_pcglasso, label = id)) +
  geom_text(size = 3) +                         # <- labels instead of points
  geom_abline(slope = 1, intercept = 0) +       # 1–1 line
  labs(title = "Diagonal entries",
       x = "Cor-GLASSO",
       y = "PCGLASSO") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("diag_PCGlasso_Corr.pdf",p,width=5,height=5)

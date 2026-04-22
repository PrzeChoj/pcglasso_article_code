source("./experiments/utils.R")

gamma <- 0.5
n.lambda <- 100
lambda_min_ratio <- 0.3
data(X.Sanger, package = "pcglassoFast")
Sigma.Sanger <- cov(X.Sanger)
n <- dim(X.Sanger)[1]

lam_max <- max(abs(Sigma.Sanger - diag(diag(Sigma.Sanger))))
lam_min <- lambda_min_ratio * lam_max
lambdas <- exp(seq(log(lam_max), log(lam_min), length.out = n.lambda))

glasso.res <- estimator_glasso(Sigma.Sanger, n, lambdas, gamma = gamma)
space.res <- estimator_space(Sigma.Sanger, n, lambdas, X.Sanger, gamma = gamma)

lam_max <- max(abs(cov2cor(Sigma.Sanger) - diag(diag(cov2cor(Sigma.Sanger)))))
lam_min <- lambda_min_ratio * lam_max
lambdas_corr <- seq(lam_max, lam_min, length.out = n.lambda)

pcglasso.res <- estimator_pcglasso(Sigma.Sanger, n, lambdas_corr, alpha_grid = 0, gamma = gamma)
Optim.pcglasso <- get_optimal_matrix(pcglasso.res$path, pcglasso.res$path.loss)

corglasso.res <- estimator_corglasso(Sigma.Sanger, n, lambdas_corr, gamma = gamma)
Optim.corglasso <- get_optimal_matrix(corglasso.res$path, corglasso.res$path.loss)

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
df_all <- rbind(df_glasso, df_corglasso, df_space, df_pcglasso)

colors_named <- c(
  "GLASSO" = "#1b9e77",
  "Cor-GLASSO" = "#7570b3",
  "SPACE" = "#e7298a",
  "PCGLASSO" = "#d95f02"
)

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
  "./outputs/Figure_2_BIC_fig.png",
  plot = fig, width = 7, height = 4
)


library(patchwork) # or library(gridExtra)

Optim.nEdge.space <- get_optimal_matrix(space.res$path, space.res$path.loss, max_edges = Optim.pcglasso$nEdges)
Optim.nEdge.corrglasso <- get_optimal_matrix(corglasso.res$path, corglasso.res$path.loss, max_edges = Optim.pcglasso$nEdges)
Optim.nEdge.glasso <- get_optimal_matrix(glasso.res$path, glasso.res$path.loss, max_edges = Optim.pcglasso$nEdges)

p2_glasso <- make_plot_matrix(Optim.nEdge.glasso$Theta_opt, "GLASSO", x_lab = "", y_lab = "")
p2_corglasso <- make_plot_matrix(Optim.nEdge.corrglasso$Theta_opt, "Cor-GLASSO", x_lab = "", y_lab = "")
p2_pcglasso <- make_plot_matrix(Optim.pcglasso$Theta_opt, "PCGLASSO", x_lab = "", y_lab = "")
p2_space <- make_plot_matrix(Optim.nEdge.space$Theta_opt, "SPACE", x_lab = "", y_lab = "")

fig <- ((p2_glasso | p2_corglasso) /
          (p2_pcglasso | p2_space))
print(fig)
ggsave(
  "./outputs/Figure_3_matrices2.png",
  plot = fig, width = 7, height = 4
)




theta.pcglasso <- Optim.pcglasso$Theta_opt
theta.glasso <- Optim.corglasso$Theta_opt
library(igraph)
library(ggraph)
glasso.corr.ind <- which(theta.pcglasso[124, -124] != 0)
res.pcglasso.ind <- which(theta.glasso[124, -124] != 0)
joint.ind <- intersect(glasso.corr.ind, res.pcglasso.ind)
glasso.corr.ind <- setdiff(glasso.corr.ind, joint.ind)
res.pcglasso.ind <- setdiff(res.pcglasso.ind, joint.ind)
sub.graph <- c(joint.ind, res.pcglasso.ind, glasso.corr.ind, 124)
sub.graph.col <- c(
  "lightgreen",
  rep("yellow", length(res.pcglasso.ind)),
  rep("lightblue", length(glasso.corr.ind)),
  "red"
)
ppcglasso <- cov2cor(theta.pcglasso[sub.graph, sub.graph])
vertex.label <- sub.graph
vertex.label[length(vertex.label)] <- "CCT8"
diag(ppcglasso) <- 0
infered_graph.pcgl <- graph_from_adjacency_matrix(ppcglasso != 0, mode = "undirected", weighted = FALSE)
V(infered_graph.pcgl)$node_color <- sub.graph.col
V(infered_graph.pcgl)$label <- vertex.label
lay <- layout_with_fr(infered_graph.pcgl)

fig_pcgl_sub <- ggraph(infered_graph.pcgl, layout = "manual", x = lay[, 1], y = lay[, 2]) +
  geom_edge_link(color = "grey70") +
  geom_node_point(aes(color = node_color), size = 10) +
  geom_node_text(aes(label = label), size = 3, family = "sans") +
  scale_color_identity() +
  labs(title = "PCGLASSO") +
  theme_graph()

print(fig_pcgl_sub)
ggsave("./outputs/Figure_5_PCGLASSO_sub.png", plot = fig_pcgl_sub, width = 7, height = 7, dpi = 300)


glasso_cor <- cov2cor(theta.glasso[sub.graph, sub.graph])
diag(glasso_cor) <- 0
infered_graph.gl_cor <- graph_from_adjacency_matrix(glasso_cor != 0, mode = "undirected", weighted = FALSE)
V(infered_graph.gl_cor)$node_color <- sub.graph.col
V(infered_graph.gl_cor)$label <- vertex.label

fig_gl_cor_sub <- ggraph(infered_graph.gl_cor, layout = "manual", x = lay[, 1], y = lay[, 2]) +
  geom_edge_link(color = "grey70") +
  geom_node_point(aes(color = node_color), size = 10) +
  geom_node_text(aes(label = label), size = 3, family = "sans") +
  scale_color_identity() +
  labs(title = "Cor-GLASSO") +
  theme_graph()

print(fig_gl_cor_sub)
ggsave("./outputs/Figure_5_PCGLASSO_corr_sub.png", plot = fig_gl_cor_sub, width = 7, height = 7, dpi = 300)



# Diagonal analysis
D_corr <- diag(Optim.corglasso$Theta_opt)
D_pcglasso <- diag(Optim.pcglasso$Theta_opt)

df <- data.frame(
  id          = seq_along(D_corr),
  D_corr      = as.numeric(D_corr),
  D_pcglasso  = as.numeric(D_pcglasso)
)

# remove any non-finite rows (safer plotting)
df <- subset(df, is.finite(D_corr) & is.finite(D_pcglasso))



p <- ggplot(df, aes(D_corr, D_pcglasso, label = id)) +
  geom_text(size = 3) + # <- labels instead of points
  geom_abline(slope = 1, intercept = 0) + # 1–1 line
  labs(
    title = "Diagonal entries",
    x = "Cor-GLASSO",
    y = "PCGLASSO"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5))
print(p)
ggsave("./outputs/Figure_4_diag_PCGlasso_Corr.pdf", p, width = 5, height = 5)

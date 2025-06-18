source("./expermient/sangerdata/analysis.R")

library(ggplot2)
library(dplyr)

# Combine data
data <- bind_rows(
  data.frame(method = "glasso", nEdges = res.loss.diag.glasso$nEdges, BIC_gamma = res.loss.diag.glasso$BIC_gamma),
  data.frame(method = "diag_glasso", nEdges = res.loss.dlasso$nEdges, BIC_gamma = res.loss.dlasso$BIC_gamma),
  data.frame(method = "pcglasso", nEdges = res.loss.pcglasso$nEdges, BIC_gamma = res.loss.pcglasso$BIC_gamma),
  data.frame(method = "glasso_sd", nEdges = res.loss.glasso.corr$nEdges, BIC_gamma = res.loss.glasso.corr$BIC_gamma),
  data.frame(method = "diag_glasso_sd", nEdges = res.loss.glasso.corr_diag$nEdges, BIC_gamma = res.loss.glasso.corr_diag$BIC_gamma)
)

p1 <- ggplot(data, aes(x = nEdges, y = BIC_gamma, color = method, shape = as.factor(method))) +
  geom_point(size = 3) +
  scale_color_manual(values = c("red", "black", "red", "black", "blue")) +
  scale_shape_manual(values = c(2, 5, 1, 4, 3)) +
  labs(
    title = "Comparison of BIC_gamma Across Methods",
    x = "Number of Edges",
    y = "BIC_gamma",
    color = "Method",
    shape = "Shape"
  ) +
  xlim(0, 1000) +
  ylim(-90, -40) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

ggsave(
  "./expermient/sangerdata/img/BIC_fig.png",
  plot = p1, width = 7, height = 4
)

#####


m.pcglasso <- which.min(res.loss.pcglasso$BIC_gamma)
m.glasso <- which.min(res.loss.diag.glasso$BIC_gamma)
m.glasso.coor <- which.min(res.loss.glasso.corr_diag$BIC_gamma)

make_plot_matrix <- function(my_matrix, my_title) {
  matrix_data <- my_matrix != 0
  df_matrix <- as.data.frame(as.table(matrix_data))
  colnames(df_matrix) <- c("Row", "Column", "Value")

  df_matrix$Row <- as.numeric(df_matrix$Row)
  df_matrix$Column <- as.numeric(df_matrix$Column)
  df_matrix$Value <- as.numeric(df_matrix$Value)

  nnz <- sum(matrix_data)

  ggplot(df_matrix, aes(x = Column, y = Row, fill = Value)) +
    geom_tile(color = "white") +
    scale_fill_gradient(low = "white", high = "blue", name = "Non-Zero") +
    labs(
      title = paste(my_title, ", nnz =", nnz),
      x = NULL,
      y = NULL
    ) +
    scale_x_continuous(breaks = c(20, 40, 60, 80, 100, 120)) +
    scale_y_continuous(breaks = c(20, 40, 60, 80, 100, 120)) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.ticks = element_line(),
      legend.position = "none",
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    )
}

p2 <- make_plot_matrix(
  my_matrix = res.pcglasso[, , m.pcglasso], my_title = "PCGLASSO"
)

ggsave(
  "./expermient/sangerdata/img/image_pcglasso.png",
  plot = p2, width = 4, height = 4
)
#####

p3 <- make_plot_matrix(
  my_matrix = res.diag.glasso[, , m.glasso], my_title = "GLASSO"
)

ggsave(
  "./expermient/sangerdata/img/image_glasso.png",
  plot = p3, width = 4, height = 4
)

p4 <- make_plot_matrix(
  my_matrix = res.glasso.corr_diag[, , m.glasso.coor], my_title = "GLASSO corr"
)

ggsave(
  "./expermient/sangerdata/img/image_glasso_corr.png",
  plot = p4, width = 4, height = 4
)


Rpc <- rowSums(res.pcglasso[, , m.pcglasso] != 0)
cbind(Rpc[order(Rpc, decreasing = T)], order(Rpc, decreasing = T))

# cbind(PCGLASSO:::unit.diag.mat(res.pcglasso[,,m.pcglasso])[124,res.pcglasso[124,,m.pcglasso]!=0],ind[which(res.pcglasso[124,,m.pcglasso]!=0)])
cbind(PCGLASSO:::unit.diag.mat(res.glasso.corr[, , m.glasso.coor])[124, res.glasso.corr[124, , m.glasso.coor] != 0], ind[which(res.glasso.corr[124, , m.glasso.coor] != 0)])

library(igraph)
par(mfrow = c(1, 1))
glasso.corr.ind <- which(res.glasso.corr[124, -124, m.glasso.coor] != 0)
res.pcglasso.ind <- (which(res.pcglasso[124, -124, m.pcglasso] != 0))
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
ppcglasso <- cov2cor(res.pcglasso[sub.graph, sub.graph, m.pcglasso])
vertex.label <- sub.graph
vertex.label[length(vertex.label)] <- "CCT8"
diag(ppcglasso) <- 0
infered_graph.pcgl <- graph_from_adjacency_matrix(ppcglasso != 0, mode = "undirected", weighted = FALSE)
V(infered_graph.pcgl)$color <- sub.graph.col
V(infered_graph.pcgl)$size <- 22
pdf("./expermient/sangerdata/img/PCGLASSO_sub.pdf", width = 7, height = 7)
lay <- layout_with_fr(infered_graph.pcgl)
plot(infered_graph.pcgl, layout = lay, vertex.label = vertex.label, main = "PCGLASSO")
dev.off()

glasso_cor <- cov2cor(res.glasso.corr_diag[sub.graph, sub.graph, m.glasso.coor])
diag(glasso_cor) <- 0
infered_graph.gl_cor <- graph_from_adjacency_matrix(glasso_cor != 0, mode = "undirected", weighted = FALSE)
V(infered_graph.gl_cor)$color <- sub.graph.col
V(infered_graph.gl_cor)$size <- 22
pdf("./expermient/sangerdata/img/PCGLASSO_corr_sub.pdf", width = 7, height = 7)
plot(infered_graph.gl_cor, layout = lay, vertex.label = vertex.label, main = "GLASSO corr")
dev.off()


# diagonal

pcglasso_sub_corr.diag <- diag(res.pcglasso[sub.graph, sub.graph, m.pcglasso])
glasso_sub_corr.diag <- diag(res.glasso[sub.graph, sub.graph, m.glasso])

# Prepare data
data_diag <- data.frame(
  glasso_corr = glasso_sub_corr.diag,
  pcglasso_corr = pcglasso_sub_corr.diag,
  labels = c(sub.graph[-length(sub.graph)], "CCT8")
)

# Generate the plot
p5 <- ggplot(data_diag, aes(x = glasso_corr, y = pcglasso_corr)) +
  geom_point(color = "black", alpha = 0.8) + # Add points
  geom_text(aes(label = labels), hjust = 0.5, vjust = -0.5, size = 3) + # Add labels
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") + # Add diagonal line
  labs(
    title = "Diagonal entries",
    x = "glasso corr",
    y = "pcglasso"
  ) +
  ylim(0, 130) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

ggsave(
  "./expermient/sangerdata/img/diagonal_sub.png",
  plot = p5, width = 4, height = 4
)

names_genes <- c(colnames(xx)[ind], "CCT8")

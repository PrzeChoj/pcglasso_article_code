# Section 4.2 — Prostate cancer patients
#
# Fits GLASSO, Cor-GLASSO, and PCGLASSO to the prostate cancer RNA-seq dataset
# (Gómez et al., 2025). Compares EBIC paths and hub diagnostics based on
# partial-correlation scores vs. precision-matrix scores.
#
# Outputs (saved to outputs/):
#   Section_4_2_adj_matrix_pcglasso.png  — adjacency matrix with EEF2 highlighted
#   Appendix_B_matrices_optimal.png     — 2x2 grid of optimal precision matrices
#   Appendix_B_BIC_fig.png              — EBIC paths for all methods
#   Section_4_2_alpha2_pcglasso.png      — sorted alpha^2(K-hat)
#   Section_4_2_alpha1_pcglasso.png      — sorted alpha^1(R-hat)
#   Section_4_2_d_pcglasso.png           — sorted D(K-hat) diagonal entries
#   Appendix_B_D_rows_cancer.png        — normalized precision rows for top-diag genes

library(patchwork)
library(ggplot2)
library(ggrepel)
source("./experiments/utils.R")

gamma    <- 0.5
n.lambda <- 50

df_reduced <- read.csv("./experiments/Section_4_2/312_df_reduced.csv")
df_reduced <- df_reduced[, -1]
df_genes   <- df_reduced[, 7 + 1:200]
gene.names <- colnames(df_genes)

Sigma.est <- cov(df_genes)
Corr.est  <- cov2cor(Sigma.est)
n         <- nrow(df_genes)
Prec.est  <- solve(Sigma.est)
PC.est    <- cov2cor(Prec.est)

alpha_2          <- apply(Prec.est, 2, function(x) sum(abs(x)^2) - 1)
alpha            <- apply(PC.est,   2, function(x) sum(abs(x)) - 1)
p_alpha_emperical <- make_alpha_plot(alpha, "Empirical estimated alpha")


######################
# GLASSO
######################
lam_max        <- 0.2  * max(abs(Sigma.est - diag(diag(Sigma.est))))
lam_min        <- 0.035 * lam_max
lambdas.glasso <- exp(seq(log(lam_min), log(lam_max), length.out = n.lambda))

glasso.res   <- estimator_glasso(Sigma.est, n, lambdas.glasso, gamma = gamma)
Optim.glasso <- get_optimal_matrix(glasso.res$path, glasso.res$path.loss)

p_glasso <- make_plot_matrix(Optim.glasso$Theta_opt, "GLASSO",
                             base_size = 6, title_size = 8,
                             axis_title_size = 6, axis_text_size = 5,
                             tick_length_pt = 1, y_lab = "", x_lab = "")
p_alpha_glasso <- make_alpha_plot(Optim.glasso$alpha, "GLASSO")


######################
# Cor-GLASSO
######################
lam_max         <- 0.3  * max(abs(Corr.est - diag(diag(Corr.est))))
lam_min         <- 0.05 * lam_max
lambdas.cglasso <- exp(seq(log(lam_min), log(lam_max), length.out = n.lambda))

cglasso.res   <- estimator_corglasso(Corr.est, n, lambdas.cglasso, gamma = 0.5)
Optim.cglasso <- get_optimal_matrix(cglasso.res$path, cglasso.res$path.loss)

p_corglasso <- make_plot_matrix(Optim.cglasso$Theta_opt, "corr-GLASSO",
                                base_size = 6, title_size = 8,
                                axis_title_size = 6, axis_text_size = 5,
                                tick_length_pt = 1, y_lab = "", x_lab = "")
p_alpha_corglasso <- make_alpha_plot(Optim.cglasso$alpha, "Corr-GLASSO")


######################
# PCGLASSO
######################
lambdas.pcglasso <- seq(1, 0.05, length.out = n.lambda)
R0 <- diag(nrow = nrow(Sigma.est))

# ~6 minutes on Apple M2
res.pgclasso <- estimator_pcglasso(Sigma.est, n, R_start = R0,
                                   lambdas = lambdas.pcglasso,
                                   verbose = 1, max_iter = 1000,
                                   solver_R = "primal")

Optim.pcglasso <- get_optimal_matrix(res.pgclasso$path, res.pgclasso$path.loss)

p_pcglasso <- make_plot_matrix(Optim.pcglasso$Theta_opt, "PCGLASSO",
                               base_size = 6, title_size = 8,
                               axis_title_size = 6, axis_text_size = 5,
                               tick_length_pt = 1, y_lab = "", x_lab = "")
p_alpha_pcglasso <- make_alpha_plot(Optim.pcglasso$alpha, "PCGLASSO")

p_pcglasso2 <- make_plot_matrix_binary_highlight_rc(
  cov2cor(Optim.pcglasso$Theta_opt), "PCGLASSO",
  highlight_index = 53, highlight_label = "EEF2"
)
ggsave("./outputs/Section_4_2_adj_matrix_pcglasso.png",
       plot = p_pcglasso2, width = 3, height = 3)


######################
# Thresholded empirical partial correlations (reference)
######################
q90           <- quantile(abs(PC.est - diag(diag(PC.est))), 0.90)
PC.est.thres  <- PC.est
PC.est.thres[abs(PC.est) < q90] <- 0

alpha             <- apply(PC.est.thres, 2, function(x) sum(abs(x)) - 1)
p_alpha_emperical <- make_alpha_plot(alpha, "Thresholded")
p_emp <- make_plot_matrix(PC.est.thres, "Thresholded",
                          base_size = 6, title_size = 8,
                          axis_title_size = 6, axis_text_size = 5,
                          tick_length_pt = 1, y_lab = "", x_lab = "")

fig <- (p_glasso | p_corglasso) / (p_pcglasso | p_emp)
print(fig)
ggsave("./outputs/Appendix_B_matrices_optimal.png",
       plot = fig, width = 8, height = 8)


######################
# EBIC path comparison
######################
df_all <- rbind(
  data.frame(Edges  = glasso.res$path.loss$nEdges,
             BIC    = glasso.res$path.loss$BIC_gamma,
             Method = "GLASSO"),
  data.frame(Edges  = cglasso.res$path.loss$nEdges,
             BIC    = cglasso.res$path.loss$BIC_gamma,
             Method = "Cor-GLASSO"),
  data.frame(Edges  = res.pgclasso$path.loss$nEdges,
             BIC    = res.pgclasso$path.loss$BIC_gamma,
             Method = "PCGLASSO")
)

colors_named <- c("GLASSO" = "#1b9e77", "Cor-GLASSO" = "#7570b3", "PCGLASSO" = "#d95f02")

fig <- ggplot(df_all, aes(x = Edges, y = BIC, colour = Method)) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = colors_named) +
  scale_x_continuous(
    name   = "Number of edges",
    breaks = seq(1000, max(df_all$Edges), by = 500),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_y_continuous(
    name   = expression(EBIC(gamma == 0.5)),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  coord_cartesian(
    xlim = c(1000, max(df_all$Edges)),
    ylim = c(min(df_all$BIC), max(df_all[df_all$Edges >= 1000, "BIC"]))
  ) +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE)) +
  theme_bw(base_size = 11) +
  theme(
    legend.position      = "bottom",
    legend.title         = element_blank(),
    legend.text          = element_text(size = 8),
    legend.key.width     = grid::unit(1.4, "lines"),
    legend.box           = "horizontal",
    panel.grid.minor     = element_blank(),
    panel.grid.major.x   = element_blank(),
    panel.border         = element_rect(linewidth = 0.4),
    plot.background      = element_rect(fill = "white", colour = NA),
    panel.background     = element_rect(fill = "white", colour = NA)
  )
print(fig)
ggsave("./outputs/Appendix_B_BIC_fig.png",
       plot = fig, width = 5.2, height = 3, dpi = 300, bg = "white")


######################
# Hub diagnostics
######################
cat("HUBS identified by Hub method: SCARNA7, MIR3609, SEMG1, SEMG2, RN7SK\n")
cat("Largest diagonal =", gene.names[order(diag(Optim.pcglasso$Theta_opt), decreasing = TRUE)[1:5]], "\n")

# alpha^2(K-hat): squared L2 row norm of precision matrix
alpha2_precision <- apply(Optim.pcglasso$Theta_opt, 2, function(x) sum(x^2))
cat("alpha^2 top genes =", gene.names[order(alpha2_precision, decreasing = TRUE)[1:5]], "\n")

fig <- make_alpha_plot(alpha2_precision,
                       "$\\textbf{Sorted}\\;\\alpha^2(\\hat{K})$",
                       latex = TRUE, names = gene.names, print.number.names = 5)
ggsave("./outputs/Section_4_2_alpha2_pcglasso.png", fig, width = 4, height = 3)

# alpha^1(R-hat): absolute row sum of partial-correlation matrix
alpha1_precision <- apply(cov2cor(Optim.pcglasso$Theta_opt), 2, function(x) sum(abs(x)))
fig <- make_alpha_plot(alpha1_precision,
                       "$\\textbf{Sorted}\\;\\alpha^1(\\hat{R})$",
                       latex = TRUE, names = gene.names, print.number.names = 1)
ggsave("./outputs/Section_4_2_alpha1_pcglasso.png", fig, width = 4, height = 3)

# D(K-hat): square root of diagonal entries
p <- nrow(Optim.pcglasso$Theta_opt)
d <- sqrt(diag(Optim.pcglasso$Theta_opt))
fig <- make_alpha_plot(d, "$\\textbf{Sorted}\\;D(\\hat{K})$",
                       latex = TRUE, names = gene.names, print.number.names = 5,
                       ylab = expression(D))
ggsave("./outputs/Section_4_2_d_pcglasso.png", fig, width = 4, height = 3)

# Collinearity check: removing SEMG2 (index 12) and re-ranking alpha^2
Sigma.est_1   <- solve(Optim.pcglasso$Theta_opt)
Sigma.est_1   <- Sigma.est_1[-12, -12]
Theta_opt_1   <- solve(Sigma.est_1)
alpha2_no_semg2 <- apply(Theta_opt_1, 2, function(x) sum(x^2))
cat("alpha^2 top genes (SEMG2 removed) =",
    gene.names[-12][order(alpha2_no_semg2, decreasing = TRUE)[1:5]], "\n")


######################
# Normalized precision rows for top-diagonal genes
######################
K            <- 5
top_diag_idx <- order(diag(Optim.pcglasso$Theta_opt), decreasing = TRUE)[1:K]
cols         <- seq_len(ncol(Prec.est))

df <- do.call(rbind, lapply(top_diag_idx, function(i) {
  data.frame(
    row_idx  = i,
    row_gene = gene.names[i],
    col_idx  = cols,
    col_gene = gene.names[cols],
    value    = as.numeric(Prec.est[i, ])
  )
}))
df        <- subset(df, row_idx != col_idx)
df$max_abs <- ave(abs(df$value), df$row_idx, FUN = max, na.rm = TRUE)
df$norm   <- ifelse(df$max_abs > 0, df$value / df$max_abs, 0)
df$label  <- ifelse(df$col_idx %in% top_diag_idx & abs(df$norm) >= 0.2, df$col_gene, "")
df$panel  <- paste0("Row: ", df$row_gene)
df$pos    <- ave(df$col_idx, df$panel, FUN = seq_along)

g <- ggplot(df, aes(x = pos, y = norm)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_point(color = "black", size = 1.7) +
  geom_text_repel(
    data = subset(df, label != ""),
    aes(label = label),
    min.segment.length = 0, box.padding = 0.2, max.overlaps = Inf, size = 3
  ) +
  facet_wrap(~ panel, ncol = 3) +
  coord_cartesian(ylim = c(-1, 0.5)) +
  labs(x = "Columns", y = "Normalized values") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text   = element_text(face = "bold")
  )

ggsave("./outputs/Appendix_B_D_rows_cancer.png", g, width = 4, height = 3)

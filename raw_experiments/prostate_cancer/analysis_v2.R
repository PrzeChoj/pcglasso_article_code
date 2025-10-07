#source("raw_experiments/estimation_function.R")
library(PCGLASSO)
library(ggplot2)
source("../hub_methods/013_Method_IPCHD.R")
source("../hub_methods/012_Method_VariableScreening.R")
source("../hub_methods/011_Method_MatrixThresholding.R")
source("../estimation_function.R")

gamma <- 0.5
#df_reduced <- read.csv("./raw_experiments/prostate_cancer/312_df_reduced.csv")
df_reduced <- read.csv("312_df_reduced.csv")
n.lambda <- 20
df_reduced <- df_reduced[, -1]
df_genes <- df_reduced[, 7 + 1:200]
gene.names <- colnames(df_genes)
Sigma.est <- cov(df_genes)
Corr.est <- cov2cor(Sigma.est)
n <- dim(df_genes)[1]
Prec.est <- solve(Sigma.est)
PC.est <- cov2cor(Prec.est)
##
# begin by looking at the raw data.
##
alpha <- apply(PC.est,2, function(x){sum(abs(x))-1})
p_alpha_emperical <- make_alpha_plot(alpha,"Emperical estimated alpha")




######################
#glasso method
#######################
lam_max <- 0.2*max(abs(Sigma.est - diag(diag(Sigma.est))))
lam_min <- 0.067 * lam_max
lambdas.glasso <- exp(seq(log(lam_min),log(lam_max), length.out = n.lambda))
glasso.res   <- estimator_glasso(Sigma.est,n, lambdas.glasso,gamma = gamma)
Optim.glasso <- get_optimal_matrix(glasso.res$path, glasso.res$path.loss)
p_glasso    <- make_plot_matrix(Optim.glasso$Theta_opt, "GLasso",
                                base_size = 6, title_size = 8,
                                axis_title_size = 6, axis_text_size = 5,
                                tick_length_pt = 1)
p_alpha_glasso <- make_alpha_plot(Optim.glasso$alpha,"GLasso")



######################
#corrglasso method
#######################
lam_max <- 0.3*max(abs(Corr.est - diag(diag(Corr.est))))
lam_min <- 0.09 * lam_max
lambdas.cglasso <- exp(seq( log(lam_min),log(lam_max), length.out = n.lambda))
cglasso.res   <- estimator_corglasso(Corr.est,n, lambdas.cglasso,gamma = 0.5)
Optim.cglasso <- get_optimal_matrix(cglasso.res$path, cglasso.res$path.loss)
p_corglasso    <- make_plot_matrix(Optim.cglasso$Theta_opt, "corrGLasso",
                                   base_size = 6, title_size = 8,
                                   axis_title_size = 6, axis_text_size = 5,
                                   tick_length_pt = 1)
p_alpha_corglasso <- make_alpha_plot(Optim.cglasso$alpha,"corrGLasso")


######################
#pcglasso method
#######################
lambdas.pcglasso <- seq(1,0.1, length.out = n.lambda)
pcglasso.est <- lambda_grid(Sigma.est,
                            0,
                            lambdas= lambdas.pcglasso,
                            max.iter=500,
                            Q_init=PC.est,
                            verbose = T)
precision_array <- array(NA, dim=c(dim(Sigma.est)[1], dim(Sigma.est)[2], length(lambdas.pcglasso)))
for(i in 1:length(lambdas.pcglasso)){
  precision_array[,,i] <- pcglasso.est$solutions[[i]]$Sinv

}

loss_path <- evaluate_objective_path(precision_array, Sigma.est, n, gamma = gamma)

Optim.pcglasso <- get_optimal_matrix(precision_array, loss_path)
p_pcglasso    <- make_plot_matrix(Optim.pcglasso$Theta_opt, "PCGLasso",
                                  base_size = 6, title_size = 8,
                                   axis_title_size = 6, axis_text_size = 5,
                                   tick_length_pt = 1)
p_alpha_pcglasso <- make_alpha_plot(Optim.pcglasso$alpha,"PCGLasso")




#screening method threshold method
# does not screen since n suff big compared to p
Res <- screening_vars_SMZL2024(X=df_genes,
                               mat_type="cov",
                               mat= Sigma.est,
                               var_inds=1:dim(df_genes)[2],
                               method="max")
Res <- sta_thresholding_perc(X = df_genes,
                      mat_type = "cov",
                      mat = Sigma.est,
                      var_inds = 1:dim(df_genes)[2],
                      true_mat=NULL, perc = 0.05)
Res.ipchd <- sta_ipchd(X = df_genes,
          mat_type = "cov",
          mat = Res$mat,
          var_inds = 1:dim(df_genes)[2],
          overest_type = "frac")
# does not give reasonable results
#space.res <- estimator_space(Sigma.est,n,lambdas.pcglasso2, data=df_genes, gamma= gamma, min_scale=log(0.2),max_scale = log(1))
q90 <- quantile(abs(PC.est-diag(diag(PC.est))),.90)
ind <- abs(PC.est) < q90
PC.est.thres <- PC.est
PC.est.thres[ind] <- 0

alpha <- apply(PC.est.thres,2, function(x){sum(abs(x))-1})
p_alpha_emperical <- make_alpha_plot(alpha,"Thresholded")
p_emp <- make_plot_matrix(PC.est.thres, "Thresholded",
                          base_size = 6, title_size = 8,
                          axis_title_size = 6, axis_text_size = 5,
                          tick_length_pt = 1)

fig <- ((p_glasso | p_corglasso) /
          (p_pcglasso |p_emp ))
print(fig)
ggsave(
  "matrices_optimal_cancer.png",
  plot = fig, width = 7, height = 4
)

alphas <- list(
  "GLasso"     = Optim.glasso$alpha,
  "Cor-GLasso" = Optim.cglasso$alpha,
  "PC-GLasso"   = Optim.pcglasso$alpha,
  "Threshold"      = alpha
)
fig_alpha <- make_alpha_grid(alphas, ncol = 2, common_y = FALSE)
print(fig_alpha)
ggsave("alphas_grid_caner.png", plot = fig_alpha, width = 7, height = 4, dpi = 300)




df_glasso <- data.frame(
  Edges  = glasso.res$path.loss$nEdges,
  BIC    = glasso.res$path.loss$BIC_gamma,
  Method = "GLasso"
)

df_corglasso <- data.frame(
  Edges  = cglasso.res$path.loss$nEdges,
  BIC    = cglasso.res$path.loss$BIC_gamma,
  Method = "Cor-GLasso"
)


df_pcglasso <- data.frame(
  Edges  = loss_path$nEdges,
  BIC    = loss_path$BIC_gamma,
  Method = "PC-GLasso"
)
# Combine all
df_all <- rbind(df_glasso, df_corglasso, df_pcglasso)#,df_pcglasso_alpha)

# Plot
# Dynamic method label
#label_pcglasso_opt <- paste("PC-GLasso", " alpha=", round(optimal_alpha, 3), sep = "")


colors_named <- c(
  "GLasso"     = "#1b9e77",
  "Cor-GLasso" = "#7570b3",
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
  coord_cartesian( xlim = c(0, 4000))

print(fig)
ggsave(
  "BIC_fig_cancer.png",
  plot = fig, width = 7, height = 4
)
plot(sort(diag(Optim.pcglasso$Theta_opt), decreasing = T),ylab='Diagonal values',xlab='')
cat('HUBS identified by Hub method: SCARNA7, MIR3609, SEMG1, SEMG2, RN7SK \n')
cat('largest diagonal = ',gene.names[order(diag(Optim.pcglasso$Theta_opt),decreasing = T)[1:5]],'\n')


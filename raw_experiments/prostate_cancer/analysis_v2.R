

source("../estimation_function.R")

gamma <- 0.5
df_reduced <- read.csv("312_df_reduced.csv")
n.lambda <- 10

df_reduced <- df_reduced[, -1]
df_genes <- df_reduced[, 7 + 1:200]
Sigma.est <- cov(df_genes)
Corr.est <- cov2cor(Sigma.est)
n <- dim(df_genes)[1]
Prec.est <- solve(Sigma.est)
PC.est <- cov2cor(Prec.est)
##
# begin by looking at the raw data.
##
alpha <- apply(PC.est,2, function(x){sum((abs(x)^(3/4)))-1})
plot(sort(alpha))
cat(order(alpha,decreasing = TRUE)[1:2], "\n")


lam_max <- 0.2*max(abs(Sigma.est - diag(diag(Sigma.est))))
lam_min <- 0.067 * lam_max
lambdas.glasso <- exp(seq(log(lam_min),log(lam_max), length.out = n.lambda))
glasso.res   <- estimator_glasso(Sigma.est,n, lambdas.glasso,gamma = gamma)
plot(glasso.res$path.loss$nEdges,glasso.res$path.loss$BIC_gamma)
alpha.glasso <- apply(cov2cor(glasso.res$path[,,20]),2, function(x){sum((abs(x)^(3/4)))-1})
glasso.best <- which.min(glasso.res$path.loss$BIC_gamma)
alpha.glasso <- apply(cov2cor(glasso.res$path[,,glasso.best]),2, function(x){sum((abs(x)^(3/4)))-1})
cat(order(alpha.glasso,decreasing = TRUE)[1:2], "\n")


lam_max <- 0.3*max(abs(Corr.est - diag(diag(Corr.est))))
lam_min <- 0.09 * lam_max
lambdas.cglasso <- exp(seq( log(lam_min),log(lam_max), length.out = n.lambda))
cglasso.res   <- estimator_corglasso(Corr.est,n, lambdas.cglasso,gamma = 0.5)
plot(cglasso.res$path.loss$nEdges,cglasso.res$path.loss$BIC_gamma)
cglasso.best <- which.min(cglasso.res$path.loss$BIC_gamma)
alpha.cglasso <- apply(cov2cor(cglasso.res$path[,,cglasso.best]),2, function(x){sum((abs(x)^(3/4)))-1})
cat(order(alpha.cglasso,decreasing = TRUE)[1:2], "\n")


lam_max <- 0.2*max(abs(Sigma.est - diag(diag(Sigma.est))))
lam_min <- 0.4 * lam_max
lambdas.pcglasso <- 2*exp(seq(log(lam_max),log(lam_min), length.out = n.lambda))
pcglasso.est <- estimator_pcglasso(Sigma.est, n, lambdas.pcglasso,gamma=gamma)

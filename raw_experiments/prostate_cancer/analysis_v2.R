#source("raw_experiments/estimation_function.R")
source("../estimation_function.R")

gamma <- 0.5
#df_reduced <- read.csv("./raw_experiments/prostate_cancer/312_df_reduced.csv")
df_reduced <- read.csv("312_df_reduced.csv")
n.lambda <- 20

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
alpha <- apply(PC.est,2, function(x){sum(abs(x))-1})
plot(sort(alpha))
pc.est.max.alpha <- order(alpha,decreasing = TRUE)[1:2]
cat(colnames(PC.est)[order(alpha,decreasing = TRUE)[1:2]], "\n")


lam_max <- 0.2*max(abs(Sigma.est - diag(diag(Sigma.est))))
lam_min <- 0.067 * lam_max
lambdas.glasso <- exp(seq(log(lam_min),log(lam_max), length.out = n.lambda))
glasso.res   <- estimator_glasso(Sigma.est,n, lambdas.glasso,gamma = gamma)
plot(glasso.res$path.loss$nEdges,glasso.res$path.loss$BIC_gamma)
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


lam_max <- 0.1*max(abs(Sigma.est - diag(diag(Sigma.est))))
lam_min <- 0.3 * lam_max
lambdas.pcglasso <- 2*exp(seq(log(lam_max),log(lam_min), length.out = n.lambda))
pcglasso.est <- estimator_pcglasso(Sigma.est, n, lambdas.pcglasso,gamma=gamma, verbose = 1) # 55 minutes
plot(pcglasso.est$path.loss$nEdges,pcglasso.est$path.loss$BIC_gamma)
print(rbind(
  pcglasso.est$path.loss$nEdges,
  pcglasso.est$path.loss$loglik
))
pcglasso.est <- lambda_grid(Sigma.est,0,lambdas= lambdas.pcglasso, max.iter=500,verbose = T)
# [1,]       0.0    939.00   1001.00   1154.00   1248.00   1449.00   1611.00   1748.00   1939.00   2113.00
# [2,] -162611.9 -91257.29 -88302.49 -83933.69 -81556.08 -77028.26 -75028.52 -73155.04 -71503.01 -69978.31
precision_array <- array(NA, dim=c(dim(Sigma.est)[1], dim(Sigma.est)[2], length(lambdas.pcglasso)))
for(i in 1:length(lambdas.pcglasso)){
  precision_array[,,i] <- pcglasso.est$solutions[[i]]$Sinv
  alpha.temp <- apply(cov2cor(precision_array[,,i]),
                      2, function(x){sum(abs(x))-1})
  cat(colnames(PC.est)[order(alpha.temp,decreasing = TRUE)[1:3]], "\n")
}
loss_path <- evaluate_objective_path(precision_array, Sigma.est, n, gamma = gamma)
plot(loss_path$nEdges,loss_path$BIC_gamma)


alpha.pc <- apply(cov2cor(precision_array[,,which.min(loss_path$BIC_gamma)]),
                  2, function(x){sum(abs(x))-1})
cat(colnames(PC.est)[order(alpha.pc,decreasing = TRUE)[1:3]], "\n")
pc.lasso.est.max.alpha <- order(alpha.pc,decreasing = TRUE)[1:2]
lambdas.pcglasso2 <- seq(1,0.1, length.out = n.lambda)
pcglasso.est2 <- lambda_grid(Sigma.est,0,lambdas= lambdas.pcglasso2, max.iter=500,Q_init=PC.est,verbose = T)
precision_array2 <- array(NA, dim=c(dim(Sigma.est)[1], dim(Sigma.est)[2], length(lambdas.pcglasso2)))
for(i in 1:length(lambdas.pcglasso2)){
  precision_array2[,,i] <- pcglasso.est2$solutions[[i]]$Sinv
  alpha.temp <- apply(cov2cor(precision_array2[,,i]),
                                  2, function(x){sum(abs(x))-1})
  cat(colnames(PC.est)[order(alpha.temp,decreasing = TRUE)[1:3]], "\n")
}

loss_path2 <- evaluate_objective_path(precision_array2, Sigma.est, n, gamma = gamma)

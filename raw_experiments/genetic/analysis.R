library(glmnet)
library(PCGLASSO)
library(glasso)

n.alpha <- 1
n.lambda <- 5

source("./expermient/sangerdata/load_S_train.R")

alpha.grid <- 0 # seq(-0.1, 0.1, length.out = n.alpha)
S_max <- max(abs(S.train[row(S.train) != col(S.train)]))
lambda.max <- 4 * S_max
lambdas <- c(exp(seq(log(0.3), log(1.36), length.out = 2 * n.lambda)))
hyperparam <- c()
for (i in 1:n.alpha) {
  alphas <- rep(alpha.grid[i], n.lambda)
  hyperparam <- rbind(hyperparam, cbind(lambdas, alphas))
}


hyperparam.2 <- c()
lambdas.2 <- c(exp(seq(log(0.04), log(0.2), length.out= 2*n.lambda))) #removed 18
for (i in 1:n.alpha) {
  alphas <- rep(alpha.grid[i], n.lambda)
  hyperparam.2 <- rbind(hyperparam.2, cbind(lambdas.2, alphas))
}

# glasso
res.glasso <- PCGLASSO:::optim.wrapper(hyperparam.2, PCGLASSO:::glasso.optim.func, S.train, Precision.hat = NULL)
res.loss.dlasso <- PCGLASSO:::loss.evaluation(res.glasso, S.train, n)
res.diag.glasso <- PCGLASSO:::optim.wrapper(hyperparam.2, PCGLASSO:::glasso.plus.diag.optim.func, S.train, Precision.hat = NULL)
res.loss.diag.glasso <- PCGLASSO:::loss.evaluation(res.diag.glasso, S.train, n)

# correlation pglasso
C.train <- PCGLASSO:::unit.diag.mat(S.train)
hyperparam.2.corr <- hyperparam.2
hyperparam.2.corr[, 1] <-  c(exp(seq(log(0.2), log(0.8), length.out= 2*n.lambda))) #-18 removed
res.glasso.corr <- PCGLASSO:::optim.wrapper(hyperparam.2.corr, PCGLASSO:::glasso.optim.func, C.train, Precision.hat = NULL)
res.glasso.corr_diag <- res.glasso.corr
for (i in 1:dim(res.glasso.corr_diag)[3]) {
  Diag <- diag(res.glasso.corr_diag[, , i])
  Par <- cov2cor(res.glasso.corr_diag[, , i])
  resD <- coord.diag(
    A = Par * C.train,
    starting_point = sqrt(Diag),
    tol = 1e-6,
    alpha = 0
  )
  res.glasso.corr_diag[, , i] <- diag(resD$D) %*% Par %*% diag(resD$D)
}


for (i in 1:dim(res.glasso.corr)[3]) {
  res.glasso.corr[, , i] <- diag(sqrt(1 / diag(S.train))) %*% res.glasso.corr[, , i] %*% diag(sqrt(1 / diag(S.train)))
  res.glasso.corr_diag[, , i] <- diag(sqrt(1 / diag(S.train))) %*% res.glasso.corr_diag[, , i] %*% diag(sqrt(1 / diag(S.train)))
}
res.loss.glasso.corr <- PCGLASSO:::loss.evaluation(res.glasso.corr, S.train, n)
res.loss.glasso.corr_diag <- PCGLASSO:::loss.evaluation(res.glasso.corr_diag, S.train, n)



# pcglasso
res.pcglasso <- PCGLASSO:::optim.wrapper(hyperparam, PCGLASSO:::pcglasso.optim.func, S.train, Precision.hat = diag(1 / diag(S.train)), ordered = TRUE)
res.loss.pcglasso <- PCGLASSO:::loss.evaluation(res.pcglasso, S.train, n)

# space
library(space)
n <- dim(X)[1]
p <- dim(X)[2]
l1 <- 1/sqrt(n)*qnorm(1-1/(2*p^2))
scale_test <- seq(0.5, 2, length.out= 2*n.lambda)
res.space <- array(dim=c(p,p,length(scale_test)))
for(i in 1:length(scale_test)){
  result.space=space.joint(as.matrix(scale(X)), lam1=n.train*l1*scale_test[i], lam2=0,iter=10)
  result.space$sig.fit <- attr(scale(X),'scaled:scale')^2*result.space$sig.fit
  res.space[,,i] <-  diag(sqrt(1/result.space$sig.fit))%*%result.space$ParCor%*%diag(sqrt(1/result.space$sig.fit))
}
make_plot_matrix(my_matrix=res.space[,,7], my_title = "space")
res.loss.space <- PCGLASSO:::loss.evaluation(res.space, S.train, n)

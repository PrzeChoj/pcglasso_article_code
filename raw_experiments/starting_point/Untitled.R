library(mvtnorm)
library(pcglassoFast)
library(microbenchmark)

draw_normal <- function(n, K) {
  p <- nrow(K)
  sigma <- solve(K)
  rmvnorm(n, mean = rep(0, p), sigma = sigma)
}

is_good_pattern <- function(S_est, pattern, tol = 0.0001) {
  p <- nrow(S_est)
  for (i in 1:p) {
    for (j in 1:p) {
      if (pattern[i, j] == 0 && abs(S_est[i, j]) > tol){
        return(FALSE)
      }
      if (pattern[i, j] == 1 && S_est[i, j] < -tol){
        return(FALSE)
      }
      if (pattern[i, j] == -1 && S_est[i, j] > tol){
        return(FALSE)
      }
    }
  }

  return(TRUE)
}

mean_is_good_pattern <- function(S_est_list, pattern, tol = 0.0001) {
  mean(sapply(S_est_list, function(S) {is_good_pattern(S, pattern, tol)}))
}

KAR <- function(p, a, b, c) {
  K <- matrix(0, nrow = p, ncol = p)
  diag(K) <- b
  K[1, 1] <- a
  K[p, p] <- a

  for (i in 2:p) {
    K[i-1, i] <- K[i, i-1] <- c
  }

  K
}
Khub <- function(p, a, b, c) {
  K <- matrix(0, nrow = p, ncol = p)
  K[, 1] <- K[1, ] <- c
  diag(K) <- b
  K[1, 1] <- a
  K
}

is_K_positive_definite <- function(K) {
  min(eigen(K, TRUE, TRUE)$values) > 0.00001
}

which_K <- 1 # 1 -> Hub, 2 -> AR
n <- 400
p <- 150
a <- c(250, 3)[which_K]
b <- 1
c_val <- c(1.2, 0.46)[which_K]

K <- list(Khub, KAR)[[which_K]](p, a, b, c_val)
stopifnot(is_K_positive_definite(K))

K_pattern <- sign(K)

X <- draw_normal(n, K)
S <- cov(X)


alpha <- 4/20
nlambda <- 100
lam_max <- max(abs(S - diag(diag(S))))
lam_min <- 0.0001 * lam_max
lambdas <- exp(seq(log(lam_max), log(lam_min), length.out = nlambda))

# paths:
ud <- function(){
  R0_big_lambda <- diag(nrow(S))
  R0_inv_big_lambda <- solve(R0_big_lambda)

  pcglassoPath(S, alpha, lambdas = lambdas, R0 = R0_big_lambda, R0_inv = R0_inv_big_lambda)
}
du <- function(){
  R0_inv_small_lambda <- cov2cor(S)
  R0_small_lambda <- solve(R0_inv_small_lambda)

  pcglassoPath(S, alpha, lambdas = rev(lambdas), R0 = R0_small_lambda, R0_inv = R0_inv_small_lambda)
}
du2 <- function(){
  R0 <- cov2cor(solve(S))
  R0_inv <- solve(R0)

  pcglassoPath(S, alpha, lambdas = rev(lambdas), R0 = R0, R0_inv = R0_inv)
}
du3 <- function(){
  R0 <- cov2cor(S)
  R0_inv <- solve(R0)

  pcglassoPath(S, alpha, lambdas = rev(lambdas), R0 = R0, R0_inv = R0_inv)
}
du4 <- function(){
  R0 <- diag(nrow(S)) # identity

  pcglassoPath(S, alpha, lambdas = rev(lambdas), R0 = R0, R0_inv = R0)
}

# 6 minutes
microbenchmark(
  du(), du2(), du3(), du4(),
  times = 1
)
# Unit: seconds
# expr       min       lq     mean   median       uq      max neval
# du()  1.602924 1.712574 1.739019 1.721601 1.823328 1.831777    10
# du2() 2.319419 2.330720 2.491350 2.536469 2.577352 2.650607    10
# du3() 2.250306 2.353569 2.511989 2.469458 2.686703 2.900138    10
# du4() 2.345895 2.449619 2.540990 2.515586 2.664263 2.784367    10

# 2 minutes
microbenchmark(
  du(), ud(),
  times = 10
)
# Unit: seconds
# expr      min       lq     mean   median       uq      max neval
# du() 1.598145 1.602527 1.791756 1.746490 1.872572 2.270280    10
# ud() 2.065285 2.075605 2.220882 2.178658 2.383253 2.503978    10

udu <- function(){
  R0_big_lambda <- diag(nrow(S))
  R0_inv_big_lambda <- solve(R0_big_lambda)

  ud_path <- pcglassoPath(S, alpha, lambdas = lambdas, R0 = R0_big_lambda, R0_inv = R0_inv_big_lambda)
  pcglassoPath(S, alpha, lambdas = lambdas, R0 = ud_path$R_path[[nlambda]], R0_inv = ud_path$Ri_path[[nlambda]])
}
dud <- function(){
  R0_inv_small_lambda <- cov2cor(S)
  R0_small_lambda <- solve(R0_inv_small_lambda)

  du_path <- pcglassoPath(S, alpha, lambdas = rev(lambdas), R0 = R0_small_lambda, R0_inv = R0_inv_small_lambda)
  pcglassoPath(S, alpha, lambdas = rev(lambdas), R0 = du_path$R_path[[nlambda]], R0_inv = du_path$Ri_path[[nlambda]])
}

# 6 minutes
microbenchmark(
  udu(), dud(), {ud(); du()},
  times = 10
)
# Unit: seconds
# expr                min       lq     mean   median       uq      max neval
# udu()          4.206750 4.281121 4.503893 4.475127 4.739928 4.861764    10
# dud()          4.005379 4.182431 4.376780 4.413087 4.553913 4.832931    10
# { ud(); du() } 3.661175 3.896683 4.070793 4.048479 4.207639 4.409770    10




#####
#
# mean_is_good_pattern(ans_ud$W_path, K_pattern)
# mean_is_good_pattern(ans_udu$W_path, K_pattern)
# mean_is_good_pattern(ans_du$W_path, K_pattern)
#
#
# for (i in 1:length(ans_du$lambdas)) {
#   lambda <- ans_du$lambdas[i]
#   S_hat_ud <- ans_ud$W_path[[length(ans_du$lambdas)-i+1]]
#   f_val_ud <- pcglassoFast:::function_to_optimize(cov2cor(S_hat_ud), sqrt(diag(S_hat_ud)), S, lambda, alpha)
#   S_hat_udu <- ans_udu$W_path[[i]]
#   f_val_udu <- pcglassoFast:::function_to_optimize(cov2cor(S_hat_udu), sqrt(diag(S_hat_udu)), S, lambda, alpha)
#   S_hat_du <- ans_du$W_path[[i]]
#   f_val_du <- pcglassoFast:::function_to_optimize(cov2cor(S_hat_du), sqrt(diag(S_hat_du)), S, lambda, alpha)
#
#   if (abs(f_val_du - f_val_udu) > 0.0002) {
#     ic(i)
#     #ic(sum((sign(S_hat_du) - sign(S_hat_ud))^2))
#     ic(sum((sign(S_hat_du) - sign(S_hat_udu))^2))
#     #ic(f_val_du - f_val_ud)
#     ic(f_val_du - f_val_udu)
#   }
# }
#
# i <- 94
# is_optimum(S, lambda, alpha, sqrt(diag(S_hat_ud)), cov2cor(S_hat_ud))
# sign(cov2cor(S_hat_ud))
# is_optimum(S, lambda, alpha, sqrt(diag(S_hat_du)), cov2cor(S_hat_du))
# is_optimum(S, lambda, alpha, sqrt(diag(S_hat_udu)), cov2cor(S_hat_udu))
#
# my_res <- pcglassoFast(S, lambda, alpha, tolerance = 1e-16, max_iter = 1000)
# round(is_optimum_BIG_PI(S, lambda, alpha, my_res$D, my_res$R), 2)

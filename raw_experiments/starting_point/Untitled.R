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

Khub <- function(p, a, b, c) {
  K <- matrix(0, nrow = p, ncol = p)
  K[, 1] <- K[1, ] <- c
  diag(K) <- b
  K[1, 1] <- a
  K
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

is_K_positive_definite <- function(K) {
  min(eigen(K, TRUE, TRUE)$values) > 0.00001
}

which_experiment <- 3 # 1 -> Sanger, 2 -> Hub, 3 -> AR,
S <- if(which_experiment == 1) {
  data(Sanger)
  cov(Sanger)
} else {
  n <- 400
  p <- 150
  a <- c(NA, 250, 3)[which_experiment]
  b <- 1
  c_val <- c(NA, 1.2, 0.46)[which_experiment]

  K <- list(NA, Khub, KAR)[[which_experiment]](p, a, b, c_val)
  stopifnot(is_K_positive_definite(K))

  K_pattern <- sign(K)

  X <- draw_normal(n, K)
  cov(X)
}
times <- c(1, 10, 10)[which_experiment]


alpha <- 0.2
nlambda <- 100
lam_max <- max(abs(S - diag(diag(S))))
lam_min <- 0.0001 * lam_max
lambdas <- exp(seq(log(lam_max), log(lam_min), length.out = nlambda))



source("./raw_experiments/starting_point/paths_functions.R")

# some minutes
microbenchmark(
  du(S, alpha, lambdas),
  du2(S, alpha, lambdas),
  du3(S, alpha, lambdas),
  du4(S, alpha, lambdas),
  times = times
)
# which_experiment == 1, Sanger data:
# Unit: seconds
# expr       min       lq     mean   median       uq      max neval
# du()  24.14990 24.14990 24.14990 24.14990 24.14990 24.14990     1
# du2() 65.84658 65.84658 65.84658 65.84658 65.84658 65.84658     1
# du3() 66.34271 66.34271 66.34271 66.34271 66.34271 66.34271     1
# du4() 69.49258 69.49258 69.49258 69.49258 69.49258 69.49258     1

# which_experiment == 2, Hub:
# Unit: seconds
# expr       min       lq     mean   median       uq      max neval
# du()  5.021417 5.134240 5.227172 5.235393 5.288931 5.378933    10
# du2() 6.823953 6.846655 6.952709 6.910798 7.016368 7.184668    10
# du3() 6.840238 6.996360 7.109834 7.177190 7.185860 7.362256    10
# du4() 6.963243 7.066019 7.150748 7.178869 7.235403 7.347942    10

# which_experiment == 3, AR:
# Unit: seconds
# expr       min       lq     mean   median       uq      max neval
# du()  5.851331 5.918153 6.200424 6.106685 6.440266 6.757789    10
# du2() 5.523438 5.591263 5.988948 5.934354 6.304790 6.672629    10
# du3() 5.925678 6.056490 6.278462 6.161392 6.549268 6.820546    10
# du4() 5.761309 5.960301 6.075738 6.100656 6.191135 6.321865    10


#####
# some minutes
microbenchmark(
  du(S, alpha, lambdas),
  ud(S, alpha, lambdas),
  times = times
)
# which_experiment == 1, Sanger data:
# Unit: seconds
# expr      min       lq     mean   median       uq      max neval
# du() 23.92645 23.92645 23.92645 23.92645 23.92645 23.92645     1
# ud() 54.82631 54.82631 54.82631 54.82631 54.82631 54.82631     1

# which_experiment == 2, Hub:
# Unit: seconds
# expr      min       lq     mean   median       uq      max neval
# du() 5.017217 5.117969 5.236935 5.290921 5.325416 5.382931    10
# ud() 6.218553 6.372678 6.551896 6.513973 6.708809 6.919347    10

# which_experiment == 3, AR:
# Unit: seconds
# expr      min       lq     mean   median       uq      max neval
# du() 5.551907 5.678631 5.895191 5.742018 6.198952 6.322515    10
# ud() 5.277002 5.289476 5.462313 5.418437 5.571002 5.903453    10

#####


# some minutes
microbenchmark(
  { ud(S, alpha, lambdas); du(S, alpha, lambdas) },
  dud(S, alpha, lambdas),
  udu(S, alpha, lambdas),
  times = times
)
# which_experiment == 1, Sanger data:
# Unit: seconds
# expr                min        lq      mean    median        uq       max neval
# { ud(); du() } 78.16442  78.16442  78.16442  78.16442  78.16442  78.16442     1
# dud()          94.20995  94.20995  94.20995  94.20995  94.20995  94.20995     1
# udu()         107.68915 107.68915 107.68915 107.68915 107.68915 107.68915     1

# which_experiment == 2, Hub:
# Unit: seconds
# expr                min       lq     mean   median       uq      max neval
# { ud(); du() } 11.10811 11.28072 13.84569 11.38007 11.56512 27.68010    10
# dud()          11.78827 12.15967 13.90151 12.33848 12.78654 26.01385    10
# udu()          12.43207 12.72229 13.98820 12.89214 12.99030 24.43083    10

# which_experiment == 3, AR:
# Unit: seconds
# expr                min       lq     mean   median       uq      max neval
# { ud(); du() } 10.76450 10.91750 11.13102 11.12650 11.27459 11.51221    10
# dud()          11.47326 11.53124 11.83025 11.82397 11.97280 12.38797    10
# udu()          10.33583 10.70178 10.85227 10.77036 11.11216 11.27507    10


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

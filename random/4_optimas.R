library(pcglassoFast)
S <- structure(c(
  33.8033029696703, 34.5199249714435, 39.3327613721754,
  34.5199249714435, 36.8145049198519, 41.1258632065145,
  39.3327613721754, 41.1258632065145, 47.4505779256602
), dim = c(3L, 3L))
p <- 3
lambda <- 2.2
alpha <- 0

#####
R_rand <- matrix(rnorm(p*p), ncol = 3); R_rand <- cov2cor(R_rand %*% t(R_rand))
pcglassoFast_ans <- pcglassoFast(
  S = S, lambda = lambda, alpha = alpha,
  max_iter = 10000, tolerance = 1e-8,
  tol_R = 1e-8,
  max_iter_R_inner = 10000, max_iter_R_outer = 10000,
  tol_D = 1e-8, max_iter_D_newton = 10000, max_iter_D_ls = 1000,
  R = R_rand
)
pcglassoFast_ans$R
pcglassoFast_ans$D

pcglassoFast_ans2 <- pcglassoFast(
  S = S, lambda = lambda, alpha = alpha,
  max_iter = 10000, tolerance = 1e-18,
  tol_R = 1e-8,
  max_iter_R_inner = 10000, max_iter_R_outer = 10000,
  tol_D = 1e-8, max_iter_D_newton = 10000, max_iter_D_ls = 1000,
  R = pcglassoFast_ans$R, R_inv = pcglassoFast_ans$R_inv, D = pcglassoFast_ans$D
)
pcglassoFast_ans2$R
pcglassoFast_ans2$D

#####
R_opt_0 <- structure(c(
  1, 0, 0,
  0, 1, 0,
  0, 0, 1),
  dim = c(3L, 3L))
d_opt_0 <- c(0.171996825393186, 0.164812640650261, 0.145170792220918)

R_opt_1 <- structure(c(
  1.0000000000000000, 0.0000000000000000, -0.599002266176439,
  0.0000000000000000, 1.0000000000000000, -0.789197120377116,
  -0.599002266176439, -0.789197120377116, 1.0000000000000000),
  dim = c(3L, 3L))
d_opt_1 <- c(0.804092872623089, 0.997787837569754, 1.10088717620153)

R_opt_2 <- structure(c(
  1.0000000000000000, -0.548595049305909, -0.823295706792424,
  -0.548595049305909, 1.0000000000000000, 0.0000000000000000,
  -0.823295706792424, 0.0000000000000000, 1.0000000000000000
), dim = c(3L, 3L))
d_opt_2 <- c(1.2177391619501, 0.667124810979768, 0.855671143807889)

R_opt_3 <- structure(c(
  1.0000000000000000, -0.497965471387785, 0.0000000000000000,
  -0.497965471387785, 1.0000000000000000, -0.855983362747829,
  0.0000000000000000, -0.855983362747829, 1.0000000000000000
), dim = c(3L, 3L))
d_opt_3 <- c(0.664551051101129, 1.21928858456384, 0.927303187183783)


#####
is_optimum_BIG_PI <- function(S, lambda, alpha, D, R, R_inv = solve(R)) {
  p <- nrow(S)
  C_hat <- cov2cor(S)

  D <- sqrt(diag(S))*D

  J_prim <- matrix(rep(1, p*p), nrow = p) - diag(p)

  BIG_PI <- (R_inv - diag(D) %*% C_hat %*% diag(D) - diag(alpha, p)) / lambda + diag(diag(J_prim %*% abs(R)))

  BIG_PI
}

#####
is_optimum_BIG_PI(S, lambda, alpha, D = d_opt_0, R = R_opt_0)
is_optimum_BIG_PI(S, lambda, alpha, D = d_opt_1, R = R_opt_1)
is_optimum_BIG_PI(S, lambda, alpha, D = d_opt_2, R = R_opt_2)
is_optimum_BIG_PI(S, lambda, alpha, D = d_opt_3, R = R_opt_3)

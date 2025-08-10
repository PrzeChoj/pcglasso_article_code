path_up_down <- function(S, alpha, lambda) {
  nlambda <- 50

  lam_max <- max(abs(S - diag(diag(S))))
  lam_min <- 0.0001 * lam_max
  lambdas <- exp(seq(log(lam_max), log(lam_min), length.out = nlambda))

  lambdas <- lambdas[lambdas > lambda]

  res <- if (length(lambdas) >= 1) {
    R0 <- diag(nrow(S))
    path_sol <- pcglassoFast::pcglassoPath(S, alpha, lambdas = lambdas, R0 = R0, R0_inv = R0)

    R <- path_sol$W_path[[length(lambdas)]]
    R_inv <- path_sol$Wi_path[[length(lambdas)]]
    pcglassoFast::pcglassoFast(S, lambda, alpha, R = R, R_inv = R_inv)
  } else {
    pcglassoFast::pcglassoFast(S, lambda, alpha)
  }

  res$loss[length(res$loss)]
}

path_up_down_up <- function(S, alpha, lambda) {
  nlambda <- 50

  lam_max <- max(abs(S - diag(diag(S))))
  lam_min <- 0.0001 * lam_max
  lambdas <- exp(seq(log(lam_max), log(lam_min), length.out = nlambda))

  # path up down
  R0 <- diag(nrow(S))
  path_ud <- pcglassoFast::pcglassoPath(S, alpha, lambdas = lambdas, R0 = R0, R0_inv = R0)

  closest_lambda_index <- which.min(abs(lambdas - lambda))

  R <- path_ud$W_path[[closest_lambda_index]]
  R_inv <- path_ud$Wi_path[[closest_lambda_index]]

  sol_ud <- pcglassoFast::pcglassoFast(S, lambda, alpha, R = R, R_inv = R_inv)

  # path down up
  R <- path_ud$W_path[[length(lambdas)]]
  R_inv <- path_ud$Wi_path[[length(lambdas)]]

  lambdas <- rev(lambdas)
  lambdas <- lambdas[lambdas < lambda]

  path_du <- pcglassoFast::pcglassoPath(S, alpha, lambdas = lambdas, R0 = R, R0_inv = R_inv)

  R <- path_du$W_path[[length(lambdas)]]
  R_inv <- path_du$Wi_path[[length(lambdas)]]
  sol_du <- pcglassoFast::pcglassoFast(S, lambda, alpha, R = R, R_inv = R_inv)

  max(sol_du$loss[length(sol_du$loss)], sol_ud$loss[length(sol_ud$loss)])
}

start_I <- function(S, alpha, lambda) {
  res <- pcglassoFast::pcglassoFast(S, lambda, alpha)

  res$loss[length(res$loss)]
}

start_cor <- function(S, alpha, lambda) {
  res <- pcglassoFast::pcglassoFast(S, lambda, alpha, R = cov2cor(solve(S)))

  res$loss[length(res$loss)]
}

start_glasso <- function(S, alpha, lambda) {
  glasso_sol <- glassoFast::glassoFast(S, rho = lambda)
  d <- sqrt(diag(glasso_sol$wi))
  d_inv <- 1/d
  R <- (glasso_sol$wi) * (d_inv %o% d_inv)
  R_inv <- (glasso_sol$w) * (d %o% d)

  res <- pcglassoFast::pcglassoFast(S, lambda, alpha, R = R, R_inv = R_inv)

  res$loss[length(res$loss)]
}

start_L2 <- function(S, alpha, lambda) {
  P_est <- rags2ridges::ridgeP(S, lambda = 0.1)

  res <- pcglassoFast::pcglassoFast(S, lambda, alpha, R = cov2cor(P_est))

  res$loss[length(res$loss)]
}

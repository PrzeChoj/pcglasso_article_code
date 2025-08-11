path_up_down_part <- function(S, alpha, lambda) {
  nlambda <- 50

  lam_max <- max(abs(cov2cor(S) - diag(nrow(S)))) + 0.001
  lam_min <- 0.0001 * lam_max
  lambdas <- exp(seq(log(lam_max), log(lam_min), length.out = nlambda))

  lambdas <- lambdas[lambdas > lambda]
  if (length(lambdas) == 0){
    return(list(
      path_ud = NA,
      pcg_sol = pcglassoFast::pcglassoFast(S, lambda, alpha)
    ))
  }

  R0 <- diag(nrow(S))
  path_ud <- pcglassoFast::pcglassoPath(S, alpha, lambdas = lambdas, R0 = R0, R0_inv = R0)

  R <- path_ud$R_path[[length(lambdas)]]
  R_inv <- path_ud$Ri_path[[length(lambdas)]]
  pcg_sol <- pcglassoFast::pcglassoFast(S, lambda, alpha, R = R, R_inv = R_inv)

  list(
    path_ud = path_ud,
    pcg_sol = pcg_sol
  )
}

path_part_up_down_full_down_up <- function(S, alpha, lambda, path_ud_part) {
  nlambda <- 50
  lam_max <- max(abs(cov2cor(S) - diag(nrow(S)))) + 0.001
  lam_min <- 0.0001 * lam_max
  lambdas <- exp(seq(log(lam_max), log(lam_min), length.out = nlambda))
  lambdas <- lambdas[lambdas <= lambda]

  # I know that length(lambdas) > 1

  path_ud <- if (length(lambdas) == 50) { # equivalent, is.na(path_ud_part)
    R0 <- diag(nrow(S))
    pcglassoFast::pcglassoPath(S, alpha, lambdas = lambdas, R0 = R0, R0_inv = R0)
  } else {
    # 1 <= length(lambdas) < 50
    R0 <- path_ud_part$R_path[[length(path_ud_part$W_path)]]
    R0_inv <- path_ud_part$Ri_path[[length(path_ud_part$Wi_path)]]
    pcglassoFast::pcglassoPath(S, alpha, lambdas = lambdas, R0 = R0, R0_inv = R0_inv)
  }

  nlambda <- 50
  lam_max <- max(abs(cov2cor(S) - diag(nrow(S)))) + 0.001
  lam_min <- 0.0001 * lam_max
  lambdas <- rev(exp(seq(log(lam_max), log(lam_min), length.out = nlambda))) # down up
  lambdas <- lambdas[lambdas <= lambda]

  # now I know that 1 <= length(lambdas)
  R0 <- path_ud$R_path[[length(path_ud$R_path)]]
  R0_inv <- path_ud$Ri_path[[length(path_ud$Ri_path)]]
  path_du <- pcglassoFast::pcglassoPath(S, alpha, lambdas = lambdas, R0 = R0, R0_inv = R0_inv)

  R0 <- path_du$R_path[[length(path_du$R_path)]]
  R0_inv <- path_du$Ri_path[[length(path_du$Ri_path)]]
  pcg_sol_2 <- pcglassoFast::pcglassoFast(S, lambda, alpha, R = R0, R_inv = R0_inv)

  pcg_sol_2$loss[length(pcg_sol_2$loss)]
}

path_up_down_up <- function(S, alpha, lambda) {
  time_ud <- system.time(
    sol_path_up_down_part <- path_up_down_part(S, alpha, lambda)
  )[["elapsed"]]

  path_ud_part <- sol_path_up_down_part$path_ud
  pcg_sol_1 <- sol_path_up_down_part$pcg_sol
  loss_ud <- pcg_sol_1$loss[length(pcg_sol_1$loss)]

  if (isFALSE(is.na(path_ud_part)) && (length(path_ud_part$iters) == 50)) {
    # lambda is smaller than all lambdas in path
    return(list(
      time_ud = time_ud,
      res_ud = loss_ud,
      time_udu = time_ud,
      res_udu = loss_ud
    ))
  }

  time_udu <- system.time(
    loss_udu <- path_part_up_down_full_down_up(S, alpha, lambda, path_ud_part)
  )[["elapsed"]]

  list(
    time_ud = time_ud,
    res_ud = loss_ud,
    time_udu = time_ud + time_udu,
    res_udu = max(loss_ud, loss_udu)
  )
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
  P_est <- rags2ridges::ridgeP(S, lambda = lambda * mean(diag(S))) # lambda is penalty for cor matrix

  res <- pcglassoFast::pcglassoFast(S, lambda, alpha, R = cov2cor(P_est))

  res$loss[length(res$loss)]
}

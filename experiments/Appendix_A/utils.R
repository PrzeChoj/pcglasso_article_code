# devtools::install_github("JackStorrorCarter/PCGLASSO")
# devtools::install_github("PrzeChoj/pcglassoFast")

library(PCGLASSO)
library(pcglassoFast)
library(MASS)

pcglasso_goal_function <- function(S, lambda, alpha, Sinv) {
  delta_matrix <- cov2cor(Sinv)
  theta_diag   <- sqrt(diag(Sinv))
  p            <- nrow(delta_matrix)

  theta_mat <- diag(theta_diag)
  log_det <- as.numeric(determinant(delta_matrix, logarithm = TRUE)$modulus)
  quad_term <- sum(diag(S %*% theta_mat %*% delta_matrix %*% theta_mat))
  l1_pen    <- lambda * sum(abs(delta_matrix - diag(p)))

  -log_det - 2 * (1 - alpha) * sum(log(theta_diag)) + quad_term + l1_pen
}

build_K_star <- function(p, cor_modifier, K_structure = c("hub","line")) {
  K_structure <- match.arg(K_structure)

  K_star <- diag(1, p)

  if (K_structure == "hub") {
    K_star[1, 2:p] <- K_star[2:p, 1] <- -cor_modifier / sqrt(p)

  } else if (K_structure == "line") {
    stopifnot(cor_modifier < 1, p > 1)
    for (i in 2:p) {
      K_star[i-1, i] <- K_star[i, i-1] <- cor_modifier / (2 * cos(pi/(p+1))) # singular for cor_modifier = 1
    }
  }

  e_min <- min(eigen(K_star, symmetric = TRUE, only.values = TRUE)$values)
  stopifnot(e_min > 1e-8)

  K_star
}

S_from_K_star <- function(K_star, n) {
  p <- nrow(K_star)
  Z <- mvrnorm(n = n, mu = rep(0, p), Sigma = solve(K_star))
  S <- t(Z) %*% Z / n
  S <- cov2cor(S)

  S
}

compute_best_value <- function(
    p,
    cor_modifier,
    lambda,
    alpha,
    best_method,
    tolerance_best,
    K_structure = c("hub","line"),
    pcglasso_tolerance_modifier = 100,
    seed = 1234) {
  set.seed(seed)
  n <- 2 * p

  K_structure <- match.arg(K_structure)
  K_star <- build_K_star(p, cor_modifier, K_structure)
  S <- S_from_K_star(K_star, n)

  tol_pcglasso <- tolerance_best * pcglasso_tolerance_modifier
  c_parameter <- 1 - alpha

  Sinv <- switch(
    best_method,
    pcglasso_I = pcglasso(
      S, lambda, c_parameter,
      Theta_start = diag(p),
      threshold   = tol_pcglasso
    ),
    pcglasso_C = pcglasso(
      S, lambda, c_parameter,
      Theta_start = solve(S),
      threshold = tol_pcglasso
    ),
    pcglassoFast_I = pcglassoFast(
      S, lambda = lambda, alpha = alpha,
      R = diag(p),
      solver_R = "fortran",
      tolerance = tolerance_best
    )$Sinv,
    pcglassoFast_C = pcglassoFast(
      S, lambda = lambda, alpha = alpha,
      R = cov2cor(solve(S)),
      solver_R = "fortran",
      tolerance = tolerance_best
    )$Sinv,
    pcglasso_cpp_I = pcglassoFast(
      S, lambda = lambda, alpha = alpha,
      R = diag(p),
      solver_R = "cpp",
      tolerance = tolerance_best
    )$Sinv,
    pcglasso_cpp_C = pcglassoFast(
      S, lambda = lambda, alpha = alpha,
      R = cov2cor(solve(S)),
      solver_R = "cpp",
      tolerance = tolerance_best
    )$Sinv,
    stop("Unknown best_method: ", best_method)
  )

  pcglasso_goal_function(S, lambda, alpha, Sinv)
}

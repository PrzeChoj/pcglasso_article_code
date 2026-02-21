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
  stopifnot(p < n) # important for future analysis

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

compute_function_value <- function(
    p,
    lambda,
    alpha,
    method,
    tolerance,
    S,
    seed = 1234,
    verbose = 0) {
  set.seed(seed)
  c_parameter <- 1 - alpha

  Sinv <- switch(
    method,
    pcglasso_I = pcglasso(
      S, lambda, c_parameter,
      Theta_start = diag(p),
      threshold   = tolerance
    ),
    pcglasso_C = pcglasso(
      S, lambda, c_parameter,
      Theta_start = solve(S),
      threshold = tolerance
    ),
    pcglassoFast_I = pcglassoFast(
      S, lambda = lambda, alpha = alpha,
      R = diag(p),
      solver_R = "fortran",
      tolerance = tolerance,
      max_iter = 10000,
      verbose = verbose
    )$Sinv,
    pcglassoFast_C = pcglassoFast(
      S, lambda = lambda, alpha = alpha,
      R = cov2cor(solve(S)),
      solver_R = "fortran",
      tolerance = tolerance,
      max_iter = 10000,
      verbose = verbose
    )$Sinv,
    pcglasso_cpp_I = pcglassoFast(
      S, lambda = lambda, alpha = alpha,
      R = diag(p),
      solver_R = "cpp",
      tolerance = tolerance,
      max_iter = 10000,
      verbose = verbose
    )$Sinv,
    pcglasso_cpp_C = pcglassoFast(
      S, lambda = lambda, alpha = alpha,
      R = cov2cor(solve(S)),
      solver_R = "cpp",
      tolerance = tolerance,
      max_iter = 10000,
      verbose = verbose
    )$Sinv,
    stop("Unknown method: ", method)
  )

  pcglasso_goal_function(S, lambda, alpha, Sinv)
}

# best method table
{
  hub_methods <- c(
    "pcglasso_C",
    "pcglasso_cpp_C",
    "pcglasso_cpp_C",
    "pcglasso_cpp_C",
    "pcglasso_cpp_I",
    "pcglasso_cpp_I",
    "pcglasso_cpp_C",
    "pcglassoFast_C",
    "pcglasso_C",
    "pcglasso_cpp_I",
    "pcglassoFast_C",
    "pcglassoFast_C",
    "pcglasso_cpp_I",
    "pcglassoFast_I",
    "pcglassoFast_I",
    "pcglassoFast_C",
    "pcglasso_cpp_I",
    "pcglasso_cpp_C",
    "pcglassoFast_C",
    "pcglassoFast_I",
    "pcglasso_cpp_C",
    "pcglassoFast_I",
    "pcglassoFast_I",
    "pcglassoFast_I",
    "pcglasso_cpp_C",
    "pcglassoFast_I",
    "pcglassoFast_C",
    "pcglassoFast_C",
    "pcglasso_cpp_C",
    "pcglassoFast_I",
    "pcglassoFast_I",
    "pcglassoFast_I"
  )
  stopifnot(length(hub_methods) == 32)
  grid_hub <- expand.grid(
    p            = c(10, 50, 100, 150),
    cor_modifier = c(1.0, 0.9),
    lambda       = c(0.1, 0.2),
    alpha        = c(0.0, 0.5),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  grid_line <- expand.grid(
    p            = c(10, 50, 100, 150),
    cor_modifier = c(0.8, 0.9),
    lambda       = c(0.1, 0.2),
    alpha        = c(0.0, 0.5),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  best_method_table <- rbind(
    transform(grid_line,
              K_structure = "line",
              best_method = "pcglassoFast_C"),
    transform(grid_hub,
              K_structure = "hub",
              best_method = hub_methods)
  )
}

get_best_method <- function(p, graph_structure, lambda, alpha) {
  cor_mod <- cor_modifier_map[[graph_structure]]

  row_best <- best_method_table[
    best_method_table$p == p &
      best_method_table$cor_modifier == cor_mod &
      best_method_table$lambda == lambda &
      best_method_table$alpha == alpha &
      best_method_table$K_structure == graph_structure,
    ,
    drop = FALSE
  ]

  if (nrow(row_best) == 1L) {
    return(row_best$best_method[1])
  }

  warning(
    "best_method not found or not unique for: ",
    "p=", p,
    " cor=", cor_mod,
    " lambda=", lambda,
    " alpha=", alpha,
    " structure=", graph_structure,
    " (fallback: pcglassoFast_C)"
  )
  "pcglassoFast_C"
}

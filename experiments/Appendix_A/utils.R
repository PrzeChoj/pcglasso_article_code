# devtools::install_github("JackStorrorCarter/PCGLASSO")
# library(devtools)
# install_local("/path/to/pcglassoFast")

library(PCGLASSO)
library(pcglassoFast)
library(MASS)

pcglasso_goal_function <- function(S, lambda, alpha, Sinv) {
  delta_matrix <- cov2cor(Sinv)
  theta_diag <- sqrt(diag(Sinv))
  p <- nrow(delta_matrix)

  theta_mat <- diag(theta_diag)
  log_det <- as.numeric(determinant(delta_matrix, logarithm = TRUE)$modulus)
  quad_term <- sum(diag(S %*% theta_mat %*% delta_matrix %*% theta_mat))
  l1_pen <- lambda * sum(abs(delta_matrix - diag(p)))

  -log_det - 2 * (1 - alpha) * sum(log(theta_diag)) + quad_term + l1_pen
}

build_K_star <- function(p, K_structure = c("hub_1", "hub_09", "AR2", "random")) {
  K_structure <- match.arg(K_structure)

  K_star <- diag(1, p)
  switch(K_structure,
    hub_1 = {
      K_star[1, 2:p] <- K_star[2:p, 1] <- -1 / sqrt(p)
    },
    hub_09 = {
      K_star[1, 2:p] <- K_star[2:p, 1] <- -0.9 / sqrt(p)
    },
    AR2 = {
      stopifnot(p > 2)
      for (i in 2:p) {
        K_star[i - 1, i] <- K_star[i, i - 1] <- 1 / 2
      }
      for (i in 3:p) {
        K_star[i - 2, i] <- K_star[i, i - 2] <- 1 / 4
      }
    },
    random = {
      # Note: In Carter's description there is no this `denominator > 0` if-statement.
      # Without it, we could have a singular K_star.
      repeat {
        K_star <- diag(1, p)
        while(sum(K_star - diag(1, p) != 0) < 3 * p / 2) {
          matrix_indexes <- sample(p, 2)
          K_star[matrix_indexes[1], matrix_indexes[2]] <- runif(1, 0.4, 1) * sample(c(-1, 1), 1)
        }
        for (i in 1:p) {
          # Note: In Carter's description there is no this `denominator > 0` if-statement.
          # Without it, we could have the 0/0 situation.
          denominator <- sum(abs(K_star[-i, i])) * 1.1
          if (denominator > 0) {
            K_star[-i, i] <- K_star[-i, i] / denominator
          }
        }

        K_star <- (K_star + t(K_star)) / 2

        if (min(eigen(K_star, symmetric = TRUE, only.values = TRUE)$values) > 1e-8) {
          break
        }
      }
    },
    stop("Unknown K_structure")
  )

  e_min <- min(eigen(K_star, symmetric = TRUE, only.values = TRUE)$values)
  stopifnot(e_min > 1e-8)

  K_star
}

S_from_K_star <- function(K_star, n) {
  p <- nrow(K_star)
  stopifnot(p < n) # important for future analysis

  Z <- mvrnorm(n = n, mu = rep(0, p), Sigma = solve(K_star))
  S <- t(Z) %*% Z / n
  S <- cov2cor(S) # as in Carter's description

  S
}

value_after_optimization <- function(S, solver, start, tol, lambda, alpha) {
  p <- nrow(S)
  c_parameter <- 1 - alpha

  starting_matrix <- switch(start,
    I = diag(p),
    C = solve(S),
    stop("Unknown start: ", start)
  )

  Sinv <- switch(solver,
    pcglasso = pcglasso(
      S, lambda, c_parameter,
      Theta_start = starting_matrix,
      threshold = tol
    ),
    pcglassoFast_Dual = pcglassoFast(
      S,
      lambda = lambda, alpha = alpha,
      R0 = cov2cor(starting_matrix),
      solver_R = "Dual",
      tolerance = tol,
      max_iter = 10000
    )$Sinv,
    pcglassoFast_Primal = pcglassoFast(
      S,
      lambda = lambda, alpha = alpha,
      R0 = cov2cor(starting_matrix),
      solver_R = "Primal",
      tolerance = tol,
      max_iter = 10000
    )$Sinv,
    stop("Unknown solver: ", solver)
  )

  pcglasso_goal_function(S, lambda, alpha, Sinv)
}

# best method table
{
  methods <- c(
    "pcglassoFast_Primal_C",
    "pcglassoFast_Primal_C",
    "pcglassoFast_Dual_I",
    "pcglassoFast_Primal_I",
    "pcglassoFast_Dual_C",
    "pcglassoFast_Primal_C",
    "pcglassoFast_Dual_I",
    "pcglassoFast_Dual_C",
    "pcglassoFast_Dual_C",
    "pcglassoFast_Primal_C",
    "pcglassoFast_Dual_C",
    "pcglassoFast_Dual_C",
    "pcglassoFast_Dual_I",
    "pcglassoFast_Primal_I",
    "pcglassoFast_Dual_C",
    "pcglassoFast_Dual_C",
    "pcglassoFast_Dual_C",
    "pcglassoFast_Primal_C",
    "pcglassoFast_Dual_I",
    "pcglassoFast_Dual_C",
    "pcglassoFast_Dual_I",
    "pcglassoFast_Dual_C",
    "pcglassoFast_Dual_I",
    "pcglassoFast_Dual_C",
    "pcglassoFast_Dual_C",
    "pcglassoFast_Dual_I",
    "pcglassoFast_Dual_I",
    "pcglassoFast_Dual_I",
    "pcglassoFast_Dual_I",
    "pcglassoFast_Dual_I",
    "pcglassoFast_Dual_I",
    "pcglassoFast_Dual_C",
    rep("pcglassoFast_Dual_I", 16), # AR2
    rep("pcglassoFast_Dual_I", 16) # random
  )
  stopifnot(length(methods) == 64)
  grid <- expand.grid(
    alpha = c(0.5, 0.0),
    lambda = c(0.1, 0.2),
    p = c(50, 100, 150, 200),
    K_structure = c("hub_1", "hub_09", "AR2", "random"),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  best_method_table <- transform(grid, best_method = methods)
}

get_best_method <- function(p, graph_structure, lambda, alpha, cor_mod = NULL) {
  row_best <- best_method_table[
    best_method_table$p == p &
      best_method_table$lambda == lambda &
      best_method_table$alpha == alpha &
      best_method_table$K_structure == graph_structure, ,
    drop = FALSE
  ]

  if (nrow(row_best) == 1L) {
    return(row_best$best_method[1])
  }

  warning("best_method not found or not unique")
  "pcglassoFast_Dual_C"
}

get_best_value <- function(S, p, graph_structure, lambda, alpha) {
  best_method <- get_best_method(p, graph_structure, lambda, alpha)
  value_after_optimization(
    S,
    substr(best_method, 1, nchar(best_method) - 2), substr(best_method, nchar(best_method), nchar(best_method)),
    1e-13, lambda, alpha
  )
}

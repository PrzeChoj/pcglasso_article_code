get_baseline <- function(p, structure) {
  row <- baseline_best_value[
    baseline_best_value$p == p & baseline_best_value$structure == structure,
    ,
    drop = FALSE
  ]
  if (nrow(row) != 1L) stop("Baseline not found/unique")
  row$best_value[1]
}

value_for_tolerance <- function(S, solver, start, tol) {
  p <- nrow(S)
  starting_matrix <- switch(
    start,
    I = diag(p),
    C = solve(S),
    stop("Unknown start: ", start)
  )

  Sinv <- switch(
    solver,
    pcglasso = pcglasso(
      S, lambda, c_parameter,
      Theta_start = starting_matrix,
      threshold   = tol
    ),
    pcglasso_fortran = pcglassoFast(
      S, lambda = lambda, alpha = alpha,
      R = cov2cor(starting_matrix),
      solver_R = "fortran",
      tolerance = tol,
      max_iter = 10000
    )$Sinv,
    pcglasso_cpp = pcglassoFast(
      S, lambda = lambda, alpha = alpha,
      R = cov2cor(starting_matrix),
      solver_R = "cpp",
      tolerance = tol,
      max_iter = 10000
    )$Sinv,
    stop("Unknown solver: ", solver)
  )

  pcglasso_goal_function(S, lambda, alpha, Sinv)
}

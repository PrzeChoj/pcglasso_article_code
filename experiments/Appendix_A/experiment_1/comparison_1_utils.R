source("./experiments/Appendix_A/utils.R")

simulate_pcglasso <- function(
    M, p, cor_modifier, lambda, alpha, tolerance_list,
    K_structure = c("hub", "line"),
    pcglasso_tolerance_modifier = 100, seed = 1234, show_progress_bar = FALSE) {
  if (!is.null(seed)) set.seed(seed)

  n       <- 2 * p
  rho     <- lambda
  c_parameter <- 1 - alpha
  n_tol   <- length(tolerance_list)

  # get data
  K_structure <- match.arg(K_structure)
  K_star <- build_K_star(p, cor_modifier, K_structure)
  S <- S_from_K_star(K_star, n)

  run_with_obj <- function(fun) {
    t0 <- proc.time()[["elapsed"]]
    Sinv <- fun()
    t1 <- proc.time()[["elapsed"]]

    list(
      time = unname(t1 - t0),
      value = pcglasso_goal_function(S, lambda, alpha, Sinv)
    )
  }

  time_pcglasso_C_mat        <- matrix(NA, nrow = M, ncol = n_tol)
  res_val_pcglasso_C_mat     <- matrix(NA, nrow = M, ncol = n_tol)
  time_pcglasso_I_mat        <- matrix(NA, nrow = M, ncol = n_tol)
  res_val_pcglasso_I_mat     <- matrix(NA, nrow = M, ncol = n_tol)
  time_pcglassoFast_I_mat    <- matrix(NA, nrow = M, ncol = n_tol)
  res_val_pcglassoFast_I_mat <- matrix(NA, nrow = M, ncol = n_tol)
  time_pcglassoFast_C_mat    <- matrix(NA, nrow = M, ncol = n_tol)
  res_val_pcglassoFast_C_mat <- matrix(NA, nrow = M, ncol = n_tol)
  time_pcglasso_cpp_I_mat    <- matrix(NA, nrow = M, ncol = n_tol)
  res_val_pcglasso_cpp_I_mat <- matrix(NA, nrow = M, ncol = n_tol)
  time_pcglasso_cpp_C_mat    <- matrix(NA, nrow = M, ncol = n_tol)
  res_val_pcglasso_cpp_C_mat <- matrix(NA, nrow = M, ncol = n_tol)

  if (show_progress_bar) {
    counter <- 0
    pb <- txtProgressBar(min = 0, max = n_tol * M, style = 3)
  }

  for (m in seq_len(M)) {
    for (i in seq_along(tolerance_list)) {
      tolerance    <- tolerance_list[i]
      tol_pcglasso <- tolerance * pcglasso_tolerance_modifier

      if (show_progress_bar) {
        counter <- counter + 1
        setTxtProgressBar(pb, counter)
      }

      # pcglasso_C
      res <- run_with_obj(function() {
        pcglasso(
          S, lambda, c_parameter,
          Theta_start = solve(S),
          threshold = tol_pcglasso
        )
      })
      time_pcglasso_C_mat[m, i]    <- res$time
      res_val_pcglasso_C_mat[m, i] <- res$value

      # pcglasso_I
      res <- run_with_obj(function() {
        pcglasso(
          S, lambda, c_parameter,
          Theta_start = diag(nrow(S)),
          threshold   = tol_pcglasso
        )
      })
      time_pcglasso_I_mat[m, i]    <- res$time
      res_val_pcglasso_I_mat[m, i] <- res$value

      # pcglassoFast_I
      res <- run_with_obj(function() {
        pcglassoFast(
          S, lambda = lambda, alpha = alpha,
          R = diag(nrow(S)),
          tolerance = tolerance
        )$Sinv
      })
      time_pcglassoFast_I_mat[m, i]    <- res$time
      res_val_pcglassoFast_I_mat[m, i] <- res$value

      # pcglassoFast_C
      res <- run_with_obj(function() {
        pcglassoFast(
          S, lambda = lambda, alpha = alpha,
          R = cov2cor(solve(S)),
          tolerance = tolerance
        )$Sinv
      })
      time_pcglassoFast_C_mat[m, i]    <- res$time
      res_val_pcglassoFast_C_mat[m, i] <- res$value

      # pcglasso_cpp_I
      res <- run_with_obj(function() {
        pcglassoFast(
          S, lambda = lambda, alpha = alpha,
          solver_R = "cpp",
          R = diag(nrow(S)),
          tolerance = tolerance
        )$Sinv
      })
      time_pcglasso_cpp_I_mat[m, i]    <- res$time
      res_val_pcglasso_cpp_I_mat[m, i] <- res$value

      # pcglasso_cpp_C
      res <- run_with_obj(function() {
        pcglassoFast(
          S, lambda = lambda, alpha = alpha,
          solver_R = "cpp",
          R = cov2cor(solve(S)),
          tolerance = tolerance
        )$Sinv
      })
      time_pcglasso_cpp_C_mat[m, i]    <- res$time
      res_val_pcglasso_cpp_C_mat[m, i] <- res$value
    }
  }
  if (show_progress_bar) {
    close(pb)
  }

  trim <- if (M >= 10) 0.1 else 0
  mean_trim <- function(mat) apply(mat, 2, mean, trim = trim)

  time_pcglasso_C        <- mean_trim(time_pcglasso_C_mat)
  res_val_pcglasso_C     <- mean_trim(res_val_pcglasso_C_mat)
  time_pcglasso_I        <- mean_trim(time_pcglasso_I_mat)
  res_val_pcglasso_I     <- mean_trim(res_val_pcglasso_I_mat)
  time_pcglassoFast_I    <- mean_trim(time_pcglassoFast_I_mat)
  res_val_pcglassoFast_I <- mean_trim(res_val_pcglassoFast_I_mat)
  time_pcglassoFast_C    <- mean_trim(time_pcglassoFast_C_mat)
  res_val_pcglassoFast_C <- mean_trim(res_val_pcglassoFast_C_mat)
  time_pcglasso_cpp_I    <- mean_trim(time_pcglasso_cpp_I_mat)
  res_val_pcglasso_cpp_I <- mean_trim(res_val_pcglasso_cpp_I_mat)
  time_pcglasso_cpp_C    <- mean_trim(time_pcglasso_cpp_C_mat)
  res_val_pcglasso_cpp_C <- mean_trim(res_val_pcglasso_cpp_C_mat)

  df <- data.frame(
    time  = c(time_pcglasso_C, time_pcglasso_I,
              time_pcglassoFast_I, time_pcglassoFast_C,
              time_pcglasso_cpp_I, time_pcglasso_cpp_C),
    value = c(res_val_pcglasso_C, res_val_pcglasso_I,
              res_val_pcglassoFast_I, res_val_pcglassoFast_C,
              res_val_pcglasso_cpp_I, res_val_pcglasso_cpp_C),
    which = factor(c(
      rep("pcglasso start C", n_tol),
      rep("pcglasso start I", n_tol),
      rep("pcglassoFast start I", n_tol),
      rep("pcglassoFast start C", n_tol),
      rep("pcglasso_cpp start I", n_tol),
      rep("pcglasso_cpp start C", n_tol)
    )),
    tolerance = rep(tolerance_list, 6L),
    M         = M,
    p         = p,
    cor_modifier = cor_modifier,
    rho       = rho
  )

  df$alg  <- sub(" (start C|start I)$", "", df$which)
  df$init <- ifelse(grepl("start C", df$which), "C", "I")

  df
}

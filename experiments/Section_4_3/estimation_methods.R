# estimation_methods.R
#
# Wrapper functions for all estimators used in the simulation study.
# Each wrapper takes the same arguments and returns a named list:
#
#   list(
#     METHOD_bic = list(Q = <p×p matrix>, timing = <seconds>, alpha = <value or NA>),
#     METHOD_cv  = list(Q = <p×p matrix>, timing = <seconds>, alpha = <value or NA>)
#   )
#
# Model selection:
#   _bic  -- BIC_gamma (gamma = 0.5) on the full training sample
#   _cv   -- Gaussian log-likelihood on a held-out test set
#
# PC-GLasso variants (6):
#   estimator_pcglasso_I_dual     identity init, dual (Fortran) solver
#   estimator_pcglasso_I_primal   identity init, primal (C++) solver
#   estimator_pcglasso_C_dual     correlation init, dual solver
#   estimator_pcglasso_C_primal   correlation init, primal solver
#   estimator_pcglasso_I_Carter   identity init, Carter algorithm
#   estimator_pcglasso_C_Carter   correlation init, Carter algorithm
#
# Baselines (3):
#   estimator_glasso              Graphical LASSO
#   estimator_corglasso           Correlation-GLASSO
#   estimator_space               SPACE


# ---- PC-GLasso: shared core -------------------------------------------------

pcglasso_estimator_core <- function(
    S_full, S_train, S_test,
    n, n_train, n_test,
    lambdas, alpha_grid, pcglasso_tolerance,
    R0_full_fun  = NULL,
    R0_train_fun = NULL,
    solver_R     = c("dual", "primal")
) {
  solver_R <- match.arg(solver_R)

  # BIC selection on full data
  t_bic <- system.time({
    best_bic <- list(bic = Inf)
    for (a in alpha_grid) {
      R0   <- if (is.null(R0_full_fun)) diag(nrow(S_full)) else R0_full_fun()
      path <- pcglassoPath(
        S_full,
        alpha             = a,
        max_edge_fraction = 0.3,
        max_iter          = 1000,
        lambdas           = lambdas,
        tolerance         = pcglasso_tolerance,
        R0                = R0,
        solver_R          = solver_R
      )
      loss <- evaluate_objective_path(path, Sigma = S_full, n = n, gamma = 0.5)
      i    <- which.min(loss$BIC_gamma)
      if (loss$BIC_gamma[i] < best_bic$bic) {
        best_bic <- list(alpha = a, bic = loss$BIC_gamma[i], W = path$W_path[[i]])
      }
    }
    Q_bic <- (best_bic$W + t(best_bic$W)) / 2
  })

  # CV selection: fit on train, evaluate on test
  t_cv <- system.time({
    best_cv <- list(loglik = -Inf)
    for (a in alpha_grid) {
      R0   <- if (is.null(R0_train_fun)) diag(nrow(S_train)) else R0_train_fun()
      path <- pcglassoPath(
        S_train,
        alpha             = a,
        max_iter          = 2000,
        max_edge_fraction = 0.3,
        lambdas           = lambdas,
        tolerance         = pcglasso_tolerance,
        R0                = R0,
        solver_R          = solver_R
      )
      loss <- evaluate_objective_path(path, Sigma = S_test, n = n_test, gamma = 0.5)
      j    <- which.max(loss$loglik)
      if (loss$loglik[j] > best_cv$loglik) {
        best_cv <- list(alpha = a, loglik = loss$loglik[j], W = path$W_path[[j]])
      }
    }
    Q_cv <- (best_cv$W + t(best_cv$W)) / 2
  })

  list(
    PCGL_bic = list(Q = Q_bic, timing = as.numeric(t_bic["elapsed"]), alpha = best_bic$alpha),
    PCGL_cv  = list(Q = Q_cv,  timing = as.numeric(t_cv["elapsed"]),  alpha = best_cv$alpha)
  )
}


# ---- PC-GLasso variants -----------------------------------------------------

estimator_pcglasso_I_dual <- function(
    S_full, S_train, S_test, n, n_train, n_test,
    lambdas, alpha_grid, pcglasso_tolerance, ...
) {
  pcglasso_estimator_core(
    S_full, S_train, S_test, n, n_train, n_test,
    lambdas, alpha_grid, pcglasso_tolerance,
    R0_full_fun = NULL, R0_train_fun = NULL,
    solver_R = "dual"
  )
}

estimator_pcglasso_I_primal <- function(
    S_full, S_train, S_test, n, n_train, n_test,
    lambdas, alpha_grid, pcglasso_tolerance, ...
) {
  pcglasso_estimator_core(
    S_full, S_train, S_test, n, n_train, n_test,
    lambdas, alpha_grid, pcglasso_tolerance,
    R0_full_fun = NULL, R0_train_fun = NULL,
    solver_R = "primal"
  )
}

estimator_pcglasso_C_dual <- function(
    S_full, S_train, S_test, n, n_train, n_test,
    lambdas, alpha_grid, pcglasso_tolerance, ...
) {
  pcglasso_estimator_core(
    S_full, S_train, S_test, n, n_train, n_test,
    lambdas, alpha_grid, pcglasso_tolerance,
    R0_full_fun  = function() cov2cor(MASS::ginv(S_full)),
    R0_train_fun = function() cov2cor(MASS::ginv(S_train)),
    solver_R = "dual"
  )
}

estimator_pcglasso_C_primal <- function(
    S_full, S_train, S_test, n, n_train, n_test,
    lambdas, alpha_grid, pcglasso_tolerance, ...
) {
  pcglasso_estimator_core(
    S_full, S_train, S_test, n, n_train, n_test,
    lambdas, alpha_grid, pcglasso_tolerance,
    R0_full_fun  = function() cov2cor(MASS::ginv(S_full)),
    R0_train_fun = function() cov2cor(MASS::ginv(S_train)),
    solver_R = "primal"
  )
}


# ---- PC-GLasso (Carter algorithm) -------------------------------------------

# Shared path builder using Carter's (PCGLASSO package) solver
pcglasso_path_carter <- function(S, lambdas, tolerance,
                                 pcglasso_tolerance_modifier = 100,
                                 alpha = 1, R0 = NULL) {
  p               <- nrow(S)
  c_parameter     <- 1 - alpha
  precision_array <- array(0, dim = c(p, p, length(lambdas)))

  for (i in seq_along(lambdas)) {
    Theta_start <- if (i == 1) R0 else cov2cor(precision_array[, , i - 1])
    Theta_start <- (Theta_start + t(Theta_start)) / 2
    precision_array[, , i] <- pcglasso(
      S, lambdas[i],
      c           = c_parameter,
      threshold   = pcglasso_tolerance_modifier * tolerance,
      Theta_start = Theta_start
    )
  }

  list(
    lambdas              = lambdas,
    R_path               = lapply(seq_along(lambdas), function(i) cov2cor(precision_array[, , i])),
    Ri_path              = NA,
    D_path               = lapply(seq_along(lambdas), function(i) sqrt(diag(precision_array[, , i]))),
    W_path               = lapply(seq_along(lambdas), function(i) precision_array[, , i]),
    Wi_path              = NA,
    objective            = NA,
    iters                = NA,
    path_optimization_time = NA
  )
}

estimator_pcglasso_C_Carter <- function(
    S_full, S_train, S_test, n, n_train, n_test,
    lambdas, alpha_grid, pcglasso_tolerance, ...
) {
  t_bic <- system.time({
    best_bic <- list(bic = Inf)
    for (a in alpha_grid) {
      path <- pcglasso_path_carter(S_full, lambdas, pcglasso_tolerance,
                                   alpha = a, R0 = cov2cor(MASS::ginv(S_full)))
      loss <- evaluate_objective_path(path, Sigma = S_full, n = n, gamma = 0.5)
      i    <- which.min(loss$BIC_gamma)
      if (loss$BIC_gamma[i] < best_bic$bic)
        best_bic <- list(alpha = a, bic = loss$BIC_gamma[i], W = path$W_path[[i]])
    }
    Q_bic <- (best_bic$W + t(best_bic$W)) / 2
  })
  t_cv <- system.time({
    best_cv <- list(loglik = -Inf)
    for (a in alpha_grid) {
      path <- pcglasso_path_carter(S_train, lambdas, pcglasso_tolerance,
                                   alpha = a, R0 = cov2cor(MASS::ginv(S_train)))
      loss <- evaluate_objective_path(path, Sigma = S_test, n = n_test, gamma = 0.5)
      j    <- which.max(loss$loglik)
      if (loss$loglik[j] > best_cv$loglik)
        best_cv <- list(alpha = a, loglik = loss$loglik[j], W = path$W_path[[j]])
    }
    Q_cv <- (best_cv$W + t(best_cv$W)) / 2
  })
  list(
    PCGL_bic = list(Q = Q_bic, timing = as.numeric(t_bic["elapsed"]), alpha = best_bic$alpha),
    PCGL_cv  = list(Q = Q_cv,  timing = as.numeric(t_cv["elapsed"]),  alpha = best_cv$alpha)
  )
}

estimator_pcglasso_I_Carter <- function(
    S_full, S_train, S_test, n, n_train, n_test,
    lambdas, alpha_grid, pcglasso_tolerance, ...
) {
  t_bic <- system.time({
    best_bic <- list(bic = Inf)
    for (a in alpha_grid) {
      path <- pcglasso_path_carter(S_full, lambdas, pcglasso_tolerance,
                                   alpha = a, R0 = diag(nrow(S_full)))
      loss <- evaluate_objective_path(path, Sigma = S_full, n = n, gamma = 0.5)
      i    <- which.min(loss$BIC_gamma)
      if (loss$BIC_gamma[i] < best_bic$bic)
        best_bic <- list(alpha = a, bic = loss$BIC_gamma[i], W = path$W_path[[i]])
    }
    Q_bic <- (best_bic$W + t(best_bic$W)) / 2
  })
  t_cv <- system.time({
    best_cv <- list(loglik = -Inf)
    for (a in alpha_grid) {
      path <- pcglasso_path_carter(S_train, lambdas, pcglasso_tolerance,
                                   alpha = a, R0 = diag(nrow(S_train)))
      loss <- evaluate_objective_path(path, Sigma = S_test, n = n_test, gamma = 0.5)
      j    <- which.max(loss$loglik)
      if (loss$loglik[j] > best_cv$loglik)
        best_cv <- list(alpha = a, loglik = loss$loglik[j], W = path$W_path[[j]])
    }
    Q_cv <- (best_cv$W + t(best_cv$W)) / 2
  })
  list(
    PCGL_bic = list(Q = Q_bic, timing = as.numeric(t_bic["elapsed"]), alpha = best_bic$alpha),
    PCGL_cv  = list(Q = Q_cv,  timing = as.numeric(t_cv["elapsed"]),  alpha = best_cv$alpha)
  )
}


# ---- Graphical LASSO --------------------------------------------------------

estimator_glasso <- function(S_full, S_train, S_test, n, n_train, n_test,
                             lambdas, ...) {
  t_bic <- system.time({
    path_full   <- glasso::glassopath(S_full,  rholist = lambdas, penalize.diagonal = FALSE)
    loss_full   <- evaluate_objective_path(path_full$wi, Sigma = S_full, n = n, gamma = 0.5)
    i           <- which.min(loss_full$BIC_gamma)
    Q_bic       <- (path_full$wi[, , i] + t(path_full$wi[, , i])) / 2
  })
  t_cv <- system.time({
    path_train  <- glasso::glassopath(S_train, rholist = lambdas, penalize.diagonal = FALSE)
    loss_cv     <- evaluate_objective_path(path_train$wi, Sigma = S_test, n = n_test, gamma = 0.5)
    j           <- which.max(loss_cv$loglik)
    Q_cv        <- (path_train$wi[, , j] + t(path_train$wi[, , j])) / 2
  })
  list(
    GL_bic = list(Q = Q_bic, timing = as.numeric(t_bic["elapsed"]), alpha = NA),
    GL_cv  = list(Q = Q_cv,  timing = as.numeric(t_cv["elapsed"]),  alpha = NA)
  )
}


# ---- Correlation-GLASSO -----------------------------------------------------

estimator_corglasso <- function(S_full, S_train, S_test, n, n_train, n_test,
                                lambdas, ...) {
  t_bic <- system.time({
    C_full      <- cov2cor(S_full)
    path_full   <- glasso::glassopath(C_full,  rholist = lambdas, penalize.diagonal = FALSE)
    vars_full   <- diag(S_full)
    # Back-transform from correlation-scale precision to covariance-scale precision
    loss_full   <- evaluate_objective_path(
      cov2cor_inv(path_full$wi, 1 / vars_full),
      Sigma = S_full, n = n, gamma = 0.5
    )
    i           <- which.min(loss_full$BIC_gamma)
    Theta_bic   <- cov2cor_inv(path_full$wi[, , i], 1 / vars_full)
    Q_bic       <- (Theta_bic + t(Theta_bic)) / 2
  })
  t_cv <- system.time({
    C_train     <- cov2cor(S_train)
    path_train  <- glasso::glassopath(C_train, rholist = lambdas, penalize.diagonal = FALSE)
    vars_train  <- diag(S_train)
    loss_cv     <- evaluate_objective_path(
      cov2cor_inv(path_train$wi, 1 / vars_train),
      Sigma = S_test, n = n_test, gamma = 0.5
    )
    j           <- which.max(loss_cv$loglik)
    Theta_cv    <- cov2cor_inv(path_train$wi[, , j], 1 / vars_train)
    Q_cv        <- (Theta_cv + t(Theta_cv)) / 2
  })
  list(
    CorGL_bic = list(Q = Q_bic, timing = as.numeric(t_bic["elapsed"]), alpha = NA),
    CorGL_cv  = list(Q = Q_cv,  timing = as.numeric(t_cv["elapsed"]),  alpha = NA)
  )
}


# ---- SPACE ------------------------------------------------------------------

estimator_space <- function(S_full, S_train, S_test, n, n_train, n_test,
                            lambdas, data, train, test, ...) {
  # BIC selection on full data
  t_bic <- system.time({
    p <- ncol(S_full)
    l1_full     <- sqrt(n) * qnorm(1 - 0.1 / (2 * p^2))   # note: typo in package docs
    scale_full  <- exp(seq(log(2), log(0.8), length.out = length(lambdas)))
    path_full   <- array(0, dim = c(p, p, length(scale_full)))
    vars_full   <- diag(S_full)
    data_scaled <- as.matrix(scale(data))

    for (i in seq_along(scale_full)) {
      invisible(capture.output({
        sp <- space.joint(data_scaled, lam1 = l1_full * scale_full[i],
                          lam2 = 0, weight = 2, iter = 3)
      }))
      Theta <- -sp$ParCor
      diag(Theta) <- 1
      Theta <- cov2cor_inv(Theta, sp$sig.fit)
      Theta <- cov2cor_inv(Theta, 1 / vars_full)
      path_full[, , i] <- (Theta + t(Theta)) / 2
    }
    loss_full   <- evaluate_objective_path(path_full, Sigma = S_full, n = n, gamma = 0.5)
    i_bic       <- which.min(loss_full$BIC_gamma)
    Q_bic       <- path_full[, , i_bic]
    scale_bic   <- scale_full[i_bic]
  })

  # CV selection: fit on train, evaluate on test
  t_cv <- system.time({
    p <- ncol(S_full)
    l1_train    <- sqrt(n_train) * qnorm(1 - 0.1 / (2 * p^2))   # note: typo in package docs
    scale_train <- exp(seq(log(2), log(0.8), length.out = length(lambdas)))
    path_train  <- array(0, dim = c(p, p, length(scale_train)))
    vars_train  <- diag(S_train)
    train_scaled <- as.matrix(scale(train))

    for (i in seq_along(scale_train)) {
      invisible(capture.output({
        sp <- space.joint(train_scaled, lam1 = l1_train * scale_train[i],
                          lam2 = 0, weight = 2, iter = 3)
      }))
      Theta <- -sp$ParCor
      diag(Theta) <- 1
      Theta <- cov2cor_inv(Theta, sp$sig.fit)
      Theta <- cov2cor_inv(Theta, 1 / vars_train)
      path_train[, , i] <- (Theta + t(Theta)) / 2
    }
    loss_cv     <- evaluate_objective_path(path_train, Sigma = S_test, n = nrow(test), gamma = 0.5)
    j_cv        <- which.max(loss_cv$loglik)
    Q_cv        <- path_train[, , j_cv]
    scale_cv    <- scale_train[j_cv]
  })

  list(
    Space_bic = list(Q = Q_bic, timing = as.numeric(t_bic["elapsed"]),
                     alpha = NA, scale = scale_bic),
    Space_cv  = list(Q = Q_cv,  timing = as.numeric(t_cv["elapsed"]),
                     alpha = NA, scale = scale_cv)
  )
}

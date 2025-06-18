# ---- PC-GLasso Estimator ----
estimator_pcglasso <- function(S_full, S_train, S_test, n, n_train, n_test, lambdas, alpha.grid, ...) {
  # BIC selection
  t_bic <- system.time({
    best_bic <- list(bic = Inf)
    for (a in alpha.grid) {
      path <- pcglassoPath(S_full, alpha = a, max.edge.fraction = 0.3,
                           lambda.min.ratio = min(lambdas) / max(lambdas), nlambda = length(lambdas))
      loss <- loss.evaluation(path, Sigma = S_full, n = n, gamma = 0)
      i <- which.min(loss$BIC_gamma)
      if (loss$BIC_gamma[i] < best_bic$bic) {
        best_bic <- list(alpha = a, lambda = path$lambda[i], bic = loss$BIC_gamma[i], W = path$W_path[[i]])
      }
    }
    Q_pc_bic <- (best_bic$W + t(best_bic$W))/2
  })
  # CV selection
  t_cv <- system.time({
    best_cv <- list(loglik = -Inf)
    for (a in alpha.grid) {
      path <- pcglassoPath(S_train, alpha = a, max.edge.fraction = 0.3,
                           lambda.min.ratio = min(lambdas) / max(lambdas), nlambda = length(lambdas))
      loss <- loss.evaluation(path, Sigma = S_test, n = n_test, gamma = 0)
      j <- which.max(loss$loglik)
      if (loss$loglik[j] > best_cv$loglik) {
        best_cv <- list(alpha = a, lambda = path$lambda[j], loglik = loss$loglik[j], W = path$W_path[[j]])
      }
    }
    Q_pc_cv <- (best_cv$W + t(best_cv$W))/2
  })
  list(
    PCGL_bic = list(Q = Q_pc_bic, timing = as.numeric(t_bic["elapsed"]), alpha = best_bic$alpha),
    PCGL_cv  = list(Q = Q_pc_cv,  timing = as.numeric(t_cv["elapsed"]),  alpha = best_cv$alpha)
  )
}

# ---- Glasso Estimator ----
estimator_glasso <- function(S_full, S_train, S_test, n, n_train, n_test, lambdas, ...) {
  t_bic <- system.time({
    gl_full_path <- glasso::glassopath(S_full, rholist = lambdas, penalize.diagonal = FALSE)
    loss_gl_full <- loss.evaluation(gl_full_path$wi, Sigma = S_full, n = n, gamma = 0.)
    idx_gl_bic   <- which.min(loss_gl_full$BIC)
    Q_gl_bic     <- (gl_full_path$wi[,,idx_gl_bic] + t(gl_full_path$wi[,,idx_gl_bic])) / 2
  })
  t_cv <- system.time({
    gl_tr_path   <- glasso::glassopath(S_train, rholist = lambdas, penalize.diagonal = FALSE)
    loss_gl_cv   <- loss.evaluation(gl_tr_path$wi, Sigma = S_test, n = n_test, gamma = 0.)
    idx_gl_cv    <- which.max(loss_gl_cv$loglik)
    Q_gl_cv      <- (gl_tr_path$wi[,,idx_gl_cv] + t(gl_tr_path$wi[,,idx_gl_cv])) / 2
  })
  list(
    GL_bic = list(Q = Q_gl_bic, timing = as.numeric(t_bic["elapsed"]), alpha = NA),
    GL_cv  = list(Q = Q_gl_cv,  timing = as.numeric(t_cv["elapsed"]),  alpha = NA)
  )
}

# ---- Correlation-Glasso Estimator ----
estimator_corglasso <- function(S_full, S_train, S_test, n, n_train, n_test, lambdas, ...) {
  # BIC
  t_bic <- system.time({
    C_full       <- cov2cor(S_full)
    cg_full_path <- glasso::glassopath(C_full, rholist = lambdas, penalize.diagonal = FALSE)
    vars_full    <- diag(S_full)
    loss_cg_full <- loss.evaluation(cov2cor.inv(cg_full_path$wi, 1/vars_full), Sigma = S_full,
                                    n = n, gamma = 0.)
    idx_cg_bic   <- which.min(loss_cg_full$BIC)
    Theta_cg_bic <- cov2cor.inv(cg_full_path$wi[,,idx_cg_bic], 1/vars_full)
    Q_cg_bic     <- (Theta_cg_bic + t(Theta_cg_bic)) / 2
  })
  # CV
  t_cv <- system.time({
    C_tr         <- cov2cor(S_train)
    cg_tr_path   <- glasso::glassopath(C_tr, rholist = lambdas, penalize.diagonal = FALSE)
    vars_tr      <- diag(S_train)
    loss_cg_cv   <- loss.evaluation(cov2cor.inv(cg_tr_path$wi, 1/vars_tr), Sigma = S_test,
                                    n = n_test, gamma = 0.)
    idx_cg_cv    <- which.max(loss_cg_cv$loglik)
    Theta_cg_cv  <- cov2cor.inv(cg_tr_path$wi[,,idx_cg_cv], 1/vars_tr)
    Q_cg_cv      <- (Theta_cg_cv + t(Theta_cg_cv)) / 2
  })
  list(
    CorGL_bic = list(Q = Q_cg_bic, timing = as.numeric(t_bic["elapsed"]), alpha = NA),
    CorGL_cv  = list(Q = Q_cg_cv,  timing = as.numeric(t_cv["elapsed"]),  alpha = NA)
  )
}

# ---- SPACE Estimator ----
estimator_space <- function(S_full, S_train, S_test, n, n_train, n_test, lambdas, data, train, test, ...) {
  # full-data BIC selection
  t_bic <- system.time({
    p <- ncol(S_full)
    l1_full    <- sqrt(n) * qnorm(1 - 0.1/(2 * p^2)) #typo in package
    scale_full <- exp(seq(log(2), log(0.8), length.out = length(lambdas)))
    res_space_f <- array(0, dim = c(p, p, length(scale_full)))
    data <- as.matrix(scale(data))

    vars_full    <- diag(S_full)
    for (i in seq_along(scale_full)) {
      invisible(
        capture.output({
      sp <- space.joint(as.matrix(scale(data)),
                        lam1 = l1_full * scale_full[i],
                        lam2 = 0,  weight=2,iter = 3)}))

      Theta <- -sp$ParCor
      diag(Theta) <- 1
      Theta <- cov2cor.inv(Theta ,sp$sig.fit)
      Theta <- cov2cor.inv(Theta, 1/vars_full)
      res_space_f[,,i] <- (Theta + t(Theta)) / 2
    }
    loss_space_full <- loss.evaluation(res_space_f, Sigma = S_full, n = n, gamma = 0)
    idx_space_bic   <- which.min(loss_space_full$BIC_gamma)
    Q_space_bic     <- res_space_f[,, idx_space_bic]
    scale_bic <- scale_full[idx_space_bic]
  })
  # train-test CV selection
  t_cv <- system.time({
    p <- ncol(S_train)
    l1_tr    <- sqrt(n_train) * qnorm(1 - 0.1/(2 * p^2)) #typo in package
    scale_tr <- exp(seq(log(2), log(0.8), length.out = length(lambdas)))
    res_space_t <- array(0, dim = c(p, p, length(scale_tr)))
    vars_full    <- diag(S_train)
    train <- as.matrix(scale(train))
    for (i in seq_along(scale_tr)) {
      invisible(
        capture.output({
      sp <- space.joint(train,
                        lam1 =  l1_tr * scale_tr[i],
                        weight=2,
                        lam2 = 0, iter = 3)}))
      #iter = 3 used in article
      Theta <- -sp$ParCor
      diag(Theta) <- 1
      Theta <- cov2cor.inv(Theta ,sp$sig.fit)
      Theta <- cov2cor.inv(Theta, 1/vars_full)
      res_space_t[,,i] <- (Theta + t(Theta)) / 2
    }
    loss_space_cv <- loss.evaluation(res_space_t, Sigma = S_test, n = nrow(test), gamma = 0)
    idx_space_cv  <- which.max(loss_space_cv$loglik)
    Q_space_cv    <- res_space_t[,, idx_space_cv]
    scale_cv <- scale_tr[idx_space_cv]
  })
  list(
    Space_bic  = list(Q = Q_space_bic,
                      timing = as.numeric(t_bic["elapsed"]),
                      alpha = NA,
                      scale = scale_bic),
    Space_cv   = list(Q = Q_space_cv,
                      timing = as.numeric(t_cv["elapsed"]),
                      alpha = NA,
                      scale = scale_cv)
  )
}

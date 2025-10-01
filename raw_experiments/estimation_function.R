library(pcglassoFast)
library(glasso)
library(space) # remotes::install_version("space", version = "0.1-1.1")
estimator_space <- function(S_full, n, lambdas, data, gamma = 0, min_scale= log(0.7), max_scale =log(4)) {
  t_full <- system.time({
    p <- ncol(S_full)
    l1_full    <- sqrt(n) * qnorm(1 - 0.1 / (2 * p^2)) # typo in paper/package
    scale_full <- exp(seq(max_scale, min_scale, length.out = length(lambdas)))
    res_space  <- array(0, dim = c(p, p, length(scale_full)))
    vars_full  <- diag(S_full)
    data       <- as.matrix(scale(data))

    for (i in seq_along(scale_full)) {
      invisible(
        capture.output({
          sp <- space.joint(data,
                            lam1   = l1_full * scale_full[i],
                            lam2   = 0,
                            weight = 2,
                            iter   = 3)
        })
      )
      Theta <- -sp$ParCor
      diag(Theta) <- 1
      Theta <- cov2cor_inv(Theta, sp$sig.fit)
      Theta <- cov2cor_inv(Theta, 1 / vars_full)
      res_space[,,i] <- (Theta + t(Theta)) / 2
    }

    loss_space <- evaluate_objective_path(res_space, Sigma = S_full, n = n, gamma = gamma)
  })

  list(
    path       = res_space,
    path.loss  = loss_space,
    timing     = as.numeric(t_full["elapsed"])
  )
}

estimator_corglasso <- function(S_full, n, lambdas, gamma =0) {
  t_full <- system.time({
    # Convert covariance to correlation matrix
    C_full <- cov2cor(S_full)
    invisible(
      capture.output({
    # Estimate correlation glasso path
    cg_full_path <- glasso::glassopath(C_full,
                                       rholist = lambdas,
                                       penalize.diagonal = FALSE)}))

    # Rescale to precision matrix scale
    vars_full <- diag(S_full)
    cg_prec_path <- cov2cor_inv(cg_full_path$wi, 1 / vars_full)

    # Evaluate loss
    loss_path <- evaluate_objective_path(cg_prec_path, Sigma = S_full, n = n, gamma = gamma)
  })

  list(
    path       = cg_prec_path,
    path.scaled = cg_full_path,
    path.loss  = loss_path,
    timing     = as.numeric(t_full["elapsed"])
  )
}

# Hub-Graphical Lasso (weighted GLASSO) on correlation scale,
# rescaled back to covariance-scale precision for evaluation.
#
# Inputs
# - S_full:  sample covariance matrix (p x p)
# - n:       sample size (passed to glasso::glasso as nobs)
# - lambdas: numeric vector of lambda values for the path
# - gamma:   extra term for your evaluate_objective_path()
# - eps_W:   ridge used in W computation; solve(C + eps_W*I)
# - penalize.diagonal, thr, maxit: passed to glasso::glasso
#
# Requirements in your env:
# - cov2cor_inv(wi_corr_array, inv_vars) -> covariance-scale precision array
# - evaluate_objective_path(prec_path, Sigma = S_full, n = n, gamma = gamma)

estimator_hubcorglasso <- function(
    S_full, n, lambdas, gamma = 0,
    eps_W = 0.1,
    penalize.diagonal = FALSE,
    thr = 1e-4,
    maxit = 1e4
) {
  stopifnot(is.matrix(S_full), nrow(S_full) == ncol(S_full))
  p <- ncol(S_full)
  if (length(lambdas) < 1L) stop("`lambdas` must be a non-empty numeric vector.")

  t_full <- system.time({
    # --- 1) Covariance -> Correlation ---------------------------------------
    C_full <- cov2cor(S_full)
    inv_mat <- solve(C_full + eps_W * diag(p))


    # Off-diagonal row sums (L1) for scaling
    ai <- rowSums(abs(inv_mat)) - abs(diag(inv_mat))
    W  <- matrix(0, p, p)
    # Fill symmetric weights; guard against zero denominators
    for (i in 2:p) {
      for (j in 1:(i - 1)) {
        aij <- abs(inv_mat[i, j])
        wij <- if (aij > 0 && ai[i] > 0 && ai[j] > 0) 1 / (aij * ai[i] * ai[j]) else 0
        W[i, j] <- wij
        W[j, i] <- wij
      }
    }
    diag(W) <- 0
    W <- W/mean(W[W!=0]) # normalize to mean 1 for interpretability
    # --- 3) Fit weighted GLASSO path on correlation scale --------------------
    K <- length(lambdas)
    OUTPUT <- vector("list", K)

    # suppress console chatter
    invisible(capture.output({
      for (k in seq_len(K)) {
        if (k == 1L) {
          OUTPUT[[k]] <- glasso::glasso(
            s = C_full, rho = lambdas[k] * W,
            nobs = n, thr = thr, maxit = maxit, approx = FALSE,
            penalize.diagonal = penalize.diagonal,
            start = "cold", w.init = NULL, wi.init = NULL,
            trace = FALSE
          )
        } else {
          OUTPUT[[k]] <- glasso::glasso(
            s = C_full, rho = lambdas[k] * W,
            nobs = n, thr = thr, maxit = maxit, approx = FALSE,
            penalize.diagonal = penalize.diagonal,
            start = "warm",
            w.init  = OUTPUT[[k - 1]]$w,
            wi.init = OUTPUT[[k - 1]]$wi,
            trace = FALSE
          )
        }
      }
    }))

    # --- 4) Collect precision path on correlation scale into an array --------
    wi_corr <- array(NA_real_, dim = c(p, p, K))
    for (k in seq_len(K)) wi_corr[, , k] <- OUTPUT[[k]]$wi

    # --- 5) Rescale to covariance-scale precision ----------------------------
    vars_full    <- diag(S_full)
    inv_vars_vec <- 1 / vars_full
    hg_prec_path <- cov2cor_inv(wi_corr, inv_vars_vec)

    # --- 6) Evaluate loss on covariance scale --------------------------------
    loss_path <- evaluate_objective_path(hg_prec_path, Sigma = S_full, n = n, gamma = gamma)
  })

  list(
    path         = hg_prec_path,        # precision on covariance scale (p x p x K)
    path.scaled  = wi_corr,             # precision on correlation scale (p x p x K)
    path.loss    = loss_path,           # numeric length K
    W            = W,                   # hub weights used
    fits         = OUTPUT,              # raw glasso fits (warm-start chain)
    timing       = as.numeric(t_full["elapsed"])
  )
}


estimator_pcglasso <- function(S_full,
                               n,
                               lambdas,
                               alpha_grid = 0,
                               gamma = 0,
                               max_edge_fraction = 0.3,
                               R_start = NULL,
                               verbose = 0) {
  t_full <- system.time({
    pc_path_list  <- list()
    pc_loss_list  <- list()
    pc_path_list_all <- list()
    if(is.null(R_start))
    {
      R_start <- diag(nrow(S_full))
    }

    for (a in alpha_grid) {
      path <- pcglassoPath(
        S_full,
        alpha = a,
        max_edge_fraction = max_edge_fraction,
        lambdas = lambdas,
        R0 = R_start,
        verbose = verbose
      )

      p <- nrow(path$W_path[[1]])
      K <- length(path$W_path)

      # Preallocate 3D array
      W <- array(0, dim = c(p, p, K))

      # Fill the array
      for (k in seq_len(K)) {
        W[,,k] <- path$W_path[[k]]
      }
      pc_path_list_all[[as.character(a)]]  <- path
      pc_path_list[[as.character(a)]] <- W
      pc_loss_list[[as.character(a)]] <- evaluate_objective_path(path, Sigma = S_full, n = n, gamma = gamma)
    }
  })
  if(length(pc_path_list) ==1)
  {
    return(list(
      path       = W,
      path.all   = pc_path_list_all[[1]],
      path.loss  = pc_loss_list[[1]],
      timing     = as.numeric(t_full["elapsed"])
    ))
  }
  list(
    path       = pc_path_list,
    path.all   = pc_path_list_all,
    path.loss  = pc_loss_list,
    alpha_grid = alpha_grid,
    timing     = as.numeric(t_full["elapsed"])
  )
}


estimator_glasso <- function(S_full, n, lambdas, gamma = 0) {
  t_full <- system.time({
    invisible(
      capture.output({
    gl_full_path <- glasso::glassopath(S_full, rholist = lambdas, penalize.diagonal = FALSE)
    }))
    loss_gl_full <- evaluate_objective_path(gl_full_path$wi, Sigma = S_full, n = n, gamma = gamma)
  })
  list(
    path       = gl_full_path$wi,
    path.loss  = loss_gl_full,
    timing     = as.numeric(t_full["elapsed"])
  )
}

#' find the optimal array or the array closest to number of edges
#' get alpha values from the precision matrix
#' @param path 3D array
#' @param path.loss EBIC, nedges
get_optimal_matrix <- function(path, path.loss, criterion = "BIC_gamma", max_edges = -1) {
  if (max_edges == -1) {
    idx_opt <- which.min(path.loss[[criterion]])
  } else {
    idx_opt <- which.min(abs(path.loss$nEdges - max_edges))
  }
  return(list(
    idx = idx_opt,
    Theta_opt = path[,,idx_opt],
    nEdges = path.loss$nEdges[idx_opt],
    criterion_value = path.loss[[criterion]][idx_opt],
    alpha = get_alpha(path[,,idx_opt])
  ))
}


#' get the alpha values from a precision matrix
#' @param Theta precision matrix
#' @return vector of alpha values
get_alpha <- function(Theta,scale=1) {
  return(apply(cov2cor(Theta), 2, function(x) {
    sum(abs(x)^scale) - 1
  }))
}

make_plot_matrix <- function(my_matrix, my_title) {
  matrix_data <- my_matrix != 0
  df_matrix <- as.data.frame(as.table(matrix_data))
  colnames(df_matrix) <- c("Row", "Column", "Value")

  df_matrix$Row    <- as.numeric(df_matrix$Row)
  df_matrix$Column <- as.numeric(df_matrix$Column)
  df_matrix$Value  <- as.numeric(df_matrix$Value)

  nnz <- sum(matrix_data)

  ggplot(df_matrix, aes(x = Column, y = Row, fill = Value)) +
    geom_tile(color = "white") +
    scale_fill_gradient(low = "white", high = "blue", name = "Non-Zero") +
    labs(
      title = paste(my_title, ", nnz =", nnz),
      x = NULL,
      y = NULL
    ) +
    scale_x_continuous(breaks = seq(0, ncol(my_matrix), by = 20)) +
    scale_y_reverse(breaks = seq(0, nrow(my_matrix), by = 20)) +  # Reverse to match matrix layout
    coord_fixed() +  # Keep aspect ratio 1:1
    theme_minimal(base_size = 12) +
    theme(
      panel.grid       = element_blank(),
      axis.ticks       = element_line(),
      legend.position  = "none",
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      plot.title       = element_text(hjust = 0.5),  # center title
      plot.title.position = "plot"
    )
}

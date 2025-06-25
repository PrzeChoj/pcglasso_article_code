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
      Theta <- cov2cor.inv(Theta, sp$sig.fit)
      Theta <- cov2cor.inv(Theta, 1 / vars_full)
      res_space[,,i] <- (Theta + t(Theta)) / 2
    }

    loss_space <- loss.evaluation(res_space, Sigma = S_full, n = n, gamma = gamma)
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
    cg_prec_path <- cov2cor.inv(cg_full_path$wi, 1 / vars_full)

    # Evaluate loss
    loss_path <- loss.evaluation(cg_prec_path, Sigma = S_full, n = n, gamma = gamma)
  })

  list(
    path       = cg_prec_path,
    path.scaled = cg_full_path,
    path.loss  = loss_path,
    timing     = as.numeric(t_full["elapsed"])
  )
}

estimator_pcglasso <- function(S_full, n, lambdas, alpha.grid = 0, gamma = 0,max.edge.fraction =0.3) {
  t_full <- system.time({
    pc_path_list  <- list()
    pc_loss_list  <- list()
    pc_path_list_all <- list()
    for (a in alpha.grid) {
      path <- pcglassoPath(S_full,
                           alpha             = a,
                           max.edge.fraction = max.edge.fraction,
                           lambdas = lambdas)

      p <- nrow(path$W[[1]])
      K <- length(path$W)

      # Preallocate 3D array
      W <- array(NA, dim = c(p, p, K))

      # Fill the array
      for (k in seq_len(K)) {
        W[,,k] <- path$W[[k]]
      }
      pc_path_list_all[[as.character(a)]]  <- path
      pc_path_list[[as.character(a)]] <- W
      pc_loss_list[[as.character(a)]] <- loss.evaluation(path, Sigma = S_full, n = n, gamma = gamma)
    }
  })
  if(length(pc_path_list) ==1)
  {
    return(list(
      path       = pc_path_list[[1]],
      path.all   =pc_path_list_all[[1]],
      path.loss  = pc_loss_list[[1]],
      timing     = as.numeric(t_full["elapsed"])
    ))
  }
  list(
    path       = pc_path_list,
    path.all   = pc_path_list_all,
    path.loss  = pc_loss_list,
    alpha.grid = alpha.grid,
    timing     = as.numeric(t_full["elapsed"])
  )
}


estimator_glasso <- function(S_full, n, lambdas, gamma = 0) {
  t_full <- system.time({
    invisible(
      capture.output({
    gl_full_path <- glasso::glassopath(S_full, rholist = lambdas, penalize.diagonal = FALSE)
    }))
    loss_gl_full <- loss.evaluation(gl_full_path$wi, Sigma = S_full, n = n, gamma = gamma)
  })
  list(
    path       = gl_full_path$wi,
    path.loss  = loss_gl_full,
    timing     = as.numeric(t_full["elapsed"])
  )
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

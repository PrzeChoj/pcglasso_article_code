library(pcglassoFast)
library(glasso)
library(space) # remotes::install_version("space", version = "0.1-1.1")

estimator_space <- function(S_full, n, lambdas, data, gamma = 0, min_scale= log(0.7), max_scale =log(4)) {
  t_full <- system.time({
    p <- ncol(S_full)
    l1_full    <- sqrt(n) * qnorm(1 - 0.1 / (2 * p^2)) # typo in SPACE paper/package
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

estimator_pcglassoFast <- function(
    S_full,
    n,
    lambdas,
    solver_R = c("primal", "dual"),
    alpha_grid = 0,
    gamma = 0,
    max_edge_fraction = 0.3,
    R_start = NULL,
    max_iter = 500,
    verbose = 0) {
  solver_R <- match.arg(solver_R)

  t_full <- system.time({
    pc_path_list  <- list()
    pc_loss_list  <- list()
    pc_path_list_all <- list()
    if(is.null(R_start))
    {
      R_start <- cov2cor(MASS::ginv(S_full))
    }

    for (a in alpha_grid) {
      path <- pcglassoPath(
        S_full,
        alpha = a,
        max_edge_fraction = max_edge_fraction,
        lambdas = lambdas,
        solver_R = solver_R,
        R0 = R_start,
        max_iter = max_iter,
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

estimator_pcglasso <- function(S_full,
                               n,
                               lambdas,
                               alpha_grid = 0,
                               gamma = 0,
                               max_edge_fraction = 0.3,
                               R_start = NULL,
                               max_iter = 500,
                               method = 'primal',
                               verbose = 0) {
  estimator_pcglassoFast(
    S_full = S_full,
    n =n,
    lambdas=  lambdas,
    solver_R = method,
    alpha_grid = alpha_grid,
    gamma = gamma,
    max_edge_fraction= max_edge_fraction,
    R_start=R_start,
    max_iter = max_iter,
    verbose = verbose
  )
}

estimator_pcglasso_primal <- function(
    S_full,
    n,
    lambdas,
    alpha_grid = 0,
    gamma = 0,
    max_edge_fraction = 0.3,
    R_start = NULL,
    verbose = 0) {
  estimator_pcglassoFast(
    S_full,
    n,
    lambdas,
    "primal",
    alpha_grid,
    gamma,
    max_edge_fraction,
    R_start,
    verbose
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
make_plot_matrix <- function(my_matrix, my_title,
                             x_lab = "Column", y_lab = "Row",
                             base_size = 6,         # overall baseline
                             title_size = 8,
                             axis_title_size = 6,
                             axis_text_size = 5,
                             tick_length_pt = 1) {

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
      title = paste(my_title, ", non-zero = ", round(100*(nnz-dim(matrix_data)[1])/(dim(matrix_data)[1]^2 - dim(matrix_data)[1]),0),'%', sep = ""),
      x = x_lab,
      y = y_lab
    ) +
    scale_x_continuous(breaks = seq(0, ncol(my_matrix), by = 20)) +
    scale_y_reverse(breaks = seq(0, nrow(my_matrix), by = 20)) +
    coord_fixed() +
    theme_minimal(base_size = base_size) +
    theme(
      panel.grid          = element_blank(),
      axis.ticks          = element_line(linewidth = 0.2),
      axis.ticks.length   = grid::unit(tick_length_pt, "pt"),
      legend.position     = "none",
      panel.background    = element_rect(fill = "white", color = NA),
      plot.background     = element_rect(fill = "white", color = NA),
      plot.title          = element_text(size = title_size, hjust = 0.5, margin = margin(b = 2)),
      axis.title.x        = element_text(size = axis_title_size, margin = margin(t = 2)),
      axis.title.y        = element_text(size = axis_title_size, margin = margin(r = 4)),
      axis.text.x         = element_text(size = axis_text_size),
      axis.text.y         = element_text(size = axis_text_size),
      plot.title.position = "plot"
    )
}
make_plot_matrix_v2 <- function(my_matrix, my_title,
                                x_lab = "Column", y_lab = "Row",
                                base_size = 6,
                                title_size = 8,
                                axis_title_size = 6,
                                axis_text_size = 5,
                                tick_length_pt = 1,
                                highlight_col = NULL,
                                highlight_label = NULL,
                                diag_color = "grey70",
                                low_color = "white",
                                high_color = "blue",
                                legend_title = "|value|",
                                show_legend = TRUE) {

  stopifnot(is.matrix(my_matrix))

  nr <- nrow(my_matrix)
  nc <- ncol(my_matrix)

  # Build plotting data
  df_matrix <- expand.grid(
    Row = seq_len(nr),
    Column = seq_len(nc)
  )
  df_matrix$Value <- as.vector(abs(my_matrix))   # magnitude
  df_matrix$IsDiagonal <- df_matrix$Row == df_matrix$Column

  # Off-diagonal non-zero percentage
  off_diag_idx <- !df_matrix$IsDiagonal
  off_diag_nnz <- sum(df_matrix$Value[off_diag_idx] != 0)
  off_diag_tot <- sum(off_diag_idx)
  nnz_pct <- if (off_diag_tot > 0) {
    round(100 * off_diag_nnz / off_diag_tot, 0)
  } else {
    NA
  }

  # Avoid degenerate scale if everything is zero
  max_mag <- max(df_matrix$Value[off_diag_idx], na.rm = TRUE)
  if (!is.finite(max_mag) || max_mag == 0) {
    max_mag <- 1
  }

  # Nice axis breaks
  x_breaks <- unique(c(1, seq(20, nc, by = 20), nc))
  y_breaks <- unique(c(1, seq(20, nr, by = 20), nr))

  p <- ggplot(df_matrix, aes(x = Column, y = Row)) +
    # First layer: magnitude heatmap
    geom_tile(aes(fill = Value), color = "white", linewidth = 0.1) +
    # Second layer: overwrite diagonal in gray
    geom_tile(
      data = subset(df_matrix, IsDiagonal),
      fill = diag_color,
      color = "white",
      linewidth = 0.1
    ) +
    scale_fill_gradient(
      low = low_color,
      high = high_color,
      limits = c(0, max_mag),
      name = legend_title
    ) +
    labs(
      title = paste0(my_title, ", non-zero = ", nnz_pct, "%"),
      x = x_lab,
      y = y_lab
    ) +
    scale_x_continuous(
      limits = c(0.5, nc + 0.5),
      breaks = x_breaks,
      expand = c(0, 0)
    ) +
    scale_y_reverse(
      limits = c(nr + 0.5, 0.5),
      breaks = y_breaks,
      expand = c(0, 0)
    ) +
    coord_fixed(clip = "off") +
    theme_minimal(base_size = base_size) +
    theme(
      panel.grid          = element_blank(),
      axis.ticks          = element_line(linewidth = 0.2),
      axis.ticks.length   = grid::unit(tick_length_pt, "pt"),
      legend.position     = if (show_legend) "right" else "none",
      panel.background    = element_rect(fill = "white", color = NA),
      plot.background     = element_rect(fill = "white", color = NA),
      plot.title          = element_text(size = title_size, hjust = 0.5,
                                         margin = margin(b = 2)),
      axis.title.x        = element_text(size = axis_title_size,
                                         margin = margin(t = 2)),
      axis.title.y        = element_text(size = axis_title_size,
                                         margin = margin(r = 4)),
      axis.text.x         = element_text(size = axis_text_size),
      axis.text.y         = element_text(size = axis_text_size),
      plot.title.position = "plot",
      plot.margin         = margin(14, 5, 5, 5)
    )

  # Optional arrow + label pointing to a chosen column
  if (!is.null(highlight_col)) {
    if (highlight_col < 1 || highlight_col > nc) {
      stop("highlight_col must be between 1 and ncol(my_matrix).")
    }

    if (is.null(highlight_label)) {
      highlight_label <- paste("Column", highlight_col)
    }

    p <- p +
      annotate(
        "segment",
        x = highlight_col, xend = highlight_col,
        y = 0.15, yend = 0.95,
        linewidth = 0.3,
        arrow = grid::arrow(length = grid::unit(1.5, "mm"), type = "closed")
      ) +
      annotate(
        "text",
        x = highlight_col,
        y = -0.05,
        label = highlight_label,
        fontface = "bold",
        vjust = 1,
        size = axis_text_size / 2.8
      )
  }

  return(p)
}
# Helper to make a single sorted-alpha plot

make_alpha_plot <- function(alpha, title, ylims = NULL,
                            latex = FALSE,
                            names = NULL,
                            print.number.names = NULL,
                            decreasing = FALSE,
                            ylab = expression(alpha)) {
  a <- alpha[is.finite(alpha)]
  stopifnot(length(a) > 0)

  # sort (ascending by default)
  ord <- order(a, decreasing = decreasing)
  a_s  <- a[ord]
  nm_s <- if (!is.null(names)) names[ord] else
    if (!is.null(names(a))) names(a)[ord] else
      as.character(ord)

  df <- data.frame(rank = seq_along(a_s), value = a_s, name = nm_s)

  # Title handling (optional LaTeX)
  title_label <- title
  if (latex) {
    if (!requireNamespace("latex2exp", quietly = TRUE))
      stop("Install 'latex2exp' to use latex=TRUE.")
    title_label <- latex2exp::TeX(title)
  }

  # Base plot
  p <- ggplot(df, aes(x = rank, y = value)) +
    geom_point(alpha = 0.6, size = 0.9) +
    geom_line(linewidth = 0.3) +
    labs(title = title_label, x = "Rank", y = ylab) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title.position = "panel",
      plot.title = element_text(hjust = 0.5),
      panel.grid.minor = element_blank()
    )

  if (!is.null(ylims)) p <- p + scale_y_continuous(limits = ylims)

  # Optional annotation of the top-K largest values with arrows + names
  if (!is.null(print.number.names) && print.number.names > 0) {
    if (!requireNamespace("ggrepel", quietly = TRUE))
      stop("Install 'ggrepel' to use print.number.names.")
    k <- min(print.number.names, nrow(df))

    # rows to annotate = largest values regardless of sorting direction
    sel <- if (!decreasing) (nrow(df) - k + 1):nrow(df) else 1:k
    labdf <- df[sel, , drop = FALSE]

    # push labels to the right; allow drawing outside panel
    nudge <- -0.2 * nrow(df)  # 6% of x-range
    p <- p +
      ggrepel::geom_text_repel(
        data = labdf,
        aes(label = name),
        nudge_x = nudge,
        direction = "y",
        min.segment.length = 0,
        box.padding = 0.3,
        segment.alpha = 0.7,
        segment.color = "grey40",
        arrow = grid::arrow(length = grid::unit(0.015, "npc"))
      ) +
      coord_cartesian(clip = "off") +
      theme(plot.margin = margin(5.5, 40, 5.5, 5.5))  # room for labels on right
  }

  p
}
# Build a grid of alpha plots from a named list
make_alpha_grid <- function(alpha_list, ncol = 2, common_y = TRUE) {
  stopifnot(length(alpha_list) > 0)
  if (is.null(names(alpha_list))) {
    names(alpha_list) <- paste0("Series ", seq_along(alpha_list))
  }
  ylims <- if (common_y) range(unlist(alpha_list), finite = TRUE, na.rm = TRUE) else NULL

  plots <- lapply(names(alpha_list), function(nm) {
    make_alpha_plot(alpha_list[[nm]], nm, ylims)
  })
  # Arrange as a grid with patchwork
  wrap_plots(plots, ncol = ncol)
}

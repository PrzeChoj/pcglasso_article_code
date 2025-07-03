library(parallel)
library(pbmcapply)
library(space)
source('estimation_methods.R')

run_single <- function(Q, n, split_train = 0.7,
                       alpha_grid = sort(unique(c(seq(-0.1, 0.1, length.out = 10), 0))),
                       nlambda = 100, lambda.min.ratio = 0.01,
                       estimators = NULL) {
  p <- ncol(Q)
  L <- Cholesky(Matrix(forceSymmetric(Q), sparse = TRUE), LDL = FALSE, perm = TRUE)
  z <- matrix(rnorm(n * p), nrow = p, ncol = n)
  x <- solve(L, solve(L, z, system = "P"), system = "Lt")
  x <- Matrix::solve(L, x, system = "Pt")
  data <- as.matrix(t(x))

  S_full  <- cov(data)
  n_train <- floor(split_train * n)
  idx     <- sample.int(n, n_train)
  train   <- data[idx, , drop = FALSE]
  test    <- data[-idx, , drop = FALSE]
  S_train <- cov(train)
  S_test  <- cov(test)
  n_test  <- nrow(test)

  lam_max <- max(abs(S_train - diag(diag(S_train))))
  lam_min <- lambda.min.ratio * lam_max
  lambdas <- exp(seq(log(lam_max), log(lam_min), length.out = nlambda))

  # Default: PC-GLasso, Glasso, CorGL, SPACE
  if (is.null(estimators)) {
    estimators <- list(
      PCGL = estimator_pcglasso,
      GL   = estimator_glasso,
      CGL  = estimator_corglasso,
      SPACE = estimator_space
    )
  }

  res_list <- list()
  for (meth in names(estimators)) {
    est <- estimators[[meth]](S_full, S_train, S_test, n, n_train, n_test, lambdas, alpha_grid=alpha_grid, data=data, train=train, test=test)
    for (sel in names(est)) {
      sel_name <- paste0(meth, "_", gsub(".*_", "", sel)) # eg "GL_bic"
      Qhat <- est[[sel]]$Q
      time <- est[[sel]]$timing
      cmp  <- compare_matrices(Q, Qhat)
      cmp$method <- sel_name
      cmp$timing <- time
      cmp$n <- n
      res_list[[sel_name]] <- cmp
    }
  }
  df <- do.call(rbind, res_list)
  rownames(df) <- NULL
  df
}
run_experiments <- function(Q, ns = c(200,500,1000), sim = 50,
                            mc_cores = parallel::detectCores(), seed=1234, estimators=NULL, ...) {
  grid <- expand.grid(n = ns, rep = seq_len(sim))
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  results <- pbmclapply(
    seq_len(nrow(grid)),
    function(i) {
      row <- grid[i, ]
      df  <- run_single(Q, n = row$n, estimators = estimators, ...)
      cbind(n = row$n, rep = row$rep, df)
    },
    mc.cores    = mc_cores,
    mc.set.seed = TRUE
  )
  do.call(rbind, results)
}
#' Summarize & plot RMSE components, error rates, timing, and true rates across sample sizes
#'
#' @param results Data.frame from run_experiments(), with columns:
#'   n, rep,
#'   rmse, rmse_diag, rmse_offdiag_zero, rmse_offdiag_nonzero,
#'   false_non0_rate, false_0_rate, true_non0_rate, true_0_rate,
#'   method, timing
#' @return A list with:
#'   - table: data.frame of mean metrics by n & method
#'   - plots: list(
#'       rmse_overall, rmse_diag,
#'       rmse_offdiag_zero, rmse_offdiag_nonzero,
#'       fp_rate, fn_rate, tp_rate, tn_rate,
#'       timing_plot,
#'       rmse_grid, rate_grid, rate_grid_true, timing_grid
#'     )
#' @import ggplot2 patchwork
#' @export
summarize_plot_results <- function(results) {
  required <- c(
    "n", "method",
    "rmse", "rmse_diag", "rmse_offdiag_zero", "rmse_offdiag_nonzero",
    "false_non0_rate", "false_0_rate", "true_non0_rate", "true_0_rate", "timing"
  )
  missing_cols <- setdiff(required, names(results))
  if (length(missing_cols)) {
    stop("Missing columns: ", paste(missing_cols, collapse = ", "))
  }

  # Standardize rate names
  results$false_pos_rate <- results$false_non0_rate
  results$false_neg_rate <- results$false_0_rate
  results$true_pos_rate  <- results$true_non0_rate
  results$true_neg_rate  <- results$true_0_rate

  # Ensure 'method' is factor for color consistency
  results$method <- as.factor(results$method)

  # Mean over reps
  summary_df <- aggregate(
    cbind(
      rmse, rmse_diag, rmse_offdiag_zero, rmse_offdiag_nonzero,
      false_pos_rate, false_neg_rate, true_pos_rate, true_neg_rate, timing
    ) ~ n + method,
    data = results, FUN = mean
  )

  table <- summary_df

  col_palette <- "Dark2" # Or Set1, Set2, Paired, etc.

  # Individual RMSE plots
  p_rmse_overall <- ggplot(summary_df, aes(x = n, y = rmse, color = method, group = method)) +
    geom_line(size = 1) + geom_point(size = 2) +
    labs(y = "RMSE Overall", x = "Sample Size (n)") +
    scale_color_brewer(palette = col_palette) +
    theme_minimal()

  p_rmse_diag <- ggplot(summary_df, aes(x = n, y = rmse_diag, color = method, group = method)) +
    geom_line(size = 1) + geom_point(size = 2) +
    labs(y = "RMSE Diagonal", x = "Sample Size (n)") +
    scale_color_brewer(palette = col_palette) +
    theme_minimal()

  p_rmse_off0 <- ggplot(summary_df, aes(x = n, y = rmse_offdiag_zero, color = method, group = method)) +
    geom_line(size = 1) + geom_point(size = 2) +
    labs(y = "RMSE Off-diag (true zero)", x = "Sample Size (n)") +
    scale_color_brewer(palette = col_palette) +
    theme_minimal()

  p_rmse_offnz <- ggplot(summary_df, aes(x = n, y = rmse_offdiag_nonzero, color = method, group = method)) +
    geom_line(size = 1) + geom_point(size = 2) +
    labs(y = "RMSE Off-diag (true non-zero)", x = "Sample Size (n)") +
    scale_color_brewer(palette = col_palette) +
    theme_minimal()

  # Error and true rate plots
  p_fp <- ggplot(summary_df, aes(x = n, y = false_pos_rate, color = method, group = method)) +
    geom_line(size = 1) + geom_point(size = 2) +
    labs(y = "False Positive Rate", x = "Sample Size (n)") +
    scale_color_brewer(palette = col_palette) +
    theme_minimal()

  p_fn <- ggplot(summary_df, aes(x = n, y = false_neg_rate, color = method, group = method)) +
    geom_line(size = 1) + geom_point(size = 2) +
    labs(y = "False Negative Rate", x = "Sample Size (n)") +
    scale_color_brewer(palette = col_palette) +
    theme_minimal()

  p_tp <- ggplot(summary_df, aes(x = n, y = true_pos_rate, color = method, group = method)) +
    geom_line(size = 1) + geom_point(size = 2) +
    labs(y = "True Positive Rate", x = "Sample Size (n)") +
    scale_color_brewer(palette = col_palette) +
    theme_minimal()

  p_tn <- ggplot(summary_df, aes(x = n, y = true_neg_rate, color = method, group = method)) +
    geom_line(size = 1) + geom_point(size = 2) +
    labs(y = "True Negative Rate", x = "Sample Size (n)") +
    scale_color_brewer(palette = col_palette) +
    theme_minimal()

  # Timing plot
  p_timing <- ggplot(summary_df, aes(x = n, y = timing, color = method, group = method)) +
    geom_line(size = 1) + geom_point(size = 2) +
    labs(y = "Timing (seconds)", x = "Sample Size (n)") +
    scale_color_brewer(palette = col_palette) +
    theme_minimal()

  # Combined grids via patchwork
  rmse_grid <- patchwork::wrap_plots(
    p_rmse_overall, p_rmse_diag,
    p_rmse_off0, p_rmse_offnz,
    ncol = 2
  )

  rate_grid <- patchwork::wrap_plots(
    p_fp, p_fn,
    ncol = 2
  )

  rate_grid_true <- patchwork::wrap_plots(
    p_tp, p_tn,
    ncol = 2
  )

  timing_grid <- p_timing

  list(
    table = table,
    plots = list(
      rmse_overall         = p_rmse_overall,
      rmse_diag            = p_rmse_diag,
      rmse_offdiag_zero    = p_rmse_off0,
      rmse_offdiag_nonzero = p_rmse_offnz,
      fp_rate              = p_fp,
      fn_rate              = p_fn,
      tp_rate              = p_tp,
      tn_rate              = p_tn,
      timing               = p_timing,
      rmse_grid            = rmse_grid,
      rate_grid            = rate_grid,
      rate_grid_true       = rate_grid_true,
      timing_grid          = timing_grid
    )
  )
}

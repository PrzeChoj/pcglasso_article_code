# simulation_functions.R
#
# Core simulation engine for the Section 4.3 / Appendix B.2 precision matrix
# estimation study.  Contains three functions:
#
#   run_single()             -- one simulation replicate for a given Q and n
#   run_experiments()        -- orchestrates replications across sample sizes
#   summarize_plot_results() -- aggregates results and produces diagnostic plots
#
# All packages are loaded in run_simulation.R; do not load them here.


# Run one simulation replicate ------------------------------------------------

run_single <- function(Q, n,
                       split_train        = 0.7,
                       alpha_grid         = sort(unique(c(seq(-0.1, 0.1, length.out = 10), 0))),
                       nlambda            = 100,
                       lambda.min.ratio   = 0.01,
                       pcglasso_tolerance = 1e-5,
                       estimators         = NULL) {

  p <- ncol(Q)

  # Generate data from N(0, Q^{-1}) via Cholesky factorisation
  L    <- Cholesky(Matrix(forceSymmetric(Q), sparse = TRUE), LDL = FALSE, perm = TRUE)
  z    <- matrix(rnorm(n * p), nrow = p, ncol = n)
  x    <- solve(L, solve(L, z, system = "P"), system = "Lt")
  x    <- Matrix::solve(L, x, system = "Pt")
  data <- as.matrix(t(x))

  # Train / test split
  S_full  <- cov(data)
  n_train <- floor(split_train * n)
  idx     <- sample.int(n, n_train)
  train   <- data[idx,  , drop = FALSE]
  test    <- data[-idx, , drop = FALSE]
  S_train <- cov(train)
  S_test  <- cov(test)
  n_test  <- nrow(test)

  # Lambda grid (on the covariance scale)
  lam_max <- max(abs(S_train - diag(diag(S_train))))
  lam_min <- lambda.min.ratio * lam_max
  lambdas <- exp(seq(log(lam_max), log(lam_min), length.out = nlambda))

  # Default estimator set: four PC-GLasso variants + three baselines
  if (is.null(estimators)) {
    estimators <- list(
      PCGLcpp_I  = estimator_pcglasso_I_primal,
      PCGLFor_I  = estimator_pcglasso_I_dual,
      GL         = estimator_glasso,
      SPACE      = estimator_space,
      CGL        = estimator_corglasso
    )
  }

  res_list <- list()
  for (meth in names(estimators)) {
    est <- tryCatch(
      estimators[[meth]](
        S_full, S_train, S_test, n, n_train, n_test,
        lambdas,
        alpha_grid         = alpha_grid,
        data               = data,
        train              = train,
        test               = test,
        pcglasso_tolerance = pcglasso_tolerance
      ),
      error = function(e) {
        warning("Estimator '", meth, "' failed: ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(est)) next

    for (sel in names(est)) {
      sel_name <- paste0(meth, "_", gsub(".*_", "", sel))   # e.g. "GL_bic"
      cmp <- tryCatch({
        Qhat      <- est[[sel]]$Q
        cmp       <- compare_matrices(Q, Qhat)
        cmp$method <- sel_name
        cmp$timing <- est[[sel]]$timing
        cmp$n      <- n
        cmp
      }, error = function(e) {
        warning("compare_matrices() failed for '", sel_name, "': ", conditionMessage(e))
        NULL
      })
      if (!is.null(cmp)) res_list[[sel_name]] <- cmp
    }
  }

  df <- if (length(res_list)) do.call(rbind, res_list) else data.frame()
  rownames(df) <- NULL
  df
}


# Run replications across sample sizes ----------------------------------------

run_experiments <- function(Q,
                            ns          = c(200, 500, 1000),
                            sim         = 50,
                            mc_cores    = parallel::detectCores(),
                            seed        = 1234,
                            estimators  = NULL,
                            pcglasso_tolerance = 0.001,
                            ...) {

  dots <- list(...)
  grid <- expand.grid(n = ns, rep = seq_len(sim))

  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)

  results <- parallel::mclapply(
    seq_len(nrow(grid)),
    function(i) {
      row <- grid[i, ]
      df  <- do.call(
        run_single,
        c(
          list(
            Q                  = Q,
            n                  = row$n,
            estimators         = estimators,
            pcglasso_tolerance = pcglasso_tolerance
          ),
          dots
        )
      )
      if (nrow(df) == 0) return(data.frame())
      cbind(n = row$n, rep = row$rep, df)
    },
    mc.cores    = mc_cores,
    mc.set.seed = TRUE
  )

  results <- Filter(function(x) !is.null(x) && nrow(x) > 0, results)
  if (!length(results)) return(data.frame())
  do.call(rbind, results)
}


# Summarise results and produce diagnostic plots ------------------------------

#' @param results  Data frame from run_experiments(), with columns:
#'   n, rep, rmse, rmse_diag, rmse_offdiag_zero, rmse_offdiag_nonzero,
#'   false_non0_rate, false_0_rate, true_non0_rate, true_0_rate, method, timing
#' @return List with:
#'   $table  -- data frame of means by (n, method)
#'   $plots  -- named list of ggplot2 / patchwork objects
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

  results$false_pos_rate <- results$false_non0_rate
  results$false_neg_rate <- results$false_0_rate
  results$true_pos_rate  <- results$true_non0_rate
  results$true_neg_rate  <- results$true_0_rate

  results$method <- as.factor(results$method)

  summary_df <- aggregate(
    cbind(
      rmse, rmse_diag, rmse_offdiag_zero, rmse_offdiag_nonzero,
      false_pos_rate, false_neg_rate, true_pos_rate, true_neg_rate, timing
    ) ~ n + method,
    data = results, FUN = mean
  )

  col_palette <- "Dark2"

  make_line <- function(y_var, y_label) {
    ggplot2::ggplot(
      summary_df,
      ggplot2::aes(x = n, y = .data[[y_var]], color = method, group = method)
    ) +
      ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::geom_point(size = 2) +
      ggplot2::labs(y = y_label, x = "Sample size (n)") +
      ggplot2::scale_color_brewer(palette = col_palette) +
      ggplot2::theme_minimal()
  }

  p_rmse_overall <- make_line("rmse",                "RMSE overall")
  p_rmse_diag    <- make_line("rmse_diag",            "RMSE diagonal")
  p_rmse_off0    <- make_line("rmse_offdiag_zero",    "RMSE off-diag (true zero)")
  p_rmse_offnz   <- make_line("rmse_offdiag_nonzero", "RMSE off-diag (true non-zero)")
  p_fp           <- make_line("false_pos_rate",       "False positive rate")
  p_fn           <- make_line("false_neg_rate",       "False negative rate")
  p_tp           <- make_line("true_pos_rate",        "True positive rate")
  p_tn           <- make_line("true_neg_rate",        "True negative rate")
  p_timing       <- make_line("timing",               "Timing (seconds)")

  rmse_grid      <- patchwork::wrap_plots(p_rmse_overall, p_rmse_diag,
                                          p_rmse_off0, p_rmse_offnz, ncol = 2)
  rate_grid      <- patchwork::wrap_plots(p_fp, p_fn, ncol = 2)
  rate_grid_true <- patchwork::wrap_plots(p_tp, p_tn, ncol = 2)

  list(
    table = summary_df,
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
      rate_grid_true       = rate_grid_true
    )
  )
}

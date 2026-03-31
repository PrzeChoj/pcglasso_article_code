library(parallel)


select_ns_part <- function(ns, part = NULL, n_parts = NULL) {
  if (is.null(part)) {
    return(ns)
  }
  if (is.null(n_parts)) {
    stop("If 'part' is provided, 'n_parts' must also be provided.")
  }
  part <- as.integer(part)
  n_parts <- as.integer(n_parts)
  if (is.na(part) || is.na(n_parts) || part < 1 || n_parts < 1) {
    stop("'part' and 'n_parts' must be positive integers.")
  }
  idx <- split(seq_along(ns), cut(seq_along(ns), breaks = n_parts, labels = FALSE))
  if (part > length(idx) || length(idx[[part]]) == 0) {
    stop("Requested part has no entries in 'ns'.")
  }
  ns[idx[[part]]]
}

solver_from_method <- function(method) {
  if (is.na(method) || is.null(method)) {
    return(NA_character_)
  }
  m <- tolower(method)
  if (grepl("cpp", m)) {
    return("cpp")
  }
  if (grepl("for", m)) {
    return("fortran")
  }
  NA_character_
}

run_single <- function(Q, n, split_train = 0.7,
                       alpha_grid = sort(unique(c(seq(-0.1, 0.1, length.out = 10), 0))),
                       nlambda = 100, lambda.min.ratio = 0.01, pcglasso_tolerance = 0.001,
                       estimators = NULL,
                       collect_errors = FALSE,
                       seed = NA_character_,
                       iter = NA_integer_,
                       rep = NA_integer_) {
  p <- ncol(Q)
  pcglasso_tolerance <- 200*pcglasso_tolerance/n
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
      PCGLcpp_I = estimator_pcglasso_I_cpp,
      PCGLFor_C = estimator_pcglasso_C_fortran,
      PCGLFor_I = estimator_pcglasso_I_fortran,
      PCGLcpp_C = estimator_pcglasso_C_cpp,
      PCGLcart_C = estimator_pcglasso_C_Carter,
      PCGLcart_I = estimator_pcglasso_I_Carter,
      GL        = estimator_glasso,
      CGL       = estimator_corglasso,
      SPACE     = estimator_space
    )
  }

  res_list <- list()
  errors <- list()
  add_error <- function(stage, method, selector, err) {
    if (!collect_errors) {
      return(invisible(NULL))
    }
    errors[[length(errors) + 1]] <<- data.frame(
      iter = iter,
      n = n,
      rep = rep,
      method = method,
      selector = selector,
      solver = solver_from_method(method),
      seed = seed,
      stage = stage,
      error = conditionMessage(err),
      stringsAsFactors = FALSE
    )
    invisible(NULL)
  }
  for (meth in names(estimators)) {
    est <- tryCatch(
      estimators[[meth]](
        S_full, S_train, S_test, n, n_train, n_test,
        lambdas, alpha_grid = alpha_grid, data = data, train = train, test = test,
        pcglasso_tolerance = pcglasso_tolerance
      ),
      error = function(e) {
        add_error("estimator", meth, NA_character_, e)
        NULL
      }
    )
    if (is.null(est)) {
      next
    }
    for (sel in names(est)) {
      sel_name <- paste0(meth, "_", gsub(".*_", "", sel)) # eg "GL_bic"
      cmp <- tryCatch({
        Qhat <- est[[sel]]$Q
        time <- est[[sel]]$timing
        cmp  <- compare_matrices(Q, Qhat)
        cmp$method <- sel_name
        cmp$timing <- time
        cmp$n <- n
        cmp
      }, error = function(e) {
        add_error("compare", meth, sel, e)
        NULL
      })
      if (!is.null(cmp)) {
        res_list[[sel_name]] <- cmp
      }
    }
  }
  df <- if (length(res_list)) do.call(rbind, res_list) else data.frame()
  rownames(df) <- NULL
  if (collect_errors) {
    err_df <- if (length(errors)) {
      do.call(rbind, errors)
    } else {
      data.frame(
        iter = integer(),
        n = integer(),
        rep = integer(),
        method = character(),
        selector = character(),
        solver = character(),
        seed = character(),
        stage = character(),
        error = character(),
        stringsAsFactors = FALSE
      )
    }
    list(results = df, errors = err_df)
  } else {
    df
  }
}
run_experiments <- function(Q, ns = c(200,500,1000), sim = 50,
                            mc_cores = parallel::detectCores(), seed=1234, estimators=NULL, pcglasso_tolerance = 0.001,
                            rep_ids = NULL, return_errors = FALSE, ...) {
  if (is.null(rep_ids)) {
    rep_ids <- seq_len(sim)
  }
  dots <- list(...)
  dots$return_errors <- NULL
  grid <- expand.grid(n = ns, rep = rep_ids)
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  results <- parallel::mclapply(
    seq_len(nrow(grid)),
    function(i) {
      row <- grid[i, ]
      seed_val <- paste(.Random.seed, collapse = ",")
      if (!return_errors) {
        df <- do.call(
          run_single,
          c(
            list(Q = Q, n = row$n, estimators = estimators, pcglasso_tolerance = pcglasso_tolerance),
            dots
          )
        )
        if (nrow(df) == 0) {
          return(data.frame())
        }
        return(cbind(n = row$n, rep = row$rep, df))
      }

      res <- tryCatch(
        do.call(
          run_single,
          c(
            list(
              Q = Q, n = row$n, estimators = estimators, pcglasso_tolerance = pcglasso_tolerance,
              collect_errors = TRUE, seed = seed_val, iter = i, rep = row$rep
            ),
            dots
          )
        ),
        error = function(e) {
          err <- data.frame(
            iter = i,
            n = row$n,
            rep = row$rep,
            method = NA_character_,
            selector = NA_character_,
            solver = NA_character_,
            seed = seed_val,
            stage = "run_single",
            error = conditionMessage(e),
            stringsAsFactors = FALSE
          )
          list(results = data.frame(), errors = err)
        }
      )

      res$results <- if (nrow(res$results) == 0) {
        res$results
      } else {
        cbind(n = row$n, rep = row$rep, res$results)
      }
      res
    },
    mc.cores    = mc_cores,
    mc.set.seed = TRUE
  )
  if (!return_errors) {
    results <- Filter(function(x) !is.null(x) && nrow(x) > 0, results)
    if (!length(results)) {
      return(data.frame())
    }
    return(do.call(rbind, results))
  }

  res_list <- lapply(results, function(x) x$results)
  err_list <- lapply(results, function(x) x$errors)
  res_list <- Filter(function(x) !is.null(x) && nrow(x) > 0, res_list)
  err_list <- Filter(function(x) !is.null(x), err_list)
  res_df <- if (length(res_list)) do.call(rbind, res_list) else data.frame()
  err_df <- if (length(err_list)) do.call(rbind, err_list) else data.frame()
  list(results = res_df, errors = err_df)
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

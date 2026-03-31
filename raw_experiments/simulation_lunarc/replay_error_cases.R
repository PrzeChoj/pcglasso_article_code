#!/usr/bin/env Rscript
#
# replay_error_cases.R
#
# Re-run error cases sequentially (no parallel) using errors_*.rds file(s).
#
# Usage:
#   Rscript replay_error_cases.R
#   Rscript replay_error_cases.R errors_file.rds [row_ids] [out_file]
#
# Examples:
#   Rscript replay_error_cases.R
#   Rscript replay_error_cases.R errors_hub_big_n5000_part1of7.rds
#   Rscript replay_error_cases.R errors_hub_big_n5000_part1of7.rds 3
#   Rscript replay_error_cases.R errors_hub_big_n5000_part1of7.rds 1,5,9 replay_out.rds
#

suppressPackageStartupMessages({
  library(space)
  library(PCGLASSOcpp)
  library(PCGLASSO)
  library(pcglassoFast)
  source("estimation_methods.R")
  source("simulation_functions.R")
  library(glasso)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 3) {
  stop("Too many arguments. Expected: errors_file [row_ids] [out_file].")
}

auto_errors <- function() {
  files <- list.files(pattern = "^errors_.*\\.rds$")
  if (length(files) == 0) {
    stop("No errors_*.rds files found in the current directory.")
  }
  files
}

parse_rows <- function(arg, n_rows) {
  if (tolower(arg) == "all") {
    return(seq_len(n_rows))
  }
  parts <- unlist(strsplit(arg, ",", fixed = TRUE))
  parts <- parts[parts != ""]
  idx <- as.integer(parts)
  if (length(idx) == 0 || any(is.na(idx))) {
    stop("row_ids must be 'all' or a comma-separated list of integers.")
  }
  idx
}

infer_graph_type <- function(path) {
  base <- tolower(basename(path))
  if (grepl("non_hub", base)) return("non_hub")
  if (grepl("hub", base)) return("hub")
  NA_character_
}

data(Q_simulated_glasso)
data(Q_simulated_pcglasso)
Q_map <- list(
  hub = Q_simulated_pcglasso,
  non_hub = Q_simulated_glasso
)

get_estimators <- function() {
  list(
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

set_seed_from_string <- function(seed_str) {
  if (is.na(seed_str) || seed_str == "") {
    stop("Missing seed in error row.")
  }
  parts <- unlist(strsplit(seed_str, ",", fixed = TRUE))
  parts <- trimws(parts)
  seed_vec <- as.integer(parts)
  if (any(is.na(seed_vec))) {
    stop("Seed string could not be parsed into integers.")
  }
  RNGkind("L'Ecuyer-CMRG")
  .Random.seed <<- seed_vec
}

build_data <- function(Q, n, split_train = 0.7,
                       nlambda = 50, lambda.min.ratio = 0.01) {
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

  list(
    S_full = S_full,
    S_train = S_train,
    S_test = S_test,
    n_train = n_train,
    n_test = n_test,
    lambdas = lambdas,
    data = data,
    train = train,
    test = test
  )
}

run_case <- function(err_row, idx, errors_path) {
  graph_type <- if ("graph_type" %in% names(err_row) && !is.na(err_row$graph_type)) {
    err_row$graph_type
  } else {
    infer_graph_type(errors_path)
  }
  if (is.na(graph_type) || !graph_type %in% names(Q_map)) {
    stop("graph_type not found; add graph_type column or use hub/non_hub in filename.")
  }
  Q <- Q_map[[graph_type]]
  Q <- (Q + t(Q)) / 2

  n <- as.integer(err_row$n)
  rep_id <- as.integer(err_row$rep)
  method <- as.character(err_row$method)
  selector <- as.character(err_row$selector)
  stage <- as.character(err_row$stage)

  set_seed_from_string(err_row$seed)

  d <- build_data(Q, n = n, split_train = 0.7, nlambda = 50, lambda.min.ratio = 0.01)
  estimators <- get_estimators()
  if (!is.na(method) && method %in% names(estimators)) {
    estimators <- estimators[method]
  }

  if (length(estimators) == 0) {
    stop("Estimator not found for method: ", method)
  }

  est <- estimators[[1]](
    d$S_full, d$S_train, d$S_test, n, d$n_train, d$n_test,
    d$lambdas, alpha_grid = 0, data = d$data, train = d$train, test = d$test,
    pcglasso_tolerance = 0.001
  )

  if (stage == "estimator") {
    return(list(
      row = idx,
      n = n,
      rep = rep_id,
      method = method,
      stage = stage,
      selector = selector,
      ok = TRUE,
      selections = names(est)
    ))
  }

  sels <- names(est)
  if (!is.na(selector) && selector %in% names(est)) {
    sels <- selector
  }

  results <- list()
  for (sel in sels) {
    Qhat <- est[[sel]]$Q
    cmp <- compare_matrices(Q, Qhat)
    cmp$method <- paste0(method, "_", gsub(".*_", "", sel))
    cmp$n <- n
    results[[sel]] <- cmp
  }

  list(
    row = idx,
    n = n,
    rep = rep_id,
    method = method,
    stage = stage,
    selector = selector,
    ok = TRUE,
    results = results
  )
}

errors_paths <- if (length(args) >= 1) {
  args[[1]]
} else {
  auto_errors()
}

run_for_file <- function(errors_path, row_arg = "all", out_path = NULL) {
  errors <- readRDS(errors_path)
  if (!is.data.frame(errors)) {
    stop("errors_file must contain a data.frame: ", errors_path)
  }
  if (nrow(errors) == 0) {
    message("No errors found in ", errors_path)
    return(invisible(NULL))
  }

  rows <- parse_rows(row_arg, nrow(errors))
  replays <- list()
  for (i in rows) {
    err_row <- errors[i, , drop = FALSE]
    message("Replaying row ", i, " (n=", err_row$n, ", rep=", err_row$rep, ", method=", err_row$method, ")")
    res <- tryCatch(
      run_case(err_row, i, errors_path),
      error = function(e) {
        list(
          row = i,
          n = err_row$n,
          rep = err_row$rep,
          method = err_row$method,
          stage = err_row$stage,
          selector = err_row$selector,
          ok = FALSE,
          error = conditionMessage(e)
        )
      }
    )
    replays[[length(replays) + 1]] <- res
  }

  if (is.null(out_path) || out_path == "") {
    base <- sub("\\.rds$", "", basename(errors_path))
    out_path <- file.path(dirname(errors_path), paste0("replay_", base, ".rds"))
  }
  saveRDS(replays, file = out_path)
  message("Saved replay results to ", out_path)
}

if (length(errors_paths) == 1 && length(args) >= 2) {
  row_arg <- args[[2]]
  out_path <- if (length(args) >= 3) args[[3]] else NULL
  run_for_file(errors_paths, row_arg = row_arg, out_path = out_path)
} else {
  for (errors_path in errors_paths) {
    message("Processing ", errors_path)
    run_for_file(errors_path)
  }
}

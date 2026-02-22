# 10 minutes on 32 cores of sr

source("./experiments/Appendix_A/utils.R")
source("./experiments/Appendix_A/experiment_2/utils.R")
source("./experiments/Appendix_A/experiment_2/0_parameters.R")

load("./experiments/Appendix_A/experiment_2/res_data/instances.RData")
baseline_best_value <- read.csv("./experiments/Appendix_A/experiment_2/res_data/baseline_best_value.csv")

library(parallel)
library(pbmcapply)

n_cores <- max(1L, detectCores(logical = FALSE) - 1L)

grid <- expand.grid(
  structure = graph_structure_vec,
  p = p_vec,
  solver = solver_vec,
  start = starting_point_vec,
  stringsAsFactors = FALSE
)

passes_at_tol <- function(S, solver, start, tol, baseline, R_strict, acceptable_error) {
  for (r in seq_len(R_strict)) {
    f_end <- value_after_optimization(S, solver, start, tol, lambda, alpha)
    gap   <- f_end - baseline
    if (!is.finite(gap) || gap > acceptable_error) return(FALSE)
  }
  TRUE
}

find_tol_found <- function(S, solver, start, baseline, tol_grid, R_strict, acceptable_error) {
  lo <- 1L
  hi <- length(tol_grid)
  ans <- NA_real_

  while (lo <= hi) {
    mid <- (lo + hi) %/% 2L
    tol_mid <- tol_grid[mid]

    ok <- passes_at_tol(S, solver, start, tol_mid, baseline, R_strict, acceptable_error)

    if (ok) {
      ans <- tol_mid
      hi <- mid - 1L
    } else {
      lo <- mid + 1L
    }
  }

  ans
}

one_job <- function(i) {
  structure <- grid$structure[i]
  p <- grid$p[i]
  solver <- grid$solver[i]
  start <- grid$start[i]

  S <- instances[[structure]][[as.character(p)]]
  if (is.null(S)) stop("Missing S for structure=", structure, " p=", p)

  baseline <- get_baseline(p, structure)

  tol_found <- find_tol_found(
    S = S,
    solver = solver,
    start = start,
    baseline = baseline,
    tol_grid = tol_grid,
    R_strict = R_strict,
    acceptable_error = acceptable_error
  )

  data.frame(
    p = p,
    structure = structure,
    solver = solver,
    start = start,
    tol_found = tol_found,
    stringsAsFactors = FALSE
  )
}

res <- pbmclapply(seq_len(nrow(grid)), one_job, mc.cores = n_cores)
calibration_tol_found <- do.call(rbind, res)

out_path <- "./experiments/Appendix_A/experiment_2/res_data/calibration_tol_found.csv"
write.csv(calibration_tol_found, out_path, row.names = FALSE)

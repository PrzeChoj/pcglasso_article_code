# ?? minutes on 32 cores of sr

source("./experiments/Appendix_A/utils.R")
source("./experiments/Appendix_A/experiment_2/utils.R")
source("./experiments/Appendix_A/experiment_2/0_parameters.R")

load("./experiments/Appendix_A/experiment_2/res_data/instances.RData")
baseline_best_value <- read.csv("./experiments/Appendix_A/experiment_2/res_data/baseline_best_value.csv")
calibration_tol_found <- read.csv("./experiments/Appendix_A/experiment_2/res_data/calibration_tol_found.csv")

stopifnot(all(!is.na(calibration_tol_found$tol_found)))

library(parallel)
library(pbmcapply)
n_cores <- max(1L, detectCores(logical = FALSE) - 1L)

get_tolerance <- function(p, structure, solver, start) {
  row <- calibration_tol_found[
    calibration_tol_found$p == p &
      calibration_tol_found$structure == structure &
      calibration_tol_found$solver == solver &
      calibration_tol_found$start == start,
    ,
    drop = FALSE
  ]
  if (nrow(row) != 1L) stop("Tolerance not found/unique")
  row$tol_found[1]
}

# build grid
grid_cfg <- expand.grid(
  structure = graph_structure_vec,
  p = p_vec,
  solver = solver_vec,
  start = starting_point_vec,
  stringsAsFactors = FALSE
)

grid_cfg$tol_found <- mapply(
  get_tolerance,
  p = grid_cfg$p,
  structure = grid_cfg$structure,
  solver = grid_cfg$solver,
  start = grid_cfg$start
)

grid <- grid_cfg[rep(seq_len(nrow(grid_cfg)), each = M), , drop = FALSE]
grid$run_id <- rep(seq_len(M), times = nrow(grid_cfg))
rownames(grid) <- NULL

time_value_for_tolerance <- function(S, solver, start, tol) {
  t0 <- proc.time()[["elapsed"]]
  f_end <- value_for_tolerance(S, solver, start, tol)
  t1 <- proc.time()[["elapsed"]]
  list(time = t1 - t0, f_end = f_end)
}

one_run <- function(i) {
  structure <- grid$structure[i]
  p <- grid$p[i]
  solver <- grid$solver[i]
  start <- grid$start[i]
  tol <- grid$tol_found[i]
  run_id <- grid$run_id[i]

  S <- instances[[structure]][[as.character(p)]]
  if (is.null(S)) stop("Missing S for structure=", structure, " p=", p)

  baseline <- get_baseline(p, structure)

  out <- time_value_for_tolerance(S, solver, start, tol)
  gap <- out$f_end - baseline
  success <- is.finite(gap) && (gap <= acceptable_error)

  data.frame(
    p = p,
    structure = structure,
    solver = solver,
    start = start,
    tol_found = tol,
    run_id = run_id,
    time = out$time,
    f_end = out$f_end,
    gap = gap,
    success = success,
    stringsAsFactors = FALSE
  )
}

timing_raw_list <- pbmclapply(seq_len(nrow(grid)), one_run, mc.cores = n_cores)
timing_raw <- do.call(rbind, timing_raw_list)
rownames(timing_raw) <- NULL

# summary
cfg_key <- paste(timing_raw$p, timing_raw$structure, timing_raw$solver, timing_raw$start, sep = "|")
idx_split <- split(seq_len(nrow(timing_raw)), cfg_key)

timing_summary <- do.call(rbind, lapply(idx_split, function(idx) {
  df <- timing_raw[idx, , drop = FALSE]
  data.frame(
    p = df$p[1],
    structure = df$structure[1],
    solver = df$solver[1],
    start = df$start[1],
    tol_found = df$tol_found[1],
    time_trimmed_mean = mean(df$time, trim = 0.1),
    success_rate = mean(df$success),
    stringsAsFactors = FALSE
  )
}))

write.csv(
  timing_raw,
  "./experiments/Appendix_A/experiment_2/res_data/timing_raw.csv",
  row.names = FALSE
)

write.csv(
  timing_summary,
  "./experiments/Appendix_A/experiment_2/res_data/timing_summary.csv",
  row.names = FALSE
)

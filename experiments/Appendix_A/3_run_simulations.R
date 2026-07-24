# 5 hours of 100 cores of AMD EPYC Rome 7742

library(parallel)
library(pbmcapply)
library(dplyr)
library(readr)

source("./experiments/Appendix_A/0_parameters.R")
source("./experiments/Appendix_A/2_functions_simulations.R")

data_dir <- "./experiments/Appendix_A/res_data"
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

#n_cores <- max(1, detectCores(logical = FALSE) - 1)
n_cores <- 48

set.seed(1234)
Sys.setenv(OMP_NUM_THREADS = 1)


param_grid <- expand.grid(
  p = p_vec,
  lambda = lambda_vec,
  alpha = alpha_vec,
  K_structure = graph_structure_vec,
  solver = solver_vec,
  starting_point = starting_point_vec,
  tolerance = tolerance_list,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

param_grid_M <- param_grid[rep(seq_len(nrow(param_grid)), each = M), , drop = FALSE]
param_grid_M$m <- rep(seq_len(M), times = nrow(param_grid))
rownames(param_grid_M) <- NULL

param_grid_M <- param_grid_M %>%
  arrange(desc(p), desc(-log10(tolerance)))

nrow(param_grid_M)

start_time <- Sys.time()
raw_list <- pbmclapply(
  seq_len(nrow(param_grid_M)),
  function(i) run_one(param_grid_M[i, ]),
  mc.cores = n_cores, mc.set.seed = TRUE
)
raw <- bind_rows(raw_list)
end_time <- Sys.time()

print(end_time - start_time)
write_csv(raw, file.path(data_dir, sprintf("raw_M%d.csv", M)))
if (all(raw$status == "ok")) {
  message("Simulation Successful")
} else {
  stop("Simulation FAILED")
}
if (anyNA(raw$time) || anyNA(raw$f_end)) {
  stop("Simulation FAILED: NA in time or f_end")
}
if (nrow(raw) != nrow(param_grid_M) ||
  !all(c("status", "time", "f_end") %in% names(raw)) ||
  anyNA(raw$time) || anyNA(raw$f_end)) {
  stop("Sth strange went wrong")
}

# summary
summary_data <- raw %>%
  group_by(p, lambda, alpha, K_structure, solver, starting_point, tolerance) %>%
  summarize(
    time_median = median(time),
    f_end_best = min(f_end),
    f_end_worst = max(f_end),
    f_end = mean(f_end),
    .groups = "drop"
  )

if (all(abs(summary_data$f_end_worst - summary_data$f_end_best) < 1e-10)) {
  message("All simuations arrived at the same point")
} else {
  message("There were simuations that arrived at different points!")
}

summary_data <- dplyr::select(summary_data, !starts_with("f_end_"))

write_csv(summary_data, file.path(data_dir, sprintf("summary_M%d.csv", M)))

# 1 minute of 7 cores of Apple's M2

library(ggplot2)
library(dplyr)
library(readr)
library(future.apply)

source("./experiments/Appendix_A/utils.R")
source("./experiments/Appendix_A/experiment_1/0_parameters.R")
source("./experiments/Appendix_A/experiment_1/2_functions.R")

load("./experiments/Appendix_A/experiment_1/res_data/instances.RData")

data_dir <- "./experiments/Appendix_A/experiment_1/res_data"
plot_dir <- "./experiments/Appendix_A/experiment_1/plots"
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

n_cores <- max(1, min(7, parallel::detectCores(logical = FALSE) - 1))

# read summary
summary_path <- file.path(data_dir, sprintf("experiment_1_summary_M%d.csv", M))
df_all <- read_csv(summary_path, show_col_types = FALSE)

# expected columns
expected_columns <- c("p", "lambda", "alpha", "K_structure", "solver", "starting_point", "tolerance", "time_trimmed_mean", "f_end")
stopifnot(all(expected_columns %in% names(df_all)))

# labels used by plot
df_all <- df_all %>%
  mutate(
    time = time_trimmed_mean,
    value = f_end,
    alg = solver,
    init = starting_point
  )

group_keys <- df_all %>%
  distinct(p, lambda, alpha, K_structure)


start_time <- Sys.time()
plan(multisession, workers = n_cores)
raw_list <- future_lapply(
  seq_len(nrow(group_keys)),
  save_plot_for_i
)
end_time <- Sys.time()

message("Done.")
print(end_time - start_time)

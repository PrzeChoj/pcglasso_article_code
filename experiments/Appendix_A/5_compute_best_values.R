# 1 minute of 7 cores of Apple M2

library(dplyr)
library(readr)
library(future.apply)

source("./experiments/Appendix_A/utils.R")
source("./experiments/Appendix_A/0_parameters.R")
source("./experiments/Appendix_A/4_functions_plot.R")

data_dir <- "./experiments/Appendix_A/res_data"
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

raw_path <- file.path(data_dir, sprintf("raw_M%d.csv", M))
df_raw <- read_csv(raw_path, show_col_types = FALSE)

n_cores <- max(1, parallel::detectCores(logical = FALSE) - 1)

group_keys <- tidyr::expand_grid(
  p = p_vec,
  lambda = lambda_vec,
  alpha = alpha_vec,
  K_structure = graph_structure_vec
) %>%
  arrange(desc(p), K_structure, lambda, alpha)

start_time <- Sys.time()
plan(multisession, workers = n_cores)

best_values <- future_sapply(seq_len(nrow(group_keys)), compute_best_value_for_i)

stopifnot(all(is.finite(best_values)))

group_keys <- group_keys %>%
  mutate(best_value = as.numeric(best_values))

out_path <- file.path(data_dir, sprintf("group_keys_with_best_value_M%d.csv", M))
write_csv(group_keys, out_path)
message("Saved: ", out_path)

end_time <- Sys.time()
message("Done.")
print(end_time - start_time)

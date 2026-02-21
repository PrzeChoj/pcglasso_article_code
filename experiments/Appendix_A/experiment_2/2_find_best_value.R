# Step 1 — Baseline map (p, structure) -> best_value

source("./experiments/Appendix_A/utils.R")
source("./experiments/Appendix_A/experiment_2/0_parameters.R")

load("./experiments/Appendix_A/experiment_2/res_data/instances.RData")

library(parallel)
library(pbmcapply)

# build grid
grid <- expand.grid(
  structure = graph_structure_vec,
  p = p_vec,
  stringsAsFactors = FALSE
)

n_cores <- max(1L, min(detectCores(logical = FALSE) - 1L, nrow(grid)))

res_list <- pbmclapply(seq_len(nrow(grid)), function(i) {
  graph_structure <- grid$structure[i]
  p <- grid$p[i]

  S <- instances[[graph_structure]][[as.character(p)]]
  if (is.null(S)) stop("Missing S for structure=", graph_structure, " p=", p)

  best_method <- get_best_method(
    p = p, graph_structure = graph_structure,
    lambda = lambda, alpha = alpha
  )

  best_value <- compute_function_value(
    p = p,
    lambda = lambda,
    alpha = alpha,
    method = best_method,
    tolerance = tol_strict,
    S = S
  )

  data.frame(
    p = p,
    structure = graph_structure,
    best_value = best_value,
    stringsAsFactors = FALSE
  )
}, mc.cores = n_cores)

baseline_best_value <- do.call(rbind, res_list)

out_path <- "./experiments/Appendix_A/experiment_2/res_data/baseline_best_value.csv"
write.csv(baseline_best_value, out_path, row.names = FALSE)

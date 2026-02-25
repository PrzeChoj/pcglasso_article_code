# 2 minutes on 1 core of dgx

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

  solver <- substr(best_method, 1, nchar(best_method) - 2)
  start <- substr(best_method, nchar(best_method), nchar(best_method))
  best_value <- value_after_optimization(S, solver, start, tol_strict, lambda, alpha)

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

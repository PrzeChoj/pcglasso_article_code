source("./experiments/Appendix_A/utils.R")
source("./experiments/Appendix_A/0_parameters.R")

set.seed(1234)

instances <- setNames(vector("list", length(graph_structure_vec)), graph_structure_vec)
for (graph_structure in graph_structure_vec) {
  list_for_p <- setNames(vector("list", length(p_vec)), as.character(p_vec))
  for (p_index in seq_along(p_vec)) {
    p <- p_vec[p_index]
    n <- p * 2
    K_star <- build_K_star(p, graph_structure)
    list_for_p[[as.character(p)]] <- S_from_K_star(K_star, n)
  }

  instances[[graph_structure]] <- list_for_p
}

path <- "./experiments/Appendix_A/res_data"
dir.create(path, showWarnings = FALSE, recursive = TRUE)
save(instances, file = file.path(path, "instances.RData"))

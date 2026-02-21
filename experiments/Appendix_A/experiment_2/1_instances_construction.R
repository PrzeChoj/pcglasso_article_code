source("./experiments/Appendix_A/utils.R")
source("./experiments/Appendix_A/experiment_2/0_parameters.R")

instances <- setNames(vector("list", length(graph_structure_vec)), graph_structure_vec)

for (graph_structure in graph_structure_vec) {
  S_list <- setNames(vector("list", length(p_vec)), as.character(p_vec))

  for (i in seq_along(p_vec)) {
    p <- p_vec[[i]]
    n <- 2 * p

    K_star <- build_K_star(p, cor_modifier_map[[graph_structure]], graph_structure)
    S_list[[i]] <- S_from_K_star(K_star, n)
  }

  instances[[graph_structure]] <- S_list
}

save(instances, file = "./experiments/Appendix_A/experiment_2/res_data/instances.RData")

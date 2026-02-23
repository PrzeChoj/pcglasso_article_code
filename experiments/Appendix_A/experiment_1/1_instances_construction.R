source("./experiments/Appendix_A/experiment_1/0_parameters.R")

set.seed(1234)

instances <- setNames(vector("list", length(graph_structure_vec)), graph_structure_vec)
for (graph_structure in graph_structure_vec) {
  cor_modifier_vec <- cor_modifier_map[[graph_structure]]

  list_for_p <- setNames(vector("list", length(p_vec)), as.character(p_vec))
  for (p_index in seq_along(p_vec)) {
    p <- p_vec[p_index]
    n <- p * 2

    S_list <- vector("list", 2)
    for (i in c(1, 2)) {
      cor <- cor_modifier_vec[i]
      K_star <- build_K_star(p, cor, graph_structure)

      S_list[[as.character(cor)]] <- S_from_K_star(K_star, n)
    }

    list_for_p[[as.character(p)]] <- S_list
  }

  instances[[graph_structure]] <- list_for_p
}

save(instances, file = "./experiments/Appendix_A/experiment_1/res_data/instances.RData")

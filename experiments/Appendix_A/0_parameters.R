M <- 2
p_vec <- c(50, 100, 150, 200)
lambda_vec <- c(0.1, 0.2)
alpha_vec <- c(0, 0.5)
tolerance_list <- exp(seq(log(0.01), log(0.00000001), length.out = 12))
graph_structure_vec <- c("hub_1", "hub_09", "AR2", "random")
solver_vec <- c("pcglassoFast_Primal", "pcglassoFast_Dual", "pcglassoFast_PrimalDual", "pcglasso")
starting_point_vec <- c("C", "I")
pcglasso_tolerance_multiplier <- 1000

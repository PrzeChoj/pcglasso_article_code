M <- 2
p_vec <- c(10, 50, 100, 150, 200)
cor_modifier_map <- list(
  hub = c(1, 0.9),
  line = c(0.8, 0.9)
)
lambda_vec <- c(0.1, 0.2)
alpha_vec <- c(0, 0.5)
tolerance_list <- exp(seq(log(0.01), log(0.00000001), length.out = 12))
graph_structure_vec <- c("hub", "line")
solver_vec <- c("pcglasso_cpp", "pcglasso_fortran", "pcglasso")
starting_point_vec <- c("C", "I")
pcglasso_tolerance_multiplier <- 100

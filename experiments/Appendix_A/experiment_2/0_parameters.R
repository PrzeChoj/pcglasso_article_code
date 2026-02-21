p_vec <- c(10, 50, 100, 150, 200)
cor_modifier_map <- c(hub = 1, line = 0.9)
lambda <- 0.1
alpha <- 0
c_parameter <- 1 - alpha
graph_structure_vec <- c("hub", "line")
solver_vec <- c("pcglasso_cpp", "pcglasso_fortran", "pcglasso")
starting_point_vec <- c("C", "I")
tol_grid <- 10^((-2):(-10))
tol_strict <- 1e-13
acceptable_error <- 1e-4
R_strict <- 2
M <- 2

# devtools::install_github("JackStorrorCarter/PCGLASSO")
# devtools::install_github("PrzeChoj/pcglassoFast")
# devtools::install_github("PrzeChoj/PCGLASSOcpp")

library(PCGLASSO)
library(pcglassoFast)
library(PCGLASSOcpp)
library(MASS)
library(ggplot2)

set.seed(1234)

pcglasso_goal_function <- function(S, lambda, alpha, delta_matrix, theta_diag) {
  p <- nrow(delta_matrix)
  res_value <- -determinant(delta_matrix)$modulus - 2 * (1 - alpha) * sum(log(theta_diag)) + sum(diag(S %*% diag(theta_diag) %*% delta_matrix %*% diag(theta_diag))) + lambda * sum(abs(delta_matrix - diag(p)))
  attr(res_value, "logarithm") <- NULL

  res_value
}

M <- 5
p <- 50
n <- p*2
rho <- 0.2
lambda <- rho
c_parameter <- 1
alpha <- 1 - c_parameter
tolerance_list <- exp(seq(log(0.01), log(0.00001), length.out = 12))
pcglasso_tolerance_modifier <- 100

cor_modifier <- 1
S_star <- diag(1, p); S_star[1,2:p] <- S_star[2:p,1] <- -cor_modifier/sqrt(p); S_star[1,1] <- 1
# S_star <- diag(1, p)
# for (i in 2:p) {
#   S_star[i-1,i] <- S_star[i,i-1] <- 0.4
# }

Z <- mvrnorm(n = n, mu = rep(0, p), Sigma = solve(S_star))
S <- t(Z) %*% Z / n
S <- cov2cor(S)


counter <- 0
pb <- txtProgressBar(min = counter, max = length(tolerance_list) * M, style = 3)
time_pcglasso_C           <- numeric(length(tolerance_list))
res_val_pcglasso_C        <- numeric(length(tolerance_list))
time_pcglasso_I           <- numeric(length(tolerance_list))
res_val_pcglasso_I        <- numeric(length(tolerance_list))
time_pcglassoFast_I       <- numeric(length(tolerance_list))
res_val_pcglassoFast_I    <- numeric(length(tolerance_list))
time_pcglassoFast_C       <- numeric(length(tolerance_list))
res_val_pcglassoFast_C    <- numeric(length(tolerance_list))
time_pcglasso_cpp_I       <- numeric(length(tolerance_list))
res_val_pcglasso_cpp_I    <- numeric(length(tolerance_list))
time_pcglasso_cpp_C       <- numeric(length(tolerance_list))
res_val_pcglasso_cpp_C    <- numeric(length(tolerance_list))
time_pcglasso_C_mat        <- matrix(NA, nrow = M, ncol = length(tolerance_list))
res_val_pcglasso_C_mat     <- matrix(NA, nrow = M, ncol = length(tolerance_list))
time_pcglasso_I_mat        <- matrix(NA, nrow = M, ncol = length(tolerance_list))
res_val_pcglasso_I_mat     <- matrix(NA, nrow = M, ncol = length(tolerance_list))
time_pcglassoFast_I_mat    <- matrix(NA, nrow = M, ncol = length(tolerance_list))
res_val_pcglassoFast_I_mat <- matrix(NA, nrow = M, ncol = length(tolerance_list))
time_pcglassoFast_C_mat    <- matrix(NA, nrow = M, ncol = length(tolerance_list))
res_val_pcglassoFast_C_mat <- matrix(NA, nrow = M, ncol = length(tolerance_list))
time_pcglasso_cpp_I_mat    <- matrix(NA, nrow = M, ncol = length(tolerance_list))
res_val_pcglasso_cpp_I_mat <- matrix(NA, nrow = M, ncol = length(tolerance_list))
time_pcglasso_cpp_C_mat    <- matrix(NA, nrow = M, ncol = length(tolerance_list))
res_val_pcglasso_cpp_C_mat <- matrix(NA, nrow = M, ncol = length(tolerance_list))
for (m in 1:M) {
  for (i in 1:length(tolerance_list)) {
    tolerance <- tolerance_list[i]

    counter <- counter + 1
    setTxtProgressBar(pb, counter)

    # pcglasso_C
    start_pcglasso_C <- Sys.time()
    res_pcglasso_C <- pcglasso(
      S, lambda, c_parameter,
      threshold = tolerance * pcglasso_tolerance_modifier
    )
    end_pcglasso_C <- Sys.time()
    time_pcglasso_C_mat[m, i] <- as.numeric(difftime(end_pcglasso_C, start_pcglasso_C, units = "secs"))

    Sinv <- res_pcglasso_C
    res_val_pcglasso_C_mat[m, i] <- pcglasso_goal_function(
      S, lambda, alpha, cov2cor(Sinv), sqrt(diag(Sinv))
    )

    # pcglasso_I
    start_pcglasso_I <- Sys.time()
    res_pcglasso_I <- pcglasso(
      S, lambda, c_parameter,
      threshold = tolerance * pcglasso_tolerance_modifier
    )
    end_pcglasso_I <- Sys.time()
    time_pcglasso_I_mat[m, i] <- as.numeric(difftime(end_pcglasso_I, start_pcglasso_I, units = "secs"))

    Sinv <- res_pcglasso_I
    res_val_pcglasso_I_mat[m, i] <- pcglasso_goal_function(
      S, lambda, alpha, cov2cor(Sinv), sqrt(diag(Sinv))
    )

    # pcglassoFast_I
    start_pcglassoFast_I <- Sys.time()
    res_pcglassoFast_I <- pcglassoFast(
      S, lambda = lambda, alpha = alpha,
      tolerance = tolerance
    )
    end_pcglassoFast_I <- Sys.time()
    time_pcglassoFast_I_mat[m, i] <- as.numeric(difftime(end_pcglassoFast_I, start_pcglassoFast_I, units = "secs"))

    Sinv <- res_pcglassoFast_I$Sinv
    res_val_pcglassoFast_I_mat[m, i] <- pcglasso_goal_function(
      S, lambda, alpha, cov2cor(Sinv), sqrt(diag(Sinv))
    )

    # pcglassoFast_C
    start_pcglassoFast_C <- Sys.time()
    res_pcglassoFast_C <- pcglassoFast(
      S, lambda = lambda, alpha = alpha,
      R = cov2cor(solve(S)),
      tolerance = tolerance
    )
    end_pcglassoFast_C <- Sys.time()
    time_pcglassoFast_C_mat[m, i] <- as.numeric(difftime(end_pcglassoFast_C, start_pcglassoFast_C, units = "secs"))

    Sinv <- res_pcglassoFast_C$Sinv
    res_val_pcglassoFast_C_mat[m, i] <- pcglasso_goal_function(
      S, lambda, alpha, cov2cor(Sinv), sqrt(diag(Sinv))
    )

    # pcglasso_cpp_I
    start_pcglasso_cpp_I <- Sys.time()
    res_pcglasso_cpp_I <- blockwise_optimization(
      S, lambda, alpha,
      tolerance = tolerance
    )
    end_pcglasso_cpp_I <- Sys.time()
    time_pcglasso_cpp_I_mat[m, i] <- as.numeric(difftime(end_pcglasso_cpp_I, start_pcglasso_cpp_I, units = "secs"))

    Sinv <- res_pcglasso_cpp_I$Sinv
    res_val_pcglasso_cpp_I_mat[m, i] <- pcglasso_goal_function(
      S, lambda, alpha, cov2cor(Sinv), sqrt(diag(Sinv))
    )

    # pcglasso_cpp_C
    start_pcglasso_cpp_C <- Sys.time()
    Q0 <- cov2cor(solve(S))
    res_pcglasso_cpp_C <- blockwise_optimization(
      S, lambda, alpha, Q = Q0, Q_inv = solve(Q0),
      tolerance = tolerance
    )
    end_pcglasso_cpp_C <- Sys.time()
    time_pcglasso_cpp_C_mat[m, i] <- as.numeric(difftime(end_pcglasso_cpp_C, start_pcglasso_cpp_C, units = "secs"))

    Sinv <- res_pcglasso_cpp_C$Sinv
    res_val_pcglasso_cpp_C_mat[m, i] <- pcglasso_goal_function(
      S, lambda, alpha, cov2cor(Sinv), sqrt(diag(Sinv))
    )
  }
}
# trim is from each side
trim <- if (m >= 10) { 0.1 } else { 0 }
time_pcglasso_C        <- apply(time_pcglasso_C_mat,        2, mean, trim = trim)
res_val_pcglasso_C     <- apply(res_val_pcglasso_C_mat,     2, mean, trim = trim)
time_pcglasso_I        <- apply(time_pcglasso_I_mat,        2, mean, trim = trim)
res_val_pcglasso_I     <- apply(res_val_pcglasso_I_mat,     2, mean, trim = trim)
time_pcglassoFast_I    <- apply(time_pcglassoFast_I_mat,    2, mean, trim = trim)
res_val_pcglassoFast_I <- apply(res_val_pcglassoFast_I_mat, 2, mean, trim = trim)
time_pcglassoFast_C    <- apply(time_pcglassoFast_C_mat,    2, mean, trim = trim)
res_val_pcglassoFast_C <- apply(res_val_pcglassoFast_C_mat, 2, mean, trim = trim)
time_pcglasso_cpp_I    <- apply(time_pcglasso_cpp_I_mat,    2, mean, trim = trim)
res_val_pcglasso_cpp_I <- apply(res_val_pcglasso_cpp_I_mat, 2, mean, trim = trim)
time_pcglasso_cpp_C    <- apply(time_pcglasso_cpp_C_mat,    2, mean, trim = trim)
res_val_pcglasso_cpp_C <- apply(res_val_pcglasso_cpp_C_mat, 2, mean, trim = trim)

df <- data.frame(
  time  = c(time_pcglasso_C, time_pcglasso_I,
            time_pcglassoFast_I, time_pcglassoFast_C,
            time_pcglasso_cpp_I, time_pcglasso_cpp_C),
  value = c(res_val_pcglasso_C, res_val_pcglasso_I,
            res_val_pcglassoFast_I, res_val_pcglassoFast_C,
            res_val_pcglasso_cpp_I, res_val_pcglasso_cpp_C),
  which = factor(c(
    rep("pcglasso start C", length(time_pcglasso_C)),
    rep("pcglasso start I", length(time_pcglasso_I)),
    rep("pcglassoFast start I", length(time_pcglassoFast_I)),
    rep("pcglassoFast start C", length(time_pcglassoFast_C)),
    rep("pcglasso_cpp start I", length(time_pcglasso_cpp_I)),
    rep("pcglasso_cpp start C", length(time_pcglasso_cpp_C))
  ))
)
df$alg  <- sub(" (start C|start I)$", "", df$which)
df$init <- ifelse(grepl("start C", df$which), "C", "I")

eps <- 1e-5 # chosen for the plot to be prettier
df$value <- df$value - min(df$value) + eps

options(scipen = 2)
ggplot(df, aes(x = time, y = value, color = alg, shape = init)) +
  geom_point(size = 5) +
  scale_shape_manual(values = c("C" = 20, "I" = 8)) +
  theme_minimal(base_size = 14) +
  scale_y_log10() +
  #scale_x_log10() +
  #expand_limits(x = 0) +
  labs(
    title = "PCGLASSO vs PCGLASSOFast vs PCGLASSOcpp",
    subtitle = paste0(
      "p = ", p, ", n = ", n, ", corr = -", cor_modifier, "/sqrt(p)"
    )
  )
``

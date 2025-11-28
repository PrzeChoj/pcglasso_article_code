# devtools::install_github("JackStorrorCarter/PCGLASSO")
# devtools::install_github("PrzeChoj/pcglassoFast")
# devtools::install_github("PrzeChoj/PCGLASSOcpp")

library(PCGLASSO)
library(pcglassoFast)
library(PCGLASSOcpp)
library(MASS)
library(ggplot2)

set.seed(1234)

pcglasso_goal_function <- function(S, lambda, alpha, Sinv) {
  delta_matrix <- cov2cor(Sinv)
  theta_diag <- sqrt(diag(Sinv))
  p <- nrow(delta_matrix)

  theta_mat <- diag(theta_diag)
  log_det   <- determinant(delta_matrix)$modulus
  attr(log_det, "logarithm") <- NULL
  quad_term <- sum(diag(S %*% theta_mat %*% delta_matrix %*% theta_mat))
  l1_pen    <- lambda * sum(abs(delta_matrix - diag(p)))

  -log_det - 2 * (1 - alpha) * sum(log(theta_diag)) + quad_term + l1_pen
}

M <- 30
p <- 50
n <- p*2
rho <- 0.2
lambda <- rho
c_parameter <- 1
alpha <- 1 - c_parameter
tolerance_list <- exp(seq(log(0.01), log(0.00001), length.out = 12))
n_tol <- length(tolerance_list)
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

run_with_obj <- function(fun) {
  start <- Sys.time()
  Sinv  <- fun()
  end   <- Sys.time()

  elapsed <- as.numeric(difftime(end, start, units = "secs"))
  value <- pcglasso_goal_function(S, lambda, alpha, Sinv)
  list(time = elapsed, value = value)
}

counter <- 0
pb <- txtProgressBar(min = counter, max = n_tol * M, style = 3)
time_pcglasso_C_mat        <- matrix(NA, nrow = M, ncol = n_tol)
res_val_pcglasso_C_mat     <- matrix(NA, nrow = M, ncol = n_tol)
time_pcglasso_I_mat        <- matrix(NA, nrow = M, ncol = n_tol)
res_val_pcglasso_I_mat     <- matrix(NA, nrow = M, ncol = n_tol)
time_pcglassoFast_I_mat    <- matrix(NA, nrow = M, ncol = n_tol)
res_val_pcglassoFast_I_mat <- matrix(NA, nrow = M, ncol = n_tol)
time_pcglassoFast_C_mat    <- matrix(NA, nrow = M, ncol = n_tol)
res_val_pcglassoFast_C_mat <- matrix(NA, nrow = M, ncol = n_tol)
time_pcglasso_cpp_I_mat    <- matrix(NA, nrow = M, ncol = n_tol)
res_val_pcglasso_cpp_I_mat <- matrix(NA, nrow = M, ncol = n_tol)
time_pcglasso_cpp_C_mat    <- matrix(NA, nrow = M, ncol = n_tol)
res_val_pcglasso_cpp_C_mat <- matrix(NA, nrow = M, ncol = n_tol)
for (m in seq_len(M)) {
  for (i in seq_along(tolerance_list)) {
    tolerance <- tolerance_list[i]
    tol_pcglasso <- tolerance * pcglasso_tolerance_modifier

    counter <- counter + 1
    setTxtProgressBar(pb, counter)

    # pcglasso_C
    res <- run_with_obj(function() {
      pcglasso(
        S, lambda, c_parameter,
        threshold = tol_pcglasso
      )
    })
    time_pcglasso_C_mat[m, i]    <- res$time
    res_val_pcglasso_C_mat[m, i] <- res$value

    # pcglasso_I
    res <- run_with_obj(function() {
      pcglasso(
        S, lambda, c_parameter,
        Theta_start = diag(nrow(S)),
        threshold   = tol_pcglasso
      )
    })
    time_pcglasso_I_mat[m, i]    <- res$time
    res_val_pcglasso_I_mat[m, i] <- res$value

    # pcglassoFast_I
    res <- run_with_obj(function() {
      pcglassoFast(
        S, lambda = lambda, alpha = alpha,
        tolerance = tolerance
      )$Sinv
    })
    time_pcglassoFast_I_mat[m, i]    <- res$time
    res_val_pcglassoFast_I_mat[m, i] <- res$value

    # pcglassoFast_C
    res <- run_with_obj(function() {
      pcglassoFast(
        S, lambda = lambda, alpha = alpha,
        R = cov2cor(solve(S)),
        tolerance = tolerance
      )$Sinv
    })
    time_pcglassoFast_C_mat[m, i]    <- res$time
    res_val_pcglassoFast_C_mat[m, i] <- res$value

    # pcglasso_cpp_I
    res <- run_with_obj(function() {
      blockwise_optimization(
        S, lambda, alpha,
        tolerance = tolerance
      )$Sinv
    })
    time_pcglasso_cpp_I_mat[m, i]    <- res$time
    res_val_pcglasso_cpp_I_mat[m, i] <- res$value

    # pcglasso_cpp_C
    res <- run_with_obj(function() {
      Q0 <- cov2cor(solve(S))
      blockwise_optimization(
        S, lambda, alpha,
        Q = Q0, Q_inv = solve(Q0),
        tolerance = tolerance
      )$Sinv
    })
    time_pcglasso_cpp_C_mat[m, i]    <- res$time
    res_val_pcglasso_cpp_C_mat[m, i] <- res$value
  }
}
close(pb)
trim <- if (M >= 10) { 0.1 } else { 0 } # trim is from each side

mean_trim <- function(mat) { apply(mat, 2, mean, trim = trim) }

time_pcglasso_C        <- mean_trim(time_pcglasso_C_mat)
res_val_pcglasso_C     <- mean_trim(res_val_pcglasso_C_mat)
time_pcglasso_I        <- mean_trim(time_pcglasso_I_mat)
res_val_pcglasso_I     <- mean_trim(res_val_pcglasso_I_mat)
time_pcglassoFast_I    <- mean_trim(time_pcglassoFast_I_mat)
res_val_pcglassoFast_I <- mean_trim(res_val_pcglassoFast_I_mat)
time_pcglassoFast_C    <- mean_trim(time_pcglassoFast_C_mat)
res_val_pcglassoFast_C <- mean_trim(res_val_pcglassoFast_C_mat)
time_pcglasso_cpp_I    <- mean_trim(time_pcglasso_cpp_I_mat)
res_val_pcglasso_cpp_I <- mean_trim(res_val_pcglasso_cpp_I_mat)
time_pcglasso_cpp_C    <- mean_trim(time_pcglasso_cpp_C_mat)
res_val_pcglasso_cpp_C <- mean_trim(res_val_pcglasso_cpp_C_mat)

df <- data.frame(
  time  = c(time_pcglasso_C, time_pcglasso_I,
            time_pcglassoFast_I, time_pcglassoFast_C,
            time_pcglasso_cpp_I, time_pcglasso_cpp_C),
  value = c(res_val_pcglasso_C, res_val_pcglasso_I,
            res_val_pcglassoFast_I, res_val_pcglassoFast_C,
            res_val_pcglasso_cpp_I, res_val_pcglasso_cpp_C),
  which = factor(c(
    rep("pcglasso start C", n_tol),
    rep("pcglasso start I", n_tol),
    rep("pcglassoFast start I", n_tol),
    rep("pcglassoFast start C", n_tol),
    rep("pcglasso_cpp start I", n_tol),
    rep("pcglasso_cpp start C", n_tol)
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
  expand_limits(x = 0) +
  labs(
    title = "PCGLASSO vs PCGLASSOFast vs PCGLASSOcpp",
    subtitle = paste0(
      "p = ", p, ", n = ", n, ", corr = -", cor_modifier, "/sqrt(p)"
    )
  )

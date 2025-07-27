library(pcglassoFast)
library(parallel)

source("./raw_experiments/irrepresentability/plot_irrep_data_functions.R")

n_cores <- 7

N <- 100 # 100 -> 10 minutes

p <- 15
a <- 60
b <- 1
c <- 1.5

number_of_ns <- 20
n_list <- floor(exp(seq(floor(log(10^2)), floor(log(10^5)), length.out = number_of_ns)))

lambda_list <- (1:200)*0.0001


Khub <- function(p, a, b, c) {
  K <- matrix(0, nrow = p, ncol = p)
  K[, 1] <- K[1, ] <- c
  diag(K) <- b
  K[1, 1] <- a
  K
}
K <- Khub(p, a, b, c)

irrepGLASSO(K) # 3
irrepPCGLASSO(K) # 0.29



# 30 seconds
start_time <- Sys.time()
results_all_jobs <- perform_jobs(
  K, n_list, lambda_list, N, n_cores
)
end_time <- Sys.time()
print(end_time - start_time)

my_ggplot <- my_plot(results_all_jobs) +
  annotate(
    "text", x = 5000, y = 0.95,
    label = paste0("IRR_PCG(K) = ", round(irrepPCGLASSO(K), 2)),
    size = 4, color = "red"
  ) +
  annotate(
    "text", x = 5000, y = 0.1,
    label = paste0("IRR_G(K) = ", round(irrepGLASSO(K), 2)),
    size = 4, color = "blue"
  ) +
  scale_x_log10(breaks = c(100, 1000, 10000, 60000), labels = c("1e2", "1e3", "1e4", "6e4"))

my_ggplot




#####
# second example; GLASSO good; PCGLASSO bad
K_2 <- structure(c(
  0.01, 0.01, -0.03, 0.1, -0.09,
  0.01, 3.53, 0.59, 0.94, 0,
  -0.03, 0.59, 3.13, 0, -0.22,
  0.1, 0.94, 0, 2.22, -0.29,
  -0.09, 0, -0.22, -0.29, 1.42
), dim = c(5L, 5L))
irrepPCGLASSO(K_2) # 1.57
irrepGLASSO(K_2) # 0.81


# 20 seconds
start_time <- Sys.time()
results_all_jobs_2 <- perform_jobs(
  K_2, n_list, lambda_list, N, n_cores
)
end_time <- Sys.time()
print(end_time - start_time)

my_ggplot_2 <- my_plot(results_all_jobs_2) +
  annotate(
    "text", x = 3000, y = 0.35,
    label = paste0("IRR_PCG(K) = ", round(irrepPCGLASSO(K_2), 2)),
    size = 4, color = "red"
  ) +
  annotate(
    "text", x = 500, y = 0.8,
    label = paste0("IRR_G(K) = ", round(irrepGLASSO(K_2), 2)),
    size = 4, color = "blue"
  ) +
  scale_x_log10(breaks = c(100, 1000, 10000, 60000), labels = c("1e2", "1e3", "1e4", "6e4"))

my_ggplot_2

library(space)
library(pcglassoFast)
library(PCGLASSO)
library(glasso)
library(parallel)
library(Matrix)
library(dplyr)
library(tidyverse)

source("estimation_methods.R")
source("simulation_functions.R")
# source("./raw_experiments/simulation/estimation_methods.R")
# source("./raw_experiments/simulation/simulation_functions.R")


# 2) Common settings (you can adjust these as needed)
set.seed(1)
graphics.off()
generate.pcglasso=T
split.train        <- 0.7      # not used directly below, but could be referenced in your functions
ns                 <- c(200)
sim                <- 7#21
nlambda            <- 50
mc_cores           <- 7
alpha_grid         <- 0
#lambda.min.ratio   <- 0.1
pcglasso_tolerance <- 0.01
lambda.min.ratio = 0.01

# 3) Load the appropriate Q‐matrix depending on generate.pcglasso
if (!generate.pcglasso) {
  data(Q_simulated_glasso)      # from pcglassoFast or simulation_functions.R
  Q <- Q_simulated_glasso
} else {
  data(Q_simulated_pcglasso)
  Q <- Q_simulated_pcglasso
}
n <- ns
# ensure symmetry
Q <- (Q + t(Q)) / 2
p <- ncol(Q)
L <- Cholesky(Matrix(forceSymmetric(Q), sparse = TRUE), LDL = FALSE, perm = TRUE)
z <- matrix(rnorm(n * p), nrow = p, ncol = n)
x <- solve(L, solve(L, z, system = "P"), system = "Lt")
x <- Matrix::solve(L, x, system = "Pt")
data <- as.matrix(t(x))
S_full  <- cov(data)
split_train = 0.7
 lambda.min.ratio = 0.01
n_train <- floor(split_train * n)
idx     <- sample.int(n, n_train)
train   <- data[idx, , drop = FALSE]
test    <- data[-idx, , drop = FALSE]
S_train <- cov(train)
S_test  <- cov(test)
n_test  <- nrow(test)

lam_max <- max(abs(S_train - diag(diag(S_train))))
lam_min <- lambda.min.ratio * lam_max
lambdas <- exp(seq(log(lam_max), log(lam_min), length.out = nlambda))

carter_path <- pcglasso_path_carter(S_full,
                                    lambdas,
                                    pcglasso_tolerance
)
loss_cg_full <- evaluate_objective_path(carter_path, Sigma = S_full,
                                        n = n, gamma = 0.5)
plot(lambdas,loss_cg_full$loglik)

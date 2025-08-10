is_K_positive_definite <- function(K) {
  min(eigen(K, TRUE, TRUE)$values) > 0.00001
}

# 1 -> hub
Khub <- function(p, a, b, c_val) {
  K <- matrix(0, nrow = p, ncol = p)
  K[, 1] <- K[1, ] <- c_val
  diag(K) <- b
  K[1, 1] <- a
  K
}

# 2 -> AR
KAR <- function(p, a, b, c_val) {
  K <- matrix(0, nrow = p, ncol = p)
  diag(K) <- b
  K[1, 1] <- a
  K[p, p] <- a

  for (i in 2:p) {
    K[i-1, i] <- K[i, i-1] <- c_val
  }

  K
}

# 3 -> Sanger
data("Sanger", package = "pcglassoFast")

# 4 -> Stock Market
data("stockdata", package = "huge")
my_data <- stockdata$data
log_returns <- log(my_data[-1, ] / my_data[-nrow(my_data), ])
log_returns <- sweep(log_returns, 1, rowMeans(log_returns))


# 1 -> hub
# 2 -> AR
# 3 -> Sanger
# 4 -> Stock Market
get_S <- function(p, which_experiment) {
  stopifnot(which_experiment %in% c(1, 2, 3, 4))
  n <- 400

  switch (which_experiment,
    { # hub
      a <- 4 * (p-1) * 1.05
      b <- 1
      c_val <- 2
      K <- Khub(p, a, b, c_val)
      stopifnot(is_K_positive_definite(K))

      X <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = solve(K))
      cov(X)
    },
    { # AR
      a <- 4
      b <- 1
      c_val <- 0.43
      K <- KAR(p, a, b, c_val)
      stopifnot(is_K_positive_definite(K))

      X <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = solve(K))
      cov(X)
    },
    { # Sanger
      stopifnot(p <= 100)

      p_max <- 123
      selected <- sample(p_max, p)
      Sanger_selected <- Sanger[, selected]

      if (nrow(Sanger_selected) < n) {
        X <- Sanger_selected[sample(nrow(Sanger_selected), n, replace = TRUE), ]
      } else {
        X <- Sanger_selected[sample(nrow(Sanger_selected), n), ]
      }

      cov(X) + diag(0.0001, p)
    },
    { # Stock Market
      stopifnot(p <= 300)

      selected <- sample(ncol(log_returns), p)
      log_returns_selected <- log_returns[, selected]

      col_means <- colMeans(log_returns_selected)
      col_sds <- apply(log_returns_selected, 2, sd)
      outlier_matrix <- abs(sweep(log_returns_selected, 2, col_means)) > 5 * rep(col_sds, each = nrow(log_returns_selected))
      keep_rows <- !apply(outlier_matrix, 1, any)
      clean_data <- log_returns_selected[keep_rows, ]

      if (nrow(clean_data) < n) {
        X <- clean_data[sample(nrow(clean_data), n, replace = TRUE), ]
      } else {
        X <- clean_data[sample(nrow(clean_data), n), ]
      }

      cov(X) + diag(0.0001, p)
    }
  )
}

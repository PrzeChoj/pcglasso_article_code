ud <- function(S, alpha, lambdas){
  R0_big_lambda <- diag(nrow(S))
  R0_inv_big_lambda <- solve(R0_big_lambda)

  pcglassoPath(S, alpha, lambdas = lambdas, R0 = R0_big_lambda, R0_inv = R0_inv_big_lambda)
}
du <- function(S, alpha, lambdas){
  R0_inv_small_lambda <- cov2cor(S)
  R0_small_lambda <- solve(R0_inv_small_lambda)

  pcglassoPath(S, alpha, lambdas = rev(lambdas), R0 = R0_small_lambda, R0_inv = R0_inv_small_lambda)
}
du2 <- function(S, alpha, lambdas){
  R0 <- cov2cor(solve(S))
  R0_inv <- solve(R0)

  pcglassoPath(S, alpha, lambdas = rev(lambdas), R0 = R0, R0_inv = R0_inv)
}
du3 <- function(S, alpha, lambdas){
  R0 <- cov2cor(S)
  R0_inv <- solve(R0)

  pcglassoPath(S, alpha, lambdas = rev(lambdas), R0 = R0, R0_inv = R0_inv)
}
du4 <- function(S, alpha, lambdas){
  R0 <- diag(nrow(S)) # identity

  pcglassoPath(S, alpha, lambdas = rev(lambdas), R0 = R0, R0_inv = R0)
}

udu <- function(S, alpha, lambdas){
  R0_big_lambda <- diag(nrow(S))
  R0_inv_big_lambda <- solve(R0_big_lambda)

  ud_path <- pcglassoPath(S, alpha, lambdas = lambdas, R0 = R0_big_lambda, R0_inv = R0_inv_big_lambda)
  pcglassoPath(S, alpha, lambdas = lambdas, R0 = ud_path$R_path[[nlambda]], R0_inv = ud_path$Ri_path[[nlambda]])
}
dud <- function(S, alpha, lambdas){
  R0_inv_small_lambda <- cov2cor(S)
  R0_small_lambda <- solve(R0_inv_small_lambda)

  du_path <- pcglassoPath(S, alpha, lambdas = rev(lambdas), R0 = R0_small_lambda, R0_inv = R0_inv_small_lambda)
  pcglassoPath(S, alpha, lambdas = rev(lambdas), R0 = du_path$R_path[[nlambda]], R0_inv = du_path$Ri_path[[nlambda]])
}

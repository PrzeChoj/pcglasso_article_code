#install_github("PrzeChoj/pcglassoFast")
library(pcglassoFast)

load("Sangerdata.Rdata")
X <- data
X$TargetID <- NULL
X <- as.matrix(X)
Sigma.Sanger <- cov(X)
min(cov2cor(Sigma.Sanger)) # 0.96

(p <- nrow(Sigma.Sanger)) # 210

# Thw following will not work, as the max_iter_R_inner = 1000000000 is too small and 10000000000 is too big (NAs prodused in Fortran)
ans_Sanger <- pcglassoFast(
  Sigma.Sanger, 0.6, 0, max_iter = 10
  , max_iter_R_inner = 10, max_iter_R_outer = 10, tol_R_outer = 1e-3,
  verbose = TRUE
)
R_curr <- ans_Sanger$R
Rinv_curr <- ans_Sanger$R_inv
D_curr <- ans_Sanger$D

ans_Sanger <- pcglassoFast(
  Sigma.Sanger, 0.2, 0, max_iter = 100,R=R_curr, R_inv = Rinv_curr, D = D_curr
  , max_iter_R_inner = 10, max_iter_R_outer = 10, tol_R_outer = 1e-3,
  verbose = TRUE
)

#ans_Sanger



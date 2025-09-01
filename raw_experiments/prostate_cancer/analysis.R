

source("../estimation_function.R")

df_reduced <- read.csv("312_df_reduced.csv")
n.lambda <- 10

df_reduced <- df_reduced[, -1]
df_genes <- df_reduced[, 7 + 1:200]
S <- cov(df_genes)
n <- dim(df_genes)[1]


#lam_max <- max(abs(S - diag(diag(S))))
#lam_min <- 0.1 * lam_max
lambdas <- exp(seq(log( 2), log(0.55), length.out = n.lambda))
lambdas <- exp(seq(log( 1), log(0.55), length.out = n.lambda))

PC <- solve(S)
R0 = cov2cor(PC)
pcglasso.res <- estimator_pcglasso(S,n, lambdas,0,gamma,max_edge_fraction=0.7, R_start = R0, max_iter=400, max_iter_R_outer=400)
#pcglasso.res_0 <- estimator_pcglasso(S,n, lambdas,0,gamma,max_edge_fraction=0.7)
print(pcglasso.res$path.loss$nEdges)
print(image(Matrix(pcglasso.res$path.all$R_path[[1]],sparse=T),sparse=T))


alpha.res <- apply(pcglasso.res$path.all$R_path[[1]], 2, function(x){sum((abs(x))^(1))}) - 1
# which.max(alpha.res)

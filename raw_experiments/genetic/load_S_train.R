library(PCGLASSO)
library(glmnet)

# just some preselection
load("./expermient/sangerdata/Sangerdata.Rdata")
# load("Sangerdata.Rdata")
y <- as.numeric(data[24266, -1])
x <- t(as.matrix(data[-24266, -1]))

colnames(x) <- data[-24266, 1]
# Gosia preprocess
quan <- quantile(as.vector(x), probs = 0.25)
out1 <- which(apply(x, 2, max) < quan)
x1 <- x[, -out1]

range <- apply(x1, 2, max) - apply(x1, 2, min)
out2 <- which(range < 2)
x2 <- x1[, -out2]
xx <- x2
p <- dim(xx)[2]
n <- dim(xx)[1]
X <- as.matrix(xx)
X <- matrix(X, ncol = ncol(X))


n <- dim(X)[1]


objl <- cv.glmnet(X, y, standardize = TRUE, intercept = TRUE)
ind <- which(coefficients(objl, s = exp(-5))[-1] != 0)

X <- cbind(X[, ind], y)
X <- X[, -18] #remove outlier


S.train <- cov(X)

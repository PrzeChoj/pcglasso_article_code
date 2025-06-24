#https://osf.io/dcebh/?view_only=
load("Deidentified Data.Rdata")
graphics.off()
library(Matrix)
library(glasso)

library(qgraph)
library(bootnet)


## ---- Sample Characteristics ----
table(data$Gender)
mean(data$Age)
sd(data$Age)
range(data$Age)
table(data$Race)
table(data$Race)/length(data$Race)


networkdata <- data[, c("CES1_1", "CES1_2", "CES1_3", "CES1_4", "CES2_1", "CES2_3", "CES2_4", "PCL5A_1", "PCL5A_2", "PCL5A_3", "PCL5A_4", "PCL5B_1", "PCL5B_2", "PCL5B_3", "PCL5B_4", "PCL5C_1", "PCL5C_2", "PCL5C_3", "PCL5C_4", "PCL5D_1", "PCL5D_2", "PCL5D_3", "PCL5D_4", "PCL5E_1", "PCL5E_2", "PCL5E_3", "PCL5E_4")]
colnames(networkdata) = c("identity", "reference", "lifestory", "colored", "lifechange", "future", "turning", "memories", "dreams", "flash", "upset", "physior", "avoidint", "avoidext", "amnesia", "beliefs", "blame", "negfeel", "lossint", "distant", "numb", "anger", "reckless", "hyper", "startle", "concen", "sleep")

set.seed(5)
n <- dim(networkdata)[1]
n.train <- round(0.5*n)
ind.train <- sample(n,n.train)
ind.val  <- setdiff(1:n, c(ind.train))
data.train <- networkdata[ind.train,]
data.val <- networkdata[ind.val,]
Sigma.train <- (cov(as.matrix(data.train)))
Sigma.val <- (cov(as.matrix(data.val)))

res <- EBICglasso(cov(networkdata), dim(networkdata)[1],returnAllResults=T)

n.alpha <- 1
n.lambda <- 20
#' create grid
alpha.grid <- c(0,-0.05,0.07)#seq(-0.1,0.02, length.out=n.alpha)
alpha.grid <- c(0.)
S_max <-  max(abs(S.train[row(S.train)!= col(S.train)]))

lambdas <- c(exp(seq(log(0.01*S_max), log(0.2*S_max), length.out= n.lambda)))
hyperparam <- c()
for(i in 1:n.alpha){
  alphas <- rep(alpha.grid[i], n.lambda)
  hyperparam <- rbind(hyperparam, cbind(lambdas,alphas))
}



alpha.grid=c(-0.1,-0.05,0,0.05,0.1,0.15,0.2)
alpha.0 <- which(alpha.grid==0)
n <- dim(networkdata)[1]

EBIC_glasso <- pcglasso.path_EBIC(networkdata,
                                  lambdas,
                               alpha.grid,
                               gamma = 0.,
                               verbose=T)


alpha.max <- which.min(apply(EBIC_glasso$EBIC,2,min))
res.alpha <- matrix(rep(alpha.grid, each=n.lambda), nrow=n.lambda)
res.lambda <- matrix(rep(lambdas, length(alpha.grid)), nrow=n.lambda)
results <- data.frame( loglik=c(EBIC_glasso$EBIC), alpha=c(res.alpha), lambda=c(res.lambda))
p <- ggplot(results, aes(x = alpha, y = lambda, z = loglik)) +
  geom_contour_filled(aes(fill = after_stat(level)))  +
  labs(title = "EBIC Surface",
       x = "alpha",
       y = "lambda",
       fill = "EBIC surface")
print(p)

Sigma.train = cov(networkdata)
n.train = dim(networkdata)[1]
resultslasso <- estimateNetwork(networkdata, default = "EBICglasso", corMethod = 'cor', corArgs = list(method='spearman'), tuning = 0)
lambdas.glasso <- resultslasso$results$lambda
Precision.hat.array <- array(NA, dim = c(dim(S.train)[1],dim(S.train)[1],length(lambdas.glasso)))
ns <- ns_<- matrix(NA,nrow=length(lambdas.glasso) )
loglik <- matrix(NA,nrow=length(lambdas.glasso)  )
Precision_corr.hat.array <- array(NA, dim = c(dim(S.train)[1],dim(S.train)[1],length(lambdas.glasso)))
ns_corr <- matrix(NA,nrow=length(lambdas.glasso) )
loglik_corr <- matrix(NA,nrow=length(lambdas.glasso)  )
D<-diag(1/sqrt(diag(Sigma.train)))
for(i in 1:length(lambdas.glasso)){
  PL<-glasso(Sigma.train,lambdas.glasso[i], penalize.diagonal=F,thr=1.0e-6,nobs=n.train);
  PL_corr<-glasso(cov2cor(Sigma.train),lambdas.glasso[i], penalize.diagonal=F,thr=1.0e-6,nobs=n.train);
  wi <- PL$wi
  PL_corr$wi <- D%*%PL_corr$wi%*%D
  Precision.hat.array[,,i] <- wi
  Precision_corr.hat.array[,,i] <- PL_corr$wi
  loglik[i] <- PL$loglik
  loglik_corr[i] <- PL_corr$loglik
  L.P <- chol(wi)
  loglik[i] <- sum(log(diag(L.P))) - 0.5 * sum(diag(wi%*%Sigma.val))
  ns[i] <- (sum(PL$wi!=0)-dim(Sigma.train)[1])/2

  ns_[i] <- (sum(resultslasso$results$results$wi[,,i]!=0)-dim(Sigma.train)[1])/2


  L.P <- chol(PL_corr$wi)
  loglik_corr[i] <- sum(log(diag(L.P))) - 0.5 * sum(diag(PL_corr$wi%*%Sigma.val))
  ns_corr[i] <- (sum(PL_corr$wi!=0)-dim(Sigma.train)[1])/2
}
gamma=0.5
EBIC <- -2*n*loglik + log(n)*ns + 4*gamma*ns*log(dim(networkdata)[2])
EBIC_corr <- -2*n*loglik_corr + log(n)*ns_corr + 4*gamma*ns_corr*log(dim(networkdata)[2])
plot(ns_corr,EBIC_corr,type='l',xlab='number of edges',ylab='EBIC',
     xlim = c(min(min(ns),min(res.ns)), 1.2*max(max(ns),max(res.ns))),ylim=c(9000,12000))
lines(EBIC_glasso$nedges[,alpha.max],EBIC_glasso$EBIC[,alpha.max],col='red')
lines(EBIC_glasso$nedges[,alpha.0 ],EBIC_glasso$EBIC[,alpha.0 ],col='blue')
legend(245, y=-10.1,c('GLASSO','alpha=0','alpha_max'), col=c('black','blue','red'),lty=c(1,1,1))

if(0){


min.GLASSO <- which.max(rowMeans(loglik))
PL<-glasso(Sigma.train,lambdas.glasso[min.GLASSO], penalize.diagonal=F,thr=1.0e-6,nobs=n.train);
sol_1 <- Precision.hat.array[,,min.GLASSO]
sol_1 <- (sol_1 + t(sol_1))/2
sol_1 <- abs(sol_1)
diag(sol_1) <- diag(sol_1)*0
newlasso <- qgraph(sol_1, layout = "spring", labels = colnames(networkdata), groups = list("CES" = 1:7, "PCL-5" = 8:27), color = c("lightblue", "sandybrown"))
min.PCGLASSO <- which.max(cross_val_glasso$loglik)

sol_2 <- cross_val_glasso$Qs[,,min.PCGLASSO]
diag(sol_2) <- diag(sol_2)*0
newPCGLASSO <- qgraph(sol_2, layout = "spring", labels = colnames(networkdata), groups = list("CES" = 1:7, "PCL-5" = 8:27), color = c("lightblue", "sandybrown"))
}


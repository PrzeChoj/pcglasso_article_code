#https://osf.io/dcebh/?view_only=
load("Deidentified Data.Rdata")
source("../estimation_function.R")
graphics.off()
library(Matrix)
library(glasso)
library(ggplot2)
library(akima)
library(qgraph)
library(bootnet)

gamma  = 0.

lambda.min.ratio = 0.05
n.lambda <- 100
## ---- Sample Characteristics ----
table(data$Gender)
mean(data$Age)
sd(data$Age)
range(data$Age)
table(data$Race)
table(data$Race)/length(data$Race)


networkdata <- data[, c("CES1_1", "CES1_2", "CES1_3", "CES1_4", "CES2_1", "CES2_3", "CES2_4", "PCL5A_1", "PCL5A_2", "PCL5A_3", "PCL5A_4", "PCL5B_1", "PCL5B_2", "PCL5B_3", "PCL5B_4", "PCL5C_1", "PCL5C_2", "PCL5C_3", "PCL5C_4", "PCL5D_1", "PCL5D_2", "PCL5D_3", "PCL5D_4", "PCL5E_1", "PCL5E_2", "PCL5E_3", "PCL5E_4")]
colnames(networkdata) = c("identity", "reference", "lifestory", "colored", "lifechange", "future", "turning", "memories", "dreams", "flash", "upset", "physior", "avoidint", "avoidext", "amnesia", "beliefs", "blame", "negfeel", "lossint", "distant", "numb", "anger", "reckless", "hyper", "startle", "concen", "sleep")
X <- networkdata

covX <- cov(X)
n <- dim(X)[1]

lam_max <- 0.2*max(abs(covX - diag(diag(covX))))
lam_min <- lambda.min.ratio * lam_max
lambdas <- exp(seq(log(lam_max), log(lam_min), length.out = n.lambda))

#alpha.grid=c(-0.1,-0.05,0,0.05,0.1,0.15,0.2)
alpha.grid <- sort(unique(c(
  seq(-0.15, -0.01, length.out = 20),
  0,
  seq(0.01, 0.1, length.out = 3)
)))
pcglasso.res <- estimator_pcglasso(covX,n, lambdas,alpha.grid,gamma,max.edge.fraction=0.7)


bic_df <- do.call(rbind, lapply(names(pcglasso.res$path.loss), function(a) {
  loss <- pcglasso.res$path.loss[[a]]
  lambda <- pcglasso.res$path.all[[a]]$lambda  # extract lambda sequence
  data.frame(
    alpha = as.numeric(a),
    lambda = lambda,
    BIC = loss$BIC_gamma
  )
}))

interp_grid <- with(bic_df, interp(
  x = alpha,
  y = lambda,
  z = BIC,
  xo = seq(min(alpha), max(alpha), length = 300),
  yo = seq(min(lambda), max(lambda), length = 300),
  duplicate = "mean"
))
bic_interp_df <- expand.grid(
  alpha  = interp_grid$x,
  lambda = interp_grid$y
)
bic_interp_df$BIC <- as.vector(interp_grid$z)
bic_interp_df <- na.omit(bic_interp_df)

# Min BIC point
bic_min <- bic_df[which.min(bic_df$BIC), ]


bic_vals <- bic_interp_df$BIC
q_breaks <- quantile(bic_vals, probs = seq(0, 1, length.out = 10), na.rm = TRUE)

fig <- ggplot(bic_interp_df, aes(x = alpha, y = lambda, z = BIC)) +
  geom_contour_filled(breaks = q_breaks) +
  geom_point(data = bic_min, aes(x = alpha, y = lambda), color = "red", size = 3) +
  annotate("text", x = bic_min$alpha, y = bic_min$lambda,
           label = "min BIC", hjust = -0.1, vjust = -1, size = 4, color = "red") +
  labs(
    title = "BIC Surface over (α, λ)",
    x = expression(alpha),
    y = expression(lambda)
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    panel.grid = element_blank()
  )
print(fig)

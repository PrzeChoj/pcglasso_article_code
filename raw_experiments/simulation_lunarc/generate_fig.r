

source("../estimation_function.R")
library(pcglassoFast)



data(Q_simulated_glasso)
p_glasso   <- make_plot_matrix(Q_simulated_glasso, "GLASSO precision", x_lab="", y_lab="")
ggsave(
  "glasso_precision.png",
  plot = p_glasso, width = 7, height = 4
)
data(Q_simulated_glasso)
p_pcglasso   <- make_plot_matrix(Q_simulated_glasso, "PCGLASSO precision", x_lab="", y_lab="")
ggsave(
  "pcglasso_precision.png",
  plot = p_pcglasso, width = 7, height = 4
)

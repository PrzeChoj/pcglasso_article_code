
library(pcglassoFast)
source("../estimation_function.R")

data(Q_simulated_glasso)      # from pcglassoFast or simulation_functions.R
#Q_glasso <- Matrix(Q_simulated_glasso, sparse = TRUE)
data(Q_simulated_pcglasso)
#Q_pcglasso <- Matrix(Q_simulated_pcglasso, sparse = TRUE)

p_nonhub     <- make_plot_matrix(Q_simulated_glasso, "NON-HUB",
                                   base_size = 6, title_size = 8,
                                   axis_title_size = 6, axis_text_size = 5,
                                   tick_length_pt = 1)

p_hub   <- make_plot_matrix(Q_simulated_pcglasso, "HUB",
                                   base_size = 6, title_size = 8,
                                   axis_title_size = 6, axis_text_size = 5,
                                   tick_length_pt = 1)


ggsave(
  "adj_nonhub.png",
  plot = p_nonhub,
  width =3,
  height = 3
)

ggsave(
  "adj_hub.png",
  plot = p_hub,
  width =3,
  height = 3
)

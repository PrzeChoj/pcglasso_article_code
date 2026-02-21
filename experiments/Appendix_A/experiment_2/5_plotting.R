library(dplyr)
library(ggplot2)
library(ggpattern)

timing_raw <- read.csv("./experiments/Appendix_A/experiment_2/res_data/timing_raw.csv")
timing_summary <- read.csv("./experiments/Appendix_A/experiment_2/res_data/timing_summary.csv")

timing_raw <- timing_raw %>%
  mutate(
    p = factor(p, levels = sort(unique(p))),
    structure = factor(structure),
    solver = factor(solver),
    start = factor(start)
  )

timing_summary <- timing_summary %>%
  mutate(
    p = factor(p, levels = sort(unique(p))),
    structure = factor(structure),
    solver = factor(solver),
    start = factor(start)
  )

p1 <- ggplot(timing_summary, aes(x = p, y = time_trimmed_mean, color = solver, linetype = start, group = interaction(solver, start))) +
  geom_line() +
  geom_point() +
  #scale_y_log10() +
  facet_wrap(~ structure, scales = "free_y") +
  labs(x = "p", y = "Trimmed mean time (trim = 0.1)") +
  theme_bw()

p1
ggsave("./experiments/Appendix_A/experiment_2/res_data/plot_time_trimmed_mean_vs_p.png", p1, width = 10, height = 5, dpi = 300)

p2 <- ggplot(
  timing_raw,
  aes(
    x = p, y = time,
    fill = solver,
    pattern = start
  )
) +
  geom_boxplot_pattern(
    position = position_dodge2(preserve = "single"),
    outlier.alpha = 0.3,
    pattern_fill = "black",
    pattern_angle = 45,
    pattern_density = 0.4,
    pattern_spacing = 0.03
  ) +
  scale_pattern_manual(values = c(I = "none", C = "stripe")) +
  scale_y_log10() +
  facet_wrap(~ structure, scales = "free_y") +
  labs(x = "p", y = "Time per run") +
  theme_bw()

p2
ggsave("./experiments/Appendix_A/experiment_2/res_data/plot_time_boxplot_vs_p.png", p2, width = 10, height = 5, dpi = 300)

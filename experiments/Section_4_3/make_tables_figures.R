# make_tables_figures.R
#
# Loads saved simulation results and produces:
#
#   Tables 1–2  (Section 4.3):   hub graph (Q_simulated_pcglasso)
#   Tables 3–4  (Appendix B.2):  non-hub graph (Q_simulated_glasso)
#
# And the following figures saved to ../../outputs/:
#
#   Section_4_3_pcglasso_rmse.png      RMSE breakdown, hub graph
#   Section_4_3_pcglasso_rates.png     FP/FN rates, hub graph
#   Section_4_3_glasso_rmse.png        RMSE breakdown, non-hub graph
#   Section_4_3_glasso_rates.png       FP/FN rates, non-hub graph
#   Section_4_3_true_precision.png     True precision matrices side-by-side
#
# Run run_simulation.R first to generate results/results_pcglasso.rds
# and results/results_glasso.rds.

library(ggplot2)
library(patchwork)
library(pcglassoFast)

source("./experiments/Section_4_3/simulation_functions.R")


# ---- Helpers ----------------------------------------------------------------

check_results_exist <- function(path) {
  if (!file.exists(path)) {
    stop(
      "Results file not found: ", path, "\n",
      "Please run run_simulation.R first (or supply the .rds files from the HPC run)."
    )
  }
}

# Inline version of make_plot_matrix() from raw_experiments/estimation_function.R
plot_precision_matrix <- function(Q, title) {
  df <- as.data.frame(as.table(Q != 0))
  colnames(df) <- c("Row", "Column", "Value")
  df$Row    <- as.numeric(df$Row)
  df$Column <- as.numeric(df$Column)
  df$Value  <- as.numeric(df$Value)

  nnz_pct <- round(
    100 * (sum(Q != 0) - nrow(Q)) / (nrow(Q)^2 - nrow(Q)), 0
  )

  ggplot(df, aes(x = Column, y = Row, fill = Value)) +
    geom_tile(color = "white") +
    scale_fill_gradient(low = "white", high = "blue", name = "Non-Zero") +
    labs(
      title = paste0(title, ", non-zero = ", nnz_pct, "%"),
      x = NULL, y = NULL
    ) +
    scale_x_continuous(breaks = seq(0, ncol(Q), by = 20)) +
    scale_y_reverse(breaks = seq(0, nrow(Q), by = 20)) +
    coord_fixed() +
    theme_minimal(base_size = 6) +
    theme(
      panel.grid        = element_blank(),
      legend.position   = "none",
      panel.background  = element_rect(fill = "white", color = NA),
      plot.background   = element_rect(fill = "white", color = NA),
      plot.title        = element_text(size = 8, hjust = 0.5)
    )
}

# Print a clean summary table with means rounded to 3 significant figures
print_table <- function(summary_df, caption) {
  cat("\n", caption, "\n", strrep("-", nchar(caption)), "\n", sep = "")
  cols <- c("n", "method",
            "rmse", "rmse_diag", "rmse_offdiag_zero", "rmse_offdiag_nonzero",
            "false_pos_rate", "false_neg_rate", "timing")
  tbl <- summary_df[, intersect(cols, names(summary_df))]
  tbl[, setdiff(names(tbl), c("n", "method"))] <-
    signif(tbl[, setdiff(names(tbl), c("n", "method"))], 3)
  tbl <- tbl[order(tbl$n, tbl$method), ]
  print(tbl, row.names = FALSE)
}


# ---- Load results -----------------------------------------------------------

res_pcglasso <- readRDS("./experiments/Section_4_3/results/results_pcglasso.rds")
message("Loaded results for ", length(unique(res_pcglasso$n)),
        " sample sizes and ", length(unique(res_pcglasso$method)), " methods (hub graph).")

res_glasso   <- readRDS("./experiments/Section_4_3/results/results_glasso.rds")
message("Loaded results for ", length(unique(res_glasso$n)),
        " sample sizes and ", length(unique(res_glasso$method)), " methods (non-hub graph).")


# ---- Summarise --------------------------------------------------------------

sum_pcglasso <- summarize_plot_results(res_pcglasso)
sum_glasso   <- summarize_plot_results(res_glasso)


# ---- Tables -----------------------------------------------------------------

print_table(sum_pcglasso$table,
            "Table 1 & 2 — Hub graph (Q_simulated_pcglasso) — Section 4.3")
print_table(sum_glasso$table,
            "Table 3 & 4 — Non-hub graph (Q_simulated_glasso) — Appendix B.2")

message("\nDone.")

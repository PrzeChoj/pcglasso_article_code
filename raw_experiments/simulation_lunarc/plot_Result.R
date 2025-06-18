library(ggplot2)
library(pcglassoFast)
library(Matrix)
source("simulation_functions.R")

data(Q_simulated_glasso)      # from pcglassoFast or simulation_functions.R
Q_glasso <- Matrix(Q_simulated_glasso, sparse = TRUE)
data(Q_simulated_pcglasso)
Q_pcglasso <- Matrix(Q_simulated_pcglasso, sparse = TRUE)

pdf('NonHub.pdf')
image(Q_glasso>0, main="Non structure",sub="")
dev.off()
pdf('Hub.pdf')
image(Q_pcglasso>0, main="Hub structure",sub="")
dev.off()


res_glasso <- readRDS(file = "results_glasso.rds")
res_pcglasso <- readRDS(file = "results_pcglasso.rds")

table_glasso <- summarize_plot_results(res_glasso)

good_rows <- suppressWarnings(!is.na(as.numeric(res_pcglasso$frob_norm)))

res_pcglasso_clean <- res_pcglasso[good_rows, ]
numeric_cols <- c("n",
  "frob_norm", "rmse", "frob_diag", "rmse_diag",
  "frob_offdiag_zero", "rmse_offdiag_zero",
  "frob_offdiag_nonzero", "rmse_offdiag_nonzero",
  "true_non0_rate", "false_non0_rate",
  "false_0_rate", "true_0_rate", "timing"
)

for (col in numeric_cols) {
  res_pcglasso_clean[[col]] <- as.numeric(res_pcglasso_clean[[col]])
}
table_pcglasso <- summarize_plot_results(res_pcglasso_clean)
#' Output two compact LaTeX tables: one for RMSE metrics, one for timing, methods as subrows, with customizable header.
#'
#' @param tab Data.frame containing columns: n, method, rmse, rmse_diag, rmse_offdiag_nonzero, timing
#' @param ns Numeric vector of sample sizes to use as columns, e.g. c(100, 500, 1000)
#' @param header Table header as a character string, e.g. "Hub Structure" (set NULL for none)
#' @return None (prints LaTeX to the console)
#' @export
latex_metric_and_timing_tables <- function(tab, ns = c(100, 500, 1000), header = "Hub Structure") {
  clean_method <- function(m) {
    m <- gsub("_cv$", " CV", m)
    m <- gsub("_bic$", " BIC", m)
    m <- gsub("_", " ", m)
    m
  }
  sig2 <- function(x) {
    if (is.na(x) || !is.finite(x)) return("")
    as.character(signif(x, 2))
  }

  metrics <- c("rmse", "rmse_diag", "rmse_offdiag_nonzero")
  metric_labels <- c(
    rmse = "RMSE",
    rmse_diag = "Diag RMSE",
    rmse_offdiag_nonzero = "Off-diag (NZ) RMSE"
  )
  methods <- unique(tab$method)
  tab_sub <- tab[tab$n %in% ns, ]

  # Precompute mins for metrics
  metric_mins <- lapply(metrics, function(metric) {
    sapply(ns, function(n) {
      vals <- tab_sub[tab_sub$n == n, metric]
      min(vals, na.rm = TRUE)
    })
  })
  names(metric_mins) <- metrics

  # Precompute mins for timing
  timing_min <- sapply(ns, function(n) {
    vals <- tab_sub[tab_sub$n == n, "timing"]
    min(vals, na.rm = TRUE)
  })

  get_val_bold <- function(method, n, metric, min_vec) {
    val <- tab_sub[tab_sub$method == method & tab_sub$n == n, metric]
    if (length(val) == 0) return("")
    sval <- sig2(val)
    minval <- min_vec[as.character(n)]
    if (is.na(val) || is.na(minval)) return(sval)
    if (abs(val - minval) < 1e-8) sprintf("\\textbf{%s}", sval) else sval
  }

  n_cols <- 2 + length(ns)

  # Table 1: Metrics (not timing)
  cat("\\begin{table}[ht]\n\\centering\n\\small\n")
  cat("\\begin{tabular}{ll", paste(rep("c", length(ns)), collapse=""), "}\n", sep="")
  cat("\\toprule\n")
  if (!is.null(header)) {
    cat(sprintf("\\multicolumn{%d}{c}{\\textbf{%s}} \\\\\n", n_cols, header))
  }
  cat("Metric & Method")
  for (n in ns) cat(sprintf(" & $n=%d$", n))
  cat(" \\\\\n\\midrule\n")
  for (i in seq_along(metrics)) {
    metric <- metrics[i]
    min_vec <- setNames(metric_mins[[metric]], as.character(ns))
    first <- TRUE
    for (method in methods) {
      method_label <- clean_method(method)
      row_label <- if (first) metric_labels[[metric]] else ""
      first <- FALSE
      cat(row_label, " & ", method_label, sep="")
      for (n in ns) cat(" & ", get_val_bold(method, n, metric, min_vec), sep="")
      cat(" \\\\\n")
    }
    cat("\\midrule\n")
  }
  cat("\\bottomrule\n\\end{tabular}\n")
  cat(sprintf("\\caption{RMSE summary for each method and sample size%s. The smallest (best) value for each metric and $n$ is bolded.}\n",
              if (!is.null(header)) paste0(" (", header, ")") else ""))
  cat("\\end{table}\n\n")

  # Table 2: Timing only
  min_vec <- setNames(timing_min, as.character(ns))
  cat("\\begin{table}[ht]\n\\centering\n\\small\n")
  cat("\\begin{tabular}{l", paste(rep("c", length(ns)), collapse=""), "}\n", sep="")
  cat("\\toprule\n")
  if (!is.null(header)) {
    cat(sprintf("\\multicolumn{%d}{c}{\\textbf{%s}} \\\\\n", 1 + length(ns), header))
  }
  cat("Method")
  for (n in ns) cat(sprintf(" & $n=%d$", n))
  cat(" \\\\\n\\midrule\n")
  for (method in methods) {
    method_label <- clean_method(method)
    cat(method_label)
    for (n in ns) cat(" & ", get_val_bold(method, n, "timing", min_vec), sep="")
    cat(" \\\\\n")
  }
  cat("\\bottomrule\n\\end{tabular}\n")
  cat(sprintf("\\caption{Computation time (seconds) for each method and sample size%s. The smallest value for each $n$ is bolded.}\n",
              if (!is.null(header)) paste0(" (", header, ")") else ""))
  cat("\\end{table}\n")
}


latex_metric_and_timing_tables(table_pcglasso$table, ns = c(200, 500, 1000,5000), header = "Hub Structure")
latex_metric_and_timing_tables(table_glasso$table, ns = c(200, 500, 1000,5000), header = "Non Hub Structure")


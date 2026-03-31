#!/usr/bin/env Rscript
#
# make_tables_auto.R
#
# Load all results_*.rds and replay_*.rds in the current folder, de-duplicate
# by (graph_type, n, rep, method), average across reps, and print LaTeX tables.
#
# Usage:
#   Rscript make_tables_auto.R
#
# Example:
#   digits <- 2L
#

digits <- 2L

infer_graph_type <- function(path) {
  base <- tolower(basename(path))
  if (grepl("non_hub", base)) return("non_hub")
  if (grepl("hub", base)) return("hub")
  NA_character_
}

read_results_file <- function(path) {
  obj <- readRDS(path)
  if (!is.data.frame(obj)) {
    message("Skipping non-data.frame: ", path)
    return(NULL)
  }
  if (!"graph_type" %in% names(obj)) {
    obj$graph_type <- infer_graph_type(path)
  }
  obj
}

read_replay_file <- function(path) {
  obj <- readRDS(path)
  if (!is.list(obj)) {
    message("Skipping non-list replay: ", path)
    return(NULL)
  }
  graph_type <- infer_graph_type(path)
  out <- list()
  for (entry in obj) {
    if (!is.list(entry) || !isTRUE(entry$ok)) {
      next
    }
    if (is.null(entry$results)) {
      next
    }
    rep_id <- if (!is.null(entry$rep)) entry$rep else NA_integer_
    n_val <- if (!is.null(entry$n)) entry$n else NA_integer_
    gt <- if (!is.null(entry$graph_type) && !is.na(entry$graph_type)) entry$graph_type else graph_type
    for (sel in names(entry$results)) {
      cmp <- entry$results[[sel]]
      if (!is.data.frame(cmp)) {
        next
      }
      if (!"n" %in% names(cmp)) {
        cmp$n <- n_val
      }
      cmp$rep <- rep_id
      cmp$graph_type <- gt
      out[[length(out) + 1]] <- cmp
    }
  }
  if (length(out)) {
    do.call(rbind, out)
  } else {
    NULL
  }
}

results_files <- list.files(pattern = "^results_.*\\.rds$")
results_files <- results_files[!grepl("_part[0-9]+of[0-9]+\\.rds$", results_files)]
results_files <- results_files[!grepl("^results_.*_part", results_files)]

replay_files <- list.files(pattern = "^replay_.*\\.rds$")

results_list <- lapply(results_files, read_results_file)
replay_list <- lapply(replay_files, read_replay_file)

results_list <- Filter(Negate(is.null), results_list)
replay_list <- Filter(Negate(is.null), replay_list)

all_results <- if (length(results_list)) do.call(rbind, results_list) else data.frame()
replay_results <- if (length(replay_list)) do.call(rbind, replay_list) else data.frame()

if (!nrow(all_results) && !nrow(replay_results)) {
  stop("No usable results found in the current directory.")
}

if (!nrow(all_results)) {
  all_results <- replay_results
} else if (nrow(replay_results)) {
  all_results <- rbind(all_results, replay_results)
}

if (!"rep" %in% names(all_results)) {
  all_results$rep <- NA_integer_
}
if (!"graph_type" %in% names(all_results)) {
  all_results$graph_type <- NA_character_
}

key <- paste(all_results$graph_type, all_results$n, all_results$rep, all_results$method, sep = "|")
all_results <- all_results[!duplicated(key, fromLast = TRUE), ]

required <- c("graph_type", "n", "method", "rmse", "rmse_diag", "rmse_offdiag_nonzero", "timing")
missing <- setdiff(required, names(all_results))
if (length(missing)) {
  stop("Missing required columns: ", paste(missing, collapse = ", "))
}

fmt_num <- function(x, digits) {
  if (is.na(x) || !is.finite(x)) return("")
  val <- round(x, digits)
  if (digits == 0 && val == 0) {
    return("0")
  }
  formatC(val, format = "f", digits = digits)
}

clean_method <- function(m) {
  m <- gsub("_cv$", " CV", m)
  m <- gsub("_bic$", " BIC", m)
  m <- gsub("_", " ", m)
  m
}

format_n <- function(n) {
  trimws(format(n, big.mark = "\\,", scientific = FALSE))
}

latex_timing_longtable <- function(tab, ns, header, label) {
  cat("\\begin{longtable}[]{@{}l", paste(rep("c", length(ns)), collapse = ""), "@{}}\n", sep = "")
  cat(sprintf("\\caption{Computation time (seconds) for each method and sample size (%s).\\label{%s}}\\\\\n", header, label))
  cat("\\toprule\n")
  cat("Method")
  for (n in ns) cat(sprintf(" & $n=%s$", format_n(n)))
  cat(" \\\\\n\\midrule\n\\endfirsthead\n\n")
  cat(sprintf("\\multicolumn{%d}{c}{{\\tablename\\ \\thetable{} -- continued from previous page}} \\\\\n", 1 + length(ns)))
  cat("\\toprule\n")
  cat("Method")
  for (n in ns) cat(sprintf(" & $n=%s$", format_n(n)))
  cat(" \\\\\n\\midrule\n\\endhead\n\n")
  cat("\\midrule \\multicolumn{", 1 + length(ns), "}{r}{{Continued on next page}} \\\\\n", sep = "")
  cat("\\endfoot\n\n\\bottomrule\n\\endlastfoot\n\n")

  timing_min <- sapply(ns, function(n) {
    vals <- tab[tab$n == n, "timing"]
    min(vals, na.rm = TRUE)
  })
  min_vec <- setNames(timing_min, as.character(ns))
  methods <- sort(unique(tab$method))
  for (method in methods) {
    cat(clean_method(method))
    for (n in ns) {
      val <- tab[tab$method == method & tab$n == n, "timing"]
      if (length(val) == 0) {
        cat(" & ")
        next
      }
      sval <- fmt_num(val, digits)
      minval <- min_vec[as.character(n)]
      if (!is.na(val) && !is.na(minval) && abs(val - minval) < 1e-8) {
        cat(" & ", sprintf("\\textbf{%s}", sval), sep = "")
      } else {
        cat(" & ", sval, sep = "")
      }
    }
    cat(" \\\\\n")
  }
  cat("\n\\end{longtable}\n\n")
}

latex_rmse_longtable <- function(tab, ns, header, label) {
  metrics <- c("rmse", "rmse_diag", "rmse_offdiag_nonzero")
  metric_labels <- c(
    rmse = "RMSE",
    rmse_diag = "Diag RMSE",
    rmse_offdiag_nonzero = "Off-diag (NZ) RMSE"
  )

  cat("\\begin{longtable}[]{@{}ll", paste(rep("c", length(ns)), collapse = ""), "@{}}\n", sep = "")
  cat(sprintf("\\caption{RMSE summary for each method and sample size (%s).\\label{%s}}\\\\\n", header, label))
  cat("\\toprule\\noalign{}\n")
  cat("Metric & Method")
  for (n in ns) cat(sprintf(" & $n=%s$", format_n(n)))
  cat(" \\\\\n\\midrule\\noalign{}\n\\endfirsthead\n\n")
  cat(sprintf("\\multicolumn{%d}{c}{{\\tablename\\ \\thetable{} -- continued from previous page}} \\\\\n", 2 + length(ns)))
  cat("\\toprule\\noalign{}\n")
  cat("Metric & Method")
  for (n in ns) cat(sprintf(" & $n=%s$", format_n(n)))
  cat(" \\\\\n\\midrule\\noalign{}\n\\endhead\n\n")
  cat("\\midrule \\multicolumn{", 2 + length(ns), "}{r}{{Continued on next page}} \\\\\n", sep = "")
  cat("\\endfoot\n\n\\bottomrule\n\\endlastfoot\n\n")

  methods <- sort(unique(tab$method))
  for (metric in metrics) {
    metric_mins <- sapply(ns, function(n) {
      vals <- tab[tab$n == n, metric]
      min(vals, na.rm = TRUE)
    })
    min_vec <- setNames(metric_mins, as.character(ns))
    first <- TRUE
    for (method in methods) {
      row_label <- if (first) metric_labels[[metric]] else ""
      first <- FALSE
      cat(row_label, " & ", clean_method(method), sep = "")
      for (n in ns) {
        val <- tab[tab$method == method & tab$n == n, metric]
        if (length(val) == 0) {
          cat(" & ")
          next
        }
        sval <- fmt_num(val, digits)
        minval <- min_vec[as.character(n)]
        if (!is.na(val) && !is.na(minval) && abs(val - minval) < 1e-8) {
          cat(" & ", sprintf("\\textbf{%s}", sval), sep = "")
        } else {
          cat(" & ", sval, sep = "")
        }
      }
      cat(" \\\\\n")
    }
    if (metric != tail(metrics, 1)) {
      cat("\n\\midrule\\noalign{}\n")
    }
  }
  cat("\n\\end{longtable}\n\n")
}

graph_types <- sort(unique(all_results$graph_type))
graph_types <- graph_types[!is.na(graph_types)]
if (length(graph_types) == 0) {
  stop("No graph_type found in results; use filenames with hub/non_hub or add graph_type column.")
}

for (gt in graph_types) {
  sub <- all_results[all_results$graph_type == gt, ]
  num_cols <- c("n", "rmse", "rmse_diag", "rmse_offdiag_nonzero", "timing")
  for (col in num_cols) {
    sub[[col]] <- suppressWarnings(as.numeric(sub[[col]]))
  }
  sub$method <- as.character(sub$method)
  summary_df <- aggregate(
    cbind(rmse, rmse_diag, rmse_offdiag_nonzero, timing) ~ n + method,
    data = sub,
    FUN = mean,
    na.rm = TRUE
  )
  ns <- sort(unique(summary_df$n))
  header <- if (gt == "hub") "Hub Structure" else "Non Hub Structure"
  tag <- if (gt == "hub") "hub" else "nonhub"

  latex_timing_longtable(summary_df, ns, header, paste0("tab:", tag, "time"))
  latex_rmse_longtable(summary_df, ns, header, paste0("tab:", tag, "rmse"))
}

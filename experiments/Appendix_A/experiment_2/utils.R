get_baseline <- function(p, structure) {
  row <- baseline_best_value[
    baseline_best_value$p == p & baseline_best_value$structure == structure,
    ,
    drop = FALSE
  ]
  if (nrow(row) != 1L) stop("Baseline not found/unique")
  row$best_value[1]
}

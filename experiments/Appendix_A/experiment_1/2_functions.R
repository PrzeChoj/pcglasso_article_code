source("./experiments/Appendix_A/utils.R")
source("./experiments/Appendix_A/experiment_1/0_parameters.R")

load("./experiments/Appendix_A/experiment_1/res_data/instances.RData")

run_one <- function(row) {
  p <- row$p
  cor_modifier <- row$cor_modifier
  lambda <- row$lambda
  alpha <- row$alpha
  K_structure <- row$K_structure
  solver <- row$solver
  starting_point <- row$starting_point
  tolerance <- row$tolerance

  S <- instances[[K_structure]][[as.character(p)]][[as.character(cor_modifier)]]
  if (is.null(S)) stop("Missing S for: ", K_structure, " p=", p, " cor=", cor_modifier)

  res <- tryCatch(
    {
      t0 <- proc.time()[["elapsed"]]
      f_end <- value_after_optimization(S, solver, starting_point, tolerance, lambda, alpha)
      t1 <- proc.time()[["elapsed"]]
      list(status = "ok", time = t1 - t0, f_end = f_end, error = NA_character_)
    },
    error = function(e) {
      list(status = "fail", time = NA_real_, f_end = NA_real_, error = conditionMessage(e))
    }
  )

  c(as.list(row), res)
}

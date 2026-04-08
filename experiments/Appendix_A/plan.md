# Experiment: comparison (pcglasso vs pcglassoFast_Dual vs pcglassoFast_Primal)

---

## Fixed settings

### Graph structures
- `K_structure in {hub_1, hub_09, AR2, random}`

Interpretation:
- `hub_1`: hub graph with partial correlation `-1 / sqrt(p)`
- `hub_09`: hub graph with partial correlation `-0.9 / sqrt(p)`
- `AR2`: AR2 graph
- `random`: random graph

### Dimensions
- `p in {50, 100, 150, 200}`

### Data generation
For each pair `(K_structure, p)`:
- sample size: `n = 2p`
- build `K_star = build_K_star(p, K_structure)`
- construct `S = S_from_K_star(K_star, n)`

These matrices `S` are generated once, with global seed:
- `set.seed(1234)`

and saved in:
- `./experiments/Appendix_A/res_data/instances.RData`

Then the same `S` is reused across all runs for a given `(K_structure, p)`.

### Regularization
- `lambda in {0.1, 0.2}`
- `alpha in {0, 0.5}`

### Tolerance grid
- `tolerance_list = exp(seq(log(1e-2), log(1e-8), length.out = 12))`

### Repeats
- repeats per configuration: `M = 100`

### Algorithms compared
The experiment compares:
- `pcglasso`
- `pcglassoFast_Dual`
- `pcglassoFast_Primal`

Each algorithm is run with two starting points:
- `starting_point in {C, I}`

Hence there are 6 method variants in total:
1. `pcglasso / C`
2. `pcglasso / I`
3. `pcglassoFast_Dual / C`
4. `pcglassoFast_Dual / I`
5. `pcglassoFast_Primal / C`
6. `pcglassoFast_Primal / I`

### Tolerance mapping
For:
- `pcglassoFast_Dual`
- `pcglassoFast_Primal`

the tolerance is used directly.

For:
- `pcglasso`

the tolerance passed to the solver is:
- `tol_pcglasso = tolerance * pcglasso_tolerance_multiplier`

with
- `pcglasso_tolerance_multiplier = 1000`

### Objective evaluation
Each run returns the final objective value:
- `f_end = value_after_optimization(S, solver, starting_point, tolerance, lambda, alpha)`

and the elapsed runtime.

---

## Step 1 - Generate instances

Generate and save one input matrix `S` for each `(K_structure, p)` pair.

For each:
- `K_structure in {hub_1, hub_09, AR2, random}`
- `p in {50, 100, 150, 200}`

do:
1. set `n = 2p`
2. construct `K_star = build_K_star(p, K_structure)`
3. compute `S = S_from_K_star(K_star, n)`

Save all such objects in:
- `./experiments/Appendix_A/res_data/instances.RData`

---

## Step 2 - Run simulations on the full grid

Construct the full parameter grid:
- `p x lambda x alpha x K_structure x solver x starting_point x tolerance`

Then replicate each row `M = 100` times by adding repetition index:
- `m in {1, ..., M}`

For every row of this replicated grid:
1. load `S = instances[[K_structure]][[as.character(p)]]`
2. if `solver == "pcglasso"`, multiply tolerance by `1000`
3. run optimization
4. record:
   - `status`
   - `time`
   - `f_end`
   - `error`

Save all raw runs into one file:
- `./experiments/Appendix_A/res_data/raw_M100.csv`

---

## Step 3 - Aggregate results

From the raw file, aggregate over repetitions separately for each:
- `(p, lambda, alpha, K_structure, solver, starting_point, tolerance)`

Compute:
- `time_median = median(time)`
- `f_end = mean(f_end)`

For `f_end`, this averaging is only a formal aggregation step. In the simulation output, for each fixed
`(p, lambda, alpha, K_structure, solver, starting_point, tolerance)`, all `M = 100` values of `f_end`
were the same, which matches the intended deterministic behavior of the solvers on a fixed input matrix `S`.
Therefore, `mean(f_end)` is effectively identical to each individual replicated value.

and save the summary to:
- `./experiments/Appendix_A/res_data/summary_M100.csv`

Thus the summary file contains one row per:
- method variant x tolerance x `(p, lambda, alpha, K_structure)`

---

## Step 4 - Compute benchmark best values

For each:
- `(p, K_structure, lambda, alpha)`

compute a reference `best_value`.

The benchmark is obtained as:
- `best_value_calculated = get_best_value(S, p, K_structure, lambda, alpha)`
- `best_value_from_simulations = min(f_end over all raw runs for this configuration)`

and finally:
- `best_value = min(best_value_calculated, best_value_from_simulations)`

Save the resulting table in:
- `./experiments/Appendix_A/res_data/group_keys_with_best_value_M100.csv`

---

## Step 5 - Plot type 1: time vs objective gap

Using `summary_M100.csv` and `group_keys_with_best_value_M100.csv`, for each:
- `(p, K_structure, lambda, alpha)`

construct a scatter plot over all solvers, starting points, and tolerances.

For each row define:
- `time = time_median`
- `value = f_end`
- `value_shifted = pmax(value - best_value, 1e-10)`

Plot:
- `x = time`
- `y = value_shifted` on log-scale
- color by algorithm
- shape by starting point

Save one PNG per configuration:
- `./experiments/Appendix_A/plots/type_1/plot_{K_structure}_p{p}_lambda{lambda}_alpha{alpha}.png`

This results in 64 plots.

---

## Step 6 - Plot type 2: tolerance needed to reach near-best objective

Define threshold:
- `thr_f_diff_to_best = 1e-5`

Using raw results, for each:
- `(p, K_structure, lambda, alpha, solver, starting_point)`

find the largest tolerance from `tolerance_list`, such that:
- `f_end - best_value < 1e-5`

Call this:
- `tol_found`

Then keep only runs with:
- `tolerance == tol_found`

From these filtered runs:
- compute median runtime across repetitions

This produces one runtime summary per:
- `(p, K_structure, lambda, alpha, solver, starting_point)`

### Type 2a - Median time vs p
For each graph structure:
- create faceted line plots of `median_time` vs `p`
- facets correspond to `(lambda, alpha)`
- separate curves for algorithm and starting point

Save:
- `./experiments/Appendix_A/plots/type_2/median_time_vs_p_{K_structure}.png`

### Type 2b - Violin plots
For each graph structure:
- create violin plots of runtime distributions versus `p`
- facets correspond to `(lambda, alpha)`
- separate violins for algorithm / starting point combinations

Save:
- `./experiments/Appendix_A/plots/type_2/violin_time_vs_p_{K_structure}.png`

---

## Final outputs

The experiment produces:

### Data files
- `instances.RData`
- `raw_M100.csv`
- `summary_M100.csv`
- `group_keys_with_best_value_M100.csv`

### Plots
- type 1: one scatter plot per `(p, K_structure, lambda, alpha)`
- type 2: one median-time plot and one violin plot per `K_structure`

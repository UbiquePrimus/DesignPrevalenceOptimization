# Spatial surveillance simulations

R scripts for comparing spatial sampling designs and design-prevalence choices in simulated host populations. Infection status is obtained by thresholding a Matérn Gaussian field. The simulations summarize prevalence estimation, coverage of positive hosts, and nearest-neighbour classification accuracy.

This README describes the refactored generalized simulation and helper scripts, together with the supplied `exemple_pr_Francois.R`. The example needs the small setup changes described below before it can run independently.

## Files

| File | Purpose |
| --- | --- |
| `Correlation_PrevEstim_Efficiency_Generalized.R` | Configures and runs simulation scenarios, produces comparisons, and saves results. |
| `functions_Correlation_PrevEstim_Efficiency_Generalized.R` | Computes excursion-set quantities, EFSA sample sizes, and surveillance metrics. |
| `exemple_pr_Francois.R` | Illustrates fixed and perimeter-corrected prevalence using one Gaussian field; exports two PNG figures. |

Keep both generalized scripts in the same directory. The main script locates its helper relative to the script location; the output directory defaults to the current working directory.

## Requirements

Install R and the packages used by the generalized simulation:

```r
install.packages(c("GpGp", "spatstat.geom", "spatstat.random",
                  "ggplot2", "gridExtra"))
```

The generalized script uses explicit package namespaces and does not install or attach packages automatically. Its plotting code uses `geom_line(linewidth = ...)`, so use ggplot2 3.4.0 or later.

The supplied example also calls `library(pdftools)`, although none of its computations use that package. Either install `pdftools` or remove that unused import.

## Run the generalized simulation

### Default run

From the directory containing the scripts:

```bash
Rscript Correlation_PrevEstim_Efficiency_Generalized.R
```

This runs all nine configured scenarios, saves `Generalized_Simulation_Results.RData`, and draws three comparison layouts. In noninteractive execution, graphics output follows R's graphics-device settings; the script does not explicitly export named plot files.

The default is **two replicates per prevalence/range combination**, intended for a preliminary run. It still includes 30 range values and nine scenarios. Increase the replicate count for analysis.

### Configure a run in R

Sourcing the refactored main script defines functions without starting simulations:

```r
source("Correlation_PrevEstim_Efficiency_Generalized.R")

config <- simulation_config()
config$seed <- 123
config$N_sim_test <- 100
config$output_dir <- "results"

results <- main(config, make_plots = FALSE)
```

For a smaller initial run:

```r
config <- simulation_config()
config$seed <- 123
config$N <- 400
config$rho_test <- c(10, 50, 100)
config$prev <- c(0.01, 0.05)
config$N_sim_test <- 2
config$output_dir <- "trial_results"

results <- main(config, make_plots = FALSE)
```

`main()` still runs all nine scenarios. Set `config$delay_gpsim <- 1` to generate a separate Gaussian field for every replicate within a scenario. Larger values reuse fields within blocks.

If the helper is elsewhere, pass its location explicitly:

```r
results <- main(config,
                helper_file = "path/to/functions_Correlation_PrevEstim_Efficiency_Generalized.R",
                make_plots = FALSE)
```

### Main parameters

| Parameter | Default | Meaning |
| --- | --- | --- |
| `N` | 10,000 | Nominal host population used for pattern generation and EFSA budgets. |
| `W_size` | 10,000 | Side length of the square observation window. |
| `prev` | 0.0005, 0.005, 0.01, 0.05 | Target true prevalences for empirical field thresholding. |
| `rho_test` | 30 values from 1.5 to 100 | True Gaussian-field ranges. |
| `nu` | 1.5 | Matérn smoothness; must exceed 1 for the spectral-moment calculation. |
| `cl` | 0.9 | Confidence parameter in the sample-size approximation. |
| `ms` | 0.8 | Detection-sensitivity parameter in that approximation. |
| `dp_base` | 0.01 | Fixed design prevalence for the standard EFSA comparison. |
| `rho_agg`, `dp_agg` | 100, 0.01 | Reference range and exceedance probability defining the target expected perimeter. |
| `N_sim_test` | 2 | Replicates per prevalence/range combination. |
| `delay_gpsim` | 50 | Number of replicates sharing a Gaussian field. |
| `strata_grid` | 4 | Number of strata per spatial dimension. |
| `risk_weight` | 0.9 | Mixture weight on normalized risk scores; the remainder is uniform. |
| `stage1_fraction` | 0.6 | Fraction of the adaptive sampling budget allocated to initial SRS. |
| `tau_noise` | 0.4 | Standard deviation of log-range noise in the optional noisy scenario. |
| `seed` | `NULL` | Optional seed set once at the beginning of `main()`. |
| `output_dir` | `"."` | Directory for saved results. |
| `verbose` | `TRUE` | Whether to print simulation progress. |

## Run one scenario

Load dependencies and helpers before using the simulation function directly:

```r
source("Correlation_PrevEstim_Efficiency_Generalized.R")
check_dependencies()
load_simulation_helpers(
  "functions_Correlation_PrevEstim_Efficiency_Generalized.R"
)

config <- simulation_config()
config$rho_test <- c(10, 50, 100)
config$prev <- c(0.01, 0.05)
set.seed(123)  # Direct function calls do not set config$seed automatically.

design_table <- make_design_table(config)
results_one <- compute_cor_prev_estim_efficiency_generalized(
  dp_df = design_table,
  N_sim = 100,
  host_pattern = "regular",
  sampling_design = "stratified",
  rho_error_scenario = "oracle",
  delay_gpsim = 1,
  config = config
)
```

For direct calls, `N_sim`, `strata_grid`, `tau_noise`, and `delay_gpsim` are function arguments. `main()` forwards their corresponding configuration values; direct calls use the function defaults unless explicitly overridden.

### Available options

| Component | Options |
| --- | --- |
| Hosts | `regular`, `poisson`, `clustered`, `repulsive` |
| Sampling | `srs`, `stratified`, `risk_based`, `two_stage` |
| Range error | `oracle`, `minus_50`, `minus_20`, `plus_20`, `plus_50`, `noisy` |

The nine default scenarios are regular/SRS/oracle; Poisson/SRS/oracle; clustered/SRS/oracle; regular/stratified/oracle; regular/risk-based/oracle; regular/two-stage/oracle; regular/SRS with −50% range error; regular/SRS with +50% range error; and regular/SRS using fixed design prevalence.

Other combinations can be run through `compute_cor_prev_estim_efficiency_generalized()`. The noisy scenario uses `exp(rnorm(1, log(rho_true), tau_noise))`: its median range is `rho_true`, and its mean is larger when the noise is nonzero.

## Helper calculations

| Function | Purpose |
| --- | --- |
| `l2(nu, rho)` | Second spectral moment: `1 / (2 * rho^2 * (nu - 1))`. |
| `E_area(W_vol, u)` | Expected excursion area: `W_vol * (1 - pnorm(u))`. |
| `E_perimeter(W_bound_vol, l2, W_vol, u)` | Original expected-perimeter expression, including the window-boundary term. |
| `compute_threshold_and_dp(...)` | Solves for a threshold in `[0, 100]` matching the target expected perimeter; returns `threshold` and `dp`. |
| `compute_efsa_sampsize(p, cl, N, ms)` | Original EFSA sample-size approximation, rounded upward. |
| `compute_min_radius_with_sample(mpp, sample_indices)` | Evaluates an already selected sample. |
| `compute_min_radius(mpp, sample_size)` | Selects an SRS sample and evaluates it. |

The perimeter expression is:

```r
W_bound_vol * (1 - pnorm(u)) +
  sqrt(l2) / 4 * W_vol * exp(-u^2 / 2)
```

For the square window, `W_bound_vol = 4 * W_size` and `W_vol = W_size^2`. The reference parameters define one target perimeter; changing the range changes the threshold and design prevalence needed to match it. Targets outside the attainable range on `[0, 100]` cause an error.

The helper retains legacy defaults that look for `cl`, `N`, `ms`, `nu`, and `perimeter_agg` in the calling environment. Prefer explicit arguments, as the refactored main script does.

## Results and plots

`main()` returns a named list containing the nine scenario results. Each scenario contains one data frame per true prevalence, named using the form `prev_0_01`. These same nine scenario names are saved as separate objects in the RData file, together with `simulation_settings`.

Each data frame has these columns:

| Column | Meaning |
| --- | --- |
| `rho_true`, `rho_used` | True range and range used to determine corrected design prevalence. |
| `efficiency` | Mean normalized positive-host coverage score. |
| `prev_estim` | Mean prevalence estimate returned by the helper. |
| `prev_estim_sd` | Sample standard deviation of replicate prevalence estimates; `NA` for a single replicate. |
| `balanced_accuracy` | Mean helper classification score. |
| `true_prevalence` | Target prevalence used for empirical thresholding. |
| `host_pattern`, `sampling_design`, `rho_error` | Scenario identifiers. |

Inspect or reload results:

```r
results$results_baseline[["prev_0_01"]]

saved <- new.env()
load("results/Generalized_Simulation_Results.RData", envir = saved)
saved$results_baseline[["prev_0_01"]]
saved$simulation_settings
```

`plot_all_comparisons(results)` draws comparisons of host patterns, sampling designs, and range misspecification. Each comparison has three vertically arranged plots, faceted over all configured true prevalences. The prevalence plot includes a dashed target-prevalence reference.

To export these layouts explicitly:

```r
pdf("comparisons.pdf", width = 12, height = 14)
plot_all_comparisons(results)
dev.off()
```

## Run the illustration script

The supplied `exemple_pr_Francois.R` uses 90,000 grid hosts, `W_size = 10000`, `rho_test = 10`, and `nu = 1.5`. It exports `rho10_0dot01.png` and `rho10_0dot01_correc.png`, both at 3000 × 3000 pixels and 300 dpi.

Its original helper filename is `functions_Correlation_PrevEstim_Efficiency.R`, which differs from the generalized helper supplied here. Replace that source line with:

```r
source("functions_Correlation_PrevEstim_Efficiency_Generalized.R")
```

After setting `W_size` and `nu`, and before computing the corrected prevalence, add the missing reference quantities:

```r
W_bound_vol <- 4 * W_size
W_vol <- W_size^2
rho_agg <- 100
dp_agg <- 0.01
perimeter_agg <- E_perimeter(
  W_bound_vol, l2(nu, rho_agg), W_vol, -qnorm(dp_agg)
)
```

Then run the example from the directory containing the helper:

```r
set.seed(123)
source("exemple_pr_Francois.R")
```

In this illustration, the corrected probability determines the fraction of positive marks assigned to the same Gaussian field. In the generalized simulations, true prevalence is specified separately and the corrected **design prevalence determines the sampling budget**. These are different uses of the correction.

## Interpretation and limitations

- **Host counts:** the regular pattern has `round(sqrt(N))^2` points. The Poisson-labelled pattern has exactly `N` uniform points, conditional on its count. Thomas and Strauss patterns have random counts. EFSA budgets use nominal `N`; samples are capped at the realized host count.
- **Spatial structure:** the Thomas model uses mean cluster size 10 and displacement scale `W_size / 20`. The Strauss model uses `gamma = 0.1`, which is soft inhibition rather than a strict hard core.
- **Sampling:** stratified sampling allocates approximately equally across grid strata and fills shortfalls uniformly. Risk-based sampling mixes normalized risk scores with a uniform component. Two-stage sampling starts with SRS and then prefers hosts near detected positives; without detections, the second stage is uniform.
- **Prevalence:** empirical quantile thresholding approximates the requested prevalence, with discretization and ties affecting the realized fraction. Unweighted sample prevalence is generally unsuitable as a design-unbiased estimator under unequal-probability/adaptive sampling. Selection weights are not inclusion probabilities.
- **Field reuse:** Gaussian fields are shared across prevalences within a scenario and across replicates within each block. Replicate standard deviations with shared fields do not capture variability from independent outbreaks. Scenarios generate their own host patterns and fields; comparisons are not paired across scenarios.
- **Range error:** for a noisy range scenario, one range draw is used for all replicates at a given prevalence/range combination.
- **Coverage score:** `radius = 1 - min_radius / sqrt(2 * window_area)`, where `min_radius` is the largest distance from a true positive to its nearest sampled positive. This is a simulation metric requiring the true population marks. The denominator is the diagonal length for square windows.
- **Classification:** nearest-neighbour predictions are scored only on unsampled hosts. Ties are resolved in sample order. An absent class contributes zero to its sensitivity or specificity.
- **Special cases retained:** no sampled positives gives efficiency 0, prevalence estimate 0, and balanced accuracy 0.5. A census or a sample containing every positive gives perfect coverage/classification scores and population prevalence, bypassing the usual prediction calculation. The census branch returns `disease_found = 1` even in an all-negative population.
- **Sensitivity:** `ms` affects only the sample-size approximation; sampled marks are observed without simulated diagnostic errors.
- **Resources:** Gaussian simulation and distance matrices can be expensive at large host/sample counts. The illustration uses 90,000 hosts, and the main simulation stores its Gaussian fields in memory.

## Validation and refactoring notes

The refactor removes repeated scenario and plotting code, adds explicit configuration and validation, and preserves the original helper formulas and nine saved scenario names. It corrects the Gaussian-field block boundary and singleton-index sampling behavior, and adds uniform fallbacks for degenerate risk/distance weights. It is not guaranteed to reproduce the original random-number sequence or bit-for-bit outputs.

Static checks were performed, but no R runtime was available during refactoring. R parsing, package API compatibility, and end-to-end numerical execution have not been verified. Run a small experiment locally before using the scripts for a large simulation study.

# ============================================================================
# Spatial surveillance simulation helpers
# ============================================================================
# Source this file; it does not attach packages or run a simulation.
# Spatial functions require spatstat.geom. Existing function names and result
# fields are retained so that current simulation scripts remain compatible.
#
# Parameters previously read from the calling environment can now be supplied
# explicitly. Legacy defaults still look for cl, N, ms, nu and perimeter_agg.
# Example:
#   compute_efsa_sampsize(p = 0.01, cl = 0.95, N = 90000, ms = 1)
#   compute_threshold_and_dp(40000, 1e8, rho = 10,
#                            nu = 1.5, perimeter_agg = 1000)
#
# Methodological conventions retained from the original:
# * The perimeter formula and EFSA approximation are unchanged.
# * `radius` is an efficiency score, not a distance. Its normalization uses
#   sqrt(2 * window_area), the diagonal length for a square window only.
# * The no-positive-sample branch returns balanced_accuracy = 0.5.
# * A census or a sample containing every positive returns perfect scores;
#   this bypasses nearest-neighbour classification. In the census branch,
#   disease_found remains 1 even if the population contains no positives.
# * In the every-positive branch, prev_estim is the population prevalence;
#   otherwise it is the unweighted sample prevalence. Unequal-probability
#   designs therefore require a separate weighted prevalence estimator.

# ============================================================================
# Input validation
# ============================================================================

.validate_scalar <- function(x, name, lower = -Inf, upper = Inf,
                             lower_open = FALSE, upper_open = FALSE,
                             integer = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) {
    stop(name, " must be one finite numeric value.", call. = FALSE)
  }
  below <- if (lower_open) x <= lower else x < lower
  above <- if (upper_open) x >= upper else x > upper
  if (below || above || (integer && x != floor(x))) {
    stop(name, " is outside its permitted range",
         if (integer) " or is not an integer" else "", ".", call. = FALSE)
  }
  invisible(x)
}

.validate_pattern <- function(mpp) {
  if (!requireNamespace("spatstat.geom", quietly = TRUE)) {
    stop("Install 'spatstat.geom' to use the spatial sampling functions.",
         call. = FALSE)
  }
  if (!inherits(mpp, "ppp") || mpp$n < 1L) {
    stop("mpp must be a nonempty spatstat point pattern (ppp).", call. = FALSE)
  }
  marks <- mpp$marks
  if (!(is.numeric(marks) || is.logical(marks) || is.factor(marks)) ||
      !is.null(dim(marks)) || length(marks) != mpp$n || anyNA(marks) ||
      !all(marks %in% c(0, 1))) {
    stop("mpp must have one binary mark (0 or 1) per point, without NA.",
         call. = FALSE)
  }
  if (spatstat.geom::area.owin(mpp$window) <= 0) {
    stop("The observation window must have positive area.", call. = FALSE)
  }
  invisible(mpp)
}

.validate_indices <- function(sample_indices, population_size) {
  if (!is.numeric(sample_indices) || anyNA(sample_indices) ||
      any(!is.finite(sample_indices)) ||
      any(sample_indices != floor(sample_indices)) ||
      any(sample_indices < 1 | sample_indices > population_size) ||
      anyDuplicated(sample_indices) > 0L) {
    stop("sample_indices must contain distinct integer indices in 1:mpp$n.",
         call. = FALSE)
  }
  invisible(sample_indices)
}

# ============================================================================
# Excursion-set calculations
# ============================================================================

# Second spectral moment under the original Matern parameterization.
# A finite value requires smoothness nu > 1 and range rho > 0.
l2 <- function(nu, rho) {
  .validate_scalar(nu, "nu", lower = 1, lower_open = TRUE)
  .validate_scalar(rho, "rho", lower = 0, lower_open = TRUE)
  1 / (2 * rho^2 * (nu - 1))
}

# Expected area above threshold u for a standardized Gaussian field.
# u may be a vector. All areas and lengths must use consistent spatial units.
E_area <- function(W_vol, u) {
  .validate_scalar(W_vol, "W_vol", lower = 0)
  W_vol * (1 - stats::pnorm(u))
}

# Original expected-perimeter expression, including the window-boundary term.
# `l2` here is the spectral moment value, not the function with the same name.
E_perimeter <- function(W_bound_vol, l2, W_vol, u) {
  .validate_scalar(W_bound_vol, "W_bound_vol", lower = 0)
  .validate_scalar(l2, "l2", lower = 0)
  .validate_scalar(W_vol, "W_vol", lower = 0)
  W_bound_vol * (1 - stats::pnorm(u)) +
    sqrt(l2) / 4 * W_vol * exp(-u^2 / 2)
}

# EFSA sample-size approximation:
# p = design prevalence, cl = confidence level, N = population size,
# ms = detection sensitivity. Retains the original ceiling and does not cap
# the result at N. This function is scalar, as are the original use cases.
compute_efsa_sampsize <- function(
    p, cl = get("cl", envir = parent.frame()),
    N = get("N", envir = parent.frame()),
    ms = get("ms", envir = parent.frame())) {
  .validate_scalar(p, "p", lower = 0, upper = 1, lower_open = TRUE)
  .validate_scalar(cl, "cl", lower = 0, upper = 1,
                   lower_open = TRUE, upper_open = TRUE)
  .validate_scalar(N, "N", lower = 1, integer = TRUE)
  .validate_scalar(ms, "ms", lower = 0, upper = 1, lower_open = TRUE)
  ceiling(((1 - (1 - cl)^(1 / (N * p))) *
             (N - 0.5 * (N * p * ms - 1))) / ms)
}

# Choose a nonnegative Gaussian threshold giving a target expected perimeter.
# The returned dp is the theoretical exceedance probability at that threshold.
# On [0, 100], expected perimeter is nonincreasing; unattainable targets produce
# an informative error instead of an opaque uniroot failure.
compute_threshold_and_dp <- function(
    W_bound_vol, W_vol, rho,
    nu = get("nu", envir = parent.frame()),
    perimeter_agg = get("perimeter_agg", envir = parent.frame())) {
  .validate_scalar(perimeter_agg, "perimeter_agg", lower = 0,
                   lower_open = TRUE)
  spectral_moment <- l2(nu, rho)
  perimeter_residual <- function(u) {
    E_perimeter(W_bound_vol, spectral_moment, W_vol, u) / perimeter_agg - 1
  }
  interval <- c(0, 100)
  residuals <- vapply(interval, perimeter_residual, numeric(1))
  if (any(!is.finite(residuals)) || all(residuals > 0) || all(residuals < 0)) {
    stop("Target perimeter cannot be attained by a threshold in [0, 100].",
         call. = FALSE)
  }
  threshold <- stats::uniroot(perimeter_residual, interval = interval)$root
  list(threshold = threshold, dp = 1 - stats::pnorm(threshold))
}

# ============================================================================
# Sampling and surveillance metrics
# ============================================================================

# Balanced accuracy on unsampled individuals, with the original convention
# of assigning zero sensitivity/specificity when the respective class is absent.
.balanced_accuracy <- function(predictions, true_labels) {
  positive <- true_labels == 1
  negative <- true_labels == 0
  sensitivity <- if (any(positive)) mean(predictions[positive] == 1) else 0
  specificity <- if (any(negative)) mean(predictions[negative] == 0) else 0
  (sensitivity + specificity) / 2
}

# Shared result constructor keeps all branches in the same field order.
.surveillance_result <- function(radius, prev_estim, disease_found,
                                 balanced_accuracy) {
  list(radius = radius, prev_estim = prev_estim,
       disease_found = disease_found, balanced_accuracy = balanced_accuracy)
}

# Evaluate a fixed sample. Indices must be unique, valid population indices.
# Empty samples retain the original no-detection return values.
# For partial samples, radius measures coverage of the full positive population
# by sampled positives. It uses true population marks and is a simulation metric.
compute_min_radius_with_sample <- function(mpp, sample_indices) {
  .validate_pattern(mpp)
  .validate_indices(sample_indices, mpp$n)
  population_size <- mpp$n
  positive_indices <- which(mpp$marks == 1)
  population_prevalence <- length(positive_indices) / population_size

  if (length(sample_indices) >= population_size) {
    return(.surveillance_result(1, population_prevalence, 1, 1))
  }

  sampled_positive_indices <- sample_indices[mpp$marks[sample_indices] == 1]
  if (length(sampled_positive_indices) == 0L) {
    return(.surveillance_result(0, 0, 0, 0.5))
  }
  if (length(sampled_positive_indices) == length(positive_indices)) {
    return(.surveillance_result(1, population_prevalence, 1, 1))
  }

  # Rows correspond to all positives; columns to sampled positives.
  # The largest row minimum is the smallest radius covering every positive.
  positive_distances <- spatstat.geom::crossdist(
    mpp$x[positive_indices], mpp$y[positive_indices],
    mpp$x[sampled_positive_indices], mpp$y[sampled_positive_indices]
  )
  min_radius <- max(apply(positive_distances, 1, min))

  # Nearest-neighbour classification is evaluated only on unsampled points.
  # which.min resolves distance ties in sample order, as in the original.
  unsampled_indices <- setdiff(seq_len(population_size), sample_indices)
  neighbour_distances <- spatstat.geom::crossdist(
    mpp$x[unsampled_indices], mpp$y[unsampled_indices],
    mpp$x[sample_indices], mpp$y[sample_indices]
  )
  nearest_indices <- apply(neighbour_distances, 1, which.min)
  predictions <- mpp$marks[sample_indices[nearest_indices]]
  accuracy <- .balanced_accuracy(predictions, mpp$marks[unsampled_indices])

  .surveillance_result(
    radius = 1 - min_radius / sqrt(2 * spatstat.geom::area.owin(mpp$window)),
    prev_estim = length(sampled_positive_indices) / length(sample_indices),
    disease_found = TRUE,
    balanced_accuracy = accuracy
  )
}

# Simple random sampling without replacement. Set the seed in the calling
# script for reproducibility; this function does not reset the RNG.
# Values above population size still request a census, as in the original.
compute_min_radius <- function(mpp, sample_size) {
  .validate_pattern(mpp)
  .validate_scalar(sample_size, "sample_size", lower = 0, integer = TRUE)
  sample_indices <- if (sample_size >= mpp$n) {
    seq_len(mpp$n)
  } else {
    sample(seq_len(mpp$n), size = sample_size, replace = FALSE)
  }
  compute_min_radius_with_sample(mpp, sample_indices)
}

# ============================================================================
# Generalized spatial surveillance simulation
# ============================================================================
# Run with Rscript to execute the configured scenarios, draw comparisons, and
# save Generalized_Simulation_Results.RData. When sourced, only functions are
# defined; call main() to run the experiment.
#
# Place the refactored functions_Correlation_PrevEstim_Efficiency_Generalized.R
# beside this file, or pass its path to main(helper_file = ...).
# To customize a run:
#   source("Correlation_PrevEstim_Efficiency_Generalized.R")
#   config <- simulation_config()
#   config$N_sim_test <- 100
#   config$seed <- 123
#   results <- main(config, make_plots = FALSE)
#
# Required packages: GpGp, spatstat.geom, spatstat.random, ggplot2, gridExtra.
# No packages are attached to the search path or installed automatically.
#
# Statistical conventions:
# * "poisson" denotes N uniform points, i.e. a homogeneous Poisson process
#   conditional on its count. Thomas and Strauss counts remain random.
# * EFSA budgets use nominal N, even when the realized host count differs.
# * True infection marks are obtained by empirical quantile thresholding.
# * Gaussian fields are reused across prevalences and within replicate blocks.
#   Replicates sharing a field are not independent epidemic realizations.
# * Noisy rho is drawn once per prevalence/rho combination, not per replicate.
# * Prevalence estimates remain unweighted, including under adaptive designs.
# * Sampling weights in sample(..., prob = ...) are sequential selection weights,
#   not the final first-order inclusion probabilities of the selected sample.
# * Metric conventions (including special-case scores) reside in the helper.
#
# Two implementation corrections: field blocks now start at replicates 1,
# delay_gpsim + 1, ...; and sampling a singleton index uses sample.int to avoid
# R's sample(x) interpretation of a single numeric x as the sequence 1:x.

# Resolve the script directory before the sourcing/execution frame disappears.
.simulation_script_dir <- local({
  source_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (is.null(source_file)) {
    file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (length(file_arg)) source_file <- sub("^--file=", "", file_arg[1])
  }
  if (is.null(source_file)) getwd() else dirname(normalizePath(source_file))
})

check_dependencies <- function() {
  packages <- c("GpGp", "spatstat.geom", "spatstat.random", "ggplot2", "gridExtra")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Missing required packages: ", paste(missing, collapse = ", "),
         ". Install them before running the simulation.", call. = FALSE)
  }
  invisible(TRUE)
}

load_simulation_helpers <- function(helper_file) {
  if (!file.exists(helper_file)) {
    stop("Helper file not found: ", helper_file, call. = FALSE)
  }
  target <- environment(load_simulation_helpers)
  sys.source(helper_file, envir = target)
  required <- c("l2", "E_perimeter", "compute_threshold_and_dp",
                "compute_efsa_sampsize", "compute_min_radius_with_sample")
  available <- vapply(required, exists, logical(1), envir = target,
                      mode = "function", inherits = TRUE)
  if (!all(available)) stop("The helper file is missing required functions.")
  # The companion refactor accepts these parameters explicitly.
  if (!all(c("nu", "perimeter_agg") %in%
           names(formals(get("compute_threshold_and_dp", envir = target))))) {
    stop("Use the refactored companion helper file with explicit parameters.",
         call. = FALSE)
  }
  invisible(TRUE)
}

# ============================================================================
# Configuration and design-prevalence tables
# ============================================================================

simulation_config <- function() {
  list(
    cl = 0.9, N = 10000L, ms = 0.8,
    prev = c(0.0005, 0.005, 0.01, 0.05),
    dp_base = 0.01,
    rho_test = c(seq(1.5, 9.5, length.out = 10),
                 seq(10, 100, length.out = 20)),
    nu = 1.5, rho_agg = 100, dp_agg = 0.01, W_size = 10000,
    # Two replicates retain the original smoke-run default; increase for analysis.
    N_sim_test = 2L, delay_gpsim = 50L, strata_grid = 4L, tau_noise = 0.4,
    risk_weight = 0.9, stage1_fraction = 0.6,
    seed = NULL,                 # Set an integer for a reproducible experiment.
    output_dir = ".", verbose = TRUE
  )
}

check_number <- function(x, name, lower = 0, upper = Inf,
                         integer = FALSE, positive = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x < lower || x > upper || (positive && x <= 0) ||
      (integer && x != floor(x))) {
    stop("Invalid value for ", name, ".", call. = FALSE)
  }
  invisible(x)
}

validate_config <- function(config) {
  check_number(config$N, "N", lower = 1, integer = TRUE)
  check_number(config$W_size, "W_size", positive = TRUE)
  check_number(config$nu, "nu", lower = 1)
  if (config$nu <= 1) stop("nu must exceed 1 for the spectral moment.")
  check_number(config$cl, "cl", positive = TRUE, upper = 1)
  if (config$cl >= 1) stop("cl must be less than 1.")
  check_number(config$ms, "ms", positive = TRUE, upper = 1)
  for (name in c("rho_agg", "dp_agg", "dp_base")) {
    check_number(config[[name]], name, positive = TRUE,
                 upper = if (name == "rho_agg") Inf else 1)
  }
  for (name in c("prev", "rho_test")) {
    values <- config[[name]]
    if (!is.numeric(values) || !length(values) || any(!is.finite(values)) ||
        any(values <= 0) || (name == "prev" && any(values >= 1))) {
      stop(name, " must contain valid positive values.", call. = FALSE)
    }
  }
  check_number(config$N_sim_test, "N_sim_test", lower = 1, integer = TRUE)
  check_number(config$delay_gpsim, "delay_gpsim", lower = 1, integer = TRUE)
  check_number(config$strata_grid, "strata_grid", lower = 1, integer = TRUE)
  check_number(config$tau_noise, "tau_noise")
  check_number(config$risk_weight, "risk_weight", upper = 1)
  check_number(config$stage1_fraction, "stage1_fraction", upper = 1)
  invisible(config)
}

perimeter_parameters <- function(config) {
  boundary <- 4 * config$W_size
  area <- config$W_size^2
  target <- E_perimeter(boundary, l2(config$nu, config$rho_agg), area,
                        -stats::qnorm(config$dp_agg))
  list(W_bound_vol = boundary, W_vol = area, nu = config$nu,
       perimeter_agg = target)
}

threshold_for_range <- function(rho, config) {
  do.call(compute_threshold_and_dp,
          c(perimeter_parameters(config), list(rho = rho)))
}

# Construct numeric columns directly; avoid transposing lists into list columns.
make_design_table <- function(config, fixed_dp = NULL) {
  values <- lapply(config$rho_test, threshold_for_range, config = config)
  table <- data.frame(
    threshold = vapply(values, `[[`, numeric(1), "threshold"),
    dp = vapply(values, `[[`, numeric(1), "dp"),
    rho = config$rho_test
  )
  # For fixed-dp designs, threshold remains the corrected theoretical threshold
  # for reference; only dp determines the budget, as in the original script.
  if (!is.null(fixed_dp)) table$dp <- fixed_dp
  table$n <- vapply(table$dp, compute_efsa_sampsize, numeric(1),
                    cl = config$cl, N = config$N, ms = config$ms)
  table
}

# ============================================================================
# Host patterns and sampling designs
# ============================================================================

generate_host_locations <- function(pattern_type, N, W_size) {
  pattern_type <- match.arg(pattern_type, c("regular", "poisson", "clustered", "repulsive"))
  check_number(N, "N", lower = 1, integer = TRUE)
  check_number(W_size, "W_size", positive = TRUE)
  window <- spatstat.geom::owin(c(0, W_size), c(0, W_size))
  if (pattern_type == "regular") {
    # Retain the original grid endpoints (1 and W_size) and rounded grid size.
    side <- seq(to = W_size, length.out = round(sqrt(N)))
    return(as.matrix(expand.grid(side, side)))
  }
  pattern <- switch(
    pattern_type,
    poisson = spatstat.random::runifpoint(N, win = window),
    clustered = spatstat.random::rThomas(
      kappa = N / (W_size^2 * 10), scale = W_size / 20, mu = 10, win = window
    ),
    # gamma = 0.1 is soft inhibition, rather than a strict hard-core process.
    repulsive = spatstat.random::rStrauss(
      beta = 2 * N / W_size^2, gamma = 0.1,
      R = 0.8 * W_size / sqrt(N), W = window
    )
  )
  cbind(pattern$x, pattern$y)
}

# Subset by sampled positions to handle length-one numeric vectors safely.
sample_indices_from <- function(indices, size, prob = NULL) {
  if (size == 0L) return(integer(0))
  indices[sample.int(length(indices), size = size, replace = FALSE, prob = prob)]
}

stratified_sample <- function(mpp, sample_size, strata_grid) {
  x_stratum <- cut(mpp$x, seq(mpp$window$xrange[1], mpp$window$xrange[2],
                            length.out = strata_grid + 1L),
                   labels = FALSE, include.lowest = TRUE)
  y_stratum <- cut(mpp$y, seq(mpp$window$yrange[1], mpp$window$yrange[2],
                            length.out = strata_grid + 1L),
                   labels = FALSE, include.lowest = TRUE)
  strata <- (x_stratum - 1L) * strata_grid + y_stratum
  n_strata <- strata_grid^2
  allocation <- rep(floor(sample_size / n_strata), n_strata)
  remainder <- sample_size %% n_strata
  if (remainder > 0L) allocation[seq_len(remainder)] <- allocation[seq_len(remainder)] + 1L
  # Equal allocation across spatial strata (not proportional allocation).
  selected <- unlist(lapply(seq_len(n_strata), function(s) {
    indices <- which(strata == s)
    sample_indices_from(indices, min(allocation[s], length(indices)))
  }), use.names = FALSE)
  # Redistribute shortfalls from empty/small strata by SRS over remaining hosts.
  remaining <- setdiff(seq_len(mpp$n), selected)
  c(selected, sample_indices_from(remaining, sample_size - length(selected)))
}

risk_sample <- function(mpp, sample_size, risk_surface, risk_weight) {
  if (is.null(risk_surface)) {
    distances <- sqrt((mpp$x - mean(mpp$window$xrange))^2 +
                      (mpp$y - mean(mpp$window$yrange))^2)
    scores <- if (max(distances) == 0) rep(1, mpp$n) else
      1 - (distances / max(distances))^2
  } else {
    scores <- risk_surface(mpp$x, mpp$y)
  }
  if (!is.numeric(scores) || length(scores) != mpp$n ||
      any(!is.finite(scores)) || any(scores < 0)) {
    stop("The risk surface must return one finite nonnegative score per host.")
  }
  # A flat zero surface has no ranking information; use uniform weights.
  weights <- if (sum(scores) == 0) rep(1 / mpp$n, mpp$n) else
    risk_weight * scores / sum(scores) + (1 - risk_weight) / mpp$n
  if (sum(weights > 0) < sample_size) {
    stop("Too few hosts have positive risk weights for the requested sample.")
  }
  sample_indices_from(seq_len(mpp$n), sample_size, weights)
}

two_stage_sample <- function(mpp, sample_size, stage1_fraction) {
  stage1 <- sample_indices_from(seq_len(mpp$n), floor(sample_size * stage1_fraction))
  positives <- stage1[mpp$marks[stage1] == 1]
  remaining <- setdiff(seq_len(mpp$n), stage1)
  weights <- NULL
  if (length(positives) && length(remaining)) {
    distances <- spatstat.geom::crossdist(
      mpp$x[remaining], mpp$y[remaining], mpp$x[positives], mpp$y[positives]
    )
    min_distance <- apply(distances, 1, min)
    # Preserve the original inverse-distance preference. If every remaining
    # point coincides with a positive, all are equally close: sample uniformly.
    if (max(min_distance) > 0) {
      weights <- 1 / (min_distance + 0.1 * max(min_distance))
    }
  }
  stage2 <- sample_indices_from(remaining, sample_size - length(stage1), weights)
  c(stage1, stage2)
}

sample_points <- function(mpp, sample_size, design_type = "srs",
                          risk_surface = NULL, strata_grid = 4,
                          detected_positives = NULL, risk_weight = 0.9,
                          stage1_fraction = 0.6) {
  design_type <- match.arg(design_type, c("srs", "stratified", "risk_based", "two_stage"))
  check_number(sample_size, "sample_size", integer = TRUE)
  check_number(strata_grid, "strata_grid", lower = 1, integer = TRUE)
  check_number(risk_weight, "risk_weight", upper = 1)
  check_number(stage1_fraction, "stage1_fraction", upper = 1)
  if (!inherits(mpp, "ppp") || mpp$n < 1L) stop("mpp must be a nonempty ppp.")
  # The original accepted this argument but never used it. Reject it explicitly
  # rather than implying that prior detections influence the implemented design.
  if (!is.null(detected_positives)) stop("detected_positives is not implemented.")
  if (sample_size >= mpp$n) return(seq_len(mpp$n))
  if (sample_size == 0L) return(integer(0))
  switch(design_type,
         srs = sample_indices_from(seq_len(mpp$n), sample_size),
         stratified = stratified_sample(mpp, sample_size, strata_grid),
         risk_based = risk_sample(mpp, sample_size, risk_surface, risk_weight),
         two_stage = two_stage_sample(mpp, sample_size, stage1_fraction))
}

# ============================================================================
# Range misspecification and simulation
# ============================================================================

compute_dp_with_error <- function(rho_true, error_type, tau = 0.4,
                                  config = simulation_config()) {
  error_type <- match.arg(error_type, c("oracle", "minus_50", "minus_20",
                                      "plus_20", "plus_50", "noisy"))
  check_number(rho_true, "rho_true", positive = TRUE)
  check_number(tau, "tau")
  rho_used <- switch(error_type, oracle = rho_true,
                     minus_50 = rho_true * 0.5, minus_20 = rho_true * 0.8,
                     plus_20 = rho_true * 1.2, plus_50 = rho_true * 1.5,
                     noisy = exp(stats::rnorm(1, log(rho_true), tau)))
  list(rho_used = rho_used, dp = threshold_for_range(rho_used, config)$dp)
}

compute_cor_prev_estim_efficiency_generalized <- function(
    dp_df, N_sim = 100, host_pattern = "regular", sampling_design = "srs",
    rho_error_scenario = "oracle", risk_surface_fn = NULL, strata_grid = 4,
    tau_noise = 0.4, delay_gpsim = 50, config = simulation_config()) {
  validate_config(config)
  check_number(N_sim, "N_sim", lower = 1, integer = TRUE)
  check_number(delay_gpsim, "delay_gpsim", lower = 1, integer = TRUE)
  rho_error_scenario <- match.arg(rho_error_scenario,
    c("oracle", "minus_50", "minus_20", "plus_20", "plus_50", "noisy"))
  if (!is.data.frame(dp_df) || !all(c("rho", "dp") %in% names(dp_df)) ||
      !nrow(dp_df) || !is.numeric(dp_df$rho) || !is.numeric(dp_df$dp) ||
      any(!is.finite(dp_df$rho)) || any(dp_df$rho <= 0) ||
      any(!is.finite(dp_df$dp)) || any(dp_df$dp <= 0 | dp_df$dp > 1)) {
    stop("dp_df must have numeric positive rho and dp columns (dp <= 1).")
  }
  log_progress <- function(...) if (isTRUE(config$verbose)) message(...)
  log_progress("Running: ", host_pattern, " / ", sampling_design, " / ", rho_error_scenario)
  locations <- generate_host_locations(host_pattern, config$N, config$W_size)
  if (!nrow(locations)) stop("The generated host pattern is empty.")
  log_progress("Generated ", nrow(locations), " hosts (nominal N = ", config$N, ").")
  window <- spatstat.geom::owin(c(0, config$W_size), c(0, config$W_size))

  # Pre-generate every field before sampling, retaining shared fields across
  # prevalences. Compact block storage replaces the original sparse lists.
  # Using (sim - 1) fixes the old boundary shift and supports delay_gpsim = 1.
  block_index <- (seq_len(N_sim) - 1L) %/% delay_gpsim + 1L
  n_blocks <- max(block_index)
  fields <- lapply(dp_df$rho, function(rho) {
    lapply(seq_len(n_blocks), function(block) {
      GpGp::fast_Gp_sim(c(1, rho, config$nu, 0), locs = locations)
    })
  })

  results <- vector("list", length(config$prev))
  names(results) <- paste0("prev_", gsub("\\.", "_", as.character(config$prev)))
  for (p_idx in seq_along(config$prev)) {
    prevalence <- config$prev[p_idx]
    log_progress("Prevalence: ", prevalence)
    rows <- vector("list", nrow(dp_df))
    for (i in seq_len(nrow(dp_df))) {
      rho <- dp_df$rho[i]
      design <- if (rho_error_scenario == "oracle") {
        list(rho_used = rho, dp = dp_df$dp[i])
      } else {
        compute_dp_with_error(rho, rho_error_scenario, tau_noise, config)
      }
      budget <- compute_efsa_sampsize(design$dp, cl = config$cl,
                                      N = config$N, ms = config$ms)
      metrics <- matrix(NA_real_, nrow = N_sim, ncol = 3L,
                        dimnames = list(NULL, c("efficiency", "prev_estim", "balanced_accuracy")))
      for (sim in seq_len(N_sim)) {
        field <- fields[[i]][[block_index[sim]]]
        threshold <- stats::quantile(field, probs = 1 - prevalence)
        marks <- as.integer(field > threshold)
        pattern <- spatstat.geom::ppp(locations[, 1], locations[, 2],
                                       window = window, marks = marks)
        indices <- sample_points(pattern, budget, sampling_design,
                                  risk_surface = risk_surface_fn, strata_grid = strata_grid,
                                  risk_weight = config$risk_weight,
                                  stage1_fraction = config$stage1_fraction)
        score <- compute_min_radius_with_sample(pattern, indices)
        metrics[sim, ] <- c(score$radius, score$prev_estim, score$balanced_accuracy)
      }
      rows[[i]] <- data.frame(
        rho_true = rho, rho_used = design$rho_used,
        efficiency = mean(metrics[, "efficiency"]),
        prev_estim = mean(metrics[, "prev_estim"]),
        prev_estim_sd = stats::sd(metrics[, "prev_estim"]),
        balanced_accuracy = mean(metrics[, "balanced_accuracy"]),
        true_prevalence = prevalence, host_pattern = host_pattern,
        sampling_design = sampling_design, rho_error = rho_error_scenario
      )
    }
    results[[p_idx]] <- do.call(rbind, rows)
  }
  results
}

# ============================================================================
# Scenario orchestration, comparisons and saving
# ============================================================================

make_risk_surface <- function(W_size) {
  function(x, y) {
    distances <- sqrt((x - W_size / 2)^2 + (y - W_size / 2)^2)
    pmax(1 - (distances / (sqrt(2) * W_size / 2))^2, 0.1)
  }
}

# A table replaces repeated near-identical calls; result object names are kept
# so existing downstream analyses can still load the same RData objects.
scenario_table <- function() {
  data.frame(
    result_name = c("results_baseline", "results_poisson", "results_clustered",
                    "results_stratified", "results_risk_based", "results_two_stage",
                    "results_rho_minus50", "results_rho_plus50", "results_efsa"),
    host = c("regular", "poisson", "clustered", rep("regular", 6)),
    design = c(rep("srs", 3), "stratified", "risk_based", "two_stage", rep("srs", 3)),
    rho_error = c(rep("oracle", 6), "minus_50", "plus_50", "oracle"),
    fixed_dp = c(rep(FALSE, 8), TRUE)
  )
}

run_scenarios <- function(config) {
  corrected <- make_design_table(config)
  fixed <- make_design_table(config, fixed_dp = config$dp_base)
  scenarios <- scenario_table()
  results <- setNames(vector("list", nrow(scenarios)), scenarios$result_name)
  for (i in seq_len(nrow(scenarios))) {
    results[[i]] <- compute_cor_prev_estim_efficiency_generalized(
      dp_df = if (scenarios$fixed_dp[i]) fixed else corrected,
      N_sim = config$N_sim_test, host_pattern = scenarios$host[i],
      sampling_design = scenarios$design[i], rho_error_scenario = scenarios$rho_error[i],
      risk_surface_fn = if (scenarios$design[i] == "risk_based")
        make_risk_surface(config$W_size) else NULL,
      strata_grid = config$strata_grid, tau_noise = config$tau_noise,
      delay_gpsim = config$delay_gpsim, config = config
    )
  }
  results
}

# Accepts named scenario data frames; each may contain one or all prevalences.
plot_comparison <- function(results_list, title_suffix = "") {
  plot_data <- do.call(rbind, lapply(names(results_list), function(name) {
    result <- results_list[[name]]
    result$scenario <- name
    result
  }))
  make_plot <- function(column, title, y_label, show_truth = FALSE) {
    plot <- ggplot2::ggplot(
      plot_data, ggplot2::aes(x = rho_true, y = .data[[column]],
                              color = scenario, linetype = scenario)
    ) +
      ggplot2::geom_line(linewidth = 1) + ggplot2::geom_point(size = 2) +
      ggplot2::facet_wrap(~ true_prevalence, scales = "free_y",
                          labeller = ggplot2::label_both) +
      ggplot2::labs(title = paste(title, title_suffix),
                    x = "True Rho (Spatial Range)", y = y_label,
                    color = "Scenario", linetype = "Scenario") +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(legend.position = "bottom") +
      ggplot2::scale_color_brewer(palette = "Set1")
    if (show_truth) {
      plot <- plot + ggplot2::geom_hline(
        data = unique(plot_data["true_prevalence"]),
        ggplot2::aes(yintercept = true_prevalence),
        inherit.aes = FALSE, linetype = "dashed", color = "black"
      )
    }
    plot
  }
  plots <- list(
    make_plot("balanced_accuracy", "Balanced Accuracy vs Rho", "Balanced Accuracy"),
    make_plot("prev_estim", "Prevalence Estimation vs Rho", "Estimated Prevalence", TRUE),
    make_plot("efficiency", "Detection Efficiency vs Rho", "Radius Efficiency")
  )
  invisible(gridExtra::grid.arrange(grobs = plots, ncol = 1))
}

plot_all_comparisons <- function(results) {
  combine <- function(names) do.call(rbind, results[[names]])
  # Facet over every configured prevalence instead of repeating the same plot
  # for each prevalence and then plotting the first prevalence a second time.
  plot_comparison(list(regular = combine("results_baseline"),
                       poisson = combine("results_poisson"),
                       clustered = combine("results_clustered")), "- Host Patterns")
  plot_comparison(list(SRS = combine("results_baseline"),
                       stratified = combine("results_stratified"),
                       risk_based = combine("results_risk_based"),
                       two_stage = combine("results_two_stage")), "- Sampling Designs")
  plot_comparison(list(oracle = combine("results_baseline"),
                       "rho -50%" = combine("results_rho_minus50"),
                       "rho +50%" = combine("results_rho_plus50"),
                       EFSA_fixed = combine("results_efsa")), "- Rho Misspecification")
  invisible(NULL)
}

main <- function(config = simulation_config(),
                 helper_file = file.path(.simulation_script_dir,
                   "functions_Correlation_PrevEstim_Efficiency_Generalized.R"),
                 make_plots = TRUE) {
  validate_config(config)
  check_dependencies()
  load_simulation_helpers(helper_file)
  if (!is.null(config$seed)) set.seed(config$seed)
  results <- run_scenarios(config)
  if (!dir.exists(config$output_dir) &&
      !dir.create(config$output_dir, recursive = TRUE)) {
    stop("Could not create output directory: ", config$output_dir)
  }
  output_file <- file.path(config$output_dir, "Generalized_Simulation_Results.RData")
  # Preserve the original nine named objects and also record the configuration.
  output_env <- list2env(c(results, list(simulation_settings = config)), parent = emptyenv())
  save(list = ls(output_env), envir = output_env, file = output_file)
  # Save before plotting so a graphics failure does not discard computed results.
  message("Results saved to: ", normalizePath(output_file))
  if (make_plots) plot_all_comparisons(results)
  invisible(results)
}

# Rscript executes the experiment; source() leaves execution to the caller.
if (sys.nframe() == 0L) main()

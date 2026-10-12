# -----------------------------------------------------------------------------
# Gaussian random field: prevalence scenarios and point-pattern illustrations
# -----------------------------------------------------------------------------
# Simulate one Matern Gaussian field on a regular grid and threshold the same
# realization at two prevalence levels:
#   1. A fixed target prevalence.
#   2. A prevalence corrected using compute_threshold_and_dp().
# Both binary fields are exported as high-resolution PNG images.
#
# Requirements:
#   - R package GpGp
#   - functions_Correlation_PrevEstim_Efficiency.R in the working directory
#   - W_bound_vol and W_vol, defined by the user or the sourced script.
#
# Run from the directory containing the helper script, or update helper_file.
# -----------------------------------------------------------------------------

# ---- Configuration -----------------------------------------------------------

n_points <- 90000L                # Requested number of grid points (perfect square)
window_size <- 10000              # Side length of the square observation window
rho <- 10                         # Matern correlation/range parameter
nu <- 1.5                         # Matern smoothness parameter
fixed_prevalence <- 0.01          # Proportion of points classified as positive
point_cex <- 1                    # Point size in the plots
image_size <- 3000L               # PNG width and height in pixels
image_resolution <- 300L          # PNG resolution in dots per inch
output_dir <- "."                 # Folder for exported images
helper_file <- "functions_Correlation_PrevEstim_Efficiency.R"
# set.seed(123)                    # Uncomment to reproduce the random field

# ---- Dependencies ------------------------------------------------------------

if (!requireNamespace("GpGp", quietly = TRUE)) {
  stop("Please install the 'GpGp' package: install.packages('GpGp').")
}
if (!file.exists(helper_file)) {
  stop("Missing helper script: ", helper_file,
       ". Adjust 'helper_file' or the working directory.")
}
source(helper_file)

# ---- Utility functions -------------------------------------------------------

# Generate a square regular lattice with exactly n_points locations.
make_regular_locations <- function(n_points, window_size) {
  side_points <- sqrt(n_points)
  if (!is.finite(side_points) || side_points < 2 ||
      side_points != floor(side_points)) {
    stop("'n_points' must be a perfect square of at least 4.")
  }
  if (!is.finite(window_size) || window_size <= 0) {
    stop("'window_size' must be strictly positive.")
  }

  axis <- seq(from = 1, to = window_size, length.out = as.integer(side_points))
  as.matrix(expand.grid(x = axis, y = axis))
}

# Map binary infection marks to consistent, color-blind-friendly colors.
get_mark_colors <- function(marks) {
  colors <- rep(NA_character_, length(marks))
  colors[!is.na(marks) & marks == 1] <- "#E69F00"  # Positive
  colors[!is.na(marks) & marks == 0] <- "#0072B2"  # Negative
  colors
}

# Classify the highest values of a Gaussian field as positive.
# Strict '>' matches the original script. Ties may yield fewer positives
# than the desired proportion (normally rare for a continuous Gaussian field).
threshold_field <- function(field, prevalence) {
  if (length(prevalence) != 1L || is.na(prevalence) ||
      !is.finite(prevalence) || prevalence < 0 || prevalence > 1) {
    stop("'prevalence' must be a single number in [0, 1].")
  }
  if (!is.numeric(field) || anyNA(field) || any(!is.finite(field))) {
    stop("'field' must contain finite numeric values without missing entries.")
  }
  cutoff <- stats::quantile(field, probs = 1 - prevalence, names = FALSE)
  as.integer(field > cutoff)
}

# Draw a point-based infection map; optionally overlay points on an existing plot.
plot_point_pattern <- function(locs, marks, cex = 1, main = NULL,
                               add = FALSE) {
  if (nrow(locs) != length(marks) || ncol(locs) != 2L) {
    stop("'locs' must have two columns and one row per mark.")
  }
  colors <- get_mark_colors(marks)

  if (add) {
    graphics::points(locs, col = colors, pch = 15, cex = cex)
  } else {
    old_par <- graphics::par(pty = "s")
    on.exit(graphics::par(old_par), add = TRUE)
    graphics::plot(locs, col = colors, pch = 15, cex = cex,
                   main = main, xlab = "", ylab = "",
                   xaxt = "n", yaxt = "n", ann = FALSE)
  }
  invisible(NULL)
}

# Open and reliably close an image device, even if plotting fails.
save_point_pattern <- function(filename, locs, marks, cex = 1,
                               size = 3000L, resolution = 300L) {
  grDevices::png(filename = filename, width = size, height = size,
                 res = resolution)
  on.exit(grDevices::dev.off(), add = TRUE)
  plot_point_pattern(locs, marks, cex = cex)
  invisible(filename)
}

# ---- Simulate one Gaussian field ---------------------------------------------

locations <- make_regular_locations(n_points, window_size)

# GpGp's Matern covariance parameterization: variance, range, smoothness, nugget.
gaussian_field <- GpGp::fast_Gp_sim(
  c(1, rho, nu, 0), locs = locations
)

# ---- Scenario 1: fixed prevalence --------------------------------------------

marks_fixed <- threshold_field(gaussian_field, fixed_prevalence)
plot_point_pattern(locations, marks_fixed, cex = point_cex)

# Preserve the original output names for the default configuration.
fixed_png <- file.path(output_dir, "rho10_0dot01.png")
save_point_pattern(fixed_png, locations, marks_fixed, cex = point_cex,
                   size = image_size, resolution = image_resolution)

# ---- Scenario 2: prevalence corrected for correlation range ------------------

# These quantities are inputs to the external function; they are not defined
# by the original example. Supply them in the helper script or above this block.
required_inputs <- c("W_bound_vol", "W_vol")
missing_inputs <- required_inputs[!vapply(
  required_inputs, exists, logical(1), inherits = TRUE
)]
if (length(missing_inputs) > 0L) {
  stop("Corrected scenario requires the following objects: ",
       paste(missing_inputs, collapse = ", "),
       ". Define them before computing corrected prevalence.")
}
if (!exists("compute_threshold_and_dp", mode = "function")) {
  stop("Helper script must define 'compute_threshold_and_dp()'.")
}

corrected_prevalence <- compute_threshold_and_dp(W_bound_vol, W_vol, rho)$dp
marks_corrected <- threshold_field(gaussian_field, corrected_prevalence)
plot_point_pattern(locations, marks_corrected, cex = point_cex)

corrected_png <- file.path(output_dir, "rho10_0dot01_correc.png")
save_point_pattern(corrected_png, locations, marks_corrected,
                   cex = point_cex, size = image_size,
                   resolution = image_resolution)

message("Fixed prevalence: ", fixed_prevalence,
        "; corrected prevalence: ", signif(corrected_prevalence, 4))
message("Images saved to: ", normalizePath(output_dir))

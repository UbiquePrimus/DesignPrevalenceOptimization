# Second spectral moment of Matern covariance function
l2 <- function(nu,rho){
  return(1/(2*rho^2*(nu-1)))
}

# Expected area of excursion set above u
E_area <- function(W_vol,u){
  return(W_vol*(1-pnorm(u)))
}

# Expected perimeter of excursion set above u
E_perimeter <- function(W_bound_vol, l2, W_vol, u){
  return(W_bound_vol*(1-pnorm(u)) + sqrt(l2)/4*W_vol*exp(-u^2/2) )
}

compute_efsa_sampsize <- function(p){
  return(ceiling(((1-(1-cl)^(1/(N*p))) * (N-(0.5*(N*p*ms -1))))/ms))
}

# Function to compute threshold and dp for other cases
compute_threshold_and_dp <- function(W_bound_vol, W_vol, rho){
  l2_rho <- l2(nu,rho)
  perimeter_u <- function(u){return(E_perimeter(W_bound_vol, l2_rho, W_vol, u)/perimeter_agg-1)}
  u_rho <- uniroot(perimeter_u, c(0,100))$root
  
  return(list(threshold = u_rho,
              dp = 1-pnorm(u_rho)))
}

# compute_min_radius <- function(mpp, sample_size) {
#   # Uniformly sample points
#   if(sample_size >= mpp$n) return(list(
#     radius = 1,
#     prev_estim = sum(mpp$marks==1)/mpp$n,
#     disease_found = 1
#   ))
#   sample_indices <- sample(1:mpp$n, sample_size)
#   sampled_points <- ppp(
#     x = mpp$x[sample_indices],
#     y = mpp$y[sample_indices],
#     window = mpp$window,
#     marks = mpp$marks[sample_indices]
#   )
#   
#   # Get positive points from sample and full population
#   sample_pos <- sampled_points[sampled_points$marks == 1]
#   all_pos <- mpp[mpp$marks == 1]
#   if(sample_pos$n == 0) return(list(
#     radius = 0,
#     prev_estim = 0,
#     disease_found = 0
#   ))
#   if(sample_pos$n == all_pos$n) return(list(
#     radius = 1,
#     prev_estim = sum(mpp$marks==1)/mpp$n,
#     disease_found = 1
#   ))
#   
#   # Create distance matrix between all positive points and sampled positive points
#   dmat <- crossdist(all_pos$x, all_pos$y, sample_pos$x, sample_pos$y)
#   # Take the maximum of minimum distances between a non-sampled positive and a sampled positive
#   min_radius <- max(apply(dmat, 1, min))
#   
#   return(list(
#     radius = 1-min_radius/sqrt(area.owin(mpp$window)*2),
#     prev_estim = sample_pos$n/sampled_points$n,
#     disease_found = sum(sampled_points$marks)>0
#   ))
# }

compute_min_radius <- function(mpp, sample_size) {
  # Uniformly sample points
  if(sample_size >= mpp$n) return(list(
    radius = 1,
    prev_estim = sum(mpp$marks==1)/mpp$n,
    disease_found = 1,
    balanced_accuracy = 1
  ))
  
  sample_indices <- sample(1:mpp$n, sample_size)
  sampled_points <- ppp(
    x = mpp$x[sample_indices],
    y = mpp$y[sample_indices],
    window = mpp$window,
    marks = mpp$marks[sample_indices]
  )
  
  # Get positive points from sample and full population
  sample_pos <- sampled_points[sampled_points$marks == 1]
  all_pos <- mpp[mpp$marks == 1]
  
  if(sample_pos$n == 0) return(list(
    radius = 0,
    prev_estim = 0,
    disease_found = 0,
    balanced_accuracy = 0.5
  ))
  
  if(sample_pos$n == all_pos$n) return(list(
    radius = 1,
    prev_estim = sum(mpp$marks==1)/mpp$n,
    disease_found = 1,
    balanced_accuracy = 1
  ))
  
  # Create distance matrix between all positive points and sampled positive points
  dmat <- crossdist(all_pos$x, all_pos$y, sample_pos$x, sample_pos$y)
  # Take the maximum of minimum distances between a non-sampled positive and a sampled positive
  min_radius <- max(apply(dmat, 1, min))
  
  # Nearest neighbor prediction for non-sampled points
  non_sample_indices <- setdiff(1:mpp$n, sample_indices)
  
  if(length(non_sample_indices) > 0) {
    # Create distance matrix between non-sampled points and all sampled points
    nn_dmat <- crossdist(mpp$x[non_sample_indices], mpp$y[non_sample_indices], 
                         sampled_points$x, sampled_points$y)
    
    # Find nearest neighbor for each non-sampled point
    nn_indices <- apply(nn_dmat, 1, which.min)
    
    # Predict based on nearest neighbor's mark
    predictions <- sampled_points$marks[nn_indices]
    
    # Get true labels for non-sampled points
    true_labels <- mpp$marks[non_sample_indices]
    
    # Calculate balanced accuracy
    # Confusion matrix components
    tp <- sum(predictions == 1 & true_labels == 1)  # True positives
    tn <- sum(predictions == 0 & true_labels == 0)  # True negatives
    fp <- sum(predictions == 1 & true_labels == 0)  # False positives
    fn <- sum(predictions == 0 & true_labels == 1)  # False negatives
    
    # Sensitivity (True Positive Rate) and Specificity (True Negative Rate)
    sensitivity <- ifelse(tp + fn > 0, tp / (tp + fn), 0)
    specificity <- ifelse(tn + fp > 0, tn / (tn + fp), 0)
    
    # Balanced accuracy = (Sensitivity + Specificity) / 2
    balanced_accuracy <- (sensitivity + specificity) / 2
  } else {
    balanced_accuracy <- 1
  }
  
  return(list(
    radius = 1 - min_radius/sqrt(area.owin(mpp$window)*2),
    prev_estim = sample_pos$n/sampled_points$n,
    disease_found = sum(sampled_points$marks) > 0,
    balanced_accuracy = balanced_accuracy
  ))
}

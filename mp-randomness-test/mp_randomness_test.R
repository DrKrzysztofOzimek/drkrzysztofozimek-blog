mp_randomness_test <- function(data, B = 1000) {
  # ---------------------------------------------------------------------------
  # Marčenko–Pastur Simulation Test (Left-Tailed)
  #
  # Args:
  #   data : A numeric matrix or data frame (rows = observations, columns = variables).
  #   B    : Number of simulations (default = 1000).
  #
  # Returns:
  #   p-value: Probability of observing a proportion of eigenvalues within
  #            MP bounds as low as or lower than the observed, assuming matrix is random.
  # ---------------------------------------------------------------------------
  
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("Input 'data' must be a data frame or matrix.")
  }
  
  data <- as.matrix(data)
  T_obs <- nrow(data)
  N_vars <- ncol(data)
  Q <- T_obs / N_vars
  
  if (Q < 1) stop("Q = T/N must be ≥ 1 for MP distribution to apply.")
  
  # Marčenko–Pastur bounds
  lambda_min <- (1 - sqrt(1 / Q))^2
  lambda_max <- (1 + sqrt(1 / Q))^2
  
  # Observed matrix
  cov_obs <- cov(scale(data, center = TRUE, scale = TRUE))
  eig_obs <- eigen(cov_obs, only.values = TRUE)$values
  prop_obs <- mean(eig_obs >= lambda_min & eig_obs <= lambda_max)
  
  # Simulate null distribution
  prop_null <- numeric(B)
  for (b in 1:B) {
    Z <- matrix(rnorm(T_obs * N_vars), nrow = T_obs)
    cov_null <- cov(scale(Z, center = TRUE, scale = TRUE))
    eig_null <- eigen(cov_null, only.values = TRUE)$values
    prop_null[b] <- mean(eig_null >= lambda_min & eig_null <= lambda_max)
  }
  
  # Left-tailed p-value: how rare is prop_obs under null?
  p_value <- mean(prop_null <= prop_obs)
  p_value
}

# Function to calculate portfolio-specific risk with rho constraint
specific_risk <- function(sigma2, rho, gamma, n) {
  if (n <= 1) {
    stop("n must be greater than 1.")
  }
  
  rho_min <- -1 / (n - 1)
  if (rho < rho_min || rho > 1) {
    warning(sprintf("rho must be in [%.4f, 1]. Given rho = %.4f", rho_min, rho))
  }
  
  term1 <- (1 - rho) * (gamma^2 + ((1 - gamma)^2) / (n - 1))
  risk <- sigma2 * (term1 + rho)
  return(risk)
}

# Example usage
sigma2 <- 0.04
rho <- 0.1
gamma <- 0.2
n <- 10

specific_risk(sigma2, rho, gamma, n)
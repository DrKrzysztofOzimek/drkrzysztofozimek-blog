library(copula)
library(kde1d)

jp <- function(data, bounds) {
  # Estimate joint probability over a rectangular region using smoothed pseudo-observations
  # and a best-fit copula selected by AIC.
  #
  # Parameters:
  # - data: R DataFrame (n x d) with continuous numeric data
  # - bounds: R DataFrame (2 x d) with lower and upper bounds
  #
  # Returns:
  # - Estimated joint probability (float)

  # Helper 1: Fit a copula to smoothed pseudo-observations and return the best-fit copula object
  #
  # This function:
  # - Applies univariate kernel density estimation (via kde1d) to each margin
  # - Transforms each margin to the unit interval [0,1] using smoothed empirical CDFs
  # - Fits several copula families (Gaussian, t, Clayton, Gumbel, Frank, Joe)
  # - Selects the best copula based on AIC
  # - Returns the fitted copula object (with corrected structure for t-copula if needed)

  copulaCDF <- function(data) {
    stopifnot(ncol(data) >= 2)
    # Step 1: Helper to infer dimension from correlation vector length
    rho2dim <- function(len) {
      d <- (1 + sqrt(1 + 8 * len)) / 2
      if (d != floor(d)) stop("Invalid rho length for tCopula.")
      return(as.integer(d))
    }

    # Step 2: Compute smoothed pseudo-observations via kde1d
    kde_list <- lapply(1:ncol(data), function(i) kde1d::kde1d(data[, i]))
    u_data <- as.data.frame(
      mapply(kde1d::pkde1d, as.list(data), kde_list)
    )
    colnames(u_data) <- paste0("u", 1:ncol(u_data))

    # Step 3: Define candidate copulas
    d <- ncol(data)
    copulas <- list(
      Gaussian = copula::normalCopula(dim = d, dispstr = "un"),
      t        = copula::tCopula(dim = d, dispstr = "un"),
      Clayton  = copula::claytonCopula(dim = d),
      Gumbel   = copula::gumbelCopula(dim = d),
      Frank    = copula::frankCopula(dim = d),
      Joe      = copula::joeCopula(dim = d)
    )

    # Step 4: Fit each copula
    fits <- lapply(copulas, function(cop) {
      tryCatch(
        suppressWarnings(copula::fitCopula(cop, data = u_data, method = "ml")),
        error = function(e) NULL
      )
    })
    names(fits) <- names(copulas)

    # Step 5: Select best fit by AIC
    criteria <- lapply(names(fits), function(name) {
      fit <- fits[[name]]
      if (!is.null(fit)) {
        data.frame(
          family = name,
          logLik = logLik(fit)[1],
          AIC = AIC(fit),
          BIC = BIC(fit)
        )
      }
    })
    results <- do.call(rbind, criteria)
    best_name <- results$family[which.min(results$AIC)]
    best_fit <- fits[[best_name]]
    best_cop <- best_fit@copula

    # Step 6: Fix tCopula (if needed)
    if (inherits(best_cop, "tCopula")) {
      est_params <- coef(best_fit)
      rho <- est_params[-length(est_params)]
      df <- round(est_params[length(est_params)])
      dim_cop <- rho2dim(length(rho))
      best_cop <- copula::tCopula(param = rho, df = df, dim = dim_cop, dispstr = "un")
    }

    return(best_cop)
  }

  # Helper 2: Evaluate smoothed marginal CDFs at a given point
  #
  # This function:
  # - Takes raw multivariate data and a point in the original data space (x_point)
  # - Fits a univariate kernel density estimate (via kde1d) to each margin
  # - Evaluates the corresponding smoothed marginal CDFs at x_point
  # - Returns the vector of pseudo-observations (u1, ..., ud) in [0,1]^d

  get_pseudo_obs <- function(data, x_point) {
    # Step 1: Ensure all columns are numeric
    if (!all(sapply(data, is.numeric))) stop("All columns in data must be numeric.")
    if (length(x_point) != ncol(data)) stop("x_point must match number of columns in data.")

    # Step 2: Fit kde1d to each column
    kde_list <- lapply(1:ncol(data), function(i) kde1d::kde1d(data[, i]))

    # Step 3: Evaluate the smoothed CDFs at each x_point value
    u_vec <- mapply(kde1d::pkde1d, x_point, kde_list)

    # Step 4: Name the result
    names(u_vec) <- paste0("u", seq_along(u_vec))
    return(u_vec)
  }

  # Helper 3: Generate all corner combinations (in [0,1]^d) from smoothed bounds
  #
  # This Helper:
  # - Applies get_pseudo_obs() to each row of the bounds matrix (lower and upper)
  # - Forms a list of pseudo-observed bounds for each variable
  # - Uses expand.grid() to generate all 2^d combinations of corners in [0,1]^d
  bounds_u <- as.data.frame(apply(bounds, 2, function(i) get_pseudo_obs(data, i)))
  bounds_u_list <- mapply(function(lo, up) c(lo, up), bounds_u[, 1], bounds_u[, 2], SIMPLIFY = FALSE)
  variations_from_bounds_u <- setNames(as.data.frame(expand.grid(bounds_u_list)), paste0("u", 1:ncol(data)))

  # Helper 4: Compute inclusion-exclusion signs for each corner
  #
  # This step:
  # - Uses expand.grid to create all binary corner index combinations (1:lower, 2:upper)
  # - Assigns sign = (-1)^(number of lower bounds used)
  variations <- expand.grid(rep(list(1:2), ncol(data)))
  signs <- apply(variations, 1, function(row) {
    num_lower <- sum(row == 1)
    (-1)^num_lower
  })

  # Main Computation: Apply inclusion-exclusion formula using the best-fit copula
  best_copula <- copulaCDF(data)
  sum(apply(variations_from_bounds_u, 1, function(i) pCopula(u = i, best_copula)) * signs)
}
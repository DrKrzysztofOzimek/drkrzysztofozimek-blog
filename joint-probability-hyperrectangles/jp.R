library(copula)
library(kde1d)

jp <- function(data, bounds) {
  # Estimate joint probability over a rectangular region using smoothed pseudo-observations
  # and a best-fit copula selected by AIC.
  #
  # Parameters:
  # - data: R DataFrame (n x d) with continuous numeric data
  # - bounds: R DataFrame (d x 2) with lower and upper bounds
  #            row i = interval for variable i: [l_i, u_i]
  #            col 1 = lower bounds (l), col 2 = upper bounds (u)
  #
  # Returns:
  # - Estimated joint probability (float)
  
  # Helper 1: Fit a copula to smoothed pseudo-observations and return the best-fit copula object
  copulaCDF <- function(data) {
    stopifnot(ncol(data) >= 2)
    
    rho2dim <- function(len) {
      d <- (1 + sqrt(1 + 8 * len)) / 2
      if (d != floor(d)) stop("Invalid rho length for tCopula.")
      return(as.integer(d))
    }
    
    kde_list <- lapply(1:ncol(data), function(i) kde1d::kde1d(data[, i]))
    u_data <- as.data.frame(
      mapply(kde1d::pkde1d, as.list(data), kde_list)
    )
    colnames(u_data) <- paste0("u", 1:ncol(u_data))
    
    d <- ncol(data)
    copulas <- list(
      Gaussian = copula::normalCopula(dim = d, dispstr = "un"),
      t        = copula::tCopula(dim = d, dispstr = "un"),
      Clayton  = copula::claytonCopula(dim = d),
      Gumbel   = copula::gumbelCopula(dim = d),
      Frank    = copula::frankCopula(dim = d),
      Joe      = copula::joeCopula(dim = d)
    )
    
    fits <- lapply(copulas, function(cop) {
      tryCatch(
        suppressWarnings(copula::fitCopula(cop, data = u_data, method = "ml")),
        error = function(e) NULL
      )
    })
    names(fits) <- names(copulas)
    
    criteria <- lapply(names(fits), function(name) {
      fit <- fits[[name]]
      if (!is.null(fit)) {
        data.frame(
          family = name,
          logLik = logLik(fit)[1],
          AIC    = AIC(fit),
          BIC    = BIC(fit)
        )
      }
    })
    results  <- do.call(rbind, criteria)
    best_name <- results$family[which.min(results$AIC)]
    best_fit  <- fits[[best_name]]
    best_cop  <- best_fit@copula
    
    if (inherits(best_cop, "tCopula")) {
      est_params <- coef(best_fit)
      rho        <- est_params[-length(est_params)]
      df         <- round(est_params[length(est_params)])
      dim_cop    <- rho2dim(length(rho))
      best_cop   <- copula::tCopula(param = rho, df = df, dim = dim_cop, dispstr = "un")
    }
    
    return(best_cop)
  }
  
  # Helper 2: Fit marginal KDEs to each column of data and return the kde1d objects
  fit_kdes <- function(data) {
    if (!all(sapply(data, is.numeric))) stop("All columns in data must be numeric.")
    lapply(1:ncol(data), function(i) kde1d::kde1d(data[, i]))
  }
  
  # Helper 3: Generate all 2^d corner combinations in [0,1]^d from bounds (d x 2).
  # For each variable i, transform l_i and u_i through the i-th marginal CDF only.
  # -Inf maps to 0 and +Inf maps to 1 (boundary values of any CDF).
  kde_list <- fit_kdes(data)
  d <- ncol(data)
  
  bounds_u_list <- lapply(1:d, function(i) {
    lo <- bounds[i, 1]
    up <- bounds[i, 2]
    u_lo <- if (is.infinite(lo) && lo < 0) 0 else kde1d::pkde1d(lo, kde_list[[i]])
    u_up <- if (is.infinite(up) && up > 0) 1 else kde1d::pkde1d(up, kde_list[[i]])
    c(u_lo, u_up)
  })
  variations_from_bounds_u <- setNames(
    as.data.frame(expand.grid(bounds_u_list)),
    paste0("u", 1:d)
  )
  
  # Helper 4: Inclusion-exclusion signs: (-1)^(number of lower bounds used)
  variations <- expand.grid(rep(list(1:2), d))
  signs <- apply(variations, 1, function(row) (-1)^sum(row == 1))
  
  # Main: sum copula CDF evaluations at all corners with inclusion-exclusion signs
  best_copula <- copulaCDF(data)
  sum(
    apply(variations_from_bounds_u, 1,
          function(i) copula::pCopula(u = as.numeric(i), copula = best_copula)) * signs
  )
}

jp(data = X, bounds = nom_bounds)/jp(data = X, bounds = den_bounds)

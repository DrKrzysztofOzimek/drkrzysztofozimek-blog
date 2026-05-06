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

  # Helper 2: Evaluate smoothed marginal CDFs at a given point
  get_pseudo_obs <- function(data, x_point) {
    if (!all(sapply(data, is.numeric))) stop("All columns in data must be numeric.")
    if (length(x_point) != ncol(data))  stop("x_point must match number of columns in data.")

    kde_list <- lapply(1:ncol(data), function(i) kde1d::kde1d(data[, i]))
    u_vec    <- mapply(kde1d::pkde1d, x_point, kde_list)
    names(u_vec) <- paste0("u", seq_along(u_vec))
    return(u_vec)
  }

  # Helper 3: Generate all 2^d corner combinations in [0,1]^d from smoothed bounds
  # FIX: iterate over ROWS of bounds (each row is a d-dimensional point)
  bounds_u <- do.call(rbind, lapply(1:nrow(bounds), function(i) {
    get_pseudo_obs(data, unlist(bounds[i, ]))
  }))
  bounds_u <- as.data.frame(bounds_u)
  colnames(bounds_u) <- paste0("u", 1:ncol(data))

  bounds_u_list <- mapply(
    function(lo, up) c(lo, up),
    bounds_u[1, ],   # lower pseudo-obs row
    bounds_u[2, ],   # upper pseudo-obs row
    SIMPLIFY = FALSE
  )
  variations_from_bounds_u <- setNames(
    as.data.frame(expand.grid(bounds_u_list)),
    paste0("u", 1:ncol(data))
  )

  # Helper 4: Inclusion-exclusion signs: (-1)^(number of lower bounds used)
  variations <- expand.grid(rep(list(1:2), ncol(data)))
  signs <- apply(variations, 1, function(row) (-1)^sum(row == 1))

  # Main: sum copula CDF evaluations at all corners with inclusion-exclusion signs
  best_copula <- copulaCDF(data)
  sum(
    apply(variations_from_bounds_u, 1,
          function(i) copula::pCopula(u = as.numeric(i), copula = best_copula)) * signs
  )
}

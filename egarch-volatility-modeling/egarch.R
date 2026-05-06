# ==============================================================================
# egarch.R  –  EGARCH(1,1) via rugarch
# Data     :  returns.csv  –  column "Returns"
# ==============================================================================
library(rugarch)
library(xts)

# ------------------------------------------------------------------------------
# 1. Load data
# ------------------------------------------------------------------------------
data    <- read.csv("returns.csv")
returns <- as.numeric(data$Returns)

# ------------------------------------------------------------------------------
# 2. Artificial daily date index
# ------------------------------------------------------------------------------
time_index  <- seq.Date(from = as.Date("2015-01-01"),
                        by   = "days",
                        length.out = length(returns))
returns_xts <- xts(returns, order.by = time_index)

# ------------------------------------------------------------------------------
# 3. Specify EGARCH(1,1) with Student-t innovations
# ------------------------------------------------------------------------------
spec <- ugarchspec(
  variance.model = list(model      = "eGARCH",
                        garchOrder = c(1, 1)),
  mean.model     = list(armaOrder  = c(0, 0),
                        include.mean = TRUE),
  distribution.model = "std"          # Student-t for heavy tails
)

# ------------------------------------------------------------------------------
# 4. Fit
# ------------------------------------------------------------------------------
fit_r <- ugarchfit(spec = spec, data = returns_xts)
show(fit_r)

# ------------------------------------------------------------------------------
# 5. One-step-ahead forecast
# ------------------------------------------------------------------------------
fc_r <- ugarchforecast(fit_r, n.ahead = 1)
cat("\n[R] Forecasted mean return (next 1 day):\n")
print(fc_r@forecast$seriesFor)
cat("\n[R] Forecasted sigma (next 1 day):\n")
print(fc_r@forecast$sigmaFor)

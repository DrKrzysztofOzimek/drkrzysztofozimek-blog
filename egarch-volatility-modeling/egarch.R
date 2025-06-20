library(rugarch)
library(xts)

# ------------------------------
# 1. Load return data from CSV
# ------------------------------
data <- read.csv("returns.csv") # or full path
returns <- as.numeric(data$Returns)

# ------------------------------
# 2. Set artificial daily date index
# ------------------------------
time_index <- seq.Date(from = as.Date("2015-01-01"), by = "days", length.out = length(returns))
returns <- xts(returns, order.by = time_index)

# ------------------------------
# 3. Specify EGARCH(1,1) model
# ------------------------------
spec <- ugarchspec(
  variance.model = list(model = "eGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(0, 0), include.mean = TRUE),
  distribution.model = "std" # Student-t for heavy tails
)

# ------------------------------
# 4. Fit the model
# ------------------------------
fit <- ugarchfit(spec = spec, data = returns)
show(fit)

# ------------------------------
# 5. Forecast
# ------------------------------
forecast <- ugarchforecast(fit, n.ahead = 1) # 10-step ahead forecast
forecast@forecast$seriesFor
forecast@forecast$sigmaFor

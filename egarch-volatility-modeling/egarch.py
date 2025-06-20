import pandas as pd
import numpy as np
from arch import arch_model

# ------------------------------
# 1. Load return data from CSV
# ------------------------------
data = pd.read_csv("returns.csv")  # Assumes a column named 'Returns'
returns = data["Returns"].dropna().astype(float)

# ------------------------------
# 2. Set artificial daily date index
# ------------------------------
start_date = pd.to_datetime("2015-01-01")
returns.index = pd.date_range(start=start_date, periods=len(returns), freq="D")

# ------------------------------
# 3. Specify EGARCH(1,1) model
# ------------------------------
model = arch_model(
    returns,
    vol="EGARCH",
    p=1,
    o=1,  # 'o' captures asymmetry (leverage effect)
    q=1,
    mean="Constant",
    dist="t",  # Student-t for heavy tails
)

# ------------------------------
# 4. Fit the model
# ------------------------------
fit = model.fit(disp="off")
print(fit.summary())

# ------------------------------
# 5. Forecast
# ------------------------------
horizon = 1
forecasts = fit.forecast(horizon=horizon)
forecast_mean = forecasts.mean.iloc[-1].values
forecast_variance = forecasts.variance.iloc[-1].values
forecast_sigma = np.sqrt(forecast_variance)

print(f"\nForecasted mean returns (next {horizon} days):")
print(forecast_mean)
print(f"\nForecasted sigma (next {horizon} days):")
print(forecast_sigma)

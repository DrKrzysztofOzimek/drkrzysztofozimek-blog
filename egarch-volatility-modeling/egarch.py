# ==============================================================================
# egarch.py  –  EGARCH(1,1) via arch
# Data      :  returns.csv  –  column "Returns"
# ==============================================================================
import pandas as pd
import numpy as np
from arch import arch_model

# ------------------------------------------------------------------------------
# 1. Load data
# ------------------------------------------------------------------------------
data    = pd.read_csv('returns.csv')
returns = data['Returns'].dropna().astype(float)

# ------------------------------------------------------------------------------
# 2. Artificial daily date index
# ------------------------------------------------------------------------------
start_date    = pd.to_datetime('2015-01-01')
returns.index = pd.date_range(start=start_date, periods=len(returns), freq='D')

# ------------------------------------------------------------------------------
# 3. Specify EGARCH(1,1) with Student-t innovations
# ------------------------------------------------------------------------------
model = arch_model(
    returns,
    vol  = 'EGARCH',
    p    = 1,
    o    = 1,   # asymmetry / leverage
    q    = 1,
    mean = 'Constant',
    dist = 't'  # Student-t for heavy tails
)

# ------------------------------------------------------------------------------
# 4. Fit
# ------------------------------------------------------------------------------
fit_py = model.fit(disp='off')
print(fit_py.summary())

# ------------------------------------------------------------------------------
# 5. One-step-ahead forecast
# ------------------------------------------------------------------------------
fc_py          = fit_py.forecast(horizon=1)
forecast_mean  = fc_py.mean.iloc[-1].values
forecast_sigma = np.sqrt(fc_py.variance.iloc[-1].values)

print()
print('[Python] Forecasted mean return (next 1 day):')
print(forecast_mean)
print('[Python] Forecasted sigma (next 1 day):')
print(forecast_sigma)

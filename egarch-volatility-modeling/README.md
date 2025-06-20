# 📊 Modeling Asymmetric Volatility with EGARCH

This folder contains R and Python code supporting the blog post on modeling time-varying volatility using the Exponential GARCH (EGARCH) model, which captures both volatility clustering and the leverage effect in financial return series.

🔗 **Read the full article:**  
[https://www.drkrzysztofozimek.com/egarch-asymmetric-volatility/](https://www.drkrzysztofozimek.com/egarch-asymmetric-volatility/)

---

## 💻 Code files included:

- `egarch.R` – R script to estimate an EGARCH(1,1) model using the `rugarch` package.
- `egarch.py` – Python script for EGARCH(1,1) estimation using the `arch` package.
- `returns.csv` – Sample return series data for model input and forecasting.

---

## 🧠 Summary

This code enables researchers, analysts, and students to estimate EGARCH models and explore how volatility responds asymmetrically to past shocks. It’s particularly useful for forecasting volatility in risk management, derivatives pricing, and volatility arbitrage strategies.

---

**Author:** [Dr Krzysztof Ozimek](https://www.drkrzysztofozimek.com)  
**License:** MIT

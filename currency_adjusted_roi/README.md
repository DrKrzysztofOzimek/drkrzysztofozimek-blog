# Currency-Adjusted ROI Calculator

This repository provides a method for calculating the return on investment (ROI) in a domestic currency when the investment is made in a foreign currency. It includes both R and Python implementations, a sample calculation, and an illustrative formula.

## 📈 Formula

rX = (rY - RXY) / (1 + RXY)

Where:
- rY = return in foreign currency
- RXY = exchange rate change (domestic per foreign unit)

## 📂 Contents

- `r_X.R` – ROI function in R
- `r_X.py` – ROI function in Python
- `example_calculation.csv` – Example with rY = 8%, e_XY = 4.00, e'_XY = 4.10
- `formula.png` – Final formula visualization
- `docs/blog_post_summary.md` – Blog post content (optional)

## ▶️ Example (Python)

```python
from r_X import r_X
r_X(0.08, 4.00, 4.10)
# Output: 0.05366

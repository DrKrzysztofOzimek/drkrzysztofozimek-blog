def hedge_coefficients(Pi, x0, h):
    # Central differences
    Pi_plus = Pi(x0 + h)
    Pi_minus = Pi(x0 - h)
    Pi_0 = Pi(x0)

    # Derivatives
    dPi_dx = (Pi_plus - Pi_minus) / (2 * h)
    d2Pi_dx2 = (Pi_plus - 2 * Pi_0 + Pi_minus) / (h**2)

    # Hedging coefficients
    m = -0.5 * d2Pi_dx2
    n = -dPi_dx - 2 * m * x0

    return n, m


# Example usage
Pi = lambda x: 0.5 * x**2 + 2 * x + 5
n, m = hedge_coefficients(Pi, x0=5, h=0.01)
print("n:", n, "m:", m)

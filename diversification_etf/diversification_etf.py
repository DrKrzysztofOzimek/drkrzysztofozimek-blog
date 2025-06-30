def specific_risk(sigma2, rho, gamma, n):
    if n <= 1:
        raise ValueError("n must be greater than 1.")

    rho_min = -1 / (n - 1)
    if rho < rho_min or rho > 1:
        print(f"Warning: rho should be in [{rho_min:.4f}, 1]. Given rho = {rho:.4f}")

    term1 = (1 - rho) * (gamma**2 + ((1 - gamma) ** 2) / (n - 1))
    risk = sigma2 * (term1 + rho)
    return risk


# Example usage
sigma2 = 0.04
rho = 0.1
gamma = 0.2
n = 10

print(specific_risk(sigma2, rho, gamma, n))

# =============================================================================
# Author:  Krzysztof Ozimek, PRM  (c) 2025
# License: Personal, non-commercial use only.  Redistribution or any
#          commercial use requires prior written permission from the author.
# =============================================================================

import numpy as np
import pandas as pd

def mp_randomness_test(data, B=1000):
    # ---------------------------------------------------------------------------
    # Marčenko–Pastur Simulation Test (Left-Tailed)
    #
    # Args:
    #   data : A numeric matrix or data frame (rows = observations, columns = variables).
    #   B    : Number of simualtions (default = 1000).
    #
    # Returns:
    #   p-value: Probability of observing a proportion of eigenvalues within
    #            MP bounds as low as or lower than the observed, assuming matrix is random.
    # ---------------------------------------------------------------------------

    if isinstance(data, pd.DataFrame):
        data = data.to_numpy()
    elif not isinstance(data, (np.ndarray, list)):
        raise ValueError("Input 'data' must be a data frame or matrix.")

    data = np.asarray(data)
    T_obs, N_vars = data.shape
    Q = T_obs / N_vars

    if Q < 1:
        raise ValueError("Q = T/N must be ≥ 1 for MP distribution to apply.")

    # Marčenko–Pastur bounds
    lambda_min = (1 - np.sqrt(1 / Q)) ** 2
    lambda_max = (1 + np.sqrt(1 / Q)) ** 2

    # Observed matrix
    data_std = (data - np.mean(data, axis=0)) / np.std(data, axis=0, ddof=1)
    cov_obs = np.cov(data_std, rowvar=False)
    eig_obs = np.linalg.eigvalsh(cov_obs)
    prop_obs = np.mean((eig_obs >= lambda_min) & (eig_obs <= lambda_max))

    # Simulate null distribution
    prop_null = np.empty(B)
    for b in range(B):
        Z = np.random.randn(T_obs, N_vars)
        Z_std = (Z - np.mean(Z, axis=0)) / np.std(Z, axis=0, ddof=1)
        cov_null = np.cov(Z_std, rowvar=False)
        eig_null = np.linalg.eigvalsh(cov_null)
        prop_null[b] = np.mean((eig_null >= lambda_min) & (eig_null <= lambda_max))

    # Left-tailed p-value: how rare is prop_obs under null?
    p_value = np.mean(prop_null <= prop_obs)
    return p_value

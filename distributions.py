"""Simulate straddle implied vol curves under Student-t, normal, lognormal, logistic, hyperbolic secant, generalized error, log-logistic, log hyperbolic secant, and log generalized error terminal distributions."""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats, optimize

from black_scholes import implied_vol_from_straddle

def student_t_scale(std: float, df: float) -> float:
    """Return Student-t scale parameter for given std and df > 2."""
    if df <= 2:
        raise ValueError("Student-t variance is finite only when df > 2.")
    if std <= 0:
        raise ValueError("Target standard deviation must be positive.")
    variance_factor = df / (df - 2)
    return std / np.sqrt(variance_factor)

def lognormal_parameters(mean: float, std: float) -> tuple[float, float]:
    """Return (mu, sigma) of lognormal with given mean and std."""
    if mean <= 0:
        raise ValueError("Lognormal mean must be positive.")
    if std <= 0:
        raise ValueError("Lognormal standard deviation must be positive.")
    variance_ratio = (std / mean) ** 2
    sigma_sq = np.log1p(variance_ratio)
    mu = np.log(mean) - 0.5 * sigma_sq
    sigma = np.sqrt(sigma_sq)
    return mu, sigma

def loglogistic_parameters(mean: float, std: float) -> tuple[float, float]:
    """Return (shape, scale) of log-logistic with given mean and std."""
    if mean <= 0:
        raise ValueError("Log-logistic mean must be positive.")
    if std <= 0:
        raise ValueError("Log-logistic standard deviation must be positive.")
    variance_ratio = (std / mean) ** 2
    target = 1.0 + variance_ratio
    def equation(scale: float) -> float:
        x = np.pi * scale
        return np.sin(x) / (x * np.cos(x)) - target
    lower = 1e-9
    upper = 0.5 - 1e-9
    scale = optimize.brentq(equation, lower, upper)
    shape = 1.0 / scale
    beta = mean * np.sin(np.pi * scale) / (np.pi * scale)
    return shape, beta

def log_hyperbolic_secant_parameters(mean: float, std: float) -> tuple[float, float]:
    """Return (mu, scale) of log hyperbolic secant with given mean and std."""
    if mean <= 0:
        raise ValueError("Log hyperbolic secant mean must be positive.")
    if std <= 0:
        raise ValueError("Log hyperbolic secant standard deviation must be positive.")
    variance_ratio = (std / mean) ** 2
    target = 1.0 + variance_ratio
    def equation(scale: float) -> float:
        half_pi_scale = 0.5 * np.pi * scale
        return (np.cos(half_pi_scale) ** 2) / np.cos(np.pi * scale) - target
    lower = 1e-9
    upper = 0.5 - 1e-9
    scale = optimize.brentq(equation, lower, upper)
    mu = np.log(mean * np.cos(0.5 * np.pi * scale))
    return mu, scale

def generalized_error_scale(std: float, power: float) -> float:
    """Return generalized error scale for given std and power > 0."""
    if std <= 0:
        raise ValueError("Generalized error standard deviation must be positive.")
    if power <= 0:
        raise ValueError("Generalized error power must be positive.")
    variance_factor = stats.gennorm.var(power, loc=0.0, scale=1.0)
    if not np.isfinite(variance_factor) or variance_factor <= 0:
        raise ValueError("Unable to compute a finite variance factor for the requested power.")
    return std / np.sqrt(variance_factor)

def log_generalized_error_parameters(mean: float, std: float, power: float) -> tuple[float, float]:
    """Return (mu, scale) of log generalized error with given mean, std, and power > 1."""
    if mean <= 0:
        raise ValueError("Log generalized error mean must be positive.")
    if std <= 0:
        raise ValueError("Log generalized error standard deviation must be positive.")
    if power <= 0:
        raise ValueError("Log generalized error power must be positive.")
    if power <= 1:
        raise ValueError("Log generalized error power must exceed 1 for finite exponential moments.")
    target_variance = std ** 2
    moment_cache: dict[float, float] = {}
    def exp_moment(scale_factor: float) -> float:
        if scale_factor not in moment_cache:
            moment = stats.gennorm.expect(lambda z: np.exp(scale_factor * z), args=(power,), loc=0.0, scale=1.0)
            if not np.isfinite(moment):
                moment_cache[scale_factor] = np.inf
            else:
                moment_cache[scale_factor] = moment
        return moment_cache[scale_factor]
    def variance_from_scale(scale: float) -> float:
        if scale <= 0:
            return 0.0
        first_moment = exp_moment(scale)
        if not np.isfinite(first_moment) or first_moment <= 0:
            return np.inf
        mu_local = np.log(mean) - np.log(first_moment)
        second_moment = exp_moment(2.0 * scale)
        if not np.isfinite(second_moment):
            return np.inf
        raw_second = np.exp(2.0 * mu_local) * second_moment
        return raw_second - mean ** 2
    lower = 1e-8
    lower_variance = variance_from_scale(lower)
    if not np.isfinite(lower_variance):
        raise RuntimeError("Failed to evaluate log generalized error variance near zero scale.")
    upper = 1e-2
    attempts = 0
    upper_variance = variance_from_scale(upper)
    while attempts < 80:
        if np.isfinite(upper_variance) and upper_variance > target_variance:
            break
        if not np.isfinite(upper_variance):
            upper = 0.5 * (lower + upper)
            if not upper > lower:
                raise RuntimeError("Unable to bracket finite variance for log generalized error calibration.")
        else:
            lower = upper
            lower_variance = upper_variance
            upper *= 2.0
            if upper > 10.0:
                raise RuntimeError("Failed to bracket target variance for log generalized error calibration.")
        upper_variance = variance_from_scale(upper)
        attempts += 1
    if not (np.isfinite(upper_variance) and upper_variance > target_variance and lower_variance <= target_variance):
        raise RuntimeError("Failed to locate a valid variance bracket for log generalized error calibration.")
    def objective(scale: float) -> float:
        return variance_from_scale(scale) - target_variance
    scale = optimize.brentq(objective, lower, upper)
    mu = np.log(mean) - np.log(exp_moment(scale))
    return mu, scale

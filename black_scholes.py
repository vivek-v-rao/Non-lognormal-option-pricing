import numpy as np
from scipy import stats
from scipy.optimize import brentq


def black_scholes_straddle_price(spot: float, strike: float,
    rate: float, time: float, vol: float) -> float:
    """Return Black-Scholes price of a straddle (call+put)"""
    if vol <= 0:
        return np.nan
    sqrt_time = np.sqrt(time)
    d1 = (np.log(spot / strike) + (rate + 0.5 * vol**2) * time) / (vol * sqrt_time)
    d2 = d1 - vol * sqrt_time
    call_price = spot * stats.norm.cdf(d1) - strike * np.exp(-rate * time) * stats.norm.cdf(d2)
    put_price = strike * np.exp(-rate * time) * stats.norm.cdf(-d2) - spot * stats.norm.cdf(-d1)
    return call_price + put_price

def implied_vol_from_straddle(price: float, spot: float,
    strike: float, rate: float, time: float) -> float:
    """Return implied volatility from a straddle price using Brent's method"""
    def objective(vol: float) -> float:
        return black_scholes_straddle_price(spot, strike, rate, time, vol) - price
    try:
        return brentq(objective, 1e-6, 5.0)
    except ValueError:
        return np.nan

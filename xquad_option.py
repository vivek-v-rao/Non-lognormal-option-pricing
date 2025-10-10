"""Simulate straddle implied volatility curves via quadrature under various terminal distributions."""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable, Dict, Iterable, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import integrate, optimize, stats

from black_scholes import implied_vol_from_straddle
from distributions import (
    generalized_error_scale,
    log_generalized_error_parameters,
    log_hyperbolic_secant_parameters,
    loglogistic_parameters,
    lognormal_parameters,
    student_t_scale,
)


QuadFunc = Callable[[float], float]


@dataclass
class QuadDistribution:
    label: str
    pdf: QuadFunc
    cdf_at_zero: float
    quantile_fn: Callable[[float], float]
    quad_points: tuple[float, ...] = ()
    _moment_cache: Dict[int, float] = field(default_factory=dict, init=False)

    def _quad(self, func: QuadFunc, lower: float, upper: float) -> float:
        """Integrate the provided function with optional break points."""
        opts = dict(limit=500, epsabs=1e-8, epsrel=1e-6)
        if self.quad_points:
            in_interval = [p for p in self.quad_points if lower < p < upper]
            if in_interval:
                opts["points"] = in_interval
        result, _ = integrate.quad(func, lower, upper, **opts)
        return result

    def expected_straddle_payoff(self, strike: float) -> float:
        """Return E[|S-K|] evaluated via quadrature."""
        return self.expected_call_payoff(strike) + self.expected_put_payoff(strike)

    def expected_call_payoff(self, strike: float) -> float:
        """Return E[(S-K)+] including any point mass at zero."""
        payoff = 0.0
        if strike < 0.0:
            payoff += self.cdf_at_zero * (-strike)
            payoff += self._quad(lambda y: (y - strike) * self.pdf(y), 0.0, np.inf)
        else:
            payoff += self._quad(lambda y: (y - strike) * self.pdf(y), max(strike, 0.0), np.inf)
        return payoff

    def expected_put_payoff(self, strike: float) -> float:
        """Return E[(K-S)+] including any point mass at zero."""
        if strike <= 0.0:
            return 0.0
        payoff = strike * self.cdf_at_zero
        payoff += self._quad(lambda y: (strike - y) * self.pdf(y), 0.0, strike)
        return payoff

    def raw_moment(self, order: int) -> float:
        """Return the non-negative raw moment using cached quadrature."""
        if order in self._moment_cache:
            return self._moment_cache[order]
        if order == 0:
            return 1.0
        value = self._quad(lambda y: (y ** order) * self.pdf(y), 0.0, np.inf)
        self._moment_cache[order] = value
        return value

    def quantile(self, prob: float) -> float:
        """Return the truncated quantile accounting for the point mass at zero."""
        if prob <= self.cdf_at_zero:
            return 0.0
        value = self.quantile_fn(prob)
        return 0.0 if not np.isfinite(value) else max(0.0, value)

    def summary(self, quantiles: Sequence[float]) -> dict[str, float | str]:
        """Summarize moments and requested quantiles for reporting."""
        mean = self.raw_moment(1)
        second = self.raw_moment(2)
        variance = second - mean ** 2
        std = math.sqrt(max(variance, 0.0))
        third = self.raw_moment(3)
        fourth = self.raw_moment(4)
        if std > 0:
            skew = (third - 3 * mean * second + 2 * mean ** 3) / std ** 3
            kurtosis = (fourth - 4 * mean * third + 6 * mean ** 2 * second - 3 * mean ** 4) / std ** 4 - 3.0
        else:
            skew = float("nan")
            kurtosis = float("nan")
        def format_label(q: float) -> str:
            text = f"q_{q:.3f}"
            text = text.rstrip('0').rstrip('.')
            return text if text else 'q_0'

        quantile_values: dict[str, float] = {}
        for q in sorted({float(q) for q in quantiles}):
            if not 0.0 <= q <= 1.0:
                continue
            label = format_label(q)
            value = self.quantile(q)
            if not np.isfinite(value):
                if q <= self.cdf_at_zero:
                    value = 0.0
                elif std > 0:
                    spread = 6.0 * std
                    value = max(0.0, mean - spread) if q < 0.5 else mean + spread
                else:
                    value = mean
            quantile_values[label] = value

        result: dict[str, float | str] = {
            "distribution": self.label,
            "mean": mean,
            "std": std,
            "skew": skew,
            "kurtosis": kurtosis,
        }
        result.update(quantile_values)
        result["prob_0"] = self.cdf_at_zero
        return result


def truncated_quantile(base_ppf: Callable[[float], float], cdf_zero: float) -> Callable[[float], float]:
    """Build a quantile function that clips probabilities at zero."""
    def _quantile(prob: float) -> float:
        if prob <= cdf_zero:
            return 0.0
        value = base_ppf(prob)
        return 0.0 if not np.isfinite(value) else max(0.0, value)

    return _quantile


def _truncated_mean(
    builder: Callable[[float], tuple[QuadFunc, float, Callable[[float], float], tuple[float, ...]]],
    loc: float,
) -> float:
    pdf, _cdf_zero, _quantile_fn, quad_points = builder(loc)
    opts = dict(limit=500, epsabs=1e-8, epsrel=1e-6)
    if quad_points:
        points = [p for p in quad_points if 0.0 < p < np.inf]
        if points:
            opts["points"] = points
    result, _ = integrate.quad(lambda y: y * pdf(y), 0.0, np.inf, **opts)
    return result


def calibrate_truncated_location(
    mean_target: float,
    initial_loc: float,
    scale_hint: float,
    builder: Callable[[float], tuple[QuadFunc, float, Callable[[float], float], tuple[float, ...]]],
) -> float:
    """Return a location shift ensuring the truncated mean matches mean_target."""
    initial_mean = _truncated_mean(builder, initial_loc)
    if abs(initial_mean - mean_target) <= 1e-8:
        return initial_loc

    step = max(abs(scale_hint), 1.0)
    lower = initial_loc
    upper = initial_loc
    f_lower = initial_mean - mean_target
    f_upper = f_lower

    attempts = 0
    while f_lower > 0.0 and attempts < 50:
        lower -= step
        f_lower = _truncated_mean(builder, lower) - mean_target
        attempts += 1

    attempts = 0
    while f_upper < 0.0 and attempts < 50:
        upper += step
        f_upper = _truncated_mean(builder, upper) - mean_target
        attempts += 1

    if f_lower > 0.0 or f_upper < 0.0:
        return initial_loc

    loc = optimize.brentq(
        lambda value: _truncated_mean(builder, value) - mean_target,
        lower,
        upper,
        xtol=1e-10,
        rtol=1e-10,
        maxiter=200,
    )
    return loc


def build_truncated_distribution(
    label: str,
    mean_target: float,
    initial_loc: float,
    scale_hint: float,
    builder: Callable[[float], tuple[QuadFunc, float, Callable[[float], float], tuple[float, ...]]],
) -> QuadDistribution:
    loc = calibrate_truncated_location(mean_target, initial_loc, scale_hint, builder)
    pdf, cdf_zero, quantile_fn, quad_points = builder(loc)
    return QuadDistribution(label, pdf, cdf_zero, quantile_fn, quad_points)


def integrate_straddle_curve(
    distribution: QuadDistribution,
    strikes: np.ndarray,
    discount_factor: float,
    spot: float,
    rate: float,
    time_to_maturity: float,
) -> list[dict[str, float | str]]:
    """Integrate straddle prices across strikes for one distribution."""
    records: list[dict[str, float | str]] = []
    for strike in strikes:
        call_payoff = distribution.expected_call_payoff(strike)
        put_payoff = distribution.expected_put_payoff(strike)
        straddle_payoff = call_payoff + put_payoff
        call_price = discount_factor * call_payoff
        put_price = discount_factor * put_payoff
        price = discount_factor * straddle_payoff
        implied_vol = implied_vol_from_straddle(price, spot, strike, rate, time_to_maturity)
        records.append(
            {
                "distribution": distribution.label,
                "strike": strike,
                "call_price": call_price,
                "put_price": put_price,
                "price": price,
                "implied_vol": implied_vol,
            }
        )
    return records


def main(
    dfs: Sequence[float] = (5.0,),
    include_normal: bool = True,
    include_lognormal: bool = True,
    include_logistic: bool = False,
    include_hyperbolic_secant: bool = True,
    include_log_logistic: bool = False,
    include_log_hyperbolic_secant: bool = True,
    generalized_error_powers: Sequence[float] = (),
    log_generalized_error_powers: Sequence[float] = (),
    mean_terminal: float = 100.0,
    std_terminal: float = 10.0,
    strike_start: float = 80.0,
    strike_end: float = 120.0,
    strike_step: float = 1.0,
    spot: float = 100.0,
    rate: float = 0.0,
    time_to_maturity: float = 1.0,
    enforce_risk_neutral_mean: bool = True,
    print_call_put: bool = False, # print prices of calls and puts separately
    print_straddles: bool = False,
    plot_terminal_distributions: bool = True,
    plot_implied_vol_markers: bool = False,
    zmin_dist: float = -4.0,
    zmax_dist: float =  4.0,
    summary_quantiles: Sequence[float] = (0.001, 0.01, 0.5, 0.99, 0.999),
) -> None:
    """Drive the quadrature-based straddle valuation experiment."""
    if any(df <= 2 for df in dfs):
        raise ValueError("All degrees of freedom must exceed 2 for finite variance.")
    if any(power <= 0 for power in generalized_error_powers):
        raise ValueError("Generalized error powers must be positive.")
    if any(power <= 1 for power in log_generalized_error_powers):
        raise ValueError("Log generalized error powers must exceed 1.")

    if not summary_quantiles:
        raise ValueError("summary_quantiles cannot be empty.")
    quantile_list = tuple(sorted({float(q) for q in summary_quantiles}))
    if any(q < 0.0 or q > 1.0 for q in quantile_list):
        raise ValueError("All summary quantiles must lie in [0, 1].")
    summary_quantiles = quantile_list

    strikes = np.arange(strike_start, strike_end + strike_step, strike_step)
    discount_factor = math.exp(-rate * time_to_maturity)
    if discount_factor <= 0.0:
        raise ValueError("Discount factor must be positive. Check rate and time_to_maturity.")
    risk_neutral_mean = spot / discount_factor
    mean_target = mean_terminal
    if enforce_risk_neutral_mean:
        mean_target = risk_neutral_mean
    if mean_target <= 0.0:
        raise ValueError("Resulting mean_target must be positive.")

    params = [
        ("dfs", ", ".join(f"{df:g}" for df in dfs)),
        ("include_normal", include_normal),
        ("include_lognormal", include_lognormal),
        ("include_logistic", include_logistic),
        ("include_hyperbolic_secant", include_hyperbolic_secant),
        ("include_log_logistic", include_log_logistic),
        ("include_log_hyperbolic_secant", include_log_hyperbolic_secant),
        ("generalized_error_powers", ", ".join(f"{p:g}" for p in generalized_error_powers) or "none"),
        ("log_generalized_error_powers", ", ".join(f"{p:g}" for p in log_generalized_error_powers) or "none"),
        ("mean_terminal", mean_terminal),
        ("enforce_risk_neutral_mean", enforce_risk_neutral_mean),
        ("risk_neutral_mean", risk_neutral_mean),
        ("mean_target", mean_target),
        ("std_terminal", std_terminal),
        ("spot", spot),
        ("rate", rate),
        ("time_to_maturity", time_to_maturity),
        ("strike_start", strike_start),
        ("strike_end", strike_end),
        ("strike_step", strike_step),
        ("print_call_put", print_call_put),
        ("plot_terminal_distributions", plot_terminal_distributions),
        ("plot_implied_vol_markers", plot_implied_vol_markers),
        ("zmin_dist", zmin_dist),
        ("zmax_dist", zmax_dist),
        ("summary_quantiles", ", ".join(f"{q:.3f}" for q in summary_quantiles)),
    ]
    params_df = pd.DataFrame(params, columns=["parameter", "value"])
    print("Simulation parameters:")
    print(params_df.to_string(index=False))

    distributions: list[QuadDistribution] = []

    for df in dfs:
        scale = student_t_scale(std_terminal, df)

        def builder(loc: float, df=df, scale=scale) -> tuple[QuadFunc, float, Callable[[float], float], tuple[float, ...]]:
            pdf = lambda x, df=df, loc=loc, scale=scale: 0.0 if x <= 0 else stats.t.pdf(x, df, loc=loc, scale=scale)
            cdf_zero = stats.t.cdf(0.0, df, loc=loc, scale=scale)
            quantile_fn = truncated_quantile(
                lambda p, df=df, loc=loc, scale=scale: stats.t.ppf(p, df, loc=loc, scale=scale),
                cdf_zero,
            )
            return pdf, cdf_zero, quantile_fn, ()

        distributions.append(
            build_truncated_distribution(f"Student-t df={df:g}", mean_target, mean_target, std_terminal, builder)
        )

    if include_normal:
        def builder(loc: float) -> tuple[QuadFunc, float, Callable[[float], float], tuple[float, ...]]:
            pdf = lambda x, loc=loc: 0.0 if x <= 0 else stats.norm.pdf(x, loc=loc, scale=std_terminal)
            cdf_zero = stats.norm.cdf(0.0, loc=loc, scale=std_terminal)
            quantile_fn = truncated_quantile(
                lambda p, loc=loc: stats.norm.ppf(p, loc=loc, scale=std_terminal),
                cdf_zero,
            )
            return pdf, cdf_zero, quantile_fn, ()

        distributions.append(
            build_truncated_distribution("Normal", mean_target, mean_target, std_terminal, builder)
        )

    if include_lognormal:
        mu, sigma = lognormal_parameters(mean_target, std_terminal)
        scale_param = math.exp(mu)
        pdf = lambda x, sigma=sigma, scale_param=scale_param: 0.0 if x <= 0 else stats.lognorm.pdf(x, s=sigma, scale=scale_param)
        quantile_fn = lambda p, sigma=sigma, scale_param=scale_param: stats.lognorm.ppf(p, s=sigma, scale=scale_param)
        distributions.append(QuadDistribution("Lognormal", pdf, 0.0, quantile_fn, quad_points=(0.0,)))

    if include_logistic:
        logistic_scale = std_terminal * math.sqrt(3.0) / math.pi
        def builder(loc: float, scale=logistic_scale) -> tuple[QuadFunc, float, Callable[[float], float], tuple[float, ...]]:
            pdf = lambda x, loc=loc, scale=scale: 0.0 if x <= 0 else stats.logistic.pdf(x, loc=loc, scale=scale)
            cdf_zero = stats.logistic.cdf(0.0, loc=loc, scale=scale)
            quantile_fn = truncated_quantile(
                lambda p, loc=loc, scale=scale: stats.logistic.ppf(p, loc=loc, scale=scale),
                cdf_zero,
            )
            return pdf, cdf_zero, quantile_fn, ()

        distributions.append(
            build_truncated_distribution("Logistic", mean_target, mean_target, std_terminal, builder)
        )

    if include_hyperbolic_secant:
        hyp_scale = std_terminal * 2.0 / math.pi
        def builder(loc: float, scale=hyp_scale) -> tuple[QuadFunc, float, Callable[[float], float], tuple[float, ...]]:
            pdf = lambda x, loc=loc, scale=scale: 0.0 if x <= 0 else stats.hypsecant.pdf(x, loc=loc, scale=scale)
            cdf_zero = stats.hypsecant.cdf(0.0, loc=loc, scale=scale)
            quantile_fn = truncated_quantile(
                lambda p, loc=loc, scale=scale: stats.hypsecant.ppf(p, loc=loc, scale=scale),
                cdf_zero,
            )
            return pdf, cdf_zero, quantile_fn, ()

        distributions.append(
            build_truncated_distribution("Hyperbolic Secant", mean_target, mean_target, std_terminal, builder)
        )

    for power in generalized_error_powers:
        scale = generalized_error_scale(std_terminal, power)
        def builder(loc: float, power=power, scale=scale) -> tuple[QuadFunc, float, Callable[[float], float], tuple[float, ...]]:
            pdf = lambda x, power=power, loc=loc, scale=scale: 0.0 if x <= 0 else stats.gennorm.pdf(x, power, loc=loc, scale=scale)
            cdf_zero = stats.gennorm.cdf(0.0, power, loc=loc, scale=scale)
            quantile_fn = truncated_quantile(
                lambda p, power=power, loc=loc, scale=scale: stats.gennorm.ppf(p, power, loc=loc, scale=scale),
                cdf_zero,
            )
            return pdf, cdf_zero, quantile_fn, ()

        label = f"Generalized Error power={power:g}"
        distributions.append(
            build_truncated_distribution(label, mean_target, mean_target, std_terminal, builder)
        )

    if include_log_logistic:
        shape, scale_param = loglogistic_parameters(mean_target, std_terminal)
        pdf = lambda x, shape=shape, scale_param=scale_param: 0.0 if x <= 0 else stats.fisk.pdf(x, shape, loc=0.0, scale=scale_param)
        quantile_fn = lambda p, shape=shape, scale_param=scale_param: stats.fisk.ppf(p, shape, loc=0.0, scale=scale_param)
        distributions.append(QuadDistribution("Log-Logistic", pdf, 0.0, quantile_fn, quad_points=(0.0,)))

    if include_log_hyperbolic_secant:
        mu, scale = log_hyperbolic_secant_parameters(mean_target, std_terminal)
        pdf = lambda x, mu=mu, scale=scale: 0.0 if x <= 0 else stats.hypsecant.pdf(math.log(x), loc=mu, scale=scale) / x
        quantile_fn = lambda p, mu=mu, scale=scale: math.exp(stats.hypsecant.ppf(p, loc=mu, scale=scale))
        distributions.append(QuadDistribution("Log Hyperbolic Secant", pdf, 0.0, quantile_fn, quad_points=(0.0,)))

    for power in log_generalized_error_powers:
        mu, scale = log_generalized_error_parameters(mean_target, std_terminal, power)
        pdf = lambda x, power=power, mu=mu, scale=scale: 0.0 if x <= 0 else stats.gennorm.pdf(math.log(x), power, loc=mu, scale=scale) / x
        quantile_fn = lambda p, power=power, mu=mu, scale=scale: math.exp(stats.gennorm.ppf(p, power, loc=mu, scale=scale))
        distributions.append(QuadDistribution(f"Log Generalized Error power={power:g}", pdf, 0.0, quantile_fn, quad_points=(0.0,)))

    summary_records = [dist.summary(summary_quantiles) for dist in distributions]
    summary_df = pd.DataFrame(summary_records)
    print("\nDistribution summary:")
    print(summary_df.to_string(index=False))

    curve_records: list[dict[str, float | str]] = []
    for dist in distributions:
        curve_records.extend(integrate_straddle_curve(dist, strikes, discount_factor, spot, rate, time_to_maturity))

    curve_df = pd.DataFrame(curve_records)
    curve_df["implied_stock"] = curve_df["call_price"] - curve_df["put_price"] + discount_factor * curve_df["strike"]
    if print_straddles:
        print("\nStraddle results:")
        straddle_view = curve_df[["distribution", "strike", "price", "implied_vol"]]
        print(straddle_view.to_string(index=False))

    if print_call_put:
        print("\nCall/Put results:")
        call_put_view = curve_df[
            [
                "distribution",
                "strike",
                "implied_vol",
                "implied_stock",
                "call_price",
                "put_price",
                "price",
            ]
        ].rename(columns={"price": "straddle_price"})
        print(call_put_view.to_string(index=False))

    plt.figure(figsize=(9, 5))
    marker = "o" if plot_implied_vol_markers else None
    for label, group in curve_df.groupby("distribution"):
        plt.plot(group["strike"], group["implied_vol"], marker=marker, linestyle="-", label=label)
    plt.xlabel("Strike")
    plt.ylabel("Implied Volatility")
    plt.title("Straddle Implied Volatility vs Strike (Quadrature)")
    plt.grid(True)
    plt.legend(title="Distribution")
    plt.tight_layout()

    if plot_terminal_distributions:
        plt.figure(figsize=(9, 5))
        overall_left: float | None = None
        overall_right: float | None = None
        for dist, summary in zip(distributions, summary_records):
            std = summary["std"]
            mean = summary["mean"]
            left = mean + (zmin_dist * std if (zmin_dist is not None and np.isfinite(std)) else 0.0)
            right = mean + (zmax_dist * std if (zmax_dist is not None and np.isfinite(std)) else dist.quantile(0.995))
            if zmin_dist is None:
                left = max(left, 0.0)
            left = max(0.0, left)
            if zmax_dist is None:
                if (not np.isfinite(right)) or right <= 0.0:
                    right = mean + 6.0 * std if std > 0 else mean
                right = max(right, strike_end * 1.1, mean + 3.0 * std if std > 0 else right)
            overall_left = left if overall_left is None else min(overall_left, left)
            overall_right = right if overall_right is None else max(overall_right, right)
            grid = np.linspace(left, right, 400)
            density = [dist.pdf(x) for x in grid]
            (line,) = plt.plot(grid, density, linewidth=1.5, label=dist.label)
            if dist.cdf_at_zero > 0.0:
                plt.vlines(
                    0.0,
                    0.0,
                    dist.cdf_at_zero,
                    colors=line.get_color(),
                    linestyles="--",
                    linewidth=1.2,
                )
        ax = plt.gca()
        if overall_left is not None and overall_right is not None:
            ax.set_xlim(overall_left, overall_right)
        plt.xlabel("Terminal Price")
        plt.ylabel("Density")
        plt.title("Terminal Price Densities (Quadrature)")
        plt.grid(True)
        plt.legend(title="Distribution")
        plt.tight_layout()

    plt.show()


if __name__ == "__main__":

    main()

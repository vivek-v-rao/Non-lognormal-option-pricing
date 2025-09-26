# Terminal Distributions

The following distributions are supported by xquad_option.py. Symbols below:

- Location: \(\mu\)
- Scale: \(\sigma\) (or named explicitly)
- Shape parameters are noted individually
- \(\Gamma(\cdot)\) is the gamma function, and \(\operatorname{sech}(x) = 1/\cosh(x)\)

## Student-t (location–scale)
- **PDF**: \[ f(x) = \frac{\Gamma\left(\frac{\nu+1}{2}\right)}{\Gamma\left(\frac{\nu}{2}\right) \sqrt{\nu \pi}\, \sigma} \left(1 + \frac{(x-\mu)^2}{\nu \sigma^2}\right)^{-\frac{\nu+1}{2}} \]
- **Mean**: \(\mu\) for \(\nu>1\)
- **Standard deviation**: \(\sigma \sqrt{\nu/(\nu-2)}\) for \(\nu>2\)
- **Skew**: 0 (\(\nu>3\))
- **Excess kurtosis**: \(6/(\nu-4)\) (\(\nu>4\))

## Normal (Gaussian)
- **PDF**: \[ f(x) = \frac{1}{\sqrt{2 \pi} \sigma} \exp\left(-\frac{(x-\mu)^2}{2 \sigma^2}\right) \]
- **Mean**: \(\mu\)
- **Standard deviation**: \(\sigma\)
- **Skew**: 0
- **Excess kurtosis**: 0

## Lognormal
Let \(\log X \sim \mathcal{N}(\mu, \sigma^2)\).
- **PDF**: \[ f(x) = \frac{1}{x \sigma \sqrt{2\pi}} \exp\left(-\frac{(\ln x - \mu)^2}{2 \sigma^2}\right),\; x>0 \]
- **Mean**: \(e^{\mu + \sigma^2/2}\)
- **Standard deviation**: \(\sqrt{(e^{\sigma^2}-1) e^{2\mu+\sigma^2}}\)
- **Skew**: \((e^{\sigma^2}+2) \sqrt{e^{\sigma^2}-1}\)
- **Excess kurtosis**: \(e^{4\sigma^2}+2e^{3\sigma^2}+3e^{2\sigma^2}-6\)

## Logistic
Scale parameter \(s>0\).
- **PDF**: \[ f(x) = \frac{\exp(-(x-\mu)/s)}{s \left(1+\exp(-(x-\mu)/s)\right)^2} \]
- **Mean**: \(\mu\)
- **Standard deviation**: \(s \pi / \sqrt{3}\)
- **Skew**: 0
- **Excess kurtosis**: \(6/5\)

## Hyperbolic Secant
Location \(\mu\), scale \(s>0\).
- **PDF**: \[ f(x) = \frac{1}{2 s} \operatorname{sech}\left(\frac{\pi (x-\mu)}{2 s}\right) \]
- **Mean**: \(\mu\)
- **Standard deviation**: \(s\)
- **Skew**: 0
- **Excess kurtosis**: 2

## Generalized Error (Exponential Power)
Shape \(\beta>0\), scale \(\alpha>0\).
- **PDF**: \[ f(x) = \frac{\beta}{2 \alpha \Gamma(1/\beta)} \exp\left(-\left|\frac{x-\mu}{\alpha}\right|^{\beta}\right) \]
- **Mean**: \(\mu\)
- **Variance**: \(\alpha^2 \Gamma(3/\beta)/\Gamma(1/\beta)\)
- **Standard deviation**: \(\alpha \sqrt{\Gamma(3/\beta)/\Gamma(1/\beta)}\)
- **Skew**: 0
- **Excess kurtosis**: \(\Gamma(5/\beta)\Gamma(1/\beta)/\Gamma(3/\beta)^2 - 3\)

## Log-Logistic (Fisk)
Shape \(\alpha>0\), scale \(\beta>0\).
- **PDF**: \[ f(x) = \frac{\alpha}{\beta} \left(\frac{x}{\beta}\right)^{\alpha-1} \left(1 + \left(\frac{x}{\beta}\right)^{\alpha}\right)^{-2},\; x>0 \]
- **Mean**: \(\beta \frac{\pi / \alpha}{\sin(\pi/\alpha)}\) for \(\alpha>1\)
- **Variance**: \(\beta^2 \left[\frac{2\pi/\alpha}{\sin(2\pi/\alpha)} - \left(\frac{\pi/\alpha}{\sin(\pi/\alpha)}\right)^2\right]\) for \(\alpha>2\)

## Log Hyperbolic Secant
Defined by \(Y = \exp(X)\) with \(X\) hyperbolic secant (location \(\mu\), scale \(s\)).
- **PDF**: \[ f(y) = \frac{1}{2 s y} \operatorname{sech}\left(\frac{\pi (\ln y - \mu)}{2 s}\right),\; y>0 \]
- **Mean**: \(e^{\mu} \operatorname{sech}(s)\) for \(|s|<\pi/2\)
- **Variance**: \(e^{2\mu} [\operatorname{sech}(2s) - \operatorname{sech}^2(s)]\)

## Log Generalized Error
Defined by \(Y = \exp(X)\), \(X\) generalized error.
- **PDF**: \[ f(y) = \frac{\beta}{2 \alpha y \Gamma(1/\beta)} \exp\left(-\left|\frac{\ln y - \mu}{\alpha}\right|^{\beta}\right),\; y>0 \]
- Higher moments depend on \(\beta\) and are omitted.

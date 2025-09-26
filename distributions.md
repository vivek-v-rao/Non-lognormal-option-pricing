# Terminal Distributions

The following distributions are supported by xquad_option.py. Symbols below:

- Location: mu
- Scale: sigma (or the scale parameter named explicitly)
- Shape parameters are noted individually
- Gamma(k) denotes the gamma function, and sech(x) = 1/cosh(x)

## Student-t (location-scale)
- PDF: f(x) = Gamma((nu+1)/2) / [ Gamma(nu/2) sqrt(nu*pi) * sigma ] * ( 1 + (x-mu)^2 / (nu*sigma^2) )^(-(nu+1)/2)
- Mean: mu (nu > 1)
- Standard deviation: sigma * sqrt(nu/(nu-2)) (nu > 2)
- Skew: 0 (nu > 3)
- Excess kurtosis: 6/(nu-4) (nu > 4)

## Normal (Gaussian)
- PDF: f(x) = (1/(sqrt(2*pi) * sigma)) * exp( -(x-mu)^2 / (2*sigma^2) )
- Mean: mu
- Standard deviation: sigma
- Skew: 0
- Excess kurtosis: 0

## Lognormal
Let log X ~ Normal(mu, sigma^2).
- PDF: f(x) = 1/(x * sigma * sqrt(2*pi)) * exp( -(ln x - mu)^2 / (2*sigma^2) ), for x > 0
- Mean: exp(mu + sigma^2/2)
- Standard deviation: sqrt( (exp(sigma^2) - 1) * exp(2*mu + sigma^2) )
- Skew: (exp(sigma^2) + 2) * sqrt(exp(sigma^2) - 1)
- Excess kurtosis: exp(4*sigma^2) + 2*exp(3*sigma^2) + 3*exp(2*sigma^2) - 6

## Logistic
Scale parameter s > 0.
- PDF: f(x) = exp(-(x-mu)/s) / [ s * (1 + exp(-(x-mu)/s))^2 ]
- Mean: mu
- Standard deviation: s * pi / sqrt(3)
- Skew: 0
- Excess kurtosis: 6/5

## Hyperbolic Secant
Location mu, scale s > 0 (SciPy convention).
- PDF: f(x) = 1/(2*s) * sech( pi*(x-mu) / (2*s) )
- Mean: mu
- Standard deviation: s
- Skew: 0
- Excess kurtosis: 2

## Generalized Error (Exponential Power)
Shape beta > 0, scale alpha > 0.
- PDF: f(x) = beta / (2*alpha*Gamma(1/beta)) * exp( -| (x-mu)/alpha |**beta )
- Mean: mu
- Variance: alpha**2 * Gamma(3/beta) / Gamma(1/beta)
- Standard deviation: alpha * sqrt( Gamma(3/beta) / Gamma(1/beta) )
- Skew: 0
- Excess kurtosis: Gamma(5/beta) * Gamma(1/beta) / Gamma(3/beta)**2 - 3

## Log-Logistic (Fisk)
Shape alpha > 0, scale beta > 0.
- PDF: f(x) = (alpha/beta) * (x/beta)**(alpha-1) * [ 1 + (x/beta)**alpha ]**(-2), for x > 0
- Mean: beta * (pi/alpha) / sin(pi/alpha)  (alpha > 1)
- Variance: beta**2 * [ (2*pi/alpha) / sin(2*pi/alpha) - ( (pi/alpha) / sin(pi/alpha) )**2 ]  (alpha > 2)

## Log Hyperbolic Secant
Y = exp(X) with X hyperbolic secant (location mu, scale s).
- PDF: f(y) = 1/(2*s*y) * sech( pi*(ln y - mu) / (2*s) ), for y > 0
- Mean: exp(mu) * sech(s)  (|s| < pi/2)
- Variance: exp(2*mu) * [ sech(2*s) - sech(s)**2 ]

## Log Generalized Error
Y = exp(X) with X generalized error (shape beta, scale alpha, location mu).
- PDF: f(y) = beta / (2*alpha*y*Gamma(1/beta)) * exp( -| (ln y - mu)/alpha |**beta ), for y > 0
- Higher moments depend on beta and are omitted.

# Terminal Distributions

The following distributions are supported by `xquad_option.py`. Symbols below:

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
- PDF: f(x) = (1/(âˆš(2Ï€) Â· Ïƒ)) Â· exp( -(x-Î¼)Â² / (2ÏƒÂ²) )
- Mean: Î¼
- Standard deviation: Ïƒ
- Skew: 0
- Excess kurtosis: 0

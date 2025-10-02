# Non-lognormal option pricing
Option prices and implied volatilities for various distributions of the terminal price, such as the [lognormal](https://en.wikipedia.org/wiki/Log-normal_distribution), [normal](https://en.wikipedia.org/wiki/Normal_distribution), [Student's t](https://en.wikipedia.org/wiki/Student%27s_t-distribution), [logistic](https://en.wikipedia.org/wiki/Logistic_distribution), [hyperbolic secant](https://en.wikipedia.org/wiki/Hyperbolic_secant_distribution), [log-logistic](https://en.wikipedia.org/wiki/Log-logistic_distribution), and log hyperbolic secant. The Student's t, hyperbolic, and logistic distributions all have excess kurtosis. If the terminal stock prices follow one of these distributions, the implied vol curve has a "smirk" resembling that of stock index options, rising in both tails but more sharply in the left tail. If the <b>log</b> prices follow these distributions, the vol curve is an almost symmetric smile. Run with<br>`python xquad_option.py`.

![Alt text](/implied_vol.png)
![Alt text](/densities.png)

```
Simulation parameters:
                    parameter                             value
                          dfs                                 5
               include_normal                              True
            include_lognormal                              True
             include_logistic                             False
    include_hyperbolic_secant                              True
         include_log_logistic                             False
include_log_hyperbolic_secant                              True
     generalized_error_powers                              none
 log_generalized_error_powers                              none
                mean_terminal                             100.0
                 std_terminal                              10.0
                         spot                             100.0
                         rate                               0.0
             time_to_maturity                               1.0
                 strike_start                              80.0
                   strike_end                             120.0
                  strike_step                               1.0
  plot_terminal_distributions                              True
     plot_implied_vol_markers                             False
                    zmin_dist                              -4.0
                    zmax_dist                               4.0
            summary_quantiles 0.001, 0.010, 0.500, 0.990, 0.999

Distribution summary:
         distribution       mean       std          skew      kurtosis   q_0.001    q_0.01      q_0.5     q_0.99    q_0.999
       Student-t df=5 100.000634  9.991504  3.836602e-02  4.990625e+00 54.349691 73.935364 100.000000 126.064636 145.650309
               Normal 100.000000 10.000000 -6.984919e-13 -1.389969e-08 69.097677 76.736521 100.000000 123.263479 130.902323
            Lognormal 100.000000 10.000000  3.010000e-01  1.615060e-01 73.108175 78.896644  99.503719 125.493172 135.429316
    Hyperbolic Secant 100.000001  9.999994  2.062134e-05  1.999716e+00 58.898734 73.557964 100.000000 126.442036 141.101266
Log Hyperbolic Secant 100.000000 10.000000  6.122449e-01  2.981234e+00 66.195530 76.554214  99.508597 129.345731 149.586547
```

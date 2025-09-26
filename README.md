# Non-lognormal option pricing
Option prices and implied volatilities for terminal distributions other than lognormal

![Alt text](/implied_vol.png)
![Alt text](/densities.png)

```
Simulation parameters:
                    parameter  value
                          dfs      6
               include_normal   True
            include_lognormal   True
             include_logistic  False
    include_hyperbolic_secant   True
         include_log_logistic  False
include_log_hyperbolic_secant   True
     generalized_error_powers   none
 log_generalized_error_powers   none
                mean_terminal  100.0
                 std_terminal   10.0
                         spot  100.0
                         rate    0.0
             time_to_maturity    1.0
                 strike_start   80.0
                   strike_end  120.0
                  strike_step    1.0
  plot_terminal_distributions   True
     plot_implied_vol_markers  False
                    zmin_dist   -4.0
                    zmax_dist    4.0

Distribution summary:
         distribution       mean       std          skew      kurtosis  min        max
       Student-t df=6 100.000186  9.997666  9.351761e-03  2.815445e+00  0.0 142.520090
               Normal 100.000000 10.000000 -6.984919e-13 -1.389969e-08  0.0 130.902323
            Lognormal 100.000000 10.000000  3.010000e-01  1.615060e-01  0.0 135.429316
    Hyperbolic Secant 100.000001  9.999994  2.062134e-05  1.999716e+00  0.0 141.101266
Log Hyperbolic Secant 100.000000 10.000000  6.122449e-01  2.981234e+00  0.0 149.586547
```

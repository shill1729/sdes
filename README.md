
# sdes

<!-- badges: start -->
<!-- badges: end -->

Solvers in R and c++ for SDEs via Euler-Maruyama and Runge-Kutta methods.

## Installation

You can install the current GitHub version via devtools

``` r
devtools::install_github("shill1729/sdes")
```


## Solving an SDE with jumps
We can simulate sample-paths of geometric Ito-Levy processes
by solving a SDE with jumps.
```r
library(sdes)
# Initial price and time-horizon
region <- c(spot = 10, time = 1)
# Infinitesimal drift and volatility coefficient functions
dynamics <- list(function(t, x) 0.05,
                 function(t, x) 0.1
)
# The jump dynamics: mean rate of jump, jump-size distribution
jumps <- list(lambda = function(t, x) 10,
              distr = "norm", 
              param = list(mean = -0.05, sd = 0.01)
)
s <- samplePathItoLevy(region, dynamics, jumps, n = 5000)
plot(s, type = "l")
```


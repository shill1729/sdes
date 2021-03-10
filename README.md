
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
by solving a SDE with jumps. Here is an example under the Merton model, together with a verification, visually, that the empirical density of the log-increments matches the exact model denisty.
```r
#============================================================
# Verification of empirical PDF of log-increments from
# a sample-path of a geometric Ito-Levy process matching
# the exact PDF, at least for the Merton model
#============================================================
library(sdes) # for generating sample-paths
library(findistr) # for PDFs
spot <- 100
maturity <- 1
n <- 10000
param <- c(0.1, 0.2, 15, -0.08, 0.01)
k <- maturity/n
# Initial price and time-horizon
region <- c(spot, maturity)
# Infinitesimal drift and volatility coefficient functions
dynamics <- list(function(t, x) param[1],
                 function(t, x) param[2]
)
# The jump dynamics: mean rate of jump, jump-size distribution
jumps <- list(lambda = function(t, x) param[3],
              distr = "norm",
              param = list(mean = param[4], sd = param[5])
)
s <- samplePathItoLevy(region, dynamics, jumps, n = n)
# Compute log-increments and densities
x <- diff(log(s$X))
epdf <- density(x)
mpdf <- dmerton(epdf$x, k, param)

# Plot sample-path and empirical vs exact density
par(mfrow = c(2, 1))
plot(s, type = "l")
plot(epdf$x, epdf$y, type = "l", ylim = c(0, max(epdf$y, mpdf)))
lines(epdf$x, mpdf, col = "blue", lty = "dashed")
```



# sdes

<!-- badges: start -->
<!-- badges: end -->

Solvers in R for SDEs via Euler-Maruyama and Runge-Kutta methods.

## Installation

You can install the current GitHub version via devtools

``` r
devtools::install_github("shill1729/sdes")
```
## Multiple models
The following models are currently available: geometric Brownian motion, Merton's jump diffusion, Heston stochastic volatility, and log-normal mixture diffusion
```r
# Simulating multiple models from the same point
spot <- 100
maturity <- 1
# Number of time-steps
n <- 5000
# Initial spot for Heston
v0 <- 0.5
# Mixture components
probs <- c(0.39, 0.01, 0.6)
mus <- c(0.02, -0.5, 0.1)
sigmas <- c(0.1, 1.5, 0.5)
# Drift and volatility for GBM and Merton
mu <- 0.9
volat <- 0.5
# Merton parameters
lambda <- 30
jm <- 0
jv <- 0.15
# Heston parameters: correlation, mean-reversion speed
# reversion level, vol-of-vol
rho <- -0.2
kappa <- 1.1
theta <- 0.6
xi <- sqrt(2*kappa*theta)/1.2
# Gathering together all the models
gbm <- list(name = "gbm", param = c(mu, volat))
merton <- list(name = "merton", param = c(mu, volat, lambda, jm, jv))
heston <- list(name = "heston", param = c(mu, rho, kappa, theta, xi, v0))
mixture <- list(name = "mixture", param = rbind(probs, mus, sigmas))
models <- list(gbm, merton, heston, mixture)
paths <- lapply(models, function(X) sdeSolve(spot, maturity, X, n, method = "em"))
par(mfrow = c(1, 1))
logRets <- lapply(paths, function(x) log(x$X/spot))
ub <- max(unlist(logRets))
lb <- min(unlist(logRets))
plot(paths[[1]]$t, logRets[[1]], type = "l", ylim = c(lb, ub))
for(i in 2:4)
{
  lines(paths[[i]]$t, logRets[[i]], col = i)
}
legend(x = "topleft", legend = c("gbm", "merton", "heston", "mixture"),
       col = c("black", "red", "green", "blue"), lty = 1, cex = 0.6)

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


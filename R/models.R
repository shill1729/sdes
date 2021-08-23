#' Simulate a geometric Brownian motion
#'
#' @param x0 initial value
#' @param tt time horizon
#' @param param vector of mean drift and volatility
#' @param n number of time sub-intervals
#' @param method "em" or "rk2"
#'
#' @description {Simulate a geometric Brownian motion sample path using
#' a numerical stochastic integrator.}
#' @return data.frame numeric
#' @export sde_gbm
sde_gbm <- function(x0, tt, param, n = 1000, method = "em")
{
  mu <- param[1]
  volat <- param[2]
  drift <- function(t, x) mu-0.5*volat^2
  diffusion <- function(t, x) volat
  dynamics <- list(drift, diffusion)
  s <- sde(0, tt, dynamics, n, method, plotG = FALSE)
  s$x <- x0*exp(s$x)
  plot(s, type  = "l", main = "GBM")
  return(s)
}

#' Simulate an OU process
#'
#' @param x0 initial value
#' @param tt time horizon
#' @param param vector of mean-reversion rate, level, and vol-of-vol
#' @param n number of time sub-intervals
#' @param method "em" or "rk2"
#'
#' @description {Simulate a OU process sample path using
#' a numerical stochastic integrator.}
#' @return data.frame numeric
#' @export sde_ou
sde_ou <- function(x0, tt, param, n = 1000, method = "em")
{
  kappa <- param[1]
  theta <- param[2]
  xi <- param[3]
  drift <- function(t, x) kappa*(theta-x)
  diffusion <- function(t, x) xi*sqrt(x)
  dynamics <- list(drift, diffusion)
  s <- sde(x0, tt, dynamics, n, method, plotG = FALSE)
  s$x <- sqrt(s$x)
  plot(s, type  = "l", main = "OU")
  return(s)
}

#' Simulate a log-normal mixture diffusion
#'
#' @param x0 initial value
#' @param tt time horizon
#' @param param matrix whose first row is probabilities, second drifts, third volatilities.
#' @param n number of time sub-intervals
#' @param method "em" or "rk2"
#'
#' @description {Simulate a log-normal mixture diffusion process sample path using
#' a numerical stochastic integrator.}
#' @return data.frame numeric
#' @export sde_mixture
sde_mixture <- function(x0, tt, param, n = 1000, method = "em")
{
  drift <- function(t, x) findistr::mixture_drift(t, x, param)-0.5*findistr::mixture_vol(t, x, param)^2
  diffusion <- function(t, x) findistr::mixture_vol(t, x, param)
  dynamics <- list(drift, diffusion)
  s <- sde(0, tt, dynamics, n, method, plotG = FALSE)
  s$x <- x0*exp(s$x)
  plot(s, type = "l", main = "Mixture")
  return(s)
}

#' Simulate a Merton jump diffusion
#'
#' @param x0 initial value
#' @param tt time horizon
#' @param param vector of mean drift, volatility,
#' mean-rate of jumps, mean jump size and jump-size sd.
#' @param n number of time sub-intervals
#'
#' @description {Simulate a Merton jump diffusion sample path using
#' a numerical stochastic integrator.}
#' @return data.frame numeric
#' @export sde_merton
sde_merton <- function(x0, tt, param, n = 1000)
{
  mu <- param[1]
  volat <- param[2]
  lambda <- param[3]
  alpha <- param[4]
  beta <- param[5]
  jumps <- list(distr = "norm", param = list(mean = alpha, sd = beta))
  eta <- exp(alpha+0.5*beta^2)-1
  drift <- function(t, x) mu-0.5*volat^2-lambda*eta
  diffusion <- function(t, x) volat
  lambdaf <- function(t, x) lambda

  s <- sde_levy(x0, tt, drift, diffusion, lambdaf, jumps, n)
  s$x <- x0*exp(s$x)
  return(s)
}

#' Simulate a Kou jump diffusion
#'
#' @param x0 initial value
#' @param tt time horizon
#' @param param vector of mean drift, volatility,
#' mean rate of jumps, probability of jump up, mean jump up size, and mean jump down size
#' @param n number of time sub-intervals
#'
#' @description {Simulate a Kou jump diffusion sample path using
#' a numerical stochastic integrator.}
#' @return data.frame numeric
#' @export sde_kou
sde_kou <- function(x0, tt, param, n = 1000)
{
  mu <- param[1]
  volat <- param[2]
  lambda <- param[3]
  p <- param[4]
  alpha <- param[5]
  beta <- param[6]
  jumps <- list(distr = "kou", param = list(p = p, alpha = alpha, beta = beta))
  eta <- meanJumpSize(jumps, 20000)
  drift <- function(t, x) mu-0.5*volat^2-lambda*eta
  diffusion <- function(t, x) volat
  lambdaf <- function(t, x) lambda
  s <- sde_levy(x0, tt, drift, diffusion, lambdaf, jumps, n)
  s$x <- x0*exp(s$x)
  return(s)
}

#' Simulate a displaced Kou jump diffusion
#'
#' @param x0 initial value
#' @param tt time horizon
#' @param param vector of mean drift, volatility,
#' mean rate of jumps, probability of jump up, mean jump up size, and mean jump down size, and displacements
#' @param n number of time sub-intervals
#'
#' @description {Simulate a displaced Kou jump diffusion sample path using
#' a numerical stochastic integrator.}
#' @return data.frame numeric
#' @export sde_dkou
sde_dkou <- function(x0, tt, param, n = 1000)
{
  mu <- param[1]
  volat <- param[2]
  lambda <- param[3]
  p <- param[4]
  alpha <- param[5]
  beta <- param[6]
  ku <- param[7]
  kd <- param[8]
  jumps <- list(distr = "dkou", param = list(p = p, alpha = alpha, beta = beta, ku = ku, kd = kd))
  eta <- meanJumpSize(jumps, 20000)
  drift <- function(t, x) mu-0.5*volat^2-lambda*eta
  diffusion <- function(t, x) volat
  lambdaf <- function(t, x) lambda
  s <- sde_levy(x0, tt, drift, diffusion, lambdaf, jumps, n)
  s$x <- x0*exp(s$x)
  return(s)
}

#' Simulate a uniform jump diffusion
#'
#' @param x0 initial value
#' @param tt time horizon
#' @param param vector of \code{(mu, volat, lambda, a, b)}
#' @param n number of time sub-intervals
#'
#' @description {Simulate a displaced uniform jump diffusion sample path using
#' a numerical stochastic integrator.}
#' @return data.frame numeric
#' @export sde_unijumps
sde_unijumps <- function(x0, tt, param, n = 1000)
{
  mu <- param[1]
  volat <- param[2]
  lambda <- param[3]
  a <- param[4]
  b <- param[5]

  jumps <- list(distr = "unif", param = list(min = a, max = b))
  eta <- meanJumpSize(jumps, 20000)
  drift <- function(t, x) mu-0.5*volat^2-lambda*eta
  diffusion <- function(t, x) volat
  lambdaf <- function(t, x) lambda
  s <- sde_levy(x0, tt, drift, diffusion, lambdaf, jumps, n)
  s$x <- x0*exp(s$x)
  return(s)
}

#' Simulate a Heston stochastic volatility model
#'
#' @param s0 initial price
#' @param v0 initial volatility
#' @param tt time horizon
#' @param param vector of mean drift, correlation, and Heston parameters
#' @param n number of time sub-intervals
#'
#' @description {Simulate a Heston stochastic volatility model sample path using
#' a numerical stochastic integrator.}
#' @return data.frame numeric
#' @export sde_heston
sde_heston <- function(s0, v0, tt, param, n = 1000)
{

  kappa <- param[1]
  theta <- param[2]
  xi <- param[3]
  mu <- param[4]
  rho <- param[5]
  drifts <- list(function(t, s, v) mu*s,
                 function(t, s, v) kappa*(theta-v)
                 )
  diffusions <- list(function(t, s, v) sqrt(v)*s,
                     function(t, s, v) xi*sqrt(v)
  )
  s <- sde_2d(IC = list(s = s0, v = v0), 0, tt, drifts, diffusions, rho, n, plotPhaseSpace = FALSE)
  plot(s$t, s$s, type = "l", main = "Heston SV")
  return(s)
}



#' Runge-Kutta solver for Ito SDEs
#'
#' @param x0 initial value of the process
#' @param t length to simulate over
#' @param drift the drift coefficient function of space-time
#' @param volat the volatility coefficient function of space-time
#' @param n number of sub-intervals in time-discretization
#'
#' @description {The stochastic Runge-Kutta for SDEs. See wikipedia for mathematical details}
#' @return data.frame
sde_rk2 <- function(x0, t, drift, volat, n = 50)
{
  h <- t/n
  tt <- seq(0, t, by = h)
  y <- matrix(0, nrow = length(tt))
  y[1] <- x0
  for(i in 1:(length(tt)-1))
  {
    dB <- sqrt(h)*stats::rnorm(1)
    gamma <- y[i]+drift(tt[i], y[i])*h+volat(tt[i], y[i])*sqrt(h)
    y[i+1] <- y[i]+drift(tt[i], y[i])*h+volat(tt[i], y[i])*dB+0.5*(volat(tt[i], gamma)-volat(tt[i], y[i]))*(dB^2-h)/sqrt(h)

  }
  return(data.frame(t = tt, X = y))
}

#' Euler-Maruyama solver for Ito SDEs
#'
#' @param x0 initial value of the process
#' @param t length to simulate over
#' @param drift the drift coefficient function of space-time
#' @param volat the volatility coefficient function of space-time
#' @param n number of sub-intervals in time-discretization
#'
#' @description {The stochastic Euler-Maruyama solver for SDEs. See wikipedia for details}
#' @return data.frame
sde_em <- function(x0, t, drift, volat, n = 50)
{
  h <- t/n
  tt <- seq(0, t, by = h)
  y <- matrix(0, nrow = length(tt))
  y[1] <- x0
  for(i in 1:(length(tt)-1))
  {

    y[i+1] <- y[i]+drift(tt[i], y[i])*h+volat(tt[i], y[i])*sqrt(h)*stats::rnorm(1)
  }
  return(data.frame(t = tt, X = y))
}


#' Simulate a sample path of an exponential Ito-Levy process with the Euler-Maruyama scheme
#'
#' @param spot the initial price of the share-price process.
#' @param t the terminal time to simulate until
#' @param mu the infinitesimal drift coefficient function of \code{(t, x)}
#' @param volat the infinitesimal volatility coefficient function of \code{(t, x)}
#' @param jumps a list of objects defining the jump-size distribution, can be NULL for purely continuous models
#' @param n number of sub-intervals in time-grid
#' @param M number of samples to use in Monte-Carlo estimate of \code{eta} for the jump-diffusion drift correction.
#'
#' @details { while \code{jumps} should contain
#' \itemize{
#' \item \code{lambda} the mean-rate of jumps function as a function of \code{(t, x)}.
#' \item \code{distr} the name of the distribution of the jump-sizes: "norm", "unif", "kou"
#' \item \code{param} named list of parameters for the distribution matching the input in \code{rdistr} for a given "distr".
#' }}
#' @return data.frame
jump_sde_em <- function(spot, t, mu, volat, jumps = NULL, n = 1000, M = 20000)
{
  lambda <- 0
  if(!is.null(jumps))
  {
    lambda <- jumps[[1]]
  }
  eta <- meanJumpSize(jumps, samples = M)
  # Solution vector to populate
  x <- matrix(0, nrow = n+1)
  x[1] <- 0
  # Fixed time-step size
  k <- t/n
  # Time grid
  tt <- seq(0, t, k)
  # Loop through time grid
  for(i in 2:(n+1))
  {
    if(!is.null(jumps))
    {
      # Extract jump-size distribution
      distr <- jumps$distr
      # Write function for calling rdistr(n)
      rdistr_name <- paste("r", distr, sep = "")
      rdistr <- function(y) do.call(what = rdistr_name, args = c(y, jumps$param))
      # Simulate number of jumps in single time-step
      m <- stats::rpois(1, lambda = lambda(tt[i-1], x[i-1])*k)
      # Compound Poisson sum of the random jump-amplitudes
      cps_pois <- sum(rdistr(m))
    } else{
      cps_pois <- 0
    }
    x[i] <- x[i-1] + (mu(tt[i-1], x[i-1])-0.5*volat(tt[i-1], x[i-1])^2-eta*lambda(tt[i-1], x[i-1]))*k+volat(tt[i-1], x[i-1])*stats::rnorm(1, 0, sd = sqrt(k))+cps_pois

  }
  X <- data.frame(t = tt, X = spot*exp(x))
  return(X)
}




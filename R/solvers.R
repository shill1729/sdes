#' Runge-Kutta solver for Ito SDEs
#'
#' @param x0 initial value of the process
#' @param t length to simulate over
#' @param drift the drift coefficient function of space-time
#' @param volat the volatility coefficient function of space-time
#' @param n number of sub-intervals in time-discretization
#'
#' @description {The stochastic Runge-Kutta for SDEs. See wikipedia for details}
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





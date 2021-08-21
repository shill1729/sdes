#' Runge-Kutta solver for Ito SDEs
#'
#' @param x0 initial value of the process
#' @param tt length to simulate over
#' @param drift the drift coefficient function of space-time in the order \code{(t, x)}
#' @param diffusion the diffusion coefficient function of space-time, in the order \code{(t, x)}
#' @param n number of sub-intervals in time-discretization
#'
#' @description {The stochastic Runge-Kutta for SDEs. See wikipedia for mathematical details.}
#' @return data.frame
sde_rk2 <- function(x0, tt, drift, diffusion, n = 1000)
{
  h <- tt/n
  timeGrid <- seq(0, tt, by = h)
  y <- matrix(0, nrow = length(tt))
  y[1] <- x0
  for(i in 1:(length(tt)-1))
  {
    dB <- sqrt(h)*stats::rnorm(1)
    gamma <- y[i]+drift(timeGrid[i], y[i])*h+diffusion(timeGrid[i], y[i])*sqrt(h)
    y[i+1] <- y[i]+drift(timeGrid[i], y[i])*h+diffusion(timeGrid[i], y[i])*dB+0.5*(diffusion(timeGrid[i], gamma)-diffusion(timeGrid[i], y[i]))*(dB^2-h)/sqrt(h)

  }
  return(data.frame(t = timeGrid, x = y))
}

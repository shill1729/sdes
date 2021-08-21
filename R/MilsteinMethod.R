#' Milstein method for solving autonomous SDEs
#'
#' @param x0 the initial point of the process
#' @param tt the time horizon to solve over
#' @param drift the infinitesimal drift coefficient, a single-variable function of space
#' @param diffusion the infinitesimal diffusion coefficient, a single-variable function of space
#' @param diff1 the derivative of the diffusion coefficient function, can be \code{NULL} for numerical approximations.
#' @param n number of sub-intervals to use in time-interval partition
#'
#' @description {A straight forward implementation of the Milstein method for solving SDEs via
#' stochastic integration. The Milstein scheme has both weak and strong order of convergence on
#' the order of the time-step, which is superior to the Euler-Maruyama scheme's square root of time order for
#' strong convergence.}
#'
#' @details {The derivative of the diffusivity coefficient function may be
#' \code{NULL} and approximated numerically.}
#' @return data.frame
sde_milstein <- function(x0, tt, drift, diffusion, diff1 = NULL, n = 1000)
{
  k <- tt/n
  if(is.null(diff1))
  {
    diff1 <- function(x) (diffusion(x+k)-diffusion(x-k))/(2*k)
  }
  y <- matrix(0, nrow = n+1)
  y[1] <- x0
  for(i in 2:(n+1))
  {
    dw <- stats::rnorm(1, 0, sqrt(k))
    y[i] <- y[i-1]+drift(y[i-1])*k+diffusion(y[i-1])*dw+0.5*diffusion(y[i-1])*diff1(y[i-1])*(dw^2-k)
  }
  return(data.frame(t = (0:n)*k, x = y))
}

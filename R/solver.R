#' Solve a one-dimensional SDE for arbitrary coefficient functions under three different numerical schemes
#'
#' @param x0 initial point of the process
#' @param tt the time horizon to solve over
#' @param dynamics list of functions defining dynamics of the stochastic process, see details
#' @param n number of time sub-intervals in the partition
#' @param method the numerical integration method to use: "em", "rk2", or "milstein", see details
#' @param plotG whether to plot the sample path in the function body
#'
#' @description {Generate sample paths of processes solving one-dimensional SDEs with near arbitrary
#' coefficient functions. Three schemes are available: Euler-Maruyama, RK2 and Milstein's method.}
#' @details {For \code{method = "milstein"} one must pass single-variable coefficient functions and
#' one can optionally pass the derivative of the volatility function if it is known, otherwise use a numeric
#' approximation, as one of the two is required in the Milstein scheme. For \code{method = "em"} or "rk2", the
#' coefficient functions ought to be functions of \code{(t, x)} but wrappers are made if this is not done.
#'
#' The \code{dynamics} list must be a list of either 2 or 3 functions defining the infinitesimal drift
#' function and the infinitesimal volatility function and optionally its derivative for the Milstein method.}
#' @return numeric/data.frame
#' @export sde
sde <- function(x0, tt, dynamics, n = 1000, method = "milstein", plotG = TRUE)
{
  drift <- 0
  diffusion <- 0
  if(method == "em" || method == "rk2")
  {
    if(length(dynamics) < 2)
    {
      stop("'dynamics' must be a list of two functions defining the drift and volatility")
    }
    if(length(formals(dynamics[[1]])) == 1 || length(formals(dynamics[[1]])) == 1 )
    {
      warning("try to use milstein for autonomous SDEs; creating wrappers for convenience")
      drift <- function(t, x) dynamics[[1]](x)
      diffusion <- function(t, x) dynamics[[2]](x)
    } else
    {
      drift <- dynamics[[1]]
      diffusion <- dynamics[[2]]
    }
    input <- list(x0 = x0, tt = tt, drift = drift, diffusion = diffusion, n = n)
  } else if(method == "milstein")
  {

    if(length(formals(dynamics[[1]])) > 1 || length(formals(dynamics[[1]])) > 1 )
    {
      stop("Use em/rk2 for time-inhomogeneous SDEs; milstein is only for autonomous SDEs here")
    }
    drift <- dynamics[[1]]
    diffusion <- dynamics[[2]]
    # If only drift and volatility are passed for milstein method, assume numerically approx
    # to derivative of volatility:
    if(length(dynamics) == 2)
    {
      diff1 <- NULL

    } else if(length(dynamics) == 3)
    {
      diff1 <- dynamics[[3]]
    }
    input <- list(x0 = x0, tt = tt, drift = drift, diffusion = diffusion, diff1 = diff1, n = n)
  }
  x <- do.call(paste("sde_", method, sep = ""), args = input)
  # graphics::par(mfrow = c(1, 1))
  if(plotG)
  {
    plot(x, type  = "l")
  }
  return(x)
}






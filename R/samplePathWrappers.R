#' Simulate sample path of Ito-process by solving a SDE
#'
#' @param region initial point and length of time as a vector
#' @param dynamics list of drift and volatility functions of \code{(t, x)}, use \code{NULL}
#' for Brownian motion default dynamics
#' @param n number of variates
#' @param method either "em" for Euler-Maruyama scheme or rk2 for Runge-Kutta scheme of order 2.
#'
#' @description {Simulate a sample path of a solution to an SDE.
#' Wrapper to internal solvers.}
#' @return data.frame of time variable and space variable
#' @export samplePathIto
samplePathIto <- function(region, dynamics = NULL, n = 1000, method = "rk2")
{

  if(method != "rk2" && method != "em")
  {
    stop("Input 'method' must be 'rk2' or 'em'")
  }
  # Default Brownian motion
  if(is.null(dynamics))
  {
    dynamics <- list(function(t, x) 0,
                     function(t, x) 1
                     )
  }
  s <- do.call(what = paste("sde_", method, sep =""),
               list(x0 = region[1],
                    t = region[2],
                    drift = dynamics[[1]],
                    volat = dynamics[[2]],
                    n = n
               )
  )
  return(s)
}

#' Simulate a sample-path of an exponential Ito-Levy process by solving a SDE with jumps
#'
#' @param region vector of initial price and time-horizon
#' @param dynamics list of drift and volatility coefficient functions
#' @param jumps the jump-specification list
#' @param n number of variates to discretize time by
#'
#' @description {Solve an SDE with jumps via the Euler-Maruyama discretization and hence
#' simulate a sample-path of the stochastic process.}
#' @details { The SDE solved is for the log-dynamics and
#' then the process is exponentiated and returned. Therefore the drift and volatility coefficents passed
#' must be those figuring in the SDE (without jumps for the sake of convenience, here)
#' \eqn{dS_t=\mu(t, S_t) S_t dt+\sigma(t, S_t) S_t dB_t}.
#'
#' Further, the argument \code{jumps} must be a named list containing
#' \itemize{
#' \item \code{lambda} the mean-rate of jumps function as a function of \code{(t, x)}.
#' \item \code{distr} the name of the distribution of the (log-)jump-sizes
#' \item \code{param} a named list of the parameters of the distribution matching
#' the names and order of the arguments in e.g. \code{runif}, etc.}
#' }
#' @return data.frame of time and price
#' @export samplePathItoLevy
samplePathItoLevy <- function(region, dynamics, jumps, n)
{
  spot <- region[1]
  t <- region[2]
  mu <- dynamics[[1]]
  volat <- dynamics[[2]]
  w <- jump_sde_em(spot, t, mu, volat, jumps, n)
  return(w)
}

#' Solve arbitrary first order SDE systems
#'
#' @param drifts a list of drift functions with as many variables as IC, and in the same order with the names
#' @param diffusions a list of volatility functions with as many variables as IC, and in the same order with the names
#' @param IC a list of variables with the same names as the non-time input as drifts, in the same order with initial values
#' @param t0 initial time
#' @param tn ending time
#' @param rho the correlation factor
#' @param n number of sub-intervals in time-grid
#'
#' @description {A general Euler scheme implementation for arbitrary first-order SDE systems with two state variables.}
#' @details {The list of functions must have syntax \code{h(t,x,y)} for each element matching the elements of \code{IC} as \code{list(x = y0, y = y0)} for example.}
#' @return data.frame
#' @export samplePathTwoStateSystem
samplePathTwoStateSystem <- function(drifts, diffusions, IC, t0 = 0, tn = 1, rho = 0, n = 1000)
{
  # Time step size
  h <- (tn-t0)/n
  # Time grid
  tt <- seq(t0, tn, length.out = n+1)
  # Number of state variables
  m <- length(IC)
  # state variables
  states <- list()
  for(i in 1:m)
  {
    states[[i]] <- matrix(0, n+1)
    # Initialize ICs
    states[[i]][1] <- IC[[i]]
  }
  # If processes are driven by the same Brownian motion, simulate these beyond
  # the scope of the loop
  z <- findistr::rcornorm(n, rho)
  # Over every time-node
  for(j in 2:(n+1))
  {
    # For each state variable
    for(i in 1:m)
    {

      # Gather input
      input <- list()
      input[[1]] <- tt[j-1]
      for(l in 1:m)
      {
        input[[l+1]] <- states[[l]][j-1]
      }
      names(input) <- c("t", names(IC))
      states[[i]][j] <- states[[i]][j-1]+h*do.call(what = drifts[[i]], args = input)+sqrt(h)*z[j-1, i]*do.call(what = diffusions[[i]], args = input)
    }
  }
  states <- data.frame(time = tt, do.call(cbind, states))
  names(states) <- c("time", names(IC))
  return(states)
}

#' Solve arbitrary first order SDE systems
#'
#' @param drifts a list of drift functions with as many variables as IC, and in the same order with the names
#' @param diffusions a list of volatility functions with as many variables as IC, and in the same order with the names
#' @param IC a list of variables with the same names as the non-time input as drifts, in the same order with initial values
#' @param t0 initial time
#' @param tn ending time
#' @param independent boolean to simulate processes driven by one Brownian motion (false) or as many independent Brownian motions as there are state variables
#' @param n number of sub-intervals in time-grid
#'
#' @description {A general Euler scheme implementation for arbitrary first-order SDE systems of finite dimension.}
#' @details {The list of functions must have syntax \code{h(t,x,y,z)} for each element matching the elements of \code{IC} as \code{list(x = y0, y = y0, z = z0)} for example.}
#' @return data.frame
#' @export samplePathSystem
samplePathSystem <- function(drifts, diffusions, IC, t0 = 0, tn = 1, independent = FALSE, n = 1000)
{
  # Time step size
  h <- (tn-t0)/n
  # Time grid
  tt <- seq(t0, tn, length.out = n+1)
  # Number of state variables
  m <- length(IC)
  # state variables
  states <- list()
  for(i in 1:m)
  {
    states[[i]] <- matrix(0, n+1)
    # Initialize ICs
    states[[i]][1] <- IC[[i]]
  }
  # If processes are driven by the same Brownian motion, simulate these beyond
  # the scope of the loop
  if(!independent)
  {
    z <- stats::rnorm(n)
  }
  # Over every time-node
  for(j in 2:(n+1))
  {
    # For each state variable
    for(i in 1:m)
    {
      # If simulating independently driven processes,
      # use a new Gaussian step each state
      if(independent)
      {
        z <- stats::rnorm(n)
      }
      # Gather input
      input <- list()
      input[[1]] <- tt[j-1]
      for(l in 1:m)
      {
        input[[l+1]] <- states[[l]][j-1]
      }
      names(input) <- c("t", names(IC))
      states[[i]][j] <- states[[i]][j-1]+h*do.call(what = drifts[[i]], args = input)+sqrt(h)*z[j-1]*do.call(what = diffusions[[i]], args = input)
    }
  }
  states <- data.frame(time = tt, do.call(cbind, states))
  names(states) <- c("time", names(IC))
  return(states)
}


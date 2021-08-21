#' Euler-Maruyama solver for Ito SDEs
#'
#' @param x0 initial value of the process
#' @param tt length to simulate over
#' @param drift the drift coefficient function of space-time in the order \code{(t, x)}
#' @param diffusion the diffusion coefficient function of space-time, in the order \code{(t, x)}
#' @param n number of sub-intervals in time-discretization
#'
#' @description {The stochastic Euler-Maruyama solver for SDEs. See wikipedia for details}
#' @return data.frame
sde_em <- function(x0, tt, drift, diffusion, n = 1000)
{
  h <- tt/n
  timeGrid <- seq(0, tt, by = h)
  y <- matrix(0, nrow = n+1)
  y[1] <- x0
  for(i in 1:(length(timeGrid)-1))
  {

    y[i+1] <- y[i]+drift(timeGrid[i], y[i])*h+diffusion(timeGrid[i], y[i])*sqrt(h)*stats::rnorm(1)
  }
  return(data.frame(t = timeGrid, x = y))
}

#' Simulate a sample path of an exponential Ito-Levy process with the Euler-Maruyama scheme
#'
#' @param spot the initial price of the share-price process.
#' @param tt the terminal time to simulate until
#' @param drift the infinitesimal drift coefficient function of \code{(t, x)}
#' @param diffusion the infinitesimal volatility coefficient function of \code{(t, x)}
#' @param lambda the mean jump-rate function of \code{(t, x)}
#' @param jumps a list of objects defining the jump-size distribution
#' @param n number of sub-intervals in time-grid
#' @param M number of samples to use in Monte-Carlo estimate of \code{eta} for the jump-diffusion drift correction.
#'
#' @details { while \code{jumps} should contain
#' \itemize{
#' \item \code{distr} the name of the distribution of the jump-sizes: "norm", "unif", "kou"
#' \item \code{param} named list of parameters for the distribution matching the input in \code{rdistr} for a given "distr".
#' }}
#' @return data.frame
#' @export sde_levy
sde_levy <- function(spot, tt, drift, diffusion, lambda, jumps, n = 1000, M = 20000)
{
  eta <- meanJumpSize(jumps, samples = M)
  # Solution vector to populate
  x <- matrix(0, nrow = n+1)
  x[1] <- 0
  # Fixed time-step size
  k <- tt/n
  # Time grid
  timeGrid <- seq(0, tt, k)
  # Loop through time grid
  for(i in 2:(n+1))
  {
    # Extract jump-size distribution
    distr <- jumps$distr
    # Write function for calling rdistr(n)
    rdistr_name <- paste("r", distr, sep = "")
    rdistr <- function(y) do.call(what = rdistr_name, args = c(y, jumps$param))
    # Simulate number of jumps in single time-step
    m <- stats::rpois(1, lambda = lambda(timeGrid[i-1], x[i-1])*k)
    # Compound Poisson sum of the random jump-amplitudes
    cps_pois <- sum(rdistr(m))
    x[i] <- x[i-1] + (drift(timeGrid[i-1], x[i-1])-0.5*diffusion(timeGrid[i-1], x[i-1])^2-eta*lambda(timeGrid[i-1], x[i-1]))*k+diffusion(timeGrid[i-1], x[i-1])*stats::rnorm(1, 0, sd = sqrt(k))+cps_pois

  }
  if(distr == "norm")
  {
    samplePathName <- "Merton"
  } else if(distr == "kou")
  {
    samplePathName <- "Kou"
  } else if(distr == "dkou")
  {
    samplePathName <- "Displaced Kou"
  } else if(distr == "unif")
  {
    samplePathName <- "Uniform jumps"
  }
  X <- data.frame(t = timeGrid, x = spot*exp(x))
  plot(X, type  = "l", main = samplePathName)
  return(X)
}

#' Solve arbitrary first order SDE systems
#'
#' @param IC a list of variables with the same names as the non-time input as drifts, in the same order with initial values
#' @param t0 initial time
#' @param tn ending time
#' @param drifts a list of drift functions with as many variables as IC, and in the same order with the names
#' @param diffusions a list of volatility functions with as many variables as IC, and in the same order with the names
#' @param rho the correlation factor
#' @param n number of sub-intervals in time-grid
#' @param plotPhaseSpace whether to plot the phase-space graph in the function body
#'
#' @description {A general Euler scheme implementation for arbitrary first-order SDE systems with two state variables.}
#' @details {The list of functions must have syntax \code{h(t,x,y)} for each element matching the elements of \code{IC} as \code{list(x = y0, y = y0)} for example.}
#' @return numeric data.frame
#' @export sde_2d
sde_2d <- function(IC, t0 = 0, tn = 1, drifts, diffusions, rho = 0, n = 1000, plotPhaseSpace = TRUE)
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
  # Generate correlated BMs
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
  states <- data.frame(t = tt, do.call(cbind, states))
  names(states) <- c("t", names(IC))
  if(plotPhaseSpace)
  {
    plot(states[,-1], type = "l", main = "Phase-space")
  }
  return(states)
}



#' Solve arbitrary first order SDE systems
#'
#' @param IC a list of variables with the same names as the non-time input as drifts, in the same order with initial values
#' @param t0 initial time
#' @param tn ending time
#' @param drifts a list of drift functions with as many variables as IC, and in the same order with the names
#' @param diffusions a list of volatility functions with as many variables as IC, and in the same order with the names
#' @param independent boolean to simulate processes driven by one Brownian motion (false) or as many independent Brownian motions as there are state variables
#' @param n number of sub-intervals in time-grid
#'
#' @description {A general Eule-Maruyama scheme implementation for arbitrary first-order SDE systems of finite dimension.}
#' @details {The list of functions must have syntax \code{h(t,x,y,z)} for each element matching the elements of \code{IC} as \code{list(x = y0, y = y0, z = z0)} for example.}
#' @return data.frame
#' @export sde_system
sde_system <- function(IC, t0 = 0, tn = 1, drifts, diffusions, independent = FALSE, n = 1000)
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


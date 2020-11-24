#' Simulate sample path of Ito-process
#'
#' @param region initial point and length of time as a vector
#' @param dynamics list of drift and volatility functions of \code{(t, x)}
#' @param control list of variates \code{n}, method and engine
#'
#' @description {Simulate a sample path of a solution to an SDE.
#' Wrapper to internal solvers.}
#' @details {The argument \code{control} must be a named list of
#' \itemize{
#' \item \code{n} number of variates in simulation, default 1000
#' \item \code{method} either "em" for Euler-Maruyama or "rk2" for Runge-Kutta order 2
#' \item \code{engine} either "r" or "cpp"}}
#' @return data.frame of time variable and space variable
#' @export samplePathIto
samplePathIto <- function(region, dynamics = NULL, control = list(n = 1000, method = "em", engine = "cpp"))
{
  if(!all.equal(names(control), c("n", "method", "engine")))
  {
    stop("control must be a named list")
  }
  n <- control$n
  method <- control$method
  engine <- control$engine


  if(method != "rk2" && method != "em")
  {
    stop("Input 'method' must be 'rk2' or 'em'")
  }
  if(engine != "r" && engine != "cpp")
  {
    stop("Input 'engine' must be 'r' or 'cpp'")
  }
  # Default Brownian motion
  if(is.null(dynamics) && engine == "r")
  {
    dynamics <- list(function(t, x) 0,
                     function(t, x) 1
                     )
  }
  if(engine == "r")
  {
    s <- do.call(what = paste("sde_", engine, sep =""),
                 list(x0 = region[1],
                      t = region[2],
                      drift = dynamics[[1]],
                      volat = dynamics[[2]],
                      n = n
                 )
    )
  } else if(engine == "cpp")
  {
    if(!is.numeric(dynamics))
    {
      stop("dynamics must be a vector of drift and volatility for cpp engine")
    }
    if(method != "em")
    {
      stop("Only Euler-Maruyama is available in cpp")
    }
    # Some day hope to implement RK2 for functions in cpp...so we'll leave this in
    s <- do.call(what = paste(method, "_ito", sep = ""), list(x0 = region[1],
                                                              t = region[2],
                                                              drift = dynamics[1],
                                                              volat = dynamics[2],
                                                              n = n)
                 )
  }
  return(s)
}

#' Solve arbitrary first order SDE systems
#'
#' @param f a list of drift functions with as many variables as IC, and in the same order with the names
#' @param g a list of volatility functions with as many variables as IC, and in the same order with the names
#' @param IC a list of variables with the same names as the non-time input as f, in the same order with initial values
#' @param t0 initial time
#' @param tn ending time
#' @param independent boolean to simulate processes driven by one Brownian motion (false) or as many independent Brownian motions as there are state variables
#' @param n number of sub-intervals in time-grid
#'
#' @description {A general Euler scheme implementation for arbitrary first-order SDE systems of finite dimension.}
#' @details {The list of functions must have syntax \code{f(t,x,y,z)} for each element matching the elements of \code{IC} as \code{list(x = y0, y = y0, z = z0)} for example.}
#' @return data.frame
#' @export samplePathSystem
samplePathSystem <- function(f, g, IC, t0 = 0, tn = 1, independent = FALSE, n = 1000)
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
      states[[i]][j] <- states[[i]][j-1]+h*do.call(what = f[[i]], args = input)+sqrt(h)*z[j-1]*do.call(what = g[[i]], args = input)
    }
  }
  states <- data.frame(time = tt, do.call(cbind, states))
  names(states) <- c("time", names(IC))
  return(states)
}


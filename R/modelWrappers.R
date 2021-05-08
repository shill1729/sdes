#' Numerically solve an SDE given a model/parameter set
#'
#' @param spot initial price
#' @param maturity time horizon to simulate over
#' @param model a list defining the model, see details
#' @param n number of variates to use
#' @param method numerical scheme either "em" or "rk2", only applicable to pure diffusions.
#'
#' @description {Numerically solve an SDE and generate a sample path under a variety of models:
#' "gbm", "mixture", "merton", "heston", for example. See details for further information.}
#' @details {The argument \code{model} must be a named list with elements \code{name} which is
#' eitehr "gbm", "mixture", "merton", "kou", "dkou", or "heston", and a vector
#' \code{param} containing the numeric parameters defining the aforementioned model.
#' \itemize{
#' \item for "gbm", \code{param} should be a vector of \eqn{(\mu, \sigma)}
#' \item for "mixture", \code{param} should be a \eqn{3 x d} matrix }
#' }
#' @return data.frame of time and price and possibly volatility for two-factor models
#' @export sdeSolve
sdeSolve <- function(spot, maturity, model, n, method = "em")
{
  param <- model$param
  s <- 0
  if(model$name == "gbm")
  {
    mu <- param[1]
    volat <- param[2]
    region <- c(0, maturity)
    dynamics <- list(drift = function(t, x) mu-0.5*volat^2,
                     diffusion = function(t, x) volat
    )
    s <- samplePathIto(region, dynamics, n, method)
    s$X <- spot*exp(s$X)
  } else if(model$name == "mixture")
  {
    mu <- function(t, x) findistr::mixture_drift(t, x, param)
    volat <- function(t, x) findistr::mixture_vol(t, x, param)
    region <- c(0, maturity)
    dynamics <- list(drift = function(t, x) mu(t,x)-0.5*volat(t,x)^2,
                     diffusion = function(t, x) volat(t,x)
    )
    s <- samplePathIto(region, dynamics, n, method)
    s$X <- spot*exp(s$X)
  } else if(model$name == "merton")
  {
    mu <- param[1]
    volat <- param[2]
    lambda <- param[3]
    jm <- param[4]
    jv <- param[5]
    region <- c(spot, maturity)
    dynamics <- list(drift = function(t, x) mu,
                     diffusion = function(t, x) volat
                     )
    jumps <- list(lambda = function(t, x) lambda,
                  distr =  "norm",
                  param = list(mean = jm, sd = jv)
                  )
    s <- samplePathItoLevy(region, dynamics, jumps, n)
  } else if(model$name == "heston")
  {
    # Assume param is in the form (mu, rho, kappa, theta, xi, v0)
    mu <- param[1]
    rho <- param[2]
    kappa <- param[3]
    theta <- param[4]
    xi <- param[5]
    v0 <- param[6]^2 # initial volatility
    IC <- list(V = v0, X = spot)
    drifts <- list(function(t, V, X) kappa*(theta-V),
                   function(t, V, X) mu*X
                   )
    diffusions <- list(function(t, V, X) xi*sqrt(V),
                       function(t, V, X) sqrt(V)*X
                       )
    s <- samplePathTwoStateSystem(drifts, diffusions, IC, 0, maturity, rho, n)
  }
  return(s)
}

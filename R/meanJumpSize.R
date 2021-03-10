#' Compute mean jump size given list of jump-dynamics variables
#'
#' @param jump_specs list of variables defining jump dynamics, see details
#' @param samples number of variates to use in Monte-Carlo estimation of mean jump-size
#'
#' @details the list \code{jump_specs} should contain
#' \itemize{
#' \item \code{distr} the jump-size distribution name: "kou", "norm" or "unif",
#' \item \code{param} the list of parameters for the above distribution}
#' @return numeric
#' @importFrom findistr rdkou rkou
#' @export meanJumpSize
meanJumpSize <- function(jump_specs, samples = 10000)
{
  if(!is.null(jump_specs))
  {
    # Extract jump-size distribution
    distr <- jump_specs$distr
    # Write function for calling rdistr(n)
    rdistr_name <- paste("r", distr, sep = "")
    rdistr <- function(y) do.call(what = rdistr_name, args = c(y, jump_specs$param))
    # Check if eta is passed and compute it if it is not
    if(is.null(jump_specs$eta))
    {

      # Error checking parameters for known distributions
      if(distr == "kou")
      {
        if(1/jump_specs$param$alpha < 1)
        {
          stop("mean positive jump must be < 1")
        }
        alpha <- 1/jump_specs$param$alpha
        beta <- 1/jump_specs$param$beta
        p <- jump_specs$param$p
        q <- 1-p
        eta <- q*beta/(beta+1)+p*alpha/(alpha-1)-1
      }
      if(distr == "unif")
      {
        if(jump_specs$param$min > jump_specs$param$max)
        {
          stop("min must be less than max for uniform distribution")
        }
        eta <- (exp(jump_specs$param$max)-exp(jump_specs$param$min))/(jump_specs$param$max-jump_specs$param$min)-1
      }
      if(distr == "norm")
      {
        if(jump_specs$param$sd <= 0)
        {
          stop("sd of jump-size must be positive")
        }
        eta <- exp(jump_specs$param$mean+0.5*jump_specs$param$sd^2)-1
      }
      if(distr == "dkou")
      {
        if(1/jump_specs$param$alpha < 1 || 1/jump_specs$param$beta < 1)
        {
          stop("mean jump sizes must be < 1")
        }
        if(jump_specs$param$p > 1 || jump_specs$param$p < 0)
        {
          stop("Probability of positive jump_specs must be in unit interval")
        }
        if(jump_specs$param$ku < 0 || jump_specs$param$kd > 0)
        {
          stop("Upward displacement must be positive and downward displacement must be negative")
        }
        p <- jump_specs$param$p
        q <- 1-p
        alpha <- 1/jump_specs$param$alpha
        beta <- 1/jump_specs$param$beta
        eta <- q*exp(jump_specs$param$kd)*beta/(beta+1)+p*exp(jump_specs$param$ku)*alpha/(alpha-1)-1
      }

      # If a known distr is not passed, use Monte-Carlo estimates
      if(!distr %in% c("kou", "unif", "norm", "dkou"))
      {
        eta <- mean(exp(rdistr(samples)))-1
      }

    } else
    {
      eta <- jump_specs$eta
    }
  } else
  {
    eta <- 0
  }
  return(eta)
}

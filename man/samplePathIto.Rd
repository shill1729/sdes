% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/samplePathWrappers.R
\name{samplePathIto}
\alias{samplePathIto}
\title{Simulate sample path of Ito-process by solving a SDE}
\usage{
samplePathIto(region, dynamics = NULL, n = 1000, method = "rk2")
}
\arguments{
\item{region}{initial point and length of time as a vector}

\item{dynamics}{list of drift and volatility functions of \code{(t, x)}, use \code{NULL}
for Brownian motion default dynamics}

\item{n}{number of variates}

\item{method}{either "em" for Euler-Maruyama scheme or rk2 for Runge-Kutta scheme of order 2.}
}
\value{
data.frame of time variable and space variable
}
\description{
{Simulate a sample path of a solution to an SDE.
Wrapper to internal solvers.}
}

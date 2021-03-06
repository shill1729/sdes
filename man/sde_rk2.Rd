% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/solvers.R
\name{sde_rk2}
\alias{sde_rk2}
\title{Runge-Kutta solver for Ito SDEs}
\usage{
sde_rk2(x0, t, drift, volat, n = 50)
}
\arguments{
\item{x0}{initial value of the process}

\item{t}{length to simulate over}

\item{drift}{the drift coefficient function of space-time}

\item{volat}{the volatility coefficient function of space-time}

\item{n}{number of sub-intervals in time-discretization}
}
\value{
data.frame
}
\description{
{The stochastic Runge-Kutta for SDEs. See wikipedia for mathematical details}
}

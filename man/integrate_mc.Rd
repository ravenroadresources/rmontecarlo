% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MC_functions.R
\name{integrate_mc}
\alias{integrate_mc}
\title{Aproximate Integral Solution by Monte Carlo Method}
\usage{
integrate_mc(func, interval, y_interval = NULL, n = 10000, plot = TRUE)
}
\arguments{
\item{func}{either a regular function of x, eithe a character string representing
the function to be integrated. Examples: \code{function(x) x + 2} OR \code{"x + 2"}}

\item{interval}{length 2 array, giving the x limits of the inegral.}

\item{y_interval}{length 2 array, giving the y limits. Default is NULL: in this case the function
calculate the minimum and maximum of the function whitin the \code{interval}.}

\item{n}{integer representing the number of repetition of the Monte Carlo method. Default = 1e4.}

\item{plot}{logical defining if the output should be plotted or not.}
}
\value{
a list fo 2 items: \code{$plot} which is the plot of the function and
    \code{$solution} with the value aproximating the exact solution of the integral.
}
\description{
Aproximate the solution of the integral of a function using Monte Carlo Method.
}
\examples{
\dontrun{
  ff <- "x * sin(x^2)"
  xx <- c(0, 2 * sqrt(pi))
  yy <- c(-10, 10)
  nn <- 1e5

  integrate_mc(ff, xx, n = nn)

  sol <- integrate_mc(function(x) x * sin(x^2), xx, y_interval = yy, n = nn)
  sol[[2]]
  sol$solution
}

}

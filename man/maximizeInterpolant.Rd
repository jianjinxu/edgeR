\name{maximizeInterpolant}
\alias{maximizeInterpolant}

\title{Maximize a function given a table of values by spline interpolation.}

\description{
Maximize a function given a table of values by spline interpolation.
}

\usage{
maximizeInterpolant(x, z, maxit=10, eps=1e-7, plot=FALSE)
}

\arguments{ 
\item{x}{numeric vector of the inputs of the function.}

\item{z}{numeric vector of the values of the function at the inputs given by \code{x}.}

\item{maxit}{numeric scalar giving the maximum number of iterations for the Newton-Raphson algorithm.}

\item{eps}{numeric scalar giving the convergence tolerance.}

\item{plot}{logical, whether or not to plot the function on those given points.}
}

\value{\code{maximizeInterpolant} returns a single value which maximizes the spline interpolation.
}

\details{
\code{maximizeInterpolant} calls the function \code{splinefun} to fit cubic spline interpolation given a set of points.

\code{maximizeInterpolant} uses Newton-Raphson algorithm in finding the maximum of the function performing the interpolation.
}

\author{Gordon Smyth}
\examples{
x <- seq(0,1,length=1000)
y <- rnorm(1000,1,1)
maximizeInterpolant(x,y)
}

\seealso{
\code{\link{splinefun}}
}

\keyword{interpolation}

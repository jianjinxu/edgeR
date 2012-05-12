\name{locallyWeightedMean}
\alias{locallyWeightedMean}

\title{Locally Weighted Mean}

\description{Smooth columns of matrix by non-robust loess curves of degree 0.}

\usage{
locallyWeightedMean(y, x, weights=rep(1,nrow(y)), span=0.5, cell=0.01)
}

\arguments{ 

\item{y}{numeric matrix of response variables.}

\item{x}{numeric vector, the covariate.}

\item{weights}{numeric vector of weights.}

\item{span}{width of the smoothing window, in terms of proportion of the data set.  Larger values produce smoother curves.}

\item{cell}{controls resolution of points used for the loess interpolation approximation, see \code{\link{loess.control}}.}
}

\value{matrix of fitted response.}

\details{
Given a matrix of response, \code{locallyWeightedMean} fits a loess curve with degree 0 for each column of the response and the covariate.
The smoothed column values are tricube-weighted means of the original values.

This function calls \code{stats:::simpleLoess}, the low-level version of \code{loess}.
}


\author{Yunshun Chen}

\seealso{
\code{\link{loess}}
}

\examples{
x <- 1:100
y <- matrix(rnorm(100*3), nrow=100, ncol=3)
head(y)
z <- locallyWeightedMean(y,x)
head(z)
}

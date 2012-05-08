\name{locallyWeightedMean}
\alias{locallyWeightedMean}

\title{Locally Weighted Mean}

\description{Apply faster tricube smoothing using simple loess with degree 0.}

\usage{
locallyWeightedMean(y, x, span=0.5, cell=0.01)
}

\arguments{ 

\item{y}{matrix of the response variable.}

\item{x}{vector of the covariate.}

\item{span}{The size of the neighbourhood that controls the degree of smoothing. Default is 0.5.}

\item{cell}{The control of the accuracy of the loess interpolation approximation. Default is 0.01}
}

\value{matrix of fitted response.}

\details{
Given a matrix of response, \code{locallyWeightedMean} fits a loess curve with degree 0 for each column of the response and the covariate.
}


\author{Yunshun Chen}

\examples{
x <- 1:1000
y <- matrix(rnorm(10000), nrow=1000, ncol=10)
l <- locallyWeightedMean(y,x)
}
\keyword{file}


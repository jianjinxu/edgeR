\name{thinCounts}
\alias{thinCounts}

\title{Binomial Thinning of Counts}

\description{
Reduce the size of Poisson-like counts by binomial thinning.
}

\usage{thinCounts(x, prob=0.5)}

\arguments{ 
	\item{x}{numeric vector or array of non-negative integers.}
	\item{prob}{numeric scalar or vector, the expected proportion of the counts to keep.}
}

\value{
A vector or array of the same dimensions as \code{x}, with thinned counts.
}

\details{
This function calls \code{rbinom} with \code{size=x} and \code{prob=prob} to generate the new counts.
}

\author{Gordon Smyth}

\examples{
x <- rpois(10,lambda=10)
thinCounts(x)
}

\keyword{models}

\name{cutWithMinN}
\alias{cutWithMinN}

\title{Cut numeric vector into intervals}

\description{
Divides the range of \code{x} into intervals, so that each interval contains a minimum number of values, and codes the values in \code{x} according to which interval they fall.
}

\usage{cutWithMinN(x, intervals=2, min.n=1)}

\arguments{ 
	\item{x}{numeric vector.}
	\item{intervals}{number of intervals (greater than or equal to 2).}
	\item{min.n}{minimum number of values in any interval.}
}

\value{
An integer vector of the same length as \code{x}, containing values from 1 to \code{intervals}.
}

\details{
This function calls \code{rbinom} with \code{size=x} and \code{prob=prob} to generate the new counts.
}

\author{Gordon Smyth}

\seealso{
\code{\link{cut}}, \code{\link{quantile}}.
}

\examples{
x <- c(1,2,3,4,5,6,7,100)
cutWithMinN(x,3,min.n=1)
}

\keyword{category}

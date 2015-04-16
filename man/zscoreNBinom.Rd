\name{zscoreNBinom}
\alias{zscoreNBinom}

\title{Z-score Equivalents of Negative Binomial Deviate}

\description{
Compute z-score equivalents of negative binomial random deviates.
}

\usage{
zscoreNBinom(q, size, mu) 
}

\arguments{
\item{q}{numeric vector or matrix giving negative binomial random values.}
\item{size}{negative binomial size parameter (>0).}
\item{mu}{mean of negative binomial distribution (>0).}
}

\value{
Numeric vector or matrix giving equivalent deviates from a standard normal distribution.
}

\details{
This function computes the mid-p value of \code{q}, then converts to the standard normal deviate with the same cumulative probability distribution value.

Care is taken to do the computations accurately in both tails of the distributions.
}

\author{Gordon Smyth}
\seealso{
\code{\link{pnbinom}}, \code{\link{qnorm}} in the stats package.
}
\examples{
zscoreNBinom(c(0,10,100), mu=10, size=1/10)
}
\keyword{distribution}

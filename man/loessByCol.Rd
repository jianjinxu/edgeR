\name{loessByCol}
\alias{loessByCol}

\title{Locally Weighted Mean By Column}

\description{Smooth columns of matrix by non-robust loess curves of degree 0.}

\usage{
loessByCol(y, x=NULL, weights=NULL, span=0.5, cell=0.01, method="loess")
}

\arguments{ 

\item{y}{numeric matrix to be smoothed.}

\item{x}{numeric covariate vector, defaults to equally spaced.}

\item{weights}{numeric vector of weights, defaults to all equal.}

\item{span}{the smoother span. This gives the proportion of points which influence the smooth at each value. Larger values give more smoothness.}

\item{cell}{controls resolution of points used for the loess interpolation approximation, see \code{\link{loess.control}}.}

\item{method}{character, \code{"loess"} or \code{"Rcode"}.}
}

\value{numeric matrix with smoothed columns.}

\details{
Given a matrix of response, \code{loessByCol} fits a loess curve with degree 0 to each column of the response matrix.
The smoothed column values are tricube-weighted means of the original values.

If \code{method="loess"}, this function calls \code{stats:::simpleLoess}, the low-level version of \code{loess}, for each column of \code{y}.
If \code{method="Rcode"}, then a simplified pure-R implementation is used.
}


\author{Yunshun Chen, Davis McCarthy, Gordon Smyth}

\seealso{
\code{\link{loess}}
}

\examples{
y <- matrix(rnorm(100*3), nrow=100, ncol=3)
head(y)
z <- loessByCol(y)
head(z)
}

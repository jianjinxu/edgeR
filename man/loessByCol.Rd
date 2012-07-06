\name{loessByCol}
\alias{loessByCol}

\title{Locally Weighted Mean By Column}

\description{Smooth columns of matrix by non-robust loess curves of degree 0.}

\usage{
loessByCol(y, x=NULL, span=0.5)
}

\arguments{ 

\item{y}{numeric matrix of response variables. Each row corresponds to an entry in the covariate vector.}

\item{x}{numeric covariate vector, defaults to equally spaced.}

\item{span}{width of the smoothing window, in terms of proportion of the data set. Larger values produce smoother curves.}
}

\value{numeric matrix with smoothed columns.}

\details{
Given a matrix of response, \code{loessByCol} fits a loess curve with degree 0 to each column of the response matrix.
The smoothed column values are tricube-weighted means of the original values.
}

\author{Aaron Lun, replacing earlier R code by Davis McCarthy}

\seealso{
\code{\link{loess}}
}

\examples{
y <- matrix(rnorm(100*3), nrow=100, ncol=3)
head(y)
out <- loessByCol(y)
head(out$fitted.values)
}

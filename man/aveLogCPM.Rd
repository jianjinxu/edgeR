\name{aveLogCPM}
\alias{aveLogCPM}

\title{Average Log Counts Per Million}

\description{
Compute genewise average log2 counts-per-million.
}

\usage{
aveLogCPM(y, dispersion=0.05, offset=0)
}

\arguments{
\item{y}{numeric matrix containing the negative binomial counts.  Rows for tags and columns for libraries.}

\item{dispersion}{numeric scalar or vector giving the dispersion parameter for each GLM.
Can be a scalar giving one value for all tags, or a vector of length equal to the number of tags giving tag-wise dispersions.}

\item{offset}{numeric vector or matrix giving the offset that is to be included in the log-linear model predictor.  Can be a scalar, a vector of length equal to the number of libraries, or a matrix of the same size as \code{y}.}
}

\details{
Uses \code{mglmOneGroup} to compute average cpm genewise for a set of libraries, and returns \code{log2(cpm+1)}.
}

\value{
Numeric vector.
}

\author{Gordon Smyth}

\examples{
## see help page for mglmOneGroup
}

\seealso{
\code{\link{mglmOneGroup}}
}

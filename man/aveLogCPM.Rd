\name{aveLogCPM}
\alias{aveLogCPM}

\title{Average Log Counts Per Million}

\description{
Compute average log2 counts-per-million for each row of counts.
}

\usage{
aveLogCPM(y, dispersion=0.05, offset=log(colSums(y)))
}

\arguments{
\item{y}{numeric matrix containing counts.  Rows for tags and columns for libraries.}

\item{dispersion}{numeric scalar or vector of negative-binomial dispersions.
Can be a scalar or a vector of length \code{nrow(y)}.}

\item{offset}{numeric vector or matrix giving the log(lib.size), where lib.size is the total library size that counts-per-million are computed relative to.}
}

\details{
Uses \code{mglmOneGroup} to compute average counts-per-million (cpm), as the intercept of a negative binomial glm, for each row of counts.
Then returns log2(cpm+1), the unit being added so that the resulting values are non-negative and infinite values are avoided.
}

\value{
Numeric vector giving log2(cpm+1) for each row of \code{y}. 
}

\author{Gordon Smyth}

\examples{
y <- matrix(c(0,100,30,40),2,2)
lib.size <- c(1000,10000)

# With disp=0, the function is equivalent to pooling the counts before dividing by lib.size:
aveLogCPM(y,offset=log(lib.size),dispersion=0)
cpms <- rowSums(y)/sum(lib.size)*1e6
log2( cpms + 1 )

# With disp large, the function is equivalent to row-wise averages of individual cpms:
aveLogCPM(y,offset=log(lib.size),dispersion=1e4)
cpms <- t(t(y)/lib.size)*1e6
log2( rowMeans( cpms +1 ))
}

\seealso{
See \code{\link{cpm}} for individual logCPM values, rather than tagwise averages.

The computations for \code{aveLogCPM} are done by \code{\link{mglmOneGroup}}.
}

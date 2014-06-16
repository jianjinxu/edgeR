\name{aveLogCPM}
\alias{aveLogCPM}
\alias{aveLogCPM.DGEList}
\alias{aveLogCPM.DGEGLM}
\alias{aveLogCPM.default}

\title{Average Log Counts Per Million}

\description{
Compute average log2 counts-per-million for each row of counts.
}

\usage{
\method{aveLogCPM}{DGEList}(y, normalized.lib.sizes=TRUE, prior.count=2, dispersion=NULL, \dots)
\method{aveLogCPM}{default}(y, lib.size=NULL, offset=NULL, prior.count=2, dispersion=NULL,
          weights=NULL, \dots)
}

\arguments{
\item{y}{numeric matrix containing counts.  Rows for tags and columns for libraries.}
\item{normalized.lib.sizes}{logical, use normalized library sizes?}
\item{prior.count}{average value to be added to each count, to avoid infinite values on the log-scale.}
\item{dispersion}{numeric scalar or vector of negative-binomial dispersions. Defaults to 0.05.}
\item{lib.size}{numeric vector of library sizes. Defaults to \code{colSums(y)}. Ignored if \code{offset} is not \code{NULL}.}
\item{offset}{numeric matrix of offsets for the log-linear models.}
\item{weights}{optional numeric matrix of observation weights.}
\item{\dots}{other arguments are not currently used.}
}

\details{
This function uses \code{mglmOneGroup} to compute average counts-per-million (AveCPM) for each row of counts, and returns log2(AveCPM).
An average value of \code{prior.count} is added to the counts before running \code{mglmOneGroup}.

This function is similar to

\code{rowMeans(cpm(y, log=TRUE, \dots))},

but with the refinement that larger library sizes are given more weight in the average.
The two versions will agree for large values of the dispersion.
}

\value{
Numeric vector giving log2(AveCPM) for each row of \code{y}. 
}

\author{Gordon Smyth}

\examples{
y <- matrix(c(0,100,30,40),2,2)
lib.size <- c(1000,10000)

# With disp large, the function is equivalent to row-wise averages of individual cpms:
aveLogCPM(y, dispersion=1e4)
cpm(y, log=TRUE, prior.count=2)

# With disp=0, the function is equivalent to pooling the counts before dividing by lib.size:
aveLogCPM(y,prior.count=0,dispersion=0)
cpms <- rowSums(y)/sum(lib.size)*1e6
log2(cpms)
}

\seealso{
See \code{\link{cpm}} for individual logCPM values, rather than tagwise averages.

The computations for \code{aveLogCPM} are done by \code{\link{mglmOneGroup}}.
}

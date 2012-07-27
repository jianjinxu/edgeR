\name{estimateCommonDisp}
\alias{estimateCommonDisp}

\title{Estimate Common Negative Binomial Dispersion by Conditional Maximum Likelihood}

\description{
Maximizes the negative binomial conditional common likelihood to give the estimate of the common dispersion across all tags.
}

\usage{
estimateCommonDisp(object, tol=1e-06, rowsum.filter=5, verbose=FALSE)
}

\arguments{
\item{object}{\code{DGEList} object}

\item{tol}{the desired accuracy, passed to \code{\link{optimize}}}

\item{rowsum.filter}{numeric scalar giving a value for the filtering out of low abundance tags in the estimation of the common dispersion.
Only tags with total sum of counts above this value are used in the estimation of the common dispersion.}

\item{verbose}{logical, if \code{TRUE} estimated dispersion and BCV will be printed to standard output.}
}

\value{Returns \code{object} with the following added components:
	\item{common.dispersion}{estimate of the common dispersion.}
	\item{pseudo.counts}{numeric matrix of quantile-quantile normalized counts.  These are counts adjusted so that the library sizes are equal, while preserving differences between groups and variability within each group.}
	\item{pseudo.lib.size}{the common library size to which the counts have been adjusted}
}

\details{
Implements the method of Robinson and Smyth (2008) for estimating a common dispersion parameter by conditional maximum likelihood.
The method of conditional maximum likelihood assumes that library sizes are equal, which is not true in general, so pseudocounts (counts adjusted so that the library sizes are equal) need to be calculated. The function \code{equalizeLibSizes} is called to adjust the counts using a quantile-to-quantile method, but this requires a fixed value for the common dispersion parameter. To obtain a good estimate for the common dispersion, pseudocounts are calculated under the Poisson model (dispersion is zero) and these pseudocounts are used to give an estimate of the common dispersion. This estimate of the common dispersion is then used to recalculate the pseudocounts, which are used to provide a final estimate of the common dispersion.
}

\references{
Robinson MD and Smyth GK (2008). Small-sample estimation of negative
binomial dispersion, with applications to SAGE data. \emph{Biostatistics},
9, 321-332
}

\author{Mark Robinson, Davis McCarthy, Gordon Smyth}
\examples{
# True dispersion is 1/5=0.2
y <- matrix(rnbinom(1000,mu=10,size=5),ncol=4)
d <- DGEList(counts=y,group=c(1,1,2,2),lib.size=c(1000:1003))
d <- estimateCommonDisp(d, verbose=TRUE)
}

\seealso{
\code{\link{equalizeLibSizes}}
}

\name{equalizeLibSizes}
\alias{equalizeLibSizes}

\title{Equalize Library Sizes by Quantile-to-Quantile Normalization}

\description{Adjusts counts so that the effective library sizes are equal, preserving fold-changes between groups and preserving biological variability within each group.}

\usage{
equalizeLibSizes(object, dispersion=NULL, common.lib.size)
}

\arguments{
\item{object}{\code{\link[edgeR:DGEList-class]{DGEList}} object}

\item{dispersion}{numeric vector of dispersion parameters.
By default, is extracted from \code{object} or, if \code{object} contains no dispersion information, is set to \code{0.05}.}

\item{common.lib.size}{numeric scalar, the library size to normalize to; default is the geometric mean of the original effective library sizes}
}

\value{A list with components
	\item{pseudo.counts}{numeric matrix of normalized pseudo-counts}
	\item{common.lib.size}{normalized library size}
}

\details{
Thus function implements the quantile-quantile normalization method of Robinson and Smyth (2008).
It computes normalized counts, or pseudo-counts, used by \code{exactTest} and \code{estimateCommonDisp}.

The output pseudo-counts are the counts that would have theoretically arisen had the effective library sizes been equal for all samples.
The pseudo-counts are computed in such as way as to preserve fold-change differences beween the groups defined by \code{object$samples$group} as well as biological variability within each group.
Consequently, the results will depend on how the groups are defined.

Note that the column sums of the \code{pseudo.counts} matrix will not generally be equal, because the effective library sizes are not necessarily the same as actual library sizes and because the normalized pseudo counts are not equal to expected counts.
}

\note{
This function is intended mainly for internal edgeR use.
It is not normally called directly by users.
}

\author{Mark Robinson, Davis McCarthy, Gordon Smyth}

\references{
Robinson MD and Smyth GK (2008).
Small-sample estimation of negative binomial dispersion, with applications to SAGE data.
\emph{Biostatistics}, 9, 321-332.
\url{http://biostatistics.oxfordjournals.org/content/9/2/321}
}

\seealso{
\code{\link{q2qnbinom}}
}

\examples{
ngenes <- 1000
nlibs <- 2
counts <- matrix(0,ngenes,nlibs)
colnames(counts) <- c("Sample1","Sample2")
counts[,1] <- rpois(ngenes,lambda=10)
counts[,2] <- rpois(ngenes,lambda=20)
summary(counts)
y <- DGEList(counts=counts)
out <- equalizeLibSizes(y)
summary(out$pseudo.counts)
}


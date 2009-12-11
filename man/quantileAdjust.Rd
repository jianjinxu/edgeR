\name{quantileAdjust}
\alias{quantileAdjust}

\title{Normalizes a Dataset by Using a Quantile Adjustment}

\description{The function adjusts (you might say normalizes) a dataset, creating pseudocounts that represents quantile-adjusted counts as if all samples had the same library size, while estimating the dispersion parameter.}


\usage{ 
quantileAdjust(object, N = prod(object$samples$lib.size)^(1/ncol(object$counts)), alpha = 0, null.hypothesis = FALSE, n.iter = 5, r.init = NULL, tol = 0.001, verbose=TRUE) 
}
\arguments{ 

\item{object}{list containing the raw counts with elements \code{counts} (table of counts), \code{group} (factor indicating group) and \code{lib.size} (numeric vector of library sizes)}

\item{N}{library size to normalize to; default is the geometric mean of the original library sizes}

\item{alpha}{weight to put on the individual tag's likelihood}

\item{null.hypothesis}{logical, whether to calculate the means and percentile under the null hypothesis; default is \code{FALSE}}

\item{n.iter}{number of iterations in estimating the size parameter}

\item{r.init}{initialized value of the size parameter; if \code{NULL}, then the common value on unadjusted counts is used}

\item{tol}{tolerance in estimating the size parameter}

\item{verbose}{whether to write comments, default \code{true}}
}

\value{list containing several elements used in downstream function calls.
	\item{r}{is the dispersion estimate}
	\item{pseudo}{is the quantile-adjusted pseudocounts}
	\item{ps}{is a list containing the abundance estimates}
	\item{N}{is the common library size}
	\item{p}{percentiles on which the quantile is based}
	\item{mu}{means on which the quantile is based}
}

\details{
An older function, no longer called by the functions recommended to carry out analysis of DGE data, namely \code{\link{estimateCommonDisp}}, \code{\link{estimateTagwiseDisp}} and \code{\link{exactTest}}. No longer recommended for use.
}


\author{Mark Robinson, Davis McCarthy}

\examples{
set.seed(0)
y<-matrix(rnbinom(40,size=1,mu=10),ncol=4)
d<-DGEList(counts=y,group=rep(1:2,each=2),lib.size=rep(c(1000:1001),2))
qA<-quantileAdjust(d,alpha=100)
}
\keyword{file}

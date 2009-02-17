\name{quantileAdjust}

\alias{quantileAdjust}

\title{Normalizes a dataset by using a quantile adjustment}

\description{The function adjusts (you might say normalizes) a dataset, creating pseudodata that represents quantile-adjusted data as if all samples had the same library size, while estimating the dispersion parameter.}


\usage{ 
quantileAdjust(object, N = prod(object$lib.size)^(1/ncol(object$data)), alpha = 0, null.hypothesis = FALSE, n.iter = 5, r.init = NULL, tol = 0.001, verbose=TRUE) 
}
\arguments{ 

\item{object}{list containing the raw data with elements \code{data} (table of counts), \code{group} (vector indicating group) and \code{lib.size} (vector of library sizes)}

\item{N}{library size to normalize to; default is the geometric mean of the original library sizes}

\item{alpha}{weight to put on the individual tag's likelihood}

\item{null.hypothesis}{logical, whether to calculate the means and percentile under the null hypothesis; default is \code{FALSE}}

\item{n.iter}{number of iterations in estimating the size parameter}

\item{r.init}{initialized value of the size parameter; if \code{NULL}, then the common value on unadjusted data is used}

\item{tol}{tolerance in estimating the size parameter}

\item{verbose}{whether to write comments, default \code{true}}
}

\value{list containing several elements used in downstream function calls.  \code{r} is the dispersion estimate, \code{pseudo} is the quantile-adjusted pseudodata, \code{ps} is a list containing the abundance estimates, \code{N} is the common library size and \code{p} and \code{mu} are the percentiles and means, respectively that the quantile is based on }

\author{Mark Robinson, Davis McCarthy}

\examples{
set.seed(0)
y<-matrix(rnbinom(40,size=1,mu=10),ncol=4)
d<-list(data=y,group=rep(1:2,each=2),lib.size=rep(c(1000:1001),2))
qA<-quantileAdjust(d,alpha=100)
}
\keyword{file}

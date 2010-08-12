\name{estimateDispIter}
\alias{estimateDispIter}

\title{Normalizes a Dataset Using Quantile Adjustment and Iteratively Estimates the Dispersion Parameter}

\description{The function equalizes the library sizes of a dataset (this could be understood as normalization), creating pseudocounts that represents quantile-adjusted counts as if all samples had the same library size, while using an iterative procedure to estimate the dispersion parameter.}


\usage{ 
estimateDispIter(object,N=exp(mean(log(object$samples$lib.size*object$samples$norm.factors))),prior.n=10,common.disp=FALSE, null.hypothesis=FALSE,n.iter=5,disp.init=NULL,tol=1e-6,verbose=TRUE) 

}
\arguments{ 

\item{object}{object of class \code{DGEList} containing (at least) the elements \code{counts} (table of raw counts), \code{group} (vector indicating group) and \code{lib.size} (vector of library sizes)}

\item{N}{numeric scalar giving the library size to which to normalize; default is the geometric mean of the original library sizes}

\item{prior.n}{numeric scalar; the smoothing parameter that indicates the weight to give to the common likelihood compared to the individual tag's likelihood; default value of \code{10} means that the common likelihood is given 10 times the weight of the individual tag/gene's likelihood in the estimation of the tag/genewise dispersion}

\item{common.disp}{logical, if \code{TRUE} then the common dispersion estimate is used for all tags/genes, otherwise tag/genewise dispersion parameters are estimated; default \code{FALSE}}

\item{null.hypothesis}{logical, whether to calculate the means and percentile under the null hypothesis; default is \code{FALSE}}

\item{n.iter}{number of iterations in estimating the dispersion parameter}

\item{disp.init}{numeric vector or scalar giving initialized value(s) of the dispersion parameter, can be a common value or tag/genewise values; if \code{NULL}, then the common value on Poisson-adjusted counts is used}

\item{tol}{numeric scalar, tolerance in estimating the dispersion parameter}

\item{verbose}{logical, whether to write comments, default \code{TRUE}}
}

\value{list containing the following elements.
	\item{dispersion}{numeric vector giving the estimate of the dispersion parameter for each tag/gene}
	\item{pseudo}{numeric matrix of quantile-adjusted pseudocounts}
	\item{conc}{list containing the estimates of the concentration of each tag in the underlying sample; \code{conc$p.common} gives estimates under the null hypothesis of no difference between groups; \code{conc$p.group} gives the estimate of the concentration for each tag within each group; concentration is a measure of abundance and thus expression level for the tags}
	\item{N}{numeric scalar, the common library size to which the counts have been adjusted}
	\item{mu}{numeric matrix of means that the quantile adjustment is based on}
}

\author{Mark Robinson, Davis McCarthy}

\examples{
set.seed(0)
y<-matrix(rnbinom(40,size=1,mu=10),ncol=4)
d<-DGEList(counts=y,group=rep(1:2,each=2),lib.size=rep(c(1000:1001),2))
disp.out<-estimateDispIter(d,prior.n=10)
}

\seealso{
The use of \code{\link{estimateCommonDisp}} and \code{\link{estimateTagwiseDisp}} are preferred for the calculation of the common dispersion and tagwise dispersion estimates, respectively.
}

\keyword{file}

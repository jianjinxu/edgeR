\name{estimateDispIter}
\alias{estimateDispIter}

\title{Normalizes a dataset using quantile adjustment and iteratively estimates the dispersion parameter}

\description{The function equalizes the library sizes of a dataset (this could be understood as normalization), creating pseudocounts that represents quantile-adjusted counts as if all samples had the same library size, while using an iterative procedure to estimate the dispersion parameter.}


\usage{ 
estimateDispIter(object,N=prod(object$samples$lib.size)^(1/ncol(object$counts)),prior.n=10,common.disp=FALSE, null.hypothesis=FALSE,n.iter=5,disp.init=NULL,tol=1e-6,verbose=TRUE) 

}
\arguments{ 

\item{object}{list containing the raw counts with elements \code{counts} (table of counts), \code{group} (vector indicating group) and \code{lib.size} (vector of library sizes)}

\item{N}{library size to normalize to; default is the geometric mean of the original library sizes}

\item{prior.n}{smoothing paramter that indicates the weight to put on the common likelihood compared to the individual tag's likelihood; default \code{10} means that the common likelihood is given 10 times the weight of the individual tag/gene's likelihood in the estimation of the tag/genewise dispersion}

\item{common.disp}{logical, if \code{TRUE} then the common dispersion estimate is used for all tags/genes, otherwise tag/genewise dispersion parameters are estimated; default \code{FALSE}}

\item{null.hypothesis}{logical, whether to calculate the means and percentile under the null hypothesis; default is \code{FALSE}}

\item{n.iter}{number of iterations in estimating the size parameter}

\item{disp.init}{initialized value(s) of the dispersion parameter, can be a common value or tag/genewise values; if \code{NULL}, then the common value on Poisson-adjusted counts is used}

\item{tol}{tolerance in estimating the size parameter}

\item{verbose}{whether to write comments, default \code{TRUE}}
}

\value{list containing several elements used in downstream function calls.  \code{dispersion} is the dispersion estimate, \code{pseudo} is the quantile-adjusted pseudocounts, \code{conc} is a list containing the abundance estimates (i.e. concentrations of the tags/genes in the RNA sample), \code{N} is the common library size and \code{mu} is the means that the quantile adjustment is based on }

\author{Mark Robinson, Davis McCarthy}

\examples{
set.seed(0)
y<-matrix(rnbinom(40,size=1,mu=10),ncol=4)
d<-DGEList(counts=y,group=rep(1:2,each=2),lib.size=rep(c(1000:1001),2))
disp.out<-estimateDispIter(d,prior.n=10)
}
\keyword{file}

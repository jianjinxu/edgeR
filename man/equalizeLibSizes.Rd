\name{equalizeLibSizes}
\alias{equalizeLibSizes}

\title{Quantile adjustment to equalize library sizes for a fixed value of the dispersion parameter}

\description{A function that uses a NB quantile-to-quantile method to adjust the libraries of counts so that library sizes are equal for a fixed value of the dispersion parameter. This is the function provides the necessary framework and calculations to call \code{q2qnbinom} that actually generates the pseudodata required by quantileAdjust to adjust (normalise) the library sizes and estimate the dispersion parameter. This function calls \code{estimatePs} to estimate the expression proportion for each tag, which is needed to calculate the \code{input.mean} and \code{output.mean} for each tag, which are passed to \code{q2qnbinom} along with the unadjusted counts and the fixed value(s) for the dispersion parameter. }


\usage{ 
equalizeLibSizes(object, disp=0, N=prod(object$samples$lib.size)^(1/ncol(object$counts)),null.hypothesis=FALSE)
}

\arguments{ 

\item{object}{\code{DGEList} object containing the raw counts with elements \code{counts} (table of counts), \code{group} (vector indicating group) and \code{lib.size} (vector of library sizes)}

\item{disp}{scalar, vector or matrix of \code{dispersion} parameters}

\item{N}{library size to normalize to; default is the geometric mean of the original library sizes}

\item{null.hypothesis}{logical, whether to calculate the \code{input.mean} and \code{output.mean} under the null hypothesis; default is \code{FALSE}}
}

\value{list with elements 
	\item{pseudo}{pseudocounts, i.e. adjusted counts for equalized libraries}
	\item{mu}{the estimated means for tags/genes for equalized libraries}
	\item{conc}{estimated tag/gene concentrations (measure of expression level); output from \code{estimatePs}}
	\item{N}{normalized library size}
}

\author{Mark Robinson, Davis McCarthy}

\examples{
y<-matrix(rnbinom(10000,size=2,mu=10),ncol=4)
d<-DGEList(counts=y,group=rep(1:2,each=2),lib.size=rep(c(1000,1010),2))
ps<-estimatePs(d,r=2)
N<-prod(d$lib.size)^(1/ncol(d$counts))
q2q.out<-equalizeLibSizes(d,disp=0.5,null.hypothesis=FALSE)
}
\keyword{algebra}

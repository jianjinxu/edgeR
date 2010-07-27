\name{equalizeLibSizes}
\alias{equalizeLibSizes}

\title{Quantile Adjustment to Equalize Library Sizes for a Fixed Value of the Dispersion Parameter}

\description{A function that uses a NB quantile-to-quantile method to adjust the libraries of counts so that library sizes are equal for a fixed value of the dispersion parameter.}


\usage{ 
equalizeLibSizes(object, disp=0, N=exp(mean(log(object$samples$lib.size))),null.hypothesis=FALSE)
}

\arguments{ 

\item{object}{\code{DGEList} object containing the raw counts with elements \code{counts} (table of counts), \code{group} (vector indicating group) and \code{lib.size} (vector of library sizes)}

\item{disp}{numeric scalar or vector of \code{dispersion} parameters; if a scalar, then a common dispersion parameter is used for all tags}

\item{N}{numeric scalar, the library size to normalize to; default is the geometric mean of the original library sizes}

\item{null.hypothesis}{logical, whether to calculate the \code{input.mean} and \code{output.mean} under the null hypothesis; default is \code{FALSE}}
}

\value{A list with elements 
	\item{pseudo}{numeric matrix of pseudocounts, i.e. adjusted counts for equalized libraries}
	\item{conc}{list with elements \code{conc.common} (vector giving overall proportion/concentration for each tag), and \code{conc.group} (matrix with columns giving estimates of tag/gene concentrations (proportion of total RNA for that group that that particular tag/gene contributes) for different groups); output from \code{estimatePs}}
	\item{N}{normalized library size}
}

\details{
The function \code{equalizeLibSizes} provides the necessary framework and calculations to call \code{q2qnbinom}, for given value(s) of the dispersion parameter. The function \code{q2qnbinom} actually generates the pseudocounts, the counts that have been adjusted for normalized library sizes. These pseudocounts are required to estimate the dispersion parameter, as the methods used by \code{\link{estimateCommonDisp}} and \code{\link{estimateTagwiseDisp}} rely on the assumption of equal library sizes. This function calls \code{estimatePs} to estimate the expression proportion for each tag, which is needed to calculate the \code{input.mean} and \code{output.mean} for each tag, which are passed to \code{q2qnbinom} along with the unadjusted counts and the fixed value(s) for the dispersion parameter. 
}

\author{Mark Robinson, Davis McCarthy}

\examples{
y<-matrix(rnbinom(10000,size=2,mu=10),ncol=4)
d<-DGEList(counts=y,group=rep(1:2,each=2),lib.size=rep(c(1000,1010),2))
ps<-estimatePs(d,r=2)
q2q.out<-equalizeLibSizes(d,disp=0.5,null.hypothesis=FALSE)
}
\keyword{algebra}

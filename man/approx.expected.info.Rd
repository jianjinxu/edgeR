\name{approx.expected.info}
\alias{approx.expected.info}

\title{Approximate Expected Information (Fisher Information)}

\description{Using a linear fit (for simplicity), the expected information from the conditional log likelihood of the dispersion parameter of the negative binomial is calculated over all genes.}


\usage{
approx.expected.info(object, d, pseudo, robust = FALSE) 
}
\arguments{ 

\item{object}{\code{DGEList} object containing the raw counts with (at least) elements \code{counts} (table of counts), \code{group} (vector indicating group) and \code{lib.size} (vector of library sizes)}

\item{d}{numeric vector giving the delta parameter for negative binomial - \code{ phi/(phi+1) }; either of length 1 or  of length equal to the number of tags/transcripts (i.e. number of rows of \code{object$counts}. }

\item{pseudo}{numeric matrix of pseudocounts from output of \code{estimateDispIter}}

\item{robust}{logical on whether to use a robust fit, default \code{FALSE}}

}

\value{numeric vector of approximate values of the Fisher information for each tag/transcript (with length same as the number of rows of the original counts)  }

\author{Mark Robinson}

\examples{
set.seed(0)
y<-matrix(rnbinom(40,size=1,mu=10),ncol=4)
d<-DGEList(counts=y,group=rep(1:2,each=2),lib.size=rep(c(1000:1001),2))
d<-estimateCommonDisp(d)
d<-estimateTagwiseDisp(d,prior.n=10)
exp.inf<-approx.expected.info(d,1/(1 + d$common.dispersion),d$pseudo.alt)
}

\seealso{
This function is used in the algorithm for estimating an appropriate amount of smoothing for the dispersion estimates carried out by \code{\link{estimateSmoothing}}.
}

\keyword{file}

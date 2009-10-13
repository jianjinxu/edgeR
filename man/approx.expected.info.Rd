\name{approx.expected.info}
\alias{approx.expected.info}

\title{Approximate of expected information (Fisher information)}

\description{Using a linear fit (for simplicity), the expected information from the conditional log likelihood of the dispersion parameter of the negative binomial is calculated over all genes.}


\usage{ 
approx.expected.info(object, d, pseudo, robust = FALSE) 
}
\arguments{ 

\item{object}{\code{DGEList} object containing the raw counts with elements \code{counts} (table of counts), \code{group} (vector indicating group) and \code{lib.size} (vector of library sizes)}

\item{d}{delta parameter for negative binomial - \code{ phi/(phi+1) } }

\item{pseudo}{matrix of pseudocounts from output of \code{estimateDispIter}}

\item{robust}{logical on whether to use a robust fit, default \code{FALSE}}

}

\value{ vector of Fisher information approximates (with length same as the number of rows of the original counts)  }

\author{Mark Robinson}

\examples{
set.seed(0)
y<-matrix(rnbinom(40,size=1,mu=10),ncol=4)
d<-DGEList(counts=y,group=rep(1:2,each=2),lib.size=rep(c(1000:1001),2))
qA<-estimateDispIter(d,prior.n=10)
exp.inf<-approx.expected.info(d,1/(1 + qA$dispersion[1]),qA$pseudo)
}
\keyword{file}

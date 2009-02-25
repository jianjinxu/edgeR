\name{approx.expected.info}

\alias{approx.expected.info}

\title{Approximate of expected information (Fisher information)}

\description{Using a linear fit (for simplicity), the expected information from the conditional log likelihood of the dispersion parameter of the negative binomial is calculated over all genes.}


\usage{ 
approx.expected.info(object, d, qA, robust = FALSE) 
}
\arguments{ 

\item{object}{\code{DGEList} object containing the raw data with elements \code{data} (table of counts), \code{group} (vector indicating group) and \code{lib.size} (vector of library sizes)}

\item{d}{delta parameter for negative binomial - \code{ phi/(phi+1) } }

\item{qA}{list from output of \code{quantileAdjust}}

\item{robust}{logical on whether to use a robust fit, default \code{FALSE}}

}

\value{ vector of Fisher information approximates (with length same as the number of rows of the original data)  }

\author{Mark Robinson}

\examples{
set.seed(0)
y<-matrix(rnbinom(40,size=1,mu=10),ncol=4)
d<-DGEList(data=y,group=rep(1:2,each=2),lib.size=rep(c(1000:1001),2))
qA<-quantileAdjust(d,alpha=100)
exp.inf<-approx.expected.info(d,1/(1 + qA$r[1]),qA)
}
\keyword{file}

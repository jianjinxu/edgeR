\name{alpha.approxeb}

\alias{alpha.approxeb}

\title{Estimate the prior weight, alpha}

\description{Estimate the prior weight, using an approximate empirical Bayes rule}


\usage{ 
alpha.approxeb(object, verbose=TRUE) 
}
\arguments{ 
\item{object}{\code{DGEList} object containing the raw data with elements \code{data} (table of counts), \code{group} (vector indicating group) and \code{lib.size} (vector of library sizes)}

\item{verbose}{whether to write comments, default \code{true}}
}

\value{ \code{EBList} object with elements \code{sigma2.0.est} (numeric scale \code{sigma_0^2} estimate), \code{alpha} (estimate for the prior weight, alpha), \code{scores} (likelihood scores), \code{infos} (Fisher expected information), \code{quantileAdjusted} (list from output of \code{quantileAdjust})}

\author{Mark Robinson, Davis McCarthy}

\examples{
y<-matrix(rnbinom(20,size=1,mu=10),nrow=5)
d<-DGEList(data=y,group=rep(1:2,each=2),lib.size=rep(c(1000:1001),2))
alpha<-alpha.approxeb(d)
}
\keyword{file}

\name{alpha.approxeb}

\alias{alpha.approxeb}

\title{Estimate the Prior Weight, Alpha}

\description{Estimate the prior weight, using an approximate empirical Bayes rule}


\usage{ 
alpha.approxeb(object, verbose=TRUE) 
}
\arguments{ 
\item{object}{\code{DGEList} object containing the raw counts with elements \code{counts} (table of counts), \code{group} (vector indicating group) and \code{lib.size} (vector of library sizes)}

\item{verbose}{whether to write comments, default \code{true}}
}

\value{ \code{EBList} object with elements \code{sigma2.0.est} (numeric scale \code{sigma_0^2} estimate), \code{alpha} (estimate for the prior weight, alpha), \code{scores} (likelihood scores), \code{infos} (Fisher expected information), \code{quantileAdjusted} (list from output of \code{quantileAdjust})}

\details{
An older function, no longer called by the functions recommended to carry out analysis of DGE data, namely \code{\link{estimateCommonDisp}}, \code{\link{estimateTagwiseDisp}} and \code{\link{exactTest}}. Estimation of the prior weight should now be done using \code{\link{estimateSmoothing}}, if at all.
}

\author{Mark Robinson, Davis McCarthy}

\examples{
y<-matrix(rnbinom(20,size=1,mu=10),nrow=5)
d<-DGEList(counts=y,group=rep(1:2,each=2),lib.size=rep(c(1000:1001),2))
alpha<-alpha.approxeb(d)
}
\keyword{file}

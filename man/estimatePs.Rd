\name{estimatePs}

\alias{estimatePs}

\title{Estimate expression proportions}

\description{Estimate expression proportions (maximum likelihood with size fixed) based on negative binomial for each tag and sample group (only 2 groups implemented at this point)}

\usage{ 
estimatePs(y1, y2, lib.size1, lib.size2, r, tol = 1e-10, maxit = 30) 
}
\arguments{ 

\item{y1}{matrix of counts for first group}

\item{y2}{matrix of counts for second group}

\item{lib.size1}{vector of library sizes for first group}

\item{lib.size2}{vector of library sizes for second group}

\item{r}{size parameter of negative binomial}

\item{tol}{tolerance between iterations}

\item{maxit}{maximum number of iterations}
}

\value{ list with elements \code{p} (overall proportion), \code{p1} (estimates for first group), \code{p2} (estimates for second group)}

\author{Mark Robinson}

\examples{
y1<-matrix(rnbinom(10,size=1,mu=10),nrow=5)
y2<-matrix(rnbinom(10,size=1,mu=5),nrow=5)
ps<-estimatePs(y1,y2,c(1000,1001),c(1000,1001),r=1)
}
\keyword{file}

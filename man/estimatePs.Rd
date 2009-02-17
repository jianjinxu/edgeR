\name{estimatePs}

\alias{estimatePs}

\title{Estimate expression proportions}

\description{Estimate expression proportions (maximum likelihood with size fixed) based on negative binomial for each tag and sample group (only 2 groups implemented at this point)}

\usage{ 
estimatePs(object, r, tol = 1e-10, maxit = 30) 
}
\arguments{ 

\item{object}{list containing the raw data with elements \code{data} (table of counts), \code{group} (vector indicating group) and \code{lib.size} (vector of library sizes)}

\item{r}{size parameter of negative binomial}

\item{tol}{tolerance between iterations}

\item{maxit}{maximum number of iterations}
}

\value{ list with elements \code{p.common} (vector giving overall proportion for each tag), \code{p.group} (matrix with columns giving estimates of proportions for different groups)}

\author{Mark Robinson, Davis McCarthy}

\examples{
set.seed(0)
y<-matrix(rnbinom(40,size=1,mu=10),ncol=4)
d<-list(data=y,group=rep(1:2,each=2),lib.size=rep(c(1000:1001),2))
ps<-estimatePs(d,r=1)
}
\keyword{file}

\name{commonCondLogLikDerDelta}
\alias{commonCondLogLikDerDelta}

\title{Common conditional log-likelihood in terms of delta}

\description{Common conditional log-likelihood parameterized in terms of delta (\code{phi / (phi+1)}) }

\usage{ 
commonCondLogLikDerDelta(y, delta, der = 0, doSum = FALSE) 
}
\arguments{ 

\item{y}{list with elements comprising the matrices of count data (or pseudocounts) for the different groups}

\item{delta}{delta (\code{phi / (phi+1)}) parameter of negative binomial}

\item{der}{derivative, either 0 (the function), 1 (first derivative) or 2 (second derivative)}

\item{doSum}{logical, whether to sum over samples or not (default \code{FALSE}}

}

\value{ numeric scalar of function/derivative evaluated at given delta}

\author{Davis McCarthy}

\examples{
counts<-matrix(rnbinom(20,size=1,mu=10),nrow=5)
d<-DGEList(counts=counts,group=rep(1:2,each=2),lib.size=rep(c(1000:1001),2))
y<-splitIntoGroups(d)
ll1<-commonCondLogLikDerDelta(y,delta=0.5,der=0,doSum=FALSE)
ll2<-commonCondLogLikDerDelta(y,delta=0.5,der=1)
}
\keyword{file}
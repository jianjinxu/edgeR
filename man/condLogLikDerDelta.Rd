\name{condLogLikDerDelta}

\alias{condLogLikDerDelta}

\title{Conditional log-likelihood in terms of delta}

\description{Conditional log-likelihood parameterized in terms of delta (\code{phi / (phi+1)}) }

\usage{ 
condLogLikDerDelta(y, delta, grid = TRUE, der = 1, doSum = TRUE) 
}
\arguments{ 

\item{y}{matrix with count data (or pseudodata)}

\item{delta}{delta (\code{phi / (phi+1)})parameter of negative binomial}

\item{grid}{logical, whether to calculate a grid over the values of delta}

\item{der}{derivative, either 0 (the function), 1 (first derivative) or 2 (second derivative)}

\item{doSum}{logical, whether to sum over samples or not (default \code{TRUE}}

}

\value{ vector or matrix of function/derivative evaluations}

\author{Mark Robinson, Davis McCarthy}

\examples{
y1<-matrix(rnbinom(10,size=1,mu=10),nrow=5)
v1<-seq(.1,.9,length=9)
ll1<-condLogLikDerDelta(y1,v1,grid=TRUE,der=0,doSum=FALSE)
ll2<-condLogLikDerDelta(y1,delta=.5,grid=FALSE,der=0)
}
\keyword{file}

\name{condLogLikDerSize}

\alias{condLogLikDerSize}

\title{Conditional log-likelihood in terms of size}

\description{Conditional log-likelihood paramterized in terms of size (\code{1 / phi}) }

\usage{ 
condLogLikDerSize(y, r, der=1)
}
\arguments{ 

\item{y}{matrix with count data (or pseudodata)}

\item{r}{size parameter of negative binomial distribution}

\item{der}{derivative, either 0 (the function), 1 (first derivative) or 2 (second derivative)}

}

\value{ vector or matrix of function/derivative evaluations}

\author{Mark Robinson, Davis McCarthy}

\examples{
y1<-matrix(rnbinom(10,size=1,mu=10),nrow=5)
ll2<-condLogLikDerSize(y1,r=10,der=0)
}
\keyword{file}
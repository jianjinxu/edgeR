\name{logLikDerP}

\alias{logLikDerP}

\title{Log-likelihood for proportion}

\description{Log-likelihood and derivatives for the proportion parameter of negative binomial (mean = library size * proportion)}

\usage{
logLikDerP(p, y, lib.size, r, der = 0)
}

\arguments{

\item{p}{vector of proportion parameters to be evaluated}

\item{y}{matrix of data}

\item{lib.size}{vector of library sizes}

\item{r}{size parameter of negative binomial distribution}

\item{der}{derivative, either 0 (the function), 1 (first derivative) or 2 (second derivative)}

}

\value{ vector of evaluations}


\author{Mark Robinson, Davis McCarthy}

\examples{
y<-matrix(rnbinom(20,size=1.5,mu=10),nrow=5)
d<-DGEList(data=y,group=rep(1:2,each=2),lib.size=rep(c(1000:1001),2))

this.p<-rowMeans( y/ outer(rep(1,nrow(y)),d$lib.size) )
d1p<-logLikDerP(this.p,y,d$lib.size,r=1.5,der=1)
}
\keyword{file}

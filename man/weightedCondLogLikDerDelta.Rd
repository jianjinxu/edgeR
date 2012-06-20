\name{weightedCondLogLikDerDelta}
\alias{weightedCondLogLikDerDelta}

\title{Weighted Conditional Log-Likelihood in Terms of Delta}

\description{Weighted conditional log-likelihood parameterized in terms of delta (\code{phi / (phi+1)}) for a given tag/gene - maximized to find the smoothed (moderated) estimate of the dispersion parameter}

\usage{
weightedCondLogLikDerDelta(y, delta, tag, prior.n=10, ntags=nrow(y[[1]]), der=0)
}

\arguments{
\item{y}{list with elements comprising the matrices of count data (or pseudocounts) for the different groups}

\item{delta}{delta (\code{phi / (phi+1)})parameter of negative binomial}

\item{tag}{tag/gene at which the weighted conditional log-likelihood is evaluated}

\item{prior.n}{smoothing paramter that indicates the weight to put on the common likelihood compared to the individual tag's likelihood; default \code{10} means that the common likelihood is given 10 times the weight of the individual tag/gene's likelihood in the estimation of the tag/genewise dispersion}

\item{ntags}{numeric scalar number of tags/genes in the dataset to be analysed}

\item{der}{derivative, either 0 (the function), 1 (first derivative) or 2 (second derivative)}
}

\value{ numeric scalar of function/derivative evaluated for the  given tag/gene and delta}

\details{
This function computes the weighted conditional log-likelihood for a given tag, parameterized in terms of delta. The value of delta that maximizes the weighted conditional log-likelihood is converted back to the \code{phi} scale, and this value is the estimate of the smoothed (moderated) dispersion parameter for that particular tag. The delta scale for convenience (delta is bounded between 0 and 1). 
}


\author{Mark Robinson, Davis McCarthy}

\examples{
counts<-matrix(rnbinom(20,size=1,mu=10),nrow=5)
d<-DGEList(counts=counts,group=rep(1:2,each=2),lib.size=rep(c(1000:1001),2))
y<-splitIntoGroups(d)
ll1<-weightedCondLogLikDerDelta(y,delta=0.5,tag=1,prior.n=10,der=0)
ll2<-weightedCondLogLikDerDelta(y,delta=0.5,tag=1,prior.n=10,der=1)
}
\keyword{file}

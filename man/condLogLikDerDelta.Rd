\name{condLogLikDerDelta}
\alias{condLogLikDerDelta}

\title{Conditional Log-Likelihood in Terms of Delta}

\description{Conditional negative binomial log-likelihood parameterized in terms of delta (\code{phi / (phi+1)}) }

\usage{
condLogLikDerDelta(y, delta, der = 1)
}
\arguments{ 
\item{y}{matrix with count data (or pseudocounts)}

\item{delta}{delta (\code{phi / (phi+1)})parameter of negative binomial}

\item{der}{derivative, either 0 (the function), 1 (first derivative) or 2 (second derivative)}
}

\value{ vector or matrix of function/derivative evaluations}

\details{
This function computes the individual tag conditional log-likelihood for each tag.
It is necessary for computing both the common conditional log-likelihood and the weighted conditional log-likelihood, which are used to find the common and tagwise (moderated) estimates of the dispersion parameter.
The delta scale for convenience (delta is bounded between 0 and 1). 
}

\author{Mark Robinson, Davis McCarthy}

\examples{
y1<-matrix(rnbinom(10,size=1,mu=10),nrow=5)
ll2<-condLogLikDerDelta(y1,delta=0.5,der=0)
}

\seealso{
\code{\link{commonCondLogLikDerDelta}} and \code{\link{weightedCondLogLikDerDelta}} rely on \code{condLogLikDerDelta}, and at a user level, \code{\link{estimateCommonDisp}} and \code{\link{estimateTagwiseDisp}} are used to estimate the common and (moderated) tagwise dispersion estimates, respectively. \code{condLogLikDerDelta} calls \code{condLogLikDerSize}, the function that does the mathematical calculations.
}

\keyword{file}

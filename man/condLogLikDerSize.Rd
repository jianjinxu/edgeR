\name{condLogLikDerSize}
\alias{condLogLikDerSize}

\title{Log-Likelihood of the Common Dispersion for a Single Equalized Group}

\description{Derivatives of the conditional negative-binomial log-likelihood (for each tag/transcript) with respect to the common dispersion parameter, for a single group of replicate libraries of the same size.
Parameterized in terms of size or precision (\code{1/phi}).}

\usage{condLogLikDerSize(y, r, der=1)}

\arguments{ 
\item{y}{matrix of (pseudo) count data}
\item{r}{size parameter of negative binomial distribution}
\item{der}{order of derivative required, either 0 (the function), 1 (first derivative) or 2 (second derivative)}
}

\value{vector of function/derivative evaluations, one for each transcript}

\details{The library sizes must be equalized before running this function. This function carries out the actual mathematical computations for the conditional log-likelihood and its derivatives, calculating the conditional log-likelihood for each tag/transcript.
}


\author{Mark Robinson, Davis McCarthy}

\examples{
y <- matrix(rnbinom(10,size=1,mu=10),nrow=5)
condLogLikDerSize(y,r=1,der=1)
}
\keyword{file}
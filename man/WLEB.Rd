\name{WLEB}
\alias{WLEB}

\title{Calculate Weighted Likelihood Empirical Bayes Estimates}

\description{
Estimates the parameters which maximize the given log-likelihood matrix using empirical Bayes method.
}

\usage{
WLEB(theta, loglik, prior.n, covariate, trend.method="locfit", span=NULL, overall=TRUE, trend=TRUE, individual=TRUE, m0=NULL, m0.out=FALSE)
}

\arguments{
\item{theta}{numeric vector of values of the parameter at which the log-likelihoods are calculated.}

\item{loglik}{numeric matrix of log-likelihood of all the candidates at those values of parameter.}

\item{prior.n}{numeric scaler, estimate of the prior weight, i.e. the smoothing parameter that indicates the weight to put on the common likelihood compared to the individual's likelihood.}

\item{covariate}{numeric vector of values across which a parameter trend is fitted}

\item{trend.method}{method for estimating the parameter trend. Possible values are \code{"none"}, \code{"movingave"} and \code{"loess"}.}

\item{span}{width of the smoothing window, as a proportion of the data set.}

\item{overall}{logical, should a single value of the parameter which maximizes the sum of all the log-likelihoods be estimated?}

\item{trend}{logical, should a parameter trend (against the covariate) which maximizes the local shared log-likelihoods be estimated?}

\item{individual}{logical, should individual estimates of all the candidates after applying empirical Bayes method along the trend be estimated?}

\item{m0}{numeric matrix of local shared log-likelihoods. If \code{Null}, they will be calculated using the method selected by \code{trend.method}.}

\item{m0.out}{logical, should local shared log-likelihoods be included in the output?}
}

\details{
This function is a generic function that calculates an overall estimate, trend estimates and individual estimates for each candidate given the values of the log-likelihood of all the candidates at some specified parameter values.
}

\value{
A list with the following:
	\item{overall}{the parameter estimate that maximizes the sum of all the log-likelihoods.}
	\item{trend}{the estimated trended parameters against the covariate.}
	\item{individual}{the individual estimates of all the candidates after applying empirical Bayes method along the trend.}
	\item{shared.loglik}{the estimated numeric matrix of local shared log-likelihoods}
}

\author{Yunshun Chen, Gordon Smyth}
\examples{
y <- matrix(rpois(100, lambda=10), ncol=4)
theta <- 7:14
loglik <- matrix(0,nrow=nrow(y),ncol=length(theta))
for(i in 1:nrow(y))
for(j in 1:length(theta))
	loglik[i,j] <- sum(dpois(y[i,], theta[j] ,log=TRUE))
covariate <- log(rowSums(y))
out <- WLEB(theta, loglik, prior.n=3, covariate)
out
}

\seealso{
\code{\link{locfitByCol}}, \code{\link{movingAverageByCol}} and \code{\link{loessByCol}} implement the local fit, moving average or loess smoothers.
}

\keyword{algebra}

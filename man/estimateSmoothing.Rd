\name{estimateSmoothing}
\alias{estimateSmoothing}

\title{Estimate the Prior Weight}

\description{Estimate the prior weight, prior.n, using an approximate empirical Bayes rule given the estimate of the common dispersion. The prior weight determines how much smoothing takes place to squeeze tag/genewise estimates of the dispersion closer to the estimate of the common dispersion.}


\usage{
estimateSmoothing(object,verbose=TRUE) 
}
\arguments{ 
\item{object}{\code{DGEList} object, output of \code{estimateCommonDisp}}

\item{verbose}{logical, whether to write comments, default \code{true}}
}

\value{ \code{estimateSmoothing} produces an object of class \code{DGEList} with the following components.
	\item{prior.n}{scalar; estimate of the prior weight, i.e. the smoothing parameter that indicates the weight to put on the common likelihood compared to the individual tag's likelihood; prior.n of 10 means that the common likelihood is given 10 times the weight of the individual tag/gene's likelihood in the estimation of the tag/genewise dispersion}
}

\details{We are not recommending this function for routine use at the moment, as it has given unexpected results on some deep-sequenced data sets. It should be considered experimental. We are instead recommending that \code{prior.n} be chosen by the user.
Values in the range 10-50 give good results in practice.}

\author{Mark Robinson, Davis McCarthy}

\examples{
y<-matrix(rnbinom(20,size=1,mu=10),nrow=5)
d<-DGEList(counts=y,group=rep(1:2,each=2),lib.size=rep(c(1000:1001),2))
d<-estimateCommonDisp(d)
prior.n<-estimateSmoothing(d)
}
\keyword{file}

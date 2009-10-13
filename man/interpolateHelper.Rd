\name{interpolateHelper}
\alias{interpolateHelper}

\title{Quantile Adjustment interpolator}

\description{Helper function to interpolate the quantile function. This is the function that actually generates the pseudocounts required by quantileAdjust to adjust (normalise) the library sizes and estimate the dispersion parameter. Given fixed values of the estimated mean (mu) and proportion in the library (p) for each tag, as well as a fixed (tagwise or common) value for the dispersion parameter (r), the function interpolates linearly the quantiles used as pseudocounts. If any value of r is infinite (corresponding to phi=0, the dispersion parameter for the negative binomial model), then a Poisson model is used, as setting phi=0 in the negative binomial model is equivalent to using a Poisson model. Otherwise, quantiles are calculated from the negative binomial distribution. }

\usage{ 
interpolateHelper(mu, p, r, count.max, verbose=TRUE) 
}

\arguments{ 

\item{mu}{matrix of means}

\item{p}{matrix of percentiles}

\item{r}{scalar, vector or matrix of \code{size} parameters}

\item{count.max}{vector of maximum counts for all tags}

\item{verbose}{whether to write comments, default \code{true}}
}

\value{matrix with quantile-adjusted pseudocounts }

\author{Mark Robinson, Davis McCarthy}

\examples{
y<-matrix(rnbinom(10000,size=2,mu=10),ncol=4)
d<-DGEList(counts=y,group=rep(1:2,each=2),lib.size=rep(c(1000,1010),2))
conc<-estimatePs(d,r=2)
N<-prod(d$lib.size)^(1/ncol(d$counts))
perc<-pnbinom(d$counts-1,size=2,mu=outer(conc$conc.common,d$samples$lib.size))+dnbinom(d$counts,size=2,mu=outer(conc$conc.common,d$samples$lib.size))/2
maxcounts<-apply(d$counts,1,max)
pseudo<-interpolateHelper(outer(conc$conc.common,rep(N,4)),perc,r=2,maxcounts)
}
\keyword{algebra}

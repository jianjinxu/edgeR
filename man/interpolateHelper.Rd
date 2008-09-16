\name{interpolateHelper}

\alias{interpolateHelper}

\title{Quantile Adjustment interpolator}

\description{Helper function to interpolate the quantile function }

\usage{ 
interpolateHelper(mu, p, r, d,verbose=TRUE) 
}

\arguments{ 

\item{mu}{matrix of means}

\item{p}{matrix of percentiles}

\item{r}{scalar, vector or matrix of \code{size} parameters}

\item{d}{original data matrix}

\item{verbose}{whether to write comments, default \code{true}}
}

\value{matrix with quantile-adjusted pseudo data }

\author{Mark Robinson}

\examples{
y<-matrix(rnbinom(10000,size=2,mu=10),ncol=4)
d<-list(data=y,group=rep(1:2,each=2),lib.size=rep(c(1000,1010),2))
ps<-estimatePs(d$data[,1:2],d$data[,3:4],d$lib.size[1:2],d$lib.size[3:4],r=2)
N<-prod(d$lib.size)^(1/ncol(d$data))
perc<-pnbinom(d$data-1,size=2,mu=outer(ps$p,d$lib.size))+dnbinom(d$data,size=2,mu=outer(ps$p,d$lib.size))/2
pseudo<-interpolateHelper(outer(ps$p,rep(N,4)),perc,r=2,d$data)
}
\keyword{algebra}

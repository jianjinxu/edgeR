\name{exactTestNB}

\alias{exactTestNB}

\title{An exact test for differences between two negative binomial groups}

\description{An exact test for differences between two negative binomial groups}

\usage{ 
exactTestNB(pseudo, group, pair=c(1,2), mus, r, verbose=TRUE)
}

\arguments{ 

\item{pseudo}{data (e.g. quantile adjusted pseudodata) on which to compute Fisher exact statistics}

\item{group}{group indicator, must be same length as ncol(pseudo)}

\item{pair}{pair of groups to be compared; default compares groups 1 and 2}

\item{mus}{vector of means under the null hypothesis (of no difference between groups)}

\item{r}{preset or estimated negative binomial \code{size} parameter.  If you want to run a Poisson test, set r very large (e.g. 1000)}

\item{verbose}{whether to write comments, default \code{TRUE}}
}

\value{ vector of p-values}

\author{Mark Robinson, Davis McCarthy}
\examples{
y<-matrix(rnbinom(20,mu=10,size=1.5),nrow=5)
group<-factor(c(1,1,2,2))
mus<-rep(10,5)
f<-exactTestNB(y,group,pair=c(1,2),mus,r=1.5)
}
\keyword{algebra}

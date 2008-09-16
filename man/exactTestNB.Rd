\name{exactTestNB}

\alias{exactTestNB}

\title{An exact test for differences between two negative binomial groups}

\description{An exact test for differences between two negative binomial groups}

\usage{ 
exactTestNB(y, g, mus, r, verbose=TRUE)
}

\arguments{ 

\item{y}{data (e.g. quantile adjusted pseudodata) to compute Fisher exact statistics on}

\item{g}{group indicator, must be same length as nrow(y)}

\item{mus}{vector of means under the null hypothesis (of no difference between groups)}

\item{r}{preset or estimated negative binomial \code{size} parameter.  If you want to run a Poisson test, set r very large (e.g. 1000)}

\item{verbose}{whether to write comments, default \code{true}}
}

\value{ list with elements \code{lr} (likelihood ratio test), \code{r} (estimates of 1/overdispersion), \code{ps} (list containing proportion estimates)}
\author{Mark Robinson}
\examples{
y<-matrix(rnbinom(20,mu=10,size=1.5),nrow=5)
group<-c(1,1,2,2)
mus<-rep(10,5)
f<-exactTestNB(y,group,mus,r=1.5)
}
\keyword{algebra}

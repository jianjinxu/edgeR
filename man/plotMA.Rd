\name{plotMA}

\alias{plotMA}
\alias{plotMA,deDGEList-method}

\title{MA-like plot for deDGEList objects}

\description{Plots }

\usage{ 
plotMA(object,xlab="A",ylab="M",pch=19,...)
}
\arguments{ 

\item{object}{\code{deDGEList} object, as output from \code{deDGE}}

\item{xlab}{x-axis label}

\item{ylab}{y-axis label}

\item{pch}{plot character}

\item{...}{further arguments to the \code{plot} command}

}

\value{A plot to the current device}

\seealso{\code{deDGE}}

\author{Mark Robinson}

\examples{
# generate raw data from NB, create list object
y<-matrix(rnbinom(20,size=1,mu=10),nrow=5)
d<-DGEList(data=y,group=rep(1:2,each=2),lib.size=rep(c(1000:1001),2))

# find alpha and call main procedure to find differences
alpha<-alpha.approxeb(d)
ms<-deDGE(d,alpha=alpha$alpha)

# plot it
plotMA(ms)
}

\keyword{file}

\name{plotMA}

\alias{plotMA}

\title{MA-like plot for deDGEList objects}

\description{Plots M (log-abundance ratio) against A (log-average abundance) for two groups}

\usage{ 
plotMA(object,pair=c(1,2),xlab="A",ylab="M",ylim=NULL,pch=19,...)
}
\arguments{ 

\item{object}{\code{deDGEList} object, as output from \code{deDGE}}

\item{pair}{pair of groups to be plotted; default plots groups 1 and 2}

\item{xlab}{x-axis label}

\item{ylab}{y-axis label}

\item{ylim}{limits on y-axis, if left at \code{NULL}, scaled to be symmetric about 0}

\item{pch}{plot character}

\item{...}{further arguments to the \code{plot} command}

}

\value{A plot to the current device}

\seealso{\code{deDGE}}

\author{Mark Robinson, Davis McCarthy}

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

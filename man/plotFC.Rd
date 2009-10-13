\name{plotFC}
\alias{plotFC}

\title{Plot of fold change against abundance for DGE data}

\description{Plots logFC (log-abundance ratio, i.e. fold change) against logConc (log-average abundance, i.e. average concentration across groups) for two groups}

\usage{ 
plotFC(object,xlab="log2Conc",ylab="log2FC",ylim=NULL,pch=19,cex=0.2,...)
}
\arguments{ 

\item{object}{a \code{data.frame} containing (at least) the elements \code{logConc} (log-average concentration/abundance) and \code{logFC} (log-abundance ratio, i.e. fold change) for the two groups being compared, as output from \code{exactTest}}

\item{xlab}{x-axis label}

\item{ylab}{y-axis label}

\item{ylim}{limits on y-axis, if left at \code{NULL}, scaled to be symmetric about 0}

\item{pch}{plot character}

\item{cex}{numerical value giving the amount by which plotting text and symbols should be magnified relative to the default; default \code{cex=0.2} to make the plotted points smaller}

\item{...}{further arguments to the \code{plot} command}

}

\value{A plot to the current device}

\seealso{\code{deDGE}}

\author{Mark Robinson, Davis McCarthy}

\examples{
# generate raw counts from NB, create list object
y<-matrix(rnbinom(2000,size=1,mu=10),nrow=500)
d<-DGEList(counts=y,group=rep(1:2,each=2),lib.size=rep(c(1000:1001),2))

# estimate the common dispersion and find differences
d<-estimateCommonDisp(d)
de<-exactTest(d)

# plot it
plotFC(de)
}

\keyword{file}

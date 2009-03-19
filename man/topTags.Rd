\name{topTags}

\alias{topTags}

\title{Displays the top differentially expressed tags in a table}

\description{Displays/Returns the top DE tags in a data frame}

\usage{ 
topTags(object,pair,n=10,adj.method= "BH", verbose=TRUE)
}

\arguments{ 
\item{object}{\code{deDGEList}, output from \code{deDGE}}

\item{pair}{pair of groups to be compared}

\item{n}{number of tags to display/return}

\item{adj.method}{method used to adjust P-values, using \code{p.adjust}}

\item{verbose}{whether to write comments, default \code{TRUE}}
}

\value{ Data frame containing the relative level of expression, log fold changes, unadjusted and adjusted P-values}

\author{Mark Robinson, Davis McCarthy}
\examples{
# generate raw data from NB, create list object
y<-matrix(rnbinom(80,size=1,mu=10),nrow=20)
d<-DGEList(data=y,group=rep(1:2,each=2),lib.size=rep(c(1000:1001),2))
rownames(d$data)<-paste("tagno",1:nrow(d$data),sep=".")

# find alpha and call main procedure to find differences
alpha<-alpha.approxeb(d)
ms<-deDGE(d,alpha=alpha$alpha)

# look at top 10
topTags(ms)
}

\references{
Robinson MD, Smyth GK. 'Small-sample estimation of negative binomial dispersion, with applications to SAGE data.' Biostatistics. 2008 Apr;9(2):321-32.

Robinson MD, Smyth GK. 'Moderated statistical tests for assessing differences in tag abundance.' Bioinformatics. 2007 Nov 1;23(21):2881-7.
}

\keyword{algebra}

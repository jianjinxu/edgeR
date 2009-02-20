\name{splitIntoGroupsPseudo}

\alias{splitIntoGroupsPseudo}

\title{Split pseudodata according to group}

\description{Given a pair of groups, split pseudodata for these groups}

\usage{ 
splitIntoGroupsPseudo(pseudo, group, pair)
}

\arguments{ 
\item{pseudo}{data matrix to be split (e.g. quantile adjusted pseudodata)}

\item{group}{group indicator, must be same length as ncol(pseudo)}

\item{pair}{pair of groups to be split from the data}
}

\value{ list in which each element is a matrix of count data for an individual group}

\author{Davis McCarthy}
\examples{
# generate raw data from NB, create list object
y<-matrix(rnbinom(80,size=1,mu=10),nrow=20)
d<-DGEList(data=y,group=rep(1:2,each=2),lib.size=rep(c(1000:1001),2))
rownames(d$data)<-paste("tagno",1:nrow(d$data),sep=".")
z<-splitIntoGroupsPseudo(d$data,d$group,pair=c(1,2))
}

\keyword{algebra}
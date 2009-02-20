\name{splitIntoGroups}

\alias{splitIntoGroups}

\title{Split the data from a DGEList object according to group}

\description{Split the data from a DGEList object according to group}

\usage{ 
splitIntoGroups(object)
}

\arguments{ 
\item{object}{\code{DGEList}, list containing the raw data with elements \code{data} (table of counts), \code{group} (vector indicating group) and \code{lib.size} (vector of library sizes)}
}

\value{ list in which each element is a matrix of count data for an individual group}

\author{Davis McCarthy}
\examples{
# generate raw data from NB, create list object
y<-matrix(rnbinom(80,size=1,mu=10),nrow=20)
d<-DGEList(data=y,group=rep(1:2,each=2),lib.size=rep(c(1000:1001),2))
rownames(d$data)<-paste("tagno",1:nrow(d$data),sep=".")
z<-splitIntoGroups(d)
}

\keyword{algebra}

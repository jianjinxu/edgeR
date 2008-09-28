\name{getData}

\alias{getData}

\title{Extract data table from DGEList object}

\description{Returns the \code{data} slot of a DGEList object}

\usage{ 
getData(object) 
}
\arguments{ 

\item{object}{matrix of counts for first group}

}

\value{ matrix of data (presumably integers) }

\author{Mark Robinson}

\examples{
# generate raw data from NB, create list object
y<-matrix(rnbinom(20,size=1,mu=10),nrow=5)
d<-DGEList(data=y,group=rep(1:2,each=2),lib.size=rep(c(1000:1001),2))
# should be 5x4
print(dim(getData(d)))
}
\keyword{file}

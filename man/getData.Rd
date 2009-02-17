\name{getData}

\alias{getData}

\title{Extract data table from DGEList object}

\description{Returns the \code{data} slot of a DGEList object}

\usage{ 
getData(object) 
}
\arguments{ 

\item{object}{list containing the raw data with elements \code{data} (table of counts), \code{group} (vector indicating group) and \code{lib.size} (vector of library sizes)}
}

\value{ matrix of data (presumably integers) }

\author{Mark Robinson, Davis McCarthy}

\examples{
# generate raw data from NB, create list object
y<-matrix(rnbinom(20,size=1,mu=10),nrow=5)
d<-DGEList(data=y,group=rep(1:2,each=2),lib.size=rep(c(1000:1001),2))
# should be 5x4
print(dim(getData(d)))
}
\keyword{file}

\name{getCounts}
\alias{getCounts}

\title{Extract table of counts from DGEList object}

\description{Returns the \code{counts} slot of a DGEList object}

\usage{ 
getCounts(object) 
}
\arguments{ 

\item{object}{list containing the raw counts with elements \code{counts} (table of counts), \code{group} (vector indicating group) and \code{lib.size} (vector of library sizes)}
}

\value{ matrix of counts (presumably integers) }

\author{Mark Robinson, Davis McCarthy}

\examples{
# generate raw counts from NB, create list object
y<-matrix(rnbinom(20,size=1,mu=10),nrow=5)
d<-DGEList(counts=y,group=rep(1:2,each=2),lib.size=rep(c(1000:1001),2))
# should be 5x4
print(dim(getCounts(d)))
}
\keyword{file}

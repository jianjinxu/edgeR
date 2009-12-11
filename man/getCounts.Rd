\name{getCounts}
\alias{getCounts}

\title{Extract Table of Counts from DGEList Object}

\description{Returns the \code{counts} slot of a DGEList object}

\usage{ 
getCounts(object) 
}
\arguments{ 

\item{object}{\code{DGEList} object containing (at least) the elements \code{counts} (table of raw counts), \code{group} (factor indicating group) and \code{lib.size} (numeric vector of library sizes)}
}

\value{\code{getCounts} returns a matrix of counts (presumably integers) }

\author{Mark Robinson, Davis McCarthy}

\examples{
# generate raw counts from NB, create list object
y<-matrix(rnbinom(20,size=1,mu=10),nrow=5)
d<-DGEList(counts=y,group=rep(1:2,each=2),lib.size=rep(c(1000:1001),2))
# should be 5x4
print(dim(getCounts(d)))
}

\seealso{
\code{\link{DGEList}} for more information about the \code{DGEList} class.
}

\keyword{file}

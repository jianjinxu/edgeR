\name{getOffset}
\alias{getOffset}

\title{Extract Vector of Offsets from DGEList Object}

\description{Returns the \code{lib.size} component of the \code{samples} component of  DGEList object multiplied by the \code{norm.factors} component}

\usage{
getOffset(object) 
}
\arguments{ 

\item{object}{\code{DGEList} object containing (at least) the elements \code{counts} (table of raw counts), \code{group} (factor indicating group) and \code{samples}, which contains \code{lib.size} (numeric vector of library sizes) and \code{norm.factors} (numeric vector of normalization factors).}
}

\value{\code{getOffset} returns a numeric vector }

\author{Gordon Smyth, Davis McCarthy}

\examples{
# generate raw counts from NB, create list object
y<-matrix(rnbinom(20,size=1,mu=10),nrow=5)
d<-DGEList(counts=y,group=rep(1:2,each=2),lib.size=rep(c(1000:1001),2))
getOffset(d)
}

\seealso{
\code{\link{DGEList}} for more information about the \code{DGEList} class.
\code{\link{as.matrix.DGEList}}.
}

\keyword{file}

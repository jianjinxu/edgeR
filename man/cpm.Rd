\name{cpm}
\alias{cpm}

\title{Calculate Counts per Millioni from DGEList or Matrix Object}

\description{Returns counts per million from a  DGEList or matrix object by dividing raw counts by library size (which can be normalized) and multiplying by one million.}

\usage{
cpm(x, normalized.lib.sizes=FALSE)
}
\arguments{ 

\item{x}{either a matrix of counts or a \code{DGEList} object with (at least) elements \code{counts} (table of unadjusted counts) and \code{samples} (data frame containing information about experimental group, library size and normalization factor for the library size)}

\item{normalized.lib.sizes}{logical, should the library sizes (total sum of counts for each library) be normalized using the \code{norm.factors} component of the \code{DGEList} object? Ignored (with a warning) if \code{x} is a count matrix.}

}

\value{\code{getPriorN} returns a numeric scalar }

\details{
A convenience function to compute the counts per million for plotting and comparing libraries on a convenient scale. Essentially just does the calculation \code{1e06*t(t(x)/lib.size)} to produce counts per million, where \code{x} is a matrix of counts and the \code{lib.size} can be the total sum of counts in each library or a normalized version of this using TMM normalization or equivalent method.
}

\author{Davis McCarthy, Gordon Smyth}

\examples{
# generate raw counts from NB, create list object
y<-matrix(rnbinom(20,size=1,mu=10),nrow=5)
cpm(y)
d<-DGEList(counts=y,group=rep(1:2,each=2),lib.size=rep(c(1000:1001),2))
# When applied to a DGEList object, x$samples$lib.size is used
cpm(d)
# As x$samples$lib.size here is very different from colSums(y), cpm(y) and cpm(d) give very different results for the counts per million
}

\seealso{
\code{\link{DGEList}} for more information about the \code{DGEList} class.
}

\keyword{file}

\name{splitIntoGroups}
\alias{splitIntoGroups}
\alias{splitIntoGroupsPseudo}


\title{Split the Counts or Pseudocounts from a DGEList Object According To Group}

\description{Split the counts from a DGEList object according to group, creating a list where each element consists of a numeric matrix of counts for a particular experimental group. Given a pair of groups, split pseudocounts for these groups, creating a list where each element is a matrix of pseudocounts for a particular gourp.}

\usage{
splitIntoGroups(object)
splitIntoGroupsPseudo(pseudo, group, pair)
}

\arguments{ 
\item{object}{\code{DGEList}, object containing (at least) the elements \code{counts} (table of raw counts), \code{group} (factor indicating group) and \code{lib.size} (numeric vector of library sizes)}

\item{pseudo}{numeric matrix of quantile-adjusted pseudocounts to be split}

\item{group}{factor indicating group to which libraries/samples (i.e. columns of \code{pseudo} belong; must be same length as ncol(pseudo)}

\item{pair}{vector of length two stating pair of groups to be split for the pseudocounts}
}

\value{\code{splitIntoGroups} outputs a list in which each element is a matrix of count counts for an individual group. \code{splitIntoGroupsPseudo} outputs a list with two elements, in which each element is a numeric matrix of (pseudo-)count data for one of the groups specified.}

\author{Davis McCarthy}
\examples{
# generate raw counts from NB, create list object
y<-matrix(rnbinom(80,size=1,mu=10),nrow=20)
d<-DGEList(counts=y,group=rep(1:2,each=2),lib.size=rep(c(1000:1001),2))
rownames(d$counts)<-paste("tagno",1:nrow(d$counts),sep=".")
z1<-splitIntoGroups(d)

z2<-splitIntoGroupsPseudo(d$counts,d$group,pair=c(1,2))
}

\keyword{algebra}

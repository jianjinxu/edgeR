\name{estimatePs}
\alias{estimatePs}

\title{Estimate expression proportions}

\description{Estimate expression proportions (maximum likelihood with size fixed) based on negative binomial for each tag/gene and sample group. Expression proportions are used to determine overall abundance of each tag/gene and differential expression of tags/genes between groups.}

\usage{ 
estimatePs(object, r, tol = 1e-10, maxit = 30) 
}
\arguments{ 

\item{object}{list containing the raw counts with elements \code{counts} (table of counts), \code{group} (vector indicating group) and \code{lib.size} (vector of library sizes)}

\item{r}{size parameter of negative binomial model (\code{size} = \code{1/phi} where \code{phi} is the dispersion parameter in the NB model)}

\item{tol}{tolerance between iterations}

\item{maxit}{maximum number of iterations}
}

\value{ list with elements \code{conc.common} (vector giving overall proportion/concentration for each tag), \code{conc.group} (matrix with columns giving estimates of tag/gene concentrations (proportion of total RNA for that group that that particular tag/gene contributes) for different groups)}

\author{Mark Robinson, Davis McCarthy}

\examples{
set.seed(0)
y<-matrix(rnbinom(40,size=1,mu=10),ncol=4)
d<-DGEList(counts=y,group=rep(1:2,each=2),lib.size=rep(c(1000:1001),2))
conc<-estimatePs(d,r=1)
}
\keyword{file}

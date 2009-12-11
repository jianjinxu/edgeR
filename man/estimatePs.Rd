\name{estimatePs}
\alias{estimatePs}

\title{Estimate Expression Levels}

\description{Estimate expression levels (i.e. proportion of all sample mRNA corresponding to each tag; or, concentration of mRNA for each tag in sample mRNA) using maximum likelihood with dispersion parameter fixed based on the negative binomial model for each tag/gene and sample group. Expression proportions are used to determine overall abundance of each tag/gene and differential expression of tags/genes between groups.}

\usage{ 
estimatePs(object, r, tol = 1e-10, maxit = 30) 
}

\arguments{ 

\item{object}{list containing (at least) the elements \code{counts} (table of counts), \code{group} (vector or factor indicating group) and \code{lib.size} (numeric vector of library sizes)}

\item{r}{numeric vector providing the size parameter of negative binomial model (\code{size} = \code{1/phi} where \code{phi} is the dispersion parameter in the NB model)}

\item{tol}{numeric scalar, tolerance between iterations}

\item{maxit}{positive integer scalar, maximum number of iterations}
}

\value{A list with elements:
	\item{conc.common}{numeric vector giving overall proportion/concentration for each tag}
	\item{conc.group}{numeric matrix with columns giving estimates of tag/gene concentrations (proportion of total RNA for that group that that particular tag/gene contributes) for different groups)}
}

\details{
The Newton-Raphson method is used to calculate iteratively the maximum likelihood estimate of the expression level (i.e. concentration of mRNA for a particular tag in the sample mRNA) for each tag/gene.
}


\author{Mark Robinson, Davis McCarthy}

\examples{
set.seed(0)
y<-matrix(rnbinom(40,size=1,mu=10),ncol=4)
d<-DGEList(counts=y,group=rep(1:2,each=2),lib.size=rep(c(1000:1001),2))
conc<-estimatePs(d,r=1)
}
\keyword{file}

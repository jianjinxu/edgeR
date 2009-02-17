\name{findMaxD2}

\alias{findMaxD2}

\title{Maximizes the negative binomial likelihood}

\description{Maximizes the negative binomial likelihood (a weighted version using the common likelihood given weight alpha) for each tag}


\usage{ 
findMaxD2(object, alpha = 0.5, grid = TRUE, tol = 1e-05, n.iter = 10, grid.length = 200)
}

\arguments{ 

\item{object}{list containing the raw data with elements \code{data} (table of counts), \code{group} (vector indicating group) and \code{lib.size} (vector of library sizes)}

\item{alpha}{weight given to common likelihood, set to 0 for individual estimates or large (e.g. 100) for common likelihood}

\item{grid}{logical, whether to use a grid search (default = \code{TRUE}); if \code{FALSE} use Newton-Rhapson steps}

\item{tol}{if \code{grid=FALSE}, tolerance for Newton-Rhapson iterations}

\item{n.iter}{if \code{grid=FALSE}, number of Newton-Rhapson iterations}

\item{grid.length}{length of the grid over which to maximize; default \code{200}}

}

\value{ vector of the values of delta that maximize the negative binomial likelihood for each tag  (where delta = \code{phi / (phi+1)} and \code{phi} is the overdispersion parameter) }

\author{Mark Robinson, Davis McCarthy}
\examples{
y<-matrix(rnbinom(1000,mu=10,size=2),ncol=4)
d<-list(data=y,group=c(1,1,2,2),lib.size=c(1000:1003))
cml1<-findMaxD2(d,alpha=10)
cml2<-findMaxD2(d,alpha=0)
}
\keyword{algebra}

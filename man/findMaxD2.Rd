\name{findMaxD2}

\alias{findMaxD2}

\title{Maximizes the negative binomial likelihood}

\description{Maximizes the negative binomial likelihood (a weighted version using the common likelihood given weight alpha) for each tag}


\usage{ 
findMaxD2(x, alpha = 0.5, grid = TRUE, tol = 1e-05, n.iter = 5, grid.length = 200)
}

\arguments{ 

\item{x}{list with elements \code{data}, \code{lib.size} and \code{group}}

\item{alpha}{weight given to common likelihood, set to 0 for individual estimates or large (e.g. 100) for common likelihood}

\item{grid}{logical, whether to use a grid search (default = \code{TRUE}); if \code{FALSE} use Newton-Rhapson steps}

\item{tol}{if \code{grid=FALSE}, tolerance for Newton-Rhapson iterations}

\item{n.iter}{if \code{grid=FALSE}, number of Newton-Rhapson iterations}

\item{grid.length}{length of the grid to maximize over; default \code{200}}

}

\value{ list with elements \code{lr} (likelihood ratio test), \code{r} (estimates of 1/overdispersion), \code{ps} (list containing proportion estimates)}
\author{Mark Robinson}
\examples{
y<-matrix(rnbinom(1000,mu=10,size=2),ncol=4)
d<-list(data=y,group=c(1,1,2,2),lib.size=c(1000:1003))
cml1<-findMaxD2(d,alpha=10)
cml2<-findMaxD2(d,alpha=0)
}
\keyword{algebra}

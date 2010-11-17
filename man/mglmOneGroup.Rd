\name{mglmOneGroup}
\alias{mglmOneGroup}
\alias{mglmLS}

\title{Fit a one-group negative binomial generalized linear model for each transcript}

\description{Fits a negative binomial generalized linear model for each transcript (tag) with the unadjusted counts provided, a value for the dispersion parameter and, optionally, offsets and weights for different libraries or transcripts, under the assumption that there is just one experimental group. }


\usage{ 
mglmOneGroup(y, offset=0, dispersion=0, maxit=50, verbose=FALSE)
mglmLS(y,design,dispersion=0,start=NULL,offset=0,tol=1e-5,maxit=50,trace=FALSE)
}

\arguments{ 

\item{y}{matrix containing the raw counts for each library (the measure of expression level) and gene}

\item{offset}{matrix giving the offset that is to be included in the NB GLM for the transcripts. If a scalar, then this value will be used as an offset for all transcripts and libraries.}

\item{design}{matrix giving the design matrix of the GLM.}

\item{dispersion}{numeric scalar or vector providing the value for the dispersion parameter that is used in fitting the GLM for each transcript. Can be a common value for all tags, or a vector of values can provide a unique dispersion value for each tag.}

\item{maxit}{scalar giving the maximum number of iterations to do for the Fisher scoring algorithm.}

\item{verbose}{logical, whether or not to print information about the convergence of the algorithm; default is \code{FALSE}.}
}

\value{ \code{glmFit} produces a vector of length equal to the number of tags/genes (number of rows of \code{y}) providing the single coefficent from the GLM fit for each tag/gene. This can be interpreted as a measure of the 'average expression' level of the tag/gene.
	
}

\details{
Given a fixed value for the dispersion parameter, a negative binomial model can be fitted to the counts for each tag/transcript in a dataset. The function fits the simplest possible model--just an intercept term--to each tag/gene in order to obtain an estimate of the 'average expression' level for each tag/gene on the natural log scale. A vectorized Fisher scoring algorithm is used to fit the models to the tags/genes.
}

\author{Davis McCarthy, Gordon Smyth}
\examples{
y<-matrix(rnbinom(1000,mu=10,size=2),ncol=4)
dispersion <- 0.1
## Fit the NB GLM to the counts
ave.expression <- glmNBOneGroup(y, dispersion=dispersion)
head(ave.expression)
}

\seealso{
\code{\link{glmFit}}, for more complicated GLM modelling for DGE data.
}

\keyword{algebra}



\name{estimateTagwiseDisp}
\alias{estimateTagwiseDisp}

\title{Estimate Empirical Bayes Tagwise Dispersion Values}

\description{
Estimates tagwise dispersion values by an empirical Bayes method based on weighted conditional maximum likelihood.
}

\usage{
estimateTagwiseDisp(object, prior.n=getPriorN(object), trend="movingave", prop.used=0.3, method="grid", grid.length=11, grid.range=c(-6,6), tol=1e-06, verbose=FALSE)
}

\arguments{ 

\item{object}{object of class \code{DGEList} containing (at least) the elements \code{counts} (table of raw counts), \code{group} (factor indicating group), \code{lib.size} (numeric vector of library sizes) and \code{pseudo.alt} (numeric matrix of quantile-adjusted pseudocounts calculated under the alternative hypothesis of a true difference between groups; recommended to use the \code{DGEList} object provided as the output of \code{estimateCommonDisp}}

\item{prior.n}{numeric scalar, smoothing parameter that indicates the weight to give to the common likelihood compared to the individual tag's likelihood; default \code{getPriorN(object)} gives a value for \code{prior.n} that is equivalent to giving the common likelihood 20 prior degrees of freedom in the estimation of the tag/genewise dispersion.}

\item{trend}{method for allowing the prior distribution for the dispersion to be abundance-dependent.
Possible values are \code{"none"}, \code{"movingave"} and \code{"tricube"}.
\code{"none"} means no trend.
\code{"movingave"} applies a moving average smoother to the local likelihood values.
\code{"tricube"} applies tricube weighting to locally smooth the common likelihood.}

\item{prop.used}{optional scalar giving the proportion of all tags/genes to be used for the locally weighted estimation of the tagwise dispersion, allowing the dispersion estimates to vary with abundance (expression level). For each tag/gene the estimate of its dispersion is based on the closest \code{prop.used} of all of the genes to that gene, where 'closeness' is based on similarity in expression level.}

\item{method}{method for maximizing the posterior likelihood.
Possible values are \code{"grid"} for interpolation on grid points or \code{"optimize"} to call the function of the same name.
}

\item{grid.length}{for \code{method="grid"}, the number of points on which the interpolation is applied for each tag.}
\item{grid.range}{for \code{method="grid"}, the range of the grid points around the trend on a log2 scale.}

\item{tol}{for \code{method="optimize"}, the tolerance for Newton-Rhapson iterations.}

\item{verbose}{logical, if \code{TRUE} then diagnostic ouput is produced during the estimation process.}

}

\details{
Maximizes the negative binomial weighted likelihood (a weighted version using the common likelihood given weight according the the smoothing parameter prior.n and the individual tag/gene likelihood) for each tag from the pseudocounts provided (i.e. assuming library sizes are equal), to give an estimate of the dispersion parameter for each tag (i.e. tagwise dispersion estimation).

\code{"tricube"} local weighting is similar to that used by \code{lowess}.
\code{"movingave"} is much faster than \code{"tricube"} and gives similar results.

\code{"optimize"} is very slow if there is a large number of tags/genes to be analysed (i.e., more than 5000).
}

\value{
An object of class \code{DGEList} with the same components as for \code{\link{estimateCommonDisp}} plus the following:
	\item{prior.n}{estimate of the prior weight, i.e. the smoothing parameter that indicates the weight to put on the common likelihood compared to the individual tag's likelihood; prior.n of 10 means that the common likelihood is given 10 times the weight of the individual tag/gene's likelihood in the estimation of the tag/genewise dispersion}
	\item{tagwise.dispersion}{tag- or gene-wise estimates of the dispersion parameter}
}

\references{
Robinson MD and Smyth GK (2007). Moderated statistical tests for
assessing differences in tag abundance. \emph{Bioinformatics} 23, 2881-2887
}

\author{Mark Robinson, Davis McCarthy and Gordon Smyth}
\examples{
y<-matrix(rnbinom(1000,mu=10,size=2),ncol=4)
d<-DGEList(counts=y,group=c(1,1,2,2),lib.size=c(1000:1003))
d<-estimateCommonDisp(d)
tgwdisp<-estimateTagwiseDisp(d, prior.n=10)
}

\seealso{
\code{\link{estimateCommonDisp}} estimates a common value for the dispersion parameter for all tags/genes - should generally be run before \code{estimateTagwiseDisp}.
}

\keyword{algebra}

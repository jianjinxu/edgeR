\name{estimateTagwiseDisp}
\alias{estimateTagwiseDisp}

\title{Maximizes the Negative Binomial Weighted Conditional Likelihood}

\description{Maximizes the negative binomial weighted likelihood (a weighted version using the common likelihood given weight according the the smoothing parameter prior.n and the individual tag/gene likelihood) for each tag from the pseudocounts provided (i.e. assuming library sizes are equal), to give an estimate of the dispersion parameter for each tag (i.e. tagwise dispersion estimation).}


\usage{
estimateTagwiseDisp(object, prior.n=10, trend=FALSE, prop.used=NULL, tol=1e-06, grid=TRUE, grid.length=200, method="movingave", verbose=TRUE)
}

\arguments{ 

\item{object}{a \code{DGEList} object containing (at least) the elements \code{counts} (table of raw counts), \code{group} (factor indicating group), \code{lib.size} (numeric vector of library sizes) and \code{pseudo,alt} (numeric matrix of quantile-adjusted pseudocounts calculated under the alternative hypothesis of a true difference between groups; recommended to use the \code{DGEList} object provided as the output of \code{estimateCommonDisp}}

\item{prior.n}{numeric scalar, smoothing parameter that indicates the weight to give to the common likelihood compared to the individual tag's likelihood; default \code{10} means that the common likelihood is given 10 times the weight of the individual tag/gene's likelihood in the estimation of the tag/genewise dispersion}

\item{trend}{logical, whether or not to let the tagwise dispersion estimates vary with tag/gene abundance (expression level), that is, whether or not to allow a trend with tag abundance in the tagwise dispersion estimates}

\item{prop.used}{optional scalar giving the proportion of all tags/genes to be used for the locally weighted estimation of the tagwise dispersion, allowing the dispersion estimates to vary with abundance (expression level). If \code{NULL}, then a default value of 0.4 (i.e. 40 per cent of tags) are used. That means that for each tag/gene the estimate of its dispersion is based on the closest 40 per cent of all of the genes to that gene, where 'closeness' is based on similarity in expression level.}

\item{tol}{numeric scalar, if \code{grid=FALSE}, tolerance for Newton-Rhapson iterations}

\item{grid}{logical, whether to use a grid search (default = \code{TRUE}); if \code{FALSE}, uses \code{optimize}, but this is very slow if there is a large number of tags/genes to be analysed (i.e. more than 5000)}

\item{grid.length}{if \code{grid=TRUE}, the number of points at which the likelihood is evaluated for each tag, so larger values improve the accuracy of the dispersion estimates; default \code{1000}}

\item{method}{names the method to use to fit the trend (so only relevant if \code{trend=TRUE}). The default is \code{"movingave"}, which uses a moving average to locally smooth the common likelihood to apply a mean-dependent trend to the dispersions. The other option is \code{"tricube"}, which uses a tricube weighting function to locally smooth the common likelihood. The moving average method is much faster and results are similar.}

\item{verbose}{logical, whether to write comments, default \code{TRUE}}

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

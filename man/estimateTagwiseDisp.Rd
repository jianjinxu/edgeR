\name{estimateTagwiseDisp}
\alias{estimateTagwiseDisp}

\title{Estimate Empirical Bayes Tagwise Dispersion Values}

\description{
Estimates tagwise dispersion values by an empirical Bayes method based on weighted conditional maximum likelihood.
}

\usage{
estimateTagwiseDisp(object, prior.df=10, trend="movingave", span=NULL, method="grid",
          grid.length=11, grid.range=c(-6,6), tol=1e-06, verbose=FALSE)
}

\arguments{
\item{object}{object of class \code{DGEList} containing (at least) the elements \code{counts} (table of raw counts), \code{group} (factor indicating group), \code{lib.size} (numeric vector of library sizes) and \code{pseudo.alt} (numeric matrix of quantile-adjusted pseudocounts calculated under the alternative hypothesis of a true difference between groups; recommended to use the \code{DGEList} object provided as the output of \code{estimateCommonDisp}}

\item{prior.df}{prior degrees of freedom.}

\item{trend}{method for estimating dispersion trend. Possible values are \code{"none"}, \code{"movingave"} and \code{"loess"}.}

\item{span}{width of the smoothing window, as a proportion of the data set.}

\item{method}{method for maximizing the posterior likelihood.
Possible values are \code{"grid"} for interpolation on grid points or \code{"optimize"} to call the function of the same name.
}

\item{grid.length}{for \code{method="grid"}, the number of points on which the interpolation is applied for each tag.}
\item{grid.range}{for \code{method="grid"}, the range of the grid points around the trend on a log2 scale.}

\item{tol}{for \code{method="optimize"}, the tolerance for Newton-Rhapson iterations.}

\item{verbose}{logical, if \code{TRUE} then diagnostic ouput is produced during the estimation process.}

}

\details{
This function implements the empirical Bayes strategy proposed by Robinson and Smyth (2007)
for estimating the tagwise negative binomial dispersions.
The experimental design is assumed to be a oneway layout with one or more experimental groups.
The empirical Bayes posterior is implemented as a conditional likelihood with tag-specific weights.

The prior values for the dispersions are determined by a global trend.
The individual tagwise dispersions are then squeezed towards this trend.
The prior degrees of freedom determines the weight given to the prior.
The larger the prior degrees of freedom, the more the tagwise dispersions are squeezed towards the global trend.
If the number of libraries is large, the prior becomes less important and the tagwise dispersion are determined more by the individual tagwise data.

If \code{trend="none"}, then the prior dispersion is just a constant, the common dispersion.
Otherwise, the trend is determined by a moving average (\code{trend="movingave"}) or loess smoother applied to the tagwise conditional log-likelihood.
\code{method="loess"} applies a loess curve of degree 0 as implemented in \code{\link{loessByCol}}.

\code{method="optimize"} is not recommended for routine use as it is very slow.
It is included for testing purposes.
}

\value{
An object of class \code{DGEList} with the same components as for \code{\link{estimateCommonDisp}} plus the following:
	\item{prior.n}{estimate of the prior weight, i.e. the smoothing parameter that indicates the weight to put on the common likelihood compared to the individual tag's likelihood; prior.n of 10 means that the common likelihood is given 10 times the weight of the individual tag/gene's likelihood in the estimation of the tag/genewise dispersion}
	\item{tagwise.dispersion}{tag- or gene-wise estimates of the dispersion parameter}
}

\references{
Robinson, MD, and Smyth, GK (2007). Moderated statistical tests for
assessing differences in tag abundance. \emph{Bioinformatics} 23, 2881-2887.
\url{http://bioinformatics.oxfordjournals.org/content/23/21/2881}
}

\author{Mark Robinson, Davis McCarthy, Yunshun Chen and Gordon Smyth}
\examples{
# See exactTest
}

\seealso{
\code{\link{estimateCommonDisp}} is usually run before \code{estimateTagwiseDisp}.

\code{\link{movingAverageByCol}} and \code{\link{loessByCol}} implement the moving average or loess smoothers.
}

\keyword{algebra}

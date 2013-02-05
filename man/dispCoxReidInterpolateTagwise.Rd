\name{dispCoxReidInterpolateTagwise}
\alias{dispCoxReidInterpolateTagwise}

\title{Estimate Tagwise Dispersion for Negative Binomial GLMs by Cox-Reid Adjusted Profile Likelihood}

\description{
Estimate tagwise dispersion parameters across multiple negative binomial generalized linear models using weighted Cox-Reid Adjusted Profile-likelihood and cubic spline interpolation over a tagwise grid.
}

\usage{
dispCoxReidInterpolateTagwise(y, design, offset=NULL, dispersion, trend=TRUE, AveLogCPM=NULL,
          min.row.sum=5, prior.df=10, span=0.3, grid.npts=11, grid.range=c(-6,6))
}

\arguments{ 
\item{y}{numeric matrix of counts}

\item{design}{numeric matrix giving the design matrix for the GLM that is to be fit.}

\item{offset}{numeric scalar, vector or matrix giving the offset (in addition to the log of the effective library size) that is to be included in the NB GLM for the transcripts. If a scalar, then this value will be used as an offset for all transcripts and libraries. If a vector, it should be have length equal to the number of libraries, and the same vector of offsets will be used for each transcript. If a matrix, then each library for each transcript can have a unique offset, if desired. In \code{adjustedProfileLik} the \code{offset} must be a matrix with the same dimension as the table of counts.}

\item{dispersion}{numeric scalar or vector giving the dispersion(s) towards which the tagwise dispersion parameters are shrunk.}

\item{trend}{logical, whether abundance-dispersion trend is used for smoothing.}

\item{AveLogCPM}{numeric vector giving average log2 counts per million for each tag.}

\item{min.row.sum}{numeric scalar giving a value for the filtering out of low abundance tags. Only tags with total sum of counts above this value are used. Low abundance tags can adversely affect the estimation of the common dispersion, so this argument allows the user to select an appropriate filter threshold for the tag abundance.}

\item{prior.df}{numeric scalar, prior degsmoothing parameter that indicates the weight to give to the common likelihood compared to the individual tag's likelihood; default \code{getPriorN(object)} gives a value for \code{prior.n} that is equivalent to giving the common likelihood 20 prior degrees of freedom in the estimation of the tag/genewise dispersion.}

\item{span}{numeric parameter between 0 and 1 specifying proportion of data to be used in the local regression moving window. Larger numbers give smoother fits.}

\item{grid.npts}{numeric scalar, the number of points at which to place knots for the spline-based estimation of the tagwise dispersion estimates.}

\item{grid.range}{numeric vector of length 2, giving relative range, in terms of \code{log2(dispersion)}, on either side of trendline for each tag for spline grid points.}
}

\value{\code{dispCoxReidInterpolateTagwise} produces a vector of tagwise dispersions having the same length as the number of genes in the count data.
}

\details{
In the \code{edgeR} context, \code{dispCoxReidInterpolateTagwise} is a low-level function called by \code{estimateGLMTagwiseDisp}.

\code{dispCoxReidInterpolateTagwise} calls the function \code{maximizeInterpolant} to fit cubic spline interpolation over a tagwise grid. 
}


\references{
Cox, DR, and Reid, N (1987). Parameter orthogonality and approximate conditional inference. \emph{Journal of the Royal Statistical Society Series B} 49, 1-39.

McCarthy, DJ, Chen, Y, Smyth, GK (2012). Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation.
\emph{Nucleic Acids Research} 40, 4288-4297.
\url{http://nar.oxfordjournals.org/content/40/10/4288}
}

\author{Yunshun Chen, Gordon Smyth}
\examples{
y <- matrix(rnbinom(1000, mu=10, size=2), ncol=4)
design <- matrix(1, 4, 1)
dispersion <- 0.5
d <- dispCoxReidInterpolateTagwise(y, design, dispersion=dispersion)
d
}

\seealso{
\code{\link{estimateGLMTagwiseDisp}}, \code{\link{maximizeInterpolant}}
}

\keyword{algebra}

\name{dispCoxReidSplineTrend}
\alias{dispCoxReidSplineTrend}
\alias{dispCoxReidPowerTrend}

\title{Estimate Dispersion Trend for Negative Binomial GLMs}

\description{
Estimate trended common dispersion parameters across multiple negative binomial generalized linear models using Cox-Reid adjusted profile likelihood.
}

\usage{
dispCoxReidSplineTrend(y, design, offset=NULL, degree = 5, subset=10000, method.optim="Nelder-Mead", trace=0)
dispCoxReidPowerTrend(y, design, offset=NULL, subset=10000, method.optim="Nelder-Mead", trace=0)
}

\arguments{ 

\item{y}{numeric matrix of counts}

\item{design}{numeric matrix giving the design matrix for the GLM that is to be fit.}

\item{offset}{numeric scalar, vector or matrix giving the offset (in addition to the log of the effective library size) that is to be included in the NB GLM for the transcripts. If a scalar, then this value will be used as an offset for all transcripts and libraries. If a vector, it should be have length equal to the number of libraries, and the same vector of offsets will be used for each transcript. If a matrix, then each library for each transcript can have a unique offset, if desired. In \code{adjustedProfileLik} the \code{offset} must be a matrix with the same dimension as the table of counts.}

\item{degree}{integer giving the degree of the spline function.}

\item{subset}{integer, number of rows to use in the calculation.  Rows used are chosen evenly spaced by abundance.}

\item{method.optim}{the method to be used in \code{optim}. See \code{\link{optim}} for more detail.}

\item{trace}{logical, should iteration information be output?}
}

\value{
Numeric vector giving the estimated trended common dispersions. It is of the same length as the number of tags in the count data.
}

\details{
In the \code{edgeR} context, these are low-level functions called by \code{estimateGLMTrendedDisp}.

\code{dispCoxReidSplineTrend} maximizes the Cox-Reid adjusted profile likelihood (Cox and Reid, 1987) by fitting spline interpolation.
\code{dispCoxReidPowerTrend} models the dispersion trend by a power function. The parameters of the power function are estimated by maximizing the Cox-Reid adjusted profile likelihood.
}

\references{
Cox, DR, and Reid, N (1987). Parameter orthogonality and approximate conditional inference. \emph{Journal of the Royal Statistical Society Series B} 49, 1-39.
}

\author{Yunshun Chen, Davis McCarthy, Gordon Smyth}

\examples{
design <- matrix(1,4,1)
y <- matrix((rnbinom(400,mu=100,size=2)),100,4)
dispCoxReidSplineTrend(y, design, degree=3)
dispCoxReidPowerTrend(y, design)
}

\seealso{
\code{\link{estimateGLMTrendedDisp}}, \code{\link{optim}}
}

\keyword{models}

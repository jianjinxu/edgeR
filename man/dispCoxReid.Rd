\name{dispCoxReid}
\alias{dispCoxReid}
\alias{dispDeviance}
\alias{dispPearson}

\title{Estimate Common Dispersion for Negative Binomial GLMs}

\description{
Estimate a common dispersion parameter across multiple negative binomial generalized linear models.
}

\usage{
dispCoxReid(y, design, offset=NULL, interval=c(0,4), tol=1e-5, min.row.sum=5, subset=10000)
dispDeviance(y, design, offset=NULL, interval=c(0,4), tol=1e-5, min.row.sum=5, subset=10000, robust=FALSE, trace=FALSE)
dispPearson(y, design, offset=NULL, interval=c(0,4), tol=1e-5, min.row.sum=5, subset=10000, robust=FALSE, trace=FALSE)
}

\arguments{ 

\item{y}{numeric matrix of counts}

\item{design}{numeric matrix giving the design matrix for the GLM that is to be fit.
Must be of full column rank.
Defaults to a single column of ones, equivalent to treating the columns as replicate libraries.}

\item{offset}{numeric scalar, vector or matrix giving the offset (in addition to the log of the effective library size) that is to be included in the NB GLM for the transcripts. If a scalar, then this value will be used as an offset for all transcripts and libraries. If a vector, it should be have length equal to the number of libraries, and the same vector of offsets will be used for each transcript. If a matrix, then each library for each transcript can have a unique offset, if desired. In \code{adjustedProfileLik} the \code{offset} must be a matrix with the same dimension as the table of counts.}

\item{interval}{numeric vector of length 2 giving allowable values for the dispersion, passed to \code{optimize}.}

\item{tol}{the desired accuracy, see \code{optimize} or \code{uniroot}.}

\item{min.row.sum}{integer. Only rows with at least this number of counts are used.}

\item{subset}{integer, number of rows to use in the calculation.  Rows used are chosen evenly spaced by abundance.}

\item{trace}{logical, should iteration information be output?}

\item{robust}{logical, should a robust estimator be used?}
}

\value{
Numeric vector of length one giving the estimated common dispersion.
}

\details{
These are low-level (non-object-orientated) functions called by \code{estimateGLMCommonDisp}.

\code{dispCoxReid} maximizes the Cox-Reid adjusted profile likelihood (Cox and Reid, 1987).
\code{dispPearson} sets the average Pearson goodness of fit statistics to its (asymptotic) expected value.
This is also known as the \emph{pseudo-likelihood} estimator.
\code{dispDeviance} sets the average residual deviance statistic to its (asymptotic) expected values.
This is also known as the \emph{quasi-likelihood} estimator.

Robinson and Smyth (2008) showed that the Pearson (pseudo-likelihood) estimator typically under-estimates the true dispersion.
It can be seriously biased when the number of libraries (\code{ncol(y)} is small.
On the other hand, the deviance (quasi-likelihood) estimator typically over-estimates the true dispersion when the number of libraries is small.
Robinson and Smyth (2008) showed the Cox-Reid estimator to be the least biased of the three options.

\code{dispCoxReid} uses \code{optimize} to maximize the adjusted profile likelihood, while \code{dispDeviance} and \code{dispPearson} use \code{uniroot} to solve the estimating equation.
The robust options use an order statistic instead the mean statistic, and have the effect that a minority of tags with very large (outlier) dispersions should have limited influence on the estimated value.
}

\references{
Cox, DR, and Reid, N (1987). Parameter orthogonality and approximate conditional inference.
\emph{Journal of the Royal Statistical Society Series B} 49, 1-39.

Robinson MD and Smyth GK (2008). Small-sample estimation of negative
binomial dispersion, with applications to SAGE data.
\emph{Biostatistics}, 9, 321-332
}

\author{Gordon Smyth}
\examples{
ntags <- 100
nlibs <- 4
y <- matrix(rnbinom(ntags*nlibs,mu=10,size=10),nrow=ntags,ncol=nlibs)
group <- factor(c(1,1,2,2))
lib.size <- rowSums(y)
design <- model.matrix(~group) # Define the design matrix for the full model
disp <- dispCoxReid(y, design, offset=log(lib.size), subset=100)
}

\seealso{
\code{\link{estimateGLMCommonDisp}}, \code{\link{optimize}}, \code{\link{uniroot}}
}

\keyword{models}

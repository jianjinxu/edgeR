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

\item{design}{numeric matrix giving the design matrix for the GLM that is to be fit.}

\item{offset}{numeric scalar, vector or matrix giving the offset (in addition to the log of the effective library size) that is to be included in the NB GLM for the transcripts. If a scalar, then this value will be used as an offset for all transcripts and libraries. If a vector, it should be have length equal to the number of libraries, and the same vector of offsets will be used for each transcript. If a matrix, then each library for each transcript can have a unique offset, if desired. In \code{adjustedProfileLik} the \code{offset} must be a matrix with the same dimension as the table of counts.}

\item{interval}{numeric vector of length 2 giving allowable values for the dispersion, passed to \code{optimize}.}

\item{tol}{the desired accuracy, see \code{optimize}.}

\item{min.row.sum}{integer. Only rows with at least this number of counts are used.}

\item{subset}{integer, number of rows to use in the calculation.  Rows used are chosen evenly spaced by abundance.}

\item{trace}{logical, should iteration information be output?}

\item{robust}{logical, should a robust estimator be used?}
}

\value{
Numeric vector of length one giving the estimated common dispersion.
}

\details{
In the \code{edgeR} context, these are low-level functions called by \code{estimateGLMCommonDisp}.

\code{dispCoxReid} maximizes the Cox-Reid adjusted profile likelihood (Cox and Reid, 1987).
\code{dispDeviance} and \code{dispPearson} set the deviance or Pearson goodness of fit statistics to their expected values.
\code{dispCoxReid} uses \code{optimize} while \code{dispDeviance} and \code{dispPearson} use \code{uniroot}.
The robust options mean that a minority of tags with very large (outlier) dispersions will not affect the estimated value.
}

\references{
Cox, DR, and Reid, N (1987). Parameter orthogonality and approximate conditional inference. \emph{Journal of the Royal Statistical Society Series B} 49, 1-39.
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

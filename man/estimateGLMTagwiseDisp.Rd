\name{estimateGLMTagwiseDisp}
\alias{estimateGLMTagwiseDisp}
\alias{estimateGLMTagwiseDisp.DGEList}
\alias{estimateGLMTagwiseDisp.default}

\title{Estimate Empirical Bayes Tagwise Dispersions for Negative Binomial GLMs}

\description{
Estimates the dispersion parameter for a DGE dataset for general experimental designs by using Cox-Reid approximate conditional inference for a negative binomial generalized linear model for each transcript (tag) with the unadjusted counts and design matrix provided.
}

\usage{
\S3method{estimateGLMTagwiseDisp}{DGEList}(y, design, offset=NULL, trend=!is.null(y$trended.dispersion), ...)
\S3method{estimateGLMTagwiseDisp}{default}(y, design, offset=NULL, dispersion, trend=TRUE, ...)
}

\arguments{ 

\item{y}{an object that contains the raw counts for each library (the measure of expression level); it can either be a matrix of counts, or a \code{DGEList} object with (at least) elements \code{counts} (table of unadjusted counts) and \code{samples} (data frame containing information about experimental group, library size and normalization factor for the library size)}

\item{design}{numeric matrix giving the design matrix for the GLM that is to be fit.}

\item{trend}{logical, should an abundance trend be applied to the grid of dispersion values over which the tagwise dispersion estimation is done?
Generally this should be \code{TRUE} if a trended dispersion has been estimated and \code{FALSE} otherwise.}

\item{offset}{numeric scalar, vector or matrix giving the offset (in addition to the log of the effective library size) that is to be included in the NB GLM for the transcripts. If a scalar, then this value will be used as an offset for all transcripts and libraries. If a vector, it should be have length equal to the number of libraries, and the same vector of offsets will be used for each transcript. If a matrix, then each library for each transcript can have a unique offset, if desired. Default is \code{NULL}; if object is a \code{DGEList} and offset is \code{NULL} then offset will be calculated automatically from code{y$samples}.}

\item{dispersion}{vector or scalar giving the dispersion value(s) to be used to set the grip of points for computation of the tagwise dispersion in \code{dispCoxReidInterpolateTagwise}.}

\item{\ldots}{other arguments are passed to lower-level functions. See \code{\link{dispCoxReidInterpolateTagwise}} for details.}
}

\value{
\code{estimateGLMTagwiseDisp.DGEList} produces a \code{DGEList} object, which contains the tagwise dispersion parameter estimate for each tag for the negative binomial model that maximizes the Cox-Reid adjusted profile likelihood. The tagwise dispersions are simply added to the \code{DGEList} object provided as the argument to the function.

\code{estimateGLMTagwiseDisp.default} returns a vector of the tagwise dispersion estimates.
}

\details{
This generic function is essentially a wrapper for \code{dispCoxReidInterpolateTagwise}.
To obtain estimates of the tagwise dispersion parameters for negative binomial GLMs we use Cox-Reid approximate conditional inference as implemented in \code{dispCoxReidInterpolateTagwise}.
The approach is to maximize the adjusted profile likelihood over the dispersion value, for the tagwise models and use these values as the tagwise dispersion parameters for differential signal testing in downstream analysis.
}

\references{
Cox, DR, and Reid, N (1987). Parameter orthogonality and approximate conditional inference. \emph{Journal of the Royal Statistical Society Series B} 49, 1-39.
}

\author{Gordon Smyth, Davis McCarthy}
\examples{
y <- matrix(rnbinom(1000,mu=10,size=10),ncol=4)
d <- DGEList(counts=y,group=c(1,1,2,2),lib.size=c(1000:1003))
design <- model.matrix(~group, data=d$samples) # Define the design matrix for the full model
d <- estimateGLMTrendedDisp(d, design, min.n=10)
d <- estimateGLMTagwiseDisp(d, design)
summary(d$tagwise.dispersion)
}

\seealso{
\code{\link{estimateGLMCommonDisp}} for common dispersion and \code{\link{estimateGLMTrendedDisp}} for trended dispersion in the context of a generalized linear model.

\code{\link{estimateCommonDisp}} for common dispersion or \code{\link{estimateTagwiseDisp}} for tagwise dispersion in the context of a multiple group experiment (one-way layout).
}

\keyword{algebra}

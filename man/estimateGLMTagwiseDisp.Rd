\name{estimateGLMTagwiseDisp}
\alias{estimateGLMTagwiseDisp}
\alias{estimateGLMTagwiseDisp.DGEList}
\alias{estimateGLMTagwiseDisp.default}

\title{Estimate Tagwise Dispersions for Negative Binomial GLMs}

\description{
Estimates the dispersion parameter for a DGE dataset for general experimental designs by using Cox-Reid approximate conditional inference for a negative binomial generalized linear model for each transcript (tag) with the unadjusted counts and design matrix provided.
}


\usage{
\S3method{estimateGLMTagwiseDisp}{DGEList}(y, design, method="trend", ...)
\S3method{estimateGLMTagwiseDisp}{default}(y, design, offset, dispersion, ...)
}

\arguments{ 

\item{y}{an object that contains the raw counts for each library (the measure of expression level); it can either be a matrix of counts, or a \code{DGEList} object with (at least) elements \code{counts} (table of unadjusted counts) and \code{samples} (data frame containing information about experimental group, library size and normalization factor for the library size)}

\item{design}{numeric matrix giving the design matrix for the GLM that is to be fit.}

\item{method}{method for defining the set up for smoothing the tagwise dispersion estimates towards some `common' value. If \code{trend} then the function uses the trended dispersion estimates from the \code{DGEList} object (\code{y$trended.dispersion}) to set the grid of points on which the tagwise dispersions are computed. If \code{common} then the function uses the common dispersion estimate from the \code{DGEList} object, (\code{y$common.dispersion}).}

\item{offset}{numeric scalar, vector or matrix giving the offset (in addition to the log of the effective library size) that is to be included in the NB GLM for the transcripts. If a scalar, then this value will be used as an offset for all transcripts and libraries. If a vector, it should be have length equal to the number of libraries, and the same vector of offsets will be used for each transcript. If a matrix, then each library for each transcript can have a unique offset, if desired.}

\item{dispersion}{vector or scalar giving the dispersion value(s) to be used to set the grip of points for computation of the tagwise dispersion in \code{dispCoxReidInterpolateTagwise}.}

\item{\ldots}{other arguments are passed to lower-level functions.}
}

\value{
\code{estimateGLMTagwiseDisp.DGEList} produces a \code{DGEList} object, which contains the tagwise dispersion parameter estimate for each tag for the negative binomial model that maximizes the Cox-Reid adjusted profile likelihood. The tagwise dispersions are simply added to the \code{DGEList} object provided as the argument to the function.

\code{estimateGLMTagwiseDisp.default} returns a  vector of the tagwise dispersion estimates.
}

\details{
This generic function is simply a wrapper for \code{dispCoxReidInterpolateTagwise}. To obtain estimates of the tagwise (i.e., genewise) dispersion parameters for negative binomial GLMs we use Cox-Reid approximate conditional inference as implemented in \code{dispCoxReidInterpolateTagwise}. The approach is to maximize the adjusted profile likelihood over the dispersion value, for the tagwise models and use these values as the tagwise dispersion parameters for differential signal testing in downstream analysis.
}


\references{
Cox, DR, and Reid, N (1987). Parameter orthogonality and approximate conditional inference. \emph{Journal of the Royal Statistical Society Series B} 49, 1-39.
}

\author{Gordon Smyth, Davis McCarthy}
\examples{
y <- matrix(rnbinom(1000,mu=10,size=10),ncol=4)
d <- DGEList(counts=y,group=c(1,1,2,2),lib.size=c(1000:1003))
design <- model.matrix(~group, data=d$samples) # Define the design matrix for the full model
d <- estimateGLMTrendedDisp(d, design)
d <- estimateGLMTagwiseDisp(d, design)
summary(d$tagwise.dispersion)
}

\seealso{
\code{\link{estimateGLMCommonDisp}}, \code{\link{estimateGLMTrendedDisp}}

\code{\link{estimateGLMCommonDisp}}, and \code{\link{estimateGLMTrendedDisp}} can be used to estimate a common value for the dispersion parameter for the whole dataset or a value for each tag/transcript that incorporates a mean-dependent trend on the dispersion, respectively, for a GLM (complicated experimental design). \code{\link{estimateCommonDisp}} and \code{\link{estimateTagwiseDisp}} use conditional maximum likelihood methods to estimate common and tagwise dispersions respectively, but are limited to multiple-group experimental designs (one-way layouts).
}

\keyword{algebra}

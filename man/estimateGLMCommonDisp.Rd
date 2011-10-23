\name{estimateGLMCommonDisp}
\alias{estimateGLMCommonDisp}
\alias{estimateGLMCommonDisp.DGEList}
\alias{estimateGLMCommonDisp.default}

\title{Estimate Common Dispersion for Negative Binomial GLMs}

\description{
Estimates a common negative binomial dispersion parameter for a DGE dataset with a general experimental design.
}

\usage{
\S3method{estimateGLMCommonDisp}{DGEList}(y, design=NULL, offset=NULL, method="CoxReid", ...)
\S3method{estimateGLMCommonDisp}{default}(y, design=NULL, offset=NULL, method="CoxReid", ...)
}

\arguments{ 

\item{y}{an object that contains the raw counts for each library (the measure of expression level); it can either be a matrix of counts, or a \code{DGEList} object with (at least) elements \code{counts} (table of unadjusted counts) and \code{samples} (data frame containing information about experimental group, library size and normalization factor for the library size)}

\item{design}{numeric matrix giving the design matrix for the GLM that is to be fit.
Must be of full column rank.
Defaults to a single column of ones, equivalent to treating the columns as replicate libraries.}

\item{method}{method for estimating the dispersion.
Possible values are \code{"CoxReid"}, \code{"Pearson"} or \code{"deviance"}.}

\item{offset}{numeric scalar, vector or matrix giving the offsets for the log-linear models.
If a scalar, then this value will be used as an offset for all transcripts and libraries.
If a vector, it should be have length equal to the number of libraries, and the same vector of offsets will be used for each transcript.
If a matrix, then it should have the same row and column dimensions as \code{y}.

In the \code{DGEList} method, the offset is calculated by default from the library sizes and normalization factors found in \code{y$samples}.}

\item{\ldots}{other arguments are passed to lower-level functions.
See \code{\link{dispCoxReid}}, \code{\link{dispPearson}} and \code{\link{dispDeviance}}
for details.}
}

\value{
The default method returns a numeric vector of length 1 containing the estimated dispersion.

The \code{DGEList} method returns the same \code{DGEList} \code{y} as input but with \code{common.dispersion} as an added component.
}

\details{
This function calls \code{dispCoxReid}, \code{dispPearson} or \code{dispDeviance} depending on the \code{method} specified.
See \code{\link{dispCoxReid}} for details of the three methods and a discussion of their relative performance.
}


\references{
Robinson MD and Smyth GK (2008). Small-sample estimation of negative
binomial dispersion, with applications to SAGE data.
\emph{Biostatistics}, 9, 321-332
}

\author{Gordon Smyth}
\examples{
#  True dispersion is 1/size=0.1
y <- matrix(rnbinom(1000,mu=10,size=10),ncol=4)
d <- DGEList(counts=y,group=c(1,1,2,2))
design <- model.matrix(~group, data=d$samples)
d1 <- estimateGLMCommonDisp(d, design)
d1$common.disp

#  Compare with classic CML estimator:
d2 <- estimateCommonDisp(d)
d2$common.disp

#  See example(glmFit) for a different example
}

\seealso{
\code{\link{dispCoxReid}}, \code{\link{dispPearson}}, \code{\link{dispDeviance}}

\code{\link{estimateGLMTrendedDisp}} for trended dispersion and \code{\link{estimateGLMTagwiseDisp}} for tagwise dispersions in the context of a generalized linear model.

\code{\link{estimateCommonDisp}} for common dispersion or \code{\link{estimateTagwiseDisp}} for tagwise dispersion in the context of a multiple group experiment (one-way layout).
}

\keyword{models}

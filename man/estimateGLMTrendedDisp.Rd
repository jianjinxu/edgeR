\name{estimateGLMTrendedDisp}
\alias{estimateGLMTrendedDisp}
\alias{estimateGLMTrendedDisp.DGEList}
\alias{estimateGLMTrendedDisp.default}

\title{Estimate Trended Dispersion for Negative Binomial GLMs}

\description{
Estimates the dispersion parameter for each transcript (tag) with a trend that depends on the overall level of expression for the transcript for a DGE dataset for general experimental designs by using Cox-Reid approximate conditional inference for a negative binomial generalized linear model for each transcript (tag) with the unadjusted counts and design matrix provided.
}


\usage{
\S3method{estimateGLMTrendedDisp}{DGEList}(y, design=NULL, offset=NULL, method="auto", ...)
\S3method{estimateGLMTrendedDisp}{default}(y, design=NULL, offset=NULL, method="auto",  ...)
}

\arguments{ 

\item{y}{an object that contains the raw counts for each library (the measure of expression level); it can either be a matrix of counts, or a \code{DGEList} object with (at least) elements \code{counts} (table of unadjusted counts) and \code{samples} (data frame containing information about experimental group, library size and normalization factor for the library size)}

\item{design}{numeric design matrix, as for \code{\link{glmFit}}.}

\item{method}{method (low-level function) used to estimated the trended dispersions.
Possible values are \code{"auto"} (default, switch to \code{"bin.spline"} method if the number of tags is great than 200 and \code{"power"} method otherwise),\code{"bin.spline"}, \code{"bin.loess"} (which both result in a call to \code{dispBinTrend}), \code{"power"} (call to \code{dispCoxReidPowerTrend}), or \code{"spline"} (call to \code{dispCoxReidSplineTrend}).}

\item{offset}{numeric scalar, vector or matrix giving the offset (in addition to the log of the effective library size) that is to be included in the NB GLM for the transcripts. If a scalar, then this value will be used as an offset for all transcripts and libraries. If a vector, it should be have length equal to the number of libraries, and the same vector of offsets will be used for each transcript. If a matrix, then each library for each transcript can have a unique offset, if desired. In \code{adjustedProfileLik} the \code{offset} must be a matrix with the same dimension as the table of counts. Default is \code{NULL}; if object is a \code{DGEList} and offset is \code{NULL} then offset will be calculated automatically from code{y$samples}.}

\item{\ldots}{other arguments are passed to lower-level functions.
See \code{\link{dispBinTrend}}, \code{\link{dispCoxReidPowerTrend}} and \code{\link{dispCoxReidSplineTrend}} for details.}
}

\value{
When the input object is a \code{DGEList}, \code{estimateGLMTrendedDisp} produces a \code{DGEList} object, which contains the estimates of the trended dispersion parameter for the negative binomial model according to the method applied.

When the input object is a numeric matrix, the output of one of the lower-level functions \code{dispBinTrend}, \code{dispCoxReidPowerTrend} of \code{dispCoxReidSplineTrend} is returned.
}

\details{
This is a wrapper function for the lower-level functions that actually carry out the dispersion estimation calculations. Provide a convenient, object-oriented interface for users.
}


\references{
Cox, DR, and Reid, N (1987). Parameter orthogonality and approximate conditional inference. \emph{Journal of the Royal Statistical Society Series B} 49, 1-39.

McCarthy, DJ, Chen, Y, Smyth, GK (2012). Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation.
\emph{Nucleic Acids Research} 40, 4288-4297.
\url{http://nar.oxfordjournals.org/content/40/10/4288}
}

\author{Gordon Smyth, Davis McCarthy, Yunshun Chen}
\examples{
y <- matrix(rnbinom(1000,mu=10,size=10),ncol=4)
d <- DGEList(counts=y,group=c(1,1,2,2),lib.size=c(1000:1003))
design <- model.matrix(~group, data=d$samples) # Define the design matrix for the full model
disp <- estimateGLMTrendedDisp(d, design, min.n=10)
}

\seealso{
\code{\link{dispBinTrend}}, \code{\link{dispCoxReidPowerTrend}} and \code{\link{dispCoxReidSplineTrend}} for details on how the calculations are done.

\code{\link{estimateGLMCommonDisp}} for common dispersion and \code{\link{estimateGLMTagwiseDisp}} for (trended) tagwise dispersion in the context of generalized linear models.

\code{\link{estimateCommonDisp}} for common dispersion or \code{\link{estimateTagwiseDisp}} for tagwise dispersion in the context of a multiple group experiment (one-way layout).
}

\keyword{models}

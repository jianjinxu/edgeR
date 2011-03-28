\name{estimateGLMCommonDisp}
\alias{estimateGLMCommonDisp}
\alias{estimateGLMCommonDisp.DGEList}
\alias{estimateGLMCommonDisp.default}

\title{Estimate Common Dispersion for Negative Binomial GLMs}

\description{
Estimates the dispersion parameter for a DGE dataset for general experimental designs by using Cox-Reid approximate conditional inference for a negative binomial generalized linear model for each transcript (tag) with the unadjusted counts and design matrix provided.
}


\usage{
\S3method{estimateGLMCommonDisp}{DGEList}(y, design, method="CoxReid", ...)
\S3method{estimateGLMCommonDisp}{default}(y, design, method="CoxReid", offset=NULL, ...)
}

\arguments{ 

\item{y}{an object that contains the raw counts for each library (the measure of expression level); it can either be a matrix of counts, or a \code{DGEList} object with (at least) elements \code{counts} (table of unadjusted counts) and \code{samples} (data frame containing information about experimental group, library size and normalization factor for the library size)}

\item{design}{numeric matrix giving the design matrix for the GLM that is to be fit.}

\item{method}{method (low-level function) used to estimated the dispersion.
Possible values are \code{"CoxReid"}, \code{"Pearson"} or \code{"deviance"}.}

\item{offset}{numeric scalar, vector or matrix giving the offset (in addition to the log of the effective library size) that is to be included in the NB GLM for the transcripts. If a scalar, then this value will be used as an offset for all transcripts and libraries. If a vector, it should be have length equal to the number of libraries, and the same vector of offsets will be used for each transcript. If a matrix, then each library for each transcript can have a unique offset, if desired. In \code{adjustedProfileLik} the \code{offset} must be a matrix with the same dimension as the table of counts.}

\item{\ldots}{other arguments are passed to lower-level functions.}
}

\value{
The default method returns a numeric vector of length 1.

The method for \code{DGEList} objects returns a , which contains the estimate of the common dispersion parameter for the negative binomial model that maximizes the Cox-Reid adjusted profile likelihood, and also the tagwise Cox-Reid dispersion estimates.

}

\details{
To obtain estimates of the common and tagwise (i.e., genewise) dispersion parameters for negative binomial GLMs we use Cox-Reid approximate conditional inference. The approach is to maximize the adjusted profile likelihood over the dispersion value, for both the common and tagwise models and use these values as the common and tagwise dispersion parameters for differential signal testing in downstream analysis.
}


\references{
Cox, DR, and Reid, N (1987). Parameter orthogonality and approximate conditional inference. \emph{Journal of the Royal Statistical Society Series B} 49, 1-39.
}

\author{Gordon Smyth}
\examples{
y <- matrix(rnbinom(1000,mu=10,size=10),ncol=4)
d <- DGEList(counts=y,group=c(1,1,2,2),lib.size=c(1000:1003))
design <- model.matrix(~group, data=d$samples) # Define the design matrix for the full model
disp <- estimateGLMCommonDisp(d, design)
}

\seealso{
\code{\link{dispCoxReid}}, \code{\link{dispPearson}}, \code{\link{dispDeviance}}

\code{\link{estimateGLMTrendedDisp}}, \code{\link{estimateGLMTagwiseDisp}}
}

\keyword{algebra}

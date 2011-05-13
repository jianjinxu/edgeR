\name{estimateCRDisp}
\alias{estimateCRDisp}

\title{Estimate Negative Binomial Dispersion By Cox-Reid Adjusted Profile Likelihood (Obsolete)}

\description{
This function has been replaced by \code{estimateGLMCommonDisp}, \code{estimateGLMTrendedDisp} and \code{estimateGLMTagwiseDisp}.
It is kept in this release of the package for backward capatability, but will be removed in future releases.

Estimates the dispersion parameter for a DGE dataset for general experimental designs by using Cox-Reid approximate conditional inference for a negative binomial generalized linear model for each transcript (tag) with the unadjusted counts and design matrix provided.
}


\usage{
estimateCRDisp(y, design=NULL, offset=NULL, npts=10, rowsum.filter=5, subset=10000, tagwise=FALSE, prior.n=10, min.n=500, lib.size=NULL, disp.trend=NULL, method.trend="binned-spline",verbose=TRUE)
}

\arguments{ 

\item{y}{an object that contains the raw counts for each library (the measure of expression level); it can either be a matrix of counts, or a \code{DGEList} object with (at least) elements \code{counts} (table of unadjusted counts) and \code{samples} (data frame containing information about experimental group, library size and normalization factor for the library size)}

\item{design}{numeric matrix giving the design matrix for the GLM that is to be fit.}

\item{offset}{numeric scalar, vector or matrix giving the offset (in addition to the log of the effective library size) that is to be included in the NB GLM for the transcripts. If a scalar, then this value will be used as an offset for all transcripts and libraries. If a vector, it should be have length equal to the number of libraries, and the same vector of offsets will be used for each transcript. If a matrix, then each library for each transcript can have a unique offset, if desired. In \code{adjustedProfileLik} the \code{offset} must be a matrix with the same dimension as the table of counts.}

\item{npts}{scalar, the number of points at which to place knots for the spline-based estimation of the common and tagwise dispersion estimates.}

\item{rowsum.filter}{numeric scalar giving a value for the filtering out of low abundance tags in the estimation of the common dispersion. Only tags with total sum of counts above this value are used in the estimation of the common dispersion. Low abundance tags can adversely affect the estimation of the common dispersion, so this argument allows the user to select an appropriate filter threshold for the tag abundance.}

\item{subset}{integer, number of rows to use in the calculation.  Rows used are chosen evenly spaced by abundance.}

\item{tagwise}{logical scalar, if \code{FALSE} (default) then the tagwise dispersions are not calculated, if \code{TRUE} then the tagwise dispersions are calculated.}

\item{prior.n}{numeric scalar, smoothing parameter that indicates the weight to give to the common likelihood compared to the individual tag's likelihood; default \code{10} means that the common likelihood is given 10 times the weight of the individual tag/gene's likelihood in the estimation of the tag/genewise dispersion}

\item{min.n}{integer, the minimum number of genes for each bin for the calculation of trended dispersions. Passed to \code{dispBinTrend}.}

\item{lib.size}{optional vector providing the (effective) library size for each library (must have length equal to the number of columns, or libraries, in the matrix of counts). If \code{NULL}, then a default is used. If \code{y} is a \code{DGEList} object then the default for \code{lib.size} is the product of the library sizes and the normalization factors (in the \code{samples} slot of the object). If \code{y} is a simple matrix of counts, then the default for \code{lib.size} is the vector of column sums of \code{y}.}

\item{disp.trend}{optional vector providing trended dispersion estimates if these have already been computed. If left as the default value, \code{NULL}, then trended dispersions are computed in \code{estimateCRDisp}.}

\item{method.trend}{character determining the method to be used to obtain trended dispersion estimates. Only required if \code{disp.trend=NULL}. If \code{"binned-spline"} or \code{"binned-loess"} then \code{dispBinTrend} is called to compute trended dispersions using a spline fit or a loess fit, respectively, through binned common dispersion values. If \code{"power"} then \code{dispCoxReidPowerTrend} is used to obtain trended dispersion estimates.}

\item{verbose}{logical scalar, if \code{TRUE} (default) then certain notification messages are displayed in some circumstances, if \code{FALSE} then these messages are not displayed.}
}

\value{\code{estimateCRDisp} produces a \code{DGEList} object, which contains the estimate of the common dispersion parameter for the negative binomial model that maximizes the Cox-Reid adjusted profile likelihood, and also the tagwise Cox-Reid dispersion estimates.
}

\details{
To obtain estimates of the common and tagwise (i.e., genewise) dispersion parameters for negative binomial GLMs we use Cox-Reid approximate conditional inference. The approach is to maximize the adjusted profile likelihood over the dispersion value, for both the common and tagwise models and use these values as the common and tagwise dispersion parameters for differential signal testing in downstream analysis.
}


\references{
Cox, DR, and Reid, N (1987). Parameter orthogonality and approximate conditional inference. \emph{Journal of the Royal Statistical Society Series B} 49, 1-39.
}

\author{Yunshun Chen, Gordon Smyth, Davis McCarthy}
\examples{
## This function is obsolete, but code below shows usage
y<-matrix(rnbinom(1000,mu=10,size=2),ncol=4)
d<-DGEList(counts=y,group=c(1,1,2,2),lib.size=c(1000:1003))
design <- model.matrix(~d$samples$group) # Define the design matrix for the full model
d<-estimateCRDisp(d, design, min.n=20)
d
}

\seealso{
This function is now obsolete, replaced by \code{\link{estimateGLMCommonDisp}}, \code{\link{estimateGLMTrendedDisp}} and \code{\link{estimateGLMTagwiseDisp}}.
}

\keyword{algebra}

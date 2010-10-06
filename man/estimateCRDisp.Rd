\name{estimateCRDisp}
\alias{estimateCRDisp}
\alias{adjustedProfileLik}
\alias{maximize.by.interpolation}

\title{Estimate the dispersion parameter for a negative binomial model using Cox-Reid approximate conditional inference}

\description{Estimates the common dispersion parameter for a DGE dataset for general experimental designs by using Cox-Reid approximate conditional inference for a negative binomial generalized linear model for each transcript (tag) with the unadjusted counts and design matrix provided. }


\usage{ 
estimateCRDisp(y, design=NULL, offset=0, npts=10, min.disp=0, max.disp=2, nselect=200, rowsum.filter=5, tagwise=FALSE, prior.n=10, trend=FALSE, lib.size=NULL, verbose=TRUE)
adjustedProfileLik(dispersion, y, design, offset)
maximize.by.interpolation(x,z,maxit=10,eps=1e-7,plot=FALSE)
}

\arguments{ 

\item{y}{an object that contains the raw counts for each library (the measure of expression level); it can either be a matrix of counts, or a \code{DGEList} object with (at least) elements \code{counts} (table of unadjusted counts) and \code{samples} (data frame containing information about experimental group, library size and normalization factor for the library size)}

\item{design}{numeric matrix giving the design matrix for the GLM that is to be fit.}

\item{offset}{numeric scalar, vector or matrix giving the offset (in addition to the log of the effective library size) that is to be included in the NB GLM for the transcripts. If a scalar, then this value will be used as an offset for all transcripts and libraries. If a vector, it should be have length equal to the number of libraries, and the same vector of offsets will be used for each transcript. If a matrix, then each library for each transcript can have a unique offset, if desired. In \code{adjustedProfileLik} the \code{offset} must be a matrix with the same dimension as the table of counts.}

\item{npts}{scalar, the number of points at which to place knots for the spline-based estimation of the common and tagwise dispersion estimates.}

\item{min.disp}{scalar, the minimum possible value for the dispersion. May need to be set smaller (e.g. 1e-04 or less) if there is no biological variability in the data.}

\item{max.disp}{scalar, the maximum possible value for the dispersion.}

\item{nselect}{scalar, the number of genes/tags to be used to get an initial 'ballpark' estimate of the magnitude of the dispersions in the dataset. Used to finesse the calculation of the estimates using all the data.}

\item{rowsum.filter}{numeric scalar giving a value for the filtering out of low abundance tags in the estimation of the common dispersion. Only tags with total sum of counts above this value are used in the estimation of the common dispersion. Low abundance tags can adversely affect the estimation of the common dispersion, so this argument allows the user to select an appropriate filter threshold for the tag abundance.}

\item{tagwise}{logical scalar, if \code{FALSE} (default) then the tagwise dispersions are not calculated, if \code{TRUE} then the tagwise dispersions are calculated.}

\item{prior.n}{numeric scalar, smoothing parameter that indicates the weight to give to the common likelihood compared to the individual tag's likelihood; default \code{10} means that the common likelihood is given 10 times the weight of the individual tag/gene's likelihood in the estimation of the tag/genewise dispersion}

\item{trend}{logical scalar, if \code{FALSE} (default) then the abundance-dispersion trend is not considered in calculating both the common dispersion and the tagwise dispersions, if \code{TRUE} then such trend is introduced in calculating both dispersions.}

\item{lib.size}{optional vector providing the (effective) library size for each library (must have length equal to the number of columns, or libraries, in the matrix of counts). If \code{NULL}, then a default is used. If \code{y} is a \code{DGEList} object then the default for \code{lib.size} is the product of the library sizes and the normalization factors (in the \code{samples} slot of the object). If \code{y} is a simple matrix of counts, then the default for \code{lib.size} is the vector of column sums of \code{y}.}

\item{verbose}{logical scalar, if \code{TRUE} (default) then certain notification messages are displayed in some circumstances, if \code{FALSE} then these messages are not displayed.}

\item{dispersion}{numeric scalar providing the common value for the dispersion parameter (the 'size' parameter in the GLM fit is equal to \code{1/dispersion}) that is used in fitting the GLM for each transcript. Poisson GLM is fitted if \code{dispersion} is set at 0. \code{estimateCRDisp} maximizes the Cox-Reid adjusted profile likelihood over dispersion to obtain the estimate for the common dispersion.}

\item{maxit}{scalar, the maximum number of iterations to perform in the spline interpolation.}

\item{eps}{scalar, the small 'epsilon' value that provides the tolerance for the convergence of the algorithm.}

\item{x}{vector of independent variables for the interpolation, i.e. return the interpolated x value that gives the maximum z value.}

\item{z}{vector of dependent variables for the interpolation, i.e. interpolate the maximum over the z values.}

\item{plot}{logical, whether or not to produce a plot of the maximization by interpolation.}

}

\value{ \code{estimateCRDisp} produces a \code{DGEList} object, which contains the estimate of the common dispersion parameter for the negative binomial model that maximizes the Cox-Reid adjusted profile likelihood, and also the tagwise Cox-Reid dispersion estimates.

\code{adjustedProfileLik} produces a vector of the tagwise Cox-Reid adjusted profile likelihood for the given counts, dispersion  value, offset and design matrices (i.e. the APL for each gene/tag).

\code{maximize.by.interpolation} returns a scalar value, which is the value for x that maximizes z.

}

\details{
To obtain estimates of the common and tagwise (i.e. genewise) dispersion parameters for negative binomial GLMs we use Cox-Reid approximate conditional inference. The approach is to maximize the adjusted profile likelihood over the dispersion value, for both the common and tagwise models and use these values as the common and tagwise dispersion parameters for differential signal testing in downstream analysis.
}


\references{
Cox DR and Reid N (1987). Parameter orthogonality and approximate conditional inference. \emph{Journal of the Royal Statistical Society. Series B (Methodological)},
49, 1, 1--39.
}

\author{Yunshun Chen, Gordon Smyth}
\examples{
y<-matrix(rnbinom(1000,mu=10,size=2),ncol=4)
d<-DGEList(counts=y,group=c(1,1,2,2),lib.size=c(1000:1003))
design <- model.matrix(~d$samples$group) # Define the design matrix for the full model
d<-estimateCRDisp(d, design)
d
}

\seealso{
\code{\link{estimateTagwiseDisp}}, and \code{\link{estimateCommonDisp}} can be used to estimate a value for the dispersion parameter for each tag/transcript and a common dispersion value, respectively. The estimates are stabilized by squeezing the estimates towards the common value calculated by \code{estimateCommonDisp}. These functions use exact conditional methods, but are restricted to less complicated experimental designs; they can deal with multiple groups, but nothing more complicated.

}

\keyword{algebra}

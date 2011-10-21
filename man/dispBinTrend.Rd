\name{dispBinTrend}
\alias{dispBinTrend}

\title{Estimate Dispersions with an Abundance-Dependent Trend for Negative Binomial GLMs}

\description{
Estimate a dispersion parameter for each of many negative binomial generalized linear models by computing the common dispersion for genes sorted into bins based on overall abundance and then using splines or a loess fit to interpolate a dispersion value for each gene, dependent on overall abundance of the gene.
}

\usage{
dispBinTrend(y, design, offset=NULL, df = 5, span=2/3, min.n=500, method.bin="CoxReid", method.trend="spline", trace=0, ...)
}

\arguments{ 

\item{y}{numeric matrix of counts}

\item{design}{numeric matrix giving the design matrix for the GLM that is to be fit.}

\item{offset}{numeric scalar, vector or matrix giving the offset (in addition to the log of the effective library size) that is to be included in the NB GLM for the transcripts. If a scalar, then this value will be used as an offset for all transcripts and libraries. If a vector, it should be have length equal to the number of libraries, and the same vector of offsets will be used for each transcript. If a matrix, then each library for each transcript can have a unique offset, if desired. In \code{adjustedProfileLik} the \code{offset} must be a matrix with the same dimension as the table of counts.}

\item{df}{scalar, the degrees of freedom for the natural cubic splines fit, used to determine the placement of the knots (number of knots is \code{df - 1}. Default is \code{5}.}

\item{span}{scalar, passed to \code{loess} to determine the amount of smoothing for the loess fit. Default is \code{2/3}.}

\item{min.n}{scalar, minimim number of genes in each of the bins into which genes are sorted to form the basis for interpolating the dispersions. Setting a minimum value ensures that there will be sufficient genes in each bin to allow reliable estimation of the common dispersion for each bin.}

\item{method.bin}{character, passed to \code{binGLMDispersion}, to specify the method used to compute the common dispersion within each bin of genes. Default is \code{"CoxReid"}, other options are \code{"Pearson"} and \code{"deviance"}.}

\item{method.trend}{character, specifies method to produce a smooth fit through the binned common dispersions in order to interpolate the trended dispersions. Default is \code{"spline"} to use natural cubic splines, other option is \code{"loess"} to use a loess fit.}

\item{trace}{logical, should iteration information be output?}

\item{\ldots}{option arguments to be passed to lower-level function \code{binGLMDispersion}.}
}

\value{
list with the following components:
\item{abundance}{numeric vector containing the overall abundance for each gene}
\item{dispersion}{numeric vector giving the trended dispersion estimate for each gene}
\item{bin.abundance}{numeric vector of length equal to \code{nbins} giving the average (mean) abundance for each bin}
\item{bin.dispersion}{numeric vector of length equal to \code{nbins} giving the estimated common dispersion for each bin}
}

\details{
This function takes the binned common dispersion and abundance from \code{\link{binGLMDispersion}} and fits a smooth curve through these binned values using either natural cubic splines or loess. From this smooth curve it predicts the dispersion value for each gene based on the gene's overall abundance. This results in estimates for the NB dispersion parameter which have a dependence on the overall expression level of the gene, and thus have an abundance-dependent trend. This function is called by \code{\link{estimateGLMTrendedDisp}}.
}

\references{
Cox, DR, and Reid, N (1987). Parameter orthogonality and approximate conditional inference. \emph{Journal of the Royal Statistical Society Series B} 49, 1-39.
}

\author{Davis McCarthy and Gordon Smyth}
\examples{
ntags <- 1000
nlibs <- 4
means <- seq(5,10000,length.out=ntags)
y <- matrix(rnbinom(ntags*nlibs,mu=rep(means,nlibs),size=0.1*means),nrow=ntags,ncol=nlibs)
keep <- rowSums(y) > 0
y <- y[keep,]
group <- factor(c(1,1,2,2))
lib.size <- colSums(y)
design <- model.matrix(~group) # Define the design matrix for the full model
disp <- dispBinTrend(y, design, offset=log(lib.size), min.n=100, span=0.3)
plot(disp$abundance, disp$dispersion)
}

\seealso{
\code{\link{binGLMDispersion}}, \code{\link{estimateGLMTrendedDisp}}
}

\keyword{models}

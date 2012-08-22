\name{bin.dispersion}
\alias{binCMLDispersion}
\alias{binGLMDispersion}

\title{Estimate Common Dispersion for Negative Binomial GLMs in Bins of Genes Sorted by Overall Abundance}

\description{
Estimates the common dispersion parameter for each of a number of bins of data for a DGE dataset. Genes are sorted into bins based on overall expression level. For multiple-group (one-way layout) experimental designs, conditional maximum likelihood (CML) methods can be used. For general experimental designs the binned common dispersions we can use Cox-Reid approximate conditional inference, Pearson or deviance estimators for a negative binomial generalized linear model.
}

\usage{
binCMLDispersion(y, nbins=50)
binGLMDispersion(y, design, offset=NULL, min.n=100, method="CoxReid", abundance=NULL, ...)
}

\arguments{
\item{y}{an object that contains the raw counts for each library (the measure of expression level); it can either be a matrix of counts, or a \code{DGEList} object with (at least) elements \code{counts} (table of unadjusted counts) and \code{samples} (data frame containing information about experimental group, library size and normalization factor for the library size)}

\item{nbins}{scalar, the number of bins for which to compute common dispersions. Default is \code{50} bins.}

\item{design}{numeric matrix giving the design matrix for the GLM that is to be fit.}

\item{min.n}{scalar, the minimum number of genes to be included in each bin.}

\item{offset}{(optional) numeric scalar, vector or matrix giving the offset (in addition to the log of the effective library size) that is to be included in the NB GLM for the transcripts. If a scalar, then this value will be used as an offset for all transcripts and libraries. If a vector, it should be have length equal to the number of libraries, and the same vector of offsets will be used for each transcript. If a matrix, then each library for each transcript can have a unique offset, if desired. Default is \code{NULL}. If \code{NULL}, then offset is \code{log(lib.size)} if \code{y} is a matrix or \code{log(y$samples$lib.size * y$samples$norm.factors)} if \code{y} is a \code{DGEList} object.}

\item{method}{method  used to estimated the dispersion. Argument passed to \code{\link{estimateGLMCommonDisp}}, which calls the functions to do the computations. Possible values are \code{"CoxReid"}, \code{"Pearson"} or \code{"deviance"}.}

\item{abundance}{numeric vector giving abundance of each gene}

\item{\ldots}{other arguments are passed to lower-level functions.}
}

\value{
Returns a list with two components:
	\item{dispersion}{numeric vector providing the common dispersion for each bin}
	\item{abundance}{numeric vector providing the average abundance (expression level) of genes in each bin}
}

\details{
To obtain estimates of the common dispersion parameters conditional maximum likelihood (\code{\link{estimateCommonDisp}}) is used for \code{binCMLDispersion} and one of Cox-Reid approximate conditional inference (\code{\link{dispCoxReid}}), the deviance (\code{\link{dispDeviance}}) or Pearson (\code{\link{dispPearson}}) estimates are used for \code{binGLMDispersion}. Genes are assigned to bins using the \code{\link{cutWithMinN}} function to obtain bins spread over the abundance range of the genes while ensuring that each bin has a minimum number of genes, thus permitting reliable estimation of the common dispersion for each bin.

If there are fewer than \code{min.n} rows of \code{y}, then one bin is used.
The number of bins is limited to 1000.
}

\references{
Cox, DR, and Reid, N (1987). Parameter orthogonality and approximate conditional inference. \emph{Journal of the Royal Statistical Society Series B} 49, 1-39.
}

\author{Gordon Smyth, Davis McCarthy}
\examples{
y <- matrix(rnbinom(1000,mu=10,size=10),ncol=4)
d <- DGEList(counts=y,group=c(1,1,2,2),lib.size=c(1000:1003))
design <- model.matrix(~group, data=d$samples) # Define the design matrix for the full model
bindisp.CML <- binCMLDispersion(d, nbins=50)
bindisp.GLM <- binGLMDispersion(d, design, min.n=10)
}

\seealso{
\code{\link{estimateGLMCommonDisp}}, \code{\link{dispCoxReid}}, \code{\link{dispPearson}}, \code{\link{dispDeviance}}
}

\keyword{algebra}

\name{estimateDisp}
\alias{estimateDisp}

\title{Estimate Common, Trended and Tagwise Negative Binomial dispersions by weighted 
likelihood empirical Bayes}

\description{
Maximizes the negative binomial likelihood to give the estimate of the common, trended 
and tagwise dispersions across all tags.
}

\usage{
estimateDisp(y, design=NULL, prior.df=NULL, trend.method="locfit", span=NULL, 
grid.length=21, grid.range=c(-10,10), robust=FALSE, winsor.tail.p=c(0.05,0.1), tol=1e-06)
}

\arguments{
\item{y}{\code{DGEList} object}

\item{design}{numeric design matrix}

\item{prior.df}{prior degrees of freedom. It is used in calculating \code{prior.n}.}

\item{trend.method}{method for estimating dispersion trend. Possible values are 
\code{"none"}, \code{"movingave"}, \code{"loess"} and \code{"locfit"}.}

\item{span}{width of the smoothing window, as a proportion of the data set.}

\item{grid.length}{the number of points on which the interpolation is applied for each 
tag.}

\item{grid.range}{the range of the grid points around the trend on a log2 scale.}

\item{robust}{logical, should the estimation of \code{prior.df} be robustified against 
outliers?}

\item{winsor.tail.p}{numeric vector of length 1 or 2, giving left and right tail 
proportions of the deviances to Winsorize when estimating \code{prior.df}.}

\item{tol}{the desired accuracy, passed to \code{\link{optimize}}}
}

\value{Returns \code{object} with the following added components:
	\item{common.dispersion}{estimate of the common dispersion.}
	\item{trended.dispersion}{estimates of the trended dispersions.}
	\item{tagwise.dispersion}{tag- or gene-wise estimates of the dispersion 
parameter.}
	\item{logCPM}{the tag abundance in log average counts per million.}
	\item{prior.df}{prior degrees of freedom. It is a vector when robust method is 
used.}
	\item{prior.n}{estimate of the prior weight, i.e. the smoothing parameter that 
indicates the weight to put on the common likelihood compared to the individual tag's 
likelihood.}
	\item{span}{width of the smoothing window used in estimating dispersions.}
}

\details{
This function calculates a matrix of likelihoods for each gene at a set of dispersion 
grid points, and then applies weighted likelihood empirical Bayes method to obtain 
posterior dispersion estimates. If there is no design matrix, it calculates the quantile 
conditional likelihood for each gene (tag) and then maximize it. The method is same as in 
the function \code{estimateCommonDisp} and \code{estimateTagwiseDisp}. If a design matrix 
is given, it then calculates the adjusted profile log-likelihood for each gene (tag) and 
then maximize it. It is similar to the functions \code{estimateGLMCommonDisp}, 
\code{estimateGLMTrendedDisp} and \code{estimateGLMTagwiseDisp}.
}

\references{
Chen, Y, Lun, ATL, and Smyth, GK (2014).
Differential expression analysis of complex RNA-seq experiments using edgeR.
In: \emph{Statistical Analysis of Next Generation Sequence Data},
Somnath Datta and Daniel S Nettleton (eds), Springer, New York.
\url{http://www.statsci.org/smyth/pubs/edgeRChapterPreprint.pdf}
}

\author{Yunshun Chen, Gordon Smyth}
\examples{
# True dispersion is 1/5=0.2
y <- matrix(rnbinom(1000, mu=10, size=5), ncol=4)
group <- c(1,1,2,2)
design <- model.matrix(~group)
d <- DGEList(counts=y, group=group)
d1 <- estimateDisp(d)
d2 <- estimateDisp(d, design)
}

\seealso{
\code{\link{estimateCommonDisp}}, \code{\link{estimateTagwiseDisp}}, 
\code{\link{estimateGLMCommonDisp}}, \code{\link{estimateGLMTrendedDisp}}, 
\code{\link{estimateGLMTagwiseDisp}}
}
